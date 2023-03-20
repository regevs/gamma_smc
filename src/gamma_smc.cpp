#include "common.h"

#include "sys.h"
#include "io.h"
#include "flow_field.h"
#include "gamma_smc.h"
#include "data_processor.h"
#include "screenoutput.h"

int main(int argc, char** argv) {
    ScreenOutput screen;

    //
    // Print command line
    //
    screen.print_subtitle("Command line:");
    for (int i = 0; i < argc; ++i) {
        cout << argv[i] << ' ';
    }
    cout << endl << endl;


    //
    // Parse flags
    //
    po::options_description desc("Allowed options");
    desc.add_options()
        ("scaled_mutation_rate,m", po::value<float>()->required(), "Scaled mutation rate")
        ("scaled_recombination_rate,r", po::value<float>()->required(), "Scaled recombination rate")
        ("flow_field,f", po::value<string>()->default_value("resources/default_flow_field.txt")->required(), "Flow field file")
        ("input,i", po::value<string>()->required(), "Input file")
        ("mask,a", po::value<string>(), "File of a global mask (empty for no mask)")
        ("masks_per_sample,b", po::value<string>(), "File of masks filenames per sample (empty for no masks)")
        ("output_text", po::value<string>(), "Output file in text format (recommended only for small datasets!)")
        ("output,o", po::value<string>(), "Output file")
        ("samples,S", po::value<string>(), "Filename of a list of subset of samples to take")
        ("samples_against,T", po::value<string>(), "Filename of a second list of subset of samples to take, to infer against first list")
        ("only_within,w", "Apply only to haplotype pairs within each diploid")
        ("output_at_stride,s", po::value<int>()->default_value(-1), "Output at positions which are multiples of this number")
        ("output_at_hets,h", po::value<bool>()->default_value(true), "Output at segregating sites")
        ("cache_size,z", po::value<int>()->default_value(1000), "Maximum cache size in basepairs")
        ("only_forward,y", "Calculate only forward pass")
        ("only_backward,d", "Calculate only backward pass")
        ("help", "Produce help message")
    ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

    // Check for --help or --version here, before notify
    if (vm.count("help")) {
        std::cout << desc << "\n";
        exit(-1);
    }
    
    po::notify(vm);

    //
    // Validate flags
    //
    float scaled_mutation_rate = vm["scaled_mutation_rate"].as<float>();
    if (scaled_mutation_rate < 0) {
        cout << "Error: --scaled_mutation_rate must be positive." << endl;
        exit(-1);
    }

    float scaled_recombination_rate = vm["scaled_recombination_rate"].as<float>();
    if (scaled_recombination_rate < 0) {
        cout << "Error: --scaled_recombination_rate must be positive." << endl;
        exit(-1);
    }

    string flow_field_filename = vm["flow_field"].as<string>();
    if (!boost::filesystem::exists(flow_field_filename)) {
        cout << boost::format("Error: Cannot open --flow_field file: %s\n") % flow_field_filename;
        exit(-1);
    }

    string input_filename = vm["input"].as<string>();
    if (!boost::filesystem::exists(input_filename)) {
        cout << boost::format("Error: Cannot open --input file: %s\n") % input_filename;
        exit(-1);
    }

    auto output_filename = vm["output"].as<string>();
    auto output_directory = boost::filesystem::path(output_filename).parent_path();
    boost::filesystem::create_directories(output_directory);
    
    if (vm.count("only_within") && vm.count("samples_against")) {
        cout << "Error: --only_within and --samples_file_against are mutually exclusive." << endl;
        exit(-1);
    }
    if (vm.count("samples_against") && !vm.count("samples")) {
        cout << "Error: --samples_file_against requires --samples_file." << endl;
        exit(-1);
    }
    if (vm.count("mask") && vm.count("masks_per_sample")) {
        cout << "Error: --mask and --masks_per_sample are mutually exclusive." << endl;
        exit(-1);
    }

    string mask_filename;
    if (vm.count("mask")) {
        mask_filename = vm["mask"].as<string>();
        if (!boost::filesystem::exists(mask_filename)) {
            cout << boost::format("Error: Cannot open --mask file: %s\n") % mask_filename;
            exit(-1);
        }
    }

    string masks_per_sample_filename;
    if (vm.count("masks_per_filename")) {
        masks_per_sample_filename = vm["masks_per_filename"].as<string>();
        if (!boost::filesystem::exists(masks_per_sample_filename)) {
            cout << boost::format("Error: Cannot open --masks_per_filename file: %s\n") % masks_per_sample_filename;
            exit(-1);
        }
    }

    string output_text_filename;
    if (vm.count("output_text")) {
        output_text_filename = vm["output_text"].as<string>();
        auto output_text_directory = boost::filesystem::path(output_text_filename).parent_path();
        boost::filesystem::create_directories(output_text_directory);
    }

    string samples_filename;
    if (vm.count("samples")) {
        samples_filename = vm["samples"].as<string>();
        if (!boost::filesystem::exists(samples_filename)) {
            cout << boost::format("Error: Cannot open --samples file: %s\n") % samples_filename;
            exit(-1);
        }
    }

    string samples_against_filename;
    if (vm.count("samples_against")) {
        samples_against_filename = vm["samples_against"].as<string>();
        if (!boost::filesystem::exists(samples_against_filename)) {
            cout << boost::format("Error: Cannot open --samples_against file: %s\n") % samples_against_filename;
            exit(-1);
        }
    }

    int output_at_stride = vm["output_at_stride"].as<int>();
    bool output_at_hets = (vm.count("output_at_hets") > 0) && (vm["output_at_hets"].as<bool>());
    if (!output_at_hets && (output_at_stride == -1)) {
        cout << "Warning: No output flags provided.\n";
    }

    bool only_forward = (vm.count("only_forward") > 0);
    bool only_backward = (vm.count("only_backward") > 0);

    int cache_size = vm["cache_size"].as<int>();
    if (cache_size <= 0) {
        cout << boost::format("Error: --cache_size must be positive.\n");
        exit(-1);
    }

    //
    // Load input file
    //
    screen.print_subtitle("Reading input file...");

    vector<unique_ptr<SegregatingSite>> input_sites;  // TODO: Have another go at getting rid of the unique_ptr and just have vector<Seg...>
    vector<string> sample_names;    
    vector<int> samples_indices;
    vector<int> samples_against_indices;

    readVcf(
        input_filename, 
        input_sites, 
        sample_names, 
        samples_filename, 
        samples_indices,
        samples_against_filename,
        samples_against_indices
    ); 

    screen.print_item(boost::str(boost::format("Read %d samples.") % sample_names.size()));
    screen.print_item(boost::str(boost::format("Read %d segregating sites.") % input_sites.size()));
    screen.print_done();

    //
    // Load masks
    //
    screen.print_subtitle("Reading mask(s) file...");

    vector<pair<int, int>> global_mask;
    if (vm.count("mask")) {
        readMask(mask_filename, global_mask);
    } else {
        // If no global mask is given, assume no mask
        global_mask.push_back(make_pair(0, input_sites.back()->pos+1));
    }

    unordered_map<string, vector<pair<int, int>>> mask_map;
    if (vm.count("masks_per_sample")) {
        readMasks(masks_per_sample_filename, mask_map, sample_names);
        screen.print_item(boost::str(boost::format("Read %d masks.") % mask_map.size()));
    }
    
    screen.print_done();

    //
    // Process segments
    //
    screen.print_subtitle("Creating segments...");

    DataProcessor data_processor(
        input_sites,
        sample_names,
        global_mask,
        mask_map,
        output_at_stride,
        cache_size,
        output_at_hets
    );

    screen.print_item(boost::str(
        boost::format("Created %d segments, with %d output positions.") % data_processor._n_segments % data_processor._seq_length
    ));

    
    //
    // Create a list of pairs to work on
    //
    vector<pair<int, int>> haplotype_pairs;
    int n_haplotypes = sample_names.size() * 2;

    if (vm.count("only_within")) {
        for (uint i = 0; i < sample_names.size(); i++) {
            haplotype_pairs.push_back(make_pair(2*i, 2*i+1));
        }
    } else {
        if (samples_against_indices.size()) {
            for (int i : samples_indices) {
                for (int j : samples_against_indices) {
                    haplotype_pairs.push_back(make_pair(i, j));
                }
            }
        } else {
            for (int i = 0; i < n_haplotypes; i++) {
                for (int j = i+1; j < n_haplotypes; j++) {
                    haplotype_pairs.push_back(make_pair(i, j));
                }
            }
        }
    }

    screen.print_item(boost::str(
        boost::format("Applying to %d haplotype pairs") % haplotype_pairs.size()
    ));
    screen.print_done();


    //
    // Load flow field file
    //
    vector<float> mean_grid_def;
    vector<float> cv_grid_def;
    vector<float> flow_field_unravelled;

    read_flow_field_raw(
        flow_field_filename,
        mean_grid_def,
        cv_grid_def,
        flow_field_unravelled
    );

    //
    // Prepare output files
    //
    // TODO: Check errors
    ostream* out = NULL;
    ofstream* output_file = NULL;
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf;
    if (vm.count("output_text")) {        
        output_file = new ofstream(output_text_filename, ios_base::out | ios_base::binary);        
        outbuf.push(boost::iostreams::gzip_compressor());
        outbuf.push(*output_file);
        out = new ostream(&outbuf);
    }

    ostream* out_raw = NULL;
    ofstream* output_file_raw_meta = NULL;
    ofstream* output_file_raw = NULL;
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_raw;
    if (output_filename.size() > 0) {
        output_file_raw = new ofstream(output_filename, ios_base::out | ios_base::binary);        
        outbuf_raw.push(boost::iostreams::gzip_compressor());
        outbuf_raw.push(*output_file_raw);
        out_raw = new ostream(&outbuf_raw);

        output_file_raw_meta = new ofstream(output_filename + ".meta", ios_base::out);        
    }

    //
    // Construct flow field
    //
    unique_ptr<FlowField> FF(new FlowField(mean_grid_def, cv_grid_def, flow_field_unravelled));

    //
    // Construct flow field cache
    //
    screen.print_subtitle("Building flow field cache...");
    
    unique_ptr<FlowFieldCache> FFC(new FlowFieldCache(
        scaled_recombination_rate,
        scaled_mutation_rate,
        move(FF),
        cache_size,
        true                // entropy clipping
    ));

    screen.print_done();

    //
    // Run
    //
    screen.print_subtitle("Running...");
    CachedPairwiseGammaSMC PPC(
        input_sites,
        haplotype_pairs,
        scaled_recombination_rate,
        scaled_mutation_rate,
        move(FFC),
        data_processor,
        output_at_stride,
        output_at_hets,
        only_forward,
        only_backward,
        out,
        output_file_raw_meta,
        out_raw
    );

    PPC.calculate_posteriors();

    //
    // Close output files
    //
    if (out != NULL) {
        boost::iostreams::close(outbuf); // Don't forget this!
        output_file->close();
    }

    if (out_raw != NULL) {
        boost::iostreams::close(outbuf_raw); // Don't forget this!
        output_file_raw->close();
        output_file_raw_meta->close();
    }

    double total_processing_time = PPC._timer_emissions + PPC._timer_forward + PPC._timer_backward;  // Excludes output time
    double total_basepairs = data_processor._segments.back().pos + data_processor._segments.back().length - data_processor._segments.front().pos;
    double time_per_bp_per_pair = total_processing_time / total_basepairs / haplotype_pairs.size();
    double time_per_segment_per_pair = total_processing_time / data_processor._segments.size() / haplotype_pairs.size();
    

    screen.print_item(boost::str(boost::format("Emissions preparation time:\t%.3f secs") % PPC._timer_emissions));
    screen.print_item(boost::str(boost::format("Forward pass time:\t\t%.3f secs") % PPC._timer_forward));
    screen.print_item(boost::str(boost::format("Backward pass time:\t\t%.3f secs") % PPC._timer_backward));
    screen.print_item(boost::str(boost::format("Output time:\t\t%.3f secs") % PPC._timer_output));
    screen.print_done();

    screen.print_subtitle("Summary:");

    screen.print_item(boost::str(boost::format("Time per Gbp per pair\t%.3f secs") % (time_per_bp_per_pair * 1e9)));
    // screen.print_item(boost::str(boost::format("Time per 1M segments per pair\t%.3f secs") % (time_per_segment_per_pair * 1e6)));
    screen.print_item(boost::str(boost::format("Overall time:\t\t%.3f secs") % mp_cputime()));
    screen.print_item(boost::str(boost::format("Peak memory:\t\t%.3f GB") % (mp_peakrss() / 1024.0 / 1024.0 / 1024.0)));	
    screen.print_done();

    return 0;
}
