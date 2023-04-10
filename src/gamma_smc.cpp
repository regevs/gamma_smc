#include "common.h"

#include "sys.h"
#include "io.h"
#include "flow_field.h"
#include "gamma_smc.h"
#include "data_processor.h"
#include "screenoutput.h"
#include "cxxopts.hpp"

int main(int argc, char** argv) {
    ScreenOutput screen;
    screen.print_header();

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
    cxxopts::Options options("Gamma-SMC", "Fast inference of pairwise coalescence times");

    options.add_options()
        ("i,input", "Input file", cxxopts::value<std::string>())
        ("o,output", "Output file", cxxopts::value<std::string>())
        ("m,scaled_mutation_rate", "Scaled mutation rate", cxxopts::value<float>())
        ("r,scaled_recombination_rate", "Scaled recombination rate", cxxopts::value<float>())
        ("t,recombination_to_mutation_ratio", "Recombination to mutation rates ratio", cxxopts::value<float>())
        ("f,flow_field", "Flow field file", cxxopts::value<std::string>())
        ("a,mask", "File of a global mask (empty for no mask)", cxxopts::value<std::string>())
        ("b,masks_per_sample", "File of masks filenames per sample (empty for no masks)", cxxopts::value<std::string>())
        ("S,samples", "Filename of a list of subset of samples to take", cxxopts::value<std::string>())
        ("T,samples_against", "Filename of a second list of subset of samples to take, to infer against first list", cxxopts::value<std::string>())
        ("w,only_within", "Apply only to haplotype pairs within each diploid")
        ("s,output_at_stride", "Output at positions which are multiples of this number", cxxopts::value<int>()->default_value("-1"))
        ("h,output_at_hets", "Output at segregating sites", cxxopts::value<bool>()->default_value("true"))
        ("z,cache_size", "Maximum cache size in basepairs", cxxopts::value<int>()->default_value("1000"))
        ("y,only_forward", "Calculate only forward pass", cxxopts::value<bool>()->default_value("false"))
        ("d,only_backward", "Calculate only backward pass", cxxopts::value<bool>()->default_value("false"))
        ("zstd_compression_level", "zstd compression level", cxxopts::value<int>()->default_value("1"))
        ("help", "Produce help message", cxxopts::value<bool>()->default_value("false"))
    ;

    auto vm = options.parse(argc, argv);

    // Check for --help or --version here, before notify
    if (vm.count("help")) {
        std::cout << options.help() << "\n";
        exit(-1);
    }
    
    //
    // Validate flags
    //
    float scaled_mutation_rate = -1;
    if (vm.count("scaled_mutation_rate")) {
        scaled_mutation_rate = vm["scaled_mutation_rate"].as<float>();
        if (scaled_mutation_rate < 0) {
            cout << "Error: --scaled_mutation_rate must be positive." << endl;
            exit(-1);
        }
    }

    if ((vm.count("scaled_recombination_rate") > 0) && (vm.count("recombination_to_mutation_ratio") > 0)) {
        cout << "Error: --scaled_recombination_rate and --recombination_to_mutation_ratio are mutually exclusive." << endl;
        exit(-1);    
    }

    if ((vm.count("scaled_recombination_rate") == 0) && (vm.count("recombination_to_mutation_ratio") == 0)) {
        cout << "Error: Either --scaled_recombination_rate or --recombination_to_mutation_ratio must be specified." << endl;
        exit(-1);    
    }

    float scaled_recombination_rate;
    if (vm.count("scaled_recombination_rate")) {
        scaled_recombination_rate = vm["scaled_recombination_rate"].as<float>();
        if (scaled_recombination_rate < 0) {
            cout << "Error: --scaled_recombination_rate must be positive." << endl;
            exit(-1);
        }
    }

    float recombination_to_mutation_ratio = -1;
    if (vm.count("recombination_to_mutation_ratio")) {
        recombination_to_mutation_ratio = vm["recombination_to_mutation_ratio"].as<float>();
        if (recombination_to_mutation_ratio < 0) {
            cout << "Error: --recombination_to_mutation_ratio must be positive." << endl;
            exit(-1);
        }
    }

    string flow_field_filename;
    if (vm.count("flow_field")) {
        flow_field_filename = vm["flow_field"].as<string>();
        if (!std::filesystem::exists(flow_field_filename)) {
            cout << boost::format("Error: Cannot open --flow_field file: %s\n") % flow_field_filename;
            exit(-1);
        }
    }

    if (vm.count("input") == 0) {
        cout << boost::format("Error: --input required.\n");
        exit(-1);
    }

    string input_filename = vm["input"].as<string>();
    if (!std::filesystem::exists(input_filename)) {
        cout << boost::format("Error: Cannot open --input file: %s\n") % input_filename;
        exit(-1);
    }

    if (vm.count("output") == 0) {
        cout << boost::format("Error: --output required.\n");
        exit(-1);
    }
    auto output_filename = vm["output"].as<string>();
    auto output_directory = std::filesystem::path(output_filename).parent_path();
    if (!output_directory.empty()) {
        std::filesystem::create_directories(output_directory);
    }
    
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
        if (!std::filesystem::exists(mask_filename)) {
            cout << boost::format("Error: Cannot open --mask file: %s\n") % mask_filename;
            exit(-1);
        }
    }

    string masks_per_sample_filename;
    if (vm.count("masks_per_filename")) {
        masks_per_sample_filename = vm["masks_per_filename"].as<string>();
        if (!std::filesystem::exists(masks_per_sample_filename)) {
            cout << boost::format("Error: Cannot open --masks_per_filename file: %s\n") % masks_per_sample_filename;
            exit(-1);
        }
    }

    string samples_filename;
    if (vm.count("samples")) {
        samples_filename = vm["samples"].as<string>();
        if (!std::filesystem::exists(samples_filename)) {
            cout << boost::format("Error: Cannot open --samples file: %s\n") % samples_filename;
            exit(-1);
        }
    }

    string samples_against_filename;
    if (vm.count("samples_against")) {
        samples_against_filename = vm["samples_against"].as<string>();
        if (!std::filesystem::exists(samples_against_filename)) {
            cout << boost::format("Error: Cannot open --samples_against file: %s\n") % samples_against_filename;
            exit(-1);
        }
    }

    int output_at_stride = vm["output_at_stride"].as<int>();
    bool output_at_hets = vm["output_at_hets"].as<bool>();
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

    screen.print_done();

    //
    // Estimated scaled mutation rate if needed
    //
    screen.print_subtitle("Calculating rates...");
    if (vm.count("scaled_mutation_rate") == 0) {
        screen.print_item("Estimating mutation rate...");
        scaled_mutation_rate = data_processor.calculate_heterozygosity();
    }

    if (vm.count("recombination_to_mutation_ratio")) {
        scaled_recombination_rate = scaled_mutation_rate * recombination_to_mutation_ratio;
    }

    screen.print_item(boost::str(boost::format("Scaled mutation rate: %f") % scaled_mutation_rate));
    screen.print_item(boost::str(boost::format("Scaled recombination rate: %f") % scaled_recombination_rate));
    screen.print_done();        
    
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

    if (vm.count("flow_field")) {
        read_flow_field_raw(
            flow_field_filename,
            mean_grid_def,
            cv_grid_def,
            flow_field_unravelled
        );
    } else {
        read_flow_field_default(
            mean_grid_def,
            cv_grid_def,
            flow_field_unravelled
        );
    }

    //
    // Prepare output files
    //
    // TODO: Check errors
    ofstream* output_file_raw_meta = NULL;
    ofstream* output_file_raw = NULL;
    if (output_filename.size() > 0) {
        output_file_raw = new ofstream(output_filename, ios_base::out | ios_base::binary);        
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
        output_file_raw_meta,
        output_file_raw,
        vm["zstd_compression_level"].as<int>()
    );

    PPC.calculate_posteriors();

    //
    // Close output files
    //
    if (output_file_raw != NULL) {
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
