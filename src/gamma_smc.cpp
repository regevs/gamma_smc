#include "common.h"

#include "sys.h"
#include "io.h"
#include "flow_field.h"
#include "gamma_smc.h"
#include "data_processor.h"

int main(int argc, char** argv) {

    // Print command line
    cout << "Command line:" << endl;
    for (int i = 0; i < argc; ++i) {
        cout << argv[i] << ' ';
    }
    cout << endl << endl;

    //
    // Parse flags
    //
    po::options_description desc("Allowed options");
    desc.add_options()
        ("scaled_mutation_rate,m", po::value<float>()->default_value(0.000375), "Scaled mutation rate")
        ("scaled_recombination_rate,r", po::value<float>()->default_value(0.0003), "Scaled recombination rate")
        ("flow_field_file,f", po::value<string>(), "Flow field file")
        ("input_file,i", po::value<string>(), "Input file")
        ("beds_file,b", po::value<string>(), "File of bed filenames (empty for no masks)")
        ("output_text_file", po::value<string>()->default_value(""), "Output file in text format (recommended only for small datasets!)")
        ("output_file,o", po::value<string>()->default_value(""), "Output file")
        ("stride,s", po::value<int>()->default_value(-1), "Output posterior every")
        ("output_at_hets,h", po::value<bool>()->default_value(false), "Output at hets")
        ("use_cache,c", po::value<bool>()->default_value(true), "Use cache")
        // ("log_coords,l", po::value<bool>()->default_value(true), "Use log coordinates in output?")
        ("cache_size,z", po::value<int>()->default_value(1000), "Stride width in basepairs")
        ("only_forward,y", po::value<bool>()->default_value(false), "Calculate only forward")
        ("only_backward", po::value<bool>()->default_value(false), "Calculate only backward")
        //("output_n_called_file", po::value<string>()->default_value(""), "Filename of n of called")
        ("samples_file,S", po::value<string>()->default_value(""), "Filename of a list of subset of samples to take")
        //("output_samples_file", po::value<string>()->default_value(""), "Filename of a list of subset of samples used in file (in order)")
        ("entropy_clipping", po::value<bool>()->default_value(false), "Clip by entropy")
    ;

    po::positional_options_description p;
    p.add("input_file", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

    po::notify(vm);

    

    // TODO: Validate flags thoroughly

    // Load input file
    //vector<SegSite_t> input_sites;
    vector<unique_ptr<SegregatingSite>> input_sites;  // TODO: Have another go at getting rid of the unique_ptr and just have vector<Seg...>
    vector<string> sample_names;    

    // readSegSitesAll(vm["input_file"].as<string>(), input_sites);
    readVcf(vm["input_file"].as<string>(), input_sites, sample_names, vm["samples_file"].as<string>()); 

    unordered_map<string, vector<pair<int, int>>> mask_map;
    if (vm.count("beds_file")) {
        readMasks(vm["beds_file"].as<string>(), mask_map, sample_names);
        cout << boost::format("Read %d bed filenames from list...\n") % mask_map.size();    
    }
    

    DataProcessor data_processor(
        input_sites,
        sample_names,
        mask_map,
        vm["stride"].as<int>(),
        vm["cache_size"].as<int>(),
        vm["output_at_hets"].as<bool>()
    );

    cout << "Evaluating at " << data_processor._n_segments << " positions, outputting at " << data_processor._seq_length << endl;

    // vector<int> n_called_array(data_processor._n_segments);
    // data_processor.intersect_masks(0, 1, n_called_array.data());

    // exit(-1);


    // TODO: Print time benchmarks, print some stats

    // int bb = 0;
    // for (const auto& stuff: input_sites) {
    //     cout << stuff.pos << " - " << stuff.type << endl;
    //     bb ++;
    //     if (bb > 10) { break; }
    // }
    cout << boost::format("Loaded file, with %d sites\n") % input_sites.size() << endl;

    // Create a list of pairs to work on
    vector<pair<int, int>> haplotype_pairs;
    int n_haplotypes = sample_names.size() * 2;
    for (int i = 0; i < n_haplotypes; i++) {
        for (int j = i+1; j < n_haplotypes; j++) {
            haplotype_pairs.push_back(make_pair(i, j));
        }
    }
    cout << boost::format("Applying to %d haplotype pairs\n") % haplotype_pairs.size() << endl;


    // Load flow field file
    vector<float> mean_grid_def;
    vector<float> cv_grid_def;
    vector<float> flow_field_unravelled;

    read_flow_field_raw(
        vm["flow_field_file"].as<string>(),
        mean_grid_def,
        cv_grid_def,
        flow_field_unravelled
    );

    // Prepare output files
    // TODO: Create dirs, check errors
    ostream* out = NULL;
    ofstream* output_file = NULL;
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf;
    if (vm["output_file"].as<string>().size() > 0) {        
        output_file = new ofstream(vm["output_file"].as<string>(), ios_base::out | ios_base::binary);        
        outbuf.push(boost::iostreams::gzip_compressor());
        outbuf.push(*output_file);
        out = new ostream(&outbuf);
    }

    ostream* out_raw = NULL;
    ofstream* output_file_raw_meta = NULL;
    ofstream* output_file_raw = NULL;
    boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf_raw;
    if (vm["output_raw_file"].as<string>().size() > 0) {
        // std::ofstream out;
        // out.open(vm["output_raw_file"].as<string>(), std::ios::out | std::ios::binary);

        output_file_raw = new ofstream(vm["output_raw_file"].as<string>(), ios_base::out | ios_base::binary);        
        outbuf_raw.push(boost::iostreams::gzip_compressor());
        outbuf_raw.push(*output_file_raw);
        out_raw = new ostream(&outbuf_raw);

        output_file_raw_meta = new ofstream(vm["output_raw_file"].as<string>() + ".meta", ios_base::out);        
    }

    bool use_cache = vm["use_cache"].as<bool>();
    if (use_cache) {
        // Construct flow field
        unique_ptr<FlowField> FF(new FlowField(mean_grid_def, cv_grid_def, flow_field_unravelled));

        // Construct flow field cache
        cout << "Building flow field cache..." << endl;

        auto t1 = std::chrono::high_resolution_clock::now();
        
        unique_ptr<FlowFieldCache> FFC(new FlowFieldCache(
            vm["scaled_recombination_rate"].as<float>(),
            vm["scaled_mutation_rate"].as<float>(),
            move(FF),
            vm["cache_size"].as<int>(),
            vm["entropy_clipping"].as<bool>()
        ));

        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<float, std::milli> ms_float = t2 - t1;

        cout << boost::format("Done (%f sec)") % (ms_float.count()/1000) << endl;

        // ofstream missing_cache;
        // missing_cache.open("/home/rs2145/temp/cache.missing.bin", std::ofstream::binary);
        // missing_cache.write((char*)&(FFC->_cached_missing_unravelled[0]), FFC->_cached_missing_unravelled.size() * sizeof(float));
        // missing_cache.close();

        // ofstream hom_cache;
        // hom_cache.open("/home/rs2145/temp/cache.hom.bin", std::ofstream::binary);
        // hom_cache.write((char*)&(FFC->_cached_hom_unravelled[0]), FFC->_cached_hom_unravelled.size() * sizeof(float));
        // hom_cache.close();

        CachedPairwiseGammaSMC PPC(
            input_sites,
            haplotype_pairs,
            vm["scaled_recombination_rate"].as<float>(),
            vm["scaled_mutation_rate"].as<float>(),
            move(FFC),
            data_processor,
            vm["stride"].as<int>(),
            vm["output_at_hets"].as<bool>(),
            vm["only_forward"].as<bool>(),
            vm["only_backward"].as<bool>(),
            out,
            output_file_raw_meta,
            out_raw
        );

        // Run
        cout << "Running..." << endl;


        t1 = std::chrono::high_resolution_clock::now();  

        // __itt_frame_begin_v3(pD, NULL);     
        PPC.calculate_posteriors();
        // __itt_frame_end_v3(pD, NULL);

        t2 = std::chrono::high_resolution_clock::now();
        ms_float = t2 - t1;
        cout << boost::format("calculate_posteriors() TOTAL (%f sec)") % (ms_float.count()/1000.0) << endl;

     
        /*if (vm["output_n_called_file"].as<string>().size() > 0) {
            std::ofstream out;
            out.open(vm["output_n_called_file"].as<string>(), std::ios::out);
            int32_t* cur_n_called_ptr = PPC._pairwise_n_called;
            for (int n_segment = 0; n_segment < PPC._n_segments; n_segment++) {
                out << PPC._segments[n_segment].pos;
                for (int i = 0; i < PPC._n_pairs_rounded_up; i++) {
                    out << "\t" << *cur_n_called_ptr;
                    cur_n_called_ptr++;
                }
                out << "\n";
            }
            out.close();
        }   
        */     
        
    } 
    /*else {
        // Construct flow field
        unique_ptr<FlowField> FF(new FlowField(mean_grid_def, cv_grid_def, flow_field_unravelled));
        
        // Construct main object
        PairwiseGammaSMC PP(
            input_sites,
            vm["scaled_recombination_rate"].as<float>(),
            vm["scaled_mutation_rate"].as<float>(),
            move(FF),
            vm["log_coords"].as<bool>()
        );   

        

        // Run
        cout << "Running..." << endl;

        auto t1 = std::chrono::high_resolution_clock::now();
        PP.calculate_posteriors();
        
        auto t2 = std::chrono::high_resolution_clock::now();
        auto ms_float = t2 - t1;

        cout << boost::format("Done (%f sec)") % (ms_float.count()/1000.0) << endl;

        // Print
        if (vm["output_file"].as<string>().size() > 0) {
            ofstream output_file;
            output_file.open(vm["output_file"].as<string>());
            for (uint i = 0; i < PP._posteriors_alpha.size(); i++) {
                if (i % vm["stride"].as<int>() == 0) {
                    output_file << boost::format("%d\t%32.20f\t%32.20f\t%32.20f\t%32.20f\t%32.20f\t%32.20f\n")
                        % i  
                        % PP._scaled_forwards_alpha[i] % PP._scaled_forwards_beta[i] 
                        % PP._rev_scaled_forwards_alpha[i] % PP._rev_scaled_forwards_beta[i] 
                        % PP._posteriors_alpha[i] % PP._posteriors_beta[i]; 
                    // output_file << i << "\t" << PP._scaled_forwards_alpha[i] << "\t" << PP._scaled_forwards_beta[i] 
                    // << "\t" << PP._rev_scaled_forwards_alpha[i] << "\t" << PP._rev_scaled_forwards_beta[i] 
                    // << "\t" << PP._posteriors_alpha[i] << "\t" << PP._posteriors_beta[i] << endl; 
                }
            }    

            output_file.close();
        }

    }
*/

    // Close output files
    if (out != NULL) {
        boost::iostreams::close(outbuf); // Don't forget this!
        output_file->close();
    }

    if (out_raw != NULL) {
        boost::iostreams::close(outbuf_raw); // Don't forget this!
        output_file_raw->close();
        output_file_raw_meta->close();
    }

    fprintf(stderr, "\nCPU: %.3f sec; Peak RSS: %.3f GB\n",mp_cputime(), mp_peakrss() / 1024.0 / 1024.0 / 1024.0);	

    return 0;
}
