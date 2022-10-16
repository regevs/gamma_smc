#pragma once

#include "common.h"
#include "flow_field.h"
#include "data_processor.h"

// https://gist.github.com/andersx/8057b2a6fd3d715d35eb

// Approximation for EXP(x) -- very fast, but not super accurate
static inline __m256 _mm256_expfaster_ps(const __m256 &q) {

     // REGEV: Added this to make it 10**x instead of exp(x)
    const __m256 ln10 =  _mm256_set1_ps (2.30258509299f);

    const __m256 C1 = _mm256_set1_ps(1064872507.1541044f);
    const __m256 C2 = _mm256_set1_ps(12102203.161561485f);

    return _mm256_castsi256_ps(_mm256_cvttps_epi32(_mm256_fmadd_ps(C2, _mm256_mul_ps(q, ln10), C1)));
}

// --------------------------------------------------------------------------------------------------

class CachedPairwiseGammaSMC {
  public:
    const vector<unique_ptr<SegregatingSite>>& _sites;
    float _scaled_recombination_rate;
    float _scaled_mutation_rate;
    unique_ptr<FlowFieldCache> _flow_field_cache;
    DataProcessor& _data_processor;
    int _posterior_every;    
    int _flow_field_cache_n_steps;

    int _n_haplotypes;
    long _n_pairs;
    long _n_pairs_rounded_up;

    bool _output_at_hets;
    bool _only_forward;
    
    const vector<position_t>& _output_positions;
    const vector<Segment_t>& _segments;
    
    int _n_segments;
    position_t _seq_length;

    // TODO: allocate less memory
    float* _scaled_forwards_mean = NULL;
    float* _scaled_forwards_cv = NULL;
    float* _scaled_backwards_mean = NULL;
    float* _scaled_backwards_cv = NULL;
    float* _posteriors_alpha = NULL;
    float* _posteriors_beta = NULL;

    int8_t* _pairwise_segment_types = NULL;

    float* _last_processed_mean_log10 = NULL;
    float* _last_processed_cv_log10 = NULL;

    int32_t* _pairwise_n_called = NULL;

    CachedPairwiseGammaSMC(
        const vector<unique_ptr<SegregatingSite>>& sites,
        float scaled_recombination_rate,
        float scaled_mutation_rate,
        unique_ptr<FlowFieldCache> flow_field_cache,  // TODO: Do we really need unique_ptr rather than const&
        DataProcessor& data_processor,
        int posterior_every,
        bool output_at_hets,
        bool only_forward
    ) : 
        _sites(sites), 
        _scaled_recombination_rate(scaled_recombination_rate),
        _scaled_mutation_rate(scaled_mutation_rate), 
        _flow_field_cache(move(flow_field_cache)),
        _data_processor(data_processor),
        _posterior_every(posterior_every),
        _flow_field_cache_n_steps(_flow_field_cache->_n_steps),
        _n_haplotypes(_sites.back()->alleles.size()),
        _n_pairs(((long) _n_haplotypes) * (_n_haplotypes-1) / 2),
        _n_pairs_rounded_up(parallel_vector_size * (long) ceil(_n_pairs / float(parallel_vector_size))),
        _output_at_hets(output_at_hets),
        _only_forward(only_forward),
        _output_positions(data_processor._output_positions),
        _segments(data_processor._segments),
        _n_segments(data_processor._n_segments),
        _seq_length(data_processor._seq_length)
    {        
        // Allocate memory
        
        // Initialize with the stationary distribution, which is Exp(1) = Gamma(alpha=1, beta=1) 
        // = Gamma(mean=1, CV=1) - so 0 in log space is fine

        size_t n_elements = (
            _seq_length *           // # of positions to output
            _n_pairs_rounded_up     // # of pairs, rounded up to AVX vector size
            );

        // _scaled_forwards_mean = aligned_alloc_float(n_elements, false);
        // _scaled_forwards_cv = aligned_alloc_float(n_elements, false);
        // _scaled_backwards_mean = aligned_alloc_float(n_elements, false);
        // _scaled_backwards_cv = aligned_alloc_float(n_elements, false);
        _posteriors_alpha = aligned_alloc_float(n_elements, true);
        _posteriors_beta = aligned_alloc_float(n_elements, true);

        // Working vector
        _last_processed_mean_log10 = aligned_alloc_float(_n_pairs_rounded_up, false);
        _last_processed_cv_log10 = aligned_alloc_float(_n_pairs_rounded_up, false);        
    }

    virtual ~CachedPairwiseGammaSMC() {
        free(_scaled_forwards_mean);
        free(_scaled_forwards_cv);
        free(_scaled_backwards_mean);
        free(_scaled_backwards_cv);
        free(_posteriors_alpha);
        free(_posteriors_beta);
        free(_pairwise_segment_types);
        free(_last_processed_mean_log10);
        free(_last_processed_cv_log10);
        free(_pairwise_n_called);
    }

    void prepare_pairwise_emissions() {
        // Allocate _n_segments rows of _n_pairs_rounded_up each
        _pairwise_segment_types = aligned_alloc_int8(_n_segments * _n_pairs_rounded_up, false);
        _pairwise_n_called = aligned_alloc_int32(_n_segments * _n_pairs_rounded_up, false);

        // Start by figuring out the number of called positions within each segment, for each pair
        int32_t* cur_n_called_ptr = _pairwise_n_called;
        int n_pair = 0;
        for (int i = 0; i < _n_haplotypes-1; i++) {
            for (int j = i+1; j < _n_haplotypes; j++) {
                _data_processor.intersect_masks(
                    i >> 1,         // First sample ID
                    j >> 1,         // Second sample ID
                    cur_n_called_ptr + n_pair,      // Pointer to place in big array
                    _n_pairs_rounded_up             // Stride into array
                );

                n_pair++;
            }
        }
        
        // Now, for each segment, figure out the emission type. If this segments
        // ends in a segregating site, look at the alleles, and also look if it falls
        // within the mask(s). 
        // Otherwise, the emission is considered to be HOM, as every segment is 
        // actually a stretch of missing and a stretch of hom.
        int8_t* cur_ptr;        
        for (long n_segment = 0; n_segment < _n_segments; n_segment++) {
            cur_ptr = _pairwise_segment_types + _n_pairs_rounded_up * n_segment;

            int seg_site_index = _segments[n_segment].seg_site_index;

            if (seg_site_index == -1) {
                // We're in the middle of a hom stretch
                // TODO: Maybe - this should not be HOM_STRETCH_HOM_SITE as this introduces a hom site at the end
                // no matter what. Perhaps we need to extend _is_seg_site_missing to all segment ends.
                memset((void*) cur_ptr, HOM_STRETCH_HOM_SITE, _n_pairs_rounded_up);                  
            } else {        
                auto& alleles = _sites[seg_site_index]->alleles;

                for (int i = 0; i < _n_haplotypes-1; i++) {
                    bool is_missing_i = _data_processor._is_seg_site_missing[i >> 1][seg_site_index];
                    for (int j = i+1; j < _n_haplotypes; j++) {
                        bool is_missing_j = _data_processor._is_seg_site_missing[j >> 1][seg_site_index];
                        if (is_missing_i || is_missing_j) {
                            (*cur_ptr) = HOM_STRETCH;
                        } else {                            
                            (*cur_ptr) = ((alleles[i] ^ alleles[j]) ? HOM_STRETCH_HET_SITE : HOM_STRETCH_HOM_SITE);
                        }
                        cur_ptr++;
                    }
                }

                memset((void*) cur_ptr, MISSING_STRETCH_MISSING_SITE, _n_pairs_rounded_up - _n_pairs); // Doesn't matter, just a default
            }
        }
    }

    void mean_cv_to_alpha_beta_log10_vec_forward(
        float* mean_log10, 
        float* cv_log10,
        float* alpha_log10_output,
        float* beta_log10_output
        ) {

        float alpha_log10, beta_log10;

        // for (int n_pair = 0; n_pair < _n_pairs_rounded_up; n_pair++) {
        //     alpha_log10 = -2 * cv_log10[n_pair];
        //     beta_log10  = alpha_log10 - mean_log10[n_pair];

        //     alpha_log10_output[n_pair] = powFastLookup(alpha_log10, _pTable);
        //     beta_log10_output[n_pair]  = powFastLookup(beta_log10, _pTable);
        // }

        const __m256 minus_two = _mm256_set1_ps(-2);
        __m256 alpha_log10_vec, beta_log10_vec;
        for (int n_pair = 0; n_pair < _n_pairs_rounded_up; n_pair += parallel_vector_size) {
            alpha_log10_vec = _mm256_mul_ps(minus_two, _mm256_load_ps(cv_log10 + n_pair));
            beta_log10_vec = _mm256_sub_ps(alpha_log10_vec,  _mm256_load_ps(mean_log10 + n_pair));            

            _mm256_store_ps(
                alpha_log10_output + n_pair,
                _mm256_expfaster_ps(alpha_log10_vec)
            );

            _mm256_store_ps(
                beta_log10_output + n_pair,
                _mm256_expfaster_ps(beta_log10_vec)
            );
        }
    }

    void mean_cv_to_alpha_beta_log10_vec_backward(
        float* mean_log10, 
        float* cv_log10,
        float* alpha_log10_output,
        float* beta_log10_output
        ) {

        float alpha_log10, beta_log10;

        // for (int n_pair = 0; n_pair < _n_pairs_rounded_up; n_pair++) {
        //     alpha_log10 = -2 * cv_log10[n_pair];
        //     beta_log10  = alpha_log10 - mean_log10[n_pair];

        //     alpha_log10_output[n_pair] += powFastLookup(alpha_log10, _pTable) - 1;
        //     beta_log10_output[n_pair]  += powFastLookup(beta_log10, _pTable) - 1;
        // }
        const __m256 minus_two = _mm256_set1_ps(-2);
        const __m256 minus_one = _mm256_set1_ps(-1);

        __m256 alpha_log10_vec, beta_log10_vec;
        for (int n_pair = 0; n_pair < _n_pairs_rounded_up; n_pair += parallel_vector_size) {
            alpha_log10_vec = _mm256_mul_ps(minus_two, _mm256_load_ps(cv_log10 + n_pair));
            beta_log10_vec = _mm256_sub_ps(alpha_log10_vec,  _mm256_load_ps(mean_log10 + n_pair));            

            _mm256_store_ps(
                alpha_log10_output + n_pair,
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_load_ps(alpha_log10_output + n_pair),
                        _mm256_expfaster_ps(alpha_log10_vec)
                    ),
                    minus_one)
            );

            _mm256_store_ps(
                beta_log10_output + n_pair,
                _mm256_add_ps(
                    _mm256_add_ps(
                        _mm256_load_ps(beta_log10_output + n_pair),
                        _mm256_expfaster_ps(beta_log10_vec)
                    ),
                    minus_one)
            );
        }
    }    

    void forward_vectorized() {
        memset((void *) _last_processed_mean_log10, 0, _n_pairs_rounded_up * sizeof(float));
        memset((void *) _last_processed_cv_log10, 0, _n_pairs_rounded_up * sizeof(float));

        Segment_t next_segment;

        // TODO: remove underscore, those are not member variables
        segment_type* _pairwise_segment_types_ptr = (segment_type*) _pairwise_segment_types;
        float* _scaled_forwards_mean_ptr = _scaled_forwards_mean;
        float* _scaled_forwards_cv_ptr = _scaled_forwards_cv;

        float* posteriors_alpha_ptr = _posteriors_alpha;
        float* posteriors_beta_ptr = _posteriors_beta;        
        
        float* _last_processed_mean_log10_ptr;
        float* _last_processed_cv_log10_ptr;

        int32_t* n_called_ptr = _pairwise_n_called;

        for (int next_processed_segment_index = 0; next_processed_segment_index < _n_segments; next_processed_segment_index++) {
            next_segment = _segments[next_processed_segment_index];

            _last_processed_mean_log10_ptr = _last_processed_mean_log10;
            _last_processed_cv_log10_ptr = _last_processed_cv_log10;
            
            for (int n_pair = 0; n_pair < _n_pairs_rounded_up; n_pair += parallel_vector_size) {
                _flow_field_cache->at_flat_vectorized(
                    _last_processed_mean_log10_ptr, 
                    _last_processed_cv_log10_ptr, 
                    (int32_t) next_segment.length,
                    n_called_ptr,                    
                    _pairwise_segment_types_ptr,  
                    true           // forward
                );                   

                _last_processed_mean_log10_ptr += parallel_vector_size;
                _last_processed_cv_log10_ptr += parallel_vector_size;
                _pairwise_segment_types_ptr += parallel_vector_size;
                n_called_ptr += parallel_vector_size;
            }

            // If we need to output, do so
            if (next_segment.output_at_end) {                        
                // memcpy(
                //     _scaled_forwards_mean_ptr,                         
                //     _last_processed_mean_log10,
                //     _n_pairs_rounded_up * sizeof(float)
                // );                    
                
                // memcpy(     
                //     _scaled_forwards_cv_ptr,                         
                //     _last_processed_cv_log10,
                //     _n_pairs_rounded_up * sizeof(float)
                // );                    

                // _scaled_forwards_mean_ptr += _n_pairs_rounded_up;
                // _scaled_forwards_cv_ptr += _n_pairs_rounded_up;

                mean_cv_to_alpha_beta_log10_vec_forward(
                    _last_processed_mean_log10,
                    _last_processed_cv_log10,
                    posteriors_alpha_ptr,
                    posteriors_beta_ptr
                );

                posteriors_alpha_ptr += _n_pairs_rounded_up;
                posteriors_beta_ptr += _n_pairs_rounded_up;
            }

            
        }        
    }


    void backward_vectorized() {  
        memset((void *) _last_processed_mean_log10, 0, _n_pairs_rounded_up * sizeof(float));
        memset((void *) _last_processed_cv_log10, 0, _n_pairs_rounded_up * sizeof(float));

        segment_type* _pairwise_segment_types_ptr; 
            
        // float* _scaled_backwards_mean_ptr = 
        //     _scaled_backwards_mean + _n_pairs_rounded_up * (_seq_length-1);
        // float* _scaled_backwards_cv_ptr =
        //      _scaled_backwards_cv + _n_pairs_rounded_up * (_seq_length-1);
        float* _last_processed_mean_log10_ptr;
        float* _last_processed_cv_log10_ptr;

        float* posteriors_alpha_ptr = 
            _posteriors_alpha + _n_pairs_rounded_up * (_seq_length-1);
        float* posteriors_beta_ptr = 
            _posteriors_beta + _n_pairs_rounded_up * (_seq_length-1);

        int32_t* n_called_ptr;

        Segment_t next_segment;

        for (int next_processed_segment_index = _n_segments-1; next_processed_segment_index >= 0; next_processed_segment_index--) {
            next_segment = _segments[next_processed_segment_index];            

            _last_processed_mean_log10_ptr = _last_processed_mean_log10;
            _last_processed_cv_log10_ptr = _last_processed_cv_log10;

            _pairwise_segment_types_ptr = (segment_type*) 
                (_pairwise_segment_types + _n_pairs_rounded_up * next_processed_segment_index);
            n_called_ptr = (int32_t*) (_pairwise_n_called + _n_pairs_rounded_up * next_processed_segment_index);

            for (int n_pair = 0; n_pair < _n_pairs_rounded_up; n_pair += parallel_vector_size) {
                // Jump to start of segment
                _flow_field_cache->at_flat_vectorized(
                    _last_processed_mean_log10_ptr, 
                    _last_processed_cv_log10_ptr, 
                    (int32_t) next_segment.length,
                    n_called_ptr,
                    _pairwise_segment_types_ptr,  
                    false           // backward
                );

                _last_processed_mean_log10_ptr += parallel_vector_size;
                _last_processed_cv_log10_ptr += parallel_vector_size;
                _pairwise_segment_types_ptr += parallel_vector_size;
                n_called_ptr += parallel_vector_size;
            }

            if (next_segment.output_at_start) {    
                // memcpy(
                //     _scaled_backwards_mean_ptr,                         
                //     _last_processed_mean_log10,
                //     _n_pairs_rounded_up * sizeof(float)
                // );                    
                
                // memcpy(     
                //     _scaled_backwards_cv_ptr,                         
                //     _last_processed_cv_log10,
                //     _n_pairs_rounded_up * sizeof(float)
                // );                    

                // _scaled_backwards_mean_ptr -= _n_pairs_rounded_up;
                // _scaled_backwards_cv_ptr -= _n_pairs_rounded_up;

                mean_cv_to_alpha_beta_log10_vec_backward(
                    _last_processed_mean_log10, 
                    _last_processed_cv_log10,
                    posteriors_alpha_ptr,
                    posteriors_beta_ptr
                );

                posteriors_alpha_ptr -= _n_pairs_rounded_up;
                posteriors_beta_ptr -= _n_pairs_rounded_up;
            }
        }
    }

    void calculate_posteriors() {
        auto t1 = std::chrono::high_resolution_clock::now();
        prepare_pairwise_emissions();
        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<float, std::milli> ms_float = t2 - t1;
        cout << boost::format("prepare_pairwise_emissions() - Done (%f sec)") % (ms_float.count()/1000) << endl;

        // t1 = std::chrono::high_resolution_clock::now();
        // forward();        
        // t2 = std::chrono::high_resolution_clock::now();
        // ms_float = t2 - t1;
        // cout << boost::format("forward() - Done (%f sec)") % (ms_float.count()/1000) << endl;

        t1 = std::chrono::high_resolution_clock::now();
        forward_vectorized();        
        t2 = std::chrono::high_resolution_clock::now();
        ms_float = t2 - t1;
        cout << boost::format("forward_vectorized() - Done (%f sec)") % (ms_float.count()/1000) << endl;

        // t1 = std::chrono::high_resolution_clock::now();
        // forward_vectorized_no_prep();        
        // t2 = std::chrono::high_resolution_clock::now();
        // ms_float = t2 - t1;
        // cout << boost::format("forward_vectorized_no_prep() - Done (%f sec)") % (ms_float.count()/1000) << endl;


        // t1 = std::chrono::high_resolution_clock::now();
        // backward();
        // t2 = std::chrono::high_resolution_clock::now();
        // ms_float = t2 - t1;
        // cout << boost::format("backward() - Done (%f sec)") % (ms_float.count()/1000) << endl;

        if (!_only_forward) {
            t1 = std::chrono::high_resolution_clock::now();
            backward_vectorized();
            t2 = std::chrono::high_resolution_clock::now();
            ms_float = t2 - t1;
            cout << boost::format("backward_vectorized() - Done (%f sec)") % (ms_float.count()/1000) << endl;
        }

        // t1 = std::chrono::high_resolution_clock::now();
        // join_forward_backward();    
        // t2 = std::chrono::high_resolution_clock::now();
        // ms_float = t2 - t1;
        // cout << boost::format("join_forward_backward() - Done (%f sec)") % (ms_float.count()/1000) << endl;
    }

    void output(ostream& output_file) {
        // output_file << "position";
        // for (int i = 0; i < _n_haplotypes-1; i++) {
        //     for (int j = i+1; j < _n_haplotypes; j++) {
        //         output_file << 
        //             boost::format("\tforward_mean_log10_%d_%d\tforward_cv_log10_%d_%d\tbackward_mean_log10_%d_%d\tbackward_cv_log10_%d_%d\tposterior_alpha_%d_%d\tposterior_beta_%d_%d")
        //                 % i % j % i % j % i % j % i % j % i % j % i % j;
        //     }
        // }
        // output_file << "\n";
        
        // for (uint i = 0; i < _seq_length; i++) {
        //     output_file << _output_positions[i]; 
        //     for (int n_pair = 0; n_pair < _n_pairs; n_pair++) {
        //         int idx = i * _n_pairs_rounded_up + n_pair;
        //         output_file 
        //             << "\t" << _scaled_forwards_mean[idx] << "\t" << _scaled_forwards_cv[idx] 
        //             << "\t" << _scaled_backwards_mean[idx] << "\t" << _scaled_backwards_cv[idx] 
        //             << "\t" << _posteriors_alpha[idx] << "\t" << _posteriors_beta[idx]; 
        //     }  
        //     output_file << "\n";
        // }
        output_file << "position";
        for (int i = 0; i < _n_haplotypes-1; i++) {
            for (int j = i+1; j < _n_haplotypes; j++) {
                output_file << 
                    boost::format("\tposterior_alpha_%d_%d\tposterior_beta_%d_%d")
                       % i % j % i % j;
            }
        }
        output_file << "\n";
        
        for (uint i = 0; i < _seq_length; i++) {
            output_file << _output_positions[i]; 
            for (int n_pair = 0; n_pair < _n_pairs; n_pair++) {
                uint64_t idx = i * _n_pairs_rounded_up + n_pair;
                output_file                    
                    << "\t" << _posteriors_alpha[idx] << "\t" << _posteriors_beta[idx]; 
            }  
            output_file << "\n";
        }
    }

    void output_raw(ostream& output_file) {           
        // Write the sequence length
        output_file.write(reinterpret_cast<const char*>(&_seq_length), sizeof(long));     

        // Write the number of pairs
        output_file.write(reinterpret_cast<const char*>(&_n_pairs), sizeof(long));

        // Write the positions
        for (uint i = 0; i < _seq_length; i++) {
            output_file.write(reinterpret_cast<const char*>(&_output_positions[i]), sizeof(long));
        }

        // Write alphas and betas
        for (uint i = 0; i < _seq_length; i++) {
            output_file.write(
                reinterpret_cast<const char*>(_posteriors_alpha + i * _n_pairs_rounded_up), 
                sizeof(float) * _n_pairs
            );
            output_file.write(
                reinterpret_cast<const char*>(_posteriors_beta + i * _n_pairs_rounded_up), 
                sizeof(float) * _n_pairs
            );
        }
        
    }

};