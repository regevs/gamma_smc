#pragma once

#include "common.h"

class DataProcessor {
  public:

    const vector<unique_ptr<SegregatingSite>>& _sites;
    const vector<string>& _sample_names;
    const vector<pair<int, int>>& _global_mask;
    unordered_map<string, vector<pair<int, int>>>& _mask_map;
    int _posterior_every;
    int _flow_field_cache_n_steps;
    bool _output_at_every;
    bool _output_at_hets;

    bool _is_global_mask;
    vector<position_t> _output_positions;
    vector<Segment_t> _segments;
    vector<vector<bool>> _is_seg_site_missing;
    int _n_segments;
    position_t _seq_length;

    DataProcessor(
        const vector<unique_ptr<SegregatingSite>>& sites,
        const vector<string>& sample_names,
        const vector<pair<int, int>>& global_mask,
        unordered_map<string, vector<pair<int, int>>>& mask_map,
        int posterior_every,
        int flow_field_cache_n_steps,
        bool output_at_hets
    ) :
        _sites(sites),
        _sample_names(sample_names),
        _global_mask(global_mask),
        _mask_map(mask_map),
        _posterior_every(posterior_every),
        _flow_field_cache_n_steps(flow_field_cache_n_steps),
        _output_at_every(_posterior_every > 0),
        _output_at_hets(output_at_hets),
        _is_global_mask(_mask_map.size() == 0)
    {
        prepare_segments();        
        precalculate_seg_sites_in_masks();
    }

    bool is_global_mask() {
        return _is_global_mask;
    }

    void prepare_segments() {
        position_t p;
        bool prev_to_output = false;
        bool to_output;
        int n_pos = -1;
        _seq_length = 0;       

        position_t last = _sites.back()->pos;
        position_t jump_to_pos;
        position_t next_output = 0;
        int next_site_index = 0;
        segment_type seg_type;

        while (n_pos < last) {
            jump_to_pos = _sites[next_site_index]->pos;
            if (_output_at_every) {
                jump_to_pos = std::min(jump_to_pos, next_output);
            }

            // While we can't jump to the next site within the cache
            while (jump_to_pos - n_pos > _flow_field_cache_n_steps) {
                _segments.push_back({
                    n_pos + _flow_field_cache_n_steps,      // pos
                    _flow_field_cache_n_steps,              // segment length
                    HOM_STRETCH,                            // segment type - meaningless default here; TODO: Remove field?
                    prev_to_output,                         // output at start
                    false,                                  // output at end
                    -1,                                     // No corresponding site
                });
                
                n_pos += _flow_field_cache_n_steps;
                prev_to_output = false;
            }

            // output only if: (i) we output hets; or (ii) it's an output site
            to_output = (_output_at_every && (jump_to_pos == next_output)) || _output_at_hets;

            // Now we can jump
            _segments.push_back({
                jump_to_pos,                       // pos
                (int) jump_to_pos - n_pos,    // segment length
                HOM_STRETCH,                     // segment type - meaningless default here; TODO: Remove field?
                prev_to_output,     // output at start
                to_output,          // output at end
                ((jump_to_pos == _sites[next_site_index]->pos) ? next_site_index : -1)                                
            });

            if (to_output) { _output_positions.push_back(jump_to_pos); }

            // Update pointers
            n_pos = jump_to_pos;  
            prev_to_output = to_output;

            if (jump_to_pos == _sites[next_site_index]->pos) {
                next_site_index++;
            }
            if (_output_at_every && (jump_to_pos == next_output)) {
                next_output += _posterior_every;
            }  
        }

        _n_segments = _segments.size();
        _seq_length = _output_positions.size();
    }

    void precalculate_seg_sites_in_masks() {
        _is_seg_site_missing.resize(_sample_names.size(), std::vector<bool>(_sites.size(), true));

        for (uint n_sample = 0; n_sample < _sample_names.size(); n_sample++) {
            int count = 0;
            uint n_interval = 0;            
            const vector<pair<int, int>>* mask;
            if (is_global_mask()) {
                mask = &_global_mask;
            } else {
                mask = &(_mask_map.at(_sample_names[n_sample]));
            };

            for (uint n_seg_site = 0; n_seg_site < _sites.size(); n_seg_site++) {                
                while ((n_interval < mask->size()) && ((*mask)[n_interval].first <= _sites[n_seg_site]->pos)) {
                    if ((*mask)[n_interval].second > _sites[n_seg_site]->pos) {
                        _is_seg_site_missing[n_sample][n_seg_site] = false;
                        count++;
                        break;
                    }
                    n_interval++;
                }                
            }
            // cout << boost::format("# not missing in %s = %d\n") % _sample_names[n_sample] % count;
        }
    }

    float calculate_heterozygosity_for_sample(int n_sample) {
        const vector<pair<int, int>>* mask;
            if (is_global_mask()) {
                mask = &_global_mask;
            } else {
                mask = &(_mask_map.at(_sample_names[n_sample]));
            };

        // Count hets not missing and not in mask
        int n_hets = 0;
        for (uint n_seg_site = 0; n_seg_site < _sites.size(); n_seg_site++) {
            // If missing, skip
            if ((_sites[n_seg_site]->alleles[2 * n_sample] == -1) || (_sites[n_seg_site]->alleles[2 * n_sample + 1] == -1)) {
                continue;
            }

            // If hom, skip
            if (_sites[n_seg_site]->alleles[2 * n_sample] == _sites[n_seg_site]->alleles[2 * n_sample + 1]) {
                continue;
            }

            // If masked, skip
            if (_is_seg_site_missing[n_sample][n_seg_site]) {
                continue;
            }

            // Otherwise, count
            n_hets++;            
        }
        
        // Count number of sites in mask
        int n_total_sites = 0;
        for (auto& interval : *mask) {
            n_total_sites += (interval.second-interval.first);
        }

        return ((float) n_hets) / n_total_sites;
    }

    float calculate_heterozygosity() {
        float theta = 0.0;
        for (uint n_sample = 0; n_sample < _sample_names.size(); n_sample++) {
            theta += calculate_heterozygosity_for_sample(n_sample);
        }
        return theta / _sample_names.size();
    }

    void intersect_masks(
        const vector<pair<int, int>>& mask1,
        const vector<pair<int, int>>& mask2,
        int32_t* n_called_array,
        int ptr_jump = 1            // stride into n_called_array
        ) {
        int i = 0, j = 0, k = 0;
        int n = mask1.size(), m = mask2.size();

        int cur_segment_start = 0;
        int cur_segment_end = _segments[0].pos;

        int n_called = 0;
        int n_interval = 0;
    
        // Loop through all intervals unless one of the interval gets exhausted
        while ((i < n) && (j < m) && (k < _n_segments)) {
            // cout << boost::format("%d %d - %d %d - %d %d\n") 
            //     % mask1[i].first % mask1[i].second % mask2[j].first % mask2[j].second % cur_segment_start % cur_segment_end;

            // Left bound for intersecting interval
            int l = max(max(mask1[i].first, mask2[j].first), cur_segment_start);
    
            // Right bound for intersecting interval
            int r = min(min(mask1[i].second, mask2[j].second), cur_segment_end);
    
            // If interval is valid print it
            if (l < r) {
                //cout << "{" << l << ", " << r << "}\n";
                n_interval++;
                n_called += (r - l);
            }
    
            if (cur_segment_end < mask1[i].second) {
                if (cur_segment_end < mask2[j].second) {
                    *n_called_array = n_called; 
                    n_called_array += ptr_jump;                   
                    n_called = 0;                    

                    k++;
                    cur_segment_start = cur_segment_end;
                    cur_segment_end = _segments[k].pos;   
                } else {
                    j++;
                }
            } else {
                if (mask1[i].second < mask2[j].second) {
                    i++;                    
                } else {
                    j++;
                }
            }            
        }

        while (k < _n_segments) {
            *n_called_array = 0;
            n_called_array += ptr_jump;
            k++;
        }
        
    }

    void intersect_global_mask(
        int32_t* n_called_array,
        int ptr_jump = 1            // stride into n_called_array
        ) {
        intersect_masks(
            _global_mask, 
            _global_mask,
            n_called_array,
            ptr_jump
        );
    }

    void intersect_masks_at_ids(
        int sample_id_1, 
        int sample_id_2,
        int32_t* n_called_array,
        int ptr_jump = 1            // stride into n_called_array
        ) {
        auto& mask1 = _mask_map.at(_sample_names[sample_id_1]);
        auto& mask2 = _mask_map.at(_sample_names[sample_id_2]);

        intersect_masks(
            mask1, 
            mask2,
            n_called_array,
            ptr_jump
        );
    }

    
};