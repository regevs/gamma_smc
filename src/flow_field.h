#pragma once

#include "common.h"

class FlowField {
  public:  
    float _mean_min_log10;
    float _mean_max_log10;
    int _mean_n_steps;
    float _mean_step;

    float _cv_min_log10;
    float _cv_max_log10;
    int _cv_n_steps;
    float _cv_step;

    // TODO: Don't copy, use smart pointers
    vector<float> _flow_field_unravelled;

    FlowField(
            const vector<float>& mean_grid_def,
            const vector<float>& cv_grid_def,
            const vector<float>& flow_field_unravelled
        ) :
        _flow_field_unravelled(flow_field_unravelled)
        {
            // TODO: Validate arguments, test

            _mean_min_log10 = log10(mean_grid_def[0]);
            _mean_max_log10 = log10(mean_grid_def[1]);
            _mean_n_steps = mean_grid_def[2];
            cout << boost::format("means: %f %f %d\n") % _mean_min_log10 % _mean_max_log10 % _mean_n_steps;
            _mean_step = (_mean_max_log10 - _mean_min_log10) / (_mean_n_steps - 1);

            _cv_min_log10 = log10(cv_grid_def[0]);
            _cv_max_log10 = log10(cv_grid_def[1]);
            _cv_n_steps = cv_grid_def[2];
            cout << boost::format("cvs: %f %f %d\n") % _cv_min_log10 % _cv_max_log10 % _cv_n_steps;
            _cv_step = (_cv_max_log10 - _cv_min_log10) / (_cv_n_steps - 1);

        }

        // Access a particular element in the grids
        void _entry(int row, int col, float& u, float& v) {
            if ((row < 0) || (row >= _mean_n_steps) || (col < 0) || (col >= _cv_n_steps)) {
                throw out_of_range((boost::format("Flow field entry (%d, %d) out of range.") % row % col).str());
            }

            u = _flow_field_unravelled[_cv_n_steps * row + col];
            v = _flow_field_unravelled[_cv_n_steps * _mean_n_steps + _cv_n_steps * row + col];
        }

        // Interpolate, grid entry given in mean+cv log scale.
        void at(float mean_log10, float cv_log10, float& u, float& v) {
            if (std::isnan(mean_log10) || std::isnan(cv_log10)) {
                throw out_of_range((boost::format("Flow field at (%f, %f) - nan.") % mean_log10 % cv_log10).str());
            }

            float orig_mean_log10 = mean_log10;
            // TODO: Track how much clipping was needed
            // Clip if needed.
            mean_log10 = min(max(mean_log10, _mean_min_log10), _mean_max_log10);
            cv_log10 = min(max(cv_log10, _cv_min_log10), _cv_max_log10);
            
            // Get the four points of interpolation
            float m = (mean_log10 - _mean_min_log10) / _mean_step;
            float c = (cv_log10 - _cv_min_log10) / _cv_step;

            int m0 = int(floor(m));
            if (m0 == _mean_n_steps-1) {
                m0--;
            }
            int m1 = m0 + 1;
            int c0 = int(floor(c));
            if (c0 == _cv_n_steps-1) {
                c0--;
            }
            int c1 = c0 + 1;

            //cout << (boost::format("flow_field->at: %f %f %f %f %f %f %f %f\n") % mean_log10 % cv_log10 % m % m0 % m1 % c % c0 % c1);

            // Fill the four grid points
            float U00, U01, U10, U11;
            float V00, V01, V10, V11;

            try {
                _entry(m0, c0, U00, V00);
                _entry(m0, c1, U01, V01);
                _entry(m1, c0, U10, V10);
                _entry(m1, c1, U11, V11);
            } catch (const std::exception& e) {
                cout << _mean_max_log10 << endl;
                cout << orig_mean_log10 << " -> " << max(orig_mean_log10, _mean_min_log10) << " -> " 
                    << min(max(orig_mean_log10, _mean_min_log10), _mean_max_log10) << endl;
                cout << (boost::format("In: %f %f\n") % mean_log10 % cv_log10);
                cout << (boost::format("%f %f %f %f %f %f\n") % m % m0 % m1 % c % c0 % c1);
                cout << (boost::format("delta %e \n") % (orig_mean_log10-mean_log10));
                cout << (boost::format("Us: %f %f %f %f\n") % U00 % U01 % U10 % U11);
                cout << (boost::format("Vs: %f %f %f %f\n") % V00 % V01 % V10 % V11);
                throw e;
            }

            // cout << (boost::format("Us: %f %f %f %f\n") % U00 % U01 % U10 % U11);
            // cout << (boost::format("Vs: %f %f %f %f\n") % V00 % V01 % V10 % V11);
            

            float w00 = (m1-m) * (c1-c);
            float w01 = (m1-m) * (c-c0);
            float w10 = (m-m0) * (c1-c);
            float w11 = (m-m0) * (c-c0);

            u = w00 * U00 + w01 * U01 + w10 * U10 + w11 * U11;
            v = w00 * V00 + w01 * V01 + w10 * V10 + w11 * V11;
        }

        void at_alpha_beta(float alpha, float beta, float& u, float& v) {
            if ((alpha < 0) || (beta < 0)) {
                throw out_of_range((boost::format("Flow field at_alpha_beta (%f, %f) - negative entry.") % alpha % beta).str());
            }

            float alpha_log10 = log10(alpha); // icsi_log10(alpha, _LOOKUP_TABLE, _log_approx_N);
            float beta_log10 = log10(beta); //icsi_log10(beta, _LOOKUP_TABLE, _log_approx_N);

            at(alpha_log10 - beta_log10, -0.5 * alpha_log10, u, v);
        }

        void at_log10_alpha_beta(float alpha_log10, float beta_log10, float& u, float& v) {
            at(alpha_log10 - beta_log10, -0.5 * alpha_log10, u, v);
        }

    
};


//
// n_steps here are n_steps of flow field step for 'missing'
// and for 'hom' it's intertwined with hom observation update
//

inline int ftoi_sse1(float f)
{
    return _mm_cvtt_ss2si(_mm_load_ss(&f));     // SSE1 instructions for float->int
}

class FlowFieldCache {
  public:

    float _scaled_recombination_rate;
    float _scaled_mutation_rate;
    unique_ptr<FlowField> _flow_field;
    int _n_steps; 

    vector<float> _cached_missing_unravelled;
    vector<float> _cached_hom_unravelled;
    vector<float> _cached_hom_stretch_hom_site_unravelled;
    vector<float> _cached_hom_stretch_het_site_unravelled;
    vector<float> _cached_hom_site_hom_stretch_unravelled;
    vector<float> _cached_het_site_hom_stretch_unravelled;

    float* _cached_missing_flat;
    float* _cached_hom_flat;
    float* _cached_hom_stretch_hom_site_flat;
    float* _cached_hom_stretch_het_site_flat;
    float* _cached_hom_site_hom_stretch_flat;
    float* _cached_het_site_hom_stretch_flat;

    float _mean_min_log10;
    float _mean_max_log10;
    int _mean_n_steps;
    float _mean_step;

    float _mean_step_recip;

    float _cv_min_log10;
    float _cv_max_log10;
    int _cv_n_steps;
    float _cv_step;

    float _cv_step_recip;

    float _min_alpha;
    float _min_beta;

    bool _log_coords;

    // Interval variables
    float _m, _c; 
    int _m0, _c0;
    float* _lookup_forward[4];
    float* _lookup_backward[4];
    float* _flat;
    float* _ptr;
    float _U00, _U01, _U10, _U11, _V00, _V01, _V10, _V11;
    float _wm, _wc, _w00, _w01, _w10, _w11;

    __m256 _mean_step_recip_v;
    __m256 _mean_min_log10_times_mean_step_recip_v;
    __m256 _cv_step_recip_v;
    __m256 _cv_min_log10_times_cv_step_recip_v;

    __m256i _shift_n_step_coef;
    __m256i _shift_m0_coef;
    __m256i _shift_c0_coef;

    __m256 _wm_v;
    __m256 _wm_mult;
    __m256 _wm_add;

    __m256 _wc_v;
    __m256 _wc_mult;
    __m256 _wc_add;

    __m256 _wmvec;
    __m256 _wcvec;
    __m256 _halfdots;    

    FlowFieldCache(
        float scaled_recombination_rate,
        float scaled_mutation_rate,
        unique_ptr<FlowField> flow_field,
        int n_steps,
        bool log_coords = true
    ) : 
        _scaled_recombination_rate(scaled_recombination_rate),
        _scaled_mutation_rate(scaled_mutation_rate), 
        _flow_field(move(flow_field)),
        _n_steps(n_steps),
        _mean_min_log10(_flow_field->_mean_min_log10),
        _mean_max_log10(_flow_field->_mean_max_log10),
        _mean_n_steps(_flow_field->_mean_n_steps),
        _mean_step(_flow_field->_mean_step),
        _mean_step_recip(1.0 / _mean_step),
        _cv_min_log10(_flow_field->_cv_min_log10),
        _cv_max_log10(_flow_field->_cv_max_log10),
        _cv_n_steps(_flow_field->_cv_n_steps),
        _cv_step(_flow_field->_cv_step),
        _cv_step_recip(1.0 / _cv_step),
        _min_alpha(0),
        _min_beta(-10),
        _log_coords(log_coords)
    {
        int cache_size = _flow_field->_flow_field_unravelled.size() * _n_steps;

        _cached_missing_unravelled.resize(cache_size, 0);
        _cached_hom_unravelled.resize(cache_size, 0);
        _cached_hom_stretch_hom_site_unravelled.resize(cache_size, 0);
        _cached_hom_stretch_het_site_unravelled.resize(cache_size, 0);
        _cached_hom_site_hom_stretch_unravelled.resize(cache_size, 0);
        _cached_het_site_hom_stretch_unravelled.resize(cache_size, 0);

        size_t flat_n_elements = (
            _mean_n_steps * _cv_n_steps *           // Flow field grid size
            _n_steps *                              // Number of cached steps
            2 *                                     // Two grids (alpha and beta)
            4                                       // 4 grid points to interpolate from
            );

        _cached_missing_flat = aligned_alloc_float(flat_n_elements);
        _cached_hom_flat = aligned_alloc_float(flat_n_elements);
        _cached_hom_stretch_hom_site_flat = aligned_alloc_float(flat_n_elements);
        _cached_hom_stretch_het_site_flat = aligned_alloc_float(flat_n_elements);
        _cached_hom_site_hom_stretch_flat = aligned_alloc_float(flat_n_elements);
        _cached_het_site_hom_stretch_flat = aligned_alloc_float(flat_n_elements);

        _lookup_forward[0] = _cached_missing_flat;
        _lookup_forward[1] = _cached_hom_stretch_hom_site_flat;
        _lookup_forward[2] = _cached_hom_stretch_het_site_flat;
        _lookup_forward[3] = _cached_hom_flat;
        
        _lookup_backward[0] = _cached_missing_flat;
        _lookup_backward[1] = _cached_hom_site_hom_stretch_flat;
        _lookup_backward[2] = _cached_het_site_hom_stretch_flat;
        _lookup_backward[3] = _cached_hom_flat;

        _wm_mult = _mm256_setr_ps(-1, -1, 1, 1, -1, -1, 1, 1);
        _wm_add = _mm256_setr_ps(1, 1, 0, 0, 1, 1, 0, 0);
        _wc_mult = _mm256_setr_ps(-1, 1, -1, 1, -1, 1, -1, 1);
        _wc_add = _mm256_setr_ps(1, 0, 1, 0, 1, 0, 1, 0);

        _mean_step_recip_v = _mm256_set1_ps(_mean_step_recip);
        _mean_min_log10_times_mean_step_recip_v = _mm256_set1_ps(_mean_min_log10 * _mean_step_recip);
        _cv_step_recip_v = _mm256_set1_ps(_cv_step_recip);
        _cv_min_log10_times_cv_step_recip_v = _mm256_set1_ps(_cv_min_log10 * _cv_step_recip);
        
        _shift_n_step_coef = _mm256_set1_epi32(2 * 4 * _mean_n_steps * _cv_n_steps);
        _shift_m0_coef =_mm256_set1_epi32(2 * 4 * _cv_n_steps);
        _shift_c0_coef = _mm256_set1_epi32(2 * 4);

        preprocess();           
    }

    virtual ~FlowFieldCache() {
        free(_cached_missing_flat);
        free(_cached_hom_flat);
        free(_cached_hom_stretch_hom_site_flat);
        free(_cached_hom_stretch_het_site_flat);
        free(_cached_hom_site_hom_stretch_flat);
        free(_cached_het_site_hom_stretch_flat);
    }

    // TODO: Replace with a robust logaddexp, or better of, a cache
    float add_const(float x, float const_to_add) {
        return std::log10(std::pow(10, x) + const_to_add);
    }

    // void mean_cv_to_alpha_beta(
    //     float mean_log10, 
    //     float cv_log10,
    //     float& alpha,
    //     float& beta
    // ) {
    //     // mean_log10 = log10(alpha / beta) = log10(alpha) - log10(beta)
    //     // cv_log10 = log10(1 / sqrt(alpha)) = -0.5 * log10(alpha)
    //     // -->
    //     // log10(alpha) = cv_log10 / -0.5
    //     // log10(beta) = log10(alpha) - mean_log10

    //     float alpha_log10 = cv_log10 / -0.5;
    //     float beta_log10 = alpha_log10 - mean_log10;
    //     alpha = std::pow(10, alpha_log10);
    //     beta = std::pow(10, beta_log10);
    // }

    void mean_cv_to_alpha_beta_log10(
        float mean_log10, 
        float cv_log10,
        float& alpha_log10,
        float& beta_log10
    ) {
        // mean_log10 = log10(alpha / beta) = log10(alpha) - log10(beta)
        // cv_log10 = log10(1 / sqrt(alpha)) = -0.5 * log10(alpha)
        // -->
        // log10(alpha) = cv_log10 / -0.5
        // log10(beta) = log10(alpha) - mean_log10

        alpha_log10 = cv_log10 / -0.5;
        beta_log10 = alpha_log10 - mean_log10;
    }

    void alpha_beta_to_mean_cv_log10(
        float alpha_log10,
        float beta_log10,
        float& mean_log10, 
        float& cv_log10
    ) {
        // mean_log10 = log10(alpha / beta) = log10(alpha) - log10(beta)
        // cv_log10 = log10(1 / sqrt(alpha)) = -0.5 * log10(alpha)
        // -->
        // log10(alpha) = cv_log10 / -0.5
        // log10(beta) = log10(alpha) - mean_log10
        mean_log10 = alpha_log10 - beta_log10;
        cv_log10 = -0.5 * alpha_log10;
    }    

    // void clip_alpha_beta_log10(
    //     float alpha_log10, 
    //     float beta_log10,
    //     float& clipped_alpha_log10,
    //     float& clipped_beta_log10        
    // ) {
    //     float mean_log10 = alpha_log10 - beta_log10;
    //     float cv_log10 = -0.5 * alpha_log10;

    //     mean_log10 = min(max(mean_log10, _mean_min_log10), _mean_max_log10);
    //     cv_log10 = min(max(cv_log10, _cv_min_log10), _cv_max_log10);  

    //     mean_cv_to_alpha_beta_log10(mean_log10, cv_log10, clipped_alpha_log10, clipped_beta_log10);
    // }

    void clip_mean_cv_log10(
        float mean_log10, 
        float cv_log10,
        float& clipped_mean_log10,
        float& clipped_cv_log10        
    ) {
        clipped_mean_log10 = min(max(mean_log10, _mean_min_log10), _mean_max_log10);
        clipped_cv_log10 = min(max(cv_log10, _cv_min_log10), _cv_max_log10);  
    }

    int unravel_entry(int row, int col, int n_step, int component) {
        return ((2 * n_step + component) * _mean_n_steps + row) * _cv_n_steps + col;        
        
        // component should be 0 or 1
        // return (
        //     _cv_n_steps * _mean_n_steps * 2 * n_step +          // Matrices before this step
        //     _cv_n_steps * _mean_n_steps * component +           // First or second matrix for this step?
        //     _cv_n_steps * row +                                 // How many rows before the entry
        //     col
        // );
    }

    // void get_missing(int row, int col, int n_step, float& alpha, float& beta) {
    //     alpha = _cached_missing_unravelled[unravel_entry(row, col, n_step, 0)];
    //     beta = _cached_missing_unravelled[unravel_entry(row, col, n_step, 1)];
    // }

    // void get_hom(int row, int col, int n_step, float& alpha, float& beta) {
    //     alpha = _cached_hom_unravelled[unravel_entry(row, col, n_step, 0)];
    //     beta = _cached_hom_unravelled[unravel_entry(row, col, n_step, 1)];
    // }

    void get_unravelled(const vector<float>& vec, int row, int col, int n_step, float& alpha, float& beta) {
        alpha = vec[unravel_entry(row, col, n_step, 0)];
        beta = vec[unravel_entry(row, col, n_step, 1)];
    }

    void preprocess() {
        float mean_log10, cv_log10;
        float conv_mean_log10, conv_cv_log10;
        float updated_mean_log10, updated_cv_log10;
        float alpha_log10, beta_log10;

        //float alpha, beta;
        float u, v;
        //float clipped_alpha, clipped_beta;

        // The first step builds off the flow field
        mean_log10 = _mean_min_log10;
        for (int row = 0; row < _mean_n_steps; row++) {
            cv_log10 = _cv_min_log10;
            for (int col = 0; col < _cv_n_steps; col++) {

                // // Convert to alpha, beta
                // mean_cv_to_alpha_beta_log10(mean_log10, cv_log10, alpha, beta);

                // Get update coefficients
                _flow_field->_entry(row, col, u, v);

                // // Updated alpha beta
                // alpha = max(alpha + u * _scaled_recombination_rate, _min_alpha);
                // beta = max(beta + v * _scaled_recombination_rate, _min_beta);

                // clip_alpha_beta_log10(alpha, beta, clipped_alpha, clipped_beta);

                clip_mean_cv_log10(
                    mean_log10 + u * _scaled_recombination_rate,
                    cv_log10 + v * _scaled_recombination_rate,
                    updated_mean_log10,
                    updated_cv_log10
                );

                // Enter in the missing
                _cached_missing_unravelled[unravel_entry(row, col, 0, 0)] = updated_mean_log10; //clipped_alpha;
                _cached_missing_unravelled[unravel_entry(row, col, 0, 1)] = updated_cv_log10; //clipped_beta;                

                // Enter in the hom
                _cached_hom_unravelled[unravel_entry(row, col, 0, 0)] = updated_mean_log10; //clipped_alpha;
                _cached_hom_unravelled[unravel_entry(row, col, 0, 1)] = updated_cv_log10; //clipped_beta;

                cv_log10 += _cv_step;
            }            
            mean_log10 += _mean_step;
        }


        for (int n_step = 1; n_step < _n_steps; n_step++) {
            for (int row = 0; row < _mean_n_steps; row++) {
                for (int col = 0; col < _cv_n_steps; col++) {
                    // 
                    // For missing:
                    //

                    // Take the previous alpha, beta
                    // alpha = _cached_missing_unravelled[unravel_entry(row, col, n_step-1, 0)];
                    // beta  = _cached_missing_unravelled[unravel_entry(row, col, n_step-1, 1)];
                    mean_log10 = _cached_missing_unravelled[unravel_entry(row, col, n_step-1, 0)];
                    cv_log10   = _cached_missing_unravelled[unravel_entry(row, col, n_step-1, 1)];

                    // This is "missing", so no observation update, just next step
                    // Get update coefficients
                    _flow_field->at(mean_log10, cv_log10, u, v);
                    // _flow_field->at_log10_alpha_beta(alpha, beta, u, v);

                    // Updated alpha beta
                    // alpha = max(alpha + u * _scaled_recombination_rate, _min_alpha);
                    // beta = max(beta + v * _scaled_recombination_rate, _min_beta);
                    // clip_alpha_beta_log10(alpha, beta, clipped_alpha, clipped_beta);
                    clip_mean_cv_log10(
                        mean_log10 + u * _scaled_recombination_rate,
                        cv_log10 + v * _scaled_recombination_rate,
                        updated_mean_log10,
                        updated_cv_log10
                    );

                     // Enter in the missing
                    // _cached_missing_unravelled[unravel_entry(row, col, n_step, 0)] = clipped_alpha;
                    // _cached_missing_unravelled[unravel_entry(row, col, n_step, 1)] = clipped_beta;
                    _cached_missing_unravelled[unravel_entry(row, col, n_step, 0)] = updated_mean_log10;
                    _cached_missing_unravelled[unravel_entry(row, col, n_step, 1)] = updated_cv_log10;

                    //
                    // For hom
                    //

                    // Take the previous alpha, beta
                    // alpha = _cached_hom_unravelled[unravel_entry(row, col, n_step-1, 0)];
                    // beta  = _cached_hom_unravelled[unravel_entry(row, col, n_step-1, 1)];
                    mean_log10 = _cached_hom_unravelled[unravel_entry(row, col, n_step-1, 0)];
                    cv_log10   = _cached_hom_unravelled[unravel_entry(row, col, n_step-1, 1)];

                    // Convert to alpha, beta, update and back
                    mean_cv_to_alpha_beta_log10(
                        mean_log10,
                        cv_log10,
                        alpha_log10,
                        beta_log10
                    );

                    // Update hom observation
                    beta_log10 = add_const(beta_log10, 2 * _scaled_mutation_rate);

                    alpha_beta_to_mean_cv_log10(                        
                        alpha_log10,
                        beta_log10,
                        mean_log10,
                        cv_log10
                    );

                    clip_mean_cv_log10(
                        mean_log10,
                        cv_log10,
                        conv_mean_log10,
                        conv_cv_log10
                    );

                    // Get update coefficients
                    _flow_field->at(conv_mean_log10, conv_cv_log10, u, v);

                    // Updated alpha beta
                    // alpha = max(alpha + u * _scaled_recombination_rate, _min_alpha);
                    // beta = max(beta + v * _scaled_recombination_rate, _min_beta);

                    // clip_alpha_beta_log10(alpha, beta, clipped_alpha, clipped_beta);
                    clip_mean_cv_log10(
                        conv_mean_log10 + u * _scaled_recombination_rate,
                        conv_cv_log10 + v * _scaled_recombination_rate,
                        updated_mean_log10,
                        updated_cv_log10
                    );

                    // Enter in the hom
                    _cached_hom_unravelled[unravel_entry(row, col, n_step, 0)] = updated_mean_log10; //clipped_alpha;
                    _cached_hom_unravelled[unravel_entry(row, col, n_step, 1)] = updated_cv_log10; //clipped_beta;

                }                
            }
        }

        flatten(_cached_missing_unravelled, _cached_missing_flat);
        flatten(_cached_hom_unravelled,     _cached_hom_flat);

        // Now go over it again and add an emission step after
        for (int n_step = 0; n_step < _n_steps; n_step++) {
            for (int row = 0; row < _mean_n_steps; row++) {
                for (int col = 0; col < _cv_n_steps; col++) {                   
                    // alpha = _cached_hom_unravelled[unravel_entry(row, col, n_step, 0)];
                    // beta  = _cached_hom_unravelled[unravel_entry(row, col, n_step, 1)];
                    mean_log10 = _cached_hom_unravelled[unravel_entry(row, col, n_step, 0)];
                    cv_log10   = _cached_hom_unravelled[unravel_entry(row, col, n_step, 1)];

                    // Update hom observation
                    // beta = add_const(beta, _scaled_mutation_rate);
                    // Convert to alpha, beta, update and back
                    mean_cv_to_alpha_beta_log10(
                        mean_log10,
                        cv_log10,
                        alpha_log10,
                        beta_log10
                    );

                    // Update hom observation
                    beta_log10 = add_const(beta_log10, 2 * _scaled_mutation_rate);

                    alpha_beta_to_mean_cv_log10(                        
                        alpha_log10,
                        beta_log10,
                        mean_log10,
                        cv_log10
                    );

                    // clip_alpha_beta_log10(alpha, beta, clipped_alpha, clipped_beta);
                    clip_mean_cv_log10(
                        mean_log10,
                        cv_log10,
                        updated_mean_log10,
                        updated_cv_log10
                    );

                    // Enter
                    _cached_hom_stretch_hom_site_unravelled[unravel_entry(row, col, n_step, 0)] = updated_mean_log10; //clipped_alpha;
                    _cached_hom_stretch_hom_site_unravelled[unravel_entry(row, col, n_step, 1)] = updated_cv_log10; //clipped_beta;

                    // Update het observation
                    // alpha = add_const(alpha, 1);
                    alpha_log10 = add_const(alpha_log10, 1);

                    alpha_beta_to_mean_cv_log10(                        
                        alpha_log10,
                        beta_log10,
                        mean_log10,
                        cv_log10
                    );

                    //clip_alpha_beta_log10(alpha, beta, clipped_alpha, clipped_beta); 
                    clip_mean_cv_log10(
                        mean_log10,
                        cv_log10,
                        updated_mean_log10,
                        updated_cv_log10
                    );

                    // Enter
                    _cached_hom_stretch_het_site_unravelled[unravel_entry(row, col, n_step, 0)] = updated_mean_log10; //clipped_alpha;
                    _cached_hom_stretch_het_site_unravelled[unravel_entry(row, col, n_step, 1)] = updated_cv_log10; //clipped_beta;
                }                
            }
        }

        flatten(_cached_hom_stretch_hom_site_unravelled, _cached_hom_stretch_hom_site_flat);
        flatten(_cached_hom_stretch_het_site_unravelled, _cached_hom_stretch_het_site_flat);

        // Go over it again and add an emission step before
        for (int n_step = 0; n_step < _n_steps; n_step++) {
            mean_log10 = _mean_min_log10;
            for (int row = 0; row < _mean_n_steps; row++) {
                cv_log10 = _cv_min_log10;
                for (int col = 0; col < _cv_n_steps; col++) {  

                    mean_cv_to_alpha_beta_log10(mean_log10, cv_log10, alpha_log10, beta_log10);

                    // if ((n_step == 999) && (row == 5) && (col >= 8)) {
                    //     cout << boost::format("*** %d %d mean_log10, cv_log10 alpha_log10, beta_log10) = %f %f %f %f \n")
                    //         % row % col % mean_log10 % cv_log10 % alpha_log10 % beta_log10;
                    // }

                    // hom 
                    //beta = add_const(beta, _scaled_mutation_rate);
                    beta_log10 = add_const(beta_log10, 2 * _scaled_mutation_rate);

                    // if ((n_step == 999) && (row == 5) && (col >= 8)) {
                    //     cout << boost::format("*** %d %d alpha_log10, beta_log10 = %f %f\n")
                    //         % row % col % alpha_log10 % beta_log10;
                    // } 

                    //clip_alpha_beta_log10(alpha, beta, clipped_alpha, clipped_beta); 
                    alpha_beta_to_mean_cv_log10(                        
                        alpha_log10,
                        beta_log10,
                        conv_mean_log10,
                        conv_cv_log10
                    );

                    // if ((n_step == 999) && (row == 5) && (col >= 8)) {
                    //     cout << boost::format("*** conv %d %d %f %f\n")
                    //         % row % col % conv_mean_log10 % conv_cv_log10;
                    // }                 


                    //clip_alpha_beta_log10(alpha, beta, clipped_alpha, clipped_beta); 
                    clip_mean_cv_log10(
                        conv_mean_log10,
                        conv_cv_log10,
                        updated_mean_log10,
                        updated_cv_log10
                    );

                    // if ((n_step == 999) && (row == 5) && (col>= 8)) {
                    //     cout << boost::format("*** updated %d %d %f %f\n")
                    //         % row % col % updated_mean_log10 % updated_cv_log10;
                    // } 



                    // at_log10_alpha_beta(
                    //     clipped_alpha, clipped_beta, n_step, HOM_STRETCH,                         
                    //     _cached_hom_site_hom_stretch_unravelled[unravel_entry(row, col, n_step, 0)], 
                    //     _cached_hom_site_hom_stretch_unravelled[unravel_entry(row, col, n_step, 1)]
                    // );

                    // at(
                    //     updated_mean_log10, updated_cv_log10, n_step, HOM_STRETCH,                         
                    //     _cached_hom_site_hom_stretch_unravelled[unravel_entry(row, col, n_step, 0)], 
                    //     _cached_hom_site_hom_stretch_unravelled[unravel_entry(row, col, n_step, 1)]
                    // );

                    at_flat(
                        updated_mean_log10, updated_cv_log10, n_step, HOM_STRETCH, false, u, v
                    );

                    _cached_hom_site_hom_stretch_unravelled[unravel_entry(row, col, n_step, 0)] = u;
                    _cached_hom_site_hom_stretch_unravelled[unravel_entry(row, col, n_step, 1)] = v;

                    // if ((n_step == 999) && (row == 5) && (col >= 8)) {
                    //     cout << boost::format("*** %d %d %f %f at %d %d\n")
                    //         % row % col 
                    //         % _cached_hom_site_hom_stretch_unravelled[unravel_entry(row, col, n_step, 0)] 
                    //         % _cached_hom_site_hom_stretch_unravelled[unravel_entry(row, col, n_step, 1)]
                    //         % unravel_entry(row, col, n_step, 0)
                    //         % unravel_entry(row, col, n_step, 1);
                    // } 

                    // het
                    //alpha = add_const(alpha, 1);
                    alpha_log10 = add_const(alpha_log10, 1);

                    //  if ((n_step == 999) && (row == 5) && (col == 9)) {
                    //     cout << boost::format("*** %d %d %f %f\n")
                    //         % row % col % alpha_log10 % beta_log10;
                    // } 

                    alpha_beta_to_mean_cv_log10(                        
                        alpha_log10,
                        beta_log10,
                        conv_mean_log10,
                        conv_cv_log10
                    );

                    // if ((n_step == 999) && (row == 5) && (col == 9)) {
                    //     cout << boost::format("*** conv %d %d %f %f\n")
                    //         % row % col % conv_mean_log10 % conv_cv_log10;
                    // } 


                    //clip_alpha_beta_log10(alpha, beta, clipped_alpha, clipped_beta); 
                    clip_mean_cv_log10(
                        conv_mean_log10,
                        conv_cv_log10,
                        updated_mean_log10,
                        updated_cv_log10
                    );

                    // if ((n_step == 999) && (row == 5) && (col == 9)) {
                    //     cout << boost::format("*** updated %d %d %f %f\n")
                    //         % row % col % updated_mean_log10 % updated_cv_log10;
                    // } 


                    // at_log10_alpha_beta(
                    //     clipped_alpha, clipped_beta, n_step, HOM_STRETCH,                         
                    //     _cached_het_site_hom_stretch_unravelled[unravel_entry(row, col, n_step, 0)], 
                    //     _cached_het_site_hom_stretch_unravelled[unravel_entry(row, col, n_step, 1)]
                    // );
                    at_flat(
                        updated_mean_log10, updated_cv_log10, n_step, HOM_STRETCH, false, u, v
                    );                
                        
                    _cached_het_site_hom_stretch_unravelled[unravel_entry(row, col, n_step, 0)] = u;
                    _cached_het_site_hom_stretch_unravelled[unravel_entry(row, col, n_step, 1)] = v;
                    
                    cv_log10 += _cv_step;
                }                
                mean_log10 += _mean_step;
            }
        }

        flatten(_cached_hom_site_hom_stretch_unravelled, _cached_hom_site_hom_stretch_flat);
        flatten(_cached_het_site_hom_stretch_unravelled, _cached_het_site_hom_stretch_flat);
    }

    void flatten(const vector<float>& unravelled, float* flat) {
        for (int n_step = 0; n_step < _n_steps; n_step++) {
            for (int row = 0; row < _mean_n_steps; row++) {
                for (int col = 0; col < _cv_n_steps; col++) {                   

                    int m0 = row;
                    if (m0 == _mean_n_steps-1) {
                        m0--;
                    }
                    int m1 = m0 + 1;
                    int c0 = col;
                    if (c0 == _cv_n_steps-1) {
                        c0--;
                    }
                    int c1 = c0 + 1;

                    float* ptr = (flat + 
                        n_step * _mean_n_steps * _cv_n_steps * 2 * 4 +          // all parts before it
                        row * _cv_n_steps * 2 * 4 +                             // rows before it
                        col * 2 * 4                                             // cols before it
                    );

                    *(ptr + 0) = unravelled[unravel_entry(m0, c0, n_step, 0)];
                    *(ptr + 1) = unravelled[unravel_entry(m0, c1, n_step, 0)];
                    *(ptr + 2) = unravelled[unravel_entry(m1, c0, n_step, 0)];
                    *(ptr + 3) = unravelled[unravel_entry(m1, c1, n_step, 0)];
                    *(ptr + 4) = unravelled[unravel_entry(m0, c0, n_step, 1)];
                    *(ptr + 5) = unravelled[unravel_entry(m0, c1, n_step, 1)];
                    *(ptr + 6) = unravelled[unravel_entry(m1, c0, n_step, 1)];
                    *(ptr + 7) = unravelled[unravel_entry(m1, c1, n_step, 1)];                    
                }                
            }
        }
    }


    void at_flat(float mean_log10, float cv_log10, int n_step, segment_type type, bool forward, float& u, float& v, bool printme = false) {  

        // TODO: Do those two in parallel      
        _m = (mean_log10 - _mean_min_log10) * _mean_step_recip;
        _c = (cv_log10 - _cv_min_log10) * _cv_step_recip;

        // _m0 = int(floor(_m));
        // _c0 = int(floor(_c));

        // TODO: Do those two in parallel
        _m0 = ftoi_sse1(_m);
        _c0 = ftoi_sse1(_c);

        

        _flat = (forward ? _lookup_forward[type] : _lookup_backward[type]);

        _ptr = (_flat + 2 * 4 * (n_step * _mean_n_steps * _cv_n_steps + _m0 * _cv_n_steps + _c0));        

        _wm = (_m - _m0);
        _wc = (_c - _c0);

        // if (printme) { 
        //     cout << boost::format("n_step = %d, _m0 = %d, _c0 = %d\n") % n_step % _m0 % _c0; 
        //     cout << boost::format("_flat = %d, vs. _cached_hom_site_hom_stretch_flat %d\n") % _flat % _cached_hom_site_hom_stretch_flat;
        //     cout << boost::format("_flat = %d, vs. _cached_hom_flat %d\n") % _flat % _cached_hom_flat;
        //     cout << boost::format("_ptr - _flat = %d\n") % (_ptr - _flat);
        //     cout << boost::format("_ptr content %f %f %f %f %f %f %f %f\n") % 
        //         _ptr[0] % _ptr[1] % _ptr[2] % _ptr[3] % _ptr[4] % _ptr[5] % _ptr[6] % _ptr[7];
        // }

        // // Fill the four grid points
        // _U00 = *(_ptr+0);
        // _U01 = *(_ptr+1);
        // _U10 = *(_ptr+2); 
        // _U11 = *(_ptr+3);        
        // _V00 = *(_ptr+4); 
        // _V01 = *(_ptr+5);
        // _V10 = *(_ptr+6); 
        // _V11 = *(_ptr+7);


        // _w00 = (1-_wm) * (1-_wc);
        // _w01 = (1-_wm) * _wc;
        // _w10 = _wm * (1-_wc);
        // _w11 = _wm * _wc;

        // u = _w00 * _U00 + _w01 * _U01 + _w10 * _U10 + _w11 * _U11;
        // v = _w00 * _V00 + _w01 * _V01 + _w10 * _V10 + _w11 * _V11;

        // AVX2 version


        // __m256 UV = _mm256_load_ps(_ptr);
        // cout << (boost::format("Us: %f %f %f %f\n") % _U00 % _U01 % _U10 % _U11);
        // cout << (boost::format("Vs: %f %f %f %f\n") % _V00 % _V01 % _V10 % _V11);

        //  cout << UV[0] << " " << UV[1] << " " << UV[2] << " " << UV[3] << " " 
        //       << UV[4] << " " << UV[5] << " " << UV[6] << " " << UV[7] << " " << endl;


        // cout << u << " / " << v << endl;
        // cout << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << " " 
        //      << res[4] << " " << res[5] << " " << res[6] << " " << res[7] << " " << endl;

        // TODO: can we save something here with FMA?
        __m256 wmvec = _mm256_fmadd_ps(_mm256_set1_ps(_wm), _wm_mult, _wm_add);
        __m256 wcvec = _mm256_fmadd_ps(_mm256_set1_ps(_wc), _wc_mult, _wc_add);
        __m256 halfdots = _mm256_dp_ps(
            _mm256_load_ps(_ptr), 
            _mm256_mul_ps(
                wmvec, wcvec                
            ), 
            0xF1);

        float* res = (float*)&halfdots;

        // if (printme) {
        //     cout << boost::format("_wm _wc = %f %f\n") % _wm % _wc;
        //     cout << boost::format("res content %f %f %f %f %f %f %f %f\n") % 
        //         res[0] % res[1] % res[2] % res[3] % res[4] % res[5] % res[6] % res[7];
        // }

        u = res[0]; v = res[4];
        

        // cout << endl << "at_flat" << endl;
        // cout << (boost::format("In: %f %f %d %d %b\n") % mean_log10 % cv_log10 % n_step % type % forward);
        // cout << (boost::format("%f %f %f %f \n") % m % m0 % c % c0);
        // cout << (boost::format("Us: %f %f %f %f\n") % U00 % U01 % U10 % U11);
        // cout << (boost::format("Vs: %f %f %f %f\n") % V00 % V01 % V10 % V11);

        // int offset = 2 * 4 * (n_step * _mean_n_steps * _cv_n_steps + m0 * _cv_n_steps + c0);
        // cout << flat << " vs " << _cached_hom_stretch_hom_site_flat << endl;
        // cout << offset << endl;
        // cout << *(_cached_hom_stretch_hom_site_flat+offset+5) << endl;
        // get_unravelled(_cached_hom_stretch_hom_site_unravelled, 5, 9, n_step, U01, V01);
        // cout << V01 << endl;
    }

    void at_flat_vectorized(
        float* mean_log10_chunk, 
        float* cv_log10_chunk, 
        int32_t n_step, 
        int32_t* n_called_ptr,
        segment_type* type_chunk, 
        bool forward) {  

        //
        // First step: Missing segments
        //

        // Load
        __m256 mean_log10_v = _mm256_load_ps(mean_log10_chunk);
        __m256 cv_log10_v = _mm256_load_ps(cv_log10_chunk);

        // _m = (mean_log10 - _mean_min_log10) * _mean_step_recip;        
        __m256 _m = _mm256_fmsub_ps(mean_log10_v, _mean_step_recip_v, _mean_min_log10_times_mean_step_recip_v);

        // _c = (cv_log10 - _cv_min_log10) * _cv_step_recip;
        __m256 _c = _mm256_fmsub_ps(cv_log10_v, _cv_step_recip_v, _cv_min_log10_times_cv_step_recip_v);

        // _m0 = floor(_m);
        // _c0 = floor(_c);
        __m256 _m0 = _mm256_floor_ps(_m);
        __m256 _c0 = _mm256_floor_ps(_c);  

        // _wm = (_m - _m0);
        // _wc = (_c - _c0);
        _wm_v = _mm256_sub_ps(_m, _m0);
        _wc_v = _mm256_sub_ps(_c, _c0);

        __m256i n_missing = _mm256_sub_epi32(_mm256_set1_epi32(n_step), _mm256_load_si256((__m256i*) n_called_ptr));

        // shifts = 2 * 4 * (n_missing * _mean_n_steps * _cv_n_steps + _m0 * _cv_n_steps + _c0)
        __m256i shifts_missing = _mm256_add_epi32(_mm256_add_epi32(
            _mm256_mullo_epi32(_mm256_sub_epi32(n_missing, _mm256_set1_epi32(1)), _shift_n_step_coef),
            _mm256_mullo_epi32(_mm256_cvtps_epi32(_m0), _shift_m0_coef)),
            _mm256_mullo_epi32(_mm256_cvtps_epi32(_c0), _shift_c0_coef)
        );


        for (int i = 0; i < parallel_vector_size; i++) {
            // If missing steps is 0, do nothing
            if (((int32_t*) (&n_missing))[i] == 0) {
                continue;
            }

            // _flat = (forward ? _lookup_forward[type] : _lookup_backward[type]);

            // This is always a missing stretch
            _flat = _cached_missing_flat;

            //_ptr = (_flat + 2 * 4 * (n_step * _mean_n_steps * _cv_n_steps + _m0 * _cv_n_steps + _c0));
            _ptr = _flat + ((int32_t*) (&shifts_missing))[i];

            _wm = ((float*) (&_wm_v))[i];
            _wc = ((float*) (&_wc_v))[i];

            // TODO: can we save something here with FMA?
            _wmvec = _mm256_fmadd_ps(_mm256_set1_ps(_wm), _wm_mult, _wm_add);
            _wcvec = _mm256_fmadd_ps(_mm256_set1_ps(_wc), _wc_mult, _wc_add);
            _halfdots = _mm256_dp_ps(
                _mm256_load_ps(_ptr), 
                _mm256_mul_ps(
                    _wmvec, _wcvec                
                ), 
                0xF1);

            float* res = (float*)&_halfdots;

            *(mean_log10_chunk+i) = res[0]; 
            *(cv_log10_chunk+i) = res[4];
        }

        // 
        // Second step: Hom stretches, with hom/het/missing emission
        //
        mean_log10_v = _mm256_load_ps(mean_log10_chunk);
        cv_log10_v = _mm256_load_ps(cv_log10_chunk);

        // _m = (mean_log10 - _mean_min_log10) * _mean_step_recip;        
        _m = _mm256_fmsub_ps(mean_log10_v, _mean_step_recip_v, _mean_min_log10_times_mean_step_recip_v);

        // _c = (cv_log10 - _cv_min_log10) * _cv_step_recip;
        _c = _mm256_fmsub_ps(cv_log10_v, _cv_step_recip_v, _cv_min_log10_times_cv_step_recip_v);

        // _m0 = floor(_m);
        // _c0 = floor(_c);
        _m0 = _mm256_floor_ps(_m);
        _c0 = _mm256_floor_ps(_c);  

        // _wm = (_m - _m0);
        // _wc = (_c - _c0);
        _wm_v = _mm256_sub_ps(_m, _m0);
        _wc_v = _mm256_sub_ps(_c, _c0);

        __m256i n_called = _mm256_load_si256((__m256i*) n_called_ptr);

        // shifts = 2 * 4 * (n_step * _mean_n_steps * _cv_n_steps + _m0 * _cv_n_steps + _c0)
        __m256i shifts_called = _mm256_add_epi32(_mm256_add_epi32(
            _mm256_mullo_epi32(_mm256_sub_epi32(n_called, _mm256_set1_epi32(1)), _shift_n_step_coef),
            _mm256_mullo_epi32(_mm256_cvtps_epi32(_m0), _shift_m0_coef)),
            _mm256_mullo_epi32(_mm256_cvtps_epi32(_c0), _shift_c0_coef)
        );


        for (int i = 0; i < parallel_vector_size; i++) {
            // If called steps is 0, do nothing
            if (((int32_t*) (&n_called))[i] == 0) {
                continue;
            }

            // _flat = (forward ? _lookup_forward[type] : _lookup_backward[type]);

            // TODO: Make this faster
            _flat = (
                forward 
                ? _lookup_forward[*(type_chunk + i)]
                : _lookup_backward[*(type_chunk + i)]
            );

            //_ptr = (_flat + 2 * 4 * (n_step * _mean_n_steps * _cv_n_steps + _m0 * _cv_n_steps + _c0));
            _ptr = _flat + ((int32_t*) (&shifts_called))[i];

            _wm = ((float*) (&_wm_v))[i];
            _wc = ((float*) (&_wc_v))[i];

            // TODO: can we save something here with FMA?
            _wmvec = _mm256_fmadd_ps(_mm256_set1_ps(_wm), _wm_mult, _wm_add);
            _wcvec = _mm256_fmadd_ps(_mm256_set1_ps(_wc), _wc_mult, _wc_add);
            _halfdots = _mm256_dp_ps(
                _mm256_load_ps(_ptr), 
                _mm256_mul_ps(
                    _wmvec, _wcvec                
                ), 
                0xF1);

            float* res = (float*)&_halfdots;

            *(mean_log10_chunk+i) = res[0]; 
            *(cv_log10_chunk+i) = res[4];
        }
    }

};