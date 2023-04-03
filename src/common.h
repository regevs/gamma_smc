#pragma once

#include <iostream>
#include <fstream>
#include <utility>
#include <stdexcept>
#include <exception>
#include <memory>
#include <chrono>
#include <algorithm>
#include <unordered_map>
#include <filesystem>
#include <cstring>

#include <boost/assert.hpp>
#include <boost/format.hpp>

#include <math.h> 
#include <stdlib.h>

#include <zstd.h>

// TODO: AVX DEFINES
#include <immintrin.h>
#include <emmintrin.h>

#include "indicators.h"
using namespace indicators;

using namespace std;

typedef long position_t;


// TODO: Unclear if this is needed anymore so much. Now that each segment is
// a stretch of missing and of hom, the difference is just the emission.
enum segment_type : int8_t {
    MISSING_STRETCH_MISSING_SITE, 
    HOM_STRETCH_HOM_SITE, 
    HOM_STRETCH_HET_SITE,
    HOM_STRETCH                 // with missing site in the end; TODO: Change name to HOM_STRETCH_MISSING_SITE?
};

// TODO: Maybe unite those two, they are not so different
struct SegSite_t {
    position_t pos;
    segment_type type;  
    vector<int32_t> alleles;
};

struct SegregatingSite {
    position_t pos;
    vector<int8_t> alleles;
};


struct Segment_t {
    position_t pos;   // TODO: Make this start_pos and end_pos, or at least clarify which is it
    int length;    
    segment_type type;
    bool output_at_start;
    bool output_at_end;
    int seg_site_index;
};

// Memory and AVX

const int parallel_vector_size = 8;
const size_t parallel_vector_size_in_bytes_float = parallel_vector_size * sizeof(float);
const size_t memory_alignment = 64;

inline float* aligned_alloc_float(size_t n_elements, bool reset = true) {
    size_t alloc_size = n_elements * sizeof(float);
    float* ptr = (float*) aligned_alloc(memory_alignment, alloc_size);
    if (ptr == NULL) {
        cerr << "Error in allocation." << endl;
        exit(-1);
    }
    if (reset) memset((void *) ptr, 0, alloc_size);
    return ptr;
}


inline int32_t* aligned_alloc_int32(size_t n_elements, bool reset = true) {
    size_t alloc_size = n_elements * sizeof(int32_t);
    int32_t* ptr = (int32_t*) aligned_alloc(memory_alignment, alloc_size);
    if (ptr == NULL) {
        cerr << "Error in allocation." << endl;
        exit(-1);
    }
    if (reset) memset((void *) ptr, 0, alloc_size);
    return ptr;
}


inline int8_t* aligned_alloc_int8(size_t n_elements, bool reset = true) {
    size_t alloc_size = n_elements * sizeof(int8_t);
    int8_t* ptr = (int8_t*) aligned_alloc(memory_alignment, alloc_size);
    if (ptr == NULL) {
        cerr << "Error in allocation." << endl;
        exit(-1);
    }
    if (reset) memset((void *) ptr, 0, alloc_size);
    return ptr;
}