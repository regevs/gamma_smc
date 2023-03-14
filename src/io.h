#pragma once

#include "common.h"

#include <boost/range/adaptor/indexed.hpp>
#include <boost/assign.hpp>
#include <boost/algorithm/string.hpp>

#include <htslib/vcf.h>

// bool has_missing_data(const string& raw_alleles, const pair<int, int>& indices) {
//     if (!strchr("ACTG01", raw_alleles[indices.first])) {
//         return true;
//     }
//     if (!strchr("ACTG01", raw_alleles[indices.second])) {
//         return true;
//     }
//     return false;
// }

/*
// Clone of MSMC2
// TODO: Give credit
// TODO: Add support for ambiguous phasing?
void readSegSites(
    string filename, 
    const vector<pair<int, int>>& indices, 
    vector<vector<SegSite_t>>& ret
    ) {
    // format: chr position nr_calledSites [alleles] -> tab separated
    // if no alleles are given, assume M=2 and "01"
    // returns data for each pair of indices

    // Initialize the result
    ret.resize(0);
    ret.resize(indices.size(), vector<SegSite_t>());

    ifstream input(filename);

    position_t lastPos = -1;

    vector<string> fields;

    for (string line; getline(input, line);) {
        // Split into fields
        boost::algorithm::trim(line);
        boost::split(fields, line, boost::is_any_of("\t "));       
        position_t pos = stoi(fields[1]);
        long nrCalledSites = stoi(fields[2]);
        string raw_alleles = fields[3];
        
        // If this is the first row, then set the last position to 
        // be the number of called sites before current position

        // TODO: That seems wrong - it causes the first site to ignore the number of called sites?
        // But maybe it's better since it's often far from 0, and otherwise would be a strech of hom-s.
        if (lastPos == -1) {
            lastPos = pos - nrCalledSites;
        }

        BOOST_ASSERT_MSG(nrCalledSites <= pos - lastPos, "Number of called sites must be no more than the gap between consecutive positions");
        BOOST_ASSERT_MSG(nrCalledSites > 0, "Number of called sites must be positive");

        // Fill the data for each pair of indices
        int i = 0;
        for (const auto& ind : indices) {
            // If there is a missing allele at this position for this pair
            if (has_missing_data(raw_alleles, ind)) {

                // If there is missing data, encode it as a stretch of missing positions...
                if (nrCalledSites < pos - lastPos) {  // missing data
                    ret[i].push_back({pos - nrCalledSites, MISSING_STRETCH_MISSING_SITE});
                }
                // ...and then a stretch of hom positions.
                if (nrCalledSites > 1) {
                    ret[i].push_back({pos - 1, HOM_STRETCH_HOM_SITE});
                }
                // Finally, the actual observed position, which is missing.
                ret[i].push_back({pos, MISSING_STRETCH_MISSING_SITE});

                lastPos = pos;
            } 
            
            // Otherwise, no missing data at this position for this pair
            else {
                if (nrCalledSites < pos - lastPos) {  // missing data
                    ret[i].push_back({pos - nrCalledSites, MISSING_STRETCH_MISSING_SITE});
                }

                segment_type status = (raw_alleles[ind.first] == raw_alleles[ind.second] ? HOM_STRETCH_HOM_SITE : HOM_STRETCH_HET_SITE);
                ret[i].push_back({pos, status});
            }

            i++;
        }

        lastPos = pos;
    }
}*/

bool has_missing_data_all(const string& raw_alleles) {
    for (const char& c : raw_alleles) {
        if (!strchr("ACTG01", c)) {
            return true;
        }
    }
    return false;
}

void readSegSitesAll(
    string filename, 
    vector<SegSite_t>& ret
    ) {
    // format: chr position nr_calledSites [alleles] -> tab separated
    // if no alleles are given, assume M=2 and "01"
    // returns data for each pair of indices

    // Initialize the result
    ret.resize(0);

    ifstream input(filename);

    position_t lastPos = -1;

    vector<string> fields;

    for (string line; getline(input, line);) {
        // Split into fields
        boost::algorithm::trim(line);
        boost::split(fields, line, boost::is_any_of("\t "));       
        position_t pos = stoi(fields[1]);
        long nrCalledSites = stoi(fields[2]);
        string raw_alleles = fields[3];
        
        // If this is the first row, then set the last position to 
        // be the number of called sites before current position

        // TODO: That seems wrong - it causes the first site to ignore the number of called sites?
        // But maybe it's better since it's often far from 0, and otherwise would be a strech of hom-s.
        if (lastPos == -1) {
            lastPos = pos - nrCalledSites;
        }

        BOOST_ASSERT_MSG(nrCalledSites <= pos - lastPos, "Number of called sites must be no more than the gap between consecutive positions");
        BOOST_ASSERT_MSG(nrCalledSites > 0, "Number of called sites must be positive");

        vector<int32_t> int_alleles;
        for(char& c : raw_alleles) {
            if (c == '0') {
                int_alleles.push_back(0);
            } else if (c == '1') {
                int_alleles.push_back(1);
            } else {
                // TODO: Deal with this case
                //cout << "Warning: Only 0/1 alleles supported for now. Treating as 1" << endl;
                int_alleles.push_back(1);
            }
        }

        // If there is a missing allele at this position for this pair
        if (has_missing_data_all(raw_alleles)) {

            // If there is missing data, encode it as a stretch of missing positions...
            if (nrCalledSites < pos - lastPos) {  // missing data
                ret.push_back({pos - nrCalledSites, MISSING_STRETCH_MISSING_SITE, int_alleles});
            }
            // ...and then a stretch of hom positions.
            if (nrCalledSites > 1) {
                ret.push_back({pos - 1, HOM_STRETCH, int_alleles});
            }
            // Finally, the actual observed position, which is missing.
            ret.push_back({pos, MISSING_STRETCH_MISSING_SITE, int_alleles});

            lastPos = pos;
        } 
        
        // Otherwise, no missing data at this position for this pair
        else {
            if (nrCalledSites < pos - lastPos) {  // missing data
                ret.push_back({pos - nrCalledSites, MISSING_STRETCH_MISSING_SITE, int_alleles});
            }

            // We don't decide HOM_STRETCH_HOM_SITE or HOM_STRETCH_HET_SITE here
            ret.push_back({pos, HOM_STRETCH, int_alleles});
        }        

        lastPos = pos;
    }
}

// void chop_segsites(
//     const vector<SegSite_t>& segsites, 
//     long maxDistance,
//     vector<SegSite_t>& ret) 
//     {
//     ret.resize(0);

//     position_t lastPos = 0;
//     for (const auto& segsite : segsites) {
//         while (segsite.pos - lastPos > maxDistance) {
//             ret.push_back({lastPos + maxDistance, min(segsite.obs, 1)});
//             lastPos += maxDistance;
//         }
//         ret.push_back(segsite);
//         lastPos = segsite.pos;
//     }
// }

void read_flow_field_raw(
    const string& flow_field_filename,
    vector<float>& mean_grid_def,
    vector<float>& cv_grid_def,
    vector<float>& flow_field_unravelled
    ) 
    {
    
    ifstream indata;
	indata.open(flow_field_filename);
	
    string line;
	vector<float> values;
	
    uint rows = 0;
	while (getline(indata, line)) {
        stringstream lineStream(line);
        string cell;
        while (getline(lineStream, cell, ' ')) {
            if (rows == 0) {
                mean_grid_def.push_back(stof(cell));
            } else if (rows == 1) {
                cv_grid_def.push_back(stof(cell));
            } else {
                flow_field_unravelled.push_back(stof(cell));
            }
        }
        ++rows;
	}
}


// 
// Open vcf
//
void readVcf(
    string filename, 
    vector<unique_ptr<SegregatingSite>>& ret,
    vector<string>& samples_in_order,
    string samples_filename,
    vector<int>& samples_indices,
    string samples_filename_against,
    vector<int>& samples_against_indices
    ) {
    vector<string> samples_names;

    vector<string> required_samples;
    vector<string> samples_names_against;
    
    // Read samples list, if available
    if (samples_filename.size() > 0) {
        ifstream indata;
        indata.open(samples_filename);        
        string line;
        while (getline(indata, line)) {
            boost::algorithm::trim(line);
            required_samples.push_back(line);
        }            
    }

    // Read second samples list, if available
    if (samples_filename_against.size() > 0) {
        ifstream indata;
        indata.open(samples_filename_against);        
        string line;
        while (getline(indata, line)) {
            boost::algorithm::trim(line);
            required_samples.push_back(line);
            samples_names_against.push_back(line);
        }            
    }
    int n_samples_required = required_samples.size();
    if (n_samples_required > 0) {
        cout << boost::format("Only %d samples required...\n") % n_samples_required;
    }

    // Resize the result
    ret.resize(0);

    // TODO CHECK ERRORS
    auto file = hts_open(filename.c_str(), "r");
    auto header = bcf_hdr_read(file);
    auto record = bcf_init();

    int n_samples = bcf_hdr_nsamples(header);
    cout << boost::format("Reading %d samples in vcf...\n") % n_samples;

    int n_relevant_sample = 0;
    vector<bool> sample_is_required;
    for (int i = 0; i < n_samples; i++) {
        samples_names.push_back(header->samples[i]);
        if (n_samples_required > 0) {
            sample_is_required.push_back(
                std::find(required_samples.begin(), required_samples.end(), header->samples[i]) != required_samples.end()
            );
        } else {
            sample_is_required.push_back(true);
        }
        if (sample_is_required.back()) {
            samples_in_order.push_back(header->samples[i]);

            if (std::find(samples_names_against.begin(), samples_names_against.end(), header->samples[i]) != samples_names_against.end()) {
                samples_against_indices.push_back(n_relevant_sample);
            } else {
                samples_indices.push_back(n_relevant_sample);
            }
            n_relevant_sample++;
        }
    }
    if (n_samples_required == 0) {
        n_samples_required = n_samples;
    }

    int32_t *gt_arr = NULL, ngt_arr = 0;
    int ngt;

    vector<int8_t> int_alleles;

    while (true) {
        // Read record
        int retcode = bcf_read(file, header, record);
        if (retcode < 0) {
            if (retcode < -1) {
                std::cerr << "Failed to parse VCF record: " << retcode << std::endl;
                exit(-1);
            }
            break;
        }


        bcf_unpack(record, BCF_UN_ALL);

        // If not SNP, skip
        if (!bcf_is_snp(record)) {
            continue;
        }

        // Reset vector
        //int_alleles.resize(0);
        int_alleles.resize(n_samples_required * 2);
        int int_alleles_index = 0;

        // Get genotypes
        ngt = bcf_get_genotypes(header, record, &gt_arr, &ngt_arr);
        if (ngt <= 0) {
            std::cerr << "No genotypes!" << std::endl;
            return;
        }

        bool het = false;
        int max_ploidy = ngt / n_samples;
        for (int i = 0; i < n_samples; i++) {
            if (!sample_is_required[i]) {
                continue;
            }

            int32_t *ptr = gt_arr + i*max_ploidy;
            for (int j = 0; j < max_ploidy; j++)
            {
                // if true, the sample has smaller ploidy
                if (ptr[j] == bcf_int32_vector_end) break;

                // missing allele
                if (bcf_gt_is_missing(ptr[j])) {
                    //int_alleles.push_back(-1);
                    int_alleles[int_alleles_index] = -1;
                    int_alleles_index++;
                } else {
                    // the VCF 0-based allele index
                    int32_t al = bcf_gt_allele(ptr[j]);
                    //int_alleles.push_back(al);
                    int_alleles[int_alleles_index] = al;
                    int_alleles_index++;

                    if (al > 0) {
                        het = true;
                    }
                }

                // is phased?
                // int is_phased = bcf_gt_is_phased(ptr[j]);                
            }            
        }

        if (het) {
            ret.push_back(make_unique<SegregatingSite>());
            ret.back()->pos = record->pos;
            ret.back()->alleles = int_alleles;        
        }
    }

    free(gt_arr);
    bcf_destroy(record);
    bcf_hdr_destroy(header);
    hts_close(file);
}

void readMask(
    string filename,
    vector<pair<int, int>>& global_mask
) {
    fstream bedfile(filename);
    // TODO: check errors

    string bedline;
    while (getline(bedfile, bedline)) {
        vector<string> bedparts;
        boost::algorithm::split(bedparts, bedline, boost::is_any_of("\t"));
        global_mask.push_back(make_pair(stoi(bedparts[1]), stoi(bedparts[2])));  // TODO: Check errors
    }
}

void readMasks(
    string filename,
    unordered_map<string, vector<pair<int, int>>>& mask_map,
    const vector<string>& sample_names
    ) 
    {
    ifstream fin(filename);
    // TODO: check errors

    string line;
    while (getline(fin, line)) {
        // Split line into tab-separated parts
        vector<string> parts;
        boost::algorithm::split(parts, line, boost::is_any_of("\t"));

        if ((std::find(sample_names.begin(), sample_names.end(), parts[0]) != sample_names.end())) {
            mask_map.emplace(parts[0], vector<pair<int, int>>());
            readMask(parts[1], mask_map[parts[0]]);
            // cout << boost::format("Mask for %s - %d segments\n") % parts[0] % mask_map[parts[0]].size();
        }
    }
    fin.close();

}