# Getting Started

Both functionality and documentation are in progress, and will be updated prior to publication.

# Installation

## Requirements

`gamma_smc` requires the following to be installed in your system:

- [Boost](https://www.boost.org/) (tested with version 1.74)
- [htslib](https://github.com/samtools/htslib) (tested with version 1.15)
- [pandas](https://pandas.pydata.org/) (tested with version 1.4.1)

To generate a flow field yourself (see below), you also need:

- [Arb](https://arblib.org/) and its prerequisites mpfr, flint (tested with version 2.22.0)
- [GSL](https://www.gnu.org/software/gsl/) (tested with version 2.7)
- [Eigen](https://eigen.tuxfamily.org/) (tested with version 3.4.0)


## Compilation

Make sure that the installation paths are set, e.g. by appending `include` paths to the `CPATH` environmental variable, and `lib` paths to `LIBRARY_PATH` and `LD_LIBRARY_PATH` environmental variables.

Download and compile with:
```
git clone https://github.com/regevs/gamma_smc
cd gamma_smc && make bin/gamma_smc
```

# Minimal Usage

A minimal command line is:
<pre>
bin/gamma_smc 
    -m <i>scaled_mutation_rate</i> 
    -r <i>scaled_recombination_rate</i>
    -i <i>input_file.vcf</i>
    -o <i>output_file.gz</i>
</pre>

There are more options. Below we discuss each in detail.

# Full Usage
## Input
<pre>
--input, -i <i>input_file.vcf</i>
</pre>
The input file is in a `vcf` format. It can also be a `vcf.gz` or `bcf` file. 

## Specifying a subset of samples
If no option is specified, all the samples in the input file will be used.
<pre>
--samples, -S <i>samples_filename.txt</i>
</pre>
If given, this is the path to a file containing a subset of sample names. This is a text file with a single sample name per line. Only these samples will be used.

## Specifying a subset of haplotype pairs
If no option is given, posteriors will be inferred for all haplotype pairs across all selected samples.

<pre>
--only_within, -w
</pre>
If specified, only haplotype pairs within each diploid sample will be used, but not across samples.
<pre>
--samples_against, -T <i>samples_filename_against.txt</i>
</pre>
An additional file of sample names, in the same format of `--samples`. If supplied, inference will be done only for haplotypes pairs between the first list and the second list.

## Masking
Many regions of the genome are difficult to align to, may have low coverage or have other issues. It is advised to specify a mask - a list of regions where the genome calling is reliable, and on which Gamma-SMC will use the sequence as observations. If no mask is provided, then it is assumed the entire chromosome should be used.
<pre>
--mask, -a <i>mask.bed</i>
</pre>
Provides a global mask file that applies to all samples, e.g. a species-wide mappability map. Provided in `bed` format. This is a *positive* mask - regions specified are those who should be *included*.
<pre>
--masks_per_sample, -b <i>masks_filenames.txt</i>
</pre>
You may want to have a different mask per sample, for example, if samples vary significantly by their calling coverage or quality along the genome. In this case you can instead supply a tab-separated text file, each line containing first the sample name (as specified in the VCF) and second a path to a bed file of the corresponding mask.

**Note**: Specifying a separate mask per sample may slow inference, since now Gamma-SMC needs to calculate the intersection of each pair of masks. In many application a single global mask will suffice.

## Mutation rate
<pre>
--scaled_mutation_rate, -m <i>u</i>
</pre>
Gamma-SMC requires the scaled mutation rate, defined as $u = 4\cdot N_e \cdot \mu$, where:
- $\mu$ is the (unscaled) mutation rate, in units of mutations per generation per base-pair
- $N_e$ is the effective population size, given as the number of haploids

Gamma-SMC is fairly robust to mis-specification of this argument, but it should be in the same ballpark as the true value. For humans (esp. of European ancestry), $u = 0.000375$ is a good estimate. You can also estimate this using Watterson's estimator.

## Recombination rate
<pre>
--scaled_recombination_rate, -r <i>r</i>
</pre>
Gamma-SMC also requires the scaled recombination rate, defined as $r = 4\cdot N_e \cdot \rho$, where:
- $\rho$ is the (unscaled) recombination rate, in units of recombinations per generation per base-pair
- $N_e$ is the effective population size, given as the number of haploids

Gamma-SMC is also robust to mis-specification of this argument, but it should be in the same ballpark as the true value. For humans (esp. of European ancestry), $r = 0.0003$ is a good estimate. 

If you have an estimate of $\mu, \rho$ and an estimate of $u$ (perhaps from the data), you can extract an estimate of $N_e$ from it and get an estimate for $r$.

## Flow field
<pre>
--flow_field, -f <i>flow_field.txt</i>
</pre>
A path to the flow field specification. Gamma-SMC requires a flow field as input to work. We supply a default at `resources/default_flow_field.txt` which should be reasonable for most use cases. If you just want to use the default flow field, you can omit this flag which should default to the correct filename. If you want to generate a new one, you may use
```
make bin/generate_canonical_flow_field
```
and run:
```
bin/generate_canonical_flow_field ...
```

## Cache size
<pre>
--cache_size, -z <i>cache_size</i>
</pre>
Gamma-SMC caches operating on stretches of homozygosity for fast inference. Larger cache may result in faster running times but slower startup and somewhat larger memory footprint. Note that there is little advantage in increasing cache size above the typical distance between segregating sites in your sample. Defaults to 1000.

## Output
<pre>
--output, -o <i>posteriors.gz</i>
</pre>
The output file is a binary, `gzip` compressed array of the alphas and betas inferred by the algorithm. It is accompanied by a file with metadata to parse it (json format), with the same name and additional suffix `.meta` (e.g., `posteriors.gz.meta`).

The function `open_posteriors` in `src/reader.py` is used to parse the binary file using the metadata, and return a pandas dataframe.

Inside the `.meta` file there is a dictionary translating from the haplotype numbers used in the dataframe to sample names, e.g. `sample_names = {0: "sample_name.0", 1: "sample_name.1", 2: "another_sample_name.0", ...}` etc.

## Specifying output positions
Inferring the posterior TMRCA distribution at each and every position is too computationally expensive and not needed. Gamma-SMC therefore supports outputting in two schemes: (i) at segregating sites; (ii) in fixed jumps along the genome. Both schemes can be used concurrently.
<pre>
--output_at_hets, -h
</pre>
Output at segregating sites - these are assumed to be all the sites in your input file. This is the default.

<pre>
--output_at_stride, -s <i>jump_size</i>
</pre>
Output every <i>jump_size</i> basepairs. Turned off by default. A good value is `-s 1000` or `-s 100`; smaller values will probably be too long to run and generate huge files.

# Interpreting the output
Expand here.