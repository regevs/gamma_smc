# Getting Started

Both functionality and documentation are in progress, and will be updated prior to publication.

# Installation

## Requirements

`gamma_smc` requires the following to be installed in your system:

- [Boost](https://www.boost.org/) (tested with version 1.74)
- [Arb](https://arblib.org/) and its prerequisites mpfr, flint (tested with version 2.22.0)
- [htslib](https://github.com/samtools/htslib) (tested with version 1.15)
- [GSL](https://www.gnu.org/software/gsl/) (tested with version 2.7)
- [Eigen](https://eigen.tuxfamily.org/) (tested with version 3.4.0)
- [pandas](https://pandas.pydata.org/) (tested with version 1.4.1)

## Compilation

Make sure that the installation paths are set, e.g. by appending `include` paths to the `CPATH` environmental variable, and `lib` paths to `LIBRARY_PATH` and `LD_LIBRARY_PATH` environmental variables.

Download and compile with:
```
git clone https://github.com/regevs/gamma_smc
cd gamma_smc && make
```

# Usage

A minimal command line is:
<pre>
gamma_smc 
    -m <i>scaled_mutation_rate</i> 
    -r <i>scaled_recombination_rate</i>
    -f <i>flow_field_file</i>
    -i <i>input_file.vcf</i>
    -o <i>output_file.gz</i>
</pre>

There are more options. Below we discuss each in detail.

## Input
```
--input_file, -i
```
The input file is in a `vcf` format.

## Output
```
--output_file, -o
```
The output file is a binary, `gzip` compressed array of the alphas and betas inferred by the algorithm. It is accompanied by a file with metadata to parse it, with the same name and additional suffix `.meta` (json format).

The function `open_posteriors` in `src/reader.py` is used to parse the binary file using the metadata, and return a pandas dataframe.

Inside the `.meta` file there is a dictionary translating from the haplotype numbers used in the dataframe to sample names, e.g. `sample_names = {0: "sample_name.0", 1: "sample_name.1", 2: "another_sample_name.0", ...}` etc.