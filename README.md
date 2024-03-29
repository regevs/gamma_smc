# Getting Started

# Installation

## Requirements

`gamma_smc` requires the following to run:
- The linux OS
- A processor that supports AVX2 instructions (run `lscpu | grep avx2` to confirm)

The following Python packages are needed to parse the raw data:
- [zstandard](https://github.com/indygreg/python-zstandard) (tested with version 0.19.0)
- [pandas](https://pandas.pydata.org/) (tested with version 1.4.1) 

## Option 1: Singularity (no compilation, preferred)
The simplest way to run `gamma_smc` is through a container, using `singularity` (or Docker, if you have root permissions), e.g.:
<pre>
$ singularity run 
    -B <i>local_dir</i>:<i>container_dir</i> 
    docker://docker.io/regevsch/gamma_smc:v0.2 
    -i <i>container_path</i> 
    -o <i>container_path</i>
    <i>other arguments</i>
</pre>

For example,
<pre>
$ singularity run 
    -B /home/regev:/mnt
    docker://docker.io/regevsch/gamma_smc:v0.2 
    -i /mnt/input.vcf.gz 
    -o /mnt/output.zst 
    -t 1
</pre>

**Note**: Using a container is easier but is somewhat slower (e.g. due to optimizing compilation to your specific processor).

## Option 2: Compilation
`gamma_smc` requires the following to be installed in your system to compile:
- [Boost](https://www.boost.org/) (tested with version 1.74)
- [htslib](https://github.com/samtools/htslib) (tested with version 1.15)
- [zstd](https://facebook.github.io/zstd/) (tested with version 1.5.2)

### Manual installation
These commands should install the prerequisites directly. Change the directories to your preferred location, where you have permissions or where it makes sense. 
```
export INSTALL_DIR=/path/to/install/dir

# Install boost. This may take some time.
wget https://boostorg.jfrog.io/artifactory/main/release/1.81.0/source/boost_1_81_0.tar.bz2 && \
    tar xjf boost_1_81_0.tar.bz2 --checkpoint=.100 && \
    cd ./boost_1_81_0 && \
    ./bootstrap.sh --with-libraries=filesystem,iostreams,system,program_options --prefix=${INSTALL_DIR} && \
    ./b2 -j 8 install && \
    cd ..
    
# Install zstd.
wget -q https://github.com/facebook/zstd/releases/download/v1.5.4/zstd-1.5.4.tar.gz && \
    tar xvf zstd-1.5.4.tar.gz && \
    cd zstd-1.5.4 && \
    make -j 8 && \
    cd ..
    
# Install htslib
wget https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 && \
    tar xjf htslib-1.17.tar.bz2 --checkpoint=.100 && \
    cd htslib-1.17 && \
    ./configure --prefix=${INSTALL_DIR} && \
    make -j 8 && \
    make install && \
    cd ..
    
# Add paths to environmental variables. Consider adding these lines to ~/.bashrc
export GAMMA_SMC_INCLUDE_PATHS=${INSTALL_DIR}/include:${INSTALL_DIR}/zstd-1.5.4/lib/
export GAMMA_SMC_LIB_PATHS=${INSTALL_DIR}/lib:${INSTALL_DIR}/zstd-1.5.4/lib/
export CPATH=$GAMMA_SMC_INCLUDE_PATHS:$CPATH
export LIBRARY_PATH=$GAMMA_SMC_LIB_PATHS:$LIBRARY_PATH
export LD_LIBRARY_PATH=$GAMMA_SMC_LIB_PATHS:$LD_LIBRARY_PATH
```

### Conda 
You could try to manage these installations using `conda`:
```
conda create --name gamma_smc --yes python=3.11
conda activate gamma_smc
conda install --yes --name gamma_smc --channel conda-forge boost-cpp==1.81 zstd==1.5.2 zstandard==0.19.0 htslib==1.18 pandas gxx==13.2
export CPATH=$CONDA_PREFIX/include:$CPATH
export LIBRARY_PATH=$CONDA_PREFIX/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
```
This creates a new environment, installs the requirements and adds their paths to the relevant environmental variables as required below.

### How to compile
Make sure that the installation paths are set, e.g. by appending `include` paths to the `CPATH` environmental variable, and `lib` paths to `LIBRARY_PATH` and `LD_LIBRARY_PATH` environmental variables. 

To compile, you need: 
- The linux OS
- A processor that supports AVX2 instructions (run `lscpu | grep avx2` to confirm)
- A modern C++ compiler (tested with gcc 8.5.0)

Download and compile with:
```
git clone git@github.com:regevs/gamma_smc.git
cd gamma_smc && make bin/gamma_smc
```

# Minimal Usage
A minimal command line is:
<pre>
$ bin/gamma_smc 
    -t <i>recombination_to_mutation_ratio</i>
    -i <i>input_file.vcf</i>
    -o <i>output_file.zst</i>
</pre>

When finishing successfully, `output_file.zst` is a binary file which can be opened as described below.

There are more options. Below we discuss each in detail.

# Full Usage
## Input
<pre>
--input, -i <i>input_file.vcf</i>
</pre>
The input file is in `vcf`, `vcf.gz`, `bcf` or `bcf.gz` formats. 

Use a separate file per chromosome.

Positions that are not specified explicitly in the `vcf` are assumed to be homozygous (unless masked out). It is therefore better and faster to process the `vcf` such that all sites are segregating (i.e. there is at least one derived allele in the sample). 

## Specifying a subset of samples
If no option is specified, all the samples in the input file will be used.
<pre>
--samples, -S <i>samples_filename.txt</i>
</pre>
If given, this is the path to a file containing a subset of sample names. This is a text file with a single sample name per line. Only these samples will be used.

## Specifying a subset of haplotype pairs
If no option is given, posteriors will be inferred for all haplotype pairs across all selected samples. **Note**: For meaningful inference between haplotypes across different samples, your data need to be phased with reasonable accuracy. If your sampled are not/poorly phased, you can still run inference within each sample (see `-w` below).

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

**Note**: Specifying a separate mask per sample may slow down inference, since now Gamma-SMC needs to calculate the intersection of each pair of masks. For many applications, a single global mask will suffice.

## Mutation rate
If the scaled mutation rate is not specified, it will be estimated from the data as the average heterzygosity across the samples, taking only masked regions into consideration, if a mask is provided. Otherwise, it can be specified explicitly. This is necessary to compare posterior inference across datasets, or across subsets of the same dataset, etc.
<pre>
--scaled_mutation_rate, -m <i>theta</i>
</pre>
Gamma-SMC requires the scaled mutation rate, defined as $\theta = 4\cdot N_e \cdot \mu$, where:
- $\mu$ is the (unscaled) mutation rate, in units of mutations per generation per base-pair
- $N_e$ is the effective population size, given as a number of diploids

Gamma-SMC is fairly robust to mis-specification of this argument, but it should be in the same ballpark as the true value. For humans (esp. of European ancestry), $\theta = 0.00075$ is a good estimate. 

## Recombination rate
The scaled recombination rate currently cannot be estimated from the data. To specify it, you can either specify it explicitly, or specify the recombination-to-mutation rate ratio.
<pre>
--scaled_recombination_rate, -r <i>rho</i>
</pre>
Gamma-SMC also requires the scaled recombination rate, defined as $\rho = 4\cdot N_e \cdot r$, where:
- $r$ is the (unscaled) recombination rate, in units of recombinations per generation per base-pair
- $N_e$ is the effective population size, given as a number of diploids

Gamma-SMC is also robust to mis-specification of this argument, but it should be in the same ballpark as the true value. For humans (esp. of European ancestry), $\rho = 0.0006$ is a good estimate. 

If you have an estimate of $\mu, r$ and also an estimate of $\theta$ (perhaps from the data), you can extract an estimate for $\rho = \theta/\mu\cdot r$.

<pre>
--recombination_to_mutation_ratio, -t <i>rho</i>
</pre>
The ratio of the recombination rate and the mutation rate. This is convenient when the scaled mutation rate is estimated from the data.

## Flow field
<pre>
--flow_field, -f <i>flow_field.txt</i>
</pre>
A path to the flow field specification. Gamma-SMC requires a flow field as input to work. If unspecified, a default one is used (hard-coded), which is also available for reference at `resources/default_flow_field.txt`. This should be reasonable for most use cases. If you want to generate a new one, you also need to install:

- [Arb](https://arblib.org/) and its prerequisites mpfr, flint (tested with version 2.22.0)
- [GSL](https://www.gnu.org/software/gsl/) (tested with version 2.7)
- [Eigen](https://eigen.tuxfamily.org/) (tested with version 3.4.0)

You may compile with
```
$ make bin/generate_canonical_flow_field
```
and run:
```
$ bin/generate_canonical_flow_field ...
```

## Cache size
<pre>
--cache_size, -z <i>cache_size</i>
</pre>
Gamma-SMC caches operating on stretches of homozygosity or missingness for fast inference. Larger cache may result in faster running times, but in slower startup time and somewhat larger memory footprint. Note that there is little advantage in increasing cache size above the typical distance between segregating sites in your sample. Defaults to 1000 basepairs.

## Output
<pre>
--output, -o <i>posteriors.zst</i>
</pre>
The output file is a binary, `zstd` compressed array of the alphas and betas inferred by the algorithm. It is accompanied by a file with metadata to parse it (json format), with the same name and additional suffix `.meta` (e.g., `posteriors.zst.meta`).

The function `open_posteriors` in `src/reader.py` is used to parse the binary file using the metadata, see below for more.

## Specifying output positions
Inferring the posterior TMRCA distribution at each and every position is too computationally expensive and not needed. Gamma-SMC therefore supports outputting in two schemes: (i) at segregating sites; (ii) in fixed jumps along the genome. Both schemes can be used concurrently.
<pre>
--output_at_hets, -h <i>or</i> --output_at_hets=false
</pre>
Output at segregating sites - these are assumed to be all the sites in your input file. This is the default.

<pre>
--output_at_stride, -s <i>jump_size</i>
</pre>
Output every <i>jump_size</i> basepairs. Turned off by default. A good value is `-s 1000` or `-s 100`; smaller values will probably be too long to run and generate huge files.

# Interpreting the output
If you haven't already, clone the repository:
```
git clone git@github.com:regevs/gamma_smc.git
```

Use `open_posteriors` from `src/reader.py` to open the output (in Python):
```Python
import sys; sys.path.append("/path/to/git/gamma_smc/src")
import reader

alphas, betas, meta = reader.open_posteriors("/path/to/output_file.zst")
```
The returned `alphas` and `betas` are pandas dataframes, and `meta` is a dictionary describing the dataset. The posterior distribution of the TMRCA at the $i$-th position, for the $j$-th pair, is described by $\Gamma(\alpha_{i,j}, \beta_{i,j})$, where $\alpha_{i,j}$ can be obtained by `alphas.iloc[i,j]`, and same for $\beta_{i,j}$. From this one can obtain, e.g. the mean posterior TMRCA, by $\alpha_{i,j}/\beta_{i,j}$. Other quantities, like the MAP (mode) or the variance can similarly be obtained. TMRCAs are defined in units of *coalescence time*, which are equivalent to units of $2N_e$ generations. This means that, if you want to convert posteriors to # of generations, you would need to estimate $2N_e$ from the scaled mutation rate (estimated from the data) and the unscaled mutation rate, which must be obtained from another source.

In `meta`, useful properties are:
- `scaled_mutation_rate` and `scaled_recombination_rate`
- `output_positions` - a list of genomic positions at which TMRCA posteriors were inferred
- `pairs` - the order of pairs for which posteriors were inferred, using serial numbers (e.g. `0_1`)
- `sample_names` - a mapping from a serial number a the sample name; e.g. `sample_names = {0: "sample_name.0", 1: "sample_name.1", 2: "another_sample_name.0", ...}` 


