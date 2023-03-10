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

Make sure that the installation path are set, e.g. by appending `include` paths to the `CPATH` environmental variables, and `lib` paths to `LIBRARY_PATH` and `LD_LIBRARY_PATH`.

Download and compile with:
```
git clone https://github.com/regevs/gamma_smc
cd gamma_smc && make
```

