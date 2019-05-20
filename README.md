# Functional analysis of genetic mutations

The aim of this script is to produce a table from a VCF for a given gene. The table contains gene and CDS positions of the mutations, protein mutations, and allele and sample numbers.

## Prerequisites

To run the script properly, you need a VCF, a GFF of the gene of interest and a population file (optional).

## Installation

To download the latest version of the files:
```
git clone https://github.com/fdchevalier/functional_analysis
```

For convenience, the script should be accessible system-wide by either including the folder in your `$PATH` or by moving the script in a folder present in your path (e.g. `$HOME/local/bin/`). The email template must be in the same directory as the script.

## Usage

A summary of available options can be obtained using `./func-table-mut.sh -h`.

## License

This project is licensed under the [GPLv3](LICENSE).
