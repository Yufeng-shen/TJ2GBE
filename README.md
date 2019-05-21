# TJ2GBE
Reconstruct grain boundary energy from triple junction geometries

This is an implementation of the paper [Determining grain boundary energies from triple junction geometries without discretizing the five-parameter space](https://doi.org/10.1016/j.actamat.2018.12.022).

## Dependencies

scipy
matplotlib

## Installation
g++ -o ../build/TJ2GBE.out main.cpp tj.cpp subdomain.cpp config.cpp -fopenmp

## Usage

## Notes

Comparing to Morawiec's code, when we indexing the cells, we start with 0 instead of 1. And we use af=1-cos(af) instead of af=cos(af) in function "find\_symeq\_num"

## License
__BSD 3 Cluase License__
