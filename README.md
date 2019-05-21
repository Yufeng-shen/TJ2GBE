# TJ2GBE
Reconstruct grain boundary energy from triple junction geometries

This is an implementation of the paper ["Determining grain boundary energies from triple junction geometries without discretizing the five-parameter space"](https://doi.org/10.1016/j.actamat.2018.12.022).

## Dependencies

-numpy

-scipy

-matplotlib

## Usage

```shell
cd Src/
g++ -o ../bin/TJ2GBE.out main.cpp tj.cpp subdomain.cpp config.cpp -fopenmp
cd bin/
./TJ2GBE.out ../CfgFile/Cubic.config
cd ../Src/
python Reconstruction.py
```

## Notes

There are two differences between the cell indexing in Morawiec's Fortran implementation and this implementation:

- In this implementation, cell index starts with 0 instead of 1.

- In function "find\_symeq\_num", the second and fourth parameters are defined as "1-cos(af)" instead of "cos(af)".

## License
__BSD 3 Cluase License__
