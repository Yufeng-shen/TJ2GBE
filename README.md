# TJ2GBE
Reconstruct grain boundary energy from triple junction geometries [1].


## Dependencies

- numpy

- scipy

- matplotlib

## Usage
This code consists of two parts, the C++ code which finds the "neighboring grain boundaries" of every grain boundary in the data set, and the Python script solves the optimization problem described in paper [1].

```shell
cd Src/
g++ -o ../bin/TJ2GBE.out main.cpp tj.cpp subdomain.cpp config.cpp -fopenmp
cd bin/
./TJ2GBE.out ../CfgFile/Cubic.config
cd ../Src/
python Reconstruction.py
```

## Notes

There are two differences between the cell indexing in Morawiec's Fortran implementation of the paper [2] and this implementation:

- In this implementation, cell index starts with 0 instead of 1.

- In function "find\_symeq\_num", the second and fourth parameters are defined as "1-cos(af)" instead of "cos(af)".

## License
__BSD 3 Cluase License__

## References
[1]  ["Determining grain boundary energies from triple junction geometries without discretizing the five-parameter space"](https://doi.org/10.1016/j.actamat.2018.12.022).
[2]  ["Method to calculate the grain boundary energy distribution over the space of macroscopic boundary parameters from the geometry of triple junctions"](https://doi.org/10.1016/S1359-6454(00)00126-9)
