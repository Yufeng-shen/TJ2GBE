# TJ2GBE
Reconstruct grain boundary energy from triple junction geometries [1,2].

This package is an implementation of the regularization based reconstruction method which is described in paper [1]. It consists of two parts: 
  *  C++ code that finds the "similar grain boundaries" of every grain boundary in the data set
  *  Python script that solves the optimization problem described in paper [1] to reconstruct the energy of every grain boundary in the data set

## Financial Support
The development of this package was supported by the National Science Foundation of the United States of America under grant DMR-1628994.

## Dependencies

- numpy

- scipy

- matplotlib

## Usage

Step 1: Compile the C++ code.
```shell
cd Src/
g++ -o ../bin/TJ2GBE.out main.cpp tj.cpp subdomain.cpp config.cpp -fopenmp
```

Step 2: Execute the compile file.
```shell
cd ../bin/
./TJ2GBE.out ../CfgFile/Cubic.config
```
The output files "rowA.binary", "colA.binary" and "valA.binary" are going to be used by Pythoon reconstruction script, other files are for debugging and tunning the hyperparameters in the configure file.

Step 3: Run the reconstruction script.
```shell
cd ../Src/
python Reconstruction.py
```
It will output the reconstructed energy for every grain boundary in the data set and a .gbdat file for plotting.

For the example data of "TJdata/triples\_30000.dat", the reconstructed grain boundary energy function with &Sigma; 7 misorientation is shown as following:
![Sigma7](https://github.com/Yufeng-shen/TJ2GBE/blob/master/misc/Sigma7.png)

## Notes

There are two differences between the cell indexing in Morawiec's Fortran implementation of the paper [2] and this implementation:

- In this implementation, cell index starts with 0 instead of 1.

- In function "find\_symeq\_num", the second and fourth parameters are defined as "1-cos(af)" instead of "cos(af)".

## License
__BSD 3 Cluase License__

## References
[1]  ["Determining grain boundary energies from triple junction geometries without discretizing the five-parameter space"](https://doi.org/10.1016/j.actamat.2018.12.022) Acta Materialia 166, 126-134

[2]  ["Method to calculate the grain boundary energy distribution over the space of macroscopic boundary parameters from the geometry of triple junctions"](https://doi.org/10.1016/S1359-6454(00)00126-9) Acta Materialia, 48 (2000) 3525-3532
