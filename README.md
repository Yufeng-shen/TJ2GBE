# TJ2GBE
Reconstruct grain boundary energy from triple junction geometries [1,2].

This package is an implementation of the regularization based reconstruction method, which is described in paper [1]. Many functions are borrowed from Adam Morawiec's Fortran implementation of his paper [2]. Moreover, the file "Src/myLOBPCG\_new.py" is [scipy's implementation](https://github.com/scipy/scipy/blob/v1.3.0/scipy/sparse/linalg/eigen/lobpcg/lobpcg.py) with minor modification. 

This package consists of two parts: 
  *  C++ code that finds the "similar grain boundaries" of every grain boundary in the data set
  *  Python script that solves the optimization problem described in paper [1] to reconstruct the energy of every grain boundary in the data set

## Financial Support
The development of this package was supported by the National Science Foundation of the United States of America under grant DMR-1628994.

## Dependencies

- numpy
```
pip install numpy
```

- scipy
```
pip install scipy
```

## Usage

Step 0: Download the code.
```shell
git clone --depth=1 https://github.com/Yufeng-shen/TJ2GBE.git
```

Step 1: Compile the C++ code.
(Thanks to @sgbaird, we find that C++11 is required)
```shell
cd Src/Cpp/
g++ -std=c++11 -o ../../bin/TJ2GBE.out main.cpp tj.cpp subdomain.cpp config.cpp -fopenmp
```

Step 2: Run the executable file with configure file.
```shell
cd ../../bin/
./TJ2GBE.out ../CfgFile/Cubic.config
```
The output files "rowA.binary", "colA.binary" and "valA.binary" are going to be used by Pythoon reconstruction script, the "NN.txt" is used for tunning the "threshold" value in the configure file.

Step 3: Change the parameters (at the top of the script) and run the reconstruction script.
```shell
cd ../Src/Python/
python Reconstruction.py
```
It will output the reconstructed energy for every grain boundary in the data set and a .gbdat file for plotting.

For the example data of "TJdata/triples\_30000.dat", the reconstructed grain boundary energy function with &Sigma;7 misorientation is shown as following:
![Sigma7](https://github.com/Yufeng-shen/TJ2GBE/blob/master/misc/Sigma7.png)
Figures are plotted by a modified version of Krzysztof Glowinski's [GBToolbox](http://imim.pl/personal/adam.morawiec/A_Morawiec_Web_Page/S/K_Glowinski/Downloads.html).

## Notes

There are two differences between the cell indexing in Morawiec's Fortran implementation of the paper [2] and this implementation:

- In this implementation, cell index starts with 0 instead of 1.

- In function "find\_symeq\_num", the second and fourth parameters are defined as "1-cos(af)" instead of "cos(af)".

## License
__BSD 3 Cluase License__

## References
[1]  ["Determining grain boundary energies from triple junction geometries without discretizing the five-parameter space"](https://doi.org/10.1016/j.actamat.2018.12.022) Acta Materialia 166, 126-134 (2019).

[2]  ["Method to calculate the grain boundary energy distribution over the space of macroscopic boundary parameters from the geometry of triple junctions"](https://doi.org/10.1016/S1359-6454(00)00126-9) Acta Materialia, 48 (2000) 3525-3532
