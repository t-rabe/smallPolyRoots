# smallPolyRoots
### Purpose
Create a scatter plot on the complex plane which represents the roots to the characteristic polynomials of all possible matrices of a certain form.
For example, there are 33,554,432 matrices of size 5x5 that are filled with values -1 and 1. Each of these matrices has a characteristic polynomial with degree 5.
This means that there are 167,772,160 total roots (counting multiplicity). It is possible to display these roots on the complex plane. Interestingly enough,
**the resulting scatter plot is structured** (see below for example). 

![Plot of roots](https://github.com/t-rabe/smallPolyRoots/blob/main/images/plot4.png?raw=true)

### Workflow
This repository is designed to do the following in `C++`:
* Create random matrices of size NxN.
  * The default number of total elements created is 17,000,000.
  * Default element values are -1 and 1.
  * Matrix size is taken as an input parameter when running the program.
* Calculate the characteristic polynomials for each of the above matrices.
  * The coefficients are written to a .csv file to be used later.
  * Creation of the matrices and coeffiecients are both done in `driver.cpp`.
* Solve the roots of each characteristic polynomial.
  * These roots are then written to various output files as an MxM matrix where M represents the number of pixel to a side of the final image.
  * See below for a more detailed description of this portion of the program.

This repository is designed to do the following in `Python`:
* Read all output files and condense them into a single MxM array.
* Use `Matplotlib` and `Numpy` to plot the resulting array.
* Adjust color values to appropriately visualize all artifacts of the plot.
* Save the image.

### Root solving method
