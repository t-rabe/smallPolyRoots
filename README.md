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

### Root finding method
Previous algorithms for root solving typically use variations of [Newton's Method](https://web.ma.utexas.edu/users/m408n/CurrentWeb/LM4-8-2.php#:~:text=Newton's%20method%20is%20a%20technique,). There have been [recent advances](https://doi.org/10.1016/j.camwa.2010.12.070) in these methods, though they still have issues with precision for large polynomials and computational time for large quantities of polynomials.

The method used to solve roots in this codebase takes a different approach. It is built upon the idea that [Horner's Method](https://www3.nd.edu/~zxu2/acms40390F13/Lec-2.6.pdf) for polynomial evaluation is extremely efficient. Taking into account that images can only display finite precision (pixel sizes can be arbitrarily small, but they are always finite in number), it stands to reason that it is possible to evaluate each characteristic polynomial at each pixel in order to find roots. 

NOTE: I do not claim to have a better performing algorithm for root finding. This repository is simply very efficient for one type of root finding and visualization.

#### Steps
1. Construct an MxM 2D histogram filled with zeros. Each bin represents a point on the complex plane (imagine a grid on the plane).
2. Start with one polynomial and construct an MxM grid on the complex plane filled with zeros for this particular polynomial.
3. Evaluate the polynomial at every square on the grid and set the value of that square equal to the value of the polynomial at that point.
4. Search for local minima.
5. Every time a local minimum is found, increase the corresponding grid on the histogram by 1.
6. Repeat steps 2-5 for every other polynomial which was previously defined by `driver.cpp`.

There are a couple of interesting facts to note here.
1. Each MxM grid for individual polynomials will have the correct number of local minimum (excluding multiplicity). This means that all roots are found for a given polynomial, assuming that the grid is large enough to encompass all roots (no roots are outside the boundaries of the grid) and that M is sufficiently large (there are enough pixels in the image for roots to lie on their own pixel).
2. This is FASTER than other methods in some cases; large ish polynomials (degree >=5) with low ish resolution images (4K images), and comparable in most other cases.
3. The results can be tested using previous methods and are found to be correct!
