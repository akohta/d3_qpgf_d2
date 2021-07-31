# d3_qpgf_d2
This is the calculation program of quasi-periodic Green's function for the Helmholtz equations. 
The quasi-periodicity is 2-dimensions ( z component is zero ), Green's function is 3-dimensions. 
This program is used my original method.

## Definitions
- quasi-periodic Green's function  
  <img src="https://latex.codecogs.com/gif.latex?\,^q\!G(\mathbf{r})=\sum_{l_1=-\infty}^{\infty}\sum_{l_2=-\infty}^{\infty}G(\mathbf{r}+l_1\mathbf{d}_1+l_2\mathbf{d}_2)\exp\left(i\mathbf{k}\cdot(l_1\mathbf{d}_1+l_2\mathbf{d}_2)\right)">  
  <img src="https://latex.codecogs.com/gif.latex?G(\mathbf{r})=\frac{\exp(ik|\mathbf{r}|)}{4\pi|\mathbf{r}|}">  
  <img src="https://latex.codecogs.com/gif.latex?\,^q\!G(\mathbf{r})"> is quasi-periodic Green's function  
  <img src="https://latex.codecogs.com/gif.latex?G(\mathbf{r})"> is Green's function of the 3-dimensional Helmholtz equation  
  <img src="https://latex.codecogs.com/gif.latex?\mathbf{k}"> is wave number vector,
  <img src="https://latex.codecogs.com/gif.latex?|\mathbf{k}|=k">  
  <img src="https://latex.codecogs.com/gif.latex?\mathbf{d}_1"> is 1st lattice vector,
  <img src="https://latex.codecogs.com/gif.latex?\mathbf{d}_1=(d_{1x},d_{1y},0)">  
  <img src="https://latex.codecogs.com/gif.latex?\mathbf{d}_2"> is 2nd lattice vector,
  <img src="https://latex.codecogs.com/gif.latex?\mathbf{d}_2=(d_{2x},d_{2y},0)">  
  These lattice vectors must be linearly independent.

- quasi-periodic condition  
  
  <img src="https://latex.codecogs.com/gif.latex?\,^q\!G(\mathbf{r}+l_1\mathbf{d}_1+l_2\mathbf{d}_2)=\exp\left(-i\mathbf{k}\cdot(l_1\mathbf{d}_1+l_2\mathbf{d}_2)\right)\,^q\!G(\mathbf{r}),l_1\in\mathbb{Z},l_2\in\mathbb{Z}">
  
## Usage of example code
1. type 'make' command to compile
2. type './example.out' to run  

Please see src/d3_qpgf_d2.h for detail of functions, src/example.c for detail of function usages.

## References
1. Capolino, Filippo, Donald R. Wilton, and William A. Johnson. "Efficient computation of the 2-D Green's function for 1-D periodic structures using the Ewald method." IEEE Transactions on Antennas and Propagation 53.9 (2005): 2977-2984.  
27
2. Beylkin, Gregory, Christopher Kurcz, and Lucas Monzón. "Fast algorithms for Helmholtz Green's functions." Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences 464.2100 (2008): 3301-3326.
28
3. Abramowitz, Milton, Irene A. Stegun, and Robert H. Romer. "Handbook of mathematical functions with formulas, graphs, and mathematical tables." (1988): 958-958.
29
4. Faddeeva Package. http://ab-initio.mit.edu/Faddeeva
