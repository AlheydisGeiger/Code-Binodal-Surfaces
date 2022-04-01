# Code-Binodal-Surfaces
This repository contains functions in OSCAR that accompany the paper "Towards tropically counting binodal surfaces" by Madeline Brandt and Alheydis Geiger, 12/2021.
It was developed for OSCAR Version 0.8.2

## How to use
After installing Oscar in Julia (see https://oscar.computeralgebra.de/install/), start Julia and run
<code> using Oscar</code>
<code>\include("PATH/code_binodals.jl")</code>

### Examples
Examples are provided in the jupyter notebook "OSCAR investigation of binoda polytopes.ipynb".

### Short description of the functions
The function <code>binodal()</code> computes the binodal variety of polynomials with a given support.

The function <code>investigate_binodalvariety()</code> prints the following data about the binodal variety:

- dimension, 
- degree,
- whether the polytope can be binodal.

The function <code>path_mult()</code> computes the multiplicity of a specified lattice path with at most one gap for a given polynomial support.


### Disclaimer
The computations do not terminate for larger values on a standard computer.  
