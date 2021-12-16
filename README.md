# Code-Binodal-Surfaces
This repository contains functions in OSCAR that accompany the paper "Towards tropically counting binodal surfaces" by Madeline Brandt and Alheydis Geiger, 12/2021.
It was developed for OSCAR Version 0.4.0

## How to use
After installing Oscar in Julia (see https://oscar.computeralgebra.de/install/), start Julia and run
<code> using Oscar</code>
<code>\include("PATH/function_countbinodalsurfaces.jl")</code>

### Examples
Examples are provided in the jupyter notebook "OSCAR investigation of binoda polytopes.ipynb"

### Short description of the functions
The function <code>binodal()</code> computes the binodal variety of polynomials with a given support.

The function <code>investigate_binodalvariety()</code> prints the following data about the binodal variety:

- dimension, 
- degree,
- whether the binodal variety contains monomials.

The function <code>con_mult()</code> computes the multiplicity of a specified connected lattice path for a given polynomial support.

A function to compute the multiplicity of a specified disconnected lattice path for a given polynomial support will be added soon.

### Disclaimer
The computations do not terminate for larger values on a standard computer.  
