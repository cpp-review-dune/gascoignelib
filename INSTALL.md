# Installation of Gascoigne

This is just a quick guide with some minimum instructions. For details and options we refer to the webpage

[www.gascoigne.de](www.gascoigne.de)

## Requirements

For installation, we require:

1) A C++ compiler, e.g. g++ or clang++ supporting at least C++-11.


2) cmake, see
[http://www.cmake.org/HTML/Index.html](http://www.cmake.org/HTML/Index.html)
(used as Makefile generator)


3) umfpack (part of suitesparse)
[https://people.engr.tamu.edu/davis/suitesparse.html](https://people.engr.tamu.edu/davis/suitesparse.html)
(used as direct sparse solver)


4) metis
[http://www-users.cs.umn.edu/~karypis/metis/](http://www-users.cs.umn.edu/~karypis/metis/)
(used for optimisation of matrix storage)


All these packages are available for all Linux distributions and als part of the homebrew ([https://brew.sh/index_de](http://www.cise.ufl.edu/research/sparse/umfpack/)) project on a Mac. 

We have no experience with Windows systems but we would welcome a guide.

## Configuration & Compilation

We assume that the Gascoigne Library is found in a directory called

> GAS/GascoigneLib

This is the directory where you should find this file

> GAS/GascoigneLib/INSTALL.md

For compilation create a build dir, e.g. 

> GAS/build

Then, to configure, just call

> cd GAS/build
> 
> cmake GAS/GascoigneLib

And finally, for compilation, in the same directory type

> make 

You are done. The Library GascoigneLib will be placed in 

> GAS/lib

## Options

There are many options to be set. Just have a look at the webpage [www.gascoigne.de](www.gascoigne.de) and try

> cd GAS/build
> 
> ccmake ../GascoigneLib




