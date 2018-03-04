# Small-Data, Large-Scale Linear Optimization

In the spirit of reproducible research, this repository contains all the code necessary for running the experiments and creating graphs in the paper:  
> "Gupta, Vishal and Rusmevichientong, Paat, Small-Data, Large-Scale Linear Optimization (Oct 31 1, 2017). Available at SSRN: https://ssrn.com/abstract=3065655."

The full-text of the paper is available on SSRN or the author's [website](http://www-bcf.usc.edu/~guptavis/research.html).

If you find this code or the paper useful, ***please consider citing it***.

## Overview
All of the source code for computing solutions by various methods can be found in KPNormal.jl, written in Julia.  

The files:
 - TestHarness_Paper.jl
 - TestCLTHarness.jl

call workhorse functions from KPNormal.jl in a simulation set-up to generate the data presented in the paper.  

Since all of the algorithms are presented are single-threaded, they are "embarrassingly parallel."  Substantive speed-ups can be achieved in multi-threading.  The files

 - test3Part.jl
 - testPortExp.jl
 - testPOAPCLT.jl
 - testLoo.jl
 
are wrappers that call the above tests in a 4-threaded environment and compile results.  Experiments for the paper were run in this way for efficiency.

Finally the folder **plotting** contains functions used to generate plots for the paper and **ArchiveExperiments** contain experiments and files not used in the final draft.  

## Licensing

This code is available under the MIT License.  
Copyright (c) 2017 Vishal Gupta
