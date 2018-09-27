#### GEGS - GPU implementation of EGSnrc
This is a project based on GPU-parallelized implementation of EGSnrc.  
It is used to speed up the calculation of dosimestry in Radiation Treatment.  
Used Monte-Carlo method to simulate. Used CUDA C for coding language.

We are using two nVidia Tesla C2050 for computation, and later one K20c.
We develop this module on Windows 7 system and using Visual Studio 2010.

EGSnrc is a Monte Carlo package for simulating Electron and Photon transportation.
This package is running on CPU and using MORTRAN (an extention of Fortran).
We first tranplant it into C then implement it into CUDA C for running on GPU.
For some cases, we gained more than 60 times speedup compare to single core of CPU.
