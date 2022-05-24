# Overview

This repository contains the codes that are needed to reproduce our results in “**A spatially resolved stochastic model reveals the role of supercoiling in transcription regulation**”. 

It can be divided into three parts:

1) the source codes for the simulation tool Lattice Microbes 

2) the codes for generating the problem-specific models (i.e., the input file for the simulation)

3) the codes for analyzing the simulation results

The codes for part 2) and 3) is stored in the subfolder “/Codes4Simulation”. Other parts of the folder are all source codes for Lattice Microbes. 

To reproduce our results, the following steps need to be done:

1) install the prerequisite packages/softwares for Lattice Microbes

2) compile Lattice Microbes 

3) generate the input file for simulation

4) run the simulation by Lattice Microbes 

5) analyze the simulation results 

# Installing Lattice Microbes

Lattice Microbes is a software package to efficiently sample trajectories from chemical master equations and reaction diffusion master equations using Gillespie algorithm. The input is a binary file that defines species (with their initial numbers) and the reactions between them (with the propensity function and the associated kinetic rates). The output is a binary file that samples the number of species at a given time interval. 

To use Lattice Microbes, several external packages should be installed:

- Python 2.7, with python development tools and the following packages
    - numpy
    - scipy
    - matplotlib
    - protobuf
    - scikit-image
    - h5py
- Git: [https://github.com/git-guides/install-git](https://github.com/git-guides/install-git)
- g++
- Cmake: [https://cmake.org/install/](https://cmake.org/install/)
- libSBML-5.18.0: [https://sourceforge.net/projects/sbml/files/libsbml/5.18.0/stable/](https://sourceforge.net/projects/sbml/files/libsbml/5.18.0/stable/)
- hdf5:
- snappy
- libxml2
- openssl
- protobuf
- ffmpeg
- Protocol Buffers v2.4.1
- lxml 3.2.4+

Most of the packages could be installed via brew/get-apt. 

When all the dependencies are ready, run the following commands:

`sudo mkdir build/`

`cd build/`

`sudo cmake ..`

`sudo make -j 8`

# Generate problem-specific models

In the manuscript “**A spatially resolved stochastic model reveals the role of supercoiling in transcription regulation**”, we built a master equation-based model to investigate the interplay between transcription and supercoiling. We mainly focused on three phenomena:

1) the cooperation of RNAP molecules at the single-gene transcription level

2) the transcriptional bursting at the topological domain level

3) the interaction between the two genes at the multi-gene transcription level

For each phenomenon, we simulated a specific system, described by a specific set of master equations. The codes that could be used to reproduce our work is stored in the “/Codes4Simulation” subfolder. 
