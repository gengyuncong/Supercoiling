# Overview

This repository contains the codes that are needed to reproduce our results in “**A spatially resolved stochastic model reveals the role of supercoiling in transcription regulation**”. 

It can be divided into three parts:

1) the source codes for the simulation tool Lattice Microbes 

2) the codes for generating the problem-specific models (i.e., the input file for the simulation)

3) the codes for analyzing the simulation results

The codes for part 2) and 3) are stored in the folder “/Codes4Simulation”. Other parts of the folder are all source codes for Lattice Microbes. 

To reproduce our results, the following steps need to be done:

1) install the prerequisite packages/softwares for Lattice Microbes and compile Lattice Microbes 

2) generate the input model in SBML format

3) convet the input file from SBML to binary format, and run the simulation by Lattice Microbes

4) analyze the simulation results in python

We will introduce each step in detail. 

# 1. Compiling Lattice Microbes

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
- g++: [https://gcc.gnu.org/](https://gcc.gnu.org/)
- Cmake: [https://cmake.org/install/](https://cmake.org/install/)
- libSBML-5.18.0: [https://sourceforge.net/projects/sbml/files/libsbml/5.18.0/stable/](https://sourceforge.net/projects/sbml/files/libsbml/5.18.0/stable/)
- hdf5
- snappy
- libxml2
- openssl
- protobuf
- ffmpeg

Most of the packages should be readily installed via brew/get-apt without specific notation. 

When all the dependencies are ready, run the following commands:

`sudo mkdir build/`

`cd build/`

`sudo cmake ..`

`sudo make -j 8`

After the compilation is done, the following executables should be found under build/utils/c:

- **build/utils/c/lm_sbml_import**: the tool to convert SBML model file to binary file (which is understandable by Lattice Microbes)

- **build/utils/c/lm_setp**: the tool to specify the simulation time and the time interval to record the trajectories

- **build/lmes**: the main program to perform stochastic simulation

It will be good to create a symbolic link of those executables under the working directory. 

# 2. Generate problem-specific models

In the manuscript “**A spatially resolved stochastic model reveals the role of supercoiling in transcription regulation**”, we built a master equation-based model to investigate the interplay between transcription and supercoiling. We mainly focused on three phenomena:

1) the cooperation of RNAP molecules at the single-gene transcription level

2) the transcriptional bursting at the topological domain level

3) the interaction between the two genes at the multi-gene transcription level

For each phenomenon, we build a specific model, described by a specific set of master equations. The model for each system is stored in a separate subfolder under “/Codes4Simulation”. 

# 3. Running the simulation

Given a specific model, the following step needs to be done:

 1) convert the model file in smbl format (say, "myfile.sbml") to binary format (say, "myfile.lm"):

`./lm_sbml_import myfile.lm myfile.sbml`
 
 2) specify the time to simulate (for example, 2000) and sampling interval (for example, 1), both in sec:

`./lm_setp myfile.lm maxTime=2000 writeInterval=$interval`
 
 3) specify the number of replicates and run:

`./lmes -r 1-1000 -f myfile.lm `


# 4. Analyzing results

The result analysis is performed in python. The "h5py" module is needed to read the binary output file from Lattice Microbe. The codes for result analysis is simply deposited in the corresponding model 
