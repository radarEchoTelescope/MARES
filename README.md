This branch will implement some new output types for MARES. We aim to have an output file with peak voltage and the corresponding peak time; a file with positions of the segments; phase coherence measure and possibly more.  
This will involve implementing a Single output file per event, following the format of a NuRadio event file: top level properties, MonteCarlo-truth simulated signal and/or the complete voltage at the receiver with noise, raytracing, multiple showers etc. 

We also aim to add signal frequency spectra as an output option directly calculated by MARES, to simplify (and standardise) analyses. This can make use of the FFTW package (https://fftw.org/), following the corresponding tutorial here (https://fftw.org/fftw3_doc/Tutorial.html#Tutorial). 

Extra: 
This event file will also include functionality to develop time evolution for flavour studies. The idea is that the stochastic energy losses of the lepton produced in the CC interaction of the neutrino will be tracked in a single MARES run and combined to produce one voltage output file with the radar echo from all the energy dumps. 
