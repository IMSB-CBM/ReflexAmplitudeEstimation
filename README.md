### Description
This repository contains code and data used for the publication *"Schmid et al. (2025), The variability of motor neuron reflex amplitude estimates in motor unit pools depends on the phenotype distribution and discharge statistics"*

### Installation
1. You need a working MATLAB installation

2. Clone the repository via:
```bash
git clone https://github.com/IMSB-CBM/ReflexAmplitudeEstimation.git
```

### Requirements
1. MATLAB 2021a
2. Signal Processing Toolbox
3. Parallel Computing Toolbox (only required for running the software on more than one core)
4. MATLAB Parallel Server (only required for running the software on more than one core)


### Repository Structure 
```
ReflexAmplitudeEstimation/
├── simulation/     % This folder contains the code to run the simulations used in the associated publication.
├── LICENSE
└── README.md
```

### How to use 
The scripts required to run the simulation:
- Run_PSC2spiketrain_parallel_rebound.m
    . This is the main script, which has to be executed.
    . Simulation settings are defined here.
    . By default runs on a single core but can be run on more cores for speed-up. The perturbation cycles will be parallelized. 
- MotorNeuronPool.m 
    . Class file for the motoneuron pool.
- CorticalInput.m
    . Helper function
    . Creates the drive to the motoneuron.
- pot2freq.m
    . Helper function
    . Determines the firing times of the motoneuron from the time course of the soma membrane potential.
- PeriStimAnalysis.m
    . Helper function
    . Applies peristimulus analysis to the firing times to create the peristimulus frequencygram (PSF) and the peristimulus time histogram (PSTH).
- reflex_amplitude.m
    . Helper function
    . Calculates the reflex amplitudes from the cumulative sum of PSF and PSTH.

