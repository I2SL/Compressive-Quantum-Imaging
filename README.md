# Compressive Quantum Imaging
Author: Nicolas Deshler

Please see [the manuscript](https://www.overleaf.com/read/prtpjfbtbpjy).

This package unifies ideas from compressive sensing and quantum parameter estimation to passively image incoherent distributed scenes at multiple resolution levels extending beyond the diffraction limit. We take inspiration from [[1]](https://iopscience.iop.org/article/10.1088/1367-2630/aa60ee) and use an adaptive bayesian approach to estimate parameters of the scene. Since natural images are generally compressible in a wavelet basis, this algorithm adaptively estimates wavelet coefficients while enforcing a sparsity prior. The bayesian framework for quantum parameter estimation is detailed in [[2]](https://ieeexplore.ieee.org/document/1054643).

This algorithm employs Spatial Mode Demultiplexing (SPADE) detailed in [[3]](https://iopscience.iop.org/article/10.1088/1367-2630/aa60ee) to decompose the optical field at the focal plane of the imaging system into an orthogonal modal basis. Our quantum measurement consists of counting the number of photons that appear in each mode.



# Algorithm Features

- Employs photon counting measurements on transverse spatial modes to outperform direct imaging
- Parameters to be estimated are surrogates to the wavelet coefficients of the image related by a linear transformation that preserve the trace 1 norm of the density operator and the non-negativity requirement of the object intensity distribution.
- Uses adaptive Bayesian framework founded on Personick quantum parameter estimation theory to update joint measurement operator
- Sparsity prior imposed on the surrogate coefficients
- Markov-Chain-Monte-Carlo methods used to sample from the posterior distribution

# System Requirements, Software Versions, Package Extensions
- 8GB RAM
- matlab 2017 or higher
- [mtimesx](https://www.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support) fast matrix multiplication extension for matlab

# Installation and Setup

To download the repository in a desired local directory `my/dir/path/` open a Git terminal run
```
cd my/dir/path/
git clone https://github.com/I2SL/Compressive-Quantum-Imaging.git
```
Enter the project directory 
```
cd Compressive-Quantum-Imaging
```
and install the 'mtimesx' package (a fast matrix-multiply package implemented in C) 
```
matlab setup_mtimesx.m
``` 
If the last command doesn't work then matlab does not have an environment variable PATH on your machine. In this case, just open matlab for your system and run `setup_mtimesx.m` from the matlab terminal.


# Running the code
This section is set up to demonstrate our compressive quantum superresolution imaging algorithm in simulation. The target scene (a square graysacle image) is the only required input parameter for the program. There are two assumptions that this program makes:
- The target scene exists entirely within the sub-Rayleigh regime. That is, the angular extent of the image is equal to the rayleigh limit. 
- The dimensions of the image are a power of 2
```
matlab main()
```
There are a collection of optional arguments that a user may choose to include as well for better performance.


# References
1) K. K. Lee, S. Guha, and A. Ashok, "Quantum-inspired Optical Super-resolution Adaptive Imaging," In: OSA Imaging and Applied Optics Congress (2021)

2) S Personick. “Application of quantum estimation theory to analog commu-
nication over quantum channels”. In: IEEE Transactions on Information
Theory 17.3 (1971), pp. 240–246.


3) Mankei Tsang. “Subdiffraction incoherent optical imaging via spatial-mode
demultiplexing”. In: New Journal of Physics 19.2 (Feb. 2017), p. 023054.
doi: 10.1088/1367- 2630/aa60ee. url: https://doi.org/10.1088/
1367-2630/aa60ee.

4) Jesus Rubio and Jacob Dunningham. “Bayesian multiparameter quantum
metrology with limited data”. In: Physical Review A 101.3 (Mar. 2020).
doi: 10.1103/physreva.101.032114. url: https://doi.org/10.1103%
2Fphysreva.101.032114.
