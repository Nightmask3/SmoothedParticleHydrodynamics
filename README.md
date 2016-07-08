# SmoothedParticleHydrodynamics
Smoothed Particle Hydrodynamics (SPH) to simulate fluid media

A simple 2D SPH implementation that is based on the work of M. Muller, D. Charypar, and M. Gross et al. 
Uses the kernel equations as described here:
http://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf

Right now this implementation is unoptimized, but runs entirely on the GPU. 
Will start working towards implementing Spatial Grid Hashing in order to reduce complexity of all-pairs distance calculation.
