# SmoothedParticleHydrodynamics
Smoothed Particle Hydrodynamics (SPH) to simulate fluid media

A simple 2D (for now) SPH implementation that is based on the work of M. Muller, D. Charypar, and M. Gross et al. 

Uses OpenCL 1.2 with the C++ bindings. 

SPH_v1 -> Unoptimized but runs entirely on the GPU, can reach approximately 4K particles with 10 iterations of SPH solver. (Completed)
SPH_v2 -> Optimization using Spatial Grid Hashing. (In Progress)


