# SmoothedParticleHydrodynamics
Smoothed Particle Hydrodynamics (SPH) to simulate fluid media

A simple 2D SPH implementation that is based on the work of M. Muller, D. Charypar, and M. Gross et al. 
Uses the kernel equations as described here:
http://www.cs.cornell.edu/~bindel/class/cs5220-f11/code/sph.pdf

Right now this implementation is unoptimized and runs entirely on the CPU. Will immediately begin work on writing an OpenCL kernel that will allow offloading of the 
interaction force computation onto the GPU.

Have begun writing a OpenCL abstraction called 'CLWrapper' which will eventually contain all the funtionality required in order to initialize, build and compile a kernel, 
request, read and write buffers etc.
