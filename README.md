PhaseField
==========

The project is intended for implementing phase field simulations on different platforms.
Phase field approaches are mesoscale models that are used widely in materials science. It can model microstructure evolution in engineering materials during processing such as heat treatments. a great amount of information can be obtained from the model, e.g. volume fractions of various phases, solute concentration distribution and microstructures.
The phase field models are usually coupled with solute diffusion equations describing solute diffision in multiple phases.
Both models are systems of nonliner partial differential equations. A explicit finite-difference scheme is used to solve the equations numerically on a uniform grid.
Since the numerical computing is both compute- and data-intensive, various high performance computing (HPC) techniques are used to develop the codes, e.g. the multithreading library OpenMP, the Message Passing Interface (MPI), and the GPU computing (NVidia CUDA). Thus, the code can be compiled and run on different platforms.
