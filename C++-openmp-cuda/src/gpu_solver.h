#ifndef __GPU_SOLVER
#define __GPU_SOLVER
#include "real.h"
/*This is a collection of functions used to solve multi-phase diffusion using Cuda.
 *
 * GPU_Init: This function should be invoke first using the computational domain shapes before calling the other functions.
 * Return 1 if GPU computing is feasible;
 * Return 0 if GPU computing is unfeasible, because of no GPU avaible or insufficient memory.
 * Inside this function, GPU memory is allocated for a few device pointers.
 * GPU_Memcpy: This function is to copy data to GPU and is called inside Diffusion_solver().
 * GPU_Diffusion_solver: This is the function called to solve diffusion for a given period
 * (the period is usually the timestep for phase field equation, and is longer than tiemstep for diffusion).
 * GPU_Clearup: This function is to deallocate the device pointers used in the Diffusion_solver. It should be invoked at the end of the program.
 *
 *
 *
 */
extern int GPU_Init(const int x,const int y, const int z);//x,y,z are the extents of the domain.
extern void GPU_Diffusion_solver(real*hdiffusivity, real*hpartition_k,real*hpartition_c,real*hvalue_phase,real*htotal_conc,int x,int y, int z,real  dt,real dx);//
//hdiffusivity[max_num_phase] is solute diffusivity; hpartition_k and hpartition_c are the thermodynamic paritioning coefficients; hvalue_phase is the phase-fraction array(4D);
// htotal_conc is the concentration array(3D); x,y,z are the domain extents; dt is the time interval; dx is the grid spacing.

extern void GPU_Cleanup();

#endif