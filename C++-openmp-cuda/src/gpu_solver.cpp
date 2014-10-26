#include <algorithm>

#include "constant.h"
#include <iostream>
#include "real.h"
using namespace std;

#ifdef __NVCC__

#include "cuda.h"
#include "cuda_runtime.h"

static dim3 block,grid;//The GPU kernel configuration, e.g. kernel<<<grid,block>>>().
static real *dvalue_phase,*dtotal_conc,*dvalue_conc;//Device pointers for phase-fraction, solute concentration and solute phase-concentration.
static int size_domain,size_phase;//The numbers of elements in total_conc array and phase-fration array;
static int dimensions;
__constant__ real dpara_partition_k[max_num_phase*max_num_phase],dpara_partition_c[max_num_phase*max_num_phase];//thermodynamic partitioning coefficients
__constant__ real ddiffusivity[max_num_phase];// solute diffusivity

/*This function does the following work:
 * 1) Check the Cuda GPU avaibility;
 * 2) Set a GPU if available and display its properties;
 * 3) Allocate memory for the device pointers;
 * 4) determine the kernel configuration parameters, i.e. grid and block.
 */
int GPU_Init(const int x,const int y, const int z) {
    size_domain=(x)*(y)*(z)*sizeof(real);
    size_phase=max_num_phase*size_domain;
    int total_byte=size_domain+size_phase*2;
    int count;
    cout<<"gput init"<<endl;
    count=0;
    cudaGetDeviceCount(&count);

    if(count==0) {
        cout<<"THERE ARE NO CUDA GPU AVAILABLE"<<endl;
        return 0;
        //exit(EXIT_SUCCESS);

    } else
        cout<<"There are "<<count<<" GPUs."<<endl;

    cudaSetDevice(0);
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop,0);
    cout<<"Device major compute capability: "<<prop.major<<endl;
    cout<<"Device minor compute capability: "<<prop.minor<<endl;
    cout<<"Global Memory(GBytes) is "<<prop.totalGlobalMem/1024/1024/1024<<endl;

    if(total_byte>prop.totalGlobalMem) {
        cout<<"GLOBAL MEMORY IS INSUFFICENT FOR CURRENT SIMULATION. GPU COMPUTING IS INHIBITED."<<endl;
        return 0;
    }

    cout<<"Async Engine Count: "<<prop.asyncEngineCount<<endl;

    if(prop.concurrentKernels>0)
        cout<<"Device can execute multiple kernels concurrently."<<endl;

    else
        cout<<"Device cannot execute multiple kernels concurrently."<<endl;

    if(prop.deviceOverlap)
        cout<<"Device can concurrently copy memory and execute a kernel."<<endl;

    else
        cout<<"Device cannot concurrently copy memory and execute a kernel."<<endl;

    if(prop.integrated)
        cout<<"Device is integrated."<<endl;

    else
        cout<<"Device is not integrated."<<endl;

    cout<<"L2 cache size (Kbytes): "<<prop.l2CacheSize/1024<<endl;
    cudaMalloc(&dtotal_conc,size_domain);
    cudaMalloc(&dvalue_conc,size_phase);
    cudaMalloc(&dvalue_phase,size_phase);
    block.x=256;
    block.y=1;
    grid.x=x/block.x;

    if(x%block.x>0)grid.x++;

    grid.y=y/block.y;

    if(y%block.y>0)grid.y++;

    grid.z=z/block.z;

    if(z%block.z>0)grid.z++;
    
    if(z>1)dimensions=3;
    if(z==1) dimensions=2;
    

    cout<<"GPU kernel Config: Grid<"<<grid.x<<","<<grid.y<<","<<grid.z<<"> Block<"<<block.x<<","<<block.y<<","<<block.z<<">"<<endl;
    return 1;
}

/***********************************************************
 * This function is to copy data from host to GPU.
 *
 *
 */
void GPU_Memcpy(real*hdiffusivity, real*hpartition_k,real*hpartition_c,real*hvalue_phase,real*htotal_conc) {
    cudaMemcpyToSymbol(dpara_partition_k,hpartition_k,sizeof(dpara_partition_k));
    cudaMemcpyToSymbol(dpara_partition_c,hpartition_c,sizeof(dpara_partition_c));
    cudaMemcpyToSymbol(ddiffusivity,hdiffusivity,sizeof(ddiffusivity));
    cudaMemcpy(dvalue_phase,hvalue_phase,size_phase,cudaMemcpyHostToDevice);
    cudaMemcpy(dtotal_conc,htotal_conc,size_domain,cudaMemcpyHostToDevice);
}


__device__ inline int maxI(int a, int b) {
    return a>b?a:b;
}

__device__ inline int minI(int a,int b) {
    return a<b?a:b;
}

__host__ __device__ inline int Get_neighbor(int i,int x,int n) {
    int ii=(i+n)>=x?(i+n)-x:i+n;
    return (ii<0?(ii+x):ii);
}

/*******************************************************************************
 * This Kernel function is to partition the solute-concentration to phase-concentration.
 */
__global__ void GPU_Partitioning(const real  * __restrict val_phase,const real* __restrict total_conc,real*val_conc,int x,int y,int z) {
    int row=blockIdx.y*blockDim.y+threadIdx.y;
    int slice=blockIdx.z*blockDim.z+threadIdx.z;
    int col=blockIdx.x*blockDim.x+threadIdx.x;
    int k=(col+x*(row+y*slice));
    real t1,t2;

    if(col<x&&row<y&&slice<z) { //Threads out of the domain will not implement the following code.
#pragma unroll
        for(int i=0; i<max_num_phase; i++) {
            t1=0.;
            t2=0.;
#pragma unroll

            for (int j=0; j<max_num_phase; j++) {
                t1+=val_phase[max_num_phase*k+j]*dpara_partition_k[i+j*max_num_phase];
                t2+=val_phase[max_num_phase*k+j]*dpara_partition_c[i+j*max_num_phase];
                //cout<<t1<<" "<<thermodynamic.para_partition_k(i,j)<<" ";
            }

            val_conc[max_num_phase*k+i]=(total_conc[k]-t2)/t1;
        }
    }
}

/*******************************************************************************
 * This Kernel function is to solve diffusion equation using explict finite difference.
 * 7-stencil scheme is used, i.e. using the nearest 6 neigbors to calculate Laplacian.
 */
__global__ void Diff3D(const real  * __restrict val_phase,real*total_conc,const real* __restrict val_conc,const int x,const int y,const int z,const real dt_dxdx,const int num_loop) {
    int row=blockIdx.y*blockDim.y+threadIdx.y;
    int slice=blockIdx.z*blockDim.z+threadIdx.z;
    int col=blockIdx.x*blockDim.x+threadIdx.x;
    int k=(col+x*(row+y*slice));
    int neighbors[6]= {Get_neighbor(col,x,-1)+x*(row+y*slice),
                       Get_neighbor(col,x,+1)+x*(row+y*slice),
                       (col+x*(Get_neighbor(row,y,-1)+y*slice)),
                       (col+x*(Get_neighbor(row,y,+1)+y*slice)),
                       (col+x*(row+y*Get_neighbor(slice,z,-1))),
                       (col+x*(row+y*Get_neighbor(slice,z,+1)))
                      };//Get the indices of neigbors in the 1D array.
    real del_conc=0.;
    int loc;

    if(col<x&&row<y&&slice<z) { //Threads out of the domain will not implement the following code.
#pragma unroll
        for(int ii=0; ii<6; ii++) { // Calculate the flux on each side
            loc= neighbors[ii]*max_num_phase;
#pragma unroll

            for(int jj=0; jj<max_num_phase; jj++) //For one side, calculate the flux in each phase
                del_conc+=sqrt( val_phase[loc+jj]*val_phase[k*max_num_phase+jj] )*
                          ( val_conc[loc+jj]-val_conc[k*max_num_phase+jj] )*ddiffusivity[jj];
        }

        total_conc[k]+=(del_conc*dt_dxdx);
    }
}

/**************************************************************************
 * This function is to coordinate the GPU for calculations.
 * First, copy necessary data to GPU;
 * Second, call kernel functions for data crunching;
 * Finaly, copy back the total concontration to Host.
 * hdiffusivity: diffusivity array on host;
 * hpartition_k: partitioning coefficients on host;
 * hpartition_c: partitioning concentrations on host;
 * hvalue_phase: phase-fraction array on host;
 * htotal_conc: total-concentration array on host;
 * x,y,z: domain sizes;
 * dt: time interval;
 * dx: grid spacing.
 */
void GPU_Diffusion_solver(real*hdiffusivity, real*hpartition_k,real*hpartition_c,real* hvalue_phase,real*htotal_conc,int x,int y, int z,real  dt, real dx) {
    GPU_Memcpy(hdiffusivity, hpartition_k,hpartition_c,hvalue_phase,htotal_conc);//Copy to GPU
    real max_diff=*max_element(hdiffusivity,hdiffusivity+max_num_phase);//Find out the maximum diffusivity
    real del_t=dx*dx/max_diff/dimensions*0.2f;// calculate the timestep for diffusion equations, ensuring stability and convergency.
    int num_loop=ceil(dt/del_t);// calculate the number of iterations to solve diffusion equation for a time interval of dt.
    real dt_dxdx=dt/num_loop/(dx*dx);// is dt/(dx*dx).
    cudaFuncSetCacheConfig(Diff3D,cudaFuncCachePreferL1);

    for(int n=1; n<=num_loop; n++) {
        GPU_Partitioning<<<grid,block>>>(dvalue_phase,dtotal_conc,dvalue_conc,x,y,z);
        Diff3D<<<grid,block>>>(dvalue_phase,dtotal_conc,dvalue_conc,x,y,z,dt_dxdx,num_loop);
    }

    cudaMemcpy(htotal_conc,dtotal_conc,size_domain,cudaMemcpyDeviceToHost);
}

/**************************************************************************
 * This function is to clean up the GPU memory and reset the device.
 */
void GPU_Cleanup() {
    cudaError_t freeError;
    freeError=cudaFree(dvalue_conc);

    if(freeError!=cudaSuccess)exit(0);

    freeError=cudaFree(dvalue_phase);

    if(freeError!=cudaSuccess)exit(0);

    freeError=cudaFree(dtotal_conc);

    if(freeError!=cudaSuccess)exit(0);

    cudaDeviceReset();
    cout<<"Free all allocated memory on GPU."<<endl;
}

/***************************************************************************
 * The following functions are used for debugging.
 */
static void GPU_CpyToHost_valconc(real *val_conc) {
    cudaMemcpy(val_conc,dvalue_conc,size_phase,cudaMemcpyDeviceToHost);
}
static void GPU_CpyToHost_valphase(real *val_phase) {
    cudaMemcpy(val_phase,dvalue_phase,size_phase,cudaMemcpyDeviceToHost);
}

#else
int GPU_Init(const int x,const int y, const int z){
  cout<<"NO GPU is available"<<endl;
  return 0;
}
void GPU_Memcpy(real*hdiffusivity, real*hpartition_k,real*hpartition_c,real*hvalue_phase,real*htotal_conc) {}
void GPU_Diffusion_solver(real*hdiffusivity, real*hpartition_k,real*hpartition_c,real* hvalue_phase,real*htotal_conc,int x,int y, int z,real  dt, real dx){}
void GPU_Cleanup() {}
#endif
