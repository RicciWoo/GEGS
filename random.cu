

__global__ void setup_kernel(curandState *state)
{
    //unsigned int ix = blockIdx.x*blockDim.x + threadIdx.x;
    //unsigned int iy = blockIdx.y*blockDim.y + threadIdx.y;
   
	//int ThreadIndex =ix + iy *(gridDim.x*blockDim.x);

	int BlockIndex = blockIdx.y * gridDim.x + blockIdx.x;
    int ThreadIndex =BlockIndex * blockDim.y * blockDim.x + threadIdx.y * blockDim.x + threadIdx.x;

	
	/* Each thread gets same seed, a different sequence number, no offset */
    curand_init(0, ThreadIndex , 0, &state[ThreadIndex ]);
	
}

void launch_setup_kernel ()
{
	d_devStates = (curandState**)malloc(GPUNo*sizeof(curandState*));

	for(int GPUId=0; GPUId<GPUNo; GPUId++) {
		
#ifdef USE_MULTIPLE_GPU
		cudaSetDevice(GPUId); ce(58002);
#endif

		//int size = sizeof(curandState);
		//printf("  size of curandState . . . . %d\n  total number of thread  . . %d\n" , size, SIMULATION_NUM_THREADS);
		cudaMalloc((void **)&d_devStates[GPUId], SIMULATION_NUM_THREADS*sizeof(curandState));

		// the grid and block configuration should match that of the simulation kernel.
		setup_kernel<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC, NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK * WARP_SIZE>>>(d_devStates[GPUId]); ce(11000);

		cudaMemcpyToSymbol(devStates, &d_devStates[GPUId], sizeof(curandState*)); ce(9106);

	}
}

/*
void launch_setup_kernel (curandState * state)
{
    // execute the kernel
    //dim3 block(16, 8, 1);
    //dim3 grid(mesh_width / block.x, mesh_height / block.y, 1);
    //setup_kernel<<< grid, block>>>(state);

	// the grid and block configuration should match that of the simulation kernel.
	setup_kernel<<<dim3(SIMULATION_BLOCKS_PER_MULTIPROC, NUM_MULTIPROC), SIMULATION_WARPS_PER_BLOCK * WARP_SIZE>>>(state); ce(11000);

}
*/