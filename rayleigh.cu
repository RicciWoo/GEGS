// Perform a Rayleigh interaction.
// This is the subroutine egs_rayleigh_sampling in the file egsnrc.mortran (v 1.72 2011/05/05) 
// of the EGSnrc code system.
__global__ void rayleigh()
{
   	indices idx = get_indices();
	curandState randStat=devStates[idx.p];  //copy the random number generator state to local.
	
	__syncthreads();

#ifdef DO_STEP_COUNT
    volatile uint *step_counters = step_counters_shared[idx.w];

    // reset step counts
    if (idx.t < NUM_CAT) step_counters[idx.t] = 0;
#endif

	particle_t p;
	p.status = p_empty;
	bool process;

	if(idx.p==0) simu_stack_is_empty = false;

	for(;;){

		if(p.status==p_empty && !simu_stack_is_empty)
			p.status = p_new_particle;

		process = (p.status == p_new_particle);

		if(process && !pop_stack(p,p_rayleigh))
			p.status = p_empty;

        if(simu_stack_is_empty){
            process = (p.status == p_empty);
            if(__all(process))
                break;
        }

		process = (p.status == p_rayleigh);

#ifdef DO_STEP_COUNT
        uint count_mask = __ballot(process);
		if (idx.t == 0) step_counters[p_rayleigh] += __popc(count_mask);
#endif

		if(process){

			region_data_t reg_dat = region_data[p.region];

			float edep=0.0f;
    
			p.status = p_photon_step;
			p.nRepeat = 0;    //no repeat for secondary particle
    
			//float xmax = 0.0F;
			//float pmax_val = 0.0F;

			float gle = logf(p.e);
			float2 ge_dat = ge[reg_dat.med];
			int lgle = int(ge_dat.x + ge_dat.y * gle);

			float2 pmax_dat = pmax[reg_dat.med * MXGE + lgle];
			float pmax_val = pmax_dat.x + pmax_dat.y * gle; 
			float xmax = HC_INVERSE * p.e;

			//int dwi = RAYCDFSIZE - 1;
			//int ibin = 0;
			//int ib = 0;

			//float xv = 0.0F;
			float costhe = 0.0F;
			float costhe2 = 0.0F;
			float sinthe = 0.0F;

			bool loop_done = false;

			do {
				bool inner_loop_done = loop_done;

				float xv = 0.0F;

				do {
					float r1 =curand_uniform(&randStat);
            
					if (!inner_loop_done) {
						float temp = r1 * pmax_val;
						// indexing in C starts at 0 and not 1 as in FORTRAN
						int dwi = RAYCDFSIZE - 1;
						int ibin = int(temp * (float)dwi);
						int ib = i_array[reg_dat.med * RAYCDFSIZE + ibin] - 1;
						int next_ib = i_array[reg_dat.med * RAYCDFSIZE + ibin + 1] - 1;

						if (next_ib > ib) {
							do {
								rayleigh_data_t ray_dat = rayleigh_data[reg_dat.med * MXRAYFF + ib + 1];
								if ((temp < ray_dat.fcum) || (ib >= RAYCDFSIZE - 2))
									break;
								ib++;
							} while (true);
						}

						rayleigh_data_t ray_dat = rayleigh_data[reg_dat.med * MXRAYFF + ib];
						temp = (temp - ray_dat.fcum) * ray_dat.c_array;
						xv = ray_dat.xgrid * expf(logf(1.0F + temp) * ray_dat.b_array);

						if (xv < xmax)
							inner_loop_done = true;
					}
            
				} while (!inner_loop_done);

				float r2 = curand_uniform(&randStat);
        
				if (!loop_done) {
					xv = xv / p.e;
					costhe = 1.0F - TWICE_HC2 * xv * xv;
					costhe2 = costhe * costhe;

					if (2.0F * r2 < 1.0F + costhe2)
						loop_done = true;
				}

			} while (!loop_done);

			sinthe = sqrtf(1.0F - costhe2);
			float cosphi,sinphi;
			uphi21(costhe,sinthe,cosphi,sinphi,p,&randStat);
			if(p.e<reg_dat.pcut)
			{
				edep = p.e;
				escore(p.region,edep*p.wt);
			}
			else
				push_stack(p,p_photon_step);

			p.status = p_empty;
		}
	}

	__syncthreads();

#ifdef DO_STEP_COUNT
    // combine the counters in shared memory and write them to global memory
    if (threadIdx.x < NUM_CAT) {
        uint total_count = 0;

        // step through the warps
        for (uchar i = 0; i < SIMULATION_WARPS_PER_BLOCK; i++)
            total_count += step_counters_shared[i][threadIdx.x];

        (*total_step_counts)[idx.b][threadIdx.x] = total_count;
    }
#endif

	devStates[idx.p]= randStat; //copy the random number generator states back to the array
}
