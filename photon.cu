// Transport the photon one step through the phantom and determine which (if any) interaction
// takes place next.
// This is the subroutine PHOTON in the file egsnrc.mortran (v 1.72 2011/05/05) of the EGSnrc 
// code system.

__global__ void photon_step()
{
	indices idx = get_indices();
	curandState randStat=devStates[idx.p];  //copy the random number generator state to local.
	
	__syncthreads();

#ifdef DO_STEP_COUNT
    volatile uint *step_counters = step_counters_shared[idx.w];

    // reset step counts
	if (idx.t < NUM_CAT) step_counters[idx.t] = 0;
#endif

	region_data_t reg_dat;
	particle_t p_original, p;
	p_original.nRepeat = 0;
	p.status = p_empty;
	bool process;
	
	if(idx.p==0) simu_stack_is_empty = false;

	for(;;){

		if(p.status==p_empty && (!simu_stack_is_empty||p.nRepeat>0))
			p.status = p_new_particle;

		process = (p.status == p_new_particle);
		
		if(process){
			if(p_original.nRepeat>0){
				p_original.nRepeat --;
				p = p_original;
				p.nRepeat=0;
			}
			else{
				do{
					if(pop_stack(p_original,p_photon_step)){
						if(p_original.charge!=0)
							push_stack(p_original,e_electron_step);
						else
							{
								p = p_original;
								p.nRepeat=0;
						}
					}
					else{
						p.status = p_empty;
						//process = false;
						break;
					}
				}while(p_original.charge!=0);
			}
		}

		bool newparticle = (p.status == p_new_particle);

#ifdef DO_STEP_COUNT
		uint count_mask = __ballot(newparticle);
		if (idx.t == 0) step_counters[p_new_particle] += __popc(count_mask);
#endif

		if(newparticle) p.status = p_photon_step;

        if(simu_stack_is_empty){
            process = (p.status == p_empty);
            if(__all(process))
                break;
        }

		process = (p.status == p_photon_step);

#ifdef DO_STEP_COUNT
        count_mask = __ballot(process);
		if (idx.t == 0) step_counters[p_photon_step] += __popc(count_mask);
#endif

		if(process){
			reg_dat = region_data[p.region];

			//float r1 =  curand_uniform(&randStat); 
			//float r2 =  curand_uniform(&randStat); 
			//float r3 =  curand_uniform(&randStat); 
			float edep = 0.0f;

			if (p.e <= reg_dat.pcut){
				edep = p.e;
				escore(p.region,edep*p.wt);
				p.status = p_empty;
				continue;
			}
			else if (p.wt <= 0.0F){  //when p.wt<=0.0, it is no need to score the energy deposition
				p.status = p_empty;
				continue;
			}

			float r1 =  curand_uniform(&randStat);
			float dpmfp = -logf(r1);

			float gle = logf(p.e);
			int   lgle = 0;
			//float gmfpr0 = 0.0F;
			float tstep = 0.0F;
			float gmfp_val = 0.0F;
			float cohfac = 0.0F;
			ushort old_medium = reg_dat.med;

			if (reg_dat.med == VACUUM)
				tstep = VACUUM_STEP;
			else {
				float2 ge_dat = ge[reg_dat.med];
				lgle = int(ge_dat.x + ge_dat.y * gle);

				float2 gmfp_dat = gmfp[reg_dat.med * MXGE + lgle];
				float gmfpr0 = gmfp_dat.x + gmfp_dat.y * gle;

				gmfp_val = gmfpr0 / reg_dat.rhof;

				if ((reg_dat.flags & f_rayleigh) > 0) {
					float2 cohe_dat = cohe[reg_dat.med * MXGE + lgle];
					cohfac = cohe_dat.x + cohe_dat.y * gle; 
					gmfp_val *= cohfac;
				}
				tstep = gmfp_val * dpmfp;
			}

			// HOWFAR
			uint new_region = 0;
			//int ix,iy,iz;
			/*
			if(p.leafIndex<-60 && p.leafIndex>=PASSMLC)
				//new_region = howfar_phantom(p,tstep,ix,iy,iz);
				new_region = howfar_phantom(p,tstep);
			else
				new_region = howfar_MLC(p,tstep);
			*/
			if(p.module==m_SecJawY)
				new_region = howfar_SecJawsY(p,tstep);
			else if(p.module==m_SecJawX)
				new_region = howfar_SecJawsX(p,tstep);
			else if(p.module==m_VarMLCs)
				new_region = howfar_MLC(p,tstep);
			else if(p.module==m_BlockMd)
				new_region = howfar_block(p,tstep);
			else if(p.module==m_WedgeMd)
				new_region = howfar_wedge(p,tstep);
			else if(p.module==m_Phantom)
				new_region = howfar_phantom(p,tstep);
			else{
				printf("error in module status!\n");
				p.status = p_empty;
				continue;
			}

			char idisc = 0;
			if (new_region == 0) {							//if going into VACUUM , commented by Tong Xu, Jan 2013
				if (reg_dat.med==VACUUM || p.leafIndex==PASSMLC)		// and current region is also VACUUM , or it has passed MLC and didn't hit the phantom commented by Tong Xu, Jan 2013
					idisc = 1;												// discard ! , commented by Tong Xu, Jan 2013
				else
					idisc = -1;
			}

			if (idisc > 0) {
				p.region = 0;  //0 means Vacuum
				edep = p.e;
				escore(p.region,edep*p.wt);
				p.status = p_empty;
				continue;
			}
			else {
				p.x += p.u * tstep;
				p.y += p.v * tstep;
				p.z += p.w * tstep;

				if (reg_dat.med != VACUUM) 
					dpmfp = max(0.0F, dpmfp - tstep / gmfp_val);

				old_medium = reg_dat.med;

				if (new_region != p.region) {
					p.region = new_region;
					reg_dat = region_data[new_region];
					/*
					if(p.region<=2)   //if the region belong to MLC
						reg_dat = region_data[new_region];
					else  //if the particle is in the phantom, we use texture fetching, hopefully faster
					{
						int2 tmp1=tex3D(med_flags_tex, ix,iy,iz);
						reg_dat.med = tmp1.x;
						reg_dat.flags = tmp1.y;
						float4 tmp2 = tex3D(rhof_rho_pcut_ecut_tex,ix,iy,iz);
						reg_dat.rhof = tmp2.x;
						reg_dat.rho  = tmp2.y;
						reg_dat.pcut = tmp2.z;
						reg_dat.ecut = tmp2.w;
					}
					*/
				}

				/*
				if (p.e <= reg_dat.pcut){
					edep = p.e;
					escore(p.region,edep*p.wt);
					p.status = p_empty;
					continue;
				}
				else if (idisc < 0){
					edep = p.e;
					escore(p.region,edep*p.wt);
					p.status = p_empty;
					continue;
				}
				*/
				if (p.e<=reg_dat.pcut || idisc<0){
					edep = p.e;
					escore(p.region,edep*p.wt);
					p.status = p_empty;
					continue;
				}
			}

			// determine next step if not already discarded

			if ((reg_dat.med != old_medium) || (reg_dat.med == VACUUM) || (dpmfp >= EPSGMFP)) {
				p.status = p_photon_step;
				continue;
			}

			if ((reg_dat.flags & f_rayleigh) > 0) {
				r1 =  curand_uniform(&randStat);
				//if (r2 < 1.0F - cohfac) {
				if (r1 < 1.0F - cohfac) {
					p.status = p_rayleigh;
					push_stack(p,p_rayleigh);
					p.status = p_empty;
					continue;
				}
			}

			float2 gbr1_dat = gbr1[reg_dat.med * MXGE + lgle];
			float gbr1_val = gbr1_dat.x + gbr1_dat.y * gle;
			float2 gbr2_dat = gbr2[reg_dat.med * MXGE + lgle];
			float gbr2_val = gbr2_dat.x + gbr2_dat.y * gle;

			r1 =  curand_uniform(&randStat);
			//if ((r3 <= gbr1_val) && (p.e > 2.0F * ELECTRON_REST_MASS_FLOAT)){
			if ((r1 <= gbr1_val) && (p.e > 2.0F * ELECTRON_REST_MASS_FLOAT)){
				p.status = p_pair;
				push_stack(p,p_pair);
				p.status = p_empty;
				//continue;
			}
			else {
				//if (r3 < gbr2_val){
				if (r1 < gbr2_val){
					p.status = p_compton;
					push_stack(p,p_compton);
					p.status = p_empty;
					//continue;
				}
				else{
					p.status = p_photo;
					push_stack(p,p_photo);
					p.status = p_empty;
					//continue;
				}
			}
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


