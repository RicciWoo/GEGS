/*=============================================================
  DISCRETE BHABHA SCATTERING (A CALL TO THIS ROUTINE) HAS BEEN
  ARBITRARILY DEFINED AND CALCULATED TO MEAN BHABHA SCATTERINGS
  WHICH IMPART TO THE SECONDARY ELECTRON SUFFICIENT ENERGY THAT
  IT BE TRANSPORTED DISCRETELY, I.E. E=AE OR T=TE. IT IS NOT
  GUARANTEED THAT THE FINAL POSITRON WILL HAVE THIS MUCH ENERGY
  HOWEVER. THE EXACT BHABHA DIFFERENTIAL CROSS SECTION IS USED.
  =============================================================*/



__global__ void bhabha()
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
	p.status = e_empty;
	bool process;

	if(idx.p==0) simu_stack_is_empty = false;

	for(;;){

		if(p.status==e_empty && !simu_stack_is_empty)
			p.status = e_new_particle;

		process = (p.status == e_new_particle);

		if(process && !pop_stack(p,e_bhabha))
			p.status = e_empty;

        if(simu_stack_is_empty){
            process = (p.status == e_empty);
            if(__all(process))
                break;
        }

		process = (p.status == e_bhabha);

#ifdef DO_STEP_COUNT
        uint count_mask = __ballot(process);
		if (idx.t == 0)step_counters[e_bhabha] += __popc(count_mask);
#endif

		if(process){

			region_data_t reg_dat = region_data[p.region];

			//float peip,pekin,pekse2,pese1,pese2,h1,dcosth;
			//float ekin,t0,e0,e02,yy,y2,yp,yp2,beta2,ep0,b1,b2,b3,b4,ep0c,rnno03,rnno04,br,rejf2;
			//float costhe,sinthe;
			float edep=0.0f;

			particle_t p1=p;
			ushort medium=reg_dat.med;

			float eip=p.e;
			float ekin=eip-prm;
			float t0=ekin/rm;
			float e0=t0+1.0f;
			float yy=1.0f/(t0+2.0f);
			float e02=e0*e0;
			float beta2=(e02-1.0f)/e02;
			float ep0=te[medium]/ekin;
			float ep0c=1.0f-ep0;
			float y2=yy*yy;
			float yp=1.0f-2.0f*yy;
			float yp2=yp*yp;
			float b4=yp2*yp;
			float b3=b4+yp2;
			float b2=yp*(3.0f+y2);
			float b1=2.0f-y2;

			//sample br from minimum(ep0) to 1.0
			float rnno03,rejf2;
			float br;
			do{
				rnno03=curand_uniform(&randStat);
				br=ep0/(1.0f-ep0c*rnno03);
				//apply rejection function
				//rnno04=curand_uniform(randStat);
				rnno03=curand_uniform(&randStat);
				rejf2=(1.0f-beta2*br*(b1-br*(b2-br*(b3-br*b4))));
				//}while(rnno04>rejf2);
			}while(rnno03>rejf2);
			//if e- got more than e+, move the e+ pointer and reject b

			/*PUTS E+ ON TOP OF STACK IF IT HAS LESS ENERGY*/
			if(br<0.5f)                
			{
				//p.charge remains +1
				p1.charge=-1;
			}
			else
			{
				p.charge=-1;
				p1.charge=1;
				br=1.0f-br;
			}

			//divide up the energy
			br=max(br,0.0f);
			float ekse2=br*ekin;
			float ese1=eip-ekse2;
			float ese2=ekse2+prm;
			p.e=ese1;
			p1.e=ese2;
			p.status = e_electron_step;
			p1.status = e_electron_step;
			p.nRepeat = 0;  //no repeat for secondary particle
			p1.nRepeat = 0;  //no repeat for secondary particle

			//bhabha angles are uniquely determined by kinematics
			float h1=(eip+prm)/ekin;
			//direction cosine change for 'old' electron
			float dcosth=min(1.00f,h1*(ese1-prm)/(ese1+prm));
			float sinthe=sqrtf(1.0f-dcosth);
			float costhe=sqrtf(dcosth);
			float sinphi,cosphi;
			uphi21(costhe,sinthe,cosphi,sinphi,p,&randStat);
			if(p.e<reg_dat.ecut){
				edep += p.e - prm;
				if(p.charge>0)
					annih_at_rest(p,reg_dat,&randStat);  //p.status properly set at the end of annih_at_rest
			}
			else
				push_stack(p,e_electron_step);


			dcosth=h1*(ese2-prm)/(ese2+prm);
			sinthe=-sqrtf(1.0f-dcosth);
			costhe=sqrtf(dcosth);
			uphi32(costhe,sinthe,cosphi,sinphi,p1);
			if(p1.e<reg_dat.ecut){
				edep += p1.e - prm;
				if(p1.charge>0)
					annih_at_rest(p1,reg_dat,&randStat);
			}
			else
				push_stack(p1,e_electron_step);        //push the second electron into stack

			escore(p.region,edep*p.wt);
			p.status=e_empty;
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