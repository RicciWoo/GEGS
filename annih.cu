/*===============================================================
  GAMMA SPECTRUM FOR TWO GAMMA IN-FLIGHT POSITRON ANNIHILATION. 
  USING SCHEME BASED ON HEITLER'S P269-27O FORMULAE.
  If the user requests radiative splitting (via nbr_split > 1), 
  this routine produces 2*nbr_split annihilation photons at once,
  each carying the fraction 1/nbr_split of the weight of the 
  incident positron. 
  Except for taking out the calculation of 
  LOG((1.0-EP0)/EP0) out of the sampling loop and using a 
  rejection function normalized to its maximum, the sampling 
  technique is the same as the original EGS4 implementation.
  ===============================================================*/



__global__ void annih()
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

		if(process && !pop_stack(p,e_annih))
			p.status = e_empty;

        if(simu_stack_is_empty){
            process = (p.status == e_empty);
            if(__all(process))
                break;
        }

		process = (p.status == e_annih);

#ifdef DO_STEP_COUNT
        uint count_mask = __ballot(process);
		if (idx.t == 0) step_counters[e_annih] += __popc(count_mask);
#endif

		if(process){

			region_data_t reg_dat = region_data[p.region];

			//float pavip,pesg1,pesg2;
			//float avip,a,g,t,pp,pot,ep0,wsamp,rnno01,rnno02,ep,rejf,esg1,esg2;
			//float aa,bb,cc,sinpsi,sindel,cosdel,us,vs,cphi,sphi;
			//int   ibr;
			//float xphi,xphi2,yphi,yphi2,rhophi2;
			//float costhe,sinthe;
			float edep=0.0f;

			if(nbr_split<=0)
			{
				p.status = e_empty;
				return ;
			}

			float avip=p.e+prm;   // PRECISE AVAILABLE ENERGY OF INCIDENT POSITRON,i.e. electron assumed to be at rest
			float a=avip/rm;
			float g=a-1.0f;
			float t=g-1.0f;
			float pp=sqrtf(a*t);
			float pot=pp/t;
			float ep0=1.0f/(a+pp);
			float wsamp=logf((1.0f-ep0)/ep0);

			float aa=p.u;
			float bb=p.v;
			float cc=p.w;
			float sinpsi=aa*aa+bb*bb;
			float sindel,cosdel;
			if(sinpsi>1e-20)
			{
				sinpsi=sqrtf(sinpsi);
				sindel=bb/sinpsi;
				cosdel=aa/sinpsi;
			}

			if(nbr_split>1)
			{
				p.wt=p.wt/nbr_split;
			}

			for(char ibr=1;ibr<=nbr_split;ibr++)
			{
				//nbr_split>1 means we want splitting for any radiative event	
				float rnno01,rejf,ep;
				do{
					rnno01=curand_uniform(&randStat);
					ep=ep0*exp(rnno01*wsamp);
					//now decide whether to accept
					//rnno02=curand_uniform(randStat);
					rnno01=curand_uniform(&randStat);
					rejf=1.0f-(ep*a-1.0f)*(ep*a-1.0f)/(ep*(a*a-2));
					//}while(rnno02>rejf);
				}while(rnno01>rejf);

				//set up energies
				float esg1=avip*ep;
				p.e=esg1;
				p.charge=0;                 
				p.status = p_photon_step;
				p.nRepeat = 0;  //no repeat for secondary particle

				float costhe=min(1.0f,(esg1-rm)*pot/esg1);
				float sinthe=sqrtf(1.0f-costhe*costhe);

				float xphi,yphi,rhophi2;
				do{
					xphi=curand_uniform(&randStat);
					xphi = 2*xphi - 1;
					//xphi2 = xphi*xphi;
					yphi=curand_uniform(&randStat);
					//yphi2 = yphi*yphi;
					//rhophi2 = xphi2 + yphi2;
					rhophi2 = xphi*xphi + yphi*yphi;
				} while (rhophi2>1);
				rhophi2 = 1/rhophi2;
				//cphi = (xphi2 - yphi2)*rhophi2;
				float cphi = (xphi*xphi - yphi*yphi)*rhophi2;
				float sphi = 2*xphi*yphi*rhophi2;

				if(sinpsi>=1e-10)
				{
					float us = sinthe*cphi;
					float vs = sinthe*sphi;
					p.u = cc*cosdel*us - sindel*vs + aa*costhe;        
					p.v = cc*sindel*us + cosdel*vs + bb*costhe;        
					p.w = cc*costhe - sinpsi*us;
				}
				else
				{
					p.u = sinthe*cphi;
					p.v = sinthe*sphi;
					p.w = cc*costhe;
				}
				if(p.e<reg_dat.pcut)
					edep += p.e;
				else
					push_stack(p,p_photon_step);

				float esg2 = avip-esg1;
				p.e = esg2;

				costhe=min(1.0f,(esg2-rm)*pot/esg2);
				sinthe=-sqrtf(1.0f-costhe*costhe);
				if(sinpsi>=1e-10)
				{   
					float us = sinthe*cphi;
					float vs = sinthe*sphi;
					p.u = cc*cosdel*us - sindel*vs + aa*costhe;    
					p.v = cc*sindel*us + cosdel*vs + bb*costhe;    
					p.w = cc*costhe - sinpsi*us;
				}
				else 
				{
					p.u = sinthe*cphi;     
					p.v = sinthe*sphi;
					p.w = cc*costhe;
				}
				if(p.e<reg_dat.pcut)
					edep += p.e;
				else
					push_stack(p,p_photon_step);            //因为前面已经有了 np+1的操作，所以此处的np就为原本的np+1
			}     

			escore(p.region,edep*p.wt);
			p.status = e_empty;
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