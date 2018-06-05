/*===============================================================
  DISCRETE MOLLER SCATTERING (A CALL TO THIS ROUTINE) HAS BEEN 
  ARBITRARILY DEFINED AND CALCULATED TO MEAN MOLLER SCATTERINGS 
  WHICH IMPART TO THE SECONDARY ELECTRON SUFFICIENT ENERGY THAT
  IT BE TRANSPORTED DISCRETELY. THE THRESHOLD TO TRANSPORT AN 
  ELECTRON DISCRETELY IS A TOTAL ENERGY OF AE OR A KINETIC ENERGY
  OF TE=AE-RM. SINCE THE KINETIC ENERGY TRANSFER IS ALWAYS, BY 
  DEFINITION, LESS THAN HALF OF THE INCIDENT KINETIC ENERGY, THIS
  IMPLIES THAT THE INCIDENT ENERGY, EIE, MUST BE LARGER THAN 
  THMOLL=TE*2+RM. THE REST OF THE COLLISION CONTRIBUTION IS 
  SUBTRACTED CONTINUOUSLY FROM THE ELECTRON AS IONIZATION 
  LOSS DURING TRANSPORT.
  ===============================================================*/


__global__ void moller()
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

		if(process && !pop_stack(p,e_moller))
			p.status = e_empty;

        if(simu_stack_is_empty){
            process = (p.status == e_empty);
            if(__all(process))
                break;
        }

		process = (p.status == e_moller);

#ifdef DO_STEP_COUNT
        uint count_mask = __ballot(process);
		if (idx.t == 0) step_counters[e_moller] += __popc(count_mask);
#endif

		if(process){

			region_data_t reg_dat = region_data[p.region];
			
			//float eie,ekin,t0,e0,e02,g2,g3,gmax,br,r,rejf4,rnno27,rnno28,extrae;
			//float sigm,pbrem,rsh,uj,sig_j;
			//int   lelke,iele,ish,nsh,ifirst,i,jj,iz;
			//float peie,pekse2,pese1,pese2,pekin,h1,dcosth,elke;
			//float costhe,sinthe;
			float edep=0.0f;

			particle_t p1=p;
			ushort medium=reg_dat.med;

			float eie=p.e;
			float ekin=eie-prm;
			float elke=log(ekin);
			bool  loop_done = false;

			if((eii_flag>0)&&(eii_nsh[medium]>0))
			{
				// The EII flag is set and this medium has shells for which we want to
				// simulate EII => sample if the interaction is with a EII shell
				uint lelke=eke01[medium].y*elke+eke01[medium].x;
				//float sigm=esig[medium*MXEKE+lelke-1].y*elke+esig[medium*MXEKE+lelke-1].x;
				float sigm=sig[medium*MXEKE+lelke-1].y*elke+sig[medium*MXEKE+lelke-1].x;
				float pbrem=ebr1[medium*MXEKE+lelke-1].y*elke+ebr1[medium*MXEKE+lelke-1].x;
				sigm=sigm*(1-pbrem);
				float rsh=curand_uniform(&randStat);
				rsh = sigm*rsh;
				for(char iele=0;iele<nne[medium];iele++){
					short iz=int(zelem[iele*MXMED+medium]+0.5);
					short nsh=eii_no[iele*MXMED+medium];
					if(nsh>0.0f){
						short ifirst=eii_first[iele*MXMED+medium];
						for(short ish=1;ish<=nsh;ish++){
							float uj=binding_energies[(iz-1)*MXSHELL+ish-1];
							if((ekin>uj)&&(uj>te[medium]||uj>ap[medium])){
								short jj=ifirst+ish-1;
								short i=eii[jj-1].x*elke+eii[jj-1].y+(jj-1)*N_EII_BINS;
								float sig_j=eii_xsection[i-1].x*elke+eii_xsection[i-1].y;
								sig_j=sig_j*pz[MXMED*iele+medium]*eii_cons[medium];
								rsh-=sig_j;
								if(rsh<0.0f){
									eii_sample(ish,iz,uj,p,reg_dat,&randStat,edep);
									//p.status=e_electron_step;  //we already set p.status in eii_sample
									if(p.e<reg_dat.ecut)
										edep += p.e - prm;
									else
										push_stack(p,e_electron_step);

									escore(p.region,edep*p.wt);
									p.status = e_empty;
									loop_done = true;
									break;  //ish loop
								}
							}
						}
						if(loop_done) break;  //iele loop
					}
				}
				if(loop_done) continue;  //outer for(;;) loop
			}
			
			if(ekin<=2*te[medium])
			{
				edep = p.e - prm;
				escore(p.region,edep*p.wt);
				p.status = e_empty;
				continue;
			}

			float t0=ekin/rm;
			float e0=t0+1.0f;
			float extrae=eie-thmoll[medium];
			float e02=e0*e0;
			//ep0=te[medium]/ekin;
			float g2=t0*t0/e02;
			float g3=(2.0f*t0+1.0f)/e02;

			float gmax=(1.0f+1.25*g2);
			float rnno27,rejf4;
			float br,r;
			do{
				//to retry if rejected
				rnno27=curand_uniform(&randStat);
				br=te[medium]/(ekin-extrae*rnno27);
				//use messel and crawfords rejection function
				r=br/(1.0f-br);
				//rnno28=curand_uniform(randStat);
				rnno27=curand_uniform(&randStat);
				rejf4=(1.0f+g2*br*br+r*(r-g3));
				//rnno28=gmax*rnno28;
				rnno27=gmax*rnno27;
				//}while(rnno28>rejf4);
			}while(rnno27>rejf4);

			float ekse2=br*ekin;
			float ese1=eie-ekse2;
			float ese2=ekse2+prm;
			p.e = ese1;
			p.status=e_electron_step;
			p.nRepeat = 0;  //no repeat for secondary particle
			p1.e=ese2;
			p1.status=e_electron_step;
			p1.nRepeat = 0;  //no repeat for secondary particle
			//since br .LE. 0.5, p1.e must be .LE. p.e
			//moller angles are uniquely determined by kinematics

			float h1=(eie+prm)/ekin;
			//direction cosine change for 'old' electron
			float dcosth=h1*(ese1-prm)/(ese1+prm);
			float sinthe=sqrtf(1.0f-dcosth);
			float costhe=sqrtf(dcosth);		
			float sinphi,cosphi;
			uphi21(costhe,sinthe,cosphi,sinphi,p,&randStat);
			if(p.e<reg_dat.ecut)
				edep += p.e - prm;
			else
				push_stack(p,e_electron_step);

			dcosth=h1*(ese2-prm)/(ese2+prm);
			sinthe=-sqrtf(1.0f-dcosth);
			costhe=sqrtf(dcosth);
			uphi32(costhe,sinthe,cosphi,sinphi,p1);
			if(p1.e<reg_dat.ecut)
				edep += p1.e - prm;
			else
				push_stack(p1,e_electron_step);  //push the second electron to stack

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