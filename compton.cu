/*===============================================================
  Subroutine for sampling incoherent (Compton) scattering 
  If the flag ibcmp(region) is zero, Klein-Nishina is used. 
  Otherwise scattering is modelled in the impulse approximation 
  (see R.Ribberfors and K.F.Berggren, Phys. Rev. A26 (1982) 3325) 
  As the total cross section from PEGS4 is not modified (and thus
  calculated using Klein-Nishina), all rejections leed to an 
  unscattered photon and a zero energy electron. 
  If a K, L1, L2, L3, M or N vacancy is created, the subsequent 
  atomic relaxation is treated in RELAX. This has as a 
  consequence that more than one particle can be created as a 
  result of an incoherent scattering. The user should therefore 
  check their user codes for possible inconsistencies
  ===============================================================*/


__global__ void compton()
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

		if(process && !pop_stack(p,p_compton))
			p.status = p_empty;

        if(simu_stack_is_empty){
            process = (p.status == p_empty);
            if(__all(process))
                break;
        }

		process = (p.status == p_compton);

#ifdef DO_STEP_COUNT
        uint count_mask = __ballot(process);
		if (idx.t == 0) step_counters[p_compton] += __popc(count_mask);
#endif

		if(process){

			region_data_t reg_dat = region_data[p.region];

			//float peig,pesg,pese;
			//float ko,  broi,  broi2,  bro,  bro1,  alph1,  alph2,  alpha,  rnno15;
			//float rnno16,rnno17,rnno18,rnno19,  br,  temp,  rejf3,  rejmax,  uj;
			//float jo,  br2,  fpz,fpz1, qc,  qc2,  af,  fmax1,  frej,  eta_incoh, eta;
			//float aux,aux1,aux2,aux3,aux4,pzmax,  pz,  pz2;
			//int   i,j;
			//bool  first_time=false;
			//int   ibcmpl;
			//float sinthe,costhe;
			float edep=0.0f;

			particle_t p1=p;
			ushort medium=reg_dat.med;

			float peig=p.e;
			float ko=peig/rm;
			float broi=1+2*ko;    //Needed for scattering angle sampling
			bool first_time = true;

			char ibcmpl;
			char tmp = reg_dat.flags&f_bound_compton_4;
			switch( tmp)
			{
				case f_bound_compton   : ibcmpl=1; break;
				case f_bound_compton_2 : ibcmpl=2; break;
				case f_bound_compton_3 : ibcmpl=3; break;
				default                : ibcmpl=0;
			}
	
//resample_compton:
			for(int compton_loops=0;;compton_loops++)
//			for(;;)
			{
				if(compton_loops>100){
					//printf("compton_loops>100! break! esig=%f\n",peig);
					break;
				}
				float uj,jo;
				short j;
				if(ibcmpl>0)    //User wants to take into account binding effects
				{               //first sample the shell and see whether an interaction is possible
					float rnno17=curand_uniform(&randStat);
					rnno17 = 1+rnno17*n_shell[medium];
					char i = int(rnno17);
					if(rnno17>eno_array[medium*MXMDSH+i-1])
					{
						i=eno_atbin_array[medium*MXMDSH+i-1];
					}
					j=shell_array[medium*MXMDSH+i-1];    //j is the shell number in the data list
					uj=be_array[j-1];                    //Uj is the binding energy in units of rm

					//Binding energy rejection
					if(ko<=uj)
					{
						if(ibcmpl==1)
							{			//reject this compton scatter, put the photon back to photon_step stack
										//Tong added these 2 sentences, July 2014
								p.status=p_photon_step;
								push_stack(p,p_photon_step);
								break;
							}
						else continue;//resample_compton loop
					}
					jo=jo_array[j-1];    //Jo is the Compton profile parameter
				}

				//We always sample the scattering angle from Klein-Nishina
//RESAMPLE:
				float bro,br,temp;
				float sinthe,costhe;
				int br_loops=0;
				do
				{
					if(br_loops>100)
					{
						//printf("br_loops >100! break\n"); 
						break;
					}
					if(ko>2)        //At high energies the original EGS4 method is most efficient
					{
						float alph1,alpha;
						if(first_time)
						{
							//broi2 = broi*broi;
							alph1 = logf(broi);
							bro = 1/broi;
							float alph2 = ko*(broi+1)*bro*bro;
							alpha = alph1+alph2;
						}
						float aux,rejf3,rnno19;
						int rnno19_loops=0;
						do
						{
							float rnno15=curand_uniform(&randStat);
							float rnno16=curand_uniform(&randStat);
							if((rnno15*alpha)<alph1)    //"Use 1/br part"
							{
								br = expf(alph1*rnno16)*bro;
							}
							else						//"Use the br part."
							{
								//br = sqrtf(rnno16*broi2 + (1-rnno16))*bro;
								br = sqrtf(rnno16*broi*broi + (1-rnno16))*bro;
							}
							temp = (1-br)/(ko*br);
							sinthe=max(0.0f,temp*(2-temp));
							aux = 1+br*br;
							rejf3 = aux - br*sinthe;
							rnno19=curand_uniform(&randStat);
							rnno19_loops++;
							if(rnno19_loops>100)
							{
								//printf("rnno19_loops>100! break\n");
								break;
							}
						}while(rnno19*aux>rejf3);
					}
					else            //At low energies it is faster to sample br uniformely
					{
						float rejmax;
						if(first_time)
						{
							bro = 1.0/broi;
							//bro1 = 1 - bro;
							rejmax = broi + bro;
						}
						float rnno16,rejf3;
						int rnno16_loops=0;
						do
						{
							float rnno15=curand_uniform(&randStat);
							rnno16=curand_uniform(&randStat);
							//br = bro + bro1*rnno15;
							br = bro + (1 - bro)*rnno15;
							temp = (1-br)/(ko*br);
							sinthe=max(0.0f,temp*(2-temp));
							rejf3 = 1 + br*br - br*sinthe;
							if(rnno16_loops++>100)
							{
								//printf("rnno16_loops>100! break\n");
								break;
							}
						}
						while(rnno16*br*rejmax>rejf3);
					}
					first_time=false;
					//if((br<bro)||(br>1))
					//{
						//if((br<0.99999/broi)||(br>1.00001))
						//{
							/* write(i_log,'(/a)') '***************** Warning: '
							write(i_log,*) ' sampled br outside of allowed range! ',ko,1.
							broi,br*/
						//}
						//goto RESAMPLE;
					//}
				}while(br<bro||br>1);

				costhe = 1 - temp;
				if(ibcmpl==0)    //User wants to use Klein-Nishina, so we are done
				{
					uj=0;
					//goto FINISHED_COMPTON_SAMPLING;  //change to if-else structure
				}
				else             //Check for rejection due to the limited range of pzmax
				{
					//br2 = br*br;
					float aux = ko*(ko-uj)*temp;
					float aux1 = 2*aux + uj*uj;
					float pzmax = aux - uj;
					if((pzmax<0)&&(pzmax*pzmax>=aux1))
					{
						if(ibcmpl==1)
						{		//reject this compton scatter, put the photon back to photon_step stack
										//Tong added these 2 sentences, July 2014
									p.status=p_photon_step;
									push_stack(p,p_photon_step);
									break;
						} else continue; //resample_compton loop
					}
					pzmax = pzmax/sqrtf(aux1);
					float qc2 = 1 + br*br - 2*br*costhe;
					float qc = sqrtf(qc2);

					float af,fmax1,fpz;
					if(pzmax>1)
					{
						pzmax = 1;
						af = 0;
						fmax1 = 1;
						fpz = 1;
						//goto RETRY_PZ;
					}
					else
					{
						float aux3 = 1 + 2*jo*abs(pzmax);
						float aux4 = 0.5*(1-aux3*aux3);
						fpz = 0.5*expf(aux4);
						af = qc*(1+br*(br-costhe)/qc2);

						if(af<0) 
						{
							if(pzmax>0)
							{
								fpz=1.0f-fpz;
							}
							float eta_incoh=curand_uniform(&randStat);
							if(eta_incoh>fpz)
							{
								if(ibcmpl==1)
								{		//reject this compton scatter, put the photon back to photon_step stack
										//Tong added these 2 sentences, July 2014
									p.status=p_photon_step;
									push_stack(p,p_photon_step);
									break;
								}
								else continue; //resample 
							}
							af = 0;
							fmax1 = 1;
							//goto RETRY_PZ;  //change to if-else structure
						}
						else
						{
							float fpz1;
							if(pzmax<-0.15)
							{
								fmax1 = 1-af*0.15;
								fpz1 = fpz*fmax1*jo;
							}
							else if(pzmax<0.15)
							{
								fmax1 = 1 + af*pzmax;
								aux3 = 1/(1+0.33267252734*aux3);  //0.33267252734 is p/sqrt(2), p is the parameter from Eq. 7.1.25
								                                  //of Abramowitz and Stegun, needed for approximating Erf
								aux4 = fpz*aux3*(0.3480242+aux3*(-0.0958798+aux3*0.7478556)) + erfjo_array[j-1];
								if(pzmax>0)
								{
									fpz1 = (1 - fmax1*fpz)*jo - 0.62665706866*af*aux4;  //0.62665706866 is sqrt(Pi/8)
									fpz = 1 - fpz;
								}
								else
								{
									fpz1 = fmax1*fpz*jo - 0.62665706866*af*aux4;
								}
							}
							else 
							{
								fmax1= 1 + af*0.15;
								fpz1 = (1 - fmax1*fpz)*jo;
								fpz = 1 - fpz;
							}
							float eta_incoh=curand_uniform(&randStat);
							if(eta_incoh*jo>fpz1)
							{
								if(ibcmpl==1)
								{		//reject this compton scatter, put the photon back to photon_step stack
										//Tong added these 2 sentences, July 2014
									p.status=p_photon_step;
									push_stack(p,p_photon_step);
									break;
								}
								else continue;//resample_compton loop
							}
						}  //else of if(af<0)
					} //if(pzmax>1)

					//At this point, all rejections are handled, now we need to sample pz
					//between -1 and pzmax using the Compton profile of the selected shell
					//and F(pz,cos(theta)) as a rejection function
//RETRY_PZ:
					float pz;
//					for(;;)
					for(int pz_loops=0;;pz_loops++) 
					{
						if(pz_loops>100){
							//printf("pz_loops>100! break!\n");
							break;
						}
						if(ibcmpl!=2)
						{
							float rnno18=curand_uniform(&randStat);
							rnno18 = rnno18*fpz;
							if(rnno18<0.5)
							{
								rnno18=max(1e-30,2*rnno18);
								pz = 0.5*(1-sqrt(1-2*logf(rnno18)))/jo;
							}
							else
							{
								rnno18 = 2*(1-rnno18);
								pz = 0.5*(sqrt(1-2*logf(rnno18))-1)/jo;
							}
							if(abs(pz)>1)
							{
								//Due to the non-relativistic approximation for pz, it has to be between -1 and 1
								//goto RETRY_PZ;
								continue;  //RETRY_PZ loop
							}
							if(pz<0.15)
							{
								float frej;
								if(pz<-0.15) frej=1-af*0.15;
								else         frej=1+af*pz;
								float eta=curand_uniform(&randStat);
								if(eta*fmax1>frej)
								{
									//goto RETRY_PZ;
									continue;  //RETRY_PZ loop
								}
							}
							//If pz > 0.15, F is always 1 => no need for rejection
						}
						else 
						{
							pz=0;    //no Doppler broadenning and no binding energy
							uj=0;
						}
						break;  //make sure to exit RETRY_PZ loop
					}  //for(pz_loop)

					//Calculate energy of scattered photon
					//pz2 = pz*pz;
					if(abs(pz)<0.01)
					{
						//br = br*(1 + pz*(qc + (br2-costhe)*pz));
						br = br*(1 + pz*(qc + (br*br-costhe)*pz));
					}
					else
					{
						//aux = 1 - pz2*br*costhe;
						aux = 1 - pz*pz*br*costhe;
						//aux1 = 1 - pz2*br2;
						aux1 = 1 - pz*pz*br*br;
						//aux2 = qc2 - br2*pz2*sinthe;
						float aux2 = qc2 - br*br*pz*pz*sinthe;
						if(aux2>1e-10)
						{
							br = br/aux1*(aux+pz*sqrt(aux2));
						}
					}
					uj = uj*prm;
				}

//FINISHED_COMPTON_SAMPLING:
				float pesg = br*peig;
				float pese = peig - pesg - uj + prm;
				sinthe = sqrtf(sinthe);
				p.e = pesg;
				p.status = p_photon_step;
				p.nRepeat = 0;    //no repeat for secondary particle
				float sinphi,cosphi;
				uphi21(costhe,sinthe,cosphi,sinphi,p,&randStat);
				if(p.e<reg_dat.pcut)  //if the photon energy is lower than the pcut, deposit the energy locally.
					edep += p.e;
				else
				{
#ifdef NUM_MLC_COMPTON
					if(p.region==2)p.latch+=1;
#endif
						push_stack(p,p_photon_step);
				}

				float aux = 1 + br*br - 2*br*costhe;
				if(aux>1e-8)
				{
					costhe = (1-br*costhe)/sqrt(aux);
					sinthe = (1-costhe)*(1+costhe);
					if(sinthe>0) sinthe=-sqrtf(sinthe);
					else sinthe=0;
				}
				else 
				{
					costhe = 0.0f;
					sinthe = -1.0f;
				}
				p1.e = pese;
				p1.charge = -1;
				p1.status = e_electron_step;
				p1.nRepeat = 0;  //no repeat for secondary particle
				uphi32(costhe,sinthe,cosphi,sinphi,p1);
				if(p1.e<reg_dat.ecut)
					edep += p1.e - prm;
				else
				{
#ifdef NUM_MLC_COMPTON
					if(p.region!=2)
#endif
					   push_stack(p1,e_electron_step);  //push the electron to stack, however if it inside MLC, we aboundent the electron
				}

				if((ibcmpl==1)||(ibcmpl==3))
				{
					if(uj>RELAX_CUTOFF)
					{
						relax(uj,shn_array[j-1],iz_array[j-1],edep,p,reg_dat,&randStat);
						//relax will put all particles with energies above ecut, pcut on the
						//stack, the remaining energy will be added to edep and deposited
					}
					else 
					{
						edep += uj;
					}
				}
				break;
			} //for(;;) loop
			
			escore(p.region,edep*p.wt);
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
