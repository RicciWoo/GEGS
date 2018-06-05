__global__ void photo()
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

		if(process && !pop_stack(p,p_photo))
			p.status = p_empty;

        if(simu_stack_is_empty){
            process = (p.status == p_empty);
            if(__all(process))
                break;
        }

		process = (p.status == p_photo);

#ifdef DO_STEP_COUNT
        uint count_mask = __ballot(process);
		if (idx.t == 0) step_counters[p_photo] += __popc(count_mask);
#endif

		if(process){

			region_data_t reg_dat = region_data[p.region];

			//float eelec,beta,gamma,alpha,ratio,rnpht,fkappa,xi,sinth2,rnpht2;
			//float peig;
			float probs[MXEL];
			//float br, sigma,  aux,aux1,  sigtot,  e_vac;
			short ints[MXEL];
			//short iz,j,k;
			//bool  do_relax;
			//float costhe,sinthe;
			//float gle;
			//short iedgfl,iphter;
			float edep=0.0f;

			ushort medium=reg_dat.med;

			float peig=p.e;
			ushort iedgfl = reg_dat.flags&f_atomic_relaxation;
			ushort iphter = reg_dat.flags&f_photo_electron_angular_distribution;
			if(peig<edge_energies[1])  //edge_energies[MXEDGE*(1-1)+2-1]  
			{
				p.charge = -1;
				p.e = peig + prm;
				p.status = e_electron_step;  //this is an electron
				p.nRepeat = 0;    //no repeat for secondary particle
				if(p.e<reg_dat.ecut){
					edep = p.e - prm;
					escore(p.region,edep*p.wt);
				}
				else
					push_stack(p,e_electron_step);

				p.status=p_empty;
				continue;
			}

			float gle=logf(p.e);
			//iz=iedgfl;
			bool do_relax = false;
			float e_vac;
			char iz,j,k;
			if(iedgfl!=0)           // User requested atomic relaxations "                              
			{					    // first sample the element
				if(nne[medium]==1)  //if there is only one element in current material
				{
					iz=int(zelem[medium]+0.5);
					for(j=0;j<edge_number[iz-1];j++)
					{
						if(peig>=edge_energies[MXEDGE*(iz-1)+j])  //find the first edge with the edge-energy lower than the energy of the incident photon, so the photo-elec can happen.
						{
							break;                       // the j carries the index of such edge. 
						}
					}
				}
				else
				{
					float aux = peig*peig;
					float aux1 = aux*peig;
					aux = aux*sqrtf(peig);
					float sigtot = 0;
					for(k=0;k<nne[medium];k++)
					{
						iz=int(zelem[MXMED*k+medium]+0.5);
						//we trust our initialization so the following error should not happen.
						/*if((iz<1)||(iz>MXELEMENT))
						{
							write(i_log,*) ' Error in PHOTO: '
								write(i_log,'(/a)') '***************** Error: '
								write(i_log,*) '   Atomic number of element ',k, ' in medi
								*um ',medium,' is not between 1 and ',100
								write(i_log,'(/a)') '***************** Quiting now.'
								call exit(1)
						}*/

						float sigma;
						if(peig>edge_energies[MXEDGE*(iz-1)])   //if the incident photon energy is bigger than the first (k shell, the highest energy) edge 
						{
							j=1;								//the ph-e effect can start from the first edge.
							//caculate the photo absorbtion cross-section of element iz at current energy
							sigma=(edge_a[MXEDGE*(iz-1)]+edge_b[MXEDGE*(iz-1)]/peig+edge_c[MXEDGE*(iz-1)]/aux+edge_d[MXEDGE*(iz-1)]/aux1)/peig;
						}
						else
						{
							for(j=2;j<=edge_number[iz-1];j++)   //find the first edge with the edge-energy lower than the energy of the incident photon, so the photo-elec can happen.
							{
								if(peig>=edge_energies[MXEDGE*(iz-1)+j-1])
									break;
							}
							// note that the cross-section forumlar is different if the first possible shell is not k-shell.
							sigma=edge_a[MXEDGE*(iz-1)+j-1]+gle*(edge_b[MXEDGE*(iz-1)+j-1]+gle*(edge_c[MXEDGE*(iz-1)+j-1]+gle*edge_d[MXEDGE*(iz-1)+j-1]));
							sigma=exp(sigma);
						}
						sigma = sigma * pz[MXMED*k+medium];
						sigtot = sigtot + sigma;
						probs[k] = sigma;
						ints[k] = j;
					}
					float br=curand_uniform(&randStat);
					br = br*sigtot;
					//sample the element number index k
					for(k=0;k<nne[medium];k++)
					{
						br=br-probs[k];
						if(br<=0)break;
					} 

					iz=int(zelem[MXMED*k+medium]+0.5); //get its z number
					j=ints[k];
				}

				// Now we know the atomic number (iZ) and the energy interval the
				// photon energy is in (j). It is time to sample the shell the photon
				// is interacting with.

				if(peig<=binding_energies[MXSHELL*(iz-1)+5])  //"Below N-shell -> local energy deposition "
				{
					edep += peig;
					//p.charge is 0, so this situation will not enter the sampling of algular distribution
				}
				else
				{
					float br=curand_uniform(&randStat);
					for(k=1;k<=5;k++)
					{
						if(peig>binding_energies[MXSHELL*(iz-1)+k-1])  //sample the index number k of the shell that will be interacting
						{
							float tmp=interaction_prob[MXSHELL*(iz-1)+k-1];
							if(br<tmp) break;
							br=(br-tmp) /(1-tmp);
						}
					}
					e_vac = binding_energies[MXSHELL*(iz-1)+k-1];
					p.e = peig - e_vac + prm;  //the energy of the secondary electron
					do_relax = true; 
					p.charge = -1;
					p.status = e_electron_step;
					p.nRepeat = 0;  //no repeat for secondary particle
				}
			}
			else             //iedgfl==0, not doing relaxiation
			{
				p.e = peig + prm;      
				p.charge = -1;
				p.status = e_electron_step;
				p.nRepeat = 0;  //no repeat for secondary particle
			}

			if(p.charge==-1)
			{
				if(p.e>reg_dat.ecut)
				{
					if(iphter)        //if we consider algular distribution of the secondary electron 
					{
					//eelec=p.e;
					//if(eelec>reg_dat.ecut)
						//beta=sqrtf((eelec-rm)*(eelec+rm))/eelec;
						float beta=sqrtf((p.e-rm)*(p.e+rm))/p.e;
						//gamma=eelec/rm;
						float gamma=p.e/rm;
						float alpha=0.5*gamma-0.5+1.0f/gamma;
						float ratio=beta/alpha;
						float xi,rnpht2,sinth2;
						float costhe;
						do{
							float rnpht=curand_uniform(&randStat);
							rnpht=2.0f*rnpht-1;
							if(ratio<=0.2)
							{
								float fkappa=rnpht+0.5*ratio*(1.0f-rnpht)*(1+rnpht);
								if(gamma<100)
								{
									costhe=(beta+fkappa)/(1+beta*fkappa);
								}
								else
								{
									if(fkappa>0)
									{
										costhe=1-(1-fkappa)*(gamma-3)/(2*(1+fkappa)*(gamma-1)*(gamma-1)*(gamma-1));
									}
									else
									{
										costhe=(beta+fkappa)/(1.0f+beta*fkappa);
									}
								}
								xi = (1+beta*fkappa)*gamma*gamma;
							}
							else
							{
								xi=gamma*gamma*(1+alpha*(sqrtf(1.0f+ratio*(2.0f*rnpht+ratio))-1.0f));
								costhe=(1.0f-1.0f/xi)/beta;
							}
							sinth2=max(0.0,(1-costhe)*(1+costhe));
							rnpht2=curand_uniform(&randStat);
						}while(rnpht2>(0.5*(1.0f+gamma)*sinth2*xi/gamma));

						float sinthe=sqrtf(sinth2);
						float sinphi,cosphi;
						uphi21(costhe,sinthe,cosphi,sinphi,p,&randStat);
					} //iphter
					//if ipher==false, we don't sample the electron directions, just push into stack, i.e. follow the same direction of incident gamma.
					push_stack(p,e_electron_step);
				} 
				else //eelec< ecut , just deposit the energy locally.
				{
					edep += p.e - prm;
				}
			}
			if(do_relax) 
			{
				relax(e_vac,k,iz,edep,p,reg_dat,&randStat);
			}

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

