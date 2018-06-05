/*==============================================================
  For a photon energy below 2.1 MeV, the energies of the pair 
  particles are uniformly distributed in the allowed range via 
  the default replacement for $SELECT-LOW-ENERGY-PAIR-PRODICTION;
  If the user has a better approach, modify this macro. 
  For a photon energy between 2.1 and 50 MeV the Bethe-Heitler 
  cross section is employed, above 50 MeV the Coulomb-corrected 
  Bethe-Heitler is used. 
  Modified from its original version to make compatible with the 
  changes made in BREMS.
  ==============================================================*/


__global__ void pair_production()
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

		if(process && !pop_stack(p,p_pair))
			p.status = p_empty;

        if(simu_stack_is_empty){
            process = (p.status == p_empty);
            if(__all(process))
                break;
        }

		process = (p.status == p_pair);

#ifdef DO_STEP_COUNT
        uint count_mask = __ballot(process);
		if (idx.t == 0) step_counters[p_pair] += __popc(count_mask);
#endif

		if(process){

			region_data_t reg_dat = region_data[p.region];

			//float peig,pese1,pese2;
			//float eig,ese1,ese2;
			//float rnno30,rnno31,rnno32,rnno33,rnno34;
			//float delta,rejf,rejmax,aux1,aux2,amax,bmax,del0,br,eminus,eplus,eavail;
			//int   l,l1;
			//float ese,pse,ztarg,tteig,ttese,esedei,eseder,ximin,ximid,rejmin,rejmid,rejtop,ya,xitry;
			//float galpha,gbeta,xitst,rejtst_on_rejtop,rejtst,rtest;
			//int   ichrg;
			//int   iq1,iq2,iprdst_use;
			//bool  do_nrc_pair;
			//float b1;
			//float theta;
			//float costhe,sinthe;
			float edep=0.0f;

			particle_t p1=p;
			ushort medium=reg_dat.med;

			float eig=p.e;
			bool do_nrc_pair=false;

			float ese1,ese2;
			char iq1,iq2;
			if(!do_nrc_pair)
			{
				if(eig<=2.1)                //BELOW 2.1,USE APPROXIMATION
				{
					float rnno30=curand_uniform(&randStat);
					ese2=prm+0.5*rnno30*(eig-2*prm);
					//rnno34=curand_uniform(randStat);
					rnno30=curand_uniform(&randStat);
					ese1=eig-ese2;
					//if(rnno34<0.5)
					if(rnno30<0.5)
					{
						iq1=-1;
						iq2=1;
					}
					else
					{
						iq1=1;
						iq2=-1;
					}  
					//uniform energy distribution is probably a better approximation than
					//a zero energy 'electron' for low energy pair production
				}
				else                        //Above 2.1, must sample
				{                           //Decide whether to use Bethe-Heitler or BH Coulomb corrected
					float amax,bmax,aux1;
					char l,l1;
					if(eig<50.0f)           //Use BH without Coulomb correction
					{
						l=5;
						l1=l+1;
						float delta=4*delcm[medium]/eig;
						if(delta<1)
						{
							amax=dl1[8*medium+l-1]+delta*(dl2[8*medium+l-1]+delta*dl3[8*medium+l-1]);
							bmax=dl1[8*medium+l1-1]+delta*(dl2[8*medium+l1-1]+delta*dl3[medium*8+l1-1]);
						}
						else
						{
							float aux2=logf(delta+dl6[medium*8+l-1]);
							amax=dl4[medium*8+l-1]+dl5[medium*8+l-1]*aux2;
							bmax=dl4[medium*8+l1-1]+dl5[medium*8+l1-1]*aux2;
						}
						aux1=1-rmt2/eig;
						aux1=aux1*aux1;
						aux1=aux1*amax/3;
						aux1=aux1/(bmax+aux1);
					}
					else                    //Use BH Coulomb-corrected
					{
						l=7;
						amax=dl1[medium*8+l-1];
						bmax=dl1[medium*8+l];
						aux1=bpar[medium*2+1]*(1-bpar[medium*8]*rm/eig);
					}
					float del0=eig*delcm[medium];
					float eavail=eig-rmt2;

					float rnno30,rejmax,rejf;
					float eminus;
					do{
						rnno30=curand_uniform(&randStat);
						float rnno31=curand_uniform(&randStat);
						float br;
						if(rnno30>aux1)
						{
							br=0.5*rnno31;
							rejmax=bmax;
							l1=l+1;
						}
						else
						{
							//rnno32=curand_uniform(randStat);
							//rnno33=curand_uniform(randStat);
							//b1=max(max(rnno31,rnno32),rnno33);
							rnno30=curand_uniform(&randStat);
							rnno30=max(rnno30,rnno31);
							rnno31=curand_uniform(&randStat);
							float b1=max(rnno30,rnno31);
							br=0.5*(1-b1);
							rejmax=amax;
							l1=l;
						}
						eminus=br*eavail+rm;
						float eplus=eig-eminus;
						float delta=del0/(eminus*eplus);
						if(delta<1)
						{
							rejf=dl1[medium*8+l1-1]+delta*(dl2[medium*8+l1-1]+delta*dl3[medium*8+l1-1]);
						}
						else
						{
							rejf=dl4[medium*8+l1-1]+dl5[medium*8+l1-1]*logf(delta+dl6[medium*8+l1-1]);
						}
						//rnno34=curand_uniform(randStat);
						rnno30=curand_uniform(&randStat);
						//}while((rnno34*rejmax)>rejf);
					}while((rnno30*rejmax)>rejf);

					ese2=eminus;
					ese1=eig-ese2;
					//rnno34=curand_uniform(randStat);
					rnno30=curand_uniform(&randStat);
					//if(rnno34<0.5)
					if(rnno30<0.5)
					{
						iq1=-1;
						iq2=1;
					}
					else
					{
						iq1=1;
						iq2=-1;
					}
				}
			}

			// ENERGY GOING TO LOWER SECONDARY HAS NOW BEEN DETERMINED

			//ese2=pese2;
			p.e = ese1;
			p1.e = ese2;
			p.charge = iq1;
			p1.charge = iq2;
			p.status = e_electron_step;
			p1.status = e_electron_step;
			p.nRepeat = 0;    //no repeat for secondary particle
			p1.nRepeat = 0;    //no repeat for secondary particle
			//This average angle of emission for both pair production
			//and bremsstrahlung is much smaller than the average angle
			//of multiple scattering for delta T transport=0.01 R.L.
			//The initial and final momenta are coplanar

			//Set up a new electron
			//Select the angle from the leading term of the angular distribution
			//usage: iprdst=0 => EGS4 default angle selection
			//       iprdst=1 => lowest order angular distribution
			//       iprdst=2 => Motz, Olsen and Koch (1969) eq. 3D-2003
			//                   If iprdst is non-zero and e_photon < $BHPAIR
			//                   The iprdst=1 distribution is used
			float theta,sinthe,costhe;
			if(iprdst>0)
			{
				char iprdst_use;
				if(iprdst==4)
				{
					float rtest=curand_uniform(&randStat);
					float gbeta=ese1/(ese1+10);
					if(rtest<gbeta)
					{
						iprdst_use=1;
					}
					else
					{
						iprdst_use=4;
					}
				}
				else if((iprdst==2)&&(eig<4.14))
				{
					iprdst_use=1;
				}
				else
				{
					iprdst_use=iprdst;
				}
				for(char ichrg=1;ichrg<=2;ichrg++)
				{
					float ese;
					if(ichrg==1)
					{
						ese=ese1;
					}
					else
					{
						ese=ese2;
						if(iprdst==4)
						{
							float gbeta=ese/(ese+10);
							float rtest=curand_uniform(&randStat);
							if(rtest<gbeta)
							{
								iprdst_use=1;
							}
							else
							{
								iprdst_use=4;
							}
						}
					}
					if(iprdst_use==1)
					{
						float pse=sqrtf(max(0.0f,(ese-rm)*(ese+rm)));
						costhe=curand_uniform(&randStat);
						costhe=1.0f-2.0*costhe;
						sinthe=rm*sqrtf((1.0f-costhe)*(1.0f+costhe))/(pse*costhe+ese);
						costhe=(ese*costhe+pse)/(pse*costhe+ese);
					}
					else if(iprdst_use==2)
					{
						float ztarg=zbrang[medium];
						float tteig=eig/rm;
						float ttese=ese/rm;
						//ttpse=sqrtf((ttese-1.0f)*(ttese+1.0f));
						float esedei=ttese/(tteig-ttese);
						float eseder=1.0f/esedei;
						float ximin=1.0f/(1.0f+(3.141593*ttese)*(3.141593*ttese));
						float rejmin=2.0f+3.0f*(esedei+eseder)-4.00*(esedei+eseder+1.0f-4.0f*(ximin-0.5f)*(ximin-0.5f))*
							(1.0f+0.25f*logf(((1.0f+eseder)*(1.0f+esedei)/(2.0f*tteig))*((1.0f+eseder)*(1.0f+esedei)/(2.0f*tteig))+ztarg*ximin*ximin));
						float ya=(2.0f/tteig)*(2.0f/tteig);
						float xitry=max(0.01f,max(ximin,min(0.5f,sqrtf(ya/ztarg))));
						float galpha=1.0f+0.25*logf(ya+ztarg*xitry*xitry);
						float gbeta=0.5f*ztarg*xitry/(ya+ztarg*xitry*xitry);
						galpha=galpha-gbeta*(xitry-0.5f);
						float ximid=galpha/(3.0f*gbeta);
						if(galpha>=0.0f)
						{
							ximid=0.5f-ximid+sqrtf(ximid*ximid+0.25);
						}
						else
						{
							ximid=0.5f-ximid-sqrtf(ximid*ximid+0.25);
						}
						ximid=max(0.01f,max(ximin,min(0.5f,ximid)));
						float rejmid=2.0f+3.0f*(esedei+eseder)-4.0f*(esedei+eseder+1.0f-4.0*(ximid-0.5f)*(ximid-0.5f))*
							(1.0f+0.25*logf(((1.0f+eseder)*(1.0f+esedei)/(2.0f*tteig))*((1.0f+eseder)*(1.0f+esedei)/(2.0f*tteig))+ztarg*ximid*ximid));
						float rejtop=1.02*max(rejmin,rejmid);

						float rtest,rejtst_on_rejtop;
						do{
							float xitst=curand_uniform(&randStat);
							float rejtst=2.0f+3.0f*(esedei+eseder)-4.0f*(esedei+eseder+1.0f-4.0f*(xitst-0.5)*(xitst-0.5))*
								(1.0f+0.25*logf(((1.0f+eseder)*(1.0f+esedei)/(2.0f*tteig))*((1.0f+eseder)*(1.0f+esedei)/(2.0f*tteig))+ztarg*xitst*xitst));
							rtest=curand_uniform(&randStat);
							theta=sqrtf(1.0f/xitst-1.0f)/ttese;
							rejtst_on_rejtop=rejtst/rejtop;
						}while((rtest>rejtst_on_rejtop)||(theta>=PI));

						__sincosf(theta,&sinthe,&costhe);
						//sinthe=sin(theta);
						//costhe=cos(theta);
					}
					else if(iprdst_use==3)
					{
						costhe=curand_uniform(&randStat);
						costhe=1.0f-2.0f*costhe;
						sinthe=(1-costhe)*(1+costhe);
						if(sinthe>0)
						{
							sinthe=sqrtf(sinthe);
						}
						else
						{
							sinthe=0;
						}
					} 
					else 
					{
						costhe=curand_uniform(&randStat);
						costhe=1-2*sqrtf(costhe);
						sinthe=(1-costhe)*(1+costhe);
						if(sinthe>0.0f)
						{
							sinthe=sqrtf(sinthe);
						}
						else
						{
							sinthe=0;
						}
					}
					float sinphi,cosphi;
					if(ichrg==1)
					{
						uphi21(costhe,sinthe,cosphi,sinphi,p,&randStat);
						if(p.e<reg_dat.ecut)
						{
							edep += p.e - prm;
							if(p.charge>0)
								annih_at_rest(p,reg_dat,&randStat);
						}
						else
							push_stack(p,e_electron_step);
					}
					else
					{
						sinthe=-sinthe;
						uphi32(costhe,sinthe,cosphi,sinphi,p1);
						if(p1.e<reg_dat.ecut)
						{
							edep += p1.e - prm;
							if(p1.charge>0)
								annih_at_rest(p1,reg_dat,&randStat);
						}
						else
							push_stack(p1,e_electron_step);
					}
				}
				escore(p.region,edep*p.wt);
				p.status=p_empty;
				continue;
			}
			else
			{
				theta=0.0f;
			}
			float sinphi,cosphi;
			uphi11(theta,costhe,sinthe,cosphi,sinphi,p,&randStat);
			if (p.e<reg_dat.ecut)
			{
				edep += p.e - prm;
				if(p.charge>0)
					annih_at_rest(p,reg_dat,&randStat);
			}
			else
				push_stack(p,e_electron_step);

			sinthe=-sinthe;
			uphi32(costhe,sinthe,cosphi,sinphi,p1);
			if(p1.e<reg_dat.ecut)
			{
				edep += p1.e - prm;
				if(p1.charge>0)
					annih_at_rest(p1,reg_dat,&randStat);
			}
			else
				push_stack(p1,e_electron_step);

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
