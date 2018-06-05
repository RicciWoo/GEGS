/*==============================================================
  Samples bremsstrahlung energy using "
   - Coulomb corrected Bethe-Heitler above 50 MeV
   - Bethe-Heitler below 50 MeV
  if ibr_nist=0, or 
   - the NIST bremsstrahlung cross section data base 
     (prepared in a form of an alias table for rapid sampling) 
  if ibr_nist=1 or 
   - the NRC bremsstrahlung cross section data base, which is 
     the same as the NIST database, but with corrections to 
     the electron-electron contribution, which are mostly 
     important for low Z and low k 
  if ibr_nist=2 
  and direction using
   - formula 2BS from from Koch and Motz if IBRDST=1 
   - leading term of the brems angular dsstr. if IBRDST=0 
   - photon direction=electron direction if IBRDST<0 
  This version replaces the original EGS4 implementation 
  because of a bug discovered in the EGS4 brems routine 
  In order to work properly, the parameter DL1,..,DL6 
  are re-calculated in subroutine fix_brems which is called 
  from HATCH 
  In addition, this version has the internal capability of 
  bremsstrahlung splitting. 
  To use bremsstrahlung splitting, set nbr_split (COMON/BREMPR/) 
  to the desired number > 1 (1 is the default) 
  Be aware that event-by-event energy conservation is NOT 
  guaranteed, so don't use for calculations where this is 
  important (e.g. calculation of detector response functions) 
  The result will be nbr_split photons, all with the weight 
  wt(npold)/nbr_split, and an electron with the original weight 
  and energy given by the incident energy-energy of last photon
  ==============================================================*/


__global__ void brems()
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

		if(process && !pop_stack(p,e_brems))
			p.status = e_empty;

        if(simu_stack_is_empty){
            process = (p.status == e_empty);
            if(__all(process))
                break;
        }

		process = (p.status == e_brems);

#ifdef DO_STEP_COUNT
        uint count_mask = __ballot(process);
		if (idx.t == 0) step_counters[e_brems] += __popc(count_mask);
#endif

		if(process){

			region_data_t reg_dat = region_data[p.region];

			//float peie,pesg,pese;
			//float eie,ekin,brmin,waux,aux,ajj,rnno06,rnno07,br,esg,ese,delta,phi1,phi2,rejf;
			//float a,b,c,sinpsi,sindel,cosdel,us,vs,ztarg,tteie,beta,y2max,y2maxi,ttese;
			//float rjarg1,rjarg2,rjarg3,rejmax,rejtst,esedei,y2tst,y2tst1,rtest;
			//float xphi,yphi,xphi2,yphi2,rhophi2,sphi,cphi;
			//int   l,l1,ibr,jj;
			//float z2max,z2maxi,aux1,aux3,aux4,aux5,aux2,weight;
			//float costhe,sinthe;
			//float elke;

			float edep=0.0f;
			particle_t p1=p;
			ushort medium=reg_dat.med;

			if(nbr_split<1)
			{
				p.status = e_empty;  //to kill this electron
				continue;        //the user can turn off brems production by setting nbr_split to zero!
			}

			float eie=p.e;
			//weight=p.wt/nbr_split;

			//DECIDE WHICH DISTRIBUTION TO USE (B-H COULOMB CORRECTED IS
			//USED FROM 50 TO 20000 MEV, B-H IS USED 1.5 TO 50 MEV)
			char l,l1;
			if(eie<50.0) l=1;
			else l=3;
			l1=l+1;

			float ekin=eie-prm;
			float elke=logf(ekin);
			float brmin=ap[medium]/ekin;
			float waux=elke-log_ap[medium];
			//this saves the time consuming log evaluation
			//log_ap = log(ap[medium]) is calculated in
			//fix_brems for each medium, elke is needed
			//in electr to calculate the branching ratios
			//and therefore it must be known at this point

			//inrdst >=0 means we will sample the photon emmision
			//angle from KM-2BS (ibrdst=1) or from the leadingterm 
			//(ibrdst=0). If nbr_split > 1, we can re-usethe following
			//quantities several time
			//defult is ibrdst=1
			float a,b,c,sinpsi,sindel,cosdel;
			float ztarg,tteie,y2maxi,z2maxi;
			if(ibrdst>=0)                           
			{
				a=p.u;
				b=p.v;
				c=p.w;
				sinpsi=a*a+b*b;
				if(sinpsi>1e-20)
				{
					sinpsi=sqrtf(sinpsi);
					sindel=b/sinpsi;
					cosdel=a/sinpsi;
				}
				ztarg=zbrang[medium];
				tteie=eie/rm;
				float beta=sqrtf((tteie-1)*(tteie+1))/tteie;
				float y2max=2*beta*(1+beta)*tteie*tteie;
				y2maxi=1/y2max;
				if(ibrdst==1)
				{
					float z2max=y2max+1;
					z2maxi=sqrtf(z2max);
				}
			}

			float ese;
			for (char ibr=1;ibr<=nbr_split;ibr++)  //User wants to use Bethe-Heitler
			{
				float rnno06,rejf;
				float esg;
				do{
					rnno06=curand_uniform(&randStat);
					//rnno07=curand_uniform(randStat);
					float br=brmin*expf(rnno06*waux);
					esg=ekin*br;
					ese=eie-esg;
					float delta=esg/eie/ese*delcm[medium];
					float aux=ese/eie;
					float phi1,phi2;
					if(delta<1)
					{
						phi1=dl1[medium*8+l-1]+delta*(dl2[medium*8+l-1]+delta*dl3[medium*8+l-1]);
						phi2=dl1[medium*8+l1-1]+delta*(dl2[medium*8+l1-1]+delta*dl3[medium*8+l1-1]);
					}
					else 
					{
						phi1=dl4[medium*8+l-1]+dl5[medium*8+l-1]*logf(delta+dl6[medium*8+l-1]);
						phi2=phi1;
					}
					rejf=(1+aux*aux)*phi1-2*aux*phi2/3;
					rnno06=curand_uniform(&randStat);
					//}while(rnno07>=rejf);
				}while(rnno06>=rejf);

				// SET UP THE NEW PHOTON   //The photon will inherit the direction from the electron.
				p1.e = esg;
				p1.charge = 0;
				//p1.wt = weight;
				p1.wt = p.wt/nbr_split;
				p1.status = p_photon_step;
				p1.nRepeat = 0;  //no repeat for secondary particle

				if(ibrdst>=0)
				{
					float y2tst;
					if(ibrdst==1)
					{
						//ttese=ese/rm;
						//esedei=ttese/tteie;
						float esedei=ese/rm/tteie;
						float rjarg1=1+esedei*esedei;
						float rjarg2=rjarg1+2*esedei;
						float aux=2*ese*tteie/esg;
						aux=aux*aux;
						float aux1=aux*ztarg;
						float rjarg3;
						if(aux1>10)
						{
							rjarg3=lzbrang[medium]+(1-aux1)/(aux1*aux1);
						}
						else 
						{
							rjarg3=logf(aux/(1+aux1));
						}
						float rejmax=rjarg1*rjarg3-rjarg2;

						float rtest,rejtst;
						do{
							y2tst=curand_uniform(&randStat);
							rtest=curand_uniform(&randStat);
							float aux3=z2maxi/(y2tst+(1-y2tst)*z2maxi);
							rtest=rtest*aux3*rejmax;
							y2tst=aux3*aux3-1;
							float y2tst1=esedei*y2tst/(aux3*aux3*aux3*aux3);
							float aux4=16*y2tst1-rjarg2;
							float aux5=rjarg1-4*y2tst1;
							if(rtest<(aux4+aux5*rjarg3))
							{
								break ;
							}
							float aux2=logf(aux/(1+aux1/(aux3*aux3*aux3*aux3)));
							rejtst=aux4+aux5*aux2;
						}while(rtest>=rejtst);
					}
					else 
					{
						y2tst=curand_uniform(&randStat);
						y2tst=y2tst/(1-y2tst+y2maxi);
					}
					float costhe=1-2*y2tst*y2maxi;
					float sinthe=sqrtf(max((1-costhe)*(1+costhe),0.0f));

					float xphi,yphi,rhophi2;
					do{
						xphi=curand_uniform(&randStat);
						xphi=2*xphi-1;
						//xphi2=xphi*xphi;
						yphi=curand_uniform(&randStat);
						//yphi2=yphi*yphi;
						//rhophi2=xphi2+yphi2;
						rhophi2=xphi*xphi+yphi*yphi;
					} while(rhophi2>1);
					rhophi2=1/rhophi2;
					//cphi=(xphi2-yphi2)*rhophi2;
					float cphi=(xphi*xphi-yphi*yphi)*rhophi2;
					float sphi=2*xphi*yphi*rhophi2;

					if(sinpsi>=1e-10)
					{
						float us=sinthe*cphi ;
						float vs=sinthe*sphi;
						p1.u = c*cosdel*us - sindel*vs + a*costhe;
						p1.v = c*sindel*us + cosdel*vs + b*costhe;
						p1.w = c*costhe - sinpsi*us;
					}
					else
					{
						p1.u = sinthe*cphi;
						p1.v = sinthe*sphi;
						p1.w = c*costhe;
					}
				}
				if(p1.e<reg_dat.pcut)
					edep += p1.e;
				else
					push_stack(p1,p_photon_step);
			}

			p.e = ese;
			p.status = e_electron_step;
			p.nRepeat = 0;  //no repeat for secondary particle

			if(p.e<reg_dat.ecut)
				edep += p.e - prm;
			else
				push_stack(p,e_electron_step);

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
