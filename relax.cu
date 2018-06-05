/*=================================================================
  Subroutine to fill a vacancy in shell n, element iZ
  by emitting fluorescent X-rays, Auger and Coster-Kronig electrons
  Transitions between K,L1,L2,L3,average M,average N are taken into
  account. Particles with energies above the transport cut-offs
  (ECUT and PCUT) are placed on the stack, energy of sub-threshold
  particles is stored in EDEP.
  In this version a global cut-off of 1 keV applies
  i.e. if ECUT-RM or PCUT is below 1 keV, binding energies below
  1 keV will still be absorbed locally (due to lack of data)
  =================================================================*/


__device__ void relax (float energy, short n, short iz, float &edep, particle_t p, region_data_t &reg_dat, curandState *randStat)
{
	short vac_array[MXVAC];
	//short n_vac,shell;
	//short final,final1,final2,iql;
	//short k;
	float e_array[MXVAC];
	//float ei,ef,ex,eta,e_cut,ekcut,pkcut,elcut;
	//float xphi,yphi,xphi2,yphi2,rhophi2, cphi,sphi;

	if((n<1)||(n>MXSHELL)) return;        //unknown vacancy

	iz_relax=iz;
	//ekcut=reg_dat.ecut-rm;
	//pkcut=reg_dat.pcut;
	//e_cut=min(ekcut,pkcut);
	float e_cut=min((reg_dat.ecut-rm),reg_dat.pcut);
	e_cut=max(RELAX_CUTOFF,e_cut);
	if (energy<=e_cut)
	{
		edep += energy;        //We assume that edep is zeroed or set to the appropriate value in the routine calling RELAX
		return;
	}

	// Set-up the array of vacancies for the relaxation cascade
	char n_vac=1;
	vac_array[n_vac-1]=n;
	e_array[n_vac-1]=energy;

//START:
	for(;;)
	{                            //Until no >N-shell vacancies
		short shell=vac_array[n_vac-1];
		float ei=e_array[n_vac-1];
		n_vac-=1;
		if(ei<=e_cut)                   //Below cut-off -> local absorption
		{
			edep += ei;
			if(n_vac>0) 
			{
				continue;
			}
			break;
		}
		//Set the relax_user common block variables
		ish_relax=shell;
		u_relax=ei;
		if(shell==MXSHELL)              //This is N-shell vacancy -> just produce Auger
		{ 
			if(ei>e_cut) 
			{
				p.e = ei + prm;
				p.charge = -1;
				p.status = e_electron_step;
				p.nRepeat = 0;  //no repeat for secondary particle

				float eta=curand_uniform(randStat);
				eta=2*eta-1;
				p.w=eta;
				eta=(1-eta)*(1+eta);
				if(eta>SMALL_POLAR_ANGLE_THRESHOLD)
				{
					eta=sqrtf(eta);
					float xphi,yphi,rhophi2;
					do{        // the box method to sample phi
						xphi=curand_uniform(randStat);
						xphi=2*xphi-1;
						//xphi2=xphi*xphi;
						yphi=curand_uniform(randStat);
						//yphi2=yphi*yphi;
						//rhophi2=xphi2+yphi2;
						rhophi2=xphi*xphi+yphi*yphi;
					}while (rhophi2>1); 
					rhophi2=1/rhophi2;
					//cphi=(xphi2-yphi2)*rhophi2;
					p.u=eta*(xphi*xphi-yphi*yphi)*rhophi2;
					//sphi=2*xphi*yphi*rhophi2;
					p.v=eta*2*xphi*yphi*rhophi2;

					//p.u=eta*cphi;
					//p.v=eta*sphi;
				}
				else 
				{
					p.u=0;
					p.v=0;
					p.w=1;
				}
				if(p.e<reg_dat.ecut)		 
					edep += p.e - prm;
				else
					push_stack(p,e_electron_step);
			}
			else
			{
				edep += ei;
			}
			if(n_vac>0)
			{
				continue;
			}
			break;
		}

		// Sample transition number for this vacancy

		float eta=curand_uniform(randStat);
		char k;
		for(k=first_transition[shell-1];k<=(last_transition[shell-1]-1);k++)
		{  
			eta=eta-relaxation_prob[MXTRANS*(iz-1)+k-1];
			if(eta<=0.0f) break;        //sample the transition index number k
		} 
		short final=final_state[k-1];
		//finala=final;
		char iql;
		float elcut,ex;
		if(final<100)
		{
			if(final<10)                //fuorescence
			{
				iql=0;
				//elcut=pkcut;
				elcut=reg_dat.pcut;     // set the photon energy cut as current cut
			}
			else                        //Coster-Kronig
			{
				final=final-10;
				iql=-1;
				//elcut=ekcut;
				elcut=reg_dat.ecut-rm;  // set the electron energy cut as current cut.
			}
			float ef=binding_energies[MXSHELL*(iz-1)+final-1];
			ex=ei-ef;
			n_vac+=1;
			vac_array[n_vac-1]=final;
			e_array[n_vac-1]=ef;
		}
		else                            //Auger will creat two new vac 
		{
			short final1=final/100;
			short final2=final-final1*100;
			n_vac=n_vac+1;
			vac_array[n_vac-1]=final1;
			e_array[n_vac-1]= binding_energies[MXSHELL*(iz-1)+final1-1];
			n_vac=n_vac+1;
			vac_array[n_vac-1]=final2;
			e_array[n_vac-1]= binding_energies[MXSHELL*(iz-1)+final2-1];
			iql=-1;
			ex=ei-e_array[n_vac-1]-e_array[n_vac-2];
			//elcut=ekcut;
			elcut=reg_dat.ecut-rm;
		}
		if(ex<=elcut)                   //Below cut-off
		{
			edep += ex;
		}
		else
		{
			p.charge=iql;
			if(iql==0)
			{
				p.e = ex;
				p.status = p_photon_step;
				p.nRepeat = 0;  //no repeat for secondary particle
			}
			else
			{
				p.e = ex + rm;
				p.status = e_electron_step;
				p.nRepeat = 0;  //no repeat for secondary particle
			}

			eta=curand_uniform(randStat);
			eta=2*eta-1;
			p.w=eta;   //cos(theta)
			eta=(1-eta)*(1+eta);  
			if(eta>SMALL_POLAR_ANGLE_THRESHOLD)
			{
				eta=sqrtf(eta);
				float xphi,yphi,rhophi2;
				do{
					xphi=curand_uniform(randStat);
					xphi=2*xphi-1;
					//xphi2=xphi*xphi;
					yphi=curand_uniform(randStat);
					//yphi2=yphi*yphi;
					//rhophi2=xphi2+yphi2;
					rhophi2=xphi*xphi+yphi*yphi;
				}while(rhophi2>1);
				rhophi2=1/rhophi2;
				//cphi=(xphi2-yphi2)*rhophi2;
				p.u=eta*(xphi*xphi-yphi*yphi)*rhophi2;
				//sphi=2*xphi*yphi*rhophi2;
				p.v=eta*2*xphi*yphi*rhophi2;

				//p.u=eta*cphi;
				//p.v=eta*sphi;
			}
			else 
			{
				p.u=0.0f;
				p.v=0.0f;
				p.w=1.0f;
			}
			if(iql==0)
			{
				if(p.e<reg_dat.pcut)
					edep += p.e;
				else
					push_stack(p,p_photon_step);
			}
			else
			{
				if(p.e<reg_dat.ecut)
					edep += p.e - prm;
				else
					push_stack(p,e_electron_step);
			}
		}
	}

	//escore(pregion,edep);  //we nologer depose energy in relax, the edep is returned back to the up-stream process  to deposit
	return;
} 
