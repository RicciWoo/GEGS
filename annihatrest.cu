/*============================================================
  It is handy to be able to initiate annihilation at rest from 
  places other than the electron discard section (e.g. AUSGAB) 
  Annihilation at rest takes a sufficent amount of time to not 
  have any real benefit from this code being inline in the 
  ELECTR subroutine.
  ============================================================*/


__device__ void annih_at_rest(particle_t &p, region_data_t &reg_dat, curandState *randStat)
{
	//float costhe,sinthe,cphi,sphi;
	//int ibr;
	//float xphi,xphi2,yphi,yphi2,rhophi2;
	float edep=0.0f;

	if(nbr_split>1)
	{
		p.wt=p.wt/nbr_split;
	}
	p.e=prm;
	p.charge=0;
	p.status=p_photon_step;
	p.nRepeat = 0;  //no repeat for secondary particle

	for(char ibr=1;ibr<=nbr_split;ibr++)
	{
		float costhe=curand_uniform(randStat);
		costhe=2*costhe-1;
		float sinthe=sqrtf(max(0.0,(1-costhe)*(1+costhe)));

		float xphi,yphi,rhophi2;
		do{	
			xphi=curand_uniform(randStat);
			xphi = 2*xphi - 1;
			//xphi2 = xphi*xphi;
			yphi=curand_uniform(randStat);
			//yphi2 = yphi*yphi;
			//rhophi2 = xphi2 + yphi2;
			rhophi2 = xphi*xphi + yphi*yphi;
		}while(rhophi2>1.0f);
		rhophi2 = 1/rhophi2;
		//cphi = (xphi2 - yphi2)*rhophi2;
		float cphi = (xphi*xphi - yphi*yphi)*rhophi2;
		float sphi = 2*xphi*yphi*rhophi2;

		p.u=sinthe*cphi;
		p.v=sinthe*sphi;
		p.w=costhe;
		if(p.e<reg_dat.pcut)
			edep += p.e;
		else
			push_stack(p,p_photon_step);

		p.u=-p.u;
		p.v=-p.v;
		p.w=-p.w;
		if(p.e<reg_dat.pcut)
			edep += p.e;
		else
			push_stack(p,p_photon_step);
	}

	escore(p.region,edep*p.wt);
	p.status = e_empty;
	return;
}

