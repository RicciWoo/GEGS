


__device__ void eii_sample(short ish, short iz, float uj, particle_t &p, region_data_t &reg_dat, curandState *randStat, float &edep)
{
	//float t,tau,tau1,tau12,tau2,p2,beta2,c1,c2,wmax,xmax,fm_s,fm_h,prob_s,prob;
	//float r1,r2,r3,wx,wxx,aux,frej;
	//float peie,pese1,pese2,dcosth,h1;
	//float sinthe,costhe;
	
	ushort medium=reg_dat.med;
	particle_t p1=p;

	/* calculate some useful constants */
	float eie=p.e;
	float t=eie-rm;
	float tau=t/rm;
	//tau1=tau+1;
	//tau12=tau1*tau1;
	float tau12=(tau+1)*(tau+1);
	//tau2=tau*tau;
	//p2=tau2+2*tau;
	float p2=tau*tau+2*tau;
	float beta2=p2/tau12;
	float wmax=0.5*(t+uj);
	float xmax=uj/wmax;
	float c1=(wmax/eie)*(wmax/eie);
	float c2=(2*tau+1)/tau12;
	float fm_s=logf(rmt2*p2/uj)-beta2-0.5f;
	float prob_s=0.66666667*fm_s*(1+xmax+xmax*xmax);
	float fm_h=2+c1-c2;
	if(fm_h<1)
	{
		fm_h=1;
	}
	float prob=fm_h+prob_s;

	float r1,frej;
	float wx;
	do{
		r1=curand_uniform(randStat);
		//r2=curand_uniform(randStat);
		//r3=curand_uniform(randStat);
		if((r1*prob)<fm_h)
		{
			r1=curand_uniform(randStat);
			//wx=1/(r2*xmax+1-r2);
			wx=1/(r1*xmax+1-r1);
			float wxx=wx*xmax;
			float aux=wxx/(2-wxx);
			frej=(1+aux*(aux-c2)+c1*wxx*wxx)/fm_h;
		}
		else 
		{
			r1=curand_uniform(randStat);
			//wx=1.0f/pow((r2*xmax*xmax*xmax+1-r2),0.333333333f);
			wx=1.0f/pow((r1*xmax*xmax*xmax+1-r1),0.333333333f);
			frej=1-logf(wx)/fm_s;
		}
		r1=curand_uniform(randStat);
	//}while(r3>=frej);
	}while(r1>=frej);

	wx=wx*uj;

	/* set-up new particles */
	float h1=(eie+prm)/t;
	float ese1=eie-wx;
	p.e = ese1;
	float dcosth=h1*(ese1-prm)/(ese1+prm);
	float sinthe=sqrtf(1-dcosth);
	float costhe=sqrtf(dcosth);
	float sinphi,cosphi;
	uphi21(costhe,sinthe,cosphi,sinphi,p,randStat);
	p.status = e_electron_step;  //go on tracking this electron, we will push it to stack in moller
	p.nRepeat = 0;  //no repeat for secondary particle

	float ese2=wx-uj+prm;
	if(ese2>ae[medium])
	{
		p1.e = ese2;
		dcosth=h1*(ese2-prm)/(ese2+prm);
		sinthe=-sqrtf(1-dcosth);
		costhe=sqrtf(dcosth);
	    uphi32(costhe,sinthe,cosphi,sinphi,p1);
		p1.status = e_electron_step;
		p1.nRepeat = 0;  //no repeat for secondary particle
		push_stack(p1,e_electron_step);
        edep=0;
	}
	else
	{
		edep=wx-uj;
	}
	relax(uj,ish,iz,edep,p,reg_dat,randStat);
	return;
}