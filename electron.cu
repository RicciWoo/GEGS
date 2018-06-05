/* * * * * * * * * * * * * * * * * * * * * *
 * Simulation Step Functions for Electrons *
 * * * * * * * * * * * * * * * * * * * * * */


//Determines the rejection function due to spin effects for
//  charge       qel  (=0 for e-, =1 for e+)
//  log(energy)  elke   //  1. MS moment  q1
//  speed        beta2  //  cos(theta)    cost
__device__ float spin_rejection(particle_t p, float elke, float beta2, float q1, float cost, bool &spin_index,
	                            char &sri, char &srj, bool is_single, region_data_t &reg_dat, curandState *randStat){
//$REAL elke,beta2,q1,cost;
//$INTEGER qel,medium;
//$LOGICAL spin_index,is_single;
//$declare_max_medium;

	//uint irl = p.region;
	//region_data_t reg_dat = region_data[irl];
	ushort medium = reg_dat.med;
	//char lelec = p.charge;
	//uchar qel = (1+lelec)/2;
	uchar qel = (1+p.charge)/2;

	//$REAL rnno,ai,qq1,aj,xi,ak;
	//$INTEGER i,j,k;
	//float rnno,ai,qq1,aj,xi,ak;
	//int   srk;
	//save i,j;

	if(spin_index){  //Determine the energy and q1 index
		spin_index = false;
		float ai;
		if(beta2>=b2spin_min){
			ai = (beta2 - b2spin_min)*dbeta2i;
			sri = ai;
			ai = ai - sri;
			sri = sri + MAXE_SPIN + 1;
		}
		else if(elke>espml){
			ai = (elke - espml)*dleneri;
			sri = ai;
			ai = ai - sri;
		}
		else{
			sri = 0;
			ai = -1;
		}
		float rnno = curand_uniform(randStat);
		if( rnno < ai ) sri = sri + 1;
		if( is_single ) srj = 0;
		else{
			float qq1 = 2*q1;
			qq1 = qq1/(1 + qq1);
			float aj = qq1*dqq1i;
			srj = aj;
			if( srj >= MAXQ_SPIN ) srj = MAXQ_SPIN;
			else{
				aj = aj - srj;
				rnno = curand_uniform(randStat);
				if(rnno<aj) srj = srj + 1;
			}
		}
	}
	float xi = sqrt(0.5*(1-cost));
	float ak = xi*MAXU_SPIN;
	char srk = ak;
	ak = ak - srk;
	float spin_reject = (1-ak)*spin_rej[(((srk*MAXQ_SPIN_PLUS_1+srj)*MAXE_SPI1_PLUS_1+sri)*2+qel)*MXMED+medium]
	                       + ak*spin_rej[((((srk+1)*MAXQ_SPIN_PLUS_1+srj)*MAXE_SPI1_PLUS_1+sri)*2+qel)*MXMED+medium];
						   //corrected by Wu, 20130719
	return spin_reject;
}

//Function to sample multiple electron scattering angles from the exact distribution resulting from
//elastic scattering described by the screened Rutherford cross section (spin_effects=.false.) or
//by the screened Rutherford cross times Mott correction (spin_effects=true)
__device__ void mscat(float lambda, float chia2, float q1, float elke, float beta2, particle_t p,
	                  bool spin_effects, bool &find_index, char &msi, char &msj, float &omega2, bool &spin_index,
					  float &cost, float &sint, region_data_t &reg_dat, curandState *randStat){
//$REAL lambda, chia2,q1,elke,beta2,cost,sint;
//$INTEGER qel,medium;
//$LOGICAL spin_effects,find_index,spin_index;

	//uint irl = p.region;
	//region_data_t reg_dat = region_data[irl];
	//ushort medium = reg_dat.med;
	//char lelec = p.charge;
	//uchar qel = (1+lelec)/2;

	//$declare_write_buffer;

	//$REAL sprob,explambda,wsum,wprob,xi,rejf,spin_rejection,
	//      cosz,sinz,phi,omega2,llmbda,ai,aj,ak,a,u,du,x1,rnno;
	//$INTEGER icount,i,j,k;
	//save i,j,omega2;
	//float sprob,explambda,wsum,wprob,xi,rejf,cosz,sinz,phi,llmbda,ai,aj,ak,a,u,du,x1,rnno;
	//int   icount,msk;
	char  sri=0,srj=0;  //variables for spin_rejection

	if(lambda<=13.8){
		//Test only for lambda = 13.8 implies a 1e-6 error, ie large-lambda cases
		//that contribute to the forward no-scattering amplitude.
		//sprob = get_rand(idx);
		float sprob = curand_uniform(randStat);
		float explambda = exp(-lambda);
		if(sprob<explambda){
			//It was a no scattering event
			cost = 1;
			sint = 0;
			return;
		}
		float wsum = (1+lambda)*explambda;
		if(sprob<wsum){
//RETRY_1:  //:RETRY_1:;
			float xi,rejf,rnno;
			do{
				xi = curand_uniform(randStat);
				xi = 2*chia2*xi/(1 - xi + chia2);
				cost = 1 - xi;
				if(spin_effects){
					rejf = spin_rejection(p,elke,beta2,q1,cost,spin_index,sri,srj,false,reg_dat,randStat);
					rnno = curand_uniform(randStat);
					//if(rnno>rejf) goto RETRY_1;
				}
			}while(spin_effects&&rnno>rejf);
			sint = sqrt(xi*(2 - xi));
			return;
		}
		if(lambda<=1){
			//this portion is introduced because with alternative BCAs mscat can be called with
			//lambda < 1 where there are no pre-calculated data
			float wprob = explambda;
			wsum = explambda;
			cost = 1;
			sint = 0;
			char icount = 0;
			do{
				icount = icount + 1;
				if(icount > 20) break;  //To avoid underflow if sprob very close to 1
				wprob = wprob*lambda/icount;
				wsum = wsum + wprob;
//RETRY_2:  //:RETRY_2:;
				float xi,rejf,rnno,cosz;
				do{
					xi = curand_uniform(randStat);
					xi = 2*chia2*xi/(1 - xi + chia2);
					cosz = 1 - xi;
					if(spin_effects){
						rejf = spin_rejection(p,elke,beta2,q1,cosz,spin_index,sri,srj,false,reg_dat,randStat);
						rnno = curand_uniform(randStat);
						//if(rnno>rejf) goto RETRY_2;
					}
				}while(spin_effects&&rnno>rejf);
				float sinz = xi*(2 - xi);
				if(sinz>1.0e-20){
					sinz = sqrt(sinz);
					xi = curand_uniform(randStat);
					float phi = xi*6.2831853;
					cost = cost*cosz - sint*sinz*cos(phi);
					sint = sqrt(max(0.0,(1-cost)*(1+cost)));
				}
			}while(wsum <= sprob);
			return;
		}
	}

	//It was a multiple scattering event. Sample the angle from the q^(2+) surface
	if(lambda<=LAMBMAX_MS){
		if(find_index){
			float llmbda = logf(lambda);
			//First fix lambda bin
			float ai = llmbda*dllambi;
			msi = ai;
			ai = ai - msi;
			float xi = curand_uniform(randStat);
			if(xi < ai) msi = msi + 1;
			//fix now q1 bin
			if(q1<QMIN_MS) msj = 0;
			else if(q1<QMAX_MS){
				float aj = q1*dqmsi;
				msj = aj;
				aj = aj - msj;
				xi = curand_uniform(randStat);
				if(xi<aj) msj = msj + 1;
			}
			else msj = MAXQ_MS;

			//Calculate omega2
			if(llmbda<2.2299)
				omega2 = chia2*(lambda + 4)*(1.347006 + llmbda*(0.209364 - llmbda*(0.45525 - llmbda*(0.50142 - 0.081234*llmbda))));
			else
				omega2 = chia2*(lambda + 4)*(-2.77164 + llmbda*(2.94874 - llmbda*(0.1535754 - llmbda*0.00552888)));

			find_index = false;
		}
		//If this is a re-iteration with the same lambda, then omega2, i, and k
		//should have been defined in the previous iteration

//RETRY_3:  //:RETRY_3:;
		float xi,rejf,rnno;
		do{
			xi = curand_uniform(randStat);
			float ak = xi*MAXU_MS;
			char msk = ak;
			ak = ak - msk;
			/*
			if(ak > wms_array[(msk * (MAXU_MS+1) + msj) * (MAXQ_MS+1) + msi])  //removed -1 by Wu, 20130610
				msk = ims_array[(msk * (MAXU_MS+1) + msj) * (MAXQ_MS+1) + msi];
			float a = fms_array[(msk * (MAXU_MS+1) + msj) * (MAXQ_MS+1) + msi];
			float u = ums_array[(msk * (MAXU_MS+1) + msj) * (MAXQ_MS+1) + msi];
			float du = ums_array[((msk+1) * (MAXU_MS+1) + msj) * (MAXQ_MS+1) + msi] - u;
			*/
            if(ak > wms_array[(msk * MAXQ_MS_PLUS_1 + msj) * (MAXL_MS_PLUS_1) + msi])  //removed -1 by Wu, 20130610
                    msk = ims_array[(msk * MAXQ_MS_PLUS_1 + msj) * (MAXL_MS_PLUS_1) + msi];
            int itmp=(msk * MAXQ_MS_PLUS_1 + msj) * (MAXL_MS_PLUS_1) + msi;
            float a = fms_array[itmp];
            float u = ums_array[itmp];
            float du = ums_array[((msk+1) * MAXQ_MS_PLUS_1 + msj) * MAXL_MS_PLUS_1 + msi] - u;  //corrected by Tong Xu, 20131120

			xi = curand_uniform(randStat);
			if(abs(a)<0.2){
				float x1 = 0.5*(1-xi)*a;
				u = u + xi*du*(1+x1*(1-xi*a));
			}
			else
				u = u - du/a*(1-sqrt(1+xi*a*(2+a)));
			//Sample u from the correspoding q_{SR}^{(2+)} distribution
			
			xi = omega2*u/(1 + 0.5*omega2 - u);  //xi = 2*omega2*u/(1 + omega2 - u);  ??
			if(xi>1.99999) xi = 1.99999;
			//"some machines have trouble when xi is very close to 2 in subsequent calculations.
			cost = 1 - xi;
			if(spin_effects){
				rejf = spin_rejection(p,elke,beta2,q1,cost,spin_index,sri,srj,false,reg_dat,randStat);
				rnno = curand_uniform(randStat);
				//if(rnno>rejf) goto RETRY_3;
			}
		}while(spin_effects&&rnno>rejf);
		sint = sqrt(xi*(2-xi));
		return;
	}
}

//single elastic scattering
__device__ void sscat(float chia2, float elke, float beta2, particle_t p, bool spin_effects, float &cost, float &sint,
	                  region_data_t &reg_dat, curandState *randStat){
//$REAL chia2,elke,beta2,cost,sint;
//$INTEGER qel,medium;
//$LOGICAL spin_effects;

	//uint irl = p.region;
	//region_data_t reg_dat = region_data[irl];
	//ushort medium = reg_dat.med;
	//char lelec = p.charge;
	//uchar qel = (1+lelec)/2;

	//$REAL xi,rnno,rejf,spin_rejection,qzero;
	//$LOGICAL spin_index;
	//float xi,rnno,rejf,qzero;
	//bool  spin_index;
	char  sri=0,srj=0;  //variables for spin_rejection

	bool  spin_index = true;
//RETRY_SPIN:  //:RETRY-SPIN:;
	float xi,rejf,rnno;
	do{
		xi = curand_uniform(randStat);
		xi = 2*chia2*xi/(1 - xi + chia2);
		cost = 1 - xi;
		if(spin_effects){
			//qzero=0;
			rejf = spin_rejection(p,elke,beta2,0.0,cost,spin_index,sri,srj,true,reg_dat,randStat);
			rnno = curand_uniform(randStat);
			//if(rnno>rejf) goto RETRY_SPIN;
		}
	}while(spin_effects&&rnno>rejf);
	sint = sqrt(xi*(2 - xi));
	return;
}

//This function models multiple elastic scattering and spatial deflections for a given path-length tustep.
__device__ void msdist_pII(float eloss, float tustep, bool spin_effects, particle_t &p0, char &msi, char &msj,  
	                       float &omega2, float &ustep, region_data_t &reg_dat, curandState *randStat){
//Input variables
//===============
//$REAL
	//e0,     "electron kinetic energy at the beginning of step
	//eloss,  "energy loss for this step
	//rhof,   "density scaling template (as in EGS)
	//tustep, "total pathlength of the step,
	//u0,     "x-direction cosine before scattering
	//v0,     "y-direction cosine before scattering
	//w0,     "z-direction cosine before scattering
	//x0,     "initial x-position
	//y0,     "initial y-position
	//z0;     "initial z-position
//$INTEGER
	//medium, "medium number
	//qel;    "=0 for e-, =1 for e+, needed for spin effects
//$LOGICAL
	//spin_effects;

	//uint irl = p0.region;
	//region_data_t reg_dat = region_data[irl];
	ushort medium = reg_dat.med;
	//float rhof = reg_dat.rhof;
	//char lelec = p0.charge;
	//uchar qel = (1+lelec)/2;
	float e0 = p0.e - rm;

//Output variables
//================
//$REAL
	//us,     "x-direction cosine after scattering
	//vs,     "y-direction cosine after scattering
	//ws,     "z-direction cosine after scattering
	//xf,     "final x-position after transport
	//yf,     "final y-position after transport
	//zf,     "final z-position after transport
	//ustep;  "straight line distance between the initial and final position

	float us;     //"x-direction cosine after scattering
	float vs;     //"y-direction cosine after scattering
	float ws;     //"z-direction cosine after scattering
	//float xf;     //"final x-position after transport
	//float yf;     //"final y-position after transport
	//float zf;     //"final z-position after transport

//Local variables
//===============
//$REAL
	//float blccc    = 0.0F;     //multiple scattering parameter
	//float xcccc    = 0.0F;     //multiple scattering parameter
	//float b        = 0.0F;     //substep transport distance
	//float c        = 0.0F;     //substep transport distance
	//float eta,eta1;            //randomization of the substep transport distances
	//float chia2    = 0.0F;     //screening angle, note: our chia2 is Moliere's chia2/4
	//float chilog   = 0.0F;     //log(1+1/chia2)
	//float cphi0    = 0.0F;     //cosine of the azimuthal angle of the initial particle relative to its coordinates
	//float cphi1    = 0.0F;     //cosine of the first azimuthal angle
	//float cphi2    = 0.0F;     //cosine of the second azimuthal angle
	//float w1       = 0.0F;     //cosine of the first substep polar scattering angle
	//float w2       = 0.0F;     //cosine of the second substep polar scattering angle
	//float w1v2     = 0.0F;     //w1*v2;
	//float delta    = 0.0F;     //transport parameter (see paper)
	//float e        = 0.0F;     //average kinetic energy over the step
	//float elke     = 0.0F;     //Log(e)"
	//float beta2    = 0.0F;     //speed at e in units of c, squared
	//float etap     = 0.0F;     //correction to the screening parameter derived from PWA
	//float xi_corr  = 0.0F;     //correction to the first MS moments due to spin
	//float ms_corr  = 0.0f;
	//float tau      = 0.0F;     //average kinetic energy over the step divided by electron mass
	//float tau2     = 0.0F;     //tau squared
	//float epsilon  = 0.0F;     //fractional energy loss
	//float epsilonp = 0.0F;     //fractional energy loss
	//float temp,temp1;          //auxilarity variables for energy loss corrections
	//float temp2;               //
	//float factor   = 0.0F;     //intermediate factor employed in the energy-loss calculations
	//float gamma    = 0.0F;     //q2/q1
	//float lambda   = 0.0F;     //distance in number of elastic scattering mean free paths
                               //for each sample of the multiple scattering angle
	//float p2       = 0.0F;     //average momentum over the step
	//float p2i      = 0.0F;     //inverse of ap2
	//float q1       = 0.0F;     //first moment of the single scattering cross section
	//float sint0    = 0.0F;     //sine of the initial particle relative to its coordinates
	//float sint02   = 0.0F;     //sint0**2
	//float sint0i   = 0.0F;     //1/sint0
	//float sint1    = 0.0F;     //sine of the first substep polar scattering angle
	//float sint2    = 0.0F;     //sine of the second substep polar scattering angle
	//float sphi0    = 0.0F;     //sine of the azimuthal angle of the initial particle relative to its coordinates
	//float sphi1    = 0.0F;     //sine of the first azimuthal angle
	//float sphi2    = 0.0F;     //sine of the second azimuthal angle
	//float u2p      = 0.0F;     //intermediate scatter or transport direction cosine
	//float u2       = 0.0F;     //sint2*cphi2;
	//float v2       = 0.0F;     //sint2*sphi2;
	//float ut       = 0.0F;     //x-direction cosine for transport
	//float vt       = 0.0F;     //y-direction cosine for transport
	//float wt       = 0.0F;     //z-direction cosine for transport
	//float xi       = 0.0F;     //first GS - moment
	//float rhophi2  = 0.0F;     //xphi**2 + yphi**2 or its inverse
	//float xphi     = 0.0F;     //x - used to calculated azimuthal angles
	//float xphi2    = 0.0F;     //xphi**2
	//float yphi     = 0.0F;     //y - used to calculated azimuthal angles
	//float yphi2    = 0.0F;     //yphi**2
	//bool  find_index = false;  //needed to save locating the q2 index in the 2. call to mscat"
	//bool  spin_index = false;  //saves locating the spin rejection index in 2. call to mscat
	//uint  lelke    = 0;

	//uint  count_pII_steps = 0;

	//$declare_max_medium;

	//atomicAdd(&count_pII_steps, 1);
	//count_pII_steps = count_pII_steps + 1;
	float blccc = blcc[medium];
	float xcccc = xcc[medium];

	//Commonly used factors
	float e = e0 - 0.5*eloss;
	float tau = e/0.5110034;
	float tau2 = tau*tau;
	float epsilon = eloss/e0;
	float epsilonp = eloss/e;
	//e = e * (1 - epsilonp*epsilonp*((6+tau*(10+5*tau))/(tau+1)/(tau+2))/24);
	e = e*(1 - epsilonp*epsilonp*(6+10*tau+5*tau2)/(24*tau2+48*tau+72));    //the coefficient may be not correct
	//e = e*(1 - epsilonp*epsilonp*(6+10*tau+5*tau2)/(24*tau2+72*tau+48));  //corrected
	float p2 = e*(e + rmt2);
	//p2i = 1/p2;
	float beta2 = p2/(p2 + rmsq);
	//chia2 = xcccc*p2i/(4*blccc);
	float chia2 = xcccc/p2/(4*blccc);
	//lambda = 0.5*tustep*rhof*blccc/beta2;  //The 0.5 implies a half-step
	float lambda = 0.5*tustep*reg_dat.rhof*blccc/beta2;
	float temp = epsilonp/((tau+1)*(tau+2));
	//temp2 = 0.166666*(4+tau*(6+tau*(7+tau*(4+tau))))*temp*temp;
	float temp1 = 0.166666*(4+tau*(6+tau*(7+tau*(4+tau))))*temp*temp;
	//temp2 = 0.166666*(4+tau*(6+tau*(7+tau*(4+tau))))*(epsilonp/((tau+1)*(tau+2)))**2;
	//lambda = lambda*(1 - temp2);
	lambda = lambda*(1 - temp1);

	float etap,xi_corr,gamma,ms_corr;
	float elke;
	if(spin_effects){
		elke = logf(e);
		//$SET INTERVAL elke,eke;
		float2 eke_dat = eke01[medium];
		uint lelke = eke_dat.x + eke_dat.y * elke;
		if(lelke<1){  //This should normally not happen
			lelke = 1;
			elke = (1 - eke01[medium].x)/eke01[medium].y;
		}
		/*
		float2 eta_ms_dat,q1c_ms_dat,q2c_ms_dat;
		//if(qel==0){
		if(p0.charge<0){
			//$EVALUATE etap USING etae_ms(elke);
			//float2 etae_ms_dat = etae_ms[medium * MXEKE + lelke-1];
			eta_ms_dat = etae_ms[medium * MXEKE + lelke-1];
			//etap = etae_ms_dat.x + etae_ms_dat.y * elke;
			//$EVALUATE xi_corr USING q1ce_ms(elke);
			//float2 q1ce_ms_dat = q1ce_ms[medium * MXEKE + lelke-1];
			q1c_ms_dat = q1ce_ms[medium * MXEKE + lelke-1];
			//xi_corr = q1ce_ms_dat.x + q1ce_ms_dat.y * elke;
			//$EVALUATE gamma USING q2ce_ms(elke);
			//float2 q2ce_ms_dat = q2ce_ms[medium * MXEKE + lelke-1];
			q2c_ms_dat = q2ce_ms[medium * MXEKE + lelke-1];
			//gamma = q2ce_ms_dat.x + q2ce_ms_dat.y * elke;
		}
		else{
			//$EVALUATE etap USING etap_ms(elke);
			//float2 etap_ms_dat = etap_ms[medium * MXEKE + lelke-1];
			eta_ms_dat = etap_ms[medium * MXEKE + lelke-1];
			//etap = etap_ms_dat.x + etap_ms_dat.y * elke;
			//$EVALUATE xi_corr USING q1cp_ms(elke);
			//float2 q1cp_ms_dat = q1cp_ms[medium * MXEKE + lelke-1];
			q1c_ms_dat = q1cp_ms[medium * MXEKE + lelke-1];
			//xi_corr = q1cp_ms_dat.x + q1cp_ms_dat.y * elke;
			//$EVALUATE gamma USING q2cp_ms(elke);
			//float2 q2cp_ms_dat = q2cp_ms[medium * MXEKE + lelke-1];
			q2c_ms_dat = q2cp_ms[medium * MXEKE + lelke-1];
			//gamma = q2cp_ms_dat.x + q2cp_ms_dat.y * elke;
		}
		*/
		uchar qel = (1+p0.charge)/2;
		float2 eta_ms_dat = eta_ms[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
		float2 q1c_ms_dat = q1c_ms[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
		float2 q2c_ms_dat = q2c_ms[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
		etap = eta_ms_dat.x + eta_ms_dat.y * elke;
		xi_corr = q1c_ms_dat.x + q1c_ms_dat.y * elke;
		gamma = q2c_ms_dat.x + q2c_ms_dat.y * elke;
		//$EVALUATE ms_corr USING blcce(elke);
		float2 blcce_dat = blcce[medium * MXEKE + lelke-1];
		ms_corr = blcce_dat.x + blcce_dat.y * elke;
	}
	else{
		etap = 1;
		xi_corr = 1;
		gamma = 1;
		ms_corr = 1;
	}

	chia2 = chia2*etap;
	lambda = lambda/(etap*(1+chia2))*ms_corr;
	float chilog = logf(1 + 1/chia2);
	float q1 = 2*chia2*(chilog*(1 + chia2) - 1);
	gamma = 6*chia2*(1 + chia2)*(chilog*(1 + 2*chia2) - 2)/q1*gamma;
	float xi = q1*lambda;

	//Sample first substep scattering angle
	bool find_index = true;
	bool spin_index = true;
	float w1,sint1;
	mscat(lambda,chia2,xi,elke,beta2,p0,spin_effects,find_index,msi,msj,omega2,spin_index,w1,sint1,reg_dat,randStat);
	//Macro for azimuthal angle selection using a sampling within a box method
	//Choose a point randomly within a box such that -1 <= x <= 1 and 0 <= y <= 1
	//Reject the set if it lies without the inscribed unit semicircle centered at (x,y) = (0,0)
	//once out of the loop, use the trigonimetric relations (TeX notation)
	//\cos 2\phi = (x^2 - y^2)/(x^2 + y^2)    \sin 2\phi = 2xy/(x^2 + y^2)
	//$SELECT-AZIMUTHAL-ANGLE(cphi1,sphi1);->
	float xphi,yphi,rhophi2;
	do{
		//xphi = get_rand(idx);
		xphi = curand_uniform(randStat);
		xphi = 2*xphi - 1;
		//xphi2 = xphi*xphi;
		//yphi = get_rand(idx);
		yphi = curand_uniform(randStat);
		//yphi2 = yphi*yphi;
		//rhophi2 = xphi2 + yphi2;
		rhophi2 = xphi*xphi + yphi*yphi;
	}while(rhophi2>1);
	rhophi2 = 1/rhophi2;
	//cphi1 = (xphi2 - yphi2)*rhophi2;
	float cphi1 = (xphi*xphi - yphi*yphi)*rhophi2;
	float sphi1 = 2*xphi*yphi*rhophi2;
	//end of $SELECT-AZIMUTHAL-ANGLE(cphi1,sphi1);

	//Sample second substep scattering angle
	float w2,sint2;
	mscat(lambda,chia2,xi,elke,beta2,p0,spin_effects,find_index,msi,msj,omega2,spin_index,w2,sint2,reg_dat,randStat);
	//$SELECT-AZIMUTHAL-ANGLE(cphi2,sphi2);->
	do{
		//xphi = get_rand(idx);
		xphi = curand_uniform(randStat);
		xphi = 2*xphi - 1;
		//xphi2 = xphi*xphi;
		//yphi = get_rand(idx);
		yphi = curand_uniform(randStat);
		//yphi2 = yphi*yphi;
		//rhophi2 = xphi2 + yphi2;
		rhophi2 = xphi*xphi + yphi*yphi;
	}while(rhophi2>1);
	rhophi2 = 1/rhophi2;
	//cphi2 = (xphi2 - yphi2)*rhophi2;
	float cphi2 = (xphi*xphi - yphi*yphi)*rhophi2;
	float sphi2 = 2*xphi*yphi*rhophi2;
	//end of $SELECT-AZIMUTHAL-ANGLE(cphi2,sphi2);

	//Final direction of motion, relative to z-axis motion
	float u2 = sint2*cphi2;
	float v2 = sint2*sphi2;
	float u2p = w1*u2 + sint1*w2;
	us = u2p*cphi1 - v2*sphi1;
	vs = u2p*sphi1 + v2*cphi1;
	ws = w1*w2 - sint1*u2;

	//Calculate delta, b, c
	xi = 2*xi*xi_corr;  //xi was for half step, xi_corr corrects for spin effects

	//eta = get_rand(idx);
	float eta = curand_uniform(randStat);
	eta = sqrt(eta);  //eta is a random number sampled from 2*eta*d(eta)
	float eta1 = 0.5*(1 - eta);
	float delta = 0.9082483-(0.1020621-0.0263747*gamma)*xi;

	//Correct the coefficients for energy loss
	temp1 = 2 + tau;
	temp = (2+tau*temp1)/((tau+1)*temp1);
	//Take logarithmic dependence into account as well
	temp = temp - (tau+1)/((tau+2)*(chilog*(1+chia2)-1));
	temp = temp * epsilonp;
	temp1 = 1 - temp;
	delta = delta + 0.40824829*(epsilon*(tau+1)/((tau+2)*(chilog*(1+chia2)-1)*(chilog*(1+2*chia2)-2)) - 0.25*temp*temp);
	//delta = delta + 0.40824829*epsilonp/(tau+1)/((tau+2)*(chilog*(1+chia2)-1)*(chilog*(1+2*chia2)-2));  //why not this
	        //0.40824829 is 1/Sqrt(6)
	float b = eta*delta;
	float c = eta*(1-delta);

	//Calculate transport direction cosines
	float w1v2 = w1*v2;
	float ut = b*sint1*cphi1 + c*(cphi1*u2 - sphi1*w1v2) + eta1*us*temp1;
	float vt = b*sint1*sphi1 + c*(sphi1*u2 + cphi1*w1v2) + eta1*vs*temp1;
	float wt = eta1*(1+temp) + b*w1 + c*w2 + eta1*ws*temp1;

	//Calculate transport distance
	ustep = tustep*sqrt(ut*ut + vt*vt + wt*wt);

	//Rotate into the final direction of motion and transport relative to original direction of motion
	float sint02 = p0.u*p0.u + p0.v*p0.v;
	if(sint02>1e-20){
		float sint0 = sqrt(sint02);
		//sint0i = 1/sint0;
		//cphi0 = sint0i*p0.u;
		float cphi0 = p0.u/sint0;
		//sphi0 = sint0i*p0.v;
		float sphi0 = p0.v/sint0;

		//Scattering angles
		u2p = p0.w*us + sint0*ws;
		ws = p0.w*ws - sint0*us;
		us = u2p* cphi0 - vs*sphi0;
		vs = u2p*sphi0 + vs*cphi0;

		//Transport angles
		u2p = p0.w*ut + sint0*wt;
		wt = p0.w*wt - sint0*ut;
		ut = u2p*cphi0 - vt*sphi0;
		vt = u2p*sphi0 + vt*cphi0;
	}
	else{
		wt = p0.w*wt;
		ws = p0.w*ws;
	}

	//Transport
	p0.x = p0.x + tustep*ut;
	p0.y = p0.y + tustep*vt;
	p0.z = p0.z + tustep*wt;
	p0.u = us;
	p0.v = vs;
	p0.w = ws;

	return;
}

//This subroutine models multiple elastic scattering and spatial deflections for a given path-length tustep
//resampling PRESTA-I behaviour.
__device__ void msdist_pI(float eloss, float tustep, bool spin_effects, particle_t &p0,
	                      float &ustep, region_data_t &reg_dat, curandState *randStat){
//Input variables
//===============
//$REAL
	//e0,     "electron kinetic energy at the beginning of step
	//eloss,  "energy loss for this step
	//rhof,   "density scaling template (as in EGS)
	//tustep, "total pathlength of the step,
	//u0,     "x-direction cosine before scattering
	//v0,     "y-direction cosine before scattering
	//w0,     "z-direction cosine before scattering
	//x0,     "initial x-position
	//y0,     "initial y-position
	//z0;     "initial z-position
//$INTEGER
	//medium, "medium number
	//qel;    "=0 for e-, =1 for e+, needed for spin effects
//$LOGICAL
	//spin_effects;

	//uint irl = p0.region;
	//region_data_t reg_dat = region_data[irl];
	ushort medium = reg_dat.med;
	//float rhof = reg_dat.rhof;
	//char lelec = p0.charge;
	//uchar qel = (1+lelec)/2;
	float e0 = p0.e - rm;

//Output variables
//================
//$REAL
	//us,     "x-direction cosine after scattering
	//vs,     "y-direction cosine after scattering
	//ws,     "z-direction cosine after scattering
	//xf,     "final x-position after transport
	//yf,     "final y-position after transport
	//zf,     "final z-position after transport
	//ustep;  "straight line distance between the initial and final position

	float us;     //"x-direction cosine after scattering
	float vs;     //"y-direction cosine after scattering
	float ws;     //"z-direction cosine after scattering
	//float xf;     //"final x-position after transport
	//float yf;     //"final y-position after transport
	//float zf;     //"final z-position after transport

//Local variables
//===============
//$REAL
	//float blccc    = 0.0F;     //multiple scattering parameter
	//float xcccc    = 0.0F;     //multiple scattering parameter
	//float z,r,z2,r2;           //used to calculate PLC and lateral deflection a la PRESTA-I
	//float r2max;
	//float chia2    = 0.0F;     //screening angle, note: our chia2 is Moliere's chia2/4
	//float chilog   = 0.0F;     //log(1+1/chia2)
	//float cphi0    = 0.0F;     //cosine of the azimuthal angle of the initial particle relative to its coordinates
	//float cphi     = 0.0F;     //cosine of the azimuthal scattering angle
	//float sphi     = 0.0F;     //sine of the azimuthal scattering angle
	//float e        = 0.0F;     //average kinetic energy over the step
	//float elke     = 0.0F;     //Log(e)
	//float beta2    = 0.0F;     //speed at e in units of c, squared
	//float etap     = 0.0F;     //correction to the screening parameter derived from PWA
	//float xi_corr  = 0.0F;     //correction to the first MS moments due to spin
	//float ms_corr  = 0.0f;
	//float epsilon  = 0.0F;     //fractional energy loss
	//float temp     = 0.0F;     //auxilarity variable for energy loss corrections
	//float factor   = 0.0F;     //intermediate factor employed in the energy-loss calculations
	//float gamma    = 0.0F;     //q2/q1
	//float lambda   = 0.0F;     //distance in number of elastic scattering mean free paths
	//float p2       = 0.0F;     //average momentum over the step
	//float p2i      = 0.0F;     //inverse of p2
	//float q1       = 0.0F;     //first moment of the single scattering cross section
	//float sint     = 0.0F;     //sine of the MS angle
	//float sint0    = 0.0F;     //sine of the initial particle relative to its coordinates
	//float sint02   = 0.0F;     //sint0**2
	//float sint0i   = 0.0F;     //1/sint0
	//float sphi0    = 0.0F;     //sine of the azimuthal angle of the initial particle relative to its coordinates
	//float u2p      = 0.0F;     //intermediate scatter or transport direction cosine
	//float ut       = 0.0F;     //x-direction cosine for transport
	//float vt       = 0.0F;     //y-direction cosine for transport
	//float wt       = 0.0F;     //z-direction cosine for transport
	//float xi       = 0.0F;     //first GS - moment
	//float rhophi2  = 0.0F;     //xphi**2 + yphi**2 or its inverse
	//float xphi     = 0.0F;     //x - used to calculated azimuthal angles
	//float xphi2    = 0.0F;     //xphi**2
	//float yphi     = 0.0F;     //y - used to calculated azimuthal angles
	//float yphi2    = 0.0F;     //yphi**2
	//bool  find_index = false;  //needed to save locating the q2 index in the 2. call to mscat"
	//bool  spin_index = false;  //saves locating the spin rejection index in 2. call to mscat
	//uint  lelke    = 0;
	
	float omega2;  //variables for mscat
	char  msi,msj;

	//$declare_max_medium;

	float blccc = blcc[medium];
	float xcccc = xcc[medium];

	//Commonly used factors
	float e = e0 - 0.5*eloss;
	float p2 = e*(e + rmt2);
	//p2i = 1/p2;
	//chia2 = xcccc*p2i/(4*blccc);
	float chia2 = xcccc/p2/(4*blccc);
	float beta2 = p2/(p2 + rmsq);
	//lambda = tustep*rhof*blccc/beta2;
	float lambda = tustep*reg_dat.rhof*blccc/beta2;

	//Account for energy loss in the MS distribution
	float factor = 1/(1 + 0.9784671*e);  //0.9784671 = 1/(2*rm)
	float epsilon = eloss/e0;
	epsilon= epsilon/(1-0.5*epsilon);
	float temp = 0.25*(1 - factor*(1 - 0.333333*factor))*epsilon*epsilon;
	lambda = lambda*(1 + temp);

	float etap,xi_corr,ms_corr;
	float elke;
	if(spin_effects){
		elke = logf(e);
		//$SET INTERVAL elke,eke;
		float2 eke_dat = eke01[medium];
		uint lelke = eke_dat.x + eke_dat.y * elke;
		if(lelke<1){  //This should normally not happen
			lelke = 1;
			elke = (1 - eke01[medium].x)/eke01[medium].y;
		}
		/*
		float2 eta_ms_dat,q1c_ms_dat;
		//if(qel==0){
		if(p0.charge<0){
			//$EVALUATE etap USING etae_ms(elke);
			//float2 etae_ms_dat = etae_ms[medium * MXEKE + lelke-1];
			eta_ms_dat = etae_ms[medium * MXEKE + lelke-1];
			//etap = etae_ms_dat.x + etae_ms_dat.y * elke;
			//$EVALUATE xi_corr USING q1ce_ms(elke);
			//float2 q1ce_ms_dat = q1ce_ms[medium * MXEKE + lelke-1];
			q1c_ms_dat = q1ce_ms[medium * MXEKE + lelke-1];
			//xi_corr = q1ce_ms_dat.x + q1ce_ms_dat.y * elke;
		}
		else{
			//$EVALUATE etap USING etap_ms(elke);
			//float2 etap_ms_dat = etap_ms[medium * MXEKE + lelke-1];
			eta_ms_dat = etap_ms[medium * MXEKE + lelke-1];
			//etap = etap_ms_dat.x + etap_ms_dat.y * elke;
			//$EVALUATE xi_corr USING q1cp_ms(elke);
			//float2 q1cp_ms_dat = q1cp_ms[medium * MXEKE + lelke-1];
			q1c_ms_dat = q1cp_ms[medium * MXEKE + lelke-1];
			//xi_corr = q1cp_ms_dat.x + q1cp_ms_dat.y * elke;
		}
		*/
		uchar qel = (1+p0.charge)/2;
		float2 eta_ms_dat = eta_ms[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
		float2 q1c_ms_dat = q1c_ms[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
		etap = eta_ms_dat.x + eta_ms_dat.y * elke;
		xi_corr = q1c_ms_dat.x + q1c_ms_dat.y * elke;
		//$EVALUATE ms_corr USING blcce(elke);
		float2 blcce_dat = blcce[medium * MXEKE + lelke-1];
		ms_corr = blcce_dat.x + blcce_dat.y * elke;
	}
	else{
		etap = 1;
		xi_corr = 1;
		ms_corr = 1;
	}

	//chia2 = xcccc*p2i/(4*blccc)*etap;
	chia2 = xcccc/p2/(4*blccc)*etap;
	lambda = lambda/etap/(1+chia2)*ms_corr;
	float chilog = logf(1 + 1/chia2);
	float q1 = 2*chia2*(chilog*(1 + chia2) - 1);
	float xi = q1*lambda;

	//Sample first substep scattering angle
	bool find_index = true;
	bool spin_index = true;
	float sint;
	mscat(lambda,chia2,xi,elke,beta2,p0,spin_effects,find_index,msi,msj,omega2,spin_index,ws,sint,reg_dat,randStat);
	//$SELECT-AZIMUTHAL-ANGLE(cphi,sphi);->
	float xphi,yphi,rhophi2;
	do{
		//xphi = get_rand(idx);
		xphi = curand_uniform(randStat);
		xphi = 2*xphi - 1;
		//xphi2 = xphi*xphi;
		//yphi = get_rand(idx);
		yphi = curand_uniform(randStat);
		//yphi2 = yphi*yphi;
		//rhophi2 = xphi2 + yphi2;
		rhophi2 = xphi*xphi + yphi*yphi;
	}while(rhophi2>1);
	rhophi2 = 1/rhophi2;
	//cphi = (xphi2 - yphi2)*rhophi2;
	float cphi = (xphi*xphi - yphi*yphi)*rhophi2;
	float sphi = 2*xphi*yphi*rhophi2;
	//end of $SELECT-AZIMUTHAL-ANGLE(cphi,sphi);
	us = sint*cphi;
	vs = sint*sphi;

	//Correct xi used for the PLC calc. for spin effects
	xi = xi*xi_corr;

	//Calculate PLC and lateral transport a la PRESTA-I
	//Note that we use here the exact expression for <z> because it is much simpler and
	//faster than the original PRESTA-I formulas (which are also second order approximations)
	float z,r,r2,z2,r2max;
	if(xi<0.1)
		z = 1 - xi*(0.5 - xi*(0.166666667 - 0.041666667*xi));
	else
		z = (1 - exp(-xi))/xi;
	r = 0.5*sint;
	r2 = r*r;
	z2 = z*z;
	r2max = 1 - z2;
	if(r2max<r2){
		r2 = r2max;
		r = sqrt(r2);
	}

	//Calculate final position vector
	float ut = r*cphi;
	float vt = r*sphi;
	float wt = z;

	//Calculate transport distance
	ustep = sqrt(z2 + r2)*tustep;

	//Rotate into the final direction of motion and transport relative to original direction of motion
	float sint02 = p0.u*p0.u + p0.v*p0.v;
	if(sint02>1e-20){
		float sint0 = sqrt(sint02);
		//sint0i = 1/sint0;
		//cphi0 = sint0i*p0.u;
		float cphi0 = p0.u/sint0;
		//sphi0 = sint0i*p0.v;
		float sphi0 = p0.v/sint0;

		//Scattering angles
		float u2p = p0.w*us + sint0*ws;
		ws = p0.w*ws - sint0*us;
		us = u2p*cphi0 - vs*sphi0;
		vs = u2p*sphi0 + vs*cphi0;

		//Transport angles
		u2p = p0.w*ut + sint0*wt;
		wt = p0.w*wt - sint0*ut;
		ut = u2p*cphi0 - vt*sphi0;
		vt = u2p*sphi0 + vt*cphi0;
	}
	else{
		wt = p0.w*wt;
		ws = p0.w*ws;
	}

	//Transport
	p0.x = p0.x + tustep*ut;
	p0.y = p0.y + tustep*vt;
	p0.z = p0.z + tustep*wt;
	p0.u = us;
	p0.v = vs;
	p0.w = ws;

	return;
}

// Transport the electron one step through the MLC or phantom and determine which (if any) interaction
// takes place next.
// This is the subroutine ELECTR in the file egsnrc.mortran (v 1.72 2011/05/05) of the EGSnrc 
// code system.
__global__ void electr_step()
{
	indices idx = get_indices();
	curandState randStat=devStates[idx.p];  //copy the random number generator state to local.
	
	__syncthreads();

#ifdef DO_STEP_COUNT
    volatile uint *step_counters = step_counters_shared[idx.w];

    // reset step counts
    if (idx.t < NUM_CAT) step_counters[idx.t] = 0;
#endif

	region_data_t reg_dat;
	particle_t p_original, p;
	p_original.nRepeat = 0;
	p.status = e_empty;
	bool process;
	
	if(idx.p==0) simu_stack_is_empty = false;
	
	float  elke     = 0.0F;       //(eke_val = kinetic enregy, rm = rest mass, all in units of MeV)
	float  demfp    = 0.0F;       //(dempf = differential electron mean free path)
	uint   lelke    = 0;          //index into the energy grid of tabulated functions
	float  sig0     = 0.0F;       //cross section before density scaling but before a step
	bool   compute_tstep = true;  //MFP resampled => calculate distance to the interaction in the USTEP loop
	float  total_tstep = 0.0F;    //total path-length to next discrete interaction
	bool   is_tstep = true;

	float omega2;  //variables for mscat
	char  msi,msj;

	for(;;){

		if(p.status==e_empty && (!simu_stack_is_empty||p.nRepeat>0))
			p.status = e_new_particle;

		process = (p.status == e_new_particle);

		if(process){
			is_tstep = true;
			if(p_original.nRepeat>0){
				p_original.nRepeat --;
				p = p_original;
			}
			else{
				if(pop_stack(p_original,e_electron_step))
					p = p_original;
				else{
					p.status = e_empty;
					//process = false;
				}
			}
		}

		bool newparticle = (p.status == e_new_particle);
		
#ifdef DO_STEP_COUNT
		uint count_mask = __ballot(newparticle);
		if (idx.t == 0) step_counters[e_new_particle] += __popc(count_mask);
#endif

		if(newparticle) p.status = e_electron_step;

        if(simu_stack_is_empty){
            process = (p.status == e_empty);
            if(__all(process))
                break;
        }

		process = (p.status == e_electron_step);

#ifdef DO_STEP_COUNT
        count_mask = __ballot(process);
		if (idx.t == 0) step_counters[e_electron_step] += __popc(count_mask);
#endif

		if(process){

			
			uint irold = p.region;  //irold=ir(np);//Initialize previous region, no STACK, no NP
			                        //irl is ir(np) that contains the region number that the current particle is in.
			uint irl = irold;       //region number in local variable
			reg_dat = region_data[irl];  //$start_new_particle; REPLACE WITH {medium=med(irl);}
			ushort medium = reg_dat.med;  //medium index of current region
			float rhof = reg_dat.rhof;    //mass density ratio
			float  edep1      = 0.0F;

			//if (!process) return;  //changed by Wu, 20130610
			//we've already use if(process) outside to call this electr_step

			//uchar  transport_algorithm = PRESTA_II;  //= PRESTA_II or PRESTA__I
			//bool   exact_bca  = true;

			//float  elke     = 0.0F;       //(eke_val = kinetic enregy, rm = rest mass, all in units of MeV)
			//double demfp    = 0.0F;       //(dempf = differential electron mean free path)
			//float  demfp    = 0.0F;       //(dempf = differential electron mean free path)
			//uint   lelke    = 0;          //index into the energy grid of tabulated functions
			//float  sigf     = 0.0F;       //cross section before density scaling but after a step
			//float  sig0     = 0.0F;       //cross section before density scaling but before a step
			//float  sig      = 0.0F;       //cross section after density scaling but before a step
			//float  dedx0    = 0.0F;       //stopping power before density scaling
			//float  dedx     = 0.0F;       //stopping power after density scaling
			//float  tstep    = 0.0F;       //total pathlength to the next discrete interation
			//float  ustep    = 0.0F;       //total pathlength of the electron step
			//float  tustep   = 0.0F;       //projected transport distance in the direction of motion at the start of the step
			//float  tvstep   = 0.0F;       //curved path-length calculated from TVSTEP (? VSTEP?)
			//double total_de = 0.0F;       //tatal energy loss to next discrete interaction
			//float  fedep    = 0.0F;       //fractional energy loss used in stopping power calculation
			//float  ekei     = 0.0F;       //used in $CALCULATE-TSTEP-FROM-DEMFP; and $COMPUTE-RANGE;
			//float  elkei    = 0.0F;       //Log(ekei), used in $CALCULATE-TSTEP-FROM-DEMFP;
			//float  ekef     = 0.0F;       //kinetic energy after a step
			//float  elkef    = 0.0F;       //Log(ekef)
			//uint   lelkef   = 0;          //index into the energy grid of tabulated functions
			//float  eketmp   = 0.0F;       //used to evaluate average kinetic energy of a step
			//float  elktmp   = 0.0F;       //log(eketmp)
			//uint   lelktmp  = 0;          //index into the energy grid of tabulated functions
			//float  dedxmid  = 0.0F;       //stopping power at mid-step before density scaling
			//float  aux      = 0.0F;       //aux. variable
			//float  tuss     = 0.0F;       //sampled path-length to a single scattering event
			//double total_tstep = 0.0F;    //total path-length to next discrete interaction
			//float  total_tstep = 0.0F;    //total path-length to next discrete interaction
			//float  tmxs_val = 0.0F;       //electron step-size restriction
			//float  range    = 0.0F;       //electron range
			//bool   random_tustep = false; //radomize tustep option
			//float  tperp    = 0.0F;       //perpendicular distance to the closest boundary
			//char   idisc    = 0;          //flag indicating user discard
			//float  blccl    = 0.0F;       //blcc(medium)*rhof
			//float  xccl     = 0.0F;       //xcc(medium)*rhof
			//float  rmt2     = 0.0F;       //2*electron mass in MeV
			//float  beta2    = 0.0F;       //incident positron velocity in units of c
			//float  rmsq     = 0.0F;       //electron mass squared in MeV**2
			//bool   spin_effects = false;  //if .true. electron/positron spin effects are taken into account
	                                        //in the sigle and multiple elasting scattering routines
			//float  etap     = 0.0F;       //correction to Moliere screening angle from PWA cross sections
			//float  ms_corr  = 0.0F;       //xi_corr, "correction to xi due to spin effects"
			//float  ssmfp    = 0.0F;       //distance of one single elastic scattering mean free path
			//float  skindepth = 0.0F;      //skin depth employed for PRESTA-II boundary crossing
			//float  skindepth_for_bca = 0.0F;  //distance from a boundary (in elastic MFP) to switch to one of the BCAs
			//uint   count_all_steps = 0;   //set default to 0;
			//bool   is_ch_step = false;
			//bool   callhowfar = false;  //true => BCA requires a call to howfar; false => BCA does not require a call to howfar
			//bool   domultiple = false;  //ture => inexact BCA requires multiple scattering;
			//bool   dosingle   = false;  //true => exact BCA requires single scattering;
			                              //false => exact BCA requires no single scattering
			//bool   callmsdist = false;  //true => normal condensed-history transport;
			                              //false => one of the BCA's will be invoked
			//float  de         = 0.0F;   //energy loss to dedx
			//float  lambda     = 0.0F;   //distance in number of elastic scattering mean free paths
			                              //for each sample of the multiple scattering angle
			//float  lambda_max = 0.0F;
			//float  ekems      = 0.0F;   //kinetic energy used to sample MS angle (normally midpoint)
			//float  elkems     = 0.0F;   //Log(ekems)
			//uint   lelkems    = 0;      //index into the energy grid of tabulated functions
			//float  chia2      = 0.0F;   //Multiple scattering screening angle
			//float  xi         = 0.0F;   //used for PLC calculations (first GS moment times path-length)
			//float  xi_corr    = 0.0F;   //correction to xi due to spin effects

			//uint   irnew      = 0;      //region after transport
			//float  ustep0     = 0.0F;   //temporary storage for ustep
			//uint   ierust     = 0;      //error counter for negative ustep errors
			//float  vstep      = 0.0F;   //transport distance after truncation by HOWFAR
			//float  edep       = 0.0F;   //"energy deposition in MeV
			//float  e_range    = 0.0F;   //range of electron before an iarg=0 ausgab call
			//float  p2         = 0.0F;   //electron momentum times c, squared
			//float  save_de    = 0.0F;   //de saved before $DE-FLUCTUATION
			//float  eold       = 0.0F;   //energy before deduction of energy loss
			//float  enew       = 0.0F;   //energy after deduction of energy loss
			//bool   findindex  = false;  //used for mscat
			//bool   spin_index = false;  //used for mscat with spin effects
			//float  theta      = 0.0F;   //polar scattering angle
			//float  sinthe     = 0.0F;   //sin(THETA)
			//float  costhe     = 0.0F;   //cos(THETA)
			//float  uscat      = 0.0F;   //x-axis direction cosine for scattering
			//float  vscat      = 0.0F;   //y-axis direction cosine for scattering
			//float  wscat      = 0.0F;   //z-axis direction cosine for scattering
			//float  xtrans     = 0.0F;   //final x-axis position after transport
			//float  ytrans     = 0.0F;   //final y-axis position after transport
			//float  ztrans     = 0.0F;   //final z-axis position after transport
			//float  x_final,y_final,z_final;  //position at end of step
			//float  u_final,v_final,w_final;  //direction at end of step  //only set (and relevant) for electrons
			//ushort medold     = 0;      //medium index of previous region
			//float  ekeold     = 0.0F;   //kinetic energy before a step
			//float  ebr1_val   = 0.0F;   //e- branching ratio into brem
			//float  pbr1_val   = 0.0F;   //e+ branching ratio into brem
			//float  pbr2_val   = 0.0F;   //e+ branching ratio into brem or Bhabha
			//char   eii_flag   = 0.0F;   //EII flag

            //float2 sig_dat, dedx_dat, eta_ms_dat, q1c_ms_dat;

            //float omega2;  //variables for mscat
	        //int   msi,msj;

	    //:NEWELECTRON:LOOP[
	    //Go once through this loop for each 'new' electron whose charge and energy has not been checked

            char   lelec = p.charge;   //lelec=iq(np);//save charge in local variable
			                           //-1 for electrons, 0 for photon and 1 for positron
			uchar  qel = (1+lelec)/2;  //qel=(1+lelec)/2;//0 for electrons, 1 for positions
			//double peie = p.e;         //precise energy of incident electron (double precision), no need
			//float  eie = peie;         //energy incident electron (conversion to single)
			float  eie = p.e;

		//if(process){
			if (eie <= reg_dat.ecut){         //IF(eie<=ecut(irl))[go to :ECUT-DISCARD:;]
				//p.status = e_cutoff_discard;  //Ecut is the lower transport threshold
				edep1 = p.e - prm;
				escore(p.region,edep1*p.wt);
				if(p.charge>0)
					annih_at_rest(p,reg_dat,&randStat);
				else
					p.status = e_empty;
				//return;
				continue;
			}
			if (p.wt <= 0.0F){                //IF(WT(NP)=0.0)[go to :USER-ELECTRON-DISCARD:;]
				//p.status = e_user_discard;
				//edep1 += p.e + p.charge*prm;
				//escore(p.region,edep1*p.wt);    //when p.wt<=0, it is no need to score energy deposition
				p.status = e_empty;
				//return;
				continue;
			}
		//}

		//:TSTEP:LOOP[
		//Go through this loop each time we recompute distance to an interaction.

			float  eke_val = eie - rm;    //kinetic energy will be known to user even for a vacuum step

		//if (process && (p.status == e_electron_step)){
			if(is_tstep){

				compute_tstep = true;  //MFP resampled => calculate distance to the interaction in the USTEP loop

				if (medium != VACUUM){
					//Not vacuum. Must sample to see how far to next interaction.

					//Macro for selection of the electron mean-free-path
					//$SELECT-ELECTRON-MFP;->
					//float rnne1 = get_rand(idx);         //$RANDOMSET RNNE1;
					float rnne1 = curand_uniform(&randStat);
					if(rnne1==0) rnne1 = 1.0E-30F;       //IF(RNNE1.EQ.0.0)[RNNE1=1.E-30;]
					demfp = max(-logf(rnne1),EPSEMFP);   //DEMFP=MAX(-LOG(RNNE1),$EPSEMFP);  //$EPSEMFP=1.E-5
					                                     //(demfp = differential electron mean free path)
					elke = logf(eke_val);                //(eke_val = kinetic energy, rm = rest mass, all in units of MeV)
					//$SET INTERVAL elke,eke; ->         //Prepare to approximate cross section
					//float2 eke_dat = eke01[medium];        //Lelke=eke1(MEDIUM)*elke+eke0(MEDIUM);
					//lelke = (int)(eke01[medium].x + eke01[medium].y * elke);
					lelke = (int)(eke01[medium].x + eke01[medium].y * elke);

					//Macros for the fictitious method
					//The following version uses sub-threshold energy loss as a measure of path-length
					//=>cross section is actual cross section divided by restricted stopping power
					//The global maximum of this quantity called esig_e (electrons) or psig_e (positron)
					//and is determined in HATCH
					//$EVALUATIE-SIG0;->
					if(sig_ismonotone[medium * 2 + qel]){
						/*
						float2 sig_dat,dedx_dat;
						if(lelec<0){
							//$EVALUATE-SIGF;->
							//$EVALUATE sigf USING esig(elke);->
							//float2 esig_dat = esig[medium * MXEKE + lelke-1];
							sig_dat = esig[medium * MXEKE + lelke-1];
							//sigf = esig_dat.x + esig_dat.y * elke;  //esig, "used for electron cross section interpolation"
							//$EVALUATE dedx0 USING ededx(elke);
							//float2 ededx_dat = ededx[medium * MXEKE + lelke-1];
							dedx_dat = ededx[medium * MXEKE + lelke-1];
							//dedx0 = ededx_dat.x + ededx_dat.y * elke;  //ededx, "used for electron dE/dx interpolation"
							//sigf /= dedx0;
						}
						else{
							//$EVALUATE sigf USING psig(elke);->
							//float2 psig_dat = psig[medium * MXEKE + lelke-1];
							sig_dat = psig[medium * MXEKE + lelke-1];
							//sigf = psig_dat.x + psig_dat.y * elke;  //psig, "used for positron cross section interpolation"
							//$EVALUATE dedx0 USING pdedx(elke);->
							//float2 pdedx_dat = pdedx[medium * MXEKE + lelke-1];
							dedx_dat = pdedx[medium * MXEKE + lelke-1];
							//dedx0 = pdedx_dat.x + pdedx_dat.y * elke;  //pdedx, "used for positron dE/dx interpolation"
							//sigf /= dedx0;
						}//end of $EVALUATE-SIGF;
						*/
						float2 sig_dat = sig[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
						float2 dedx_dat = dedx[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
						float sigf = sig_dat.x + sig_dat.y * elke;
						float dedx0 = dedx_dat.x + dedx_dat.y * elke;
						sigf /= dedx0;
						sig0 = sigf;        //sig0, "cross section before density scaling but before a step"
					}
					else{
						/*
						if(lelec<0)
							sig0 = esig_e[medium];  //esig_e, "maximum electron cross section per energy loss for each medium"
						else
							sig0 = psig_e[medium];  //psig_e, "maximum positron cross section per energy loss for each medium"
						*/
						sig0 = sig_e[qel*MXMED+medium];
					}//end of $EVALUATIE-SIG0;
					// The fix up of the fictitious method uses cross section per energy loss.
					// Therefor, demfp/sig is sub-thrshold energy loss until the next discrete interaction occures (see below).
					// As this quantity is a single constant for a material, $SET INTERVAL is not necessary at this point.
					// However, to not completely alter the logic of the TSTEP and USTEP loops, this is left for now.
				}//end non-vacuum test

			}

		//do{  //:USTEP:LOOP[
			//Here for each check with user geometry.
			//Compute size of maxium acceptable step, which is limited by multiple scattering or other approximations.

			float  sig_val      = 0.0F;       //cross section after density scaling but before a step
			float  dedx_val     = 0.0F;       //stopping power after density scaling
			float  tstep    = 0.0F;       //total pathlength to the next discrete interation
			float  ustep    = 0.0F;       //total pathlength of the electron step
			float  tustep   = 0.0F;       //projected transport distance in the direction of motion at the start of the step
			float  tvstep   = 0.0F;       //curved path-length calculated from TVSTEP (? VSTEP?)
			float  vstep    = 0.0F;       //transport distance after truncation by HOWFAR
			float  range    = 0.0F;       //electron range
			float  de       = 0.0F;       //energy loss to dedx

			float  blccl    = 0.0F;       //blcc(medium)*rhof
			float  xccl     = 0.0F;       //xcc(medium)*rhof

			float  beta2    = 0.0F;       //incident positron velocity in units of c
			float  etap     = 0.0F;       //correction to Moliere screening angle from PWA cross sections
			float  ekems    = 0.0F;       //kinetic energy used to sample MS angle (normally midpoint)
			float  elkems   = 0.0F;       //Log(ekems)
			uint   lelkems  = 0;          //index into the energy grid of tabulated functions
			float  chia2    = 0.0F;       //Multiple scattering screening angle
			float  xi       = 0.0F;       //used for PLC calculations (first GS moment times path-length)
			float  xi_corr  = 0.0F;       //correction to xi due to spin effects

			bool   callhowfar = false;    //true => BCA requires a call to howfar; false => BCA does not require a call to howfar
			bool   domultiple = false;    //ture => inexact BCA requires multiple scattering;
			bool   dosingle   = false;    //true => exact BCA requires single scattering; false => exact BCA requires no single scattering
			bool   callmsdist = false;    //true => normal condensed-history transport; false => one of the BCA's will be invoked

			if(medium == VACUUM){  //vacuum
				tstep  = VACUUM_STEP;  //tstep = vacdst;  //vacdst = infinity (actually 10^8)
				//tstep = total pathlength to the next discrete interaction
				ustep  = tstep;  //ustep = projected transport distance in the direction of motion at the start of the step 
				tustep = ustep;  //tustep = total pathlength of the electron step
				callhowfar = true;  //Always call HOWFAR for vacuum steps!
				//Note that tustep and ustep are modified below. these provide defaults
			}
			else{  //non-vacuum
				//Density ratio scaling macro (to over-ride density in a particular region)
				//$SET-RHOF;->RHOF=RHOR(IRL)/RHO(MEDIUM);  //This is already calculated in rhof
				//density ratio scaling template, EGS allows the density to vary continuously (user option)

				//Because the cross section is interactions per energy loss, no rhof-scalling is required
				//$SCALE-SIG0;->
				sig_val = sig0;  //sig, "cross section after density scaling but before a step"
				if(sig_val<=0){
					//This can happen if the threshold for brems, (ap + rm), is greater than ae.
					//Moller threshold is 2*ae - rm. If sig is zero, we are below the thresholds for
					//both bremsstrahlung and Moller. In this case we will just lose energy by
					//ionization loss untill we go below cut-off. Do not assume range is available, so
					//just ask for step same as vacuum. Electron transport will reduce into little steps.
					//(ae is the lower threshold for creation of a secondary Moller electron,
					// ap is the lower threshold for creation of a brem.)
					tstep = VACUUM_STEP;
					sig0  = 1.0E-15F;
				}
				else{
					//Once the sub-threshold processes energy loss to the next discrete interaction is determined,
					//the corresponding path-length has to be calculated. This is done by the macro below.
					//This macro assumes the energy at the begining to be eke_val, the logarithm of it elke,
					//lelke - the corresponding interpolation index and makes use of $COMPUTE-DRANGE(#,#,#,#,#,#)
					//$CALCULATE-TSTEP-FROM-DEMFP;->
					if(compute_tstep){
						float total_de = demfp/sig_val;
						//fedep = total_de;
						//ekef  = eke_val - fedep;
						float ekef  = eke_val - total_de;
						if(ekef <= e_array[medium * MXEKE + 0])  //corrected 2013 June 5
							tstep = VACUUM_STEP;
						else{
							float elkef = logf(ekef);
							//SET INTERVAL elkef,eke;->
							//float2 eke_dat = eke01[medium];  //Lelkef=eke1(MEDIUM)*elkef+eke0(MEDIUM)
							//lelkef = (int)(eke_dat.x + eke_dat.y * elkef);
							uint lelkef = (int)(eke01[medium].x + eke01[medium].y * elkef);
							if(lelkef == lelke){
								//initial and final energy are in the same interpolation bin

								//The following macro computes the path-length traveled while going from energy {P1}
								//to energy {P2}, both energies being in the same interpolation bin, given by {P3}.
								//{P4} and {P5} are the logarithms of {P1} and {P2}. The expression is based on
								//logarithmic interpolation as used in EGSnrc (i.e. dedx = a + b*Log(E)) and a power
								//series expansion of the ExpIntegralEi function that is the result of the integration.
								//The result is returned in {P6}.
								//$COMPUTE-DRANGE(eke_val,ekef,lelke,elke,elkef,tstep);->
								float fedep = 1 - ekef/eke_val;
								float elktmp = 0.5*(elke+elkef+0.25*fedep*fedep*(1+fedep*(1+0.875*fedep)));
								//the above evaluates the logarithm of the midpoint energy
								uint lelktmp = lelke;
								/*
								float2 dedx_dat;
								if(lelec<0){
									//$EVALUATE dedxmid USING ededx(elktmp);->
									//float2 ededx_dat = ededx[medium * MXEKE + lelktmp-1];
									dedx_dat = ededx[medium * MXEKE + lelktmp-1];
									//dedxmid = ededx_dat.x + ededx_dat.y * elktmp;
									//dedxmid = 1/dedxmid;
									//aux = ededx[medium * MXEKE + lelktmp-1].y * dedxmid;
								}
								else{
									//$EVALUATE dedxmid USING pdedx(elktmp);->
									//float2 pdedx_dat = pdedx[medium * MXEKE + lelktmp-1];
									dedx_dat = pdedx[medium * MXEKE + lelktmp-1];
									//dedxmid = pdedx_dat.x + pdedx_dat.y * elktmp;
									//dedxmid = 1/dedxmid;
									//aux = pdedx[medium * MXEKE + lelktmp-1].y * dedxmid;
								}
								*/
								float2 dedx_dat = dedx[qel*MXEKE*MXMED + medium*MXEKE + lelktmp-1];
								float dedxmid = dedx_dat.x + dedx_dat.y * elktmp;
								dedxmid = 1/dedxmid;
								float aux = dedx_dat.y * dedxmid;
								float temp = (fedep/(2-fedep));
								aux = aux*(1+2*aux)*temp*temp/6;  //aux = aux*(1+2*aux)*(fedep/(2-fedep))**2/6;
								tstep = fedep*eke_val*dedxmid*(1+aux);
								//end of $COMPUTE-DRANGE(eke_val,ekef,lelke,elke,elkef,tstep);
							}
							else{
								//initial and final energy are in different interpolation bins,
								//calc range from ekef to E(lelkef+1) and from E(lelke) to eke_val
								//and add the pre-calculated range from E(lelke+1) to E(lelke)
								float ekei = e_array[medium * MXEKE + lelke-1];
								float elkei = (lelke - eke01[medium].x)/eke01[medium].y;
								//$COMPUTE-DRANGE(eke_val,ekei,lelke,elke,elkei,tuss);->
								float fedep = 1 - ekei/eke_val;
								float elktmp = 0.5*(elke+elkei+0.25*fedep*fedep*(1+fedep*(1+0.875*fedep)));
								//the above evaluates the logarithm of the midpoint energy
								uint lelktmp = lelke;
								/*
								float2 dedx_dat;
								if(lelec<0){
									//$EVALUATE dedxmid USING ededx(elktmp);->
									//float2 ededx_dat = ededx[medium * MXEKE + lelktmp-1];
									dedx_dat = ededx[medium * MXEKE + lelktmp-1];
									//dedxmid = ededx_dat.x + ededx_dat.y * elktmp;
									//dedxmid = 1/dedxmid;
									//aux = ededx[medium * MXEKE + lelktmp-1].y * dedxmid;
								}
								else{
									//$EVALUATE dedxmid USING pdedx(elktmp);->
									//float2 pdedx_dat = pdedx[medium * MXEKE + lelktmp-1];
									dedx_dat = pdedx[medium * MXEKE + lelktmp-1];
									//dedxmid = pdedx_dat.x + pdedx_dat.y * elktmp;
									//dedxmid = 1/dedxmid;
									//aux = pdedx[medium * MXEKE + lelktmp-1].y * dedxmid;
								}
								*/
								float2 dedx_dat = dedx[qel*MXEKE*MXMED + medium*MXEKE + lelktmp-1];
								float dedxmid = dedx_dat.x + dedx_dat.y * elktmp;
								dedxmid = 1/dedxmid;
								float aux = dedx_dat.y * dedxmid;
								float temp = (fedep/(2-fedep));
								aux = aux*(1+2*aux)*temp*temp/6;
								float tuss = fedep*eke_val*dedxmid*(1+aux);
								//end of $COMPUTE-DRANGE(eke_val,ekei,lelke,elke,elkei,tuss);

								ekei = e_array[medium * MXEKE + lelkef + 1-1];
								elkei = (lelkef + 1 - eke01[medium].x)/eke01[medium].y;
								//$COMPUTE-DRANGE(ekei,ekef,lelkef,elkei,elkef,tstep);->
								fedep = 1 - ekef/ekei;
								elktmp = 0.5*(elkei+elkef+0.25*fedep*fedep*(1+fedep*(1+0.875*fedep)));
								lelktmp = lelkef;
								/*
								if(lelec<0){
									//$EVALUATE dedxmid USING ededx(elktmp);->
									//float2 ededx_dat = ededx[medium * MXEKE + lelktmp-1];
									dedx_dat = ededx[medium * MXEKE + lelktmp-1];
									//dedxmid = ededx_dat.x + ededx_dat.y * elktmp;
									//dedxmid = 1/dedxmid;
									//aux = ededx[medium * MXEKE + lelktmp-1].y * dedxmid;
								}
								else{
									//$EVALUATE dedxmid USING pdedx(elktmp);->
									//float2 pdedx_dat = pdedx[medium * MXEKE + lelktmp-1];
									dedx_dat = pdedx[medium * MXEKE + lelktmp-1];
									//dedxmid = pdedx_dat.x + pdedx_dat.y * elktmp;
									//dedxmid = 1/dedxmid;
									//aux = pdedx[medium * MXEKE + lelktmp-1].y * dedxmid;
								}
								*/
								dedx_dat = dedx[qel*MXEKE*MXMED + medium*MXEKE + lelktmp-1];
								dedxmid = dedx_dat.x + dedx_dat.y * elktmp;
								dedxmid = 1/dedxmid;
								aux = dedx_dat.y * dedxmid;
								temp = (fedep/(2-fedep));
								aux = aux*(1+2*aux)*temp*temp/6;
								tstep = fedep*ekei*dedxmid*(1+aux);
								//end of $COMPUTE-DRANGE(ekei,ekef,lelkef,elkei,elkef,tstep);
								tstep = tstep + tuss + range_ep[(medium*MXEKE+lelke-1)*2+qel]
								                     - range_ep[(medium*MXEKE+lelkef+1-1)*2+qel];
							}
						}
						total_tstep = tstep;
						compute_tstep = false;
					}
					tstep = total_tstep/rhof;  //non-default density scaling
					//end of $CALCULATE-TSTEP-FROM-DEMFP;
				}//end sig if-else

				//calculate stopping power
				/*
				float2 dedx_dat;
				if(lelec<0){  //"e-"
					//$EVALUATE dedx0 USING ededx(elke);->
					//float2 ededx_dat = ededx[medium * MXEKE + lelke-1];
					dedx_dat = ededx[medium * MXEKE + lelke-1];
					//dedx0 = ededx_dat.x + ededx_dat.y * elke;
				}
				else{         //"e+"
					//$EVALUATE dedx0 USING pdedx(elke);->
					//float2 pdedx_dat = pdedx[medium * MXEKE + lelke-1];
					dedx_dat = pdedx[medium * MXEKE + lelke-1];
					//dedx0 = pdedx_dat.x + pdedx_dat.y * elke;
				}
				*/
				float2 dedx_dat = dedx[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
				float dedx0 = dedx_dat.x + dedx_dat.y * elke;
				dedx_val = rhof * dedx0;

				//Determine maximum step-size (Formerly $SET-TUSTEP)
				float2 tmxs_dat = tmxs[medium * MXEKE + lelke-1];  //TMXS, used for maximum step-size interpolation
				float tmxs_val = tmxs_dat.x + tmxs_dat.y * elke;
				tmxs_val /= rhof;

				//Compute the range to E_min(medium) (e_min is the first energy in the table). Do not go
				//more than range. Don't replace this macro and don't override range, because the energy
				//loss evaluation below relies on the accurate (and self-consistent) evaluation of range!

				//The following macro computes the range to the minimum table energy. It uses $COMPUTE-DRANGE
				//Note that range_ep array is precomputed in subroutine mscati and gives the range
				//from the energy interval end points to AE for each medium.
				//$COMPUTE-RANGE;->
				float ekei = e_array[medium * MXEKE + lelke-1];
				float elkei = (lelke - eke01[medium].x)/eke01[medium].y;
				//$COMPUTE-DRANGE(eke_val,ekei,lelke,elke,elkei,range);->
				float fedep = 1 - ekei/eke_val;
				float elktmp = 0.5*(elke+elkei+0.25*fedep*fedep*(1+fedep*(1+0.875*fedep)));
				uint lelktmp = lelke;
				/*
				if(lelec<0){
					//$EVALUATE dedxmid USING ededx(elktmp);->
					//float2 ededx_dat = ededx[medium * MXEKE + lelktmp-1];
					dedx_dat = ededx[medium * MXEKE + lelktmp-1];
					//dedxmid = ededx_dat.x + ededx_dat.y * elktmp;
					//dedxmid = 1/dedxmid;
					//aux = ededx[medium * MXEKE + lelktmp-1].y * dedxmid;
				}
				else{
					//$EVALUATE dedxmid USING pdedx(elktmp);->
					//float2 pdedx_dat = pdedx[medium * MXEKE + lelktmp-1];
					dedx_dat = pdedx[medium * MXEKE + lelktmp-1];
					//dedxmid = pdedx_dat.x + pdedx_dat.y * elktmp;
					//dedxmid = 1/dedxmid;
					//aux = pdedx[medium * MXEKE + lelktmp-1].y * dedxmid;
				}
				*/
				dedx_dat = dedx[qel*MXEKE*MXMED + medium*MXEKE + lelktmp-1];
				float dedxmid = dedx_dat.x + dedx_dat.y * elktmp;
				dedxmid = 1/dedxmid;
				float aux = dedx_dat.y * dedxmid;
				float temp = (fedep/(2-fedep));
				aux = aux*(1+2*aux)*temp*temp/6;
				range = fedep*eke_val*dedxmid*(1+aux);
				//end of $COMPUTE-DRANGE(eke_val,ekei,lelke,elke,elkei,range);
				range = (range + range_ep[(medium * MXEKE + lelke-1) * 2 + qel]) / rhof;
				//end of $COMPUTE-RANGE;

				//The RANDOMIZE-TUSTEP option forced the electrons to approach discrete events
				//(Moller,brems etc.) only in a single scattering mode => waste of CPU time.
				bool random_tustep = false;  //switches tustep radomization off
				if(random_tustep){
					//float rnnotu = get_rand(idx);
					float rnnotu = curand_uniform(&randStat);
					tmxs_val = rnnotu*min(tmxs_val,smaxir[medium]);  //reg_dat.smaxir is smxir(irl)
					                                //check here and below, do we need -1 or not???
				}
				else
					tmxs_val = min(tmxs_val,smaxir[medium]);  //smaxir, "geom. step-size constrain for each region"
				tustep = min(tstep,min(tmxs_val,range));

				//optional tustep restriction in EM field
				//$SET-TUSTEP-EM-FIELD;

				// HOWNEAR
				//CALL-HOWNEAR(tperp);
				float tperp;
				/*
				if(p.leafIndex<-60 && p.leafIndex>=PASSMLC)
					hownear_phantom(p, tperp);
				else
					hownear_MLC(p, tperp);
				*/
				if(p.module==m_SecJawY)
					hownear_SecJawsY(p,tperp);
				else if(p.module==m_SecJawX)
					hownear_SecJawsX(p,tperp);
				else if(p.module==m_VarMLCs)
					hownear_MLC(p,tperp);
				else if(p.module==m_BlockMd)
					hownear_block(p,tperp);
				else if(p.module==m_WedgeMd)
					hownear_wedge(p,tperp);
				else if(p.module==m_Phantom)
					hownear_phantom(p,tperp);
				else{
					printf("error in module status!\n");
					p.status = e_empty;
					continue;
				}
				p.dnear = tperp;

				//optional regional range rejection for particles below e_max_rr if i_do_rr set

				//macro to do range rejection on a region by region basis if the user requests it.
				//The variables e_max_rr and i_do_rr are in COMIN ET-CONTROL. This macro is called
				//immediately after $USER-RANGE-DISCARD in ELECTR and everytime called the electrons
				//current range has been computed and stored in range and the distance to the nearest
				//boundary has just been computed and is in tperp. e_max_rr and i_do_rr
				//are initialized to zero in BLOCK DATA so range rejection is not done unless
				//Since option must be turned on by the user, it is considered a USER-ELECTRON-DISCARD.
				//Note this technique implies an approximation because the particle is not allowed to
				//create a brem particle which might escape the region. This is why e_max_rr is used,
				//to allow high energy electrons to be tracked in case they give off brem.
				//$RANGE-DISCARD;->
				if(i_do_rr[medium] == 1 && p.e < e_max_rr[medium]){    //check here, do we need reg_dat.med-1 or not
					if(tperp >= range){  //particle cannot escape local region
						edep1 = p.e - prm;
						escore(p.region,edep1*p.wt);
						if(lelec>0)
							annih_at_rest(p,reg_dat,&randStat);
						else
							p.status=e_empty;
						//return;
						continue;
					}
				}

				//This macro sets the minimum step size for a condensed history (CH) step. When the exact BCA is 
				//used, the minimum CH step is determined by efficiency consideratioins only. At about 3 elastic 
				//MFP's single scattering becomes more efficient than CH and so the algorithm switches off CH. 
				//If one of the various inexact BCA's is invoked, this macro provides a simple way to include 
				//more sophisticated decisions about the maximum acceptable approximated CH step.
				//The parameters passed to the macro in ELECTR are eke_val and elke
				//$SET-SKINDEPTH(eke_val,elke);->
				//This macro calculates the elastic scattering MFP. If spin_effects is false, the screened
				//Rutherford cross section is used, else the elastic MFP is based on PWA cross sections
				//$CALCULATE-ELASTIC-SCATTERING-MFP(ssmfp,eke_val,elke});->
				blccl = rhof * blcc[medium];
				xccl  = rhof * xcc[medium];
				float p2 = eke_val * (eke_val + rmt2);
				beta2 = p2 / (p2 + rmsq);
				if(spin_effects){
					/*
					float2 eta_ms_dat;
					if(lelec<0){
						//$EVALUATE etap USING etae_ms(elke);->
						//float2 etae_ms_dat = etae_ms[medium * MXEKE + lelke-1];
						eta_ms_dat = etae_ms[medium * MXEKE + lelke-1];
						//etap = etae_ms_dat.x + etae_ms_dat.y * elke;
						//etae_ms, "for interpolation of screening parameter (e-)"
					}
					else{
						//$EVALUATE etap USING etap_ms(elke);->
						//float2 etap_ms_dat = etap_ms[medium * MXEKE + lelke-1];
						eta_ms_dat = etap_ms[medium * MXEKE + lelke-1];
						//etap = etap_ms_dat.x + etap_ms_dat.y * elke;
						//etap_ms, "for interpolation of screening parameter (e+)"
					}
					*/
					float2 eta_ms_dat = eta_ms[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
					etap = eta_ms_dat.x + eta_ms_dat.y * elke;
					//$EVALUATE ms_corr USING blcce(elke);
					float2 blcce_dat = blcce[medium * MXEKE + lelke-1];
					float ms_corr = blcce_dat.x + blcce_dat.y * elke;
					blccl = blccl/etap/(1+0.25*etap*xccl/blccl/p2)*ms_corr;
					//blcce, "for interpolation of scattering power correction necessary to
					//account for scattering already taken into account in discrete Moller/Bhabha"
				}
				float ssmfp = beta2/blccl;
				//end of $CALCULATE-ELASTIC-SCATTERING-MFP(ssmfp,eke_val,elke});
				//float skindepth = SKIN_DEPTH_FOR_BCA * ssmfp;
				float skindepth = skindepth_for_bca * ssmfp;
				//end of $SET-SKINDEPTH(eke_val,elke);
				tustep = min(tustep, max(tperp, skindepth));

				//The transport logic below is determined by the logical variables: 
				//callhowfar, domultiple and dosingle
				// 
				//There are the following possibilities:
				// 
				//    callhowfar = .false.  This indicates that the intended step is
				//    ====================  shorter than tperp independent of BCA used
				//  - domultiple = .false. dosingle = .false. and callmsdist = .true.
				//        ==> everything has been done in msdist
				//  - domultiple = .true. and dosingle = .false.
				//        ==> should happen only if exact_bca = .false.
				//            indicates that MS remains to be done
				//  - domultiple = .false. and dosingle = .true.
				//        ==> should happen only if exact_bca = .true.
				//            sampled distance to a single scattering event is shorter than tperp
				//            ==> do single scattering at the end of the step
				//  - domultiple = .true. and dosingle = .true.
				//        ==> error condition, something with the logic is wrong!
				//
				//    callhowfar = .true.  This indicates that the intended step is longer than tperp
				//    ===================  and forces a call to howfar which returns the straight line
				//                         distance to the boundary in the intial direction of motion
				//                         (via a modification of ustep)
				//  - domultiple = .false. and dosingle = .false.
				//        ==> should happen only if exact_bca = .true.
				//            simply put the particle on the boundary
				//  - domultiple = .false. and dosingle = .true.
				//        ==> should happen only if exact_bca = .true.
				//            single elastic scattering has to be done
				//  - domultiple = .true. and dosingle = .false.
				//        ==> should happen only if exact_bca = .false.
				//            indicates that MS remains to be done
				//  - domultiple = .true. and dosingle = .true.
				//        ==> error condition, something with the logic is wrong!

				if(tustep<=tperp && (!exact_bca || tustep>skindepth)){
					//We are further way from a boundary than a skindepth, so perform a normal condensed-history step
					callhowfar = false;  //Do not call HOWFAR
					domultiple = false;  //Multiple scattering done here
					dosingle   = false;  //MS => mo single scattering
					callmsdist = true;   //Remember that msdist has been called

					//Fourth order technique for de
					//The following is a generalized version of $COMPUTE-ELOSS
					//$COMPUTE-ELOSS-G(tustep,eke_val,elke,lelke,de);->
					float tuss = range - range_ep[(medium * MXEKE + lelke-1) * 2 + qel]/rhof;
					//here tuss is the range between the initial energy
					//and the next lower energy on the interpolation grid
					if(tuss>=tustep){  //Final energy is in the same interpolation bin
						//The following macro computes the energy loss due to sub-threshold processes for a
						//path-length {P1}. The energy at the beginning of the step is {P2}, {P3}=Log({P2}),
						//{P4} is the interpolation index. The formulae are based on the logarithmic
						//interpolation for dedx used in EGSnrc. The result is returned in {P5}.
						//Assumes that initial and final energy are in the same interpolation bin.
						//$COMPUTE-ELOSS(tustep,eke_val,elke,lelke,de);->
						/*
						if(lelec<0){
							//float2 ededx_dat = ededx[medium * MXEKE + lelke-1];
							dedx_dat = ededx[medium * MXEKE + lelke-1];
							//dedxmid = ededx_dat.x + ededx_dat.y * elke;
							//aux = ededx[medium * MXEKE + lelke-1].y / dedxmid;
						}
						else{
							//float2 pdedx_dat = pdedx[medium * MXEKE + lelke-1];
							dedx_dat = pdedx[medium * MXEKE + lelke-1];
							//dedxmid = pdedx_dat.x + pdedx_dat.y * elke;
							//aux = pdedx[medium * MXEKE + lelke-1].y / dedxmid;
						}
						*/
						dedx_dat = dedx[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
						dedxmid = dedx_dat.x + dedx_dat.y * elke;
						aux = dedx_dat.y / dedxmid;
						de = dedxmid * tustep * rhof;
						fedep = de / eke_val;  //de = eke_val - Ei;  //fedep = eke_val/Ei - 1;
						de = de*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1-0.25*fedep*(2-aux*(4-aux)))));
						//my calculated result is as below
						//de = de*(1+0.5*fedep*aux*(1-0.333333*fedep*(2*aux+1-0.25*fedep*(2+6*aux*(1+aux)))));
						//end of $COMPUTE-ELOSS(tustep,eke_val,elke,lelke,de);
					}
					else{  //Must find first the table index where the step ends using pre-calculated ranges
						lelktmp = lelke;
						tuss = (range - tustep) * rhof;  //now tuss is the range of the final energy electron
						//scaled to the default mass density from PEGS4
						if(tuss<=0)
							de = eke_val - te[medium] * 0.99;
						//i.e., if the step we intend to take is longer than the particle range,
						//the particle energy goes down to the threshold, ({P2} is the initial
						//particle energy), originally the entire energy was lost,
						//but msdist_xxx is not prepared to deal with such large eloss fractions
						else{
							while(tuss<range_ep[(medium * MXEKE + lelktmp-1) * 2 + qel])
								lelktmp = lelktmp - 1;
							elktmp = (lelktmp + 1 - eke01[medium].x) / eke01[medium].y;
							float eketmp = e_array[medium * MXEKE + lelktmp + 1-1];
							tuss = (range_ep[(medium * MXEKE + lelktmp + 1-1) * 2 + qel] - tuss) / rhof;
							//$COMPUTE-ELOSS(tuss,eketmp,elktmp,lelktmp,de);->
							/*
							if(lelec<0){
								//float2 ededx_dat = ededx[medium * MXEKE + lelktmp-1];
								dedx_dat = ededx[medium * MXEKE + lelktmp-1];
								//dedxmid = ededx_dat.x + ededx_dat.y * elktmp;
								//aux = ededx[medium * MXEKE + lelktmp-1].y / dedxmid;
							}
							else{
								//float2 pdedx_dat = pdedx[medium * MXEKE + lelktmp-1];
								dedx_dat = pdedx[medium * MXEKE + lelktmp-1];
								//dedxmid = pdedx_dat.x + pdedx_dat.y * elktmp;
								//aux = pdedx[medium * MXEKE + lelktmp-1].y / dedxmid;
							}
							*/
							dedx_dat = dedx[qel*MXEKE*MXMED + medium*MXEKE + lelktmp-1];
							dedxmid = dedx_dat.x + dedx_dat.y * elktmp;
							aux = dedx_dat.y / dedxmid;
							de = dedxmid * tuss * rhof;
							fedep = de / eketmp;
							de = de*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1-0.25*fedep*(2-aux*(4-aux)))));
							//end of $COMPUTE-ELOSS(tuss,eketmp,elktmp,lelktmp,de);
							de = de + eke_val - eketmp;
						}
					}//end of $COMPUTE-ELOSS-G(tustep,eke_val,elke,lelke,de);

					tvstep = tustep;

					if(transport_algorithm==PRESTA_II)
						msdist_pII(de,tustep,spin_effects,p,msi,msj,omega2,ustep,reg_dat,&randStat);
					else
						msdist_pI(de,tustep, spin_effects,p,ustep,reg_dat,&randStat);
				}
				else{
					//We are within a skindepth from a boundary, invoke one of the various boundary-crossing algorithms
					callmsdist = false;  //Remember that msdist has not been called
					if(exact_bca){
						//Cross the boundary in a single scattering mode 
						domultiple = false;  //Do not do multiple scattering
						                     //Sample the distance to a single scattering event
						//float rnnoss = get_rand(idx);
						float rnnoss = curand_uniform(&randStat);
						float lambda = - logf(1 - rnnoss);
						float tmp =(eke_val/rm+1);
						float lambda_max = 0.5*blccl*rm/dedx_val*tmp*tmp*tmp;  //0.5 means decrease to half step of the initial ??
						if(lambda>=0 && lambda_max>0){
							float tuss;
							if(lambda<lambda_max)
								tuss = lambda*ssmfp*(1-0.5*lambda/lambda_max);
							else
								tuss = 0.5*lambda*ssmfp;  //looks like the average ??
							if(tuss<tustep){
								tustep = tuss;
								dosingle = true;
							}
							else
								dosingle = false;
						}
						else{
							//$egs_warning(*,' lambda > lambda_max: ',lambda,lambda_max,' eke dedx: ',eke,dedx,
							//    ' ir medium blcc: ',ir(np),medium,blcc(medium),' position = ',x(np),y(np),z(np));
							dosingle = false;
							p.status = e_empty;    //check here!!!!
							//return;
							continue;
						}
						ustep = tustep;
					}
					else{
						//Boundary crossing a la EGS4/PRESTA-I but using exact PLC(path-length-correction)
						dosingle = false;
						domultiple = true;
						//The following macros are used in subroutine electr in order to make
						//path length corrections and restrictions
						//$SET-USTEP;->
						ekems = eke_val - 0.5*tustep*dedx_val;  //Use mid-point energy to calculate
						                                    //energy dependent quantites
						//$CALCULATE-XI(tustep);->
						p2 = ekems*(ekems+rmt2);
						beta2 = p2/(p2 + rmsq);
						chia2 = xccl/(4*blccl*p2);  //note that our chia2 is Moliere chia2/4
						                            //note also that xcc is now old egs xcc**2
						xi = 0.5*xccl/p2/beta2*tustep;
						if(spin_effects){
							elkems = logf(ekems);
							lelkems = eke01[medium].x + eke01[medium].y * elkems;
							/*
							float2 eta_ms_dat,q1c_ms_dat;
							if(lelec<0){
								//$EVALUATE etap USING etae_ms(elkems);
								//float2 etae_ms_dat = etae_ms[medium * MXEKE + lelkems-1];
								eta_ms_dat = etae_ms[medium * MXEKE + lelkems-1];
								//etap = etae_ms_dat.x + etae_ms_dat.y * elkems;
								//$EVALUATE xi_corr USING q1ce_ms(elkems);
								//float2 q1ce_ms_dat = q1ce_ms[medium * MXEKE + lelkems-1];
								q1c_ms_dat = q1ce_ms[medium * MXEKE + lelkems-1];
								//xi_corr = q1ce_ms_dat.x + q1ce_ms_dat.y * elkems;
							}
							else{
								//$EVALUATE etap USING etap_ms(elkems);
								//float2 etap_ms_dat = etap_ms[medium * MXEKE + lelkems-1];
								eta_ms_dat = etap_ms[medium * MXEKE + lelkems-1];
								//etap = etap_ms_dat.x + etap_ms_dat.y * elkems;
								//$EVALUATE xi_corr USING q1cp_ms(elkems);
								//float2 q1cp_ms_dat = q1cp_ms[medium * MXEKE + lelkems-1];
								q1c_ms_dat = q1cp_ms[medium * MXEKE + lelkems-1];
								//xi_corr = q1cp_ms_dat.x + q1cp_ms_dat.y * elkems;
							}
							*/
							float2 eta_ms_dat = eta_ms[qel*MXEKE*MXMED + medium*MXEKE + lelkems-1];
							float2 q1c_ms_dat = q1c_ms[qel*MXEKE*MXMED + medium*MXEKE + lelkems-1];
							etap = eta_ms_dat.x + eta_ms_dat.y * elkems;
							xi_corr = q1c_ms_dat.x + q1c_ms_dat.y * elkems;
							chia2 = chia2*etap;
							xi = xi*xi_corr;
							//$EVALUATE ms_corr USING blcce(elkems);
							float2 blcce_dat = blcce[medium * MXEKE + lelkems-1];
							float ms_corr = blcce_dat.x + blcce_dat.y * elkems;
							blccl = blccl*ms_corr;
						}
						else{
							xi_corr = 1;
							etap = 1;
						}
						xi = xi*(logf(1+1.0/chia2)-1/(1+chia2));
						//end of $CALCULATE-XI(tustep);

						if(xi<0.1)
							ustep = tustep*(1 - xi*(0.5 - xi*0.166667));
						else
							ustep = tustep*(1 - expf(-xi))/xi;
						//end of $SET-USTEP;
					}

					if(ustep<tperp)
						callhowfar = false;
					else
						callhowfar = true;
				}
            }//end non-vacuum test

			uint  irnew = p.region;  //default new region is old region
			char  idisc = 0;         //default is no discard (this flag is initialized here)
			float ustep0 = ustep;    //Save the intended ustep

			//$CALL-HOWFAR-IN-ELECTR;->  //should be checked here
			//int ix,iy,iz;
			/*
			if(callhowfar || p.wt<=0){
				if(p.leafIndex<-60&& p.leafIndex >= PASSMLC)
					//irnew = howfar_phantom(p,ustep,ix,iy,iz);
					irnew = howfar_phantom(p,ustep);
				else
					irnew = howfar_MLC(p,ustep);
			}
			*/
			if(callhowfar || p.wt<=0){
				if(p.module==m_SecJawY)
					irnew = howfar_SecJawsY(p,ustep);
				else if(p.module==m_SecJawX)
					irnew = howfar_SecJawsX(p,ustep);
				else if(p.module==m_VarMLCs)
					irnew = howfar_MLC(p,ustep);
				else if(p.module==m_BlockMd)
					irnew = howfar_block(p,ustep);
				else if(p.module==m_WedgeMd)
					irnew = howfar_wedge(p,ustep);
				else if(p.module==m_Phantom)
					irnew = howfar_phantom(p,ustep);
				else{
					printf("error in module status!\n");
					p.status = e_empty;
					continue;
				}
			}

			//discard strategy is same as photon
			if (irnew == 0) {  //if going into VACUUM , commented by Tong Xu, Jan 2013
				if (reg_dat.med == VACUUM||p.leafIndex==PASSMLC)  // and current region is also VACUUM ,
					//or it has passed MLC and didn't hit the phantom commented by Tong Xu, Jan 2013
					idisc = 1;							   // discard ! , commented by Tong Xu, Jan 2013
				else
					idisc = -1;
			}

			//Now see if user requested discard
			if(idisc>0){
				p.region = 0;    //0 means VACUUM
				//p.status = e_user_discard;
				edep1 = p.e + p.charge*prm;
				escore(p.region,edep1*p.wt);
				p.status = e_empty;
				//return;
				continue;
			}

			//$CHECK-NEGATIVE-USTEP;->    //This is should be considered later
			if(ustep<=0){
				//Negative ustep--probable truncation problem at a boundary, which means we are not in region 
				//we think we are in. The default macro assumes that user has set irnew to the region we are 
				//really most likely to be in. A message is written out whenever ustep is less than -1.e-4
				/*
				if(ustep<-1e-4){
				    ierust = ierust + 1;
				    //OUTPUT ierust,ustep,dedx,e(np)-prm,ir(np),irnew,irold,x(np),y(np),z(np);
				    //    (i4,' Negative ustep = ',e12.5,' dedx=',F8.4,' ke=',F8.4,
				    //        ' ir,irnew,irold =',3i4,' x,y,z =',4e10.3);
				    if(ierust>1000){
					    //OUTPUT;(////' Called exit---too many ustep errors'///);
					    exit(1);  //$CALL_EXIT(1);
				    }
			    }
			    */
				ustep = 0;
			}//end of $CHECK-NEGATIVE-USTEP;

			if(ustep==0 || medium==VACUUM){  //corrected by Wu, 20130604
				//Do fast step
				if(ustep!=0){
					//Step in vacuum
					vstep = ustep;   //vstep is ustep truncated (possibly) by howfar
					tvstep = vstep;  //tvstep is the total curved path associated with vstep
					//edep = 0;        //edep = pzero;  no energy loss in vacuum
					//$VACUUM-TRANSPORT-EM-FIELD;  //additional vacuum transport in em field
					//e_range = VACUUM_STEP;  //e_range = vacdst;
					//$AUSCALL($TRANAUSB);

					//Transport the particle
					p.x = p.x + p.u * vstep;
					p.y = p.y + p.v * vstep;
					p.z = p.z + p.w * vstep;
					p.dnear = p.dnear - vstep;  //dnear is distance to the nearest boundary that goes along with 
					//particle stack and which the user's howfar can supply (option)
					irold = p.region;  //save previous region
					//$SET-ANGLES-EM-FIELD;  //allows for EM field deflection
				}//end of vacuum step

				//$electron_region_change;->
				p.region = irnew;
				irl = irnew;
				reg_dat = region_data[irl];
				/*
				if(p.region<=2)   //if the region belong to MLC
					reg_dat = region_data[irl];  //check here!!!!!!!!!!
				else  //if the particle is in the phantom, we use texture fetching, hopefully faster
				{
					int2 tmp1=tex3D(med_flags_tex, ix,iy,iz);
					reg_dat.med = tmp1.x;
					reg_dat.flags = tmp1.y;
					float4 tmp2 = tex3D(rhof_rho_pcut_ecut_tex,ix,iy,iz);
					reg_dat.rhof = tmp2.x;
					reg_dat.rho  = tmp2.y;
					reg_dat.pcut = tmp2.z;
					reg_dat.ecut = tmp2.w;
				}
				*/
				//medium = reg_dat.med;  //this is not necessary
				//rhof = reg_dat.rhof;

				//if(ustep!=0) $AUSCALL($TRANAUSA);
				if(eie<=reg_dat.ecut){
					//p.status = e_cutoff_discard;
					edep1 = p.e - prm;
					escore(p.region,edep1*p.wt);
					if(p.charge>0)
						annih_at_rest(p,reg_dat,&randStat);
					else
						p.status = e_empty;
					//return;
					continue;
				}
				if(ustep!=0 && idisc<0){
					//p.status = e_user_discard;
					edep1 = p.e + p.charge*prm;
					escore(p.region,edep1*p.wt);
					p.status = e_empty;
					//return;
					continue;
				}
				//p.status = e_electron_step;
				is_tstep = true;
				continue;
				//return;  //NEXT :TSTEP:  //start again at :STEP:
			}//Go try another big step in (possibly) new medium

			vstep = ustep;
			if(callhowfar){
				if(exact_bca){
					//if callhowfar = .true. and exact_bca = .true. => we are in a single scattering mode
					tvstep = vstep;
					if(tvstep!=tustep)
						//Boundary was crossed. Shut off single scattering
						dosingle = false;
				}
				else{
					//callhowfar = .true. and exact_bca = .false. =>we are doing an approximate CH step
					//calculate the average curved path-length corresponding to vstep
					//$SET-TVSTEP;->
					if(vstep<ustep0){
						ekems = eke_val - 0.5*tustep*vstep/ustep0*dedx_val;
						//This estimates the energy loss to the boundary.
						//tustep was the intended curved path-length,
						//ustep0 is the average transport distance in the initial direction
						//       resulting from tustep
						//vstep = ustep is the reduced average transport distance in the initial direction
						//              due to boundary crossing
						//$CALCULATE-XI(vstep);->
						float p2 = ekems*(ekems+rmt2);
						beta2 = p2/(p2 + rmsq);
						chia2 = xccl/(4*blccl*p2);  //Note that our chia2 is Moliere chia2/4
						//Note also that xcc is now old egs xcc**2
						xi = 0.5*xccl/p2/beta2*vstep;
						if(spin_effects){
							elkems = logf(ekems);
							//$SET INTERVAL elkems,eke;
							lelkems = eke01[medium].x + eke01[medium].y * elkems;
							/*
							float2 eta_ms_dat,q1c_ms_dat;
							if(lelec<0){
								//$EVALUATE etap USING etae_ms(elkems);
								//float2 etae_ms_dat = etae_ms[medium * MXEKE + lelkems-1];
								eta_ms_dat = etae_ms[medium * MXEKE + lelkems-1];
								//etap = etae_ms_dat.x + etae_ms_dat.y * elkems;
								//$EVALUATE xi_corr USING q1ce_ms(elkems);
								//float2 q1ce_ms_dat = q1ce_ms[medium * MXEKE + lelkems-1];
								q1c_ms_dat = q1ce_ms[medium * MXEKE + lelkems-1];
								//xi_corr = q1ce_ms_dat.x + q1ce_ms_dat.y * elkems;
							}
							else{
								//$EVALUATE etap USING etap_ms(elkems);
								//float2 etap_ms_dat = etap_ms[medium * MXEKE + lelkems-1];
								eta_ms_dat = etap_ms[medium * MXEKE + lelkems-1];
								//etap = etap_ms_dat.x + etap_ms_dat.y * elkems;
								//$EVALUATE xi_corr USING q1cp_ms(elkems);
								//float2 q1cp_ms_dat = q1cp_ms[medium * MXEKE + lelkems-1];
								q1c_ms_dat = q1cp_ms[medium * MXEKE + lelkems-1];
								//xi_corr = q1cp_ms_dat.x + q1cp_ms_dat.y * elkems;
							}
							*/
							float2 eta_ms_dat = eta_ms[qel*MXEKE*MXMED + medium*MXEKE + lelkems-1];
							float2 q1c_ms_dat = q1c_ms[qel*MXEKE*MXMED + medium*MXEKE + lelkems-1];
							etap = eta_ms_dat.x + eta_ms_dat.y * elkems;
							xi_corr = q1c_ms_dat.x + q1c_ms_dat.y * elkems;
							chia2 = chia2*etap;
							xi = xi*xi_corr;
							//$EVALUATE ms_corr USING blcce(elkems);
							float2 blcce_dat = blcce[medium * MXEKE + lelkems-1];
							float ms_corr = blcce_dat.x + blcce_dat.y * elkems;
							blccl = blccl*ms_corr;
						}
						else{
							xi_corr = 1;
							etap = 1;
						}
						xi = xi*(logf(1+1.0/chia2)-1/(1+chia2));
						//end of $CALCULATE-XI(vstep);

						if(xi<0.1)
							tvstep = vstep*(1 + xi*(0.5 + xi*0.333333));
						else{
							if(xi<0.999999)
								tvstep = - vstep*logf(1 - xi)/xi;
							else{
							    //This is an error condition because the average transition
							    //in the initial direction of motion is always smaller than 1/Q1
							    //$egs_info(*,' Stoped in SET-TVSTEP because xi > 1! ');
							    //$egs_info(*,' Medium: ',medium);
							    //$egs_info(*,' Initial energy: ',eke);
							    //$egs_info(*,' Average step energy: ',ekems);
							    //$egs_info(*,' tustep: ',tustep);
							    //$egs_info(*,' ustep0: ',ustep0);
							    //$egs_info(*,' vstep: ',vstep);
							    //$egs_info(*,' ==> xi = ',xi);
							    //$egs_fatal(*,'This is a fatal error condition');
								p.status = e_empty;
								//return;  //check here!!!!!!
								continue;
						    }
						}
					}
					else
						tvstep = tustep;
					//end of $SET-TVSTEP;
				}
				//Fourth order technique for dedx
				//Must be done for an approx. CH step or a single scattering step.

				//$COMPUTE-ELOSS-G(tvstep,eke_val,elke,lelke,de);
				float tuss = range - range_ep[(medium * MXEKE + lelke-1) * 2 + qel]/rhof;
				//here tuss is the range between the initial energy
				//and the next lower energy on the interpolation grid
				if(tuss>=tvstep){  //Final energy is in the same interpolation bin
					//$COMPUTE-ELOSS(tvstep,eke_val,elke,lelke,de);
					/*
					float2 dedx_dat;
					if(lelec<0){
						//float2 ededx_dat = ededx[medium * MXEKE + lelke-1];
						dedx_dat = ededx[medium * MXEKE + lelke-1];
						//dedxmid = ededx_dat.x + ededx_dat.y * elke;
						//aux = ededx[medium * MXEKE + lelke-1].y / dedxmid;
					}
					else{
						//float2 pdedx_dat = pdedx[medium * MXEKE + lelke-1];
						dedx_dat = pdedx[medium * MXEKE + lelke-1];
						//dedxmid = pdedx_dat.x + pdedx_dat.y * elke;
						//aux = pdedx[medium * MXEKE + lelke-1].y / dedxmid;
					}
					*/
					float2 dedx_dat = dedx[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
					float dedxmid = dedx_dat.x + dedx_dat.y * elke;
					float aux = dedx_dat.y / dedxmid;
					de = dedxmid * tvstep * rhof;
					float fedep = de / eke_val;
					de = de*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1-0.25*fedep*(2-aux*(4-aux)))));
					//end of $COMPUTE-ELOSS(tvstep,eke_val,elke,lelke,de);
				}
				else{  //Must find first the table index where the step ends using pre-calculated ranges
					uint lelktmp = lelke;
					tuss = (range - tvstep) * rhof;  //now tuss is the range of the final energy electron
					//scaled to the default mass density from PEGS4
					if(tuss<=0)
						de = eke_val - te[medium] * 0.99;
					//i.e., if the step we intend to take is longer than the particle range,
					//the particle energy goes down to the threshold, (eke_val is the initial particle energy)
					//originally the entire energy was lost, 
					//but msdist_xxx is not prepared to deal with such large eloss fractions
					else{
						while(tuss<range_ep[(medium * MXEKE + lelktmp-1) * 2 + qel])
							lelktmp = lelktmp -1;
						float elktmp = (lelktmp + 1 - eke01[medium].x) / eke01[medium].y;
						float eketmp = e_array[medium * MXEKE + lelktmp + 1-1];
						tuss = (range_ep[(medium * MXEKE + lelktmp + 1-1) * 2 + qel] - tuss) / rhof;
						//$COMPUTE-ELOSS(tuss,eketmp,elktmp,lelktmp,de);->
						/*
						float2 dedx_dat;
						if(lelec<0){
							//float2 ededx_dat = ededx[medium * MXEKE + lelktmp-1];
							dedx_dat = ededx[medium * MXEKE + lelktmp-1];
							//dedxmid = ededx_dat.x + ededx_dat.y * elktmp;
							//aux = ededx[medium * MXEKE + lelktmp-1].y / dedxmid;
						}
						else{
							//float2 pdedx_dat = pdedx[medium * MXEKE + lelktmp-1];
							dedx_dat = pdedx[medium * MXEKE + lelktmp-1];
							//dedxmid = pdedx_dat.x + pdedx_dat.y * elktmp;
							//aux = pdedx[medium * MXEKE + lelktmp-1].y / dedxmid;
						}
						*/
						float2 dedx_dat = dedx[qel*MXEKE*MXMED + medium*MXEKE + lelktmp-1];
						float dedxmid = dedx_dat.x + dedx_dat.y * elktmp;
						float aux = dedx_dat.y / dedxmid;
						de = dedxmid * tuss * rhof;
						float fedep = de / eketmp;
						de = de*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1-0.25*fedep*(2-aux*(4-aux)))));
						//end of $COMPUTE-ELOSS(tuss,eketmp,elktmp,lelktmp,de);
						de = de + eke_val - eketmp;
					}
				}//end of $COMPUTE-ELOSS-G(tvstep,eke_val,elke,lelke,de);
		    }
		    else{
				//callhowfar = .false. => step has not been reduced due to boundaries
				tvstep = tustep;
				if(!callmsdist){
					//Second order technique for dedx
					//Already done in a normal CH step with call to msdist

					//$COMPUTE-ELOSS-G(tvstep,eke_val,elke,lelke,de);
					float tuss = range - range_ep[(medium * MXEKE + lelke-1) * 2 + qel]/rhof;
					//here tuss is the range between the initial energy
					//and the next lower energy on the interpolation grid
					if(tuss>=tvstep){  //Final energy is in the same interpolation bin
						//$COMPUTE-ELOSS(tvstep,eke_val,elke,lelke,de);
						/*
						float2 dedx_dat;
						if(lelec<0){
							//float2 ededx_dat = ededx[medium * MXEKE + lelke-1];
							dedx_dat = ededx[medium * MXEKE + lelke-1];
							//dedxmid = ededx_dat.x + ededx_dat.y * elke;
							//aux = ededx[medium * MXEKE + lelke-1].y / dedxmid;
						}
						else{
							//float2 pdedx_dat = pdedx[medium * MXEKE + lelke-1];
							dedx_dat = pdedx[medium * MXEKE + lelke-1];
							//dedxmid = pdedx_dat.x + pdedx_dat.y * elke;
							//aux = pdedx[medium * MXEKE + lelke-1].y / dedxmid;
						}
						*/
						float2 dedx_dat = dedx[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
						float dedxmid = dedx_dat.x + dedx_dat.y * elke;
						float aux = dedx_dat.y / dedxmid;
						de = dedxmid * tvstep * rhof;
						float fedep = de / eke_val;
						de = de*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1-0.25*fedep*(2-aux*(4-aux)))));
						//end of $COMPUTE-ELOSS(tvstep,eke_val,elke,lelke,de);
					}
					else{  //Must find first the table index where the step ends using pre-calculated ranges
						uint lelktmp = lelke;
						tuss = (range - tvstep) * rhof;  //now tuss is the range of the final energy electron
						//scaled to the default mass density from PEGS4
						if(tuss<=0)
							de = eke_val - te[medium] * 0.99;
						//i.e., if the step we intend to take is longer than the particle range,
						//the particle energy goes down to the threshold, (eke_val is the initial particle energy)
						//originally the entire energy was lost, 
						//but msdist_xxx is not prepared to deal with such large eloss fractions
						else{
							while(tuss<range_ep[(medium * MXEKE + lelktmp-1) * 2 + qel])
								lelktmp = lelktmp -1;
							float elktmp = (lelktmp + 1 - eke01[medium].x) / eke01[medium].y;
							float eketmp = e_array[medium * MXEKE + lelktmp + 1-1];
							tuss = (range_ep[(medium * MXEKE + lelktmp + 1-1) * 2 + qel] - tuss) / rhof;
							//$COMPUTE-ELOSS(tuss,eketmp,elktmp,lelktmp,de);->
							/*
							float2 dedx_dat;
							if(lelec<0){
								//float2 ededx_dat = ededx[medium * MXEKE + lelktmp-1];
								dedx_dat = ededx[medium * MXEKE + lelktmp-1];
								//dedxmid = ededx_dat.x + ededx_dat.y * elktmp;
								//aux = ededx[medium * MXEKE + lelktmp-1].y / dedxmid;
							}
							else{
								//float2 pdedx_dat = pdedx[medium * MXEKE + lelktmp-1];
								dedx_dat = pdedx[medium * MXEKE + lelktmp-1];
								//dedxmid = pdedx_dat.x + pdedx_dat.y * elktmp;
								//aux = pdedx[medium * MXEKE + lelktmp-1].y / dedxmid;
							}
							*/
							float2 dedx_dat = dedx[qel*MXEKE*MXMED + medium*MXEKE + lelktmp-1];
							float dedxmid = dedx_dat.x + dedx_dat.y * elktmp;
							float aux = dedx_dat.y / dedxmid;
							de = dedxmid * tuss * rhof;
							float fedep = de / eketmp;
							de = de*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1-0.25*fedep*(2-aux*(4-aux)))));
							//end of $COMPUTE-ELOSS(tuss,eketmp,elktmp,lelktmp,de);
							de = de + eke_val - eketmp;
						}
					}//end of $COMPUTE-ELOSS-G(tvstep,eke_val,elke,lelke,de);
				}
			}

			//Calculates tvstep given vstep
			float save_de = de;  //the energy loss is used to calculate the number of MFP gone up to now.

			float edep = de;            //energy deposition variable for user
			float ekef = eke_val - de;  //(final kinetic energy)
			//eold = eie;           //save old value
			//enew = eold - de;     //energy at end of transport

			//Now do multiple scattering
			if(!callmsdist){  //everything done if callmsdist = .true.

				float  sinthe   = 0.0F;       //sin(THETA)
				float  costhe   = 0.0F;       //cos(THETA)

				if(domultiple){
					//Approximated CH step => do multiple scattering
					//ekems, elkems, beta2 have been set in either $SET-TUSTEP or $SET-TVSTEP
					//if spin_effects is true, they are not needed if spin_effects is false.
					//chia2, etap, xi, xi_corr are also set in the above macros.
					//qel (0 for e-, 1 for e+) and medium are now also required
					//(for the spin rejection loop)
					float lambda = blccl*tvstep/beta2/etap/(1+chia2);
					xi = xi/xi_corr;  //do not consider the q1 correction
					bool findindex = true;
					bool spin_index = true;
					mscat(lambda,chia2,xi,elkems,beta2,p,spin_effects,findindex,msi,msj,omega2,spin_index,costhe,sinthe,reg_dat,&randStat);
				}
				else{
					if(dosingle){
						//Single scattering
						ekems = max(ekef,reg_dat.ecut-rm);
						float p2 = ekems*(ekems + rmt2);
						beta2 = p2/(p2 + rmsq);
						chia2 = xcc[medium]/(4*blcc[medium]*p2);
						if(spin_effects){
							elkems = logf(ekems);
							//$SET INTERVAL elkems,eke;
							lelkems = eke01[medium].x + eke01[medium].y * elkems;
							/*
							float2 eta_ms_dat;
							if(lelec<0){
								//$EVALUATE etap USING etae_ms(elkems);
								//float2 etae_ms_dat = etae_ms[medium * MXEKE + lelkems-1];
								eta_ms_dat = etae_ms[medium * MXEKE + lelkems-1];
								//etap = etae_ms_dat.x + etae_ms_dat.y * elkems;
							}
							else{
								//$EVALUATE etap USING etap_ms(elkems);
								//float2 etap_ms_dat = etap_ms[medium * MXEKE + lelkems-1];
								eta_ms_dat = etap_ms[medium * MXEKE + lelkems-1];
								//etap = etap_ms_dat.x + etap_ms_dat.y * elkems;
							}
							*/
							float2 eta_ms_dat = eta_ms[qel*MXEKE*MXMED + medium*MXEKE + lelkems-1];
							etap = eta_ms_dat.x + eta_ms_dat.y * elkems;
							chia2 = chia2*etap;
						}
						sscat(chia2,elkems,beta2,p,spin_effects,costhe,sinthe,reg_dat,&randStat);
					}
					else{
						//theta = 0;  //No deflection in single scattering model
						sinthe = 0;
						costhe = 1;
					}
				}
			//}  //there are two if(!callmsdist) in succession, so we merge to one

			//We now know distance and amount of energy loss for this step,
			//and the angle by which the electron will be scattered. Hence,
			//it is time to call the user and inform him of this transport,
			//after which we will do it.

			//Now transport, deduct energy loss, and do multiple scatter.
			//e_range = range;

			//Put expected final position and direction in common block variables
			//so that they are available to the user for things such as scoring
			//on a grid that is different from the geometry grid
			//if(callmsdist){  //x,y,z,u,v,w are already changed in msdist
				//Deflection and scattering have been calculated/sampled in msdist
				//u_final = uscat;
				//v_final = vscat;
				//w_final = wscat;
				//x_final = xtrans;
				//y_final = ytrans;
				//z_final = ztrans;
			//}
			//else{  //here p.u, p.v and p.w are not changed
			//now we don't need if(callmsdist), so the else changed to if(!callmsdist)
			//if(!callmsdist){
				//x_final = p.x + p.u * vstep;
				//y_final = p.y + p.v * vstep;
				//z_final = p.z + p.w * vstep;
				p.x = p.x + p.u * vstep;
				p.y = p.y + p.v * vstep;
				p.z = p.z + p.w * vstep;
				if(domultiple || dosingle){
					//float u_tmp = p.u;
					//float v_tmp = p.v;
					//float w_tmp = p.w;
					float sinphi,cosphi;
					uphi21(costhe,sinthe,cosphi,sinphi,p,&randStat);
					//uphi(2,1);  //Apply the deflection, save call to uphi if
					              //no deflection in a single scattering mode
					//u_final = p.u;
					//v_final = p.v;
					//w_final = p.w;
					//p.u = u_tmp;
					//p.v = v_tmp;
					//p.w = w_tmp;
				}
				//else{  //stay unchanged
					//u_final = p.u;
					//v_final = p.v;
					//w_final = p.w;
				//}
			}

			//$AUSCALL($TRANAUSB);

			//Transport the particle
			//p.x = x_final;
			//p.y = y_final;
			//p.z = z_final;
			//p.u = u_final;
			//p.v = v_final;
			//p.w = w_final;

			p.dnear = p.dnear - vstep;
			irold = p.region;  //save previous region

			//Now done with multiple scattering, update energy and see if below cut
			//peie = peie - edep;
			//eie  = peie;
			//p.e  = peie;
			eie = eie - edep;
			p.e = eie;

			//After transport call to user scoring routine
			//Tong moved the escore for transportation edep here. July 2014
			// because the following ecut condition may result in end the current loop without score the energy
			escore(irold,edep*p.wt);

			if(irnew==irl && eie<=reg_dat.ecut){
				//p.status = e_cutoff_discard;
				edep1 = p.e - prm;
				escore(p.region,edep1*p.wt);
				if(p.charge>0)
					annih_at_rest(p,reg_dat,&randStat);
				else
					p.status = e_empty;
				//return;
				continue;
			}

			ushort medold = medium;
			if(medium!=VACUUM){  //corrected by Wu, 20130711  //if come to here, this will always be true
				//float ekeold = eke_val;
				eke_val = eie - rm;  //update kinetic energy
				elke = logf(eke_val);
				//$SET INTERVAL elke,eke;  //Get updated interval
				//float2 eke_dat = eke01[medium];
				//lelke = (int)(eke_dat.x + eke_dat.y * elke);
				lelke = (int)(eke01[medium].x + eke01[medium].y * elke);
			}

			if(irnew!=irold){
				//$electron_region_change;
				p.region = irnew;
				irl = irnew;
				reg_dat = region_data[irl];
				/*
				if(p.region<=2)   //if the region belong to MLC
					reg_dat = region_data[irl];  //check here!!!!!!!!!!
				else  //if the particle is in the phantom, we use texture fetching, hopefully faster
				{
					int2 tmp1=tex3D(med_flags_tex, ix,iy,iz);
					reg_dat.med = tmp1.x;
					reg_dat.flags = tmp1.y;
					float4 tmp2 = tex3D(rhof_rho_pcut_ecut_tex,ix,iy,iz);
					reg_dat.rhof = tmp2.x;
					reg_dat.rho  = tmp2.y;
					reg_dat.pcut = tmp2.z;
					reg_dat.ecut = tmp2.w;
				}
				*/
				medium = reg_dat.med;
				//rhof = reg_dat.rhof;  //This will affect the below $UPDATE-DEMFP;
				                        //and there the rhof should be the previous one, so we don't update the rhof
			}


			if(eie<=reg_dat.ecut){
				//p.status = e_cutoff_discard;
				edep1 = p.e - prm;
				escore(p.region,edep1*p.wt);
				if(p.charge>0)
					annih_at_rest(p,reg_dat,&randStat);
				else
					p.status = e_empty;
				//return;
				continue;
			}

			//Now check for deferred discard request.
			//May have been set by either howfar, or one of the transport ausgab calls
			if(idisc<0){
				//p.status = e_user_discard;
				edep1 = p.e + p.charge*prm;
				escore(p.region,edep1*p.wt);
				p.status = e_empty;
				//return;
				continue;
			}

			if(medium!=medold){
				is_tstep = true;
				continue;
				//return;
			}
			//NEXT :TSTEP:;

			//$UPDATE-DEMFP;->
			demfp = demfp - save_de*sig_val;
			//total_de = total_de - save_de;    //This may not be necessary, by Wu 20130417
			total_tstep = total_tstep - tvstep*rhof;
			if(total_tstep<1e-9) demfp = 0;
			//end of $UPDATE-DEMFP;

			if(demfp>=EPSEMFP){
				is_tstep = false;
				continue;
				//return;
			}
		//}while(demfp>=EPSEMFP);  //]UNTIL(demfp < $EPSEMFP);
		//end ustep loop

			//Compute final sigma to see if resample is needed.
			//this will take the energy variation of the sigma into account using the fictitious sigma method.

			//$EVALUATE-SIGF;->
			/*
			float2 sig_dat,dedx_dat;
			if(lelec<0){
				//$EVALUATE sigf USING esig(elke);
				//float2 esig_dat = esig[medium * MXEKE + lelke-1];
				sig_dat = esig[medium * MXEKE + lelke-1];
				//sigf = esig_dat.x + esig_dat.y * elke;
				//$EVALUATE dedx0 USING ededx(elke);
				//float2 ededx_dat = ededx[medium * MXEKE + lelke-1];
				dedx_dat = ededx[medium * MXEKE + lelke-1];
				//dedx0 = ededx_dat.x + ededx_dat.y * elke;
				//sigf /= dedx0;
			}
			else{
				//$EVALUATE sigf USING psig(elke);
				//float2 psig_dat = psig[medium * MXEKE + lelke-1];
				sig_dat = psig[medium * MXEKE + lelke-1];
				//sigf = psig_dat.x + psig_dat.y * elke;
				//$EVALUATE dedx0 USING pdedx(elke);
				//float2 pdedx_dat = pdedx[medium * MXEKE + lelke-1];
				dedx_dat = pdedx[medium * MXEKE + lelke-1];
				//dedx0 = pdedx_dat.x + pdedx_dat.y * elke;
				//sigf /= dedx0;
			}
			*/
			float2 sig_dat = sig[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
			float2 dedx_dat = dedx[qel*MXEKE*MXMED + medium*MXEKE + lelke-1];
			float sigf = sig_dat.x + sig_dat.y * elke;
			float dedx0 = dedx_dat.x + dedx_dat.y * elke;
			sigf /= dedx0;
			//end of $EVALUATE-SIGF;

			float sigratio = sigf/sig0;

			//float rfict = get_rand(idx);  //$RANDOMSET rfict;
			float rfict = curand_uniform(&randStat);

			if(rfict>sigratio){
				is_tstep = true;
				continue;
				//return;
			}
		//]UNTIL (rfict <= sigratio) ;
		//}//end tstep loop

			//Now sample electron interaction

			if(lelec<0){
				//e-, check branching ratio
				//$EVALUATE-EBREM-FRACTION;->
				//$EVALUATE ebr1 USING ebr1(elke);->
				float2 ebr1_dat = ebr1[medium * MXEKE + lelke-1];
				float ebr1_val = ebr1_dat.x + ebr1_dat.y * elke;
				//float rnno24 = get_rand(idx);  //$RANDOMSET rnno24;
				float rnno24 = curand_uniform(&randStat);
				if(rnno24<=ebr1_val){
					//It was bremsstrahlung
					//go to :EBREMS:;->
					//call brems;
					p.status = e_brems;
					push_stack(p,e_brems);
					p.status = e_empty;
					continue;
					//return;
				}
				else{
					//It was Moller, but first check the kinematics.
					//However, if EII is on, we should still permit an interaction even if E<moller threshold
					//as EII interaction go down to the ionization threshold which may be less than thmoll.
					if(p.e<=thmoll[medium] && eii_flag==0){  //corrected by Wu, 20130603
						//(thmoll = lower Moller threshold)
						//Not enough energy for Moller, so force it to be a bremsstrahlung---provided ok kinematically.
						if(ebr1_val<=0){
							//go to :NEWELECTRON:;
							is_tstep = true;  //Brems not allowed either
							continue;
							//return;
						}
						//go to :EBREMS:;->
						//$AUSCALL($BREMAUSB);
						p.status = e_brems;
						push_stack(p,e_brems);
						p.status = e_empty;
						continue;
						//return;
					}
					//$AUSCALL($MOLLAUSB);
					p.status = e_moller;
					push_stack(p,e_moller);
					p.status = e_empty;
					continue;
					//return;
				}
				//go to :NEWELECTRON:; //Electron is lowest energy-follow it
				//after calling electron interaction function, set particle status to be e_new_particle
			}

			//e+ interaction. pbr1 = brems/(brems + bhabha + annih
			//$EVALUATE-PBREM-FRACTION;->
			//$EVALUATE pbr1 USING pbr1(elke);->
			float2 pbr1_dat = pbr1[medium * MXEKE + lelke-1];
			float pbr1_val = pbr1_dat.x + pbr1_dat.y * elke;
			//float rnno25 = get_rand(idx);
			float rnno25 = curand_uniform(&randStat);
			if(rnno25<pbr1_val){
				//It was bremsstrahlung
				//go to :EBREMS:;->
				//$AUSCALL($BREMAUSB);
				p.status = e_brems;
				push_stack(p,e_brems);
				p.status = e_empty;
				continue;
				//return;
			}
			//Decide between bhabha and annihilation
			//pbr2 is (brems + bhabha)/(brems + bhabha + annih)
			//$EVALUATE-BHABHA-FRACTION;->
			//$EVALUATE pbr2 USING pbr2(elke);
			float2 pbr2_dat = pbr2[medium * MXEKE + lelke-1];
			float pbr2_val = pbr2_dat.x + pbr2_dat.y * elke;
			if(rnno25<pbr2_val){
				//It is bhabha
				//$AUSCALL($BHABAUSB);
				p.status = e_bhabha;
				push_stack(p,e_bhabha);
				p.status = e_empty;
				continue;
				//return;
			}
			else{
				//It is in-flight annihilation
				//$AUSCALL($ANNIHFAUSB);
				p.status = e_annih;
				push_stack(p,e_annih);
				p.status = e_empty;
				continue;
				//return;
			}

		//] REPEAT  //newelectron
		//after calling electron interaction function, set particle status to be e_new_particle


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
