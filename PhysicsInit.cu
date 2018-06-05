//#include "CUDA_EGS.h"
//#include "EGSnrc_Parameters.h"
//#include "commonBlocks host.h"

//extern "C" void logstr(const char *format, ...);

 int egs_swap_4 (char c[4])
 {
	 char tmp;
	 tmp=c[3];
	 c[3]=c[0];
	 tmp=c[2];
	 c[2]=c[1];
	 c[1]=tmp;
	 return 1 ;
 }
  
 int egs_swap_2 (char c[2])
 {
	 char tmp;
	 tmp=c[1];
	 c[1]=c[0];
	 c[0]=tmp;
	 return 1;
 }

 int set_spline(double *x,double *f,double *a,double *b,double *c,double *d,int n)
 {
	 double r,s;
	 int m, m1,m2,mr;
	 m1=2;
	 m2=n-1;
	 s=0;
	 for (m=0;m<m2;m++)
	 {
		 d[m]=x[m+1]-x[m];
		 r=(f[m+1]-f[m])/d[m];
		 c[m]=r-s;
		 s=r;
	 }
	 s=0;
	 r=0;
	 c[0]=0;
	 c[n-1]=0;
	 for(m=m1-1;m<m2;m++)
	 {
		 c[m]+=r*c[m-1];
		 b[m]=2*(x[m-1]-x[m+1])-r*s;
		 s=d[m];
		 r=s/b[m];
	 }
	 mr=m2;
	 for(m=m1-1;m<m2;m++)
	 {
		 c[mr]=(d[mr]*c[mr+1]-c[mr])/b[mr];
		 mr=mr-1;
	 }
	 for (m=0;m<m2;m++)
	 {
		 s=d[m];
		 r=c[m+1]-c[m];
		 d[m]=r/s;
		 c[m]=3*c[m];
		 b[m]=(f[m+1]-f[m])/s-(c[m]+r)*s;
		 a[m]=f[m];
	 }
	 return 0;
 }
 double spline(double  s, double *x, double *a,double *b,double *c,double *d,int  n)
 {
	 //int i_1;
	 double ret_val;
	 int m;
	 double q;
	 int direction, ml, mu, mav, m_lower, m_upper;
	 if(x[0]>x[n-1])
	 {
		 direction = 1;
        m_lower = n;
        m_upper = 0;
	 }
	 else 
	 {
		 direction = 0;
        m_lower = 0;
        m_upper = n;
	 }
	 if(s>=x[m_upper+direction-1])
	 {
		  m = m_upper + 2*direction - 1;
	 }
	 else if(s<=x[m_lower-direction])
	 {
		 m = m_lower - 2*direction + 1;
	 }
	 else 
	 {
		 ml=m_lower;
		 mu=m_upper;
L5651:
		 if(abs(mu-ml)<=1)
		 {
			 goto L5652;
		 }
		 mav=(ml+mu)/2;
		 if(s<x[mav-1])
		 {
			 mu=mav;
		 }
		 else
		 {ml=mav;}
		 goto L5651;
 L5652:
          m=mu+direction-1;
	 }
	 q=s-x[m-1];
	 ret_val = a[m-1] + q * (b[m-1] + q * (c[m-1] + q * d[m-1]));
    return ret_val;
	
 }

 double erf1(double x)
{
	//int n;
	int k;
	double ret_val;
	double a[46]={ 1.0954712997776232,
	    -.289175401126989,.1104563986337951,-.0412531882278565,
	    .0140828380706516,-.0043292954474314,.0011982719015923,
	    -2.999729623532e-4,6.83258603789e-5,-1.42469884549e-5,
	    2.7354087728e-6,-4.861912872e-7,8.03872762e-8,-1.24184183e-8,
	    1.7995326e-9,-2.454795e-10,3.16251e-11,-3.859e-12,4.472e-13,
	    -4.93e-14,5.2e-15,-5e-16,1e-16,.9750834237085559,
	    -.0240493938504146,8.204522408804e-4,-4.34293081303e-5,
	    3.018447034e-6,-2.544733193e-7,2.4858353e-8,-2.7317201e-9,
	    3.308472e-10,1.464e-13,-2.44e-14,4.2e-15,-8e-16,1e-16,0.,0.,0.,0.,
	    0.,0.,0.,0.,0. };
/*	
	{ 1.0954712997776232 , -0.2891754011269890 , 0.1104563986337951 ,
		-0.0412531882278565 , 0.0140828380706516 , -0.0043292954474314 ,
		0.0011982719015923 , -0.0002999729623532 , 0.0000683258603789
      , -0.0000142469884549 , 0.0000027354087728 , -0.0000004861912872
     , 0.0000000803872762 , -0.0000000124184183 , 0.0000000017995326 ,
     -0.0000000002454795 , 0.0000000000316251 , -0.0000000000038590 , 0.0000000000004472 ,
	 -0.0000000000000493 , 0.0000000000000052 , -0.0000000000000005 ,
	 0.0000000000000001 , 0.9750834237085559 , -0.0240493938504146 ,
	 0.0008204522408804 , -0.0000434293081303 , 0.0000030184470340 ,
	 -0.0000002544733193 , 0.0000000248583530 , -0.0000000027317201 ,
	 0.0000000003308472 , 0.0000000000001464 , -0.0000000000000244 ,
	 0.0000000000000042 , -0.0000000000000008 , 0.0000000000000001 ,
	 0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};
	 */
	int nlim[2]={22,16};
	double cons1,bn,bn1,bn2,y,fac;
	cons1=1.128379167095513;
	if(x>3)
	{
		y=3.0/x;
		k=2;

	}
	else
	{
		y=x/3.0;
		k=1;
	}
	fac=2.0*(2.0*y*y-1.0);
	bn1=0.0f;
	bn=0.0f;
	for(int n=nlim[k-1];n>=0;n--)
	{
		bn2=bn1;
		bn1=bn;
		bn=fac*bn1-bn2+a[23*(k-1)+n];
	}
	if(k==1)
	{
		ret_val=cons1*y*(bn-bn1)-1.0;  //move the "-1" from init_compton here to reduce the round-off error
	}
	else
	{
		ret_val=-cons1*exp(-(x*x))*(bn-bn2+a[23*(k-1)])/(4.0*x);  // was "1.0-cons*exp..., now remove the "1" and removed the "-1" in the init_compton
	}
	return ret_val;
}

int egs_set_defaults()
{
	//h_vacdst=1e8;
	/*	
	for(int i=0;i<MXREG;i++)
	{
		h_ecut[i]=0.7f;  // was 0.0f
		h_pcut[i]=0.01f;  //was 0.0f
		h_ibcmp[i] = 1;
        h_iedgfl[i] = 1;
        h_iphter[i] = 1;
        h_smaxir[i] = 1e10;
        h_i_do_rr[i] = 0;
        h_e_max_rr[i] = 0;
        h_med[i] = 1;
        h_rhor[i] = 0;
        h_iraylr[i] = 0;
	}
	*/

	/*
	eii_xfile = 'Off'
	xsec_out = 0
	photon_xsections = 'si'
	comp_xsections = 'default'
	*/

	for(int i=0;i<MXMED;i++)
	{
		//h_iraylm[i]=0;
	
		h_ae[i]=0;
		h_ap[i]=0;
		h_ue[i]=0;
		h_te[i]=0;
		h_thmoll[i]=0;

	}
	for(int i=0;i<MXSHELL;i++)
	{
		for(int j=0;j<MXELEMENT;j++)
		{
			h_binding_energies[MXSHELL*j+i]=0;
		}
	}

#ifdef EII_FLAG
	h_eii_flag=1;
#else
	h_eii_flag=0;
#endif

#ifdef BAS_KM
	h_ibrdst = IBRDST_KM;
#else
	h_ibrdst = IBRDST_SIMPLE;
#endif

	//h_ibr_nist = 0;
	//h_pair_nrc = 0;
	//h_itriplet = 0;
	h_iprdst = IPRDST;

	h_dunit=1.0;
	//h_rhof = 1;
	/*
	for (int i=0;i<5;i++)
	{
		h_iausfl[i]=0;
	}
	for(int i=5;i<29;i++)
	{
		h_iausfl[i]=0;
	}
	*/

	h_ximax = 0.5;
	h_estepe = 0.25;
	//h_skindepth_for_bca = 3;  //no use, reinitialized by next line
	h_skindepth_for_bca = SKIN_DEPTH_FOR_BCA;
#ifdef ESA_PII
	h_transport_algorithm = PRESTA_II;
#else
	h_transport_algorithm = PRESTA__I;
#endif
#ifdef BCA_EXACT
	h_bca_algorithm = 0;
#else
	h_bca_algorithm = 1;
#endif
	h_exact_bca = true;  //no use, will be determined by h_bca_algorithm
#ifdef SPIN_FLAG
	h_spin_effects = true;
#else
	h_spin_effects = false;
#endif
	//h_count_pII_steps = 0;
	//h_count_all_steps = 0;
	//h_radc_flag = 0;
	h_nmed = 1;
	//h_rng_seed = 999999;
	//h_latchi = 0;
	h_rm = (float)0.5110034;
	h_prm = (score_t)0.5110034;
	h_prmt2 = 2*h_prm;
	h_rmt2 = 2*h_rm;
	h_rmsq = h_rm*h_rm;
	//h_pi = (float)3.1415926;
	//h_pi = 4*datan(1d0);              //???????
	//h_twopi = 2*h_pi;
	//h_pi5d2 = 2.5*h_pi;

	h_nbr_split = 1;

	//h_i_play_rr = 0;
	//h_i_survived_rr = 0;
	//h_prob_rr = -1;
	//h_n_rr_warning = 0;
	/*
	double P=1.0;   
	for( int i=0;i<50;i++)
	{
		h_pwr2i[i]=P;
		P=P/2.0;
	}
	*/

	return 1;
}



 int init_spin()
{
	double fine=137.03604;
	double tf_constant=0.88534138;
	//int i_1,i_2;
	//double  r_1,r_2,r_3;
	//double d_1;
	//short equiv_0[1];
	//float equiv_1[1];
	double e;
	int i,j,k;
	double z,eta_array[2*MAXE_SPI1_PLUS_1];
	//char spin_file[256];
	//int spin_unit;
	double af[MAXE_SPI1_PLUS_1],bf[MAXE_SPI1_PLUS_1],cf[MAXE_SPI1_PLUS_1],df[MAXE_SPI1_PLUS_1];
    int egs_swap_4(char *c), egs_swap_2(char *c);
	int je,iq;
	double z23;
	int iz;
	//int rec_length;
	float fmax_array[MAXQ_SPIN_PLUS_1];
	int ii4;
	double  aae,g_m,g_r,eta,eil;
	int n_q;
	double  sig,tau,tmp,si1e,si2e;
	double si1p,si2p;
	int neke;
	double dedx,etap,tauc,fmax1;
	int leil,irec;
	bool swap=false;
	//extern /* Subroutine */ int exit_(integer *);
	double beta2;
	int i_ele;
	double gamma;
	int ndata;
	double dloge,eloge,sum_a,sum_z,sum_z2;
	int n_ener,medium;
	double earray[MAXE_SPI1_PLUS_1],farray[MAXE_SPI1_PLUS_1];
	double sum_pz,c_array[2*MAXE_SPI1_PLUS_1],elarray[MAXE_SPI1_PLUS_1],g_array[2*MAXE_SPI1_PLUS_1];
	//double tmp_array[MAXE_SPI1_PLUS_1],f_array[2*MAXE_SPI1_PLUS_1];
	int n_point;
	short i2_array[512];
	char data_version[32],endianess[5];
	float tmp_4;
	float dum1,dum2,dum3,aux_o;

    ndata=MAXE_SPI1+1;
	
	// 前面有一堆的输出 字符  暂且省略
	string file1 = ("\\spinms.data");
	string spin_file1 = string(data_dir) + file1;
	FILE *i_spinms = fopen(spin_file1.c_str(),"rb");
	if(i_spinms==NULL)
	{
		logstr("spin data file not found: %s\n",spin_file1);
		//exit(1);
		return 1;
	}

	fread(data_version,sizeof(char),32,i_spinms);
	fread(endianess,sizeof(char),4,i_spinms);
	fread(&h_espin_min,sizeof(float),1,i_spinms);
	fread(&h_espin_max,sizeof(float),1,i_spinms);
	fread(&h_b2spin_min,sizeof(float),1,i_spinms);
	fread(&h_b2spin_max,sizeof(float),1,i_spinms);
	endianess[4]=0; //make sure the char array has a proper ending before convert to string.
	string endianess_str (endianess);
	if(endianess_str !="1234")swap = true;
	if(swap)
	{

		tmp_4 = h_espin_min;
        egs_swap_4((char*)&tmp_4);
        h_espin_min = tmp_4;
        tmp_4 = h_espin_max;
        egs_swap_4((char*)&tmp_4);
        h_espin_max = tmp_4;
        tmp_4 = h_b2spin_min;
        egs_swap_4((char*)&tmp_4);
        h_b2spin_min = tmp_4;
        tmp_4 = h_b2spin_max;
        egs_swap_4((char*)&tmp_4);
        h_b2spin_max = tmp_4;

	}

	n_ener=MAXE_SPIN;
	n_q=MAXQ_SPIN;
	n_point=MAXU_SPIN;
	dloge=log(h_espin_max/h_espin_min)/n_ener;
	eloge=log(h_espin_min);
	earray[0]=h_espin_min;
	for (i=1;i<=n_ener;i++)
	{
		eloge=eloge+dloge;
		earray[i]=exp(eloge);
	}
	h_dbeta2=(h_b2spin_max-h_b2spin_min)/n_ener;
	beta2=h_b2spin_min;
	earray[n_ener+1]=h_espin_max;
	for (i=n_ener+2;i<=2*n_ener+1;i++)
	{
		beta2=beta2+h_dbeta2;
		if(beta2<0.999f)
		{
			earray[i]=511.0034*(1/sqrtf(1-beta2)-1);

		}
		else 
		{
			earray[i]=50585.1;
		}

	}
	h_espin_min = h_espin_min/1000;
	h_espin_max = h_espin_max/1000;
	h_dlener = log(h_espin_max/h_espin_min)/MAXE_SPIN;
	h_dleneri = 1/h_dlener;
	h_espml = log(h_espin_min);
	h_dbeta2 = (h_b2spin_max-h_b2spin_min)/MAXE_SPIN;
	h_dbeta2i = 1/h_dbeta2;
	h_dqq1 = 0.5/n_q;
	h_dqq1i = 1/h_dqq1;


	for (medium=0;medium<h_nmed;medium++)
	{
		logstr("ini_spin for medium # %d\n" ,medium);

		for (iq=0;iq<=1;iq++)
		{
			for (i=0;i<=MAXE_SPI1;i++)
			{	 
				eta_array[i*2+iq]=0;
				c_array[i*2+iq]=0;
				g_array[i*2+iq]=0;

				for(j=0;j<=MAXE_SPI1;j++)
				{
					for(k=0;k<=MAXU_SPIN;k++)
					{

						h_spin_rej[MXMED*2*MAXE_SPI1_PLUS_1*MAXQ_SPIN_PLUS_1*k+j*MXMED*2*MAXE_SPI1_PLUS_1+i*MXMED*2+iq*MXMED+medium]=0;
					}
				}
			}     

		}
		sum_z2=0;
		sum_a=0;
		sum_pz=0;
		sum_z=0;
		for (i_ele=0;i_ele<h_nne[medium];i_ele++)
		{
			z=h_zelem[i_ele*MXMED+medium];
			iz=int(z+0.5f);
			tmp=h_pz[i_ele*MXMED+medium]*z*(z+1);
			sum_z2=sum_z2+tmp;
			sum_z+=h_pz[i_ele*MXMED+medium]*z;
			sum_a+=h_pz[i_ele*MXMED+medium]*h_wa[i_ele*MXMED+medium];
			sum_pz+=h_pz[i_ele*MXMED+medium];
			z23=pow(z,(double)0.66666666666667);
			for(iq=0;iq<=1;iq++)
			{
				for (i=0;i<=MAXE_SPI1;i++)
				{
					irec=1+(iz-1)*4*(n_ener+1)+2*iq*(n_ener+1)+i+1; 

					//we are going to read the irec-th record, so we skip (i rec-1) records first. Tong Xu
					dum1=0;dum2=0;dum3=0;aux_o=0;
					int item_read;
					fseek(i_spinms,(irec-1)*sizeof(float)*276,0);
					item_read=fread(&dum1,sizeof(float),1,i_spinms); if(item_read==0){logstr("error while reading spinms.dat, exiting...\n");exit(1);}
					item_read=fread (&dum2,sizeof(float),1,i_spinms);if(item_read==0){logstr("error while reading spinms.dat, exiting...\n");exit(1);}
					item_read=fread(&dum3,sizeof(float),1,i_spinms);if(item_read==0){logstr("error while reading spinms.dat, exiting...\n");exit(1);}
					item_read=fread(&aux_o,sizeof(float),1,i_spinms);if(item_read==0){logstr("error while reading spinms.dat, exiting...\n");exit(1);}
					item_read=fread(fmax_array,sizeof(float),MAXQ_SPIN_PLUS_1,i_spinms);if(item_read!=MAXQ_SPIN_PLUS_1){logstr("error while reading spinms.dat,  fmax_array,exiting...\n");exit(1);}
					item_read=fread(i2_array,sizeof(short),512,i_spinms);if(item_read!=512){logstr("error while reading spinms.dat, i2_array,exiting...\n");exit(1);}

					if(swap)
					{
						tmp_4 = dum1;
						egs_swap_4((char*)&tmp_4);
						dum1 = tmp_4;
						tmp_4 = dum2;
						egs_swap_4((char*)&tmp_4);
						dum2 = tmp_4;
						tmp_4 = dum3;
						egs_swap_4((char*)&tmp_4);
						dum3 = tmp_4;
						tmp_4 = aux_o;
						egs_swap_4((char*)&tmp_4);
						aux_o = tmp_4;
					}
					eta_array[2*i+iq]+=tmp*log(z23*aux_o);
					tau=earray[i]/511.0034;
					beta2=tau*(tau+2)/((tau+1)*(tau+1));
					eta=z23/((fine*tf_constant)*(fine*tf_constant))*aux_o/4/tau/(tau+2);
					c_array[2*i+iq]+=tmp*(log(1+1/eta)-1/(1+eta))*dum1*dum3;
					g_array[2*i+iq]+=tmp*dum2;
					for (j=0;j<=MAXQ_SPIN;j++)
					{
						tmp_4=fmax_array[j]; 
						if(swap) egs_swap_4((char*)&tmp_4);
						for (k=0;k<=MAXU_SPIN;k++)
						{
							short ii2=i2_array[(n_point+1)*j+k];           //remove the "+1" because the C++ arrays start from zero
							if(swap)
							{
								egs_swap_2((char*)&ii2); //!!!!
							}
							ii4=ii2;
							if(ii4<0)
							{
								ii4+=65536;
							}
							dum1=ii4;
							dum1=dum1*tmp_4/65535;
							h_spin_rej[MXMED*2*MAXE_SPI1_PLUS_1*MAXQ_SPIN_PLUS_1*k+j*MXMED*2*MAXE_SPI1_PLUS_1+i*MXMED*2+iq*MXMED+medium]+=tmp*dum1;
						}

					}

				}
			}
		}
		for(iq=0;iq<=1;iq++)
		{
			for(i=0;i<=MAXE_SPI1;i++)
			{

				for(j=0;j<=MAXQ_SPIN;j++)
				{
					fmax1=0;
					for(k=0;k<=MAXU_SPIN;k++)
					{
						if(h_spin_rej[MXMED*2*MAXE_SPI1_PLUS_1*MAXQ_SPIN_PLUS_1*k+j*MXMED*2*MAXE_SPI1_PLUS_1+i*MXMED*2+iq*MXMED+medium]>fmax1)
						{
							fmax1=h_spin_rej[MXMED*2*MAXE_SPI1_PLUS_1*MAXQ_SPIN_PLUS_1*k+j*MXMED*2*MAXE_SPI1_PLUS_1+i*MXMED*2+iq*MXMED+medium];
						}
					}
					//if(medium==4&&iq==0)
					//	logstr("fmax, i,j %f, %d, %d \n",fmax1,i,j);
					for(k=0;k<=MAXU_SPIN;k++)
					{
						h_spin_rej[MXMED*2*MAXE_SPI1_PLUS_1*MAXQ_SPIN_PLUS_1*k+j*MXMED*2*MAXE_SPI1_PLUS_1+i*MXMED*2+iq*MXMED+medium]=
							h_spin_rej[MXMED*2*MAXE_SPI1_PLUS_1*MAXQ_SPIN_PLUS_1*k+j*MXMED*2*MAXE_SPI1_PLUS_1+i*MXMED*2+iq*MXMED+medium]/fmax1;
					}

				}
			}
		}
		for (i=0;i<=MAXE_SPI1;i++)
		{
			tau=earray[i]/511.0034;
			beta2=tau*(tau+2)/((tau+1)*(tau+1));
			for(iq=0;iq<=1;iq++)
			{
				aux_o=exp(eta_array[2*i+iq]/sum_z2)/((fine*tf_constant)*(fine*tf_constant));
				eta_array[2*i+iq]=0.26112447*aux_o*h_blcc[medium]/h_xcc[medium];
				eta=aux_o/4/tau/(tau+2);
				gamma=3*(1+eta)*(log(1+1/eta)*(1+2*eta)-2)/(log(1+1/eta)*(1+eta)-1);
				g_array[i*2+iq]=g_array[i*2+iq]/sum_z2/gamma;
				c_array[i*2+iq]=c_array[i*2+iq]/sum_z2/(log(1+1/eta)-1/(1+eta));
			}
		}
		/* //the following was used for debuging
		if(medium==3)
		{
			 FILE *iSum=fopen("sum_for_init_spin.txt","wt");
			 FILE *iDump;
			 double sum;
			 DUMP_1D_REAL(init_spin,eta_array,2*MAXE_SPI1_PLUS_1);
			 DUMP_1D_REAL(init_spin,g_array,2*MAXE_SPI1_PLUS_1);
			 DUMP_1D_REAL(init_spin,c_array,2*MAXE_SPI1_PLUS_1);
			 fclose(iSum);
		}
		*/
		eil=(1-h_eke01[medium].x)/h_eke01[medium].y;
		e=exp(eil);
		if(e<=h_espin_min)
		{
			si1e=eta_array[0];
			si1p=eta_array[1];
		}
		else 
		{  
			if(e<=h_espin_max)
			{
				aae=(eil-h_espml)*h_dleneri;
				je=aae;
				aae-=je;
			}
			else
			{
				tau=e/0.5110034;
				beta2=tau*(tau+2)/((tau+1)*(tau+1));
				aae=(beta2-h_b2spin_min)*h_dbeta2i;
				je=aae;
				aae-=je;
				je+=MAXE_SPIN+1;

			}
			si1e=(1-aae)*eta_array[je*2]+aae*eta_array[(je+1)*2];
			si1p=(1-aae)*eta_array[je*2+1]+aae*eta_array[(je+1)*2+1];
		}
		neke=h_meke[medium];
		for(i=0;i<(neke-1);i++)
		{
			eil=(i+2-h_eke01[medium].x)/h_eke01[medium].y;
			e=exp(eil);
			if(e<=h_espin_min)
			{
				si2e=eta_array[0];
				si2p=eta_array[1];
			}
			else 
			{
				if(e<=h_espin_max)
				{
					aae=(eil-h_espml)*h_dleneri;
					je=aae;
					aae-=je;
				}
				else
				{
					tau=e/0.5110034;
					beta2=tau*(tau+2)/((tau+1)*(tau+1));
					aae=(beta2-h_b2spin_min)*h_dbeta2i;
					je = aae;
					aae = aae - je;
					je = je + MAXE_SPIN + 1;
				}
				si2e=(1-aae)*eta_array[je*2]+aae*eta_array[(je+1)*2];
				si2p=(1-aae)*eta_array[je*2+1]+aae*eta_array[(je+1)*2+1];
			}
			//h_etae_ms[medium*MXEKE+i].y=(si2e-si1e)*h_eke01[medium].y;
			h_eta_ms[medium*MXEKE+i].y=(si2e-si1e)*h_eke01[medium].y;
			//h_etae_ms[medium*MXEKE+i].x=si2e-h_etae_ms[medium*MXEKE+i].y*eil;
			h_eta_ms[medium*MXEKE+i].x=si2e-h_eta_ms[medium*MXEKE+i].y*eil;
			//h_etap_ms[medium*MXEKE+i].y=(si2p-si1p)*h_eke01[medium].y;
			h_eta_ms[MXEKE*MXMED+medium*MXEKE+i].y=(si2p-si1p)*h_eke01[medium].y;
			//h_etap_ms[medium*MXEKE+i].x=si2p-h_etap_ms[medium*MXEKE+i].y*eil;
			h_eta_ms[MXEKE*MXMED+medium*MXEKE+i].x=si2p-h_eta_ms[MXEKE*MXMED+medium*MXEKE+i].y*eil;
			si1e=si2e;
			si1p=si2p;

		}
		//h_etae_ms[medium*MXEKE+neke-1].y=h_etae_ms[medium*MXEKE+neke-2].y;
		h_eta_ms[medium*MXEKE+neke-1].y=h_eta_ms[medium*MXEKE+neke-2].y;
		//h_etae_ms[medium*MXEKE+neke-1].x=h_etae_ms[medium*MXEKE+neke-2].x;
		h_eta_ms[medium*MXEKE+neke-1].x=h_eta_ms[medium*MXEKE+neke-2].x;
		//h_etap_ms[medium*MXEKE+neke-1].y=h_etap_ms[medium*MXEKE+neke-2].y;
		h_eta_ms[MXEKE*MXMED+medium*MXEKE+neke-1].y=h_eta_ms[MXEKE*MXMED+medium*MXEKE+neke-2].y;
		//h_etap_ms[medium*MXEKE+neke-1].x=h_etap_ms[medium*MXEKE+neke-2].x;       //!!!!!!!! 改过
		h_eta_ms[MXEKE*MXMED+medium*MXEKE+neke-1].x=h_eta_ms[MXEKE*MXMED+medium*MXEKE+neke-2].x;
		for(i=0;i<=MAXE_SPIN;i++)
		{
			elarray[i]=log(earray[i]/1000);
			farray[i]=c_array[i*2];
		}
		for(i=MAXE_SPIN+1;i<=MAXE_SPI1-1;i++)
		{
			elarray[i]=log(earray[i+1]/1000);
			farray[i]=c_array[(i+1)*2];

		}
		if(h_ue[medium]>1e5)
		{
			elarray[ndata-1]=log(h_ue[medium]);
		}
		else
		{
			elarray[ndata-1]=log(1e5);
		}
		farray[ndata-1]=1;

		set_spline(elarray,farray,af,bf,cf,df,ndata)  ;       //           一个外部函数

		eil=(1-h_eke01[medium].x)/h_eke01[medium].y;	 
		si1e=spline(eil,elarray,af,bf,cf,df,ndata);

		for (i=0;i<neke-1;i++)
		{
			eil=(i+2-h_eke01[medium].x)/h_eke01[medium].y;
			si2e = spline(eil,elarray,af,bf,cf,df,ndata)   ;          //         spline is a external function  
			//h_q1ce_ms[MXEKE*medium+i].y=(si2e-si1e)*h_eke01[medium].y;
			h_q1c_ms[MXEKE*medium+i].y=(si2e-si1e)*h_eke01[medium].y;
			//h_q1ce_ms[MXEKE*medium+i].x=si2e-h_q1ce_ms[MXEKE*medium+i].y*eil;
			h_q1c_ms[MXEKE*medium+i].x=si2e-h_q1c_ms[MXEKE*medium+i].y*eil;
			si1e=si2e;
		}
		//h_q1ce_ms[medium*MXEKE+neke-1].y=h_q1ce_ms[medium*MXEKE+neke-2].y;
		h_q1c_ms[medium*MXEKE+neke-1].y=h_q1c_ms[medium*MXEKE+neke-2].y;
		//h_q1ce_ms[medium*MXEKE+neke-1].x=h_q1ce_ms[medium*MXEKE+neke-2].x;
		h_q1c_ms[medium*MXEKE+neke-1].x=h_q1c_ms[medium*MXEKE+neke-2].x;
		for(i=0;i<=MAXE_SPIN;i++)
		{
			farray[i]=c_array[2*i+1];
		}
		for(i=(MAXE_SPIN+1);i<=(MAXE_SPI1-1);i++)
		{
			farray[i]=c_array[2*(i+1)+1];
		}
		set_spline(elarray,farray,af,bf,cf,df,ndata)  ;                    
		eil=(1-h_eke01[medium].x)/h_eke01[medium].y;
		si1e = spline(eil,elarray,af,bf,cf,df,ndata);
		for(i=0;i<neke-1;i++)
		{
			eil=(i+2-h_eke01[medium].x)/h_eke01[medium].y;
			si2e = spline(eil,elarray,af,bf,cf,df,ndata);
			//h_q1cp_ms[medium*MXEKE+i].y=(si2e-si1e)*h_eke01[medium].y;
			h_q1c_ms[MXEKE*MXMED+medium*MXEKE+i].y=(si2e-si1e)*h_eke01[medium].y;
			//h_q1cp_ms[medium*MXEKE+i].x=si2e-h_q1cp_ms[medium*MXEKE+i].y*eil;
			h_q1c_ms[MXEKE*MXMED+medium*MXEKE+i].x=si2e-h_q1c_ms[MXEKE*MXMED+medium*MXEKE+i].y*eil;
			si1e=si2e;
		}
		//h_q1cp_ms[medium*MXEKE+neke-1].y=h_q1cp_ms[medium*MXEKE+neke-2].y;
		h_q1c_ms[MXEKE*MXMED+medium*MXEKE+neke-1].y=h_q1c_ms[MXEKE*MXMED+medium*MXEKE+neke-2].y;
		//h_q1cp_ms[medium*MXEKE+neke-1].x=h_q1cp_ms[medium*MXEKE+neke-2].x;
		h_q1c_ms[MXEKE*MXMED+medium*MXEKE+neke-1].x=h_q1c_ms[MXEKE*MXMED+medium*MXEKE+neke-2].x;
		for(i=0;i<=MAXE_SPIN;i++)
		{
			farray[i]=g_array[i*2];              //!!! 改过
		}
		for(i=(MAXE_SPIN+1);i<=(MAXE_SPI1-1);i++)
		{
			farray[i]=g_array[2*(i+1)];
		}
		set_spline(elarray,farray,af,bf,cf,df,ndata);
		eil=(1-h_eke01[medium].x)/h_eke01[medium].y;
		si1e=spline(eil,elarray,af,bf,cf,df,ndata);
		logstr( " init_spin Interpolation table for q2 correction (e-)");
		for (i=0;i<neke-1;i++)
		{
			eil=(i+2-h_eke01[medium].x)/h_eke01[medium].y;
			si2e = spline(eil,elarray,af,bf,cf,df,ndata);
			//h_q2ce_ms[MXEKE*medium+i].y=(si2e-si1e)*h_eke01[medium].y;
			h_q2c_ms[MXEKE*medium+i].y=(si2e-si1e)*h_eke01[medium].y;
			//h_q2ce_ms[medium*MXEKE+i].x=si2e-h_q2ce_ms[medium*MXEKE+i].y*eil;
			h_q2c_ms[medium*MXEKE+i].x=si2e-h_q2c_ms[medium*MXEKE+i].y*eil;
			si1e=si2e;
		}
		//h_q2ce_ms[medium*MXEKE+neke-1].y=h_q2ce_ms[medium*MXEKE+neke-2].y;
		h_q2c_ms[medium*MXEKE+neke-1].y=h_q2c_ms[medium*MXEKE+neke-2].y;
		//h_q2ce_ms[medium*MXEKE+neke-1].x=h_q2ce_ms[medium*MXEKE+neke-2].x;
		h_q2c_ms[medium*MXEKE+neke-1].x=h_q2c_ms[medium*MXEKE+neke-2].x;
		for(i=0;i<=MAXE_SPIN;i++)
		{
			farray[i]=g_array[2*i+1];
		}
		for(i=(MAXE_SPIN+1);i<=(MAXE_SPI1-1);i++)
		{
			farray[i]=g_array[(i+1)*2+1];
		}
		set_spline(elarray,farray,af,bf,cf,df,ndata)  ;        //  !!!!!!!!!!!
		eil=(1-h_eke01[medium].x)/h_eke01[medium].y;
		si1e = spline(eil,elarray,af,bf,cf,df,ndata);
		for(i=0;i<(neke-1);i++)
		{
			eil=(i+2-h_eke01[medium].x)/h_eke01[medium].y;
			si2e = spline(eil,elarray,af,bf,cf,df,ndata);
			//h_q2cp_ms[medium*MXEKE+i].y=(si2e-si1e)*h_eke01[medium].y;
			h_q2c_ms[MXEKE*MXMED+medium*MXEKE+i].y=(si2e-si1e)*h_eke01[medium].y;
			//h_q2cp_ms[medium*MXEKE+i].x=si2e-h_q2cp_ms[medium*MXEKE+i].y*eil;
			h_q2c_ms[MXEKE*MXMED+medium*MXEKE+i].x=si2e-h_q2c_ms[MXEKE*MXMED+medium*MXEKE+i].y*eil;
			si1e=si2e;
		}
		//h_q2cp_ms[medium*MXEKE+neke-1].y=h_q2cp_ms[medium*MXEKE+neke-2].y;
		h_q2c_ms[MXEKE*MXMED+medium*MXEKE+neke-1].y=h_q2c_ms[MXEKE*MXMED+medium*MXEKE+neke-2].y;
		//h_q2cp_ms[medium*MXEKE+neke-1].x=h_q2cp_ms[medium*MXEKE+neke-2].x;
		h_q2c_ms[MXEKE*MXMED+medium*MXEKE+neke-1].x=h_q2c_ms[MXEKE*MXMED+medium*MXEKE+neke-2].x;
		tauc=h_te[medium]/0.5110034;
		si1e=1;
		for(i=0;i<neke-1;i++)
		{
			eil=(i+2-h_eke01[medium].x)/h_eke01[medium].y;
			e=exp(eil);
			leil=i+1+1;
			tau=e/0.5110034;
			if(tau>2*tauc)
			{
				//sig=h_esig[medium*MXEKE+leil-1].y*eil+h_esig[medium*MXEKE+leil-1].x;
				sig=h_sig[medium*MXEKE+leil-1].y*eil+h_sig[medium*MXEKE+leil-1].x;
				//dedx=h_ededx[medium*MXEKE+leil-1].y*eil+h_ededx[medium*MXEKE+leil-1].x;
				dedx=h_dedx[medium*MXEKE+leil-1].y*eil+h_dedx[medium*MXEKE+leil-1].x;
				sig=sig/dedx;
				if(sig>1e-6)
				{
					//etap=h_etae_ms[medium*MXEKE+leil-1].y*eil+h_etae_ms[medium*MXEKE+leil-1].x;
					etap=h_eta_ms[medium*MXEKE+leil-1].y*eil+h_eta_ms[medium*MXEKE+leil-1].x;
					eta=0.25*etap*h_xcc[medium]/h_blcc[medium]/tau/(tau+2);
					g_r=(1+2*eta)*log(1+1/eta)-2;
					g_m=log(0.5*tau/tauc)+ (1+((tau+2)/(tau+1))*((tau+2)/(tau+1)))*log(2*(tau
						-tauc+2)/(tau+4))- 0.25*(tau+2)*(tau+2+2*(2*tau+1)/((tau+
						1)*(tau+1)))* log((tau+4)*(tau-tauc)/tau/(tau-tauc+2))+ 0.5*(tau
						-2*tauc)*(tau+2)*(1/(tau-tauc)-1/((tau+1)*(tau+1)));
					if(g_m<g_r)
					{
						g_m=g_m/g_r;
					}
					else
					{
						g_m=1;
					}
					si2e=1-g_m*sum_z/sum_z2;
				}
				else
				{
					si2e=1;
				}
			}
			else
			{
				si2e=1;
			}

			h_blcce[medium*MXEKE+i].y=(si2e-si1e)*h_eke01[medium].y;
			h_blcce[medium*MXEKE+i].x=si2e-h_blcce[medium*MXEKE+i].y*eil;
			si1e=si2e;
		}
		h_blcce[medium*MXEKE+neke-1].y=h_blcce[medium*MXEKE+neke-2].y;
		h_blcce[medium*MXEKE+neke-1].x=h_blcce[medium*MXEKE+neke-2].x;
		//write(i_log,'(a)') ' done'                             
    }
	fclose(i_spinms);


	//write(i_log,'(/a)') '***************** Error: '
	//write(i_log,'(a,a)') 'Failed to open spin data file ',spin_file(:l
	//*nblnk1(spin_file))
	//write(i_log,'(/a)') '***************** Quiting now.'
	//call exit(1)
	//write(i_log,'(/a)') '***************** Error: '
	//write(i_log,*) 'Error while reading spin data file for element',iZ
	//write(i_log,'(/a)') '***************** Quiting now.'      
	return 1;
}



 int init_ms_sr()
{
	int i,j,k1,k2,k;
	//logstr( "Reading screened Rutherford MS data ......" );
	/*    *......... ' */
	string file1 = ("\\msnew.data");
	string ms_sr_file1 = string(data_dir) + file1; 
	FILE *i_mscat = fopen(ms_sr_file1.c_str(),"r");
	if(i_mscat==NULL)
	{
		logstr("MS scatter data file not found: %s\n",ms_sr_file1);
		//exit(1);
		return 1;
	}
	int simsarray=0;
	int buf;
	for (i=0;i<=MAXL_MS;i++)
	{
		for (j=0;j<=MAXQ_MS;j++)
		{
			for (k1=0;k1<=MAXU_MS;k1++)	fscanf (i_mscat,"%f",&h_ums_array[k1*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1+j*MAXL_MS_PLUS_1+i]);
			for (k1=0;k1<=MAXU_MS;k1++)	fscanf (i_mscat,"%f",&h_fms_array[k1*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1+j*MAXL_MS_PLUS_1+i]);
			for (k2=0;k2<=MAXU_MS-1;k2++)	fscanf (i_mscat,"%f",&h_wms_array[k2*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1+j*MAXL_MS_PLUS_1+i]);
			for (k2=0;k2<=MAXU_MS-1;k2++)
				{
					fscanf (i_mscat,"%d",&buf);   //Becareful!, h_ims_array is a short integer array, so we can't read it directly with "%d" formate, which require a int pointer. 
					h_ims_array[k2*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1+j*MAXL_MS_PLUS_1+i]=buf;
				}
			for (k=0;k<=MAXU_MS-1;k++)
			{
				h_fms_array[k*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1+j*MAXL_MS_PLUS_1+i]=h_fms_array[(k+1)*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1+j*MAXL_MS_PLUS_1+i]/
				h_fms_array[k*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1+j*MAXL_MS_PLUS_1+i] -1;
				h_ims_array[k*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1+j*MAXL_MS_PLUS_1+i]-=1;
				simsarray+=h_ims_array[k*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1+j*MAXL_MS_PLUS_1+i];
			}
			h_fms_array[MAXU_MS*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1+j*MAXL_MS_PLUS_1+i]=h_fms_array[(MAXU_MS-1)*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1+j*MAXL_MS_PLUS_1+i];
		}
	}
	//logstr("Done\n"); //, sum of ims_array = %d, [64]=%d,[128]=%d\n",simsarray,h_ims_array[64],h_ims_array[128]);
	h_llammin=log(LAMBMIN_MS);
	h_llammax=log(LAMBMAX_MS);
	h_dllamb=(h_llammax-h_llammin)/MAXL_MS;
	h_dllambi = 1.0f/h_dllamb;
	h_dqms = QMAX_MS/MAXQ_MS;
	h_dqmsi = 1.0f/h_dqms;
	return 1;
}

 int mscati(void)
{
	int i_1,i_2,i;
	float r_1,sige_old,sigp_old,p2,ei,si,rm,eke,eil,sig,aux,eip1,sip1;
	//float r_2;
	float ekef,elke;
	int neke,leil;
	float beta2,dedx0,eip1l,elkef,ededx;
	int lelke;
	float sigee,sigep,ektmp,chi_a2,tstbm;
	int leip1l,lelkef,medium;
	float ecutmn,elktmp,tstbmn,estepx;
	int lelktmp;
	bool ise_monoton,isp_monoton;
	
	if (h_bca_algorithm==0)
	{
		h_exact_bca=true;
	}
	else
	{
		h_exact_bca=false;
	}
	if((h_estepe<=0.f)||(h_estepe>=1.0f)) 
		{
			h_estepe=MAX_ELOSS;
	}
	if((h_ximax<=0.f)||(h_ximax>=1.0f)) 
	{
		if (h_exact_bca==true)
		{
			h_ximax=EXACT_BCA_XIMAX;
		}
		else
		{
			h_ximax=INEXACT_BCA_XIMAX;
		}
	}
	if(h_transport_algorithm!=PRESTA_II&&h_transport_algorithm!=PRESTA__I&&h_transport_algorithm!=VMC)
	{
		h_transport_algorithm=PRESTA_II;
	}
	if(h_skindepth_for_bca<=(1e-4f))
	{
		if(!h_exact_bca)
		{
			//write(i_log,*) ' old PRESTA calculates default min. step-size
			//*for BCA: '
			ecutmn=1e30f;
			for( i=0;i<total_number_regions;i++)
			{
				/*
				if (h_med[i]>0&&h_med[i]<=h_nmed) 
				{
					ecutmn=min(ecutmn, h_ecut[i]);
				}
				*/
				if( h_region_data[i].med>=0&& h_region_data[i].med<h_nmed)
					ecutmn = min(ecutmn, h_region_data[i].ecut);
			}
			//write(i_log,*) '     minimum ECUT found: ',ecutmn
			tstbmn=1e30f;
			i_1=h_nmed;
			for (medium=0;medium<i_1;medium++)
			{
				r_1=ecutmn;
				tstbm=(ecutmn-0.5110034f)*(ecutmn+0.5110034f)/(r_1*r_1);
				r_1=ecutmn/h_xcc[medium];
				tstbm=h_blcc[medium]*tstbm*(r_1*r_1);
				aux=log(tstbm);
				if(aux>300.f)
				{
					//write(i_log,*) 'aux > 300 ? ',aux	
				}
				tstbm=log(tstbm/aux);
				tstbmn=min(tstbmn,tstbm);


			}
			//write(i_log,*) '     default BLCMIN is: ',tstbmn
			h_skindepth_for_bca=exp(tstbmn);
			/*write(i_log,*) '     this corresponds to ',skindepth_for_bca,
			*' elastic MFPs '*/

		}
		else 
		{
			h_skindepth_for_bca=SKIN_DEPTH_FOR_BCA;
		}
	}

	init_ms_sr();                                                 //    / an undefined function
    i_1=h_nmed;
	for (medium=0;medium<i_1;medium++)
	{
		h_blcc[medium]*=1.16699413758864573f;
		r_1=h_xcc[medium];
		h_xcc[medium]=r_1*r_1;
	}
	if (h_spin_effects==true)
	{
		init_spin();                                     // an undefined function;
	}                                          
	//write(i_log,*) ' '
	h_esige_max=0.0f;
	h_psige_max=0.0f;
	i_1=h_nmed;
	for(medium=0;medium<i_1;medium++)
	{
		sigee=1e-15f;
		sigep=1e-15f;
		neke=h_meke[medium];
		ise_monoton=true;
		isp_monoton=true;
		sige_old=-1.0f;
		sigp_old=-1.0f;
		i_2=neke;
		for (i=1;i<=i_2;i++)
		{
			ei=exp(((float)i-h_eke01[medium].x)/h_eke01[medium].y);
			eil=log(ei);
			leil=i;
			//ededx=h_ededx[leil+medium*MXEKE-1].y*eil+h_ededx[leil+medium*MXEKE-1].x;
			ededx=h_dedx[leil+medium*MXEKE-1].y*eil+h_dedx[leil+medium*MXEKE-1].x;
			//sig=h_esig[leil+medium*MXEKE-1].y*eil+h_esig[leil+medium*MXEKE-1].x;
			sig=h_sig[leil+medium*MXEKE-1].y*eil+h_sig[leil+medium*MXEKE-1].x;
			sig/=ededx;
			if (sig>sigee)
			{
				sigee=sig;
			}
			if(sig<sige_old)
			{ 
				ise_monoton=false;
			}
			sige_old=sig;
			//ededx=h_pdedx[leil+medium*MXEKE-1].y*eil+h_pdedx[leil+medium*MXEKE-1].x;
			ededx=h_dedx[MXEKE*MXMED+leil+medium*MXEKE-1].y*eil+h_dedx[MXEKE*MXMED+leil+medium*MXEKE-1].x;
			//sig=h_psig[leil+medium*MXEKE-1].y*eil+h_psig[leil+medium*MXEKE-1].x;
			sig=h_sig[MXEKE*MXMED+leil+medium*MXEKE-1].y*eil+h_sig[MXEKE*MXMED+leil+medium*MXEKE-1].x;
			sig/=ededx;
			if (sig>sigep)
			{ 
				sigep=sig;
			}
			if(sig<sigp_old)
			{
				isp_monoton=false;
			}
			sigp_old=sig;
		}
		/*write(i_log,*) ' Medium ',medium,' sige = ',sigee,sigep,' monoto
		*ne = ', ise_monoton,isp_monoton*/
		h_sig_ismonotone[2*medium]=ise_monoton;             
		h_sig_ismonotone[2*medium+1]=isp_monoton;           
		//h_esig_e[medium]=sigee;
		h_sig_e[medium]=sigee;
		//h_psig_e[medium]=sigep;
		h_sig_e[MXMED+medium]=sigep;
		if (sigee>h_esige_max)
		{
			h_esige_max=sigee;
		}
		if(sigep>h_psige_max)
		{
			h_psige_max=sigep;
		}

	}
	//write(i_log,*) ' '
	//write(i_log,*) ' Initializing tmxs for estepe = ',estepe,' and xim
	//*ax = ',ximax
	//write(i_log,*) ' '
	rm=0.5110034f;
	i_1=h_nmed;
	for (medium=1;medium<=i_1;medium++)
	{
		ei=exp((1-h_eke01[medium-1].x)/h_eke01[medium-1].y);
		eil=log(ei);
		leil=1;
		h_e_array[(medium-1)*MXEKE]=ei;
		h_expeke1[medium-1]=exp(1.0f/h_eke01[medium-1].y)-1;
		h_range_ep[(medium-1)*2*MXEKE]=0;
		h_range_ep[(medium-1)*2*MXEKE+1]=0;
		neke=h_meke[medium-1];
		i_2=neke-1;
		for (i=1;i<=i_2;i++)
		{
			eip1=exp(((float)(i+1)-h_eke01[medium-1].x)/h_eke01[medium-1].y);
			h_e_array[(medium-1)*MXEKE+i]=eip1;
			eke=(eip1+ei)*0.5f;
			elke=log(eke);
			lelke=h_eke01[medium-1].y*elke+h_eke01[medium-1].x;
			//ededx=h_pdedx[(medium-1)*MXEKE+lelke-1].y*elke+h_pdedx[(medium-1)*MXEKE+lelke-1].x;
			ededx=h_dedx[MXEKE*MXMED+(medium-1)*MXEKE+lelke-1].y*elke+h_dedx[MXEKE*MXMED+(medium-1)*MXEKE+lelke-1].x;
			//aux=h_pdedx[(medium-1)*MXEKE+i-1].y/ededx;
			aux=h_dedx[MXEKE*MXMED+(medium-1)*MXEKE+i-1].y/ededx;
			r_1=(eip1-ei)/eke;
			h_range_ep[(medium-1)*2*MXEKE+i*2+1]=h_range_ep[(medium-1)*2*MXEKE+(i-1)*2+1]+(eip1-ei)/ededx*(1+aux*(aux*2+1)*(r_1*r_1)/24);//???
			//ededx=h_ededx[(medium-1)*MXEKE+lelke-1].y*elke+h_ededx[lelke+(medium-1)*MXEKE-1].x;
			ededx=h_dedx[(medium-1)*MXEKE+lelke-1].y*elke+h_dedx[lelke+(medium-1)*MXEKE-1].x;
			//aux=h_ededx[i+(medium-1)*MXEKE-1].y/ededx;
			aux=h_dedx[i+(medium-1)*MXEKE-1].y/ededx;
			r_1=(eip1-ei)/eke;
			h_range_ep[(medium-1)*2*MXEKE+i*2]=h_range_ep[(medium-1)*2*MXEKE+(i-1)*2]+(eip1-ei)/ededx*(aux*(aux*2+1)*(r_1*r_1)/24+1);
			ei=eip1;
        }
		eil=(1-h_eke01[medium-1].x)/h_eke01[medium-1].y;
		ei=exp(eil);
		leil=1;
		p2=ei*(ei+rm*2);
		beta2=p2/(p2+rm*rm);
		chi_a2=h_xcc[medium-1]/(p2*4*h_blcc[medium-1]);
		//dedx0=h_ededx[leil+(medium-1)*MXEKE-1].y*eil+h_ededx[leil+(medium-1)*MXEKE-1].x;
		dedx0=h_dedx[leil+(medium-1)*MXEKE-1].y*eil+h_dedx[leil+(medium-1)*MXEKE-1].x;
		estepx=p2*2*beta2*dedx0/ei/h_xcc[medium-1]/(log(1.0f/chi_a2+1)*(chi_a2+1)-1);
		estepx*=h_ximax;
		if(estepx>h_estepe)
		{
			estepx=h_estepe;
		}
		si=estepx*ei/dedx0;
		i_2=neke-1;
		for (i=1;i<=i_2;++i)
		{
			elke=(i+1-h_eke01[medium-1].x)/h_eke01[medium-1].y;
			eke=exp(elke);
			lelke=i+1;
			p2=eke*(eke+rm*2);
			beta2=p2/(p2+rm*rm);
			chi_a2=h_xcc[medium-1]/(p2*4*h_blcc[medium-1]);
			//ededx=h_ededx[lelke+(medium-1)*MXEKE-1].y*elke+h_ededx[lelke+(medium-1)*MXEKE-1].x;
			ededx=h_dedx[lelke+(medium-1)*MXEKE-1].y*elke+h_dedx[lelke+(medium-1)*MXEKE-1].x;
			estepx=p2*2*beta2*ededx/eke/h_xcc[medium-1]/(log(1.0f/chi_a2+1)*(chi_a2+1)-1);
			estepx*=h_ximax;
			if(estepx>h_estepe)
			{
				estepx=h_estepe;
			}
			ekef=(1-estepx)*eke;
			if(ekef<=h_e_array[(medium-1)*MXEKE])
			{
				sip1=(h_e_array[(medium-1)*MXEKE]-ekef)/dedx0;
				ekef=h_e_array[(medium-1)*MXEKE];
				elkef=(1-h_eke01[medium-1].x)/h_eke01[medium-1].y;
				lelkef=1;
			}
			else
			{
				elkef=log(ekef);
				lelkef=h_eke01[medium-1].y*elkef+h_eke01[medium-1].x;
				leip1l=lelkef+1;
				eip1l=(leip1l-h_eke01[medium-1].x)/h_eke01[medium-1].y;
				eip1=h_e_array[(medium-1)*MXEKE+leip1l-1];
				aux=(eip1-ekef)/eip1;
				elktmp = 0.5*(elkef+eip1l+0.25*aux*aux*(1+aux*(1+0.875*aux)));
				ektmp=(ekef+eip1)*0.5f;
				lelktmp=lelkef;
				//ededx=h_ededx[lelktmp+(medium-1)*MXEKE-1].y*elktmp+h_ededx[lelktmp+(medium-1)*MXEKE-1].x;
				ededx=h_dedx[lelktmp+(medium-1)*MXEKE-1].y*elktmp+h_dedx[lelktmp+(medium-1)*MXEKE-1].x;
				//aux=h_ededx[lelktmp+(medium-1)*MXEKE-1].y/ededx;
				aux=h_dedx[lelktmp+(medium-1)*MXEKE-1].y/ededx;
				r_1=(eip1-ekef)/ektmp;
				sip1=(eip1-ekef)/ededx*(1+aux*(aux*2+1)*(r_1*r_1)/24);
			}
			sip1=sip1+h_range_ep[(medium-1)*2*MXEKE+(lelke-1)*2]-h_range_ep[(medium-1)*2*MXEKE+lelkef*2];//????
			h_tmxs[i+(medium-1)*MXEKE-1].y=(sip1-si)*h_eke01[medium-1].y;
			h_tmxs[i+(medium-1)*MXEKE-1].x=sip1-h_tmxs[i-1+(medium-1)*MXEKE].y*elke;
			si=sip1;
		}
		h_tmxs[neke+(medium-1)*MXEKE-1].x=h_tmxs[neke-1+(medium-1)*MXEKE-1].x;
		h_tmxs[neke+(medium-1)*MXEKE-1].y=h_tmxs[neke-1+(medium-1)*MXEKE-1].y;
	}
	return 1;
}

int edgset(int nreglo,int nreghi)
{   
	//float a[MXTRANS][MXELEMENT];
	//bool got_data=false;
	int i_1;
	bool do_relax;
	int i,i1,j,j1,k;
	//int k1;
	int zn;
	//float binding_energies1[MXSHELL][],binding_energies2[][],interaction_prob1[][],interaction_prob2[][];
	//logical is_batch
	//integer NREGLO,NREGHI
	//integer*4 i,j,k,jj,iz
	//logical do_relax
	//logical got_data
	//save got_data
	//data got_data/.false./
	//if (got_data)
	//{
	//	return ;
	//}
	//logstr(" Output from subroutine EDGSET:");
	//*======================='

	do_relax=true;
	string file1 = ("\\photo_relax.data");
	string relax_file1 = string(data_dir) + file1;
	FILE *photo_relax = fopen(relax_file1.c_str(),"r");
	if(photo_relax==NULL)
	{
		logstr("photo_relax data file not found: %s\n",relax_file1);
		//exit(1);
		return 1;
	}

	//logstr(" Output from subroutine EDGSET:\n");
/*
	for (k=0;k<MXREG;k++)
	{
		if((h_iedgfl[k]>0)&&(h_iedgfl[k]<=MXELEMENT))
		{
			do_relax=true;
			goto L3102;
		}
	}
	*/
//L3102:

	if(!do_relax)
	{
		//logstr( " Atomic relaxations not requested!\n " );
		return 1 ;
	}
	//logstr( " Atomic relaxations requested! \n" );
	/*      write(i_log,'(a/)') ' Atomic relaxations requested! ' */
    //logstr("Reading photo-absorption data .....");
	//got_data=true;
	for(i=0;i<MXELEMENT;i++)
	{ 
		fscanf(photo_relax,"%d",&zn);
		for (k=0;k<MXSHELL;k++)
		{
			fscanf(photo_relax,"%f",&h_binding_energies[MXSHELL*i+k]); //明天再讨论数据的存储位置
		}

	}
    for(i1=0;i1<MXELEMENT;i1++)
	{
		for (j1=0;j1<MXSHELL;j1++)
		{
			h_binding_energies[MXSHELL*i1+j1]=h_binding_energies[MXSHELL*i1+j1]*1e-6;
		}
	}

//	read(i_photo_relax,*)             //空读一行
	for(i=0;i<MXELEMENT;i++)
	{   
		fscanf(photo_relax,"%d",&zn);
		for (j=0;j<(MXSHELL-1);j++)
		{
			fscanf(photo_relax,"%f",&h_interaction_prob[MXSHELL*i+j]);
		}
		h_interaction_prob[MXSHELL*i+MXSHELL-1]=1.01f;
	}
//	/*      write(i_log,'(a)') ' Done' */
// /*      write(i_log,'(/a$)') ' Reading relaxation data .....' */
// /*      read(i_photo_relax,*) */
//	read(i_photo_relax,*)                     /空读一行

	for(i=0;i<MXELEMENT;i++)
	{
		fscanf(photo_relax,"%d",&zn);
		for(k=0;k<19;k++)
		{fscanf(photo_relax,"%f",&h_relaxation_prob[MXTRANS*i+k]);}
	}

	for (i=0;i<MXELEMENT;i++)
	{
		fscanf(photo_relax,"%d",&zn);
		for(k=19;k<26;k++)
		{
			fscanf(photo_relax,"%f",&h_relaxation_prob[MXTRANS*i+k]);
		}
	}
	for(i=0;i<MXELEMENT;i++)
	{
		fscanf(photo_relax,"%d",&zn);
		for(k=26;k<32;k++)
		{fscanf(photo_relax,"%f",&h_relaxation_prob[MXTRANS*i+k]);}
	}
	for (i=0;i<MXELEMENT;i++)
	{
		fscanf(photo_relax,"%d",&zn);
		for (k=32;k<37;k++)
		{
			fscanf(photo_relax,"%f",&h_relaxation_prob[MXTRANS*i+k]);

		}
	}
	for (i=0;i<MXELEMENT;i++)
	{
		fscanf(photo_relax,"%d",&zn);
		fscanf(photo_relax,"%f",&h_relaxation_prob[MXTRANS*i+37]);
	}
	//logstr(" Done\n");
	fclose(photo_relax);
	//logstr("Reading photo cross section data .....");
	string file2 = ("\\photo_cs.data");
	string cs_file2 = string(data_dir) + file2;
	FILE *photo_cs = fopen(cs_file2.c_str(),"r");
	if(photo_cs==NULL)
	{
		logstr("photo_cs data file not found: %s\n",cs_file2);
		//exit(1);
		return 1;
	}
	
	for (i=0;i<MXELEMENT;i++)
	{
		fscanf(photo_cs,"%d %d",&zn,&h_edge_number[i]);
		i_1=h_edge_number[i];
		for(j=0;j<i_1;j++)
		{
			fscanf(photo_cs,"%f %f %f %f %f",&h_edge_a[MXEDGE*i+j],&h_edge_b[MXEDGE*i+j],&h_edge_c[MXEDGE*i+j],&h_edge_d[MXEDGE*i+j],&h_edge_energies[MXEDGE*i+j]); //note that the edges are arranged from high edge energy to low.
		}
	}
	fclose(photo_cs);
	//logstr(" Done\n");

	return 1;
}

void init_compton()
{
	int i_1,i_2,i_3;
	//float d_1;
	int i,j,iz,j_h,j_l,nsh;
	float aux;
	bool getd;
	float atav;
	float pztot;
	int medium;
	double aux_erf;
	float rm=0.5110034;
	getd=false;
	char dummy[80];

	//logstr("Output from init_compton :\n");
//we will no longer check the bound compton scatter requirment region by region.
// just set it to be always true by default
/*
for (j=0;j<MXREG;j++)
	{
		medium=h_med(j);
		if((medium>0)&&(medium<=h_nmed))
		{
			if(h_ibcmp[j]>0)
			{
				getd=true;
				goto L1782;
			}
		}
	}
L1782:
*/
	getd= true;
	if(!getd)
	{
		//logstr(" Bound Compton scattering not requested! \n");
		return;
	}
	//logstr(" Bound Compton scattering  requested! \n");
	
	string file1 = ("\\incoh.data");
	string incoh_file1 = string(data_dir) + file1;
	FILE *incoh = fopen(incoh_file1.c_str(),"r");
	if(incoh==NULL)
	{
		logstr("incoh data file not found: %s\n",incoh_file1);
		//exit(1);
		return;
	}
	for(j=0;j<18;j++)
	{
		fgets (dummy,80,incoh);  //skip first line?
	}
	for(j=0;j<MXTOTSH;j++)
	{
		fscanf(incoh,"%d %d %d %f %f",&h_iz_array[j],&h_shn_array[j],&h_ne_array[j],&h_jo_array[j],&h_be_array[j]);
		h_jo_array[j]=h_jo_array[j]*137.0;
		h_be_array[j]=h_be_array[j]*(1e-6)/rm;
		aux_erf=0.70710678119*(1.0f+0.3f*h_jo_array[j]);
		//h_erfjo_array[j]=0.82436063535*(erf1(aux_erf)-1.0f)  ;      // !!!!!!!!!!!!!!!!!!!!!!!!!!???????? 调用一个外部函数 init spline erf1
		h_erfjo_array[j]=0.82436063535*(erf1(aux_erf))  ;      // !!!!!!!!!!!!!!!!!!!!!!!!!!???????? 调用一个外部函数 init spline erf1
	}

    //logstr("Initializing Bound Compton scattering ......total number of medium: h_nmed=%d\n ",h_nmed);
    
	for(medium=0;medium<h_nmed ;medium++)
	{
		pztot=0;
		nsh=0;
		i_2=h_nne[medium];
//		logstr("medium=%d\n",medium);
		for(i=0;i<i_2;i++)
		{
			iz=int(h_zelem[i*MXMED+medium]);
//			logstr("\n looking for shells for element Z=%d\n",iz);
			for(j=0;j<MXTOTSH;j++)
			{

				if(iz==h_iz_array[j])
				{
					nsh=nsh+1;
					if(nsh>MXMDSH)
					{
						logstr("***************** Error: for medium # %d, the number of shells is larger than MXMDSH",medium);
					    exit(1);
					}
					h_shell_array[medium*MXMDSH+nsh-1]=j+1;
					aux=h_pz[i*MXMED+medium]*h_ne_array[j];
					h_eno_array[medium*MXMDSH+nsh-1]=aux;
                    pztot=pztot+aux;
				}
			}
		}
		if(nsh==0)
		{
			logstr("***************** Error: for medium # %d has zero shells!",medium);
		    exit(1);
		}
		h_n_shell[medium]=nsh;
		//logstr("medium %d has %d shells\n", medium, nsh);
		i_3=nsh;
		for (i=0;i<i_3;i++)
		{
			j=h_shell_array[medium*MXMDSH+i];
			h_eno_array[medium*MXMDSH+i]=h_eno_array[medium*MXMDSH+i]/pztot;
		//	write(i_log,'(i3,i4,i3,f9.5,e10.3,f10.3)') i,j,h_shn_array(j),en
     // *    o_array(i,medium), h_jo_array(j),h_be_array(j)*rm*1000.
			h_eno_array[MXMDSH*medium+i]=-1*h_eno_array[MXMDSH*medium+i];
			h_eno_atbin_array[medium*MXMDSH+i]=i+1;
		}
		atav=1.0f/nsh;
		i_2=nsh-1;
		for(i=0;i<i_2;i++)
		{
			for(j_h=0;j_h<i_2;j_h++)
			{
				if(h_eno_array[medium*MXMDSH+j_h]<0)
				{
					if(abs(h_eno_array[medium*MXMDSH+j_h])>atav)
					{
						goto L1862;
					}                                           // ?
				}
			}
L1862:         
			i_3=nsh-1;
			for(j_l=0;j_l<i_3;j_l++)
			{
				if(h_eno_array[MXMDSH*medium+j_l]<0)
				{
					if(abs(h_eno_array[MXMDSH*medium+j_l])<atav)
					{
						goto L8722;
					}
				}
			}
L8722:
			aux=atav-abs(h_eno_array[MXMDSH*medium+j_l]);
			h_eno_array[medium*MXMDSH+j_h]=h_eno_array[medium*MXMDSH+j_h]+aux;
			h_eno_array[medium*MXMDSH+j_l]=-h_eno_array[medium*MXMDSH+j_l]/atav+j_l+1;
			h_eno_atbin_array[medium*MXMDSH+j_l]=j_h+1;
			if(i==(nsh-2))
			{
				h_eno_array[MXMDSH*medium+j_h]=j_h+1+1;
			}
		}
		i_1=nsh;
		for(i=0;i<i_1;i++)
		{
			if(h_eno_array[medium*MXMDSH+i]<0)
		{ 
			h_eno_array[medium*MXMDSH+i]=1+i+1;
			
			}		
		}

}
	//we shoul always call edgset already. so no need to check and call as the following code. 
/*
getd=false;
	for (j=0;j<MXREG;j++)
	{
		if((h_iedgfl[j]>0)&&iedgl[j]<=MXELEMENT)
		{
			getd=true;
			goto L1892;
		}
	}
L1892:

	if(getd)
	{
		return ;
	}
	// write(i_log,'(a/,a/,a/,a//)') ' In subroutine init_compton: ', '
    //  * fluorescence not set but relaxation data are required for ', '
    // *bound Compton scattering. ', '   calling EDGSET. '
	h_iedgfl[0]=1;
	edgeset(1,1);                        //调用edgeset函数
	h_iedgfl[0]=1;
	*/
	return ;
}

float fcoulc(float z)
{
	float fine,asq;
	float ret_val;
	fine=137.03604;
	asq=z/fine;
	asq=asq*asq;
	ret_val=asq * (1.f/(asq + 1.f) + 0.20206f+asq*(asq*(asq* (-0.002f) 
	    + .0083f) - 0.0369f));
	return ret_val;

}
float xsif(float z)
{
	int iz;
	float a1440,a183,ret_val;
	float fcoulc(float);
	float alrad[4]={5.31,4.79,4.74,4.71};
	float alradp[4]={6.144,5.621,5.805,5.924};
	a1440=1194.0;
	a183=184.15;
	if (z<=4)
	{
		iz=z;
	ret_val=alradp[iz-1]/(alrad[iz-1]-fcoulc(z));
	}
	else ret_val=log(a1440*pow(z,-0.666667f))/(log(a183*pow(z,-0.33333f))-fcoulc(z));

	return ret_val;

}
int fix_brems()
{
	float fcoulc(float z);
	float xsif(float z);
	float fcoulc(float );
	float xsif(float);
	int i_1,i_2;

	int i;
	float fc,pi,zb,zf,zg,xi,zi,zt,zv,aux;
	float fmax1,fmax2;
	int medium;
	i_1=h_nmed;
	for(medium=0;medium<i_1;medium++)
	{
		
		h_log_ap[medium]=log(h_ap[medium]);
        zt=0;
		zb=0;
		zf=0;
		i_2=h_nne[medium];
		for(i=0;i<i_2;i++)
		{
			zi=h_zelem[MXMED*i+medium];
			pi=h_pz[MXMED*i+medium];
			fc=fcoulc(zi);
			xi=xsif(zi);
			aux=pi*zi*(zi+xi);
			zt=zt+aux;
			zb=zb-aux*log(zi)/3;
			zf=zf+aux*fc;
		}
		zv=(zb-zf)/zt;
		zg=zb/zt;
		fmax1=2*(20.863+4*zg)-2*(20.029+4*zg)/3;
		fmax2=2*(20.863+4*zv)-2*(20.029+4*zv)/3;
		h_dl1[medium*8+0]=(20.863+4*zg)/fmax1;
		h_dl2[medium*8+0]=-3.242/fmax1;
		h_dl3[medium*8+0]=0.625/fmax1;
		h_dl4[medium*8+0]=(21.12+4*zg)/fmax1;
		h_dl5[medium*8+0]=-4.184/fmax1;
		h_dl6[medium*8+0]=0.952;
		h_dl1[medium*8+1]=(20.029+4*zg)/fmax1;
		h_dl2[medium*8+1]=-1.93/fmax1;
		h_dl3[medium*8+1]=-0.086/fmax1;
		h_dl4[medium*8+1]=(21.12+4*zg)/fmax1;
        h_dl5[medium*8+1]=-4.184/fmax1;
		h_dl6[medium*8+1]=0.952;
		h_dl1[medium*8+2]=(20.863 + 4*zv)/fmax2;
		h_dl2[medium*8+2]=-3.242/fmax2;
		h_dl3[medium*8+2]=0.625/fmax2;
		h_dl4[medium*8+2]=(21.12+4*zv)/fmax2;
		h_dl5[medium*8+2]=-4.184/fmax2;
		h_dl6[medium*8+2]=0.952;
		h_dl1[medium*8+3]=(20.029+4*zv)/fmax2;
		h_dl2[medium*8+3]=-1.93/fmax2;
		h_dl3[medium*8+3]=-0.086/fmax2;
		h_dl4[medium*8+3]=(21.12+4*zv)/fmax2;
		h_dl5[medium*8+3]=-4.184/fmax2;
		h_dl6[medium*8+3]=0.952;
		h_dl1[medium*8+4]=(3*(20.863 + 4*zg) - (20.029 + 4*zg));
		h_dl2[medium*8+4]=(3*(-3.242) - (-1.930));
		h_dl3[medium*8+4]=(3*(0.625)-(-0.086));
		h_dl4[medium*8+4]=(2*21.12+8*zg);
		h_dl5[medium*8+4]=(2*(-4.184));
		h_dl6[medium*8+4]=0.952;
		h_dl1[medium*8+5]=(3*(20.863 + 4*zg) + (20.029 + 4*zg));
		h_dl2[medium*8+5]=(3*(-3.242) + (-1.930));
		h_dl3[medium*8+5]=(3*0.625+(-0.086));
		h_dl4[medium*8+5]=(4*21.12+16*zg);
		h_dl5[medium*8+5]=(4*(-4.184));
		h_dl6[medium*8+5]=0.952;
		h_dl1[medium*8+6]=(3*(20.863 + 4*zv) - (20.029 + 4*zv));
		h_dl2[medium*8+6]=(3*(-3.242) - (-1.930));
		h_dl3[medium*8+6]=(3*(0.625)-(-0.086));
		h_dl4[medium*8+6]= (2*21.12+8*zv);
		h_dl5[medium*8+6]=(2*(-4.184));
		h_dl6[medium*8+6]=0.952;
		h_dl1[medium*8+7]=(3*(20.863 + 4*zv) + (20.029 + 4*zv));
		h_dl2[medium*8+7]=(3*(-3.242) + (-1.930));
		h_dl3[medium*8+7]=(3*0.625+(-0.086));
		h_dl4[medium*8+7]=(4*21.12+16*zv);
		h_dl5[medium*8+7]=(4*(-4.184));
		h_dl6[medium*8+7]=0.952;
		h_bpar[medium*2+1]=h_dl1[medium*8+6]/(3*h_dl1[medium*8+7]+h_dl1[medium*8+6]);
		h_bpar[medium*2+0]=12*h_dl1[medium*8+7]/(3*h_dl1[medium*8+7]+h_dl1[medium*8+6]);
			
	}
	return 0;
}

int eii_init()
{
	int occn_numbers[4]={2,2,2,4};
//	address a_1[6];                         ??????????????????????
	int i_1,i_2,i_5;
	//int i_3[6],i_4;
	double r_1,r_2;
	//double r_3,r_4;
	//char ch_1[512];
	//char eii_file[128];
	double  dedx_old,wbrem;
	//double  sigm_old;
	//int eii_unit;
	double sum_occn,sum_dedx,e;
	int i,j,k;
	double u, e_eii_min, sigma_old, sigma_max, wbrem_old, p2, 
	    sum_sigma;
	float aux_array[N_EII_BINS];
	int tmp_array[MXELEMENT];
	double ec,de;
	int ii,jj,iz;
    double ecc,sig;
	int iii,jjj,i8;
	int ish,nsh;
	double tau,aux,uwm,sh_0,sh_1,ss_0,ss_1;
	bool is_monotone;
	double av_e;
	int imed,iele,nbin;
	float emax, fmax1;
	double loge,dedx,sigm;
	//bool getd;
	double sigo;
	//int itmp;
	double wmax,beta2, sig_j,sigma;
	int nskip;
	//int want_eii_unit;
    double sum_z, sum_a, sum_wa,sum_sh,sum_pz;

	float con_med;
	//int eii_out;
	int nsh_tot;
	for (j=0;j<MXELEMENT;j++)
	{
		h_eii_nshells[j]=0.0f;
	}
	for (j=0;j<MXMED;j++)
	{
		h_eii_nsh[j]=0;
	}

	if(h_eii_flag==0)
	{
		return 1;
	}

//edgset should always called before this. so no need to check.

	e_eii_min=1e30f;
	i_1=h_nmed;
	for(imed=0;imed<i_1;imed++)
	{
		if((h_ae[imed]-h_rm)<e_eii_min)
		{
           e_eii_min=h_ae[imed]-h_rm;
		}
		if(h_ap[imed]<e_eii_min)
		{
			e_eii_min=h_ap[imed];
		}
	}

//	write(i_log,*) ' '
    //logstr("eii_init: minimum threshold energy found: %f\n",e_eii_min);


	/*determine elements that need to load EII data*/
	i_1=h_nmed;
	for(imed=0;imed<i_1;imed++)
	{
		i_2=h_nne[imed];
		for(iele=0;iele<i_2;iele++)
		{
			iz=int(h_zelem[MXMED*iele+imed]+0.5);
			if(h_eii_nshells[iz-1]==0)
			{
				nsh=0;
				for(ish=0;ish<4;ish++)
				{
					if(h_binding_energies[MXSHELL*(iz-1)+ish]>e_eii_min)
					{
						nsh=nsh+1;
					}
				}
				h_eii_nshells[iz-1]=nsh;
			}
		}

	}
	/* total number of shells that need to be loaded */

	nsh=0;
	for(iz=0;iz<MXELEMENT;iz++)
	{
		nsh+=h_eii_nshells[iz];
	}
	if(nsh==0)
	{
		logstr("*** EII requested but no shells with binding energies above the specified threshold found, turning off EII\n");
		h_eii_flag=0;
	}
	if(nsh>MAX_EII_SHELLS)
	{
		logstr("*** Number of shells with binding energies greater than the specified thresholds is %d\n",nsh);
        logstr(" Increase the macro $MAX_EII_SHELLS and retry");
	    exit(1);
	}
    //logstr("eii_init: number of shells to simulate EII: %d\n ",nsh);


	nsh_tot=nsh;
	tmp_array[0]=0;
	for(j=1;j<MXELEMENT;j++)
	{
		tmp_array[j]=tmp_array[j-1]+h_eii_nshells[j-1];
	}
	/*Get relaxation data if necessary so that the binding energies
are available*/

//	itmp=h_iedgfl[0];
//	h_iedgfl[0]=1;
//	edgset(1,1);                                      // edgset should be called already.
//	h_iedgfl[0]=itmp;

/* set EII active shells per medium and for each element */
	for(imed=0;imed<h_nmed;imed++)
	{
		nsh=0;
		i_2=h_nne[imed];
		for(iele=0;iele<i_2;iele++)
		{
			iz=int(h_zelem[MXMED*iele+imed]+0.5f);
			h_eii_no[iele*MXMED+imed]=h_eii_nshells[iz-1];
			nsh+=h_eii_nshells[iz-1];
			if(h_eii_nshells[iz-1]>0)
			{
				h_eii_first[MXMED*iele+imed]=tmp_array[iz-1]+1;
			}
			else
			{
				h_eii_first[MXMED*iele+imed]=0.0f;
			}
		}
		h_eii_nsh[imed]=nsh;
	}

                     /* read EII data */

	string file1 = ("eii_ik.data");
	string eii_file1 = string(data_dir) + file1;
	FILE *eii_data = fopen(eii_file1.c_str(),"r");
	if(eii_data==NULL)
	{
		logstr("Eii data file not found: %s\n",eii_file1);
		//exit(1);
		return 1;
	}
	
	fscanf (eii_data,"%d\n",&nskip);
    char dummy[81];
	i_1=nskip;
	for (j=0;j<i_1;j++)
	{
		fgets(dummy,80,eii_data);
	}


	fscanf (eii_data,"%f %d",&emax,&nbin);

//      IF (( nbin .NE. 9039 )) THEN
	if (nbin!=N_EII_BINS)
       {
		   logstr("***************** Error:\n Inconsistent EII data file\n Exiting ...");
	       exit(1);
	}
        ii = 0;
		for (j=0;j<MXELEMENT;j++)
		{
			fscanf (eii_data,"%d %d",&iz,&nsh);
			if(nsh<h_eii_nshells[iz-1])
			{
				logstr("EII data file has data for %d shells for element %d\n",nsh,iz);
				logstr("but according to binding energies and thresholds %d shells are required\n",h_eii_nshells[iz-1]);
				logstr("fatal error!, exiting ...");
				 exit(1);
			}
			for (ish=1;ish<=nsh;ish++)
			{
				fscanf (eii_data,"%f",&fmax1);
			
				for (i8=0;i8<N_EII_BINS;i8++)
				{
					fscanf(eii_data,"%f",&aux_array[i8]);
				}
				if (ish<=h_eii_nshells[iz-1])
				{

				ii=ii+1;
	            h_eii_z[ii-1]=iz;
	            h_eii_sh[ii-1]=ish;
	            h_eii[ii-1].x=nbin;
	            h_eii[ii-1].x=h_eii[ii-1].x/log(emax/h_binding_energies[MXSHELL*(iz-1)+ish-1]);
	            h_eii[ii-1].y=1-h_eii[ii-1].x*log(h_binding_energies[MXSHELL*(iz-1)+ish-1]);
            	i_2=nbin;
		 	
	            for (k=1;k<=i_2;k++)
	   {
		if(k>1)
		{
			sigo=fmax1*aux_array[k-2];
		}
		else
		{
			sigo=0;
		}
		loge=(k-h_eii[ii-1].y)/h_eii[ii-1].x;
		iii=nbin*(ii-1)+k;
		h_eii_xsection[iii-1].x=(fmax1*aux_array[k-1]-sigo)*h_eii[ii-1].x;
		h_eii_xsection[iii-1].y=sigo-h_eii_xsection[iii-1].x*loge;
				}
				
				}
			}

		     
			
			if(ii==nsh_tot)
	       {
	                	goto L3882;
	             }
			}
		
L3882:

	fclose (eii_data);  


	for(imed=0;imed<h_nmed;imed++)
	{
		ec=h_ae[imed]-h_rm;
		ecc=min(ec,h_ap[imed]);
		sum_z=0.f;
		sum_pz=0.f;
		sum_a=0.f;
		sum_wa=0.f;
		for(iele=0;iele<h_nne[imed];iele++)
		{
			sum_z+=h_pz[MXMED*iele+imed]*h_zelem[MXMED*iele+imed];
			sum_pz+=h_pz[MXMED*iele+imed];
			sum_wa+=h_rhoz[MXMED*iele+imed];
			sum_a+=h_pz[MXMED*iele+imed]*h_wa[MXMED*iele+imed];
		}
		con_med=h_rho[imed]/1.6605655/sum_a;
		h_eii_cons[imed]=con_med;
		if (h_eii_nsh[imed]>0)
		{
			is_monotone=true;
			sigma_max=0;
//			i_3=h_meke[imed];
			int medium;
			for(j=0;j<h_meke[imed];j++)
			{
				loge=(j+1-h_eke01[imed].x)/h_eke01[imed].y;
				e=exp(loge);
				tau=e/h_rm;
				beta2=tau*(tau+2)/((tau+1)*(tau+1));
				p2=2*h_rm*tau*(tau+2);
				int lloge=j;
				medium=imed;
				//dedx=h_ededx[medium*MXEKE+lloge].y*loge+h_ededx[medium*MXEKE+lloge].x;
				dedx=h_dedx[medium*MXEKE+lloge].y*loge+h_dedx[medium*MXEKE+lloge].x;
				if(e>h_ap[medium]||2>2*ec)
				{
					//sig=h_esig[medium*MXEKE+lloge].y*loge+h_esig[medium*MXEKE+lloge].x;
					sig=h_sig[medium*MXEKE+lloge].y*loge+h_sig[medium*MXEKE+lloge].x;
				}
				else 
				{
					sig=0;
				}
				if(e>2*ec)
				{
					wbrem=h_ebr1[medium*MXEKE+lloge].y*loge+h_ebr1[medium*MXEKE+lloge].x;
				    sigm=sig*(1-wbrem);
//					logstr("lloge,wbrem,sig,sigm %d,%f,%f,%f\n",lloge,wbrem,sig,sigm);
				}
				else 
				{
					sigm=0;
					wbrem=1;
				     }
				sum_occn=0.f;
				sum_sigma=0.f;
				sum_dedx=0.f;
//				i_4=h_nne[imed];
				for(iele=0;iele<h_nne[imed];iele++)
				{
					iz=int(h_zelem[MXMED*iele+imed]+0.5f);
					sum_sh=0.0f;
					i_5=h_eii_no[MXMED*iele+imed];
					for(ish=0;ish<i_5;ish++)
					{
						jj=h_eii_first[MXMED*iele+imed]+ish-1 +1;
						jjj=h_eii_sh[jj-1];
						u=h_binding_energies[MXSHELL*(iz-1)+jjj-1];
						wmax=(e+u)/2;
						uwm=u/wmax;
						if((u<e)&&(u>ecc))
						{
//						logstr("iz,jj,jjj,u,%d,%d,%d,%f\n",iz,jj,jjj,u);
						sum_sh +=occn_numbers[jjj-1];
						ss_0=2*(log(p2/u)-uwm*uwm*uwm*log(p2/wmax)-(beta2+0.833333)*(1-uwm*uwm*uwm))/3/u;
						r_1=e+h_rm;
						r_2=tau+1;
						sh_0=((1-uwm)*(1+uwm/(2-uwm))+u*(wmax-u)/(r_1*r_1)-(2*tau+1)/(r_2*r_2)*uwm/2*log((2-uwm)/uwm))/u;
						ss_1=log(p2/u)-uwm*uwm*log(p2/wmax)-(beta2+1)*(1-uwm*uwm);
						sh_1=log(wmax/u/(2-uwm))+2*(wmax-u)/(2*wmax-u)+(wmax*wmax-u*u)/((e+h_rm)*(e+h_rm))/2-(2*tau+1)/((tau+1)*(tau+1))*log((2*wmax-u)/wmax);
						av_e=(ss_1+sh_1)/(ss_0+sh_0);
						i=h_eii[jjj-1].x*loge+h_eii[jjj-1].y;
						i=(jj-1)*N_EII_BINS+i;
						sig_j=h_eii_xsection[i-1].x*loge+h_eii_xsection[i-1].y;
						sig_j=sig_j*h_pz[MXMED*iele+imed]*con_med;
						sum_sigma+=sig_j;
						sum_dedx+=sig_j*av_e;
						}
					}
					sum_occn+=sum_sh*h_pz[MXMED*iele+imed];
				}
				sigm=sigm+sum_sigma;
				dedx=dedx-sum_dedx;
				aux=ec/e;
				if(e>2*ec)
				{
					r_1=tau/(tau+1);
					r_2=tau+1;
				sigo=sum_occn*0.153536f*h_rho[imed]/(beta2*ec)*((1-aux*2)*(1+aux/(1-aux)+r_1*r_1*aux/2)-(tau*2+1)/(
			    r_2* r_2)*aux*log((1-aux)/aux))/sum_a;
				de = sum_occn*0.153536f*h_rho[imed]/beta2*(log(0.25f/aux/(1-aux))+(1-aux* 2)/(1-aux)+r_1*r_1*(1-aux*4*aux)/8-(tau*2+1)/(r_2*r_2)*log((1- aux)*2))/sum_a;
				sigm-=sigo;
				dedx+=de;
//				logstr("sigm,dedx,aux,sum_occn,sigo,de  %f,%f,%f,%f,%f,%f\n",sigm,dedx,aux,sum_occn,sigo,de);
				}
				sigma=sigm+wbrem*sig;
//				logstr(" %f,%f,%f,%f,%f,%f\n",sigm,dedx,aux,sigma,wbrem,sig);
				if((sigma/dedx)>sigma_max)
				{
					sigma_max=sigma/dedx;
				}
				if(sigma>0.0f)
				{
					wbrem=wbrem*sig/sigma;
				}
				else 
				{
					wbrem=1.0f;
				}
				if(j>0)
				{
					//h_ededx[imed*MXEKE+j-1].y=(dedx-dedx_old)*h_eke01[imed].y;
					h_dedx[imed*MXEKE+j-1].y=(dedx-dedx_old)*h_eke01[imed].y;
					//h_ededx[imed*MXEKE+j-1].x=dedx-h_ededx[imed*MXEKE+j-1].y*loge;
					h_dedx[imed*MXEKE+j-1].x=dedx-h_dedx[imed*MXEKE+j-1].y*loge;
					//h_esig[imed*MXEKE+j-1].y=(sigma-sigma_old)*h_eke01[imed].y;
					h_sig[imed*MXEKE+j-1].y=(sigma-sigma_old)*h_eke01[imed].y;
					//h_esig[imed*MXEKE+j-1].x=sigma-h_esig[imed*MXEKE+j-1].y*loge;
					h_sig[imed*MXEKE+j-1].x=sigma-h_sig[imed*MXEKE+j-1].y*loge;
					h_ebr1[imed*MXEKE+j-1].y=(wbrem-wbrem_old)*h_eke01[imed].y;
					h_ebr1[imed*MXEKE+j-1].x=wbrem-h_ebr1[imed*MXEKE+j-1].y*loge;
					if(sigma/dedx<sigma_old/dedx_old)
					{
						is_monotone=false;
					}
				}
				dedx_old=dedx;
				//sigm_old=sigm;
				sigma_old=sigma;
				wbrem_old=wbrem;
			}
			//h_ededx[MXEKE*imed+h_meke[imed]-1].y=h_ededx[MXEKE*imed+h_meke[imed]-2].y;
			h_dedx[MXEKE*imed+h_meke[imed]-1].y=h_dedx[MXEKE*imed+h_meke[imed]-2].y;
			//h_ededx[MXEKE*imed+h_meke[imed]-1].x=h_ededx[MXEKE*imed+h_meke[imed]-2].x;
			h_dedx[MXEKE*imed+h_meke[imed]-1].x=h_dedx[MXEKE*imed+h_meke[imed]-2].x;
			//h_esig[MXEKE*imed+h_meke[imed]-1].y=h_esig[MXEKE*imed+h_meke[imed]-2].y;
			h_sig[MXEKE*imed+h_meke[imed]-1].y=h_sig[MXEKE*imed+h_meke[imed]-2].y;
			//h_esig[MXEKE*imed+h_meke[imed]-1].x=h_esig[MXEKE*imed+h_meke[imed]-2].x;
			h_sig[MXEKE*imed+h_meke[imed]-1].x=h_sig[MXEKE*imed+h_meke[imed]-2].x;
			h_ebr1[MXEKE*imed+h_meke[imed]-1].y=h_ebr1[MXEKE*imed+h_meke[imed]-2].y;
			h_ebr1[MXEKE*imed+h_meke[imed]-1].x=h_ebr1[MXEKE*imed+h_meke[imed]-2].x;
			h_sig_ismonotone[2*imed] = is_monotone ; 
			//h_esig_e[imed]=sigma_max;
			h_sig_e[imed]=sigma_max;
		}

	} 
	return 0;
}
