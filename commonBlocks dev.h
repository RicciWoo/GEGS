//following are the variables defined in EGSnrc common Blocks.
//The array pairs for linear interploations such as ge0[] and ge1[], are now combined into a __device__ float2 array ge[]

//#include "EGSnrc_Parameters.h"
/* Common Block Declarations */

//COMMON/BOUNDS/
    //__device__ float ecut[MXREG], pcut[MXREG], vacdst;
    
//COMMON/BREMPR/
    __device__ float dl1[8*MXMED], dl2[8*MXMED], 
	    dl3[8*MXMED], dl4[8*MXMED], 
	    dl5[8*MXMED], dl6[8*MXMED];
	//__device__ float alphi[2*MXMED];
	__device__ float bpar[2*MXMED];
	//__device__ float delpos[2*MXMED];
	__device__ float wa[MXMED*MXEL],
	    pz[MXMED*MXEL], zelem[MXMED*MXEL];
	//__device__ float rhoz[MXMED*MXEL];
	__device__ float delcm[MXMED], zbrang[MXMED], lzbrang[MXMED];
	//__device__ float pwr2i[MXPWR2I]; 
    __device__ int nne[MXMED], ibrdst, iprdst;
	//__device__ int ibr_nist, pair_nrc, itriplet;
    //__device__ char asym[MXMED*MXEL*2]	;
//#define brempr_1 brempr_

/*
//COMMON/nist_brems
    __device__ float nb_fdata[MXBRXS_PLUS_1*MXBRES*MXMED], nb_xdata[MXBRXS_PLUS_1*MXBRES*MXMED], nb_wdata[MXBRXS*MXBRES*MXMED];
    __device__ int nb_idata[MXBRXS*MXBRES*MXMED];
    __device__ float nb_emin[MXMED], nb_emax[MXMED], nb_lemin[MXMED], nb_lemax[MXMED], 
	    nb_dle[MXMED], nb_dlei[MXMED];
*/
	__device__ float log_ap[MXMED];
//#define nist_brems1 nist_brems

/*
//COMMON/nrc_pair
    __device__ float nrcp_fdata[NRC_PAIR_NXX*NRC_PAIR_NEE*MXMED], nrcp_wdata[NRC_PAIR_NXX*NRC_PAIR_NEE*MXMED]	;
    __device__ int nrcp_idata[NRC_PAIR_NXX*NRC_PAIR_NEE*MXMED];
    __device__ float nrcp_xdata[NRC_PAIR_NXX], nrcp_emin, nrcp_emax, nrcp_dle, nrcp_dlei;
//#define nrc_pair1 nrc_pair
*/

/*
//COMMON/triplet_data/
    __device__ float a_triplet[MAX_TRIPLET*MXMED], b_triplet[MAX_TRIPLET*MXMED], 
	    dl_triplet, dli_triplet, bli_triplet, log_4rm;
//#define triplet_data1 triplet_data
*/

//Common/compton_data/
    __device__ int iz_array[MXTOTSH];
    __device__ float be_array[MXTOTSH], jo_array[MXTOTSH], erfjo_array[MXTOTSH];
    __device__ int ne_array[MXTOTSH], shn_array[MXTOTSH], shell_array[MXMDSH*MXMED];
    __device__ float eno_array[MXMDSH*MXMED];
    __device__ int eno_atbin_array[MXMDSH*MXMED], n_shell[MXMED];
	//__device__ int radc_flag;
    //__device__ short ibcmp[MXREG];
//#define compton_data1 compton_data

//COMMON/edge/
	__device__ float binding_energies[MXSHELL*MXELEMENT], 
		interaction_prob[MXSHELL*MXELEMENT], 
		relaxation_prob[MXTRANS*MXELEMENT], 
		edge_energies[MXEDGE*MXELEMENT];
	__device__ int edge_number[MXELEMENT];
	__device__ float edge_a[MXEDGE*MXELEMENT], edge_b[MXEDGE*MXELEMENT],	
		edge_c[MXEDGE*MXELEMENT],
		edge_d[MXEDGE*MXELEMENT];
	//__device__ short iedgfl[MXREG], iphter[MXREG];
//#define edge_1 (edge_._1)
//#define edge_2 (edge_._2)

//COMMON/elecin/
    __device__ float esige_max, psige_max, range_ep[2*MXEKE*MXMED], e_array[MXEKE*MXMED];
	//__device__ float esig_e[MXMED];
	//__device__ float psig_e[MXMED]; 
    //__device__ float2 etae_ms[MXEKE*MXMED];
	//__device__ float2 etap_ms[MXEKE*MXMED];
  	//__device__ float2 q1ce_ms[MXEKE*MXMED];
  	//__device__ float2 q1cp_ms[MXEKE*MXMED];
  	//__device__ float2 q2ce_ms[MXEKE*MXMED];
  	//__device__ float2 q2cp_ms[MXEKE*MXMED];

	__device__ float  sig_e[MXMED*2];
	__device__ float2 eta_ms[MXEKE*MXMED*2];
	__device__ float2 q1c_ms[MXEKE*MXMED*2];
	__device__ float2 q2c_ms[MXEKE*MXMED*2];

  	__device__ float2 blcce[MXEKE*MXMED];
	__device__ float2 eke01[MXMED];

	//__device__ float xr0[MXMED], teff0[MXMED];
	__device__ float blcc[MXMED], xcc[MXMED];
	//__device__ float2 esig[MXEKE*MXMED];
	//__device__ float2 psig[MXEKE*MXMED];
	//__device__ float2 ededx[MXEKE*MXMED];
	//__device__ float2 pdedx[MXEKE*MXMED];
	__device__ float2 ebr1[MXEKE*MXMED];
	__device__ float2 pbr1[MXEKE*MXMED];
	__device__ float2 pbr2[MXEKE*MXMED]; 
	__device__ float2 tmxs[MXEKE*MXMED];

	__device__ float2 sig[MXEKE*MXMED*2];
	__device__ float2 dedx[MXEKE*MXMED*2];

	__device__ float expeke1[MXMED];
	//__device__ int iunrst[MXMED], epstfl[MXMED], iaprim[MXMED];
	__device__ bool sig_ismonotone[2*MXMED]	;
//#define elecin_1 elecin_

//COMMON/eii_data/
    __device__ float2 eii_xsection[MAX_EII_BINS];  //use to be eii_xsection_a, eii_xsection_b.
	__device__ float eii_cons[MXMED];
	__device__ float2 eii[MAX_EII_SHELLS];   //use to be eii_a, eii_b. now combine them to eii.x, eii.y.
    __device__ int eii_z[MAX_EII_SHELLS], eii_sh[MAX_EII_SHELLS], eii_nshells[MXELEMENT], eii_nsh[MXMED], 
	    eii_first[MXMED*MXEL], eii_no[MXMED*MXEL], eii_flag;
//#define eii_data1 eii_data

//COMMON/user_relax/
    __device__ float u_relax;
    __device__ int ish_relax, iz_relax;
//#define user_relax1 user_relax

//COMMON/et_control/
	//__device__ float smaxir[MXREG];
    __device__ float smaxir[MXMED];  //was smaxir[MXREG]
	__device__ float estepe, ximax, skindepth_for_bca;
    __device__ int transport_algorithm, bca_algorithm;
    __device__ bool exact_bca, spin_effects;
//#define et_control1 et_control

//COMMON/ms_data
    __device__ float ums_array[MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1*MAXU_MS_PLUS_1], fms_array[MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1*MAXU_MS_PLUS_1], wms_array[MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1*MAXU_MS_PLUS_1];
    __device__ short ims_array[MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1*MAXU_MS_PLUS_1];
    __device__ float llammin, llammax, dllamb, dllambi, dqms, dqmsi;
//#define ms_data1 ms_data

//COMMON/spin_data
    __device__ float spin_rej[MXMED*2*MAXE_SPI1_PLUS_1*MAXQ_SPIN_PLUS_1*MAXU_SPIN_PLUS_1], 
	    espin_min, espin_max, espml, b2spin_min, b2spin_max, 
	    dbeta2, dbeta2i, dlener, dleneri, dqq1, dqq1i;
    //__device__ bool fool_intel_optimizer;
//#define spin_data1 spin_data

/*
//COMMON/ch_steps
    __device__ long count_pii_steps, count_all_steps;
    __device__ bool is_ch_step;
//#define ch_steps1 ch_steps
*/

/*
//COMMON/epcont/
    __device__ float edep, tstep, tustep, ustep, tvstep, vstep, rhof, eold, enew, eke, 
	    elke, gle, e_range, x_final, y_final, z_final, u_final, 
	    v_final, w_final;
    __device__ int idisc, irold, irnew, iausfl[MXAUS];
//#define epcont_1 epcont_
*/

//COMMON/media/
    __device__ float rlc[MXMED];
	//__device__ float rldu[MXMED];
	//__device__ float rho[MXMED];
    //__device__ int msge[MXMED], mge[MXMED], mseke[MXMED], meke[MXMED], mleke[MXMED], mcmfp[MXMED], mrange[MXMED];
	//__device__ int iraylm[MXMED];
    //__device__ char media[24*MXMED], photon_xsections[16], comp_xsections[16];
    __device__ float apx, upx;
    __device__ char eii_xfile[16];
    //__device__ int nmed;
//#define media_1 media_

//COMMON/misc/
    //__device__ float dunit;
/*
    __device__ int kmpi, kmpo;
    __device__ float rhor[MXREG];
    __device__ short med[MXREG], iraylr[MXREG];
//#define misc_1 misc_
*/

//COMMON/photin/
    //__device__ float ebinda[MXMED]; 
    //__device__ float2 ge[MXMED], gmfp[MXGE*MXMED], gbr1[MXGE*MXMED], gbr2[MXGE*MXMED], cohe[MXGE*MXMED];
	//__device__ float2 rco[MXMED], rsct[MXRAYFF*MXMED];
	//__device__ float dpmfp;
    __device__ int mpgem[MXSGE*MXMED];
	//__device__ int ngr[MXMED];
//#define photin_1 photin_

/*
//COMMON/stack/
    __device__ float e[MXSTACK], x[MXSTACK], y[MXSTACK], z[MXSTACK], u[MXSTACK], v[MXSTACK], w[MXSTACK], 
	    dnear[MXSTACK], wt[MXSTACK];
    __device__ int iq[MXSTACK], ir[MXSTACK], latch[MXSTACK], latchi, np, npold;
//#define stack_1 stack_
*/

//COMMON/thresh/
    __device__ float rmt2, rmsq, ap[MXMED], ae[MXMED], up[MXMED], ue[MXMED];
	__device__ float te[MXMED], thmoll[MXMED];
//#define thresh_1 thresh_

/*
//COMMON/uphiin/
    __device__ float sinc0, sinc1; 
    __device__ float2 sin[MXSINC];
//#define uphiin_1 uphiin_

//COMMON/uphiot/
    __device__ float theta, sinthe, costhe, sinphi, cosphi, 
	__device__ float pi, twopi, pi5d2;
//#define uphiot_1 uphiot_
*/

//COMMON/useful/
    //__device__ float pzero;
	__device__ float prm, prmt2, rm;
    //__device__ int medium, medold;
//#define useful_1 useful_

//COMMON/etaly1/
    //__device__ float esum[4*MXREG*5]	;
//#define etaly1_1 etaly1_

//COMMON/ntaly1/
    //__device__ int nsum[4*MXREG*5]	;
//#define ntaly1_1 ntaly1_

//COMMON/rayleigh_inputs/
//Name of 101 of Media and file names, probably will not used in C
    //__device__ char iray_ff_media[2424], iray_ff_file[12928];
//#define rayleigh_inputs1 rayleigh_inputs

/*
//COMMON/rad_compton/
    __device__ float radc_sigs[RADC_NE_PLUS_1], radc_sigd[RADC_NE_PLUS_1], radc_frej[RADC_NE_PLUS_1*RADC_NU_PLUS_1],
        radc_x[RADC_NX], radc_fdat[RADC_NBOX], radc_smax[RADC_NBOX], 
	    radc_emin, radc_emax, radc_dw, radc_dle, radc_dlei, 
	    radc_le1;
    __device__ short radc_bins[RADC_NBOX], radc_ixmin1[RADC_NBOX], radc_ixmax1[RADC_NBOX], 
	    radc_ixmin2[RADC_NBOX], radc_ixmax2[RADC_NBOX], radc_ixmin3[RADC_NBOX], 
	    radc_ixmax3[RADC_NBOX], radc_ixmin4[RADC_NBOX], radc_ixmax4[RADC_NBOX], 
	    radc_startx[RADC_NE_PLUS_1], radc_startb[RADC_NE_PLUS_1];
//#define rad_compton1 rad_compton
*/

//COMMON/egs_vr/
    __device__ float e_max_rr[MXMED];  //the maximum energy that will do range rejection (rr)
	__device__ short i_do_rr[MXMED];   //which medium will do range rejection.

	__device__ int nbr_split;

	//__device__ float prob_rr;
	//__device__ int i_play_rr, i_survived_rr, n_rr_warning;  // the "rr" here means " Russian Roulette"
	
//#define egs_vr1 egs_vr
	__constant__ short first_transition[5]={1,20,27,33,38};                       //was static
	__constant__ short last_transition[5]={19,26,32,37,39};                       //was static
	__constant__ short final_state[39]={4,3,5,6,202,302,402,404,403,303,502,503,504,602,603,604,505,605,606,13,14,5,6,505,605,606,14,5,6,505,605,606,5,6,505,605,606,6,606};                               
