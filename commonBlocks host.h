//following are the variables defined in EGSnrc common Blocks.

//#include "f2c.h"
//#include "EGSnrc_Parameters.h"
/* Common Block Declarations */

//COMMON/BOUNDS/
    //float h_ecut[MXREG], h_pcut[MXREG], h_vacdst;
    
//COMMON/BREMPR/
    float h_dl1[8*MXMED], h_dl2[8*MXMED], 
	    h_dl3[8*MXMED], h_dl4[8*MXMED], 
	    h_dl5[8*MXMED], h_dl6[8*MXMED];
	//float h_alphi[2*MXMED];
	float h_bpar[2*MXMED];
	//float h_delpos[2*MXMED];
	float h_wa[MXMED*MXEL],
	    h_pz[MXMED*MXEL], h_zelem[MXMED*MXEL],
	    h_rhoz[MXMED*MXEL],
	    h_delcm[MXMED], h_zbrang[MXMED], h_lzbrang[MXMED];
	//float h_pwr2i[MXPWR2I];
    int h_nne[MXMED], h_ibrdst, h_iprdst;
	//int h_ibr_nist, h_pair_nrc, h_itriplet;
    //char h_asym[MXMED*MXEL*2]	;
//#define brempr_1 brempr_

/*
//COMMON/nist_brems
    float h_nb_fdata[MXBRXS_PLUS_1*MXBRES*MXMED], h_nb_xdata[MXBRXS_PLUS_1*MXBRES*MXMED], h_nb_wdata[MXBRXS*MXBRES*MXMED];
    int h_nb_idata[MXBRXS*MXBRES*MXMED]	;
    float h_nb_emin[MXMED], h_nb_emax[MXMED], h_nb_lemin[MXMED], h_nb_lemax[MXMED], 
	    h_nb_dle[MXMED], h_nb_dlei[MXMED];
*/
	float h_log_ap[MXMED];
//#define nist_brems1 nist_brems

/*
//COMMON/nrc_pair
    float h_nrcp_fdata[NRC_PAIR_NXX*NRC_PAIR_NEE*MXMED], h_nrcp_wdata[NRC_PAIR_NXX*NRC_PAIR_NEE*MXMED]	;
    int h_nrcp_idata[NRC_PAIR_NXX*NRC_PAIR_NEE*MXMED];
    float h_nrcp_xdata[NRC_PAIR_NXX], h_nrcp_emin, h_nrcp_emax, h_nrcp_dle, h_nrcp_dlei;
//#define nrc_pair1 nrc_pair
*/

/*
//COMMON/triplet_data/
    float h_a_triplet[MAX_TRIPLET*MXMED]	, h_b_triplet[MAX_TRIPLET*MXMED], 
	    h_dl_triplet, h_dli_triplet, h_bli_triplet, h_log_4rm;
//#define triplet_data1 triplet_data
*/

//Common/compton_data/
    int h_iz_array[MXTOTSH];
    float h_be_array[MXTOTSH], h_jo_array[MXTOTSH], h_erfjo_array[MXTOTSH];
    int h_ne_array[MXTOTSH], h_shn_array[MXTOTSH], h_shell_array[MXMDSH*MXMED];
    float h_eno_array[MXMDSH*MXMED];
    int h_eno_atbin_array[MXMDSH*MXMED], h_n_shell[MXMED];
	//int h_radc_flag;
    //short h_ibcmp[MXREG];
//#define compton_data1 compton_data

//COMMON/edge/
	float h_binding_energies[MXSHELL*MXELEMENT], 
		h_interaction_prob[MXSHELL*MXELEMENT], 
		h_relaxation_prob[MXTRANS*MXELEMENT], 
		h_edge_energies[MXEDGE*MXELEMENT];
	int h_edge_number[MXELEMENT];
	float h_edge_a[MXEDGE*MXELEMENT], h_edge_b[MXEDGE*MXELEMENT],	
		h_edge_c[MXEDGE*MXELEMENT],
		h_edge_d[MXEDGE*MXELEMENT];
	//short h_iedgfl[MXREG], h_iphter[MXREG];
//#define edge_1 (edge_._1)
//#define edge_2 (edge_._2)

//COMMON/elecin/
    float h_esige_max, h_psige_max, h_range_ep[2*MXEKE*MXMED], h_e_array[MXEKE*MXMED];
	//float h_esig_e[MXMED];
	//float h_psig_e[MXMED];
    //float2 h_etae_ms[MXEKE*MXMED];
	//float2 h_etap_ms[MXEKE*MXMED];
	//float2 h_q1ce_ms[MXEKE*MXMED];
	//float2 h_q1cp_ms[MXEKE*MXMED];
	//float2 h_q2ce_ms[MXEKE*MXMED];
	//float2 h_q2cp_ms[MXEKE*MXMED];

	float  h_sig_e[MXMED*2];
	float2 h_eta_ms[MXEKE*MXMED*2];
	float2 h_q1c_ms[MXEKE*MXMED*2];
	float2 h_q2c_ms[MXEKE*MXMED*2];

	float2 h_blcce[MXEKE*MXMED];
	float2 h_eke01[MXMED];

	//float h_xr0[MXMED], h_teff0[MXMED];
	float h_blcc[MXMED], h_xcc[MXMED];
	//float2 h_esig[MXEKE*MXMED];
	//float2 h_psig[MXEKE*MXMED];
	//float2 h_ededx[MXEKE*MXMED];
	//float2 h_pdedx[MXEKE*MXMED];
	float2 h_ebr1[MXEKE*MXMED];
	float2 h_pbr1[MXEKE*MXMED];
	float2 h_pbr2[MXEKE*MXMED];
	float2 h_tmxs[MXEKE*MXMED];

	float2 h_sig[MXEKE*MXMED*2];
	float2 h_dedx[MXEKE*MXMED*2];

	float h_expeke1[MXMED];
	//int h_iunrst[MXMED], h_epstfl[MXMED], h_iaprim[MXMED];
	bool h_sig_ismonotone[2*MXMED]	;
//#define elecin_1 elecin_

//COMMON/eii_data/
    float2 h_eii_xsection[MAX_EII_BINS];  //use to be h_eii_xsection_a, h_eii_xsection_b.
	float h_eii_cons[MXMED];
	float2 h_eii[MAX_EII_SHELLS];   //use to be h_eii_a, h_eii_b. now combine them to h_eii.x, h_eii.y.
    int h_eii_z[MAX_EII_SHELLS], h_eii_sh[MAX_EII_SHELLS], h_eii_nshells[MXELEMENT], h_eii_nsh[MXMED], 
	    h_eii_first[MXMED*MXEL], h_eii_no[MXMED*MXEL], h_eii_flag;
//#define eii_data1 eii_data

//COMMON/user_relax/
    float h_u_relax;
    int h_ish_relax, h_iz_relax;
//#define user_relax1 user_relax

//COMMON/et_control/
    //float h_smaxir[MXREG]
	float h_smaxir[MXMED];  //was h_smaxir[MXREG]
	float h_estepe, h_ximax, h_skindepth_for_bca;
    int h_transport_algorithm, h_bca_algorithm;
    bool h_exact_bca, h_spin_effects;
//#define et_control1 et_control

//COMMON/ms_data
    float h_ums_array[MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1*MAXU_MS_PLUS_1], h_fms_array[MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1*MAXU_MS_PLUS_1], h_wms_array[MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1*MAXU_MS_PLUS_1];
    short h_ims_array[MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1*MAXU_MS_PLUS_1];
    float h_llammin, h_llammax, h_dllamb, h_dllambi, h_dqms, h_dqmsi;
//#define ms_data1 ms_data

//COMMON/spin_data
    float h_spin_rej[MXMED*2*MAXE_SPI1_PLUS_1*MAXQ_SPIN_PLUS_1*MAXU_SPIN_PLUS_1], 
	    h_espin_min, h_espin_max, h_espml, h_b2spin_min, h_b2spin_max, 
	    h_dbeta2, h_dbeta2i, h_dlener, h_dleneri, h_dqq1, h_dqq1i;
    //bool h_fool_intel_optimizer;
//#define spin_data1 spin_data

/*
//COMMON/ch_steps
    double h_count_pii_steps, h_count_all_steps;
    bool h_is_ch_step;
//#define ch_steps1 ch_steps
*/

/*
//COMMON/epcont/
    float h_edep, h_tstep, h_tustep, h_ustep, h_tvstep, h_vstep, h_rhof, h_eold, h_enew,h_eke01, 
	    h_elke, h_gle, h_e_range, h_x_final, h_y_final, h_z_final, h_u_final, 
	    h_v_final, h_w_final;
    int h_idisc, h_irold, h_irnew, h_iausfl[MXAUS];
//#define epcont_1 epcont_
*/

//COMMON/media/
    float h_rlc[MXMED];
	//float h_rldu[MXMED];
	float h_rho[MXMED];
    //int h_msge[MXMED], h_mge[MXMED], h_mseke[MXMED];
	int h_meke[MXMED];
	//int h_mleke[MXMED], h_mcmfp[MXMED], h_mrange[MXMED];
	//int h_iraylm[MXMED];
	//char h_media[24*MXMED], h_photon_xsections[16], h_comp_xsections[16];
    float h_apx, h_upx;
    char h_eii_xfile[16];
    int h_nmed;
//#define media_1 media_

//COMMON/misc/
    float h_dunit;
/*
	int h_kmpi, h_kmpo;
	float h_rhor[MXREG];
	short h_med[MXREG], h_iraylr[MXREG];
//#define misc_1 misc_
*/

//COMMON/photin/
    //float h_ebinda[MXMED]; 
    //float2 h_ge[MXMED], h_gmfp[MXGE*MXMED],	h_gbr1[MXGE*MXMED], h_gbr2[MXGE*MXMED], h_cohe[MXGE*MXMED]; 
	//float2 h_rco[MXMED], h_rsct[MXRAYFF*MXMED];
	//float h_dpmfp;
    int h_mpgem[MXSGE*MXMED];
	//int h_ngr[MXMED];
//#define photin_1 photin_

/*
//COMMON/stack/
    float h_e[MXSTACK], h_x[MXSTACK], h_y[MXSTACK], h_z[MXSTACK], h_u[MXSTACK], h_v[MXSTACK], h_w[MXSTACK], 
	    h_dnear[MXSTACK], h_wt[MXSTACK];
    int h_iq[MXSTACK], h_ir[MXSTACK], h_latch[MXSTACK], h_latchi, h_np, h_npold;
//#define stack_1 stack_
*/

//COMMON/thresh/
    float h_rmt2, h_rmsq, h_ap[MXMED], h_ae[MXMED], h_up[MXMED];
	float h_ue[MXMED], h_te[MXMED], h_thmoll[MXMED];
//#define thresh_1 thresh_

/*
//COMMON/uphiin/
    float h_sinc0, h_sinc1; 
    float2 h_sin[MXSINC]; 
//#define uphiin_1 uphiin_

//COMMON/uphiot/
    float h_theta, h_sinthe, h_costhe, h_sinphi, h_cosphi, 
	float h_pi, h_twopi, h_pi5d2,h_dummy;
//#define uphiot_1 uphiot_
*/

//COMMON/useful/
    //float h_pzero;
	float h_prm, h_prmt2, h_rm;
    //int h_medium, h_medold;
//#define useful_1 useful_

//COMMON/etaly1/
    //float h_esum[4*MXREG*5]	;
//#define etaly1_1 etaly1_

//COMMON/ntaly1/
    //int h_nsum[4*MXREG*5]	;
//#define ntaly1_1 ntaly1_

//COMMON/rayleigh_inputs/
//Name of 101 of Media and file names, probably will not used in C
    //char h_iray_ff_media[2424], h_iray_ff_file[12928];
//#define rayleigh_inputs1 rayleigh_inputs

/*
//COMMON/rad_compton/
    float h_radc_sigs[RADC_NE_PLUS_1], h_radc_sigd[RADC_NE_PLUS_1], h_radc_frej[RADC_NE_PLUS_1*RADC_NU_PLUS_1],
        h_radc_x[RADC_NX], h_radc_fdat[RADC_NBOX], h_radc_smax[RADC_NBOX], 
	    h_radc_emin, h_radc_emax, h_radc_dw, h_radc_dle, h_radc_dlei, 
	    h_radc_le1;
    short h_radc_bins[RADC_NBOX], h_radc_ixmin1[RADC_NBOX], h_radc_ixmax1[RADC_NBOX], 
	    h_radc_ixmin2[RADC_NBOX], h_radc_ixmax2[RADC_NBOX], h_radc_ixmin3[RADC_NBOX], 
	    h_radc_ixmax3[RADC_NBOX], h_radc_ixmin4[RADC_NBOX], h_radc_ixmax4[RADC_NBOX], 
	    h_radc_startx[RADC_NE_PLUS_1], h_radc_startb[RADC_NE_PLUS_1];
//#define rad_compton1 rad_compton
*/

//COMMON/egs_vr/
	float h_e_max_rr[MXMED];  //the maximum energy that will do range rejection (rr)
	short h_i_do_rr[MXMED];   //which medium will do range rejection.

	int h_nbr_split;

	//float h_prob_rr;
	//int h_i_play_rr, h_i_survived_rr, h_n_rr_warning;  // the "rr" here means "Russian Roulette"
