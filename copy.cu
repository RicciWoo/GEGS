

int copyToDevice()
{

	cudaMemcpyToSymbol(dl1,h_dl1, sizeof(float)*8*MXMED);
	cudaMemcpyToSymbol(dl2,h_dl2, sizeof(float)*8*MXMED);
	cudaMemcpyToSymbol(dl3,h_dl3, sizeof(float)*8*MXMED);
	cudaMemcpyToSymbol(dl4,h_dl4, sizeof(float)*8*MXMED);
	cudaMemcpyToSymbol(dl5,h_dl5, sizeof(float)*8*MXMED);
	cudaMemcpyToSymbol(dl6,h_dl6, sizeof(float)*8*MXMED);

	//cudaMemcpyToSymbol(alphi,h_alphi , sizeof(float)*2*MXMED);
	cudaMemcpyToSymbol(bpar, h_bpar, sizeof(float)*2*MXMED);
	//cudaMemcpyToSymbol(delpos, h_delpos, sizeof(float)*2*MXMED);
	cudaMemcpyToSymbol(wa,h_wa , sizeof(float)*MXMED*MXEL);
	cudaMemcpyToSymbol(pz, h_pz, sizeof(float)*MXMED*MXEL);
	cudaMemcpyToSymbol(zelem, h_zelem, sizeof(float)*MXMED*MXEL);
	//cudaMemcpyToSymbol(rhoz, h_rhoz, sizeof(float)*MXMED*MXEL);
	//cudaMemcpyToSymbol(pwr2i, h_pwr2i, sizeof(float)*MXPWR2I);
	cudaMemcpyToSymbol(delcm,h_delcm , sizeof(float)*MXMED);
	cudaMemcpyToSymbol(zbrang,h_zbrang , sizeof(float)*MXMED);
	cudaMemcpyToSymbol(lzbrang, h_lzbrang, sizeof(float)*MXMED);
	cudaMemcpyToSymbol(nne, h_nne, sizeof(int)*MXMED);

	cudaMemcpyToSymbol(smaxir, h_smaxir, sizeof(float)*MXMED);

	cudaMemcpyToSymbol(ibrdst, &h_ibrdst, sizeof(int));
	cudaMemcpyToSymbol(iprdst, &h_iprdst, sizeof(int));
	//cudaMemcpyToSymbol(ibr_nist,&h_ibr_nist , sizeof(int));
	//cudaMemcpyToSymbol(pair_nrc,&h_pair_nrc , sizeof(int));
	//cudaMemcpyToSymbol(itriplet,&h_itriplet , sizeof(int));

	//cudaMemcpyToSymbol(asym,h_asym , sizeof(char)*MXMED*MXEL*2);
	//cudaMemcpyToSymbol(nb_fdata,h_nb_fdata , sizeof(float)*MXBRXS_PLUS_1*MXBRES*MXMED);
	//cudaMemcpyToSymbol(nb_xdata,h_nb_xdata , sizeof(float)*MXBRXS_PLUS_1*MXBRES*MXMED);
	//cudaMemcpyToSymbol(nb_wdata,h_nb_wdata , sizeof(float)*MXBRXS*MXBRES*MXMED);
	//cudaMemcpyToSymbol(nb_idata,h_nb_idata , sizeof(int)*MXBRXS*MXBRES*MXMED);
	//cudaMemcpyToSymbol(nb_emin,h_nb_emin , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(nb_emax,h_nb_emax , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(nb_lemin,h_nb_lemin , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(nb_lemax,h_nb_lemax , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(nb_dle,h_nb_dle , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(nb_dlei,h_nb_dlei , sizeof(float)*MXMED);
	cudaMemcpyToSymbol(log_ap,h_log_ap , sizeof(float)*MXMED);

	//cudaMemcpyToSymbol(nrcp_fdata,h_nrcp_fdata , sizeof(float)*NRC_PAIR_NXX*NRC_PAIR_NEE*MXMED);
	//cudaMemcpyToSymbol(nrcp_wdata,h_nrcp_wdata , sizeof(float)*NRC_PAIR_NXX*NRC_PAIR_NEE*MXMED);
	//cudaMemcpyToSymbol(nrcp_idata,h_nrcp_idata , sizeof(int)*NRC_PAIR_NXX*NRC_PAIR_NEE*MXMED);
	//cudaMemcpyToSymbol(nrcp_xdata,h_nrcp_xdata , sizeof(float)*NRC_PAIR_NXX);
	//cudaMemcpyToSymbol(nrcp_emin,&h_nrcp_emin , sizeof(float));
	//cudaMemcpyToSymbol(nrcp_emax,&h_nrcp_emax , sizeof(float));
	//cudaMemcpyToSymbol(nrcp_dle,&h_nrcp_dle , sizeof(float));
	//cudaMemcpyToSymbol(nrcp_dlei,&h_nrcp_dlei , sizeof(float));

	//cudaMemcpyToSymbol(a_triplet,h_a_triplet , sizeof(float)*MAX_TRIPLET*MXMED);
	//cudaMemcpyToSymbol(b_triplet,h_b_triplet , sizeof(float)*MAX_TRIPLET*MXMED);
	//cudaMemcpyToSymbol(dl_triplet,&h_dl_triplet , sizeof(float));
	//cudaMemcpyToSymbol(dli_triplet,&h_dli_triplet , sizeof(float));
	//cudaMemcpyToSymbol(bli_triplet,&h_bli_triplet , sizeof(float));
	//cudaMemcpyToSymbol(log_4rm,&h_log_4rm , sizeof(float));

	cudaMemcpyToSymbol(iz_array,h_iz_array , sizeof(int)*MXTOTSH);
	cudaMemcpyToSymbol(be_array,h_be_array , sizeof(float)*MXTOTSH);
	cudaMemcpyToSymbol(jo_array,h_jo_array , sizeof(float)*MXTOTSH);
	cudaMemcpyToSymbol(erfjo_array,h_erfjo_array , sizeof(float)*MXTOTSH);
	cudaMemcpyToSymbol(ne_array,h_ne_array , sizeof(int)*MXTOTSH);
	cudaMemcpyToSymbol(shn_array,h_shn_array , sizeof(int)*MXTOTSH);
	cudaMemcpyToSymbol(shell_array,h_shell_array , sizeof(int)*MXMDSH*MXMED);
	cudaMemcpyToSymbol(eno_array,h_eno_array , sizeof(float)*MXMDSH*MXMED);
	cudaMemcpyToSymbol(eno_atbin_array,h_eno_atbin_array , sizeof(int)*MXMDSH*MXMED);
	cudaMemcpyToSymbol(n_shell,h_n_shell , sizeof(int)*MXMED);
	//cudaMemcpyToSymbol(radc_flag,&h_radc_flag , sizeof(int));
	cudaMemcpyToSymbol(binding_energies,h_binding_energies , sizeof(float)*MXSHELL*MXELEMENT);
	cudaMemcpyToSymbol(interaction_prob,h_interaction_prob , sizeof(float)*MXSHELL*MXELEMENT);
	cudaMemcpyToSymbol(relaxation_prob,h_relaxation_prob , sizeof(float)*MXTRANS*MXELEMENT);
	cudaMemcpyToSymbol(edge_energies,h_edge_energies , sizeof(float)*MXEDGE*MXELEMENT);
	cudaMemcpyToSymbol(edge_number,h_edge_number , sizeof(int)*MXELEMENT);
	cudaMemcpyToSymbol(edge_a,h_edge_a , sizeof(float)*MXEDGE*MXELEMENT);
	cudaMemcpyToSymbol(edge_b,h_edge_b , sizeof(float)*MXEDGE*MXELEMENT);
	cudaMemcpyToSymbol(edge_c,h_edge_c , sizeof(float)*MXEDGE*MXELEMENT);
	cudaMemcpyToSymbol(edge_d,h_edge_d , sizeof(float)*MXEDGE*MXELEMENT);

	//cudaMemcpyToSymbol(esig_e,h_esig_e , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(psig_e,h_psig_e , sizeof(float)*MXMED);
	cudaMemcpyToSymbol(esige_max,&h_esige_max , sizeof(float));
	cudaMemcpyToSymbol(psige_max,&h_psige_max , sizeof(float));
	cudaMemcpyToSymbol(range_ep,h_range_ep , sizeof(float)*2*MXEKE*MXMED);
	cudaMemcpyToSymbol(e_array,h_e_array , sizeof(float)*MXEKE*MXMED);
	//cudaMemcpyToSymbol(etae_ms,h_etae_ms , sizeof(float2)*MXEKE*MXMED);
	//cudaMemcpyToSymbol(etap_ms,h_etap_ms , sizeof(float2)*MXEKE*MXMED);
	//cudaMemcpyToSymbol(q1ce_ms,h_q1ce_ms , sizeof(float2)*MXEKE*MXMED);
	//cudaMemcpyToSymbol(q1cp_ms,h_q1cp_ms , sizeof(float2)*MXEKE*MXMED);
	//cudaMemcpyToSymbol(q2ce_ms,h_q2ce_ms , sizeof(float2)*MXEKE*MXMED);
	//cudaMemcpyToSymbol(q2cp_ms,h_q2cp_ms , sizeof(float2)*MXEKE*MXMED);

	cudaMemcpyToSymbol(sig_e,h_sig_e , sizeof(float)*MXMED*2);
	cudaMemcpyToSymbol(eta_ms,h_eta_ms , sizeof(float2)*MXEKE*MXMED*2);
	cudaMemcpyToSymbol(q1c_ms,h_q1c_ms , sizeof(float2)*MXEKE*MXMED*2);
	cudaMemcpyToSymbol(q2c_ms,h_q2c_ms , sizeof(float2)*MXEKE*MXMED*2);

	cudaMemcpyToSymbol(blcce,h_blcce , sizeof(float2)*MXEKE*MXMED);
	cudaMemcpyToSymbol(eke01,h_eke01 , sizeof(float2)*MXMED);
	//cudaMemcpyToSymbol(xr0,h_xr0 , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(teff0,h_teff0 , sizeof(float)*MXMED);
	cudaMemcpyToSymbol(blcc,h_blcc , sizeof(float)*MXMED);
	cudaMemcpyToSymbol(xcc,h_xcc , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(esig,h_esig , sizeof(float2)*MXEKE*MXMED);
	//cudaMemcpyToSymbol(psig,h_psig , sizeof(float2)*MXEKE*MXMED);
	cudaMemcpyToSymbol(sig,h_sig , sizeof(float2)*MXEKE*MXMED*2);
	//cudaMemcpyToSymbol(ededx,h_ededx , sizeof(float2)*MXEKE*MXMED);
	//cudaMemcpyToSymbol(pdedx,h_pdedx , sizeof(float2)*MXEKE*MXMED);
	cudaMemcpyToSymbol(dedx,h_dedx , sizeof(float2)*MXEKE*MXMED*2);
	cudaMemcpyToSymbol(ebr1,h_ebr1 , sizeof(float2)*MXEKE*MXMED);
	cudaMemcpyToSymbol(pbr1,h_pbr1 , sizeof(float2)*MXEKE*MXMED);
	cudaMemcpyToSymbol(pbr2,h_pbr2 , sizeof(float2)*MXEKE*MXMED);
	cudaMemcpyToSymbol(tmxs,h_tmxs , sizeof(float2)*MXEKE*MXMED);
	cudaMemcpyToSymbol(expeke1,h_expeke1 , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(iunrst,h_iunrst , sizeof(int)*MXMED);
	//cudaMemcpyToSymbol(epstfl,h_epstfl , sizeof(int)*MXMED);
	//cudaMemcpyToSymbol(iaprim,h_iaprim , sizeof(int)*MXMED);
	cudaMemcpyToSymbol(sig_ismonotone,h_sig_ismonotone , sizeof(bool)*2*MXMED);
	cudaMemcpyToSymbol(eii_xsection,h_eii_xsection , sizeof(float2)*MAX_EII_BINS);
	cudaMemcpyToSymbol(eii_cons,h_eii_cons , sizeof(float)*MXMED);
	cudaMemcpyToSymbol(eii,h_eii , sizeof(float2)*MAX_EII_SHELLS);
	cudaMemcpyToSymbol(eii_z,h_eii_z , sizeof(int)*MAX_EII_SHELLS);
	cudaMemcpyToSymbol(eii_sh,h_eii_sh , sizeof(int)*MAX_EII_SHELLS);
	cudaMemcpyToSymbol(eii_nshells,h_eii_nshells , sizeof(int)*MXELEMENT);
	cudaMemcpyToSymbol(eii_nsh,h_eii_nsh , sizeof(int)*MXMED);
	cudaMemcpyToSymbol(eii_first,h_eii_first , sizeof(int)*MXMED*MXEL);
	cudaMemcpyToSymbol(eii_no,h_eii_no , sizeof(int)*MXMED*MXEL);
	cudaMemcpyToSymbol(eii_flag,&h_eii_flag , sizeof(int));
	cudaMemcpyToSymbol(u_relax,&h_u_relax , sizeof(float));
	cudaMemcpyToSymbol(ish_relax,&h_ish_relax , sizeof(int));
	cudaMemcpyToSymbol(iz_relax,&h_iz_relax , sizeof(int));
	cudaMemcpyToSymbol(transport_algorithm,&h_transport_algorithm , sizeof(int));
	cudaMemcpyToSymbol(bca_algorithm,&h_bca_algorithm , sizeof(int));
	cudaMemcpyToSymbol(exact_bca,&h_exact_bca , sizeof(bool));
	cudaMemcpyToSymbol(spin_effects,&h_spin_effects , sizeof(bool));
	cudaMemcpyToSymbol(estepe,&h_estepe , sizeof(float));
	cudaMemcpyToSymbol(ximax,&h_ximax , sizeof(float));
	cudaMemcpyToSymbol(skindepth_for_bca,&h_skindepth_for_bca , sizeof(float));
	cudaMemcpyToSymbol(ums_array,h_ums_array , sizeof(float)*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1*MAXU_MS_PLUS_1);
	cudaMemcpyToSymbol(fms_array,h_fms_array , sizeof(float)*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1*MAXU_MS_PLUS_1);
	cudaMemcpyToSymbol(wms_array,h_wms_array , sizeof(float)*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1*MAXU_MS_PLUS_1);
	cudaMemcpyToSymbol(ims_array,h_ims_array , sizeof(short)*MAXL_MS_PLUS_1*MAXQ_MS_PLUS_1*MAXU_MS_PLUS_1);

	cudaMemcpyToSymbol(spin_rej,h_spin_rej , sizeof(float)*MXMED*2*MAXE_SPI1_PLUS_1*MAXQ_SPIN_PLUS_1*MAXU_SPIN_PLUS_1);


	cudaMemcpyToSymbol(llammin,&h_llammin , sizeof(float));
	cudaMemcpyToSymbol(llammax,&h_llammax , sizeof(float));
	cudaMemcpyToSymbol(dllamb,&h_dllamb , sizeof(float));
	cudaMemcpyToSymbol(dllambi,&h_dllambi , sizeof(float));
	cudaMemcpyToSymbol(dqms,&h_dqms , sizeof(float));
	cudaMemcpyToSymbol(dqmsi,&h_dqmsi , sizeof(float));
	cudaMemcpyToSymbol(espin_min,&h_espin_min , sizeof(float));
	cudaMemcpyToSymbol(espin_max,&h_espin_max , sizeof(float));
	cudaMemcpyToSymbol(espml,&h_espml , sizeof(float));
	cudaMemcpyToSymbol(b2spin_min,&h_b2spin_min , sizeof(float));
	cudaMemcpyToSymbol(b2spin_max,&h_b2spin_max , sizeof(float));
	cudaMemcpyToSymbol(dbeta2,&h_dbeta2 , sizeof(float));
	cudaMemcpyToSymbol(dbeta2i,&h_dbeta2i , sizeof(float));
	cudaMemcpyToSymbol(dlener,&h_dlener , sizeof(float));
	cudaMemcpyToSymbol(dleneri,&h_dleneri , sizeof(float));
	cudaMemcpyToSymbol(dqq1,&h_dqq1 , sizeof(float));
	cudaMemcpyToSymbol(dqq1i,&h_dqq1i , sizeof(float));

	//cudaMemcpyToSymbol(count_pii_steps,&h_count_pii_steps , sizeof(long));
	//cudaMemcpyToSymbol(count_all_steps,&h_count_all_steps , sizeof(long));
	//cudaMemcpyToSymbol(is_ch_step,&h_is_ch_step , sizeof(bool));
	//cudaMemcpyToSymbol(edep,&h_edep , sizeof(float));
	//cudaMemcpyToSymbol(tstep,&h_tstep , sizeof(float));
	//cudaMemcpyToSymbol(tustep,&h_tustep , sizeof(float));
	//cudaMemcpyToSymbol(ustep,&h_ustep , sizeof(float));
	//cudaMemcpyToSymbol(tvstep,&h_tvstep , sizeof(float));
	//cudaMemcpyToSymbol(vstep,&h_vstep , sizeof(float));
	//cudaMemcpyToSymbol(rhof,&h_rhof , sizeof(float));
	//cudaMemcpyToSymbol(eold,&h_eold , sizeof(float));
	//cudaMemcpyToSymbol(enew,&h_enew , sizeof(float));
	//cudaMemcpyToSymbol(eke01,&h_eke01 , sizeof(float));
	//cudaMemcpyToSymbol(elke,&h_elke , sizeof(float));
	//cudaMemcpyToSymbol(gle,&h_gle , sizeof(float));
	//cudaMemcpyToSymbol(e_range,&h_e_range , sizeof(float));
	//cudaMemcpyToSymbol(x_final,&h_x_final , sizeof(float));
	//cudaMemcpyToSymbol(y_final,&h_y_final , sizeof(float));
	//cudaMemcpyToSymbol(z_final,&h_z_final , sizeof(float));
	//cudaMemcpyToSymbol(u_final,&h_u_final , sizeof(float));
	//cudaMemcpyToSymbol(v_final,&h_v_final , sizeof(float));
	//cudaMemcpyToSymbol(w_final,&h_w_final , sizeof(float));
	//cudaMemcpyToSymbol(idisc,&h_idisc , sizeof(int));
	//cudaMemcpyToSymbol(irold,&h_irold , sizeof(int));
	//cudaMemcpyToSymbol(irnew,&h_irnew , sizeof(int));
	//cudaMemcpyToSymbol(iausfl,h_iausfl , sizeof(int)*MXAUS);

	//cudaMemcpyToSymbol(rlc,h_rlc , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(rldu,h_rldu , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(rho,h_rho , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(msge,h_msge , sizeof(int)*MXMED);
	//cudaMemcpyToSymbol(mge,h_mge , sizeof(int)*MXMED);
	//cudaMemcpyToSymbol(mseke,h_mseke , sizeof(int)*MXMED);
	//cudaMemcpyToSymbol(meke,h_meke , sizeof(int)*MXMED);
	//cudaMemcpyToSymbol(mleke,h_mleke , sizeof(int)*MXMED);
	//cudaMemcpyToSymbol(mcmfp,h_mcmfp , sizeof(int)*MXMED);
	//cudaMemcpyToSymbol(mrange,h_mrange , sizeof(int)*MXMED);

	//cudaMemcpyToSymbol(iraylm,h_iraylm , sizeof(int)*MXMED);
	//cudaMemcpyToSymbol(media,h_media , sizeof(char)*24*MXMED);
	//cudaMemcpyToSymbol(photon_xsections,h_photon_xsections , sizeof(char)*16);
	cudaMemcpyToSymbol(apx,&h_apx , sizeof(float));
	cudaMemcpyToSymbol(upx,&h_upx , sizeof(float));
	//cudaMemcpyToSymbol(nmed,&h_nmed , sizeof(int));
	//cudaMemcpyToSymbol(dunit,&h_dunit , sizeof(float));
	//cudaMemcpyToSymbol(kmpi,h_kmpi , sizeof(int));
	//cudaMemcpyToSymbol(kmpo,h_kmpo , sizeof(int));
	//cudaMemcpyToSymbol(ebinda,h_ebinda , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(ge,h_ge , sizeof(float2)*MXMED);
	//cudaMemcpyToSymbol(gbr1,h_gbr1 , sizeof(float2)*MXGE*MXMED);
	//cudaMemcpyToSymbol(gmfp,h_gmfp , sizeof(float2)*MXGE*MXMED);
	//cudaMemcpyToSymbol(gbr2,h_gbr2 , sizeof(float2)*MXGE*MXMED);
	//cudaMemcpyToSymbol(rco,h_rco , sizeof(float2)*MXMED);
	//cudaMemcpyToSymbol(rsct,h_rsct , sizeof(float2)*MXRAYFF*MXMED);
	//cudaMemcpyToSymbol(cohe,h_cohe , sizeof(float2)*MXGE*MXMED);
	//cudaMemcpyToSymbol(dpmfp,&h_dpmfp , sizeof(float));
	cudaMemcpyToSymbol(mpgem,h_mpgem , sizeof(int)*MXSGE*MXMED);
	//cudaMemcpyToSymbol(ngr,h_ngr , sizeof(int)*MXMED);
	cudaMemcpyToSymbol(rmt2,&h_rmt2 , sizeof(float));
	cudaMemcpyToSymbol(rmsq,&h_rmsq , sizeof(float));
	cudaMemcpyToSymbol(ap,h_ap , sizeof(float)*MXMED);
	cudaMemcpyToSymbol(ae,h_ae , sizeof(float)*MXMED);
	cudaMemcpyToSymbol(up,h_up , sizeof(float)*MXMED);
	cudaMemcpyToSymbol(ue,h_ue , sizeof(float)*MXMED);
	cudaMemcpyToSymbol(te,h_te , sizeof(float)*MXMED);
	cudaMemcpyToSymbol(thmoll,h_thmoll , sizeof(float)*MXMED);
	//cudaMemcpyToSymbol(pzero,&h_pzero , sizeof(float));
	cudaMemcpyToSymbol(prm,&h_prm , sizeof(float));
	cudaMemcpyToSymbol(prmt2,&h_prmt2 , sizeof(float));
	cudaMemcpyToSymbol(rm,&h_rm , sizeof(float));
	//cudaMemcpyToSymbol(medium,&h_medium , sizeof(int));
	//cudaMemcpyToSymbol(medold,&h_medold , sizeof(int));

	//cudaMemcpyToSymbol(radc_emin,&h_radc_emin , sizeof(float));
	//cudaMemcpyToSymbol(radc_emax,&h_radc_emax , sizeof(float));
	//cudaMemcpyToSymbol(radc_dw,&h_radc_dw , sizeof(float));
	//cudaMemcpyToSymbol(radc_dle,&h_radc_dle , sizeof(float));
	//cudaMemcpyToSymbol(radc_dlei,&h_radc_dlei , sizeof(float));
	//cudaMemcpyToSymbol(radc_le1,&h_radc_le1 , sizeof(float));
	//cudaMemcpyToSymbol(radc_sigs,&h_radc_sigs , sizeof(float)*RADC_NE_PLUS_1);
	//cudaMemcpyToSymbol(radc_sigd,h_radc_sigd , sizeof(float)*RADC_NE_PLUS_1);
	//cudaMemcpyToSymbol(radc_frej,h_radc_frej , sizeof(float)*RADC_NE_PLUS_1*RADC_NU_PLUS_1);
	//cudaMemcpyToSymbol(radc_x,h_radc_x , sizeof(float)*RADC_NX);
	//cudaMemcpyToSymbol(radc_fdat,h_radc_fdat , sizeof(float)*RADC_NBOX);
	//cudaMemcpyToSymbol(radc_smax,h_radc_smax , sizeof(float)*RADC_NBOX);
	//cudaMemcpyToSymbol(radc_bins,h_radc_bins , sizeof(short)*RADC_NBOX);
	//cudaMemcpyToSymbol(radc_ixmin1,h_radc_ixmin1 , sizeof(short)*RADC_NBOX);
	//cudaMemcpyToSymbol(radc_ixmax1,h_radc_ixmax1 , sizeof(short)*RADC_NBOX);
	//cudaMemcpyToSymbol(radc_ixmin2,h_radc_ixmin2 , sizeof(short)*RADC_NBOX);
	//cudaMemcpyToSymbol(radc_ixmax2,h_radc_ixmax2 , sizeof(short)*RADC_NBOX);
	//cudaMemcpyToSymbol(radc_ixmin3,h_radc_ixmin3 , sizeof(short)*RADC_NBOX);
	//cudaMemcpyToSymbol(radc_ixmax3,h_radc_ixmax3 , sizeof(short)*RADC_NBOX);
	//cudaMemcpyToSymbol(radc_ixmin4,h_radc_ixmin4 , sizeof(short)*RADC_NBOX);
	//cudaMemcpyToSymbol(radc_ixmax4,h_radc_ixmax4 , sizeof(short)*RADC_NBOX);
	//cudaMemcpyToSymbol(radc_startx,h_radc_startx , sizeof(short)*RADC_NE_PLUS_1);
	//cudaMemcpyToSymbol(radc_startb,h_radc_startb , sizeof(short)*RADC_NE_PLUS_1);

	cudaMemcpyToSymbol(i_do_rr, &h_i_do_rr , sizeof(short)*MXMED);
	cudaMemcpyToSymbol(e_max_rr, &h_e_max_rr , sizeof(float)*MXMED);


	//cudaMemcpyToSymbol(prob_rr,&h_prob_rr , sizeof(float));
	cudaMemcpyToSymbol(nbr_split,&h_nbr_split , sizeof(int));
	//cudaMemcpyToSymbol(i_survived_rr,&h_i_survived_rr , sizeof(int));
	//cudaMemcpyToSymbol(n_rr_warning,&h_n_rr_warning , sizeof(int));

	return 1;
}

