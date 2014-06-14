/*
**  $Id: gsw_check_functions.c,v d7a5468a0b8c 2014/06/14 02:35:01 fmd $
**  $Version: 3.03.0 $
*/
#include <gswteos-10.h>

static void	report(char *what, double acceptable, double actual,
			double expected);

/*
** Acceptable differences from check values:
*/
#define abs_pressure_from_p_ca  2.300031483173370e-004
#define adiabatic_lapse_rate_from_ct_ca  5.699743845487985e-019
#define adiabatic_lapse_rate_from_t_ca  5.699738675609156e-019
#define adiabatic_lapse_rate_t_exact_ca  5.699738675609156e-019
#define alpha_ca  8.264713918662917e-015
#define alpha_ct_exact_ca  8.257102480598889e-015
#define alpha_ct_exact_rab_ca  8.257102480598889e-015
#define alpha_on_beta_ca 1.052907760978883e-11
#define alpha_rab_ca  8.264713918662917e-015
#define alpha_wrt_ct_t_exact_ca  8.257094010269417e-015
#define alpha_wrt_pt_t_exact_ca  8.599068316130290e-015
#define alpha_wrt_t_exact_ca  8.594856868316542e-015
#define beta_ca  1.846179459308317e-015
#define beta_const_ct_t_exact_ca  1.847372081698051e-015
#define beta_const_pt_t_exact_ca  1.805738718274608e-015
#define beta_const_t_exact_ca  1.804871356536619e-015
#define beta_ct_exact_ca  1.847372081698051e-015
#define beta_ct_exact_rab_ca  1.847372081698051e-015
#define beta_rab_ca  1.846179459308317e-015
#define brinesa_ct_ca  1.300080043620255e-010
#define brinesa_t_ca  1.300293206440983e-010
#define cabbeling_ca  1.634722766223964e-16
#define c_from_sp_ca  6.163816124171717e-010
#define chem_potential_salt_t_exact_ca  4.805315256817266e-007
#define chem_potential_t_exact_ca  1.317864928296331e-009
#define chem_potential_water_t_exact_ca  4.811790859093890e-007
#define cp_t_exact_ca  2.813976607285440e-009
#define ct_freezing_ca  2.257127817983928e-011
#define ct_from_entropy_ca  6.261107188265669e-010
#define ct_from_pt_ca  6.261107188265669e-010
#define ct_from_rho_ca  6.298552790440226e-010
#define ct_from_rho_exact_ca  6.280096442878858e-010
#define ct_from_t_ca  6.261107188265669e-010
#define ct_maxdensity_ca  6.688338771709823e-011
#define ct_maxdensity_exact_ca  5.839062566792563e-011
#define ct_pt_ca  2.964295475749168e-013
#define ct_pt_pt_ca  5.551115123125783e-014
#define ct_sa_ca  1.006231122729906e-012
#define ct_sa_pt_ca  1.457167719820518e-014
#define ct_sa_sa_ca  1.431146867680866e-014
#define deltasa_from_sp_ca  6.963318810448982e-013
#define deltasa_atlas_ca  6.945514042372425e-013
#define depth_from_z_ca  2.287223921371151e-008
#define distance_ca  4.470348358154297e-008
#define drho_dct_ca 8.315403920988729e-12
#define drho_dp_ca 1.782157321023048e-18
#define drho_dsa_ca 1.943112337698949e-12
#define dynamic_enthalpy_ca  2.288754734930485e-007
#define dynamic_enthalpy_ct_exact_ca  2.288797986693680e-007
#define dynamic_enthalpy_t_exact_ca  2.288943505845964e-007
#define enthalpy_ca  2.499356924090534e-006
#define enthalpy_ct_exact_ca  2.499349648132920e-006
#define enthalpy_diff_ca  1.154347728515859e-010
#define enthalpy_diff_ct_exact_ca  7.275957614183426e-011
#define enthalpy_t_exact_ca  2.499349648132920e-006
#define entropy_from_ct_ca  9.028163105995191e-009
#define entropy_from_pt_ca  9.028163105995191e-009
#define entropy_from_t_ca  9.028163105995191e-009
#define entropy_t_exact_ca  9.028163105995191e-009
#define eta_ct_ca  3.137579085432662e-011
#define eta_ct_ct_ca  2.381358998881922e-013
#define eta_sa_ca  4.653527563291959e-012
#define eta_sa_ct_ca  4.981922535098049e-014
#define eta_sa_sa_ca  6.995931611797346e-013
#define f_ca  1.000000000000000e-015
#define fdelta_ca  2.702939055302528e-014
#define g_ct_ca  2.854610769986721e-009
#define geo_strf_cunningham_ca  9.926097845891491e-007
#define geo_strf_dyn_height_ca  9.969444363377988e-007
#define geo_strf_dyn_height_pc_ca  4.210555459849275e-008
#define geo_strf_dyn_height_pc_p_mid_ca  1.000000000000000e-015
#define geo_strf_isopycnal_ca  4.297791695861974e-007
#define geo_strf_isopycnal_pc_ca  8.540336438045415e-009
#define geo_strf_isopycnal_pc_p_mid_ca  1.000000000000000e-015
#define geo_strf_montgomery_ca  9.915690188933013e-007
#define geo_strf_velocity_ca  8.024344404222727e-009
#define geo_strf_velocity_mid_lat_ca  1.000000000000000e-015
#define geo_strf_velocity_mid_long_ca  1.000000000000000e-015
#define grav_ca  5.329070518200751e-014
#define h_ct_ca  3.160494088660926e-010
#define h_ct_ct_ca  1.135647131889073e-011
#define helmholtz_energy_t_exact_ca  2.440137905068696e-007
#define h_p_ca  2.818925648462312e-016
#define h_sa_ca  2.371365326325758e-010
#define h_sa_ct_ca  2.188554892867956e-012
#define h_sa_sa_ca  7.264189250122399e-013
#define internal_energy_ca  2.499342372175306e-006
#define internal_energy_ct_exact_ca  2.499335096217692e-006
#define internal_energy_t_exact_ca  2.499335096217692e-006
#define ionic_strength_from_sa_ca  2.768674178810215e-012
#define ipvfn2_ca  3.474816878679121e-009
#define isochoric_heat_cap_t_exact_ca  1.614353095646948e-009
#define isopycnal_slope_ratio_ca  3.781384094736495e-010
#define kappa_ca 1.717743939542931e-21
#define kappa_const_t_exact_ca  1.697064424229105e-021
#define kappa_t_exact_ca  1.712677458291044e-021
#define latentheat_evap_ct_ca  1.455657184123993e-006
#define latentheat_evap_t_ca  1.443084329366684e-006
#define latentheat_melting_ca  6.286427378654480e-008
#define molality_from_sa_ca  4.446665258228677e-012
#define n2_ca  1.578186366313350e-014
#define ntpptct_ca  5.024869409453459e-013
#define osmotic_coefficient_t_exact_ca  3.583799923490005e-013
#define osmotic_pressure_t_exact_ca  1.023465756588848e-009
#define p_from_abs_pressure_ca  2.300066626048647e-008
#define p_from_z_ca  2.301931090187281e-008
#define p_mid_g_ct_ca  2.300021151313558e-008
#define p_mid_ipvfn2_ca  2.300021151313558e-008
#define p_mid_n2_ca  2.300021151313558e-008
#define p_mid_tursr_ca  2.300021151313558e-008
#define pot_enthalpy_from_pt_ca  2.499356924090534e-006
#define pot_rho_t_exact_ca  2.929709808086045e-010
#define pt0_from_t_ca  6.054037271496782e-010
#define pt_ct_ca  2.733369086627135e-013
#define pt_ct_ct_ca  4.440892098500626e-014
#define pt_from_ct_ca  6.054037271496782e-010
#define pt_from_entropy_ca  6.054072798633570e-010
#define pt_from_t_ca  6.054037271496782e-010
#define pt_sa_ca  9.670458878119348e-013
#define pt_sa_ct_ca  1.199127602768968e-014
#define pt_sa_sa_ca  2.081668171172169e-014
#define r_from_sp_ca  1.436317731418058e-011
#define rho_ca  2.945625965367071e-010
#define rho_ct_exact_ca  2.944489096989855e-010
#define rho_ct_exact_rab_ca  2.944489096989855e-010
#define rho_rab_ca  2.944489096989855e-010
#define rho_t_exact_ca  2.944489096989855e-010
#define rsubrho_ca  1.709803143512545e-008
#define sa_from_rho_ca  1.308677610722953e-010
#define sa_from_rho_ct_exact_ca  1.306119656874216e-010
#define sa_from_rho_t_exact_ca  1.304769625676272e-010
#define sa_from_sp_ca  1.300080043620255e-010
#define sa_from_sstar_ca  1.300222152167407e-010
#define sa_sa_sstar_from_sp_ca  1.300008989346679e-010
#define sigma0_ca  2.933120413217694e-010
#define sigma0_ct_exact_ca  2.929709808086045e-010
#define sigma0_pt0_exact_ca  2.929709808086045e-010
#define sigma1_ca  2.999058779096231e-010
#define sigma1_ct_exact_ca  2.994511305587366e-010
#define sigma2_ca  3.060449671465904e-010
#define sigma2_ct_exact_ca  3.055902197957039e-010
#define sigma3_ca  3.119566827081144e-010
#define sigma3_ct_exact_ca  3.117293090326712e-010
#define sigma4_ca  3.180957719450817e-010
#define sigma4_ct_exact_ca  3.176410245941952e-010
#define sound_speed_ca  2.596152626210824e-009
#define sound_speed_ct_exact_ca  2.590240910649300e-009
#define sound_speed_t_exact_ca  2.590240910649300e-009
#define specvol_anom_ca  2.810252031082428e-016
#define specvol_anom_ct_exact_ca  2.805915222392486e-016
#define specvol_anom_t_exact_ca  2.805915222392486e-016
#define specvol_ca  2.821094052807283e-016
#define specvol_ct_exact_ca  2.818925648462312e-016
#define specvol_t_exact_ca  2.818925648462312e-016
#define sp_from_c_ca  1.297193463756230e-010
#define sp_from_r_ca  2.681299626772216e-012
#define sp_from_sa_ca  1.297113527698457e-010
#define sp_from_sr_ca  1.297122409482654e-010
#define sp_from_sstar_ca  1.297122409482654e-010
#define sp_salinometer_ca  1.297131291266851e-010
#define sr_from_sp_ca  1.303233077010191e-010
#define sstar_from_sa_ca  1.300008989346679e-010
#define sstar_from_sp_ca  1.300008989346679e-010
#define sstar_sa_sstar_from_sp_ca  1.300008989346679e-010
#define steric_height_ca  1.017674460257467e-007
#define t90_from_t48_ca  5.997407015456702e-010
#define t90_from_t68_ca  5.998579410970706e-010
#define t_freezing_ca  2.157829470661454e-011
#define t_from_ct_ca  6.000142604989378e-010
#define t_from_rho_exact_ca  6.032974120273593e-010
#define thermobaric_ca  4.890907320111513e-23
#define t_maxdensity_exact_ca  6.274447628129565e-011
#define tu_ca  2.190718007000214e-008
#define z_from_depth_ca  2.287223921371151e-008
#define z_from_p_ca  2.287223921371151e-008

int	gsw_error_flag=0, debug;

/*
**  Main
*/
int
main(int argc, char **argv)
{
	int	nz=3;
	double	sp, sa, sstar, sr, t, ct, pt, p, p_bs, p_ref, c, rho;
	double	lon, long_bs, lat, lat_bs, saturation_fraction;
	double	cabbeling, thermobaric, alpha_on_beta;
	double	sigma0, sigma1, sigma2, sigma3, sigma4;
	double	drho_dsa, drho_dct, drho_dp, drho_dsa_error, drho_dct_error;
	double	drho_dp_error;

	double	n2[2], p_mid_n2[2], n2_error[2], p_mid_n2_error[2];
	double	ipvfn2[2],p_mid_ipvfn2[2];
	double	ipvfn2_error[2], p_mid_ipvfn2_error[2];
	double	tu[2], rsubrho[2], p_mid_tursr[2];
	double	tu_error[2], rsubrho_error[2], p_mid_tursr_error[2];
	double	sa_profile[3], ct_profile[3], p_profile[3], lat_profile[3];

	sa_profile[0] = 35.5e0;
	sa_profile[1] = 35.7e0;
	sa_profile[2] = 35.6e0;
	ct_profile[0] = 12.5e0;
	ct_profile[1] = 15e0;
	ct_profile[2] = 10e0;
	p_profile[0] = 00e0;
	p_profile[1] = 50e0;
	p_profile[2] = 100e0;
	lat_profile[0] = 10e0;
	lat_profile[1] = 10e0;
	lat_profile[2] = 10e0;

	sp	=  35.5e0;
	sa	= 35.7e0;
	sstar	= 35.5e0;
	sr	= 35.5e0;
	t	= 15e0;
	ct	= 20e0;
	pt	= 15e0;
	p	= 300e0;
	p_bs	= 50e0;
	p_ref	= 100e0;
	c	= 43.6e0;
	lon	= 260e0;
	long_bs	= 20e0;
	lat	= 16e0;
	lat_bs	= 60e0;
	saturation_fraction	= 0.5e0;

	if (argc==2 && !strcmp(argv[1],"-debug"))
	    debug	= 1;

	printf(
"============================================================================\n"
" Gibbs SeaWater (GSW) Oceanographic Toolbox of TEOS-10 version 3.02 (C)\n"
"============================================================================\n"
"\n"
" These are the check values for the subset of functions that have been \n"
" converted into C from the Gibbs SeaWater (GSW) Oceanographic Toolbox \n"
" of TEOS-10 (version 3.03).\n"
	);
	printf(
"\nPractical Salinity, PSS-78:\n\n"
	);
	report("gsw_sp_from_c", sp_from_c_ca,
		gsw_sp_from_c(c,t,p), 35.500961780774482e0);
	report("gsw_c_from_sp", c_from_sp_ca,
		gsw_c_from_sp(sp,t,p), 43.598945605280484e0);
	printf(
"\nAbsolute Salinity, Preformed Salinity and Conservative Temperature:\n\n"
	);
	report("gsw_sa_from_sp", sa_from_sp_ca,
		gsw_sa_from_sp(sp,p,lon,lat), 35.671358392019094e0);
	report("gsw_sstar_from_sp", sstar_from_sp_ca,
		gsw_sstar_from_sp(sa,p,lon,lat), 35.866946753006239e0);
	report("gsw_ct_from_t", ct_from_t_ca,
		gsw_ct_from_t(sa,t,p),  14.930280459895560e0);
	printf(
"\nOther conversions between temperatures, salinities, entropy, pressure\n"
"and height:\n\n"
	);
	report("gsw_deltasa_from_sp", deltasa_from_sp_ca,
		gsw_deltasa_from_sp(sp,p,lon,lat),  3.96067773336028495e-3);
	report("gsw_sr_from_sp", sr_from_sp_ca,
		gsw_sr_from_sp(sp),  35.667397714285734e0);
	report("gsw_sp_from_sr", sp_from_sr_ca,
		gsw_sp_from_sr(sr),  35.333387933015295e0);
	report("gsw_sp_from_sa", sp_from_sa_ca,
		gsw_sp_from_sa(sa,p,lon,lat), 35.528504019167094e0);
	report("gsw_sstar_from_sa", sstar_from_sa_ca,
		gsw_sstar_from_sa(sa,p,lon,lat), 35.694648791860907e0);
	report("gsw_sp_from_sstar", sp_from_sstar_ca,
		gsw_sp_from_sstar(sstar,p,lon,lat), 35.334761242083573e0);
	report("gsw_sa_from_sstar", sa_from_sstar_ca,
		gsw_sa_from_sstar(sstar,p,lon,lat), 35.505322027120805e0);
	report("gsw_pt_from_ct", pt_from_ct_ca,
		gsw_pt_from_ct(sa,ct), 20.023899375975017e0);
	report("gsw_t_from_ct", t_from_ct_ca,
		gsw_t_from_ct(sa,ct,p), 20.079820359223014e0);
	report("gsw_ct_from_pt", ct_from_pt_ca,
		gsw_ct_from_pt(sa,pt), 14.976021403957613e0);
	report("gsw_pt0_from_t", pt0_from_t_ca,
		gsw_pt0_from_t(sa,t,p),  14.954241363902305e0);
	report("gsw_pt_from_t", pt_from_t_ca,
		gsw_pt_from_t(sa,t,p,p_ref), 14.969381237883740e0);
	report("gsw_z_from_p", z_from_p_ca,
		gsw_z_from_p(p,lat), -2.980161553316402e2);
	report("gsw_entropy_from_t", entropy_from_t_ca,
		gsw_entropy_from_t(sa,t,p), 212.30166821093002e0);
	report("gsw_adiabatic_lapse_rate_from_ct",
		adiabatic_lapse_rate_from_ct_ca,
		gsw_adiabatic_lapse_rate_from_ct(sa,ct,p),
		1.877941744191155e-8);
	printf(
"\nDensity and enthalpy, based on the 48-term expression for density:\n\n"
	);
	report("gsw_rho", rho_ca,
		gsw_rho(sa,ct,p), 1026.4562376198473e0);
	report("gsw_alpha", alpha_ca,
		gsw_alpha(sa,ct,p), 2.62460550806784356e-4);
	report("gsw_beta", beta_ca,
		gsw_beta(sa,ct,p), 7.29314455934463365e-4 );
	report("gsw_alpha_on_beta", alpha_on_beta_ca,
		gsw_alpha_on_beta(sa,ct,p), 0.359872958325632e0 );
	gsw_rho_first_derivatives(sa,ct,p,&drho_dsa,&drho_dct,&drho_dp);
	drho_dsa_error = fabs(drho_dsa - 0.748609372480258e0);
	drho_dct_error = fabs(drho_dct + 0.269404269504765e0);
	drho_dp_error = fabs(drho_dp - 4.287533235942749e-7);
	printf("\t%-32.32s  ........  ", "gsw_rho_first_derivatives");
	if (drho_dsa_error < drho_dsa_ca && drho_dct_error < drho_dct_ca
	    && drho_dp_error < drho_dp_ca)
            printf("passed\n");
        else {
            printf("failed\n");
            gsw_error_flag      = 1;
 	}
	report("gsw_specvol", specvol_ca,
		gsw_specvol(sa,ct,p), 9.74225654586897711e-4 );
	report("gsw_specvol_anom", specvol_anom_ca,
		gsw_specvol_anom(sa,ct,p), 2.90948181201264571e-6 );
	report("gsw_sigma0", sigma0_ca,
		gsw_sigma0(sa,ct), 25.165674636323047e0);
	report("gsw_sigma1", sigma1_ca,
		gsw_sigma1(sa,ct), 29.434338510752923e0);
	report("gsw_sigma2", sigma2_ca,
		gsw_sigma2(sa,ct), 33.609842926904093e0);
	report("gsw_sigma3", sigma3_ca,
		gsw_sigma3(sa,ct), 37.695147569371784e0);
	report("gsw_sigma4", sigma4_ca,
		gsw_sigma4(sa,ct), 41.693064726656303e0);
	report("gsw_sound_speed", sound_speed_ca,
		gsw_sound_speed(sa,ct,p), 1527.2011773569989e0 );
	report("gsw_kappa", kappa_ca,
		gsw_kappa(sa,ct,p), 4.177024873349404e-010);
	report("gsw_cabbeling", cabbeling_ca,
		gsw_cabbeling(sa,ct,p), 9.463053321129075e-6);
	report("gsw_thermobaric", thermobaric_ca,
		gsw_thermobaric(sa,ct,p), 1.739078662082863e-12);
	report("gsw_internal_energy", internal_energy_ca,
		gsw_internal_energy(sa,ct,p), 79740.482561720783e0 );
	report("gsw_enthalpy", enthalpy_ca,
		gsw_enthalpy(sa,ct,p), 82761.872939932495e0 );
	report("gsw_dynamic_enthalpy", dynamic_enthalpy_ca,
		gsw_dynamic_enthalpy(sa,ct,p),  2924.5137975399025e0 );
	rho	= gsw_rho(sa,ct,p);
	report("gsw_sa_from_rho", rho_ca,
		gsw_sa_from_rho(rho,ct,p), sa);
	printf(
  "\nWater column properties, based on the 48-term expression for density:\n\n"
	);
	gsw_nsquared(sa_profile,ct_profile,p_profile,lat_profile,nz,n2,
			p_mid_n2);
	n2_error[0] = fabs(n2[0] + 0.070960392693051e-3);
	n2_error[1] = fabs(n2[1] - 0.175435821615983e-3);
	p_mid_n2_error[0] = fabs(p_mid_n2[0] - 25e0);
	p_mid_n2_error[1] = fabs(p_mid_n2[1] - 75e0);
	printf("\t%-32.32s  ........  ", "gsw_nsquared");
	if (n2_error[0] < n2_ca && n2_error[1] < n2_ca &&
	    p_mid_n2_error[0] < p_mid_n2_ca && p_mid_n2_error[1] < p_mid_n2_ca)
            printf("passed\n");
        else {
            printf("failed\n");
            gsw_error_flag      = 1;
 	}
	gsw_turner_rsubrho(sa_profile,ct_profile,p_profile,nz,tu,
				rsubrho,p_mid_tursr);
	tu_error[0] = fabs(tu[0] + 1.187243981606485e2);
	tu_error[1] = fabs(tu[1] - 0.494158257088517e2);
	rsubrho_error[0] = fabs(rsubrho[0] - 3.425146897090065e0);
	rsubrho_error[1] = fabs(rsubrho[1] - 12.949399443139164e0);
	p_mid_tursr_error[0] = fabs(p_mid_tursr[0] - 25e0);
	p_mid_tursr_error[1] = fabs(p_mid_tursr[1] - 75e0);
	printf("\t%-32.32s  ........  ", "gsw_turner_rsubrho");
	if (tu_error[0] < tu_ca && tu_error[1] < tu_ca &&  
	    rsubrho_error[0] < rsubrho_ca && rsubrho_error[1] < rsubrho_ca &&  
	    p_mid_tursr_error[0] < p_mid_tursr_ca &&
	    p_mid_tursr_error[1] < p_mid_tursr_ca)
            printf("passed\n");
        else {
            printf("failed\n");
            gsw_error_flag      = 1;
 	}
	gsw_ipv_vs_fnsquared_ratio(sa_profile,ct_profile,p_profile,nz,
					ipvfn2,p_mid_ipvfn2);
	ipvfn2_error[0] = fabs(ipvfn2[0] - 0.996783975249010e0);
	ipvfn2_error[1] = fabs(ipvfn2[1] - 0.992112251478320e0);
	p_mid_ipvfn2_error[0] = fabs(p_mid_ipvfn2[0] - 25e0);
	p_mid_ipvfn2_error[1] = fabs(p_mid_ipvfn2[1] - 75e0);
	printf("\t%-32.32s  ........  ", "gsw_ipv_vs_fnsquared_ratio");
	if (ipvfn2_error[0] < ipvfn2_ca && 
	    ipvfn2_error[1] < ipvfn2_ca && 
	    p_mid_ipvfn2_error[0] < p_mid_ipvfn2_ca && 
	    p_mid_ipvfn2_error[1] < p_mid_ipvfn2_ca)
            printf("passed\n");
        else {
            printf("failed\n");
            gsw_error_flag      = 1;
 	}
	printf(
	"\nFreezing temperatures:\n\n"
	);
	report("gsw_ct_freezing", ct_freezing_ca,
		gsw_ct_freezing(sa,p,saturation_fraction), -2.1801450326174852e0);
	report("gsw_t_freezing", t_freezing_ca,
		gsw_t_freezing(sa,p,saturation_fraction), -2.1765521998023516e0);
	printf(
"\nIsobaric melting enthalpy and isobaric evaporation enthalpy:\n\n"
	);
	report("gsw_latentheat_melting", latentheat_melting_ca,
		gsw_latentheat_melting(sa,p), 329330.54839618353e0);
	report("gsw_latentheat_evap_ct", latentheat_evap_ct_ca,
		gsw_latentheat_evap_ct(sa,ct), 2450871.0228523901e0);
	report("gsw_latentheat_evap_t", latentheat_evap_t_ca,
		gsw_latentheat_evap_t(sa,t), 2462848.2895522709e0);
	printf(
"\nBasic thermodynamic properties in terms of in-situ t, based on\n"
"the exact Gibbs function:\n\n"
	);
	report("gsw_rho_t_exact", rho_t_exact_ca,
		gsw_rho_t_exact(sa,t,p), 1027.7128170207150e0);
	report("gsw_pot_rho_t_exact", pot_rho_t_exact_ca,
		gsw_pot_rho_t_exact(sa,t,p,p_ref), 1026.8362655887486e0);
	report("gsw_alpha_wrt_t_exact", alpha_wrt_t_exact_ca,
		gsw_alpha_wrt_t_exact(sa,t,p), 2.19066952410728916e-4);
	report("gsw_beta_const_t_exact", beta_const_t_exact_ca,
		gsw_beta_const_t_exact(sa,t,p),  7.44744841648729426e-4);
	report("gsw_specvol_t_exact", specvol_t_exact_ca,
		gsw_specvol_t_exact(sa,t,p), 9.73034473676164815e-4);
	report("gsw_sound_speed_t_exact", sound_speed_t_exact_ca,
		gsw_sound_speed_t_exact(sa,t,p), 1512.2053940303056e0);
	report("gsw_kappa_t_exact", kappa_t_exact_ca,
		gsw_kappa_t_exact(sa,t,p), 4.25506953386609075e-010);
	report("gsw_enthalpy_t_exact", enthalpy_t_exact_ca,
		gsw_enthalpy_t_exact(sa,t,p), 62520.680485510929e0);
	report("gsw_cp_t_exact", cp_t_exact_ca,
		gsw_cp_t_exact(sa,t,p), 3982.7832563441461e0);
	printf(
"\nLibrary functions of the GSW toolbox:\n\n"
	);
	report("gsw_deltasa_atlas", deltasa_atlas_ca,
		gsw_deltasa_atlas(p,lon,lat), 3.87660373016291727e-3);
	report("gsw_fdelta", fdelta_ca,
		gsw_fdelta(p,lon,lat), 1.49916256924158942e-004);
	report("gsw_sa_from_sp_baltic", sa_from_sp_ca,
		gsw_sa_from_sp_baltic(sp,long_bs,lat_bs) , 35.666154857142850e0);
	report("gsw_sp_from_sa_baltic", sp_from_sa_ca,
		gsw_sp_from_sa_baltic(sa,long_bs,lat_bs), 35.533769845749660e0);

	if (gsw_error_flag)
	    printf("\nYour installation of the Gibbs SeaWater (GSW) "
		"Oceanographic Toolbox has errors !\n");
	else
	    printf("\nWell done! The gsw_check_fuctions confirms that the\n"
		"Gibbs SeaWater (GSW) Oceanographic Toolbox is "
		"installed correctly.\n");

	return (0);
}

static void
report(char *what, double acceptable, double actual, double expected)
{
	double	diff=fabs(actual-expected);

	printf("\t%-32.32s  ........  ", what);
	if (diff < acceptable)
	    printf("passed\n");
	else {
	    printf("failed\n");
	    gsw_error_flag	= 1;
	    if (debug)
		printf("actual=%g\texpected=%g\tdiff=%g\n",actual,expected,
			diff);
	}
}

/*
**  The End.
*/
