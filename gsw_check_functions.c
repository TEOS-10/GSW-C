/*
**  $Id: gsw_check_functions.c,v 0fa6ed68e79e 2011/09/25 18:18:19 fdelahoyde $
*/
#include <gswteos-10.h>

static void	report(char *what, double actual, double expected);

int
main(int argc, char **argv)
{
	double	sp, sa, sstar, sr, t, ct, pt, p, p_bs, p_ref ;
	double	lon, long_bs, lat, lat_bs, saturation_fraction;

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
	lon	= 260e0;
	long_bs	= 20e0;
	lat	= 16e0;
	lat_bs	= 60e0;
	saturation_fraction	= 0.5e0;

	printf(
"============================================================================\n"
" Gibbs SeaWater (GSW) Oceanographic Toolbox of TEOS-10 version 3.0 (C)\n"
"============================================================================\n"
"\n"
" These are the check values for the subset of functions that have been \n"
" converted into C from the Gibbs SeaWater (GSW) Oceanographic Toolbox \n"
" of TEOS-10 (version 3.0).\n\n"
	);
	printf("\t%-24.24s  %18.18s %18.18s\n", "", "Actual", "Expected");
	printf(
"\nAbsolute Salinity, Preformed Salinity and Conservative Temperature:\n\n"
	);
	report("gsw_sa_from_sp",
		gsw_sa_from_sp(sp,p,lon,lat), 35.671358392019094e0);
	report("gsw_sstar_from_sp",
		gsw_sstar_from_sp(sa,p,lon,lat), 35.866946753006239e0);
	report("gsw_ct_from_t",
		gsw_ct_from_t(sa,t,p),  14.930280459895560e0);
	printf(
"\nOther conversions between temperatures, salinities, entropy, pressure\n"
"and height:\n\n"
	);
	report("gsw_deltasa_from_sp",
		gsw_deltasa_from_sp(sp,p,lon,lat),  3.96067773336028495e-3);
	report("gsw_sr_from_sp",
		gsw_sr_from_sp(sp),  35.667397714285734e0);
	report("gsw_sp_from_sr",
		gsw_sp_from_sr(sr),  35.333387933015295e0);
	report("gsw_sp_from_sa",
		gsw_sp_from_sa(sa,p,lon,lat), 35.528504019167094e0);
	report("gsw_sstar_from_sa",
		gsw_sstar_from_sa(sa,p,lon,lat), 35.694648791860907e0);
	report("gsw_sp_from_sstar",
		gsw_sp_from_sstar(sstar,p,lon,lat), 35.334761242083573e0);
	report("gsw_sa_from_sstar",
		gsw_sa_from_sstar(sstar,p,lon,lat), 35.505322027120805e0);
	report("gsw_pt_from_ct",
		gsw_pt_from_ct(sa,ct), 20.023899375975017e0);
	report("gsw_t_from_ct",
		gsw_t_from_ct(sa,ct,p), 20.079820359223014e0);
	report("gsw_ct_from_pt",
		gsw_ct_from_pt(sa,pt), 14.976021403957613e0);
	report("gsw_pt0_from_t",
		gsw_pt0_from_t(sa,t,p),  14.954241363902305e0);
	report("gsw_pt_from_t",
		gsw_pt_from_t(sa,t,p,p_ref), 14.969381237883740e0);
	printf(
"\nDensity and enthalpy, based on the 48-term expression for density:\n\n"
	);
	report("gsw_rho",
		gsw_rho(sa,ct,p), 1026.4562376198473e0);
	report("gsw_alpha",
		gsw_alpha(sa,ct,p), 2.62460550806784356e-4);
	report("gsw_beta",
		gsw_beta(sa,ct,p), 7.29314455934463365e-4 );
	report("gsw_specvol",
		gsw_specvol(sa,ct,p), 9.74225654586897711e-4 );
	report("gsw_specvol_anom",
		gsw_specvol_anom(sa,ct,p), 2.90948181201264571e-6 );
	report("gsw_sound_speed",
		gsw_sound_speed(sa,ct,p), 1527.2011773569989e0 );
	report("gsw_internal_energy",
		gsw_internal_energy(sa,ct,p), 79740.482561720783e0 );
	report("gsw_enthalpy",
		gsw_enthalpy(sa,ct,p), 82761.872939932495e0 );
	report("gsw_dynamic_enthalpy",
		gsw_dynamic_enthalpy(sa,ct,p),  2924.5137975399025e0 );
	printf(
	"\nFreezing temperatures:\n\n"
	);
	report("gsw_ct_freezing",
		gsw_ct_freezing(sa,p,saturation_fraction), -2.1801450326174852e0);
	report("gsw_t_freezing",
		gsw_t_freezing(sa,p,saturation_fraction), -2.1765521998023516e0);
	printf(
"\nIsobaric melting enthalpy and isobaric evaporation enthalpy:\n\n"
	);
	report("gsw_latentheat_melting",
		gsw_latentheat_melting(sa,p), 329330.54839618353e0);
	report("gsw_latentheat_evap_ct",
		gsw_latentheat_evap_ct(sa,ct), 2450871.0228523901e0);
	report("gsw_latentheat_evap_t",
		gsw_latentheat_evap_t(sa,t), 2462848.2895522709e0);
	printf(
"\nBasic thermodynamic properties in terms of in-situ t, based on\n"
"the exact Gibbs function:\n\n"
	);
	report("gsw_rho_t_exact",
		gsw_rho_t_exact(sa,t,p), 1027.7128170207150e0);
	report("gsw_pot_rho_t_exact",
		gsw_pot_rho_t_exact(sa,t,p,p_ref), 1026.8362655887486e0);
	report("gsw_alpha_wrt_t_exact",
		gsw_alpha_wrt_t_exact(sa,t,p), 2.19066952410728916e-4);
	report("gsw_beta_const_t_exact",
		gsw_beta_const_t_exact(sa,t,p),  7.44744841648729426e-4);
	report("gsw_specvol_t_exact",
		gsw_specvol_t_exact(sa,t,p), 9.73034473676164815e-4);
	report("gsw_sound_speed_t_exact",
		gsw_sound_speed_t_exact(sa,t,p), 1512.2053940303056e0);
	report("gsw_kappa_t_exact",
		gsw_kappa_t_exact(sa,t,p), 4.25506953386609075e-010);
	report("gsw_enthalpy_t_exact",
		gsw_enthalpy_t_exact(sa,t,p), 62520.680485510929e0);
	report("gsw_entropy_t_exact",
		gsw_entropy_t_exact(sa,t,p), 212.30166821093002e0);
	report("gsw_cp_t_exact",
		gsw_cp_t_exact(sa,t,p), 3982.7832563441461e0);
	printf(
"\nLibrary functions of the GSW toolbox:\n\n"
	);
	report("gsw_delta_sa_ref",
		gsw_delta_sa_ref(p,lon,lat), 3.87660373016291727e-3);
	report("gsw_fdelta",
		gsw_fdelta(p,lon,lat), 1.49916256924158942e-004);
	report("gsw_sa_from_sp_baltic",
		gsw_sa_from_sp_baltic(sp,long_bs,lat_bs) , 35.666154857142850e0);
	report("gsw_sp_from_sa_baltic",
		gsw_sp_from_sa_baltic(sa,long_bs,lat_bs), 35.533769845749660e0);
	return (0);
}

static void
report(char *what, double actual, double expected)
{
	double	diff=actual-expected;
	char	stat[64];

	if (fabs(diff)>1e-14)
	    sprintf(stat,"***\n\t\t\tDiff=%g",diff);
	else
	    sprintf(stat,"OK");
	printf("\t%-24.24s: %.18g %.18g %s\n", what, actual, expected, stat);
}

/*
**  The End.
*/
