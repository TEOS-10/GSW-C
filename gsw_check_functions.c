/*
**  $Id: gsw_check_functions.c,v e3398b4e1644 2011/09/23 21:39:29 fdelahoyde $
*/
#include <teos-10.h>

int
main(int argc, char **argv)
{
double sp, sa, sstar, sr, t, ct, pt, p, p_bs, p_ref ;
double lon, long_bs, lat, lat_bs, saturation_fraction;

sp =  35.5e0;
sa = 35.7e0;
sstar = 35.5e0;
sr = 35.5e0;
t = 15e0;
ct = 20e0;
pt = 15e0;
p = 300e0;
p_bs = 50e0;
p_ref = 100e0;
lon = 260e0;
long_bs = 20e0;
lat = 16e0;
lat_bs = 60e0;
saturation_fraction = 0.5e0;

printf("============================================================================\n");
printf(" Gibbs SeaWater (GSW) Oceanographic Toolbox of TEOS-10 version 3.0 (C)\n");
printf("============================================================================\n");
printf("\n");
printf(" These are the check values for the subset of functions that have been \n");
printf(" converted into C from the Gibbs SeaWater (GSW) Oceanographic Toolbox \n");
printf(" of TEOS-10 (version 3.0).\n");
printf("\n");
printf("Absolute Salinity, Preformed Salinity and Conservative Temperature\n");
printf("gsw_sa_from_sp          : %.18g %.18g\n", gsw_sa_from_sp(sp,p,lon,lat), 35.671358392019094e0);
printf("gsw_sstar_from_sp       : %.18g %.18g\n", gsw_sstar_from_sp(sa,p,lon,lat), 35.866946753006239e0);
printf("gsw_ct_from_t           : %.18g %.18g\n", gsw_ct_from_t(sa,t,p),  14.930280459895560e0);
printf(" \n");
printf("other conversions between temperatures, salinities, entropy, pressure and height\n");
printf("gsw_deltasa_from_sp     : %.18g %.18g\n", gsw_deltasa_from_sp(sp,p,lon,lat),  3.96067773336028495e-3);
printf("gsw_sr_from_sp          : %.18g %.18g\n", gsw_sr_from_sp(sp),  35.667397714285734e0);
printf("gsw_sp_from_sr          : %.18g %.18g\n", gsw_sp_from_sr(sr),  35.333387933015295e0);
printf("gsw_sp_from_sa          : %.18g %.18g\n", gsw_sp_from_sa(sa,p,lon,lat), 35.528504019167094e0);
printf("gsw_sstar_from_sa       : %.18g %.18g\n", gsw_sstar_from_sa(sa,p,lon,lat), 35.694648791860907e0);
printf("gsw_sp_from_sstar       : %.18g %.18g\n", gsw_sp_from_sstar(sstar,p,lon,lat), 35.334761242083573e0);
printf("gsw_sa_from_sstar       : %.18g %.18g\n", gsw_sa_from_sstar(sstar,p,lon,lat), 35.505322027120805e0);
printf("gsw_pt_from_ct          : %.18g %.18g\n", gsw_pt_from_ct(sa,ct), 20.023899375975017e0);
printf("gsw_t_from_ct           : %.18g %.18g\n", gsw_t_from_ct(sa,ct,p), 20.079820359223014e0);
printf("gsw_ct_from_pt          : %.18g %.18g\n", gsw_ct_from_pt(sa,pt), 14.976021403957613e0);
printf("gsw_pt0_from_t          : %.18g %.18g\n", gsw_pt0_from_t(sa,t,p),  14.954241363902305e0);
printf("gsw_pt_from_t           : %.18g %.18g\n", gsw_pt_from_t(sa,t,p,p_ref), 14.969381237883740e0);
printf(" \n");
printf("density and enthalpy, based on the 48-term expression for density\n");
printf("gsw_rho                 : %.18g %.18g\n", gsw_rho(sa,ct,p), 1026.4562376198473e0);
printf("gsw_alpha               : %.18g %.18g\n", gsw_alpha(sa,ct,p), 2.62460550806784356e-4);
printf("gsw_beta                : %.18g %.18g\n", gsw_beta(sa,ct,p), 7.29314455934463365e-4 );
printf("gsw_specvol             : %.18g %.18g\n", gsw_specvol(sa,ct,p), 9.74225654586897711e-4 );
printf("gsw_specvol_anom        : %.18g %.18g\n", gsw_specvol_anom(sa,ct,p), 2.90948181201264571e-6 );
printf("gsw_sound_speed         : %.18g %.18g\n", gsw_sound_speed(sa,ct,p), 1527.2011773569989e0 );
printf("gsw_internal_energy     : %.18g %.18g\n", gsw_internal_energy(sa,ct,p), 79740.482561720783e0 );
printf("gsw_enthalpy            : %.18g %.18g\n", gsw_enthalpy(sa,ct,p), 82761.872939932495e0 );
printf("gsw_dynamic_enthalpy    : %.18g %.18g\n", gsw_dynamic_enthalpy(sa,ct,p),  2924.5137975399025e0 );
printf(" \n");
printf("freezing temperatures\n");
printf("gsw_ct_freezing         : %.18g %.18g\n", gsw_ct_freezing(sa,p,saturation_fraction), -2.1801450326174852e0);
printf("gsw_t_freezing          : %.18g %.18g\n", gsw_t_freezing(sa,p,saturation_fraction), -2.1765521998023516e0);
printf(" \n");
printf("isobaric melting enthalpy and isobaric evaporation enthalpy\n");
printf("gsw_latentheat_melting  : %.18g %.18g\n", gsw_latentheat_melting(sa,p), 329330.54839618353e0);
printf("gsw_latentheat_evap_ct  : %.18g %.18g\n", gsw_latentheat_evap_ct(sa,ct), 2450871.0228523901e0);
printf("gsw_latentheat_evap_t   : %.18g %.18g\n", gsw_latentheat_evap_t(sa,t), 2462848.2895522709e0);
printf(" \n");
printf("basic thermodynamic properties in terms of in-situ t, based on the exact Gibbs function\n");
printf("gsw_rho_t_exact         : %.18g %.18g\n", gsw_rho_t_exact(sa,t,p), 1027.7128170207150e0);
printf("gsw_pot_rho_t_exact     : %.18g %.18g\n", gsw_pot_rho_t_exact(sa,t,p,p_ref), 1026.8362655887486e0);
printf("gsw_alpha_wrt_t_exact   : %.18g %.18g\n", gsw_alpha_wrt_t_exact(sa,t,p), 2.19066952410728916e-4);
printf("gsw_beta_const_t_exact  : %.18g %.18g\n", gsw_beta_const_t_exact(sa,t,p),  7.44744841648729426e-4);
printf("gsw_specvol_t_exact     : %.18g %.18g\n", gsw_specvol_t_exact(sa,t,p), 9.73034473676164815e-4);
printf("gsw_sound_speed_t_exact : %.18g %.18g\n", gsw_sound_speed_t_exact(sa,t,p), 1512.2053940303056e0);
printf("gsw_kappa_t_exact       : %.18g %.18g\n", gsw_kappa_t_exact(sa,t,p), 4.25506953386609075e-010);
printf("gsw_enthalpy_t_exact    : %.18g %.18g\n", gsw_enthalpy_t_exact(sa,t,p), 62520.680485510929e0);
printf("gsw_entropy_t_exact     : %.18g %.18g\n", gsw_entropy_t_exact(sa,t,p), 212.30166821093002e0);
printf("gsw_cp_t_exact          : %.18g %.18g\n", gsw_cp_t_exact(sa,t,p), 3982.7832563441461e0);
printf(" \n");
printf("library functions of the GSW toolbox\n");
printf("gsw_delta_sa_ref        : %.18g %.18g\n", gsw_delta_sa_ref(p,lon,lat), 3.87660373016291727e-3);
printf("gsw_fdelta              : %.18g %.18g\n", gsw_fdelta(p,lon,lat), 1.49916256924158942e-004);
printf("gsw_sa_from_sp_baltic   : %.18g %.18g\n", gsw_sa_from_sp_baltic(sp,long_bs,lat_bs) , 35.666154857142850e0);
printf("gsw_sp_from_sa_baltic   : %.18g %.18g\n", gsw_sp_from_sa_baltic(sa,long_bs,lat_bs), 35.533769845749660e0);
return (0);
}

           
