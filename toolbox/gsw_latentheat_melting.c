/*
!--------------------------------------------------------------------------
! isobaric melting enthalpy and isobaric evaporation enthalpy
!--------------------------------------------------------------------------

!==========================================================================
function gsw_latentheat_melting(sa,p)  
!==========================================================================

! Calculates latent heat, or enthalpy, of melting.
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! 
! latentheat_melting : latent heat of melting              [kg/m^3]
*/
double
gsw_latentheat_melting(double sa, double p)
{
	GSW_TEOS10_CONSTANTS;
	double	tf = gsw_t_freezing(sa,p,0.0);

	return (1000.0*(gsw_chem_potential_water_t_exact(sa,tf,p)
           - (gsw_t0 + tf)*gsw_t_deriv_chem_potential_water_t_exact(sa,tf,p))
           - gsw_enthalpy_ice(tf,p));
}
