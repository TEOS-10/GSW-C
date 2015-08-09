/*
!==========================================================================
function gsw_latentheat_evap_t(sa,t)  
!==========================================================================
!
! Calculates latent heat, or enthalpy, of evaporation.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! 
! gsw_latentheat_evap_t : latent heat of evaporation       [J/kg]
*/
double
gsw_latentheat_evap_t(double sa, double t)
{

	double	ct = gsw_ct_from_pt(sa,t);

	return (gsw_latentheat_evap_ct(sa,ct));
}
