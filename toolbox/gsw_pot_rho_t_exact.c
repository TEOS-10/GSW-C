/*
!==========================================================================
function gsw_pot_rho_t_exact(sa,t,p,p_ref)  
!==========================================================================

! Calculates the potential density of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! p_ref  : reference sea pressure                          [dbar]
! 
! gsw_pot_rho_t_exact : potential density                  [kg/m^3]
*/
double
gsw_pot_rho_t_exact(double sa, double t, double p, double p_ref)
{
	double	pt = gsw_pt_from_t(sa,t,p,p_ref);

	return (gsw_rho_t_exact(sa,pt,p_ref));
}
