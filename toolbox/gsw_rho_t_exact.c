/*
!==========================================================================
function gsw_rho_t_exact(sa,t,p)  
!==========================================================================

! Calculates in-situ density of seawater from Absolute Salinity and 
! in-situ temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_rho_t_exact : in-situ density                        [kg/m^3]
*/
double
gsw_rho_t_exact(double sa, double t, double p)
{
	int	n0=0, n1=1;

	return (1.0/gsw_gibbs(n0,n0,n1,sa,t,p));
}
