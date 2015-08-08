/*
!==========================================================================
function gsw_alpha_wrt_t_exact(sa,t,p)  
!==========================================================================

! Calculates thermal expansion coefficient of seawater with respect to 
! in-situ temperature
!
! sa     : Absolute Salinity                               [g/kg]
! t      : insitu temperature                              [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_alpha_wrt_t_exact : thermal expansion coefficient    [1/K]
!                         wrt (in-situ) temperature
*/
double
gsw_alpha_wrt_t_exact(double sa, double t, double p)
{
	int	n0=0, n1=1;

	return (gsw_gibbs(n0,n1,n1,sa,t,p)/gsw_gibbs(n0,n0,n1,sa,t,p));
}
