/*
!==========================================================================
function gsw_beta_const_t_exact(sa,t,p)  
!==========================================================================

! Calculates saline (haline) contraction coefficient of seawater at 
! constant in-situ temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! beta_const_t_exact : haline contraction coefficient      [kg/g]
*/
double
gsw_beta_const_t_exact(double sa, double t, double p)
{
	int	n0=0, n1=1;

	return (-gsw_gibbs(n1,n0,n1,sa,t,p)/gsw_gibbs(n0,n0,n1,sa,t,p));
}
