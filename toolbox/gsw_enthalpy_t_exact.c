/*
!==========================================================================
function gsw_enthalpy_t_exact(sa,t,p)  
!==========================================================================

! Calculates the specific enthalpy of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_enthalpy_t_exact : specific enthalpy                 [J/kg]
*/
double
gsw_enthalpy_t_exact(double sa, double t, double p)
{
	GSW_TEOS10_CONSTANTS;
	int	n0=0, n1=1;

	return (gsw_gibbs(n0,n0,n0,sa,t,p) -
		(t+gsw_t0)*gsw_gibbs(n0,n1,n0,sa,t,p));
}
