/*
!==========================================================================
function gsw_specvol_t_exact(sa,t,p)  
!==========================================================================

! Calculates the specific volume of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! specvol_t_exact : specific volume                        [kg/m^3]
*/
double
gsw_specvol_t_exact(double sa, double t, double p)
{
	int	n0=0, n1=1;

	return (gsw_gibbs(n0,n0,n1,sa,t,p));
}
