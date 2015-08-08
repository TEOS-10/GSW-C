/*
!==========================================================================
function gsw_entropy_from_t(sa,t,p)
!==========================================================================

! Calculates the specific entropy of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_entropy_from_t : specific entropy                    [J/(kg K)]
*/
double
gsw_entropy_from_t(double sa, double t, double p)
{
	int	n0=0, n1=1;

	return (-gsw_gibbs(n0,n1,n0,sa,t,p));

}
