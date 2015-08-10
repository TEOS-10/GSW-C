/*
!==========================================================================
function gsw_cp_t_exact(sa,t,p)
!==========================================================================

! Calculates isobaric heat capacity of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_cp_t_exact : heat capacity                           [J/(kg K)]
*/
double
gsw_cp_t_exact(double sa, double t, double p)
{
	int	n0, n2;

	n0 = 0;
	n2 = 2;

	return (-(t+273.15e0)*gsw_gibbs(n0,n2,n0,sa,t,p));
}
