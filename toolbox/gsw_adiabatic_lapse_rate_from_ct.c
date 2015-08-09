/*
!==========================================================================
function gsw_adiabatic_lapse_rate_from_ct(sa,ct,p)
!==========================================================================

! Calculates the adiabatic lapse rate from Conservative Temperature
!
! sa     : Absolute Salinity                                 [g/kg]
! ct     : Conservative Temperature                          [deg C]
! p      : sea pressure                                      [dbar]
!
! gsw_adiabatic_lapse_rate_from_ct : adiabatic lapse rate    [K/Pa]
*/
double
gsw_adiabatic_lapse_rate_from_ct(double sa, double ct, double p)
{
	int	n0=0, n1=1, n2=2;
	double	pt0, pr0=0.0, t;

	pt0	= gsw_pt_from_ct(sa,ct);
	t	= gsw_pt_from_t(sa,pt0,pr0,p);

	return (-gsw_gibbs(n0,n1,n1,sa,t,p)/gsw_gibbs(n0,n2,n0,sa,t,p));

}
