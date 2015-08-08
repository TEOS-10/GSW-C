/*
!==========================================================================
function gsw_t_from_ct(sa,ct,p)  
!==========================================================================

! Calculates in-situ temperature from Conservative Temperature of seawater  
!
! sa      : Absolute Salinity                              [g/kg]
! ct      : Conservative Temperature                       [deg C]
!
! gsw_t_from_ct : in-situ temperature                      [deg C]
*/
double
gsw_t_from_ct(double sa, double ct, double p)
{
	double	pt0, p0=0.0;

	pt0	= gsw_pt_from_ct(sa,ct);
	return (gsw_pt_from_t(sa,pt0,p0,p));
}
