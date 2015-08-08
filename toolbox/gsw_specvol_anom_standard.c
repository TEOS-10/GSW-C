/*
!==========================================================================
function gsw_specvol_anom_standard(sa,ct,p)
!==========================================================================
!
!  Calculates specific volume anomaly of seawater.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
! 
! specvol_anom  :  specific volume anomaly of seawater
*/
double
gsw_specvol_anom_standard(double sa, double ct, double p)
{
	return (gsw_specvol(sa,ct,p) - gsw_specvol_sso_0(p));
}
