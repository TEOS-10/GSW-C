/*
!==========================================================================
function gsw_t_freezing(sa,p,saturation_fraction)
!==========================================================================

! Calculates the in-situ temperature at which seawater freezes
!
! sa     : Absolute Salinity                                 [g/kg]
! p      : sea pressure                                      [dbar]
!         ( i.e. absolute pressure - 10.1325 dbar )
! saturation_fraction : the saturation fraction of dissolved air
!                       in seawater
!
! t_freezing : in-situ temperature at which seawater freezes.[deg C]
*/
double
gsw_t_freezing(double sa, double p, double saturation_fraction)
{
	return (gsw_t_freezing_exact(sa,p,saturation_fraction));
}
