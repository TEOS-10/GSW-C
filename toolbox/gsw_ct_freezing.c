/*
!--------------------------------------------------------------------------
! freezing temperatures
!--------------------------------------------------------------------------

!==========================================================================
function gsw_ct_freezing(sa,p,saturation_fraction)
!==========================================================================

!  Calculates the Conservative Temperature at which seawater freezes.
!  The error of this fit ranges between -5e-4 K and 6e-4 K when compared
!  with the Conservative Temperature calculated from the exact in-situ
!  freezing temperature which is found by a Newton-Raphson iteration of the
!  equality of the chemical potentials of water in seawater and in ice.
!  Note that the Conservative temperature freezing temperature can be found
!  by this exact method using the function gsw_ct_freezing_exact.
!
! sa     : Absolute Salinity                                 [g/kg]
! p      : sea pressure                                      [dbar]
! saturation_fraction : saturation fraction of dissolved air in
|                       seawater
!
! ct_freezing : Conservative Temperature at freezing point   [deg C]
!                That is, the freezing temperature expressed in
!                terms of Conservative Temperature (ITS-90).
*/
double
gsw_ct_freezing(double sa, double p, double saturation_fraction)
{
	return (gsw_ct_freezing_exact(sa, p, saturation_fraction));
}
