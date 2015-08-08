/*
!==========================================================================
function gsw_sigma4(sa,ct) 
!==========================================================================

!  Calculates potential density anomaly with reference pressure of 4000 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  function has inputs of Absolute Salinity and Conservative Temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! sigma4  : potential density anomaly with reference pressure of 4000
*/
double
gsw_sigma4(double sa, double ct)
{
	return (gsw_rho(sa,ct,4000.0) - 1000.0);
}
