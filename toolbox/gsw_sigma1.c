/*
!==========================================================================
function gsw_sigma1(sa,ct) 
!==========================================================================

!  Calculates potential density anomaly with reference pressure of 1000 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  function has inputs of Absolute Salinity and Conservative Temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! sigma1 : potential density anomaly with reference pressure of 1000
*/
double
gsw_sigma1(double sa, double ct)
{
	return (gsw_rho(sa,ct,1000.0) - 1000.0);
}
