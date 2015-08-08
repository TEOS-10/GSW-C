/*
!==========================================================================
function gsw_sigma3(sa,ct) 
!==========================================================================

!  Calculates potential density anomaly with reference pressure of 3000 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  function has inputs of Absolute Salinity and Conservative Temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! sigma3 : potential density anomaly with reference pressure of 3000
*/
double
gsw_sigma3(double sa, double ct)
{
	return (gsw_rho(sa,ct,3000.0) - 1000.0);
}
