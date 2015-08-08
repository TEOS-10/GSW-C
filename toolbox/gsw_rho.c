/*
!--------------------------------------------------------------------------
! density and enthalpy, based on the 48-term expression for density
!--------------------------------------------------------------------------

!==========================================================================
function gsw_rho(sa,ct,p)  
!==========================================================================

!  Calculates in-situ density from Absolute Salinity and Conservative 
!  Temperature, using the computationally-efficient expression for
!  specific volume in terms of SA, CT and p (Roquet et al., 2014).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
!	   ( i.e. absolute pressure - 10.1325 dbar )
! 
! rho    : in-situ density				   [kg/m]
*/
double
gsw_rho(double sa, double ct, double p)
{
	return (1.0/gsw_specvol(sa,ct,p));
}
