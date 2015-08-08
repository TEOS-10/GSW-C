/*
!==========================================================================
function gsw_enthalpy(sa,ct,p)  
!==========================================================================

!  Calculates specific enthalpy of seawater using the computationally-
!  efficient expression for specific volume in terms of SA, CT and p
!  (Roquet et al., 2014).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
!         ( i.e. absolute pressure - 10.1325 dbar )
! 
! enthalpy  :  specific enthalpy of seawater               [J/kg]
*/
double
gsw_enthalpy(double sa, double ct, double p)
{
	GSW_TEOS10_CONSTANTS;
	return (gsw_cp0*ct + gsw_dynamic_enthalpy(sa,ct,p));
}
