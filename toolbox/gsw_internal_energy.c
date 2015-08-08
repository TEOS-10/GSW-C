/*
!==========================================================================
function gsw_internal_energy(sa,ct,p)  
!==========================================================================

!  Calculates internal energy of seawater.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
! 
! internal_energy  :  internal_energy of seawater          [J/kg]
*/
double
gsw_internal_energy(double sa, double ct, double p)
{
	GSW_TEOS10_CONSTANTS;

	return (gsw_enthalpy(sa,ct,p) - (gsw_p0 + db2pa*p)
		*gsw_specvol(sa,ct,p));
}
