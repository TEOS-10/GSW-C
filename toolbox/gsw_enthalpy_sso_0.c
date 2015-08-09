/*
!==========================================================================
function gsw_enthalpy_sso_0(p)
!==========================================================================

!  This function calculates enthalpy at the Standard Ocean Salinity, SSO,
!  and at a Conservative Temperature of zero degrees C, as a function of
!  pressure, p, in dbar, using a streamlined version of the
!  computationally-efficient expression for specific volume, that is, a
!  streamlined version of the code "gsw_enthalpy(SA,CT,p)".
!
! p      : sea pressure                                    [dbar]
!
! enthalpy_sso_0 : enthalpy(sso,0,p)
*/
double
gsw_enthalpy_sso_0(double p)
{
	GSW_TEOS10_CONSTANTS;
	GSW_SPECVOL_COEFFICIENTS;
	double	dynamic_enthalpy_sso_0_p, z;

	z = p*1.0e-4;

	dynamic_enthalpy_sso_0_p =
			z*( 9.726613854843870e-4 + z*(-2.252956605630465e-5
			+ z*( 2.376909655387404e-6 + z*(-1.664294869986011e-7
			+ z*(-5.988108894465758e-9 + z*(h006 + h007*z))))));
	return (dynamic_enthalpy_sso_0_p*db2pa*1.0e4);
}
