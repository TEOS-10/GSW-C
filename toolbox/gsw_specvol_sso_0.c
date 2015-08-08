/*
!==========================================================================
function gsw_specvol_sso_0(p) 
!==========================================================================

!  This function calculates specific volume at the Standard Ocean Salinty,
!  SSO, and at a Conservative Temperature of zero degrees C, as a function 
!  of pressure, p, in dbar, using a streamlined version of the CT version
!  of specific volume, that is, a streamlined version of the code
!  "gsw_specvol(SA,CT,p)".
!
! p      : sea pressure                                    [dbar]
! 							     3   -1
! specvol_sso_0 : specvol(sso,0,p)                         [m  kg  ]
*/
double
gsw_specvol_sso_0(double p)
{
	GSW_SPECVOL_COEFFICIENTS;
	double	z, return_value;

	z = p*1.0e-4;

	return_value = 9.726613854843870e-04 + z*(-4.505913211160929e-05
		+ z*(7.130728965927127e-06 + z*(-6.657179479768312e-07
		+ z*(-2.994054447232880e-08 + z*(v005 + v006*z)))));
	return (return_value);
}
