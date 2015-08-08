/*
!==========================================================================
elemental subroutine gsw_ct_freezing_first_derivatives_poly (sa, p, &
                          saturation_fraction, ctfreezing_sa, ctfreezing_p)
!==========================================================================
!
!  Calculates the first derivatives of the Conservative Temperature at
!  which seawater freezes, with respect to Absolute Salinity SA and
!  pressure P (in Pa) of the comptationally efficient polynomial fit of the
!  freezing temperature (McDougall et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!  saturation_fraction = the saturation fraction of dissolved air in 
!                        seawater
!
!  CTfreezing_SA = the derivative of the Conservative Temperature at
!                  freezing (ITS-90) with respect to Absolute Salinity at
!                  fixed pressure              [ K/(g/kg) ] i.e. [ K kg/g ]
!
!  CTfreezing_P  = the derivative of the Conservative Temperature at
!                  freezing (ITS-90) with respect to pressure (in Pa) at
!                  fixed Absolute Salinity                         [ K/Pa ]
!--------------------------------------------------------------------------
*/
void
gsw_ct_freezing_first_derivatives_poly(double sa, double p,
		double saturation_fraction, double *ctfreezing_sa,
		double *ctfreezing_p)
{
	GSW_TEOS10_CONSTANTS;
	GSW_FREEZING_POLY_COEFFICIENTS;
	double	p_r, sa_r, x, d = -a - a*b - 2.4*b/gsw_sso,
		e = 2.0*a*b/gsw_sso;

	sa_r = sa*1e-2;
	x = sqrt(sa_r);
	p_r = p*1e-4;

	if (ctfreezing_sa != NULL) *ctfreezing_sa = 
	    (c1 + x*(1.5*c2 + x*(2.0*c3 + x*(2.5*c4 + x*(3.0*c5
	        + 3.5*c6*x)))) + p_r*(c10 + x*(1.5*c11 + x*(2.0*c13
		+ x*(2.5*c16 + x*(3.0*c19 + 3.5*c22*x)))) 
	        + p_r*(c12 + x*(1.5*c14 + x*(2.0*c17 + 2.5*c20*x))
	        + p_r*(c15 + x*(1.5*c18 + 2.0*c21*x)))))*1e-2
		- saturation_fraction*1e-3*(d - sa*e);

	if (ctfreezing_p != NULL) *ctfreezing_p =
	    (c7 + sa_r*(c10 + x*(c11 + x*(c13 + x*(c16 + x*(c19 + c22*x)))))
	        + p_r*(2.0*c8 + sa_r*(2.0*c12 + x*(2.0*c14 + x*(2.0*c17
		+ 2.0*c20*x))) + p_r*(3.0*c9 + sa_r*(3.0*c15 + x*(3.0*c18
		+ 3.0*c21*x)))))*1e-8;
}
