/*
!==========================================================================
elemental subroutine gsw_enthalpy_first_derivatives (sa, ct, p, h_sa, h_ct)
!==========================================================================
!
!  Calculates the following two derivatives of specific enthalpy (h) of
!  seawater using the computationally-efficient expression for 
!  specific volume in terms of SA, CT and p (Roquet et al., 2014).  
!   (1) h_SA, the derivative with respect to Absolute Salinity at 
!       constant CT and p, and
!   (2) h_CT, derivative with respect to CT at constant SA and p. 
!  Note that h_P is specific volume (1/rho) it can be caclulated by calling
!  gsw_specvol(SA,CT,p).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  h_SA  =  The first derivative of specific enthalpy with respect to 
!           Absolute Salinity at constant CT and p.     
!                                            [ J/(kg (g/kg))]  i.e. [ J/g ]
!  h_CT  =  The first derivative of specific enthalpy with respect to 
!           CT at constant SA and p.                           [ J/(kg K) ]
!--------------------------------------------------------------------------
*/
void
gsw_enthalpy_first_derivatives(double sa, double ct, double p, double *h_sa,
				double *h_ct)
{
	GSW_TEOS10_CONSTANTS;
	GSW_SPECVOL_COEFFICIENTS;
	double	dynamic_h_ct_part, dynamic_h_sa_part, xs, ys, z;

	xs = sqrt(gsw_sfac*sa + offset);
	ys = ct*0.025;
	z = p*1e-4;

	if (h_sa != NULL) {

	    dynamic_h_sa_part =  z*(h101 + xs*(2.0*h201 + xs*(3.0*h301
	        + xs*(4.0*h401 + xs*(5.0*h501 + 6.0*h601*xs)))) + ys*(h111
		+ xs*(2.0*h211 + xs*(3.0*h311 + xs*(4.0*h411
		+ 5.0*h511*xs))) + ys*(h121 + xs*(2.0*h221 + xs*(3.0*h321
	        + 4.0*h421*xs)) + ys*(h131 + xs*(2.0*h231 + 3.0*h331*xs)
		+ ys*(h141 + 2.0*h241*xs + h151*ys)))) + z*(h102
		+ xs*(2.0*h202 + xs*(3.0*h302 + xs*(4.0*h402
		+ 5.0*h502*xs))) + ys*(h112 + xs*(2.0*h212 + xs*(3.0*h312
	        + 4.0*h412*xs)) + ys*(h122 + xs*(2.0*h222 + 3.0*h322*xs)
		+ ys*(h132 + 2.0*h232*xs + h142*ys ))) + z*(h103 + xs*(2.0*h203
		+ xs*(3.0*h303 + 4.0*h403*xs)) + ys*(h113 + xs*(2.0*h213
		+ 3.0*h313*xs) + ys*(h123 + 2.0*h223*xs + h133*ys))
		+ z*(h104 + 2.0*h204*xs + h114*ys + h105*z))));

	    *h_sa = 1e8*0.5*gsw_sfac*dynamic_h_sa_part/xs;

	}

	if (h_ct != NULL) {

	    dynamic_h_ct_part = z*(h011 + xs*(h111
		+ xs*(h211 + xs*(h311 + xs*(h411
	        + h511*xs)))) + ys*(2.0*(h021 + xs*(h121 + xs*(h221 + xs*(h321
	        + h421*xs)))) + ys*(3.0*(h031 + xs*(h131 + xs*(h231 + h331*xs)))
	        + ys*(4.0*(h041 + xs*(h141 + h241*xs)) + ys*(5.0*(h051
		+ h151*xs) + 6.0*h061*ys)))) + z*(h012 + xs*(h112 + xs*(h212
		+ xs*(h312 + h412*xs))) + ys*(2.0*(h022 + xs*(h122 + xs*(h222
		+ h322*xs))) + ys*(3.0*(h032 + xs*(h132 + h232*xs))
		+ ys*(4.0*(h042 + h142*xs) + 5.0*h052*ys))) + z*(h013
		+ xs*(h113 + xs*(h213 + h313*xs)) + ys*(2.0*(h023 + xs*(h123
	        + h223*xs)) + ys*(3.0*(h033 + h133*xs) + 4.0*h043*ys))
		+ z*(h014 + h114*xs + 2.0*h024*ys + h015*z ))));

	    *h_ct = gsw_cp0 + 1e8*0.025*dynamic_h_ct_part;

	}
}
