/*
!==========================================================================
elemental subroutine gsw_enthalpy_second_derivatives (sa, ct, p, h_sa_sa, &
                                                      h_sa_ct, h_ct_ct)
! =========================================================================
!
!  Calculates the following three second-order derivatives of specific
!  enthalpy (h),using the computationally-efficient expression for 
!  specific volume in terms of SA, CT and p (Roquet et al., 2014).
!   (1) h_SA_SA, second-order derivative with respect to Absolute Salinity 
!       at constant CT & p.
!   (2) h_SA_CT, second-order derivative with respect to SA & CT at 
!       constant p. 
!   (3) h_CT_CT, second-order derivative with respect to CT at constant SA 
!       and p. 
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  h_SA_SA  =  The second derivative of specific enthalpy with respect to 
!              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
!  h_SA_CT  =  The second derivative of specific enthalpy with respect to 
!              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
!  h_CT_CT  =  The second derivative of specific enthalpy with respect to 
!              CT at constant SA and p.                      [ J/(kg K^2) ]
!--------------------------------------------------------------------------
*/
void
gsw_enthalpy_second_derivatives(double sa, double ct, double p,
		double *h_sa_sa, double *h_sa_ct, double *h_ct_ct)
{
	GSW_TEOS10_CONSTANTS;
	GSW_SPECVOL_COEFFICIENTS;
	double	dynamic_h_ct_ct_part, dynamic_h_sa_ct_part,
		dynamic_h_sa_sa_part, xs, xs2, ys, z;

	xs = sqrt(gsw_sfac*sa + offset);
	ys = ct*0.025;
	z = p*1e-4;

	if (h_sa_sa != NULL) {

	    xs2 = pow(xs,2);
	    dynamic_h_sa_sa_part = z*(-h101 + xs2*(3.0*h301 + xs*(8.0*h401
	        + xs*(15.0*h501 + 24.0*h601*xs))) + ys*(- h111
		+ xs2*(3.0*h311 + xs*(8.0*h411 + 15.0*h511*xs)) + ys*(-h121
	        + xs2*(3.0*h321 + 8.0*h421*xs) + ys*(-h131 + 3.0*h331*xs2
		+ ys*(-h141 - h151*ys)))) + z*(-h102 + xs2*(3.0*h302
		+ xs*(8.0*h402 + 15.0*h502*xs)) + ys*(-h112 + xs2*(3.0*h312
		+ 8.0*h412*xs) + ys*(-h122 + 3.0*h322*xs2 + ys*(-h132
		- h142*ys ))) + z*(xs2*(8.0*h403*xs + 3.0*h313*ys)
	        + z*(-h103 + 3.0*h303*xs2 + ys*(-h113 + ys*(-h123 - h133*ys))
	        + z*(-h104 - h114*ys - h105*z)))));

	    *h_sa_sa = 1e8*0.25*gsw_sfac*gsw_sfac*dynamic_h_sa_sa_part/
			pow(xs,3);

	}

	if (h_sa_ct != NULL) {

	    dynamic_h_sa_ct_part = z*(h111 + xs*(2.0*h211 + xs*(3.0*h311
	        + xs*(4.0*h411 + 5.0*h511*xs))) + ys*(2.0*h121
		+ xs*(4.0*h221 + xs*(6.0*h321 + 8.0*h421*xs))
		+ ys*(3.0*h131 + xs*(6.0*h231 + 9.0*h331*xs)
		+ ys*(4.0*h141 + 8.0*h241*xs + 5.0*h151*ys ))) + z*(h112
		+ xs*(2.0*h212 + xs*(3.0*h312 + 4.0*h412*xs))
		+ ys*(2.0*h122 + xs*(4.0*h222 + 6.0*h322*xs)
		+ ys*(3.0*h132 + 6.0*h232*xs + 4.0*h142*ys)) + z*(h113
		+ xs*(2.0*h213 + 3.0*h313*xs) + ys*(2.0*h123
		+ 4.0*h223*xs + 3.0*h133*ys) + h114*z)));

	    *h_sa_ct = 1e8*0.025*0.5*gsw_sfac*dynamic_h_sa_ct_part/xs;

	}

	if (h_ct_ct != NULL) {

	    dynamic_h_ct_ct_part = z*(2.0*h021 + xs*(2.0*h121 + xs*(2.0*h221
	        + xs*(2.0*h321 + 2.0*h421*xs))) + ys*(6.0*h031
		+ xs*(6.0*h131 + xs*(6.0*h231 + 6.0*h331*xs))
	        + ys*(12.0*h041 + xs*(12.0*h141 + 12.0*h241*xs)
		+ ys*(20.0*h051 + 20.0*h151*xs + 30.0*h061*ys)))
		+ z*(2.0*h022 + xs*(2.0*h122 + xs*(2.0*h222
		+ 2.0*h322*xs)) + ys*(6.0*h032 + xs*(6.0*h132
		+ 6.0*h232*xs) + ys*(12.0*h042 + 12.0*h142*xs
		+ 20.0*h052*ys)) + z*(2.0*h023 + xs*(2.0*h123
		+ 2.0*h223*xs) + ys*(6.0*h133*xs + 6.0*h033
		+ 12.0*h043*ys) + 2.0*h024*z)));

	    *h_ct_ct = 1e8*6.25e-4*dynamic_h_ct_ct_part;

	}
}
