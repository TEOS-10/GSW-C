/*
!==========================================================================
function gsw_alpha_on_beta(sa,ct,p)  
!==========================================================================

!  Calculates alpha divided by beta, where alpha is the thermal expansion
!  coefficient and beta is the saline contraction coefficient of seawater 
!  from Absolute Salinity and Conservative Temperature.  This function uses
!  the computationally-efficient expression for specific volume in terms of 
!  SA, CT and p (Roquet et al., 2014).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! alpha_on_beta
!        : thermal expansion coefficient with respect to   [kg g^-1 K^-1]
!          Conservative Temperature divided by the saline
!          contraction coefficient at constant Conservative
!          Temperature
*/
double
gsw_alpha_on_beta(double sa, double ct, double p)
{
	GSW_TEOS10_CONSTANTS;
	GSW_SPECVOL_COEFFICIENTS;
	double	xs, ys, z, v_ct_part, v_sa_part;

	xs	= sqrt(gsw_sfac*sa + offset);
	ys	= ct*0.025;
	z	= p*1e-4;

	v_ct_part = a000
		+ xs*(a100 + xs*(a200 + xs*(a300 + xs*(a400 + a500*xs))))
	        + ys*(a010 + xs*(a110 + xs*(a210 + xs*(a310 + a410*xs)))
	        + ys*(a020 + xs*(a120 + xs*(a220 + a320*xs)) + ys*(a030
	        + xs*(a130 + a230*xs) + ys*(a040 + a140*xs + a050*ys ))))
	        + z*(a001 + xs*(a101 + xs*(a201 + xs*(a301 + a401*xs)))
	        + ys*(a011 + xs*(a111 + xs*(a211 + a311*xs)) + ys*(a021
	        + xs*(a121 + a221*xs) + ys*(a031 + a131*xs + a041*ys)))
	        + z*(a002 + xs*(a102 + xs*(a202 + a302*xs)) + ys*(a012
	        + xs*(a112 + a212*xs) + ys*(a022 + a122*xs + a032*ys))
	        + z*(a003 + a103*xs + a013*ys + a004*z)));

	v_sa_part = b000
		+ xs*(b100 + xs*(b200 + xs*(b300 + xs*(b400 + b500*xs))))
		+ ys*(b010 + xs*(b110 + xs*(b210 + xs*(b310 + b410*xs)))
		+ ys*(b020 + xs*(b120 + xs*(b220 + b320*xs)) + ys*(b030
		+ xs*(b130 + b230*xs) + ys*(b040 + b140*xs + b050*ys))))
		+ z*(b001 + xs*(b101 + xs*(b201 + xs*(b301 + b401*xs)))
		+ ys*(b011 + xs*(b111 + xs*(b211 + b311*xs)) + ys*(b021
		+ xs*(b121 + b221*xs) + ys*(b031 + b131*xs + b041*ys)))
		+ z*(b002 + xs*(b102 + xs*(b202 + b302*xs))+ ys*(b012
		+ xs*(b112 + b212*xs) + ys*(b022 + b122*xs + b032*ys))
		+ z*(b003 +  b103*xs + b013*ys + b004*z)));

	return (-(v_ct_part*xs)/(20.0*gsw_sfac*v_sa_part));
}
