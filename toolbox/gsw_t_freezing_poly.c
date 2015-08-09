/*
!==========================================================================
elemental function gsw_t_freezing_poly (sa, p, saturation_fraction, polynomial)
!==========================================================================
!
!  Calculates the in-situ temperature at which seawater freezes from a 
!  computationally efficient polynomial.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!  saturation_fraction = the saturation fraction of dissolved air in 
!                        seawater
!
!  t_freezing = in-situ temperature at which seawater freezes.    [ deg C ]
!               (ITS-90)                
!--------------------------------------------------------------------------
*/
double
gsw_t_freezing_poly(double sa, double p, double saturation_fraction,
			int polynomial)
{
	GSW_TEOS10_CONSTANTS;
	GSW_FREEZING_POLY_COEFFICIENTS;
	double	p_r, sa_r, x, ctf, sfrac, return_value;
	int	direct_poly=polynomial;

	if (! direct_poly) {
	    sfrac = saturation_fraction;
	    ctf = gsw_ct_freezing_poly(sa,p,sfrac);
	    return_value = gsw_t_from_ct(sa,ctf,p);
	} else {
	    /* Alternative calculation ... */
	    sa_r = sa*1e-2;
	    x = sqrt(sa_r);
	    p_r = p*1e-4;

	    return_value = t0
		+ sa_r*(t1 + x*(t2 + x*(t3 + x*(t4 + x*(t5 + t6*x)))))
		+ p_r*(t7 + p_r*(t8 + t9*p_r))
		+ sa_r*p_r*(t10 + p_r*(t12 + p_r*(t15 + t21*sa_r))
		+ sa_r*(t13 + t17*p_r + t19*sa_r)
		+ x*(t11 + p_r*(t14 + t18*p_r) + sa_r*(t16 + t20*p_r
		+ t22*sa_r)));

	    /* Adjust for the effects of dissolved air */
	    return_value = return_value -
                  saturation_fraction*(1e-3)*(2.4 - sa/(2.0*gsw_sso));
	}
	return (return_value);
}
