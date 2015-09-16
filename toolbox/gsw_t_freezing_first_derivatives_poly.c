/*
!==========================================================================
elemental subroutine gsw_t_freezing_first_derivatives_poly (sa, p, &
                            saturation_fraction, tfreezing_sa, tfreezing_p)
!==========================================================================
!
!  Calculates the first derivatives of the in-situ temperature at which 
!  seawater freezes with respect to Absolute Salinity SA and pressure P (in
!  Pa).  These expressions come from differentiating the expression that
!  defines the freezing temperature, namely the equality between the 
!  chemical potentials of water in seawater and in ice.  
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!  saturation_fraction = the saturation fraction of dissolved air in 
!                        seawater
!
!  tfreezing_SA = the derivative of the in-situ freezing temperature 
!                 (ITS-90) with respect to Absolute Salinity at fixed    
!                 pressure                     [ K/(g/kg) ] i.e. [ K kg/g ]
!
!  tfreezing_P  = the derivative of the in-situ freezing temperature  
!                 (ITS-90) with respect to pressure (in Pa) at fixed  
!                 Absolute Salinity                                [ K/Pa ]
!--------------------------------------------------------------------------
*/
void
gsw_t_freezing_first_derivatives_poly(double sa, double p,
	double saturation_fraction, double *tfreezing_sa, double *tfreezing_p)
{
	GSW_TEOS10_CONSTANTS;
	GSW_FREEZING_POLY_COEFFICIENTS;
	double	p_r, sa_r, x, c = 1e-3/(2.0*gsw_sso);

	sa_r = sa*1e-2;
	x = sqrt(sa_r);
	p_r = p*1e-4;

	if (tfreezing_sa != NULL)
	    *tfreezing_sa =
	    (t1 + x*(1.5*t2 + x*(2.0*t3 + x*(2.5*t4 + x*(3.0*t5
	        + 3.5*t6*x)))) + p_r*(t10 + x*(1.5*t11 + x*(2.0*t13
		+ x*(2.5*t16 + x*(3.0*t19 + 3.5*t22*x))))
	        + p_r*(t12 + x*(1.5*t14 + x*(2.0*t17 + 2.5*t20*x))
	        + p_r*(t15 + x*(1.5*t18 + 2.0*t21*x)))))*1e-2
	        + saturation_fraction*c;

	if (tfreezing_p != NULL)
	    *tfreezing_p =
	    (t7 + sa_r*(t10 + x*(t11 + x*(t13 + x*(t16 + x*(t19 + t22*x)))))
	        + p_r*(2.0*t8 + sa_r*(2.0*t12 + x*(2.0*t14 + x*(2.0*t17
		+ 2.0*t20*x))) + p_r*(3.0*t9 + sa_r*(3.0*t15 + x*(3.0*t18
		+ 3.0*t21*x)))))*1e-8;

	return;
}
