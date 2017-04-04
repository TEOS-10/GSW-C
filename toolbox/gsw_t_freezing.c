/*
!==========================================================================
function gsw_t_freezing(sa,p,saturation_fraction)
!==========================================================================

! Calculates the in-situ temperature at which seawater freezes
!
! sa     : Absolute Salinity                                 [g/kg]
! p      : sea pressure                                      [dbar]
!         ( i.e. absolute pressure - 10.1325 dbar )
! saturation_fraction : the saturation fraction of dissolved air
!                       in seawater
!
! t_freezing : in-situ temperature at which seawater freezes.[deg C]
*/
double
gsw_t_freezing(double sa, double p, double saturation_fraction)
{
	GSW_TEOS10_CONSTANTS;
	GSW_FREEZING_POLY_COEFFICIENTS;
	double sa_r, x, p_r;
	double	df_dt, tf, tfm, tf_old, f, return_value;

	/* The initial value of t_freezing_exact (for air-free seawater) */
	sa_r = sa*1e-2;
	x = sqrt(sa_r);
	p_r = p*1e-4;

	tf = t0
	+ sa_r*(t1 + x*(t2 + x*(t3 + x*(t4 + x*(t5 + t6*x)))))
	+ p_r*(t7 + p_r*(t8 + t9*p_r))
	+ sa_r*p_r*(t10 + p_r*(t12 + p_r*(t15 + t21*sa_r))
	+ sa_r*(t13 + t17*p_r + t19*sa_r)
	+ x*(t11 + p_r*(t14 + t18*p_r) + sa_r*(t16 + t20*p_r
	+ t22*sa_r)));

	/* Adjust for the effects of dissolved air */
	tf -= saturation_fraction*(1e-3)*(2.4 - sa/(2.0*gsw_sso));

	df_dt = 1e3*gsw_t_deriv_chem_potential_water_t_exact(sa,tf,p) -
		gsw_gibbs_ice(1,0,tf,p);
/*
! df_dt here is the initial value of the derivative of the function f whose
! zero (f = 0) we are finding (see Eqn. (3.33.2) of IOC et al (2010)).
*/

	tf_old = tf;
	f = 1e3*gsw_chem_potential_water_t_exact(sa,tf_old,p) -
		gsw_gibbs_ice(0,0,tf_old,p);
	tf = tf_old - f/df_dt;
	tfm = 0.5*(tf + tf_old);
	df_dt = 1e3*gsw_t_deriv_chem_potential_water_t_exact(sa,tfm,p) -
		gsw_gibbs_ice(1,0,tfm,p);
	tf = tf_old - f/df_dt;

	tf_old = tf;
	f = 1e3*gsw_chem_potential_water_t_exact(sa,tf_old,p) -
		gsw_gibbs_ice(0,0,tf_old,p);
	tf = tf_old - f/df_dt;

	/* Adjust for the effects of dissolved air */
	return_value = tf -
                saturation_fraction*(1e-3)*(2.4 - sa/(2.0*gsw_sso));
	return (return_value);
}
