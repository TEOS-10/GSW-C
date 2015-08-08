/*
!==========================================================================
elemental function gsw_pt0_cold_ice_poly (pot_enthalpy_ice)
!==========================================================================
!
!  Calculates an initial estimate of pt0_ice when it is less than about
!  -100 deg C. 
!
!  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
!
!  pt0_cold_ice_poly  =  initial estimate of potential temperatur 
!                        of very cold ice in dgress C (not K)     [ deg C ] 
!--------------------------------------------------------------------------
*/
double
gsw_pt0_cold_ice_poly(double pot_enthalpy_ice)
{
	GSW_TEOS10_CONSTANTS;
	double	log_abs_theta0, log_h_diff,
		/*h00 = gsw_enthalpy_ice(-gsw_t0,0)*/
		h00 = -6.320202333358860e5,

		s0 =  1.493103204647916,
		s1 =  2.372788609320607e-1,
		s2 = -2.014996002119374e-3,
		s3 =  2.640600197732682e-6,
		s4 =  3.134706016844293e-5,
		s5 =  2.733592344937913e-6,
		s6 =  4.726828010223258e-8,
		s7 = -2.735193883189589e-9,
		s8 = -8.547714991377670e-11;

	log_h_diff = log(pot_enthalpy_ice - h00);

	log_abs_theta0 = s0 + log_h_diff*(s1 + log_h_diff*(s2 + log_h_diff*(s3
	                + log_h_diff*(s4 + log_h_diff*(s5 + log_h_diff*(s6
	                + log_h_diff*(s7 + log_h_diff*s8)))))));

	return (exp(log_abs_theta0) - gsw_t0);
}
