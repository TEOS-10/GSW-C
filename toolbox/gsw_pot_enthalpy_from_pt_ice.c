/*
!==========================================================================
elemental function gsw_pot_enthalpy_from_pt_ice (pt0_ice)
!==========================================================================
!
!  Calculates the potential enthalpy of ice from potential temperature of
!  ice (whose reference sea pressure is zero dbar).  
!
!  pt0_ice  =  potential temperature of ice (ITS-90)              [ deg C ]
!
!  gsw_pot_enthalpy_ice  =  potential enthalpy of ice              [ J/kg ]
!--------------------------------------------------------------------------
*/
double
gsw_pot_enthalpy_from_pt_ice(double pt0_ice)
{
	GSW_TEOS10_CONSTANTS;
	GSW_GIBBS_ICE_COEFFICIENTS;
	double	tau;
	double complex	h0_part, sqtau_t1, sqtau_t2;

	tau = (pt0_ice + gsw_t0)*rec_tt;

	sqtau_t1 = (tau/t1)*(tau/t1);
	sqtau_t2 = (tau/t2)*(tau/t2);

	h0_part = r1*t1*(clog(1.0 - sqtau_t1) + sqtau_t1)
	          + r20*t2*(clog(1.0 - sqtau_t2) + sqtau_t2);

	return (g00 + tt*creal(h0_part));
}
