/*
!==========================================================================
elemental function gsw_pt_from_pot_enthalpy_ice_poly (pot_enthalpy_ice)
!==========================================================================
!
!  Calculates the potential temperature of ice (whose reference sea 
!  pressure is zero dbar) from the potential enthalpy of ice.  This is a
!  compuationally efficient polynomial fit to the potential enthalpy of
!  ice.
!
!  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
!
!  pt0_ice  =  potential temperature of ice (ITS-90)              [ deg C ]
!--------------------------------------------------------------------------
*/
double
gsw_pt_from_pot_enthalpy_ice_poly(double pot_enthalpy_ice)
{
	double	q0 = 2.533588268773218e2,
		q1 = 2.594351081876611e-3,
		q2 = 1.765077810213815e-8,
		q3 = 7.768070564290540e-14,
		q4 = 2.034842254277530e-19,
		q5 = 3.220014531712841e-25,
		q6 = 2.845172809636068e-31,
		q7 = 1.094005878892950e-37;
/*    
! The error of this fit ranges between -5e-5 and 2e-4 deg C over the potential 
! temperature range of -100 to 2 deg C, or the potential enthalpy range of 
! -5.7 x 10^5 to -3.3 x 10^5 J/kg. 
*/ 
	return (q0
         + pot_enthalpy_ice*(q1 + pot_enthalpy_ice*(q2 + pot_enthalpy_ice*(q3
         + pot_enthalpy_ice*(q4 + pot_enthalpy_ice*(q5 + pot_enthalpy_ice*(q6
         + pot_enthalpy_ice*q7)))))));
}
