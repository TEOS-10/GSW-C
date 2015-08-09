/*
!==========================================================================
elemental function gsw_pot_enthalpy_from_pt_ice_poly (pt0_ice)
!==========================================================================
!
!  Calculates the potential enthalpy of ice from potential temperature of
!  ice (whose reference sea pressure is zero dbar).  This is a
!  compuationally efficient polynomial fit to the potential enthalpy of
!  ice.
!   
!  pt0_ice  =  potential temperature of ice (ITS-90)              [ deg C ]
!
!  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
!--------------------------------------------------------------------------
*/
double
gsw_pot_enthalpy_from_pt_ice_poly(double pt0_ice)
{
	int	iteration;
	double	df_dt, f, pot_enthalpy_ice,
		pot_enthalpy_ice_mid, pot_enthalpy_ice_old,
		p0 = -3.333601570157700e5,
		p1 =  2.096693916810367e3,
		p2 =  3.687110754043292,
		p3 =  4.559401565980682e-4,
		p4 = -2.516011957758120e-6,
		p5 = -1.040364574632784e-8,
		p6 = -1.701786588412454e-10,
		p7 = -7.667191301635057e-13;
    
	/*initial estimate of the potential enthalpy.*/
	pot_enthalpy_ice = p0 + pt0_ice*(p1 + pt0_ice*(p2 + pt0_ice*(p3
	                   + pt0_ice*(p4 + pt0_ice*(p5 + pt0_ice*(p6
	                   + pt0_ice*p7))))));

	df_dt = gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice);

	for (iteration = 1; iteration <= 5; iteration++) {
	    pot_enthalpy_ice_old = pot_enthalpy_ice;
	    f = gsw_pt_from_pot_enthalpy_ice_poly(pot_enthalpy_ice_old)
			- pt0_ice;
	    pot_enthalpy_ice = pot_enthalpy_ice_old - f/df_dt;
	    pot_enthalpy_ice_mid = 0.5*(pot_enthalpy_ice+pot_enthalpy_ice_old);
	    df_dt = gsw_pt_from_pot_enthalpy_ice_poly_dh(pot_enthalpy_ice_mid);
	    pot_enthalpy_ice = pot_enthalpy_ice_old - f/df_dt;
	}
	/*
	! The error of this fit ranges between -6e-3 and 6e-3 J/kg over the
	| potential temperature range of -100 to 2 deg C, or the potential
	| enthalpy range of -5.7 x 10^5 to -3.3 x 10^5 J/kg. 
	*/
	return (pot_enthalpy_ice);
}
