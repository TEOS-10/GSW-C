/*
!==========================================================================
elemental function gsw_sa_freezing_from_ct_poly (ct, p, saturation_fraction)
!==========================================================================
!
!  Calculates the Absolute Salinity of seawater at the freezing temperature.  
!  That is, the output is the Absolute Salinity of seawater, with the 
!  fraction saturation_fraction of dissolved air, that is in equilibrium 
!  with ice at Conservative Temperature CT and pressure p.  If the input 
!  values are such that there is no positive value of Absolute Salinity for
!  which seawater is frozen, the output is put equal to Nan.
!
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!  saturation_fraction  =  the saturation fraction of dissolved air in 
!                          seawater
!
!  sa_freezing_from_ct  =  Absolute Salinity of seawater when it freezes,
!                 for given input values of Conservative Temperature
!                 pressure and air saturation fraction.            [ g/kg ]
!--------------------------------------------------------------------------
*/
double
gsw_sa_freezing_from_ct_poly(double ct, double p, double saturation_fraction)
{
	int	i_iter, number_of_iterations = 2;
	double	ct_freezing, ct_freezing_zero_sa, dct_dsa, sa, sa_old, sa_mean;
	/*
	! This is the band of sa within +- 2.5 g/kg of sa = 0, which we treat
	! differently in calculating the initial values of both SA and dCT_dSA. 
	*/
	double	sa_cut_off = 2.5;
	/*
	! Find CT > CT_freezing_zero_SA.  If this is the case, the input values
	! represent seawater that is not frozen (at any positive SA). 
	*/
	ct_freezing_zero_sa = gsw_ct_freezing_poly(0.0,p,saturation_fraction);
	if (ct > ct_freezing_zero_sa)
	    return (GSW_INVALID_VALUE);

	/*Form the first estimate of SA from a polynomial in CT and p */
	sa = gsw_sa_freezing_estimate(p,saturation_fraction,&ct,NULL);
	if (sa < -sa_cut_off)
	    return (GSW_INVALID_VALUE);
	/*
	! Form the first estimate of dCT_dSA, the derivative of CT with respect 
	! to SA at fixed p.  
	*/
	sa = max(sa,0.0);
	gsw_ct_freezing_first_derivatives_poly(sa,p,saturation_fraction,
	                                            &dct_dsa, NULL);
	/*
	! For -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA  
	! with one based on (CT_freezing_zero_SA - CT).
	*/
	if (fabs(sa) < sa_cut_off)
	    sa = (ct - ct_freezing_zero_sa)/dct_dsa;
	/*
	!-----------------------------------------------------------------------
	! Begin the modified Newton-Raphson method to solve the root of 
	! CT_freezing = CT for SA. 
	!-----------------------------------------------------------------------
	*/
	for (i_iter = 1; i_iter <= number_of_iterations; i_iter++) {
	    sa_old = sa;
	    ct_freezing = gsw_ct_freezing_poly(sa_old,p,saturation_fraction);
	    sa = sa_old - (ct_freezing - ct)/dct_dsa;
	    sa_mean = 0.5*(sa + sa_old);
	    gsw_ct_freezing_first_derivatives_poly(sa_mean,p,
				saturation_fraction, &dct_dsa, NULL);
	    sa = sa_old - (ct_freezing - ct)/dct_dsa;
	}

	if (gsw_sa_p_inrange(sa,p))
	    return (sa);
	return (GSW_INVALID_VALUE);
}
