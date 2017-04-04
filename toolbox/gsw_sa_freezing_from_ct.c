/*
!==========================================================================
elemental function gsw_sa_freezing_from_ct (ct, p, saturation_fraction)
!==========================================================================
!
!  Calculates the Absolute Salinity of seawater at the freezing temperature.  
!  That is, the output is the Absolute Salinity of seawater, with 
!  Conservative Temperature CT, pressure p and the fraction 
!  saturation_fraction of dissolved air, that is in equilibrium 
!  with ice at the same in situ temperature and pressure.  If the input 
!  values are such that there is no positive value of Absolute Salinity for
!  which seawater is frozen, the output is made a NaN.
!
!  CT  =  Conservative Temperature of seawater (ITS-90)           [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!  saturation_fraction  =  the saturation fraction of dissolved air in 
!                          seawater
!
!  sa_freezing_from_ct  =  Absolute Salinity of seawater when it freezes,
!                 for given input values of its Conservative Temperature,
!                 pressure and air saturation fraction.            [ g/kg ]
!--------------------------------------------------------------------------
*/
double
gsw_sa_freezing_from_ct(double ct, double p, double saturation_fraction)
{
	int	i_iter, number_of_iterations = 3;
	double	ct_freezing_zero_sa, f, ctfreezing_sa, sa, sa_mean, sa_old;
	/*
	! This is the band of sa within +- 2.5 g/kg of sa = 0, which we treat
	! differently in calculating the initial values of both SA and dCT_dSA. 
	*/
	double	sa_cut_off = 2.5;
	/*
	! Find CT > CT_freezing_zero_SA.  If this is the case, the input values
	! represent seawater that is not frozen (for any positive SA). 
	*/
	ct_freezing_zero_sa = gsw_ct_freezing(0.0,p,saturation_fraction);
	if (ct > ct_freezing_zero_sa)
	    return (GSW_INVALID_VALUE);

	/*Form the first estimate of SA from a polynomial in CT and p*/
	sa = gsw_sa_freezing_estimate(p,saturation_fraction,&ct,NULL);
	if (sa < -sa_cut_off)
	    return (GSW_INVALID_VALUE);
	/*
	! Form the first estimate of CTfreezing_SA,
	! the derivative of CT_freezing with respect to SA at fixed p.
	*/
	sa = max(sa,0.0);
	gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,
	                                       &ctfreezing_sa, NULL);
	/*
	! For -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA
	! with one based on (CT_freezing_zero_SA - CT).
	*/
	if (fabs(sa) < sa_cut_off)
	    sa = (ct - ct_freezing_zero_sa)/ctfreezing_sa;
	/*
	!-----------------------------------------------------------------------
	! Begin the modified Newton-Raphson method to solve
	! f = (CT_freezing - CT) = 0 for SA.
	!-----------------------------------------------------------------------
	*/
	for (i_iter = 1; i_iter <= number_of_iterations; i_iter++) {
	    sa_old = sa;
	    f = gsw_ct_freezing(sa,p,saturation_fraction) - ct;
	    sa = sa_old - f/ctfreezing_sa;
	    sa_mean = 0.5*(sa + sa_old);
	    gsw_ct_freezing_first_derivatives(sa_mean,p,saturation_fraction,
					&ctfreezing_sa, NULL);
	    sa = sa_old - f/ctfreezing_sa;
	}

	if (gsw_sa_p_inrange(sa,p))
	    return (sa);
	return (GSW_INVALID_VALUE);
}
