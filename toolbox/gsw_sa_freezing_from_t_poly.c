/*
!==========================================================================
elemental function gsw_sa_freezing_from_t_poly (t, p, saturation_fraction)
!==========================================================================
!
!  Calculates the Absolute Salinity of seawater at the freezing temperature.
!  That is, the output is the Absolute Salinity of seawater, with the
!  fraction saturation_fraction of dissolved air, that is in equilibrium
!  with ice at in-situ temperature t and pressure p.  If the input values
!  are such that there is no positive value of Absolute Salinity for which
!  seawater is frozen, the output is put equal to Nan.
!
!  t  =  in-situ Temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!        ( i.e. absolute pressure - 10.1325 dbar )
!  saturation_fraction = the saturation fraction of dissolved air in
!                        seawater
!
!  sa_freezing_from_t_poly  =  Absolute Salinity of seawater when it freezes,
!                for given input values of in situ temperature, pressure and
!                air saturation fraction.                          [ g/kg ]
!--------------------------------------------------------------------------
*/
double
gsw_sa_freezing_from_t_poly(double t, double p, double saturation_fraction)
{
	int	i_iter, number_of_iterations = 5;
	double	dt_dsa, sa, sa_old, sa_mean, t_freezing, t_freezing_zero_sa;
	/*
	! This is the band of sa within +- 2.5 g/kg of sa = 0, which we treat
	! differently in calculating the initial values of both SA and dCT_dSA.
	*/
	double	sa_cut_off = 2.5;
	/*
	! Find t > t_freezing_zero_SA.  If this is the case, the input values
	! represent seawater that is not frozen (at any positive SA).
	*/
	t_freezing_zero_sa = gsw_t_freezing_poly(0.0,p,saturation_fraction);
	if (t > t_freezing_zero_sa)
	    return (GSW_INVALID_VALUE);
	/*
	! This is the inital guess of SA using a purpose-built
	! polynomial in CT and p
	*/
	sa = gsw_sa_freezing_estimate(p,saturation_fraction,NULL,&t);
	if (sa < -sa_cut_off)
	    return (GSW_INVALID_VALUE);
	/*
	! Form the first estimate of dt_dSA, the derivative of t with respect
	! to SA at fixed p.
	*/
	sa = max(sa,0.0);
	gsw_t_freezing_first_derivatives_poly(sa,p,saturation_fraction,
	                                           &dt_dsa, NULL);
	/*
	! For -SA_cut_off < SA < SA_cut_off, replace the above estimate of SA
	! with one based on (t_freezing_zero_SA - t).
	*/
	if (fabs(sa) < sa_cut_off)
	    sa = (t - t_freezing_zero_sa)/dt_dsa;
	/*
	!-----------------------------------------------------------------------
	! Begin the modified Newton-Raphson method to find the root of
	! t_freezing = t for SA.
	!-----------------------------------------------------------------------
	*/
	for (i_iter = 1; i_iter <= number_of_iterations; i_iter++) {
	    sa_old = sa;
	    t_freezing = gsw_t_freezing_poly(sa_old,p,saturation_fraction);
	    sa = sa_old - (t_freezing - t)/dt_dsa;
	    sa_mean = 0.5*(sa + sa_old);
	    gsw_t_freezing_first_derivatives_poly(sa_mean,p,saturation_fraction,
	                                               &dt_dsa, NULL);
	    sa = sa_old - (t_freezing - t)/dt_dsa;
	}

	if (gsw_sa_p_inrange(sa,p))
	    return (sa);
	return (GSW_INVALID_VALUE);
}
