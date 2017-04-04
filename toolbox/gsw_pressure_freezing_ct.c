/*
!==========================================================================
elemental function gsw_pressure_freezing_ct (sa, ct, saturation_fraction)
!==========================================================================
!
!  Calculates the pressure (in dbar) of seawater at the freezing
!  temperature.  That is, the output is the pressure at which seawater,
!  with Absolute Salinity SA, Conservative Temperature CT, and with
!  saturation_fraction of dissolved air, freezes.  If the input values are
!  such that there is no value of pressure in the range between 0 dbar and
!  10,000 dbar for which seawater is at the freezing temperature, the
!  output, pressure_freezing, is put equal to NaN.
!
!  SA  =  Absolute Salinity of seawater                            [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  saturation_fraction = the saturation fraction of dissolved air in
!                        seawater
!
!  pressure_freezing = sea pressure at which the seawater freezes  [ dbar ]
!        ( i.e. absolute pressure - 10.1325 dbar )
!--------------------------------------------------------------------------
*/
double
gsw_pressure_freezing_ct(double sa, double ct, double saturation_fraction)
{
	GSW_TEOS10_CONSTANTS;
	int	i_iter, number_of_iterations = 3;
	double	ct_freezing_p0, ct_freezing_p10000, dctf_dp, f,
		pf, pfm, pf_old, ctfreezing_p;

	/*
	! rec_Pa2dbar is to have dCTf_dp in units of K/dbar rather than K/Pa

	! Find CT > CT_freezing_p0.  If this is the case, the input CT value
	! represent seawater that will not be frozen at any positive p.
	*/
	ct_freezing_p0 = gsw_ct_freezing(sa,0.0,saturation_fraction);
	if (ct > ct_freezing_p0) {
	    return (GSW_INVALID_VALUE);
	}
	/*
	! Find CT < CT_freezing_p10000.  If this is the case, the input CT value
	! represent seawater that is frozen even at p = 10,000 dbar.
	*/
	ct_freezing_p10000 = gsw_ct_freezing(sa,1e4,saturation_fraction);
	if (ct < ct_freezing_p10000) {
	    return (GSW_INVALID_VALUE);
	}
	/*
	! This is the initial (linear) guess of the freezing pressure, in dbar.
	*/
	pf = rec_pa2db*(ct_freezing_p0 - ct)/
			(ct_freezing_p0 - ct_freezing_p10000);

	gsw_ct_freezing_first_derivatives(sa,pf,saturation_fraction,
	                                       NULL,&ctfreezing_p);
	dctf_dp = rec_pa2db*ctfreezing_p;
	    /*
	    !  this dctf_dp is the initial value of the partial derivative of
	    !  ct_freezing with respect to pressure (in dbar) at fixed sa,
	    !  assuming that the saturation_fraction is zero.
	    */
	for (i_iter = 1; i_iter <= number_of_iterations; i_iter++) {
	    pf_old = pf;
	    f = gsw_ct_freezing(sa,pf_old,saturation_fraction) - ct;
	    pf = pf_old - f/dctf_dp;
	    pfm = 0.5*(pf + pf_old);
	    gsw_ct_freezing_first_derivatives(sa,pfm,saturation_fraction,
	                                           NULL,&ctfreezing_p);
	    dctf_dp = rec_pa2db*ctfreezing_p;
	    pf = pf_old - f/dctf_dp;
	}

	if (gsw_sa_p_inrange(sa,pf))
	    return (pf);
	return (GSW_INVALID_VALUE);
}
