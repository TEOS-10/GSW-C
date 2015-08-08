/*
!==========================================================================
elemental function gsw_ct_from_enthalpy (sa, h, p)
!==========================================================================
!
!  Calculates the Conservative Temperature of seawater, given the Absolute 
!  Salinity, specific enthalpy, h, and pressure p.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  h   =  specific enthalpy                                        [ J/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325d0 dbar ) 
!
!  CT  =  Conservative Temperature ( ITS-90)                      [ deg C ]
!--------------------------------------------------------------------------
*/
double
gsw_ct_from_enthalpy(double sa, double h, double p)
{
	GSW_TEOS10_CONSTANTS;
	double	ct, ct_freezing, ct_mean, ct_old, f, h_freezing, h_ct, h_40,
		ct_40 = 40.0;

	ct_freezing = gsw_ct_freezing(sa,p,0.0);

	h_freezing = gsw_enthalpy(sa,ct_freezing,p);
	if (h < (h_freezing - gsw_cp0)) {
	    /*
	    ! The input, seawater enthalpy h, is less than the enthalpy at the
	    ! freezing temperature, i.e. the water is frozen.
	    */
	    return (GSW_INVALID_VALUE);
	}

	h_40 = gsw_enthalpy(sa,ct_40,p);
	if (h > h_40) {
	    /*
	    ! The input seawater enthalpy is greater than the enthalpy
	    ! when CT is 40C
	    */
	    return (GSW_INVALID_VALUE);
	}

	/* first guess of ct */
	ct = ct_freezing + (ct_40 - ct_freezing)*(h - h_freezing)/
				(h_40 - h_freezing);
	gsw_enthalpy_first_derivatives(sa,ct,p,NULL,&h_ct);

	/*
	!------------------------------------------------------
	! Begin the modified Newton-Raphson iterative procedure 
	!------------------------------------------------------
	*/

	ct_old = ct;
	f = gsw_enthalpy(sa,ct_old,p) - h;
	ct = ct_old - f/h_ct;
	ct_mean = 0.5*(ct + ct_old);
	gsw_enthalpy_first_derivatives(sa,ct_mean,p,NULL,&h_ct);
	ct = ct_old - f/h_ct;

	ct_old = ct;
	f = gsw_enthalpy(sa,ct_old,p) - h;
	ct = ct_old - f/h_ct;
	/*
	! After 1.5d0 iterations of this modified Newton-Raphson iteration,
	! the error in CT is no larger than 4x10^-13 degrees C, which 
	! is machine precision for this calculation. 
	*/
	return (ct);
}
