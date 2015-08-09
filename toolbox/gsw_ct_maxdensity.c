/*
!==========================================================================
elemental function gsw_ct_maxdensity (sa, p)
!==========================================================================
!
!  Calculates the Conservative Temperature of maximum density of seawater. 
!  This function returns the Conservative temperature at which the density
!  of seawater is a maximum, at given Absolute Salinity, SA, and sea 
!  pressure, p (in dbar).
!
!  SA =  Absolute Salinity                                         [ g/kg ]
!  p  =  sea pressure                                              [ dbar ]
!        ( i.e. absolute pressure - 10.1325 dbar ) 
!
!  CT_maxdensity  =  Conservative Temperature at which            [ deg C ]
!                    the density of seawater is a maximum for
!                    given Absolute Salinity and pressure.
!--------------------------------------------------------------------------
*/
double
gsw_ct_maxdensity(double sa, double p)
{
	int	number_of_iterations;
	double	alpha, ct, ct_mean, ct_old, dalpha_dct,
		dct = 0.001;

	ct = 3.978 - 0.22072*sa;         /*the initial guess of ct.*/

	dalpha_dct = 1.1e-5;             /*the initial guess for dalpha_dct.*/

	for (number_of_iterations = 1; number_of_iterations <= 3;
	    number_of_iterations++) {
	    ct_old = ct;
	    alpha = gsw_alpha(sa,ct_old,p);
	    ct = ct_old - alpha/dalpha_dct;
	    ct_mean = 0.5*(ct + ct_old);
	    dalpha_dct = (gsw_alpha(sa,ct_mean+dct,p)
	                  - gsw_alpha(sa,ct_mean-dct,p))/(dct + dct);
	    ct = ct_old - alpha/dalpha_dct;
	}
	/*
	! After three iterations of this modified Newton-Raphson (McDougall and 
	! Wotherspoon, 2012) iteration, the error in CT_maxdensity is typically
	! no larger than 1x10^-15 degress C.  
	*/
	return (ct);
}
