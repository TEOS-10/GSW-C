/*
! =========================================================================
elemental function gsw_gibbs_ice (nt, np, t, p)
! =========================================================================
!
!  Ice specific Gibbs energy and derivatives up to order 2.
!
!  nt  =  order of t derivative                      [ integers 0, 1 or 2 ]
!  np  =  order of p derivative                      [ integers 0, 1 or 2 ]
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!   
!  gibbs_ice = Specific Gibbs energy of ice or its derivatives.
!            The Gibbs energy (when nt = np = 0) has units of:     [ J/kg ]
!            The temperature derivatives are output in units of: 
!                                                      [ (J/kg) (K)^(-nt) ]
!            The pressure derivatives are output in units of:
!                                                     [ (J/kg) (Pa)^(-np) ]
!            The mixed derivatives are output in units of:
!                                           [ (J/kg) (K)^(-nt) (Pa)^(-np) ]
!  Note. The derivatives are taken with respect to pressure in Pa, not
!    withstanding that the pressure input into this routine is in dbar.
!--------------------------------------------------------------------------
*/
double
gsw_gibbs_ice (int nt, int np, double t, double p)
{
	GSW_TEOS10_CONSTANTS;
	GSW_GIBBS_ICE_COEFFICIENTS;
	double	dzi, g0, g0p, g0pp, sqrec_pt;
	double complex	r2, r2p, r2pp, g, sqtau_t1, sqtau_t2, tau,
			tau_t1, tau_t2;
	double	s0 = -3.32733756492168e3;

	tau = (t + gsw_t0)*rec_tt;

	dzi = db2pa*p*rec_pt;

	if (nt == 0 && np == 0) {

	    tau_t1 = tau/t1;
	    sqtau_t1 = tau_t1*tau_t1;
	    tau_t2 = tau/t2;
	    sqtau_t2 = tau_t2*tau_t2;

	    g0 = g00 + dzi*(g01 + dzi*(g02 + dzi*(g03 + g04*dzi)));

	    r2 = r20 + dzi*(r21 + r22*dzi);

	    g = r1*(tau*clog((1.0 + tau_t1)/(1.0 - tau_t1))
	        + t1*(clog(1.0 - sqtau_t1) - sqtau_t1))
		+ r2*(tau*clog((1.0 + tau_t2)/(1.0 - tau_t2))
		+ t2*(clog(1.0 - sqtau_t2) - sqtau_t2));

	    return (g0 - tt*(s0*tau - creal(g)));

	} else if (nt == 1 && np == 0) {

	    tau_t1 = tau/t1;
	    tau_t2 = tau/t2;

	    r2 = r20 + dzi*(r21 + r22*dzi);

	    g = r1*(clog((1.0 + tau_t1)/(1.0 - tau_t1)) - 2.0*tau_t1)
	        + r2*(clog((1.0 + tau_t2)/(1.0 - tau_t2)) - 2.0*tau_t2);

	    return (-s0 + creal(g));

	} else if (nt == 0 && np == 1) {

	    tau_t2 = tau/t2;
	    sqtau_t2 = tau_t2*tau_t2;

	    g0p = rec_pt*(g01 + dzi*(2.0*g02 + dzi*(3.0*g03 + 4.0*g04*dzi)));

	    r2p = rec_pt*(r21 + 2.0*r22*dzi);

	    g = r2p*(tau*clog((1.0 + tau_t2)/(1.0 - tau_t2))
	        + t2*(clog(1.0 - sqtau_t2) - sqtau_t2));

	    return (g0p + tt*creal(g));

	} else if (nt == 1 && np == 1) {

	    tau_t2 = tau/t2;

	    r2p = rec_pt*(r21 + 2.0*r22*dzi) ;

	    g = r2p*(clog((1.0 + tau_t2)/(1.0 - tau_t2)) - 2.0*tau_t2);

	    return (creal(g));

	} else if (nt == 2 && np == 0) {

	    r2 = r20 + dzi*(r21 + r22*dzi);

	    g = r1*(1.0/(t1 - tau) + 1.0/(t1 + tau) - 2.0/t1)
	        + r2*(1.0/(t2 - tau) + 1.0/(t2 + tau) - 2.0/t2);

	    return (rec_tt*creal(g));

	} else if (nt == 0 && np == 2) {

	    sqrec_pt = rec_pt*rec_pt;

	    tau_t2 = tau/t2;
	    sqtau_t2 = tau_t2*tau_t2;

	    g0pp = sqrec_pt*(2.0*g02 + dzi*(6.0*g03 + 12.0*g04*dzi));

	    r2pp = 2.0*r22*sqrec_pt;

	    g = r2pp*(tau*clog((1.0 + tau_t2)/(1.0 - tau_t2))
	        + t2*(clog(1.0 - sqtau_t2) - sqtau_t2));

	   return (g0pp + tt*creal(g));

	} else
	   return (GSW_INVALID_VALUE);
}
