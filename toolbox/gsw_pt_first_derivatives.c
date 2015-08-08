/*
!==========================================================================
elemental subroutine gsw_pt_first_derivatives (sa, ct, pt_sa, pt_ct)
! =========================================================================
!
!  Calculates the following two partial derivatives of potential temperature 
!  (the regular potential temperature whose reference sea pressure is 0 dbar) 
!  (1) pt_SA, the derivative with respect to Absolute Salinity at 
!       constant Conservative Temperature, and
!  (2) pt_CT, the derivative with respect to Conservative Temperature at 
!       constant Absolute Salinity. 
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  pt_SA =  The derivative of potential temperature with respect to 
!           Absolute Salinity at constant Conservative Temperature. 
!                                                               [ K/(g/kg)]
!  pt_CT =  The derivative of potential temperature with respect to 
!           Conservative Temperature at constant Absolute Salinity.
!           pt_CT is dimensionless.                            [ unitless ]
!--------------------------------------------------------------------------
*/
void
gsw_pt_first_derivatives (double sa, double ct, double *pt_sa, double *pt_ct)
{
	GSW_TEOS10_CONSTANTS;
	double	abs_pt, ct_pt, ct_sa, pt, pr0 = 0.0;
	int	n0=0, n1=1, n2=2;

	pt = gsw_pt_from_ct(sa,ct);
	abs_pt = (gsw_t0 + pt);

	ct_pt = -(abs_pt*gsw_gibbs(n0,n2,n0,sa,pt,pr0))/gsw_cp0;

	if (pt_sa != NULL) {

	    ct_sa = (gsw_gibbs(n1,n0,n0,sa,pt,pr0) -
                abs_pt*gsw_gibbs(n1,n1,n0,sa,pt,pr0))/gsw_cp0;

	    *pt_sa = -ct_sa/ct_pt;

	}

	if (pt_ct != NULL)
	    *pt_ct = 1.0/ct_pt;
}
