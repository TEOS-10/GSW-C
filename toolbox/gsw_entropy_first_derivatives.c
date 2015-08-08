/*
!==========================================================================
elemental subroutine gsw_entropy_first_derivatives (sa, ct, eta_sa, eta_ct)
! =========================================================================
!
!  Calculates the following two partial derivatives of specific entropy
!  (eta) 
!   (1) eta_SA, the derivative with respect to Absolute Salinity at 
!       constant Conservative Temperature, and
!   (2) eta_CT, the derivative with respect to Conservative Temperature at 
!       constant Absolute Salinity. 
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  eta_SA =  The derivative of specific entropy with respect to 
!            Absolute Salinity (in units of g kg^-1) at constant  
!            Conservative Temperature.  
!            eta_SA has units of:         [ J/(kg K(g/kg))]  or [ J/(g K) ]
!  eta_CT =  The derivative of specific entropy with respect to 
!            Conservative Temperature at constant Absolute Salinity.
!            eta_CT has units of:                            [ J/(kg K^2) ]
!--------------------------------------------------------------------------
*/
void
gsw_entropy_first_derivatives(double sa, double ct, double *eta_sa,
	double *eta_ct)
{
	GSW_TEOS10_CONSTANTS;
	double	pt, pr0 = 0.0;
	int	n0=0, n1=1;

	pt = gsw_pt_from_ct(sa,ct);

	if (eta_sa != NULL)
	    *eta_sa = -(gsw_gibbs(n1,n0,n0,sa,pt,pr0))/(gsw_t0 + pt);

	if (eta_ct != NULL)
	    *eta_ct = gsw_cp0/(gsw_t0 + pt);
}
