/*
!==========================================================================
elemental subroutine gsw_entropy_second_derivatives (sa, ct, eta_sa_sa, &
                                                     eta_sa_ct, eta_ct_ct)
! =========================================================================
!
!  Calculates the following three second-order partial derivatives of 
!  specific entropy (eta) 
!   (1) eta_SA_SA, the second derivative with respect to Absolute  
!       Salinity at constant Conservative Temperature, and
!   (2) eta_SA_CT, the derivative with respect to Absolute Salinity and 
!       Conservative Temperature. 
!   (3) eta_CT_CT, the second derivative with respect to Conservative 
!       Temperature at constant Absolute Salinity. 
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  eta_SA_SA =  The second derivative of specific entropy with respect 
!               to Absolute Salinity (in units of g kg^-1) at constant 
!               Conservative Temperature.    
!               eta_SA_SA has units of:                 [ J/(kg K(g/kg)^2)]
!  eta_SA_CT =  The second derivative of specific entropy with respect 
!               to Conservative Temperature at constant Absolute 
!               Salinity. eta_SA_CT has units of:     [ J/(kg (g/kg) K^2) ]
!  eta_CT_CT =  The second derivative of specific entropy with respect 
!               to Conservative Temperature at constant Absolute
!               Salinity.  eta_CT_CT has units of:           [ J/(kg K^3) ]
!--------------------------------------------------------------------------
*/
void
gsw_entropy_second_derivatives(double sa, double ct,
	double *eta_sa_sa, double *eta_sa_ct, double *eta_ct_ct)
{
	GSW_TEOS10_CONSTANTS;
	double	abs_pt, ct_pt, ct_sa, pt, ct_ct, pr0 = 0.0;
	int	n0=0, n1=1, n2=2;

	pt = gsw_pt_from_ct(sa,ct);
	abs_pt = gsw_t0 + pt;
	
	ct_pt = -(abs_pt*gsw_gibbs(n0,n2,n0,sa,pt,pr0))/gsw_cp0;
	
	ct_ct = -gsw_cp0/(ct_pt*abs_pt*abs_pt);
	
	if ((eta_sa_ct != NULL) || (eta_sa_sa != NULL)) {
	
	    ct_sa = (gsw_gibbs(n1,n0,n0,sa,pt,pr0) -
	               (abs_pt*gsw_gibbs(n1,n1,n0,sa,pt,pr0)))/gsw_cp0;
	
	    if (eta_sa_ct != NULL) *eta_sa_ct = -ct_sa*ct_ct;
	
	    if (eta_sa_sa != NULL)
		*eta_sa_sa = -gsw_gibbs(n2,n0,n0,sa,pt,pr0)/abs_pt +
	                                     ct_sa*ct_sa*ct_ct;
	}
	
	if (eta_ct_ct != NULL) *eta_ct_ct = ct_ct;
}
