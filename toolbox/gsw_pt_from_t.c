/*
!==========================================================================
function gsw_pt_from_t(sa,t,p,p_ref)  
!==========================================================================
   
! Calculates potential temperature of seawater from in-situ temperature 
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! p_ref  : reference sea pressure                          [dbar]
!
! gsw_pt_from_t : potential temperature                    [deg C]
*/
double
gsw_pt_from_t(double sa, double t, double p, double p_ref)
{
        GSW_TEOS10_CONSTANTS;
	int	n0=0, n2=2, no_iter;
	double	s1, pt, ptm, pt_old, dentropy, dentropy_dt,
		true_entropy_part;

	s1	= sa/gsw_ups;
	pt	= t+(p-p_ref)*( 8.65483913395442e-6  -
			  s1 *  1.41636299744881e-6  -
		   (p+p_ref) *  7.38286467135737e-9  +
			  t  *(-8.38241357039698e-6  +
			  s1 *  2.83933368585534e-8  +
			  t  *  1.77803965218656e-8  +
		   (p+p_ref) *  1.71155619208233e-10));

	dentropy_dt	= gsw_cp0/((gsw_t0 + pt)*(1.0-0.05*(1.0 - sa/gsw_sso)));
	true_entropy_part	= gsw_entropy_part(sa,t,p);
	for (no_iter=1; no_iter <= 2; no_iter++) {
	    pt_old	= pt;
	    dentropy	= gsw_entropy_part(sa,pt_old,p_ref) - true_entropy_part;
	    pt		= pt_old - dentropy/dentropy_dt;
	    ptm		= 0.5*(pt + pt_old);
	    dentropy_dt	= -gsw_gibbs(n0,n2,n0,sa,ptm,p_ref);
	    pt		= pt_old - dentropy/dentropy_dt;
	}
	return (pt);
}
