/*
!==========================================================================
function gsw_pt0_from_t(sa,t,p)  
!==========================================================================
   
! Calculates potential temperature with reference pressure, p_ref = 0 dbar. 
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_pt0_from_t : potential temperature, p_ref = 0        [deg C]
*/
double
gsw_pt0_from_t(double sa, double t, double p)
{
	GSW_TEOS10_CONSTANTS;
	int	no_iter;
	double	pt0, pt0_old, dentropy, dentropy_dt;
	double	s1, true_entropy_part, pt0m;

	s1	= sa/gsw_ups;

	pt0	= t+p*( 8.65483913395442e-6  -
        	  s1 *  1.41636299744881e-6  -
		   p *  7.38286467135737e-9  +
		   t *(-8.38241357039698e-6  +
		  s1 *  2.83933368585534e-8  +
		   t *  1.77803965218656e-8  +
		   p *  1.71155619208233e-10));

	dentropy_dt	= gsw_cp0/((gsw_t0+pt0)*(1.0-0.05*(1.0-sa/gsw_sso)));

	true_entropy_part = gsw_entropy_part(sa,t,p);

	for (no_iter=1; no_iter <= 2; no_iter++) {
	    pt0_old	= pt0;
	    dentropy	= gsw_entropy_part_zerop(sa,pt0_old) -
			  true_entropy_part;
	    pt0		= pt0_old - dentropy/dentropy_dt;
	    pt0m	= 0.5*(pt0 + pt0_old);
	    dentropy_dt	= -gsw_gibbs_pt0_pt0(sa,pt0m);
	    pt0		= pt0_old - dentropy/dentropy_dt;
	}
	return (pt0);
}
