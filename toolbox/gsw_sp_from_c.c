/*
!--------------------------------------------------------------------------
! Practical Salinity (SP), PSS-78
!--------------------------------------------------------------------------

!==========================================================================
function gsw_sp_from_c(c,t,p)
!==========================================================================

!  Calculates Practical Salinity, SP, from conductivity, C, primarily using
!  the PSS-78 algorithm.  Note that the PSS-78 algorithm for Practical
!  Salinity is only valid in the range 2 < SP < 42.  If the PSS-78
!  algorithm produces a Practical Salinity that is less than 2 then the
!  Practical Salinity is recalculated with a modified form of the Hill et
!  al. (1986) formula.  The modification of the Hill et al. (1986)
!  expression is to ensure that it is exactly consistent with PSS-78
!  at SP = 2.  Note that the input values of conductivity need to be in
!  units of mS/cm (not S/m).
!
! c      : conductivity                                     [ mS/cm ]
! t      : in-situ temperature [ITS-90]                     [deg C]
! p      : sea pressure                                     [dbar]
!
! sp     : Practical Salinity                               [unitless]
*/
double
gsw_sp_from_c(double c, double t, double p)
{
	GSW_TEOS10_CONSTANTS;
	GSW_SP_COEFFICIENTS;
	double	sp, t68, ft68, r, rt_lc, rp, rt, rtx,
		hill_ratio, x, sqrty, part1, part2,
		sp_hill_raw;

	t68	= t*1.00024e0;
	ft68	= (t68 - 15e0)/(1e0 + k*(t68 - 15e0));
    /*
     ! The dimensionless conductivity ratio, R, is the conductivity input, C,
     ! divided by the present estimate of C(SP=35, t_68=15, p=0) which is
     ! 42.9140 mS/cm (=4.29140 S/m), (Culkin and Smith, 1980).
    */

	r = c/gsw_c3515;	/* 0.023302418791070513 = 1./42.9140 */

	/*rt_lc corresponds to rt as defined in the UNESCO 44 (1983) routines.*/
	rt_lc	= c0 + (c1 + (c2 + (c3 + c4*t68)*t68)*t68)*t68;
	rp	= 1e0 + (p*(e1 + e2*p + e3*p*p))/(1e0 + d1*t68 + d2*t68*t68 +
		  (d3 + d4*t68)*r);
	rt	= r/(rp*rt_lc);

	if (rt < 0.0) {
	    return (GSW_INVALID_VALUE);
	}

	rtx	= sqrt(rt);

	sp	= a0 + (a1 + (a2 + (a3 + (a4 + a5*rtx)*rtx)*rtx)*rtx)*rtx + 
		  ft68*(b0 + (b1 + (b2 + (b3 +
			(b4 + b5*rtx)*rtx)*rtx)*rtx)*rtx);
    /*
     ! The following section of the code is designed for SP < 2 based on the
     ! Hill et al. (1986) algorithm.  This algorithm is adjusted so that it is
     ! exactly equal to the PSS-78 algorithm at SP = 2.
    */

	if (sp < 2) {
	    hill_ratio	= gsw_hill_ratio_at_sp2(t);
	    x		= 400e0*rt;
	    sqrty	= 10e0*rtx;
	    part1	= 1e0 + x*(1.5e0 + x);
	    part2	= 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
	    sp_hill_raw	= sp - a0/part1 - b0*ft68/part2;
	    sp		= hill_ratio*sp_hill_raw;
	}

    /* This line ensures that SP is non-negative. */
	if (sp < 0.0) {
	    sp	= GSW_INVALID_VALUE;
	}

	return (sp);
}
