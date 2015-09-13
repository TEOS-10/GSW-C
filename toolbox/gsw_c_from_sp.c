/*
!==========================================================================
function gsw_c_from_sp(sp,t,p)
!==========================================================================

!  Calculates conductivity, C, from (SP,t,p) using PSS-78 in the range
!  2 < SP < 42.  If the input Practical Salinity is less than 2 then a
!  modified form of the Hill et al. (1986) fomula is used for Practical
!  Salinity.  The modification of the Hill et al. (1986) expression is to
!  ensure that it is exactly consistent with PSS-78 at SP = 2.
!
!  The conductivity ratio returned by this function is consistent with the
!  input value of Practical Salinity, SP, to 2x10^-14 psu over the full
!  range of input parameters (from pure fresh water up to SP = 42 psu).
!  This error of 2x10^-14 psu is machine precision at typical seawater
!  salinities.  This accuracy is achieved by having four different
!  polynomials for the starting value of Rtx (the square root of Rt) in
!  four different ranges of SP, and by using one and a half iterations of
!  a computationally efficient modified Newton-Raphson technique (McDougall
!  and Wotherspoon, 2012) to find the root of the equation.
!
!  Note that strictly speaking PSS-78 (Unesco, 1983) defines Practical
!  Salinity in terms of the conductivity ratio, R, without actually
!  specifying the value of C(35,15,0) (which we currently take to be
!  42.9140 mS/cm).
!
! sp     : Practical Salinity                               [unitless]
! t      : in-situ temperature [ITS-90]                     [deg C]
! p      : sea pressure                                     [dbar]
!
! c      : conductivity                                     [ mS/cm ]
*/
double
gsw_c_from_sp(double sp, double t, double p)
{
	GSW_TEOS10_CONSTANTS;
	GSW_SP_COEFFICIENTS;
	double	p0 = 4.577801212923119e-3,	p1 = 1.924049429136640e-1,
		p2 = 2.183871685127932e-5,	p3 = -7.292156330457999e-3,
		p4 = 1.568129536470258e-4,	p5 = -1.478995271680869e-6,
		p6 = 9.086442524716395e-4,	p7 = -1.949560839540487e-5,
		p8 = -3.223058111118377e-6,	p9 = 1.175871639741131e-7,
		p10 = -7.522895856600089e-5,	p11 = -2.254458513439107e-6,
		p12 = 6.179992190192848e-7,	p13 = 1.005054226996868e-8,
		p14 = -1.923745566122602e-9,	p15 = 2.259550611212616e-6,
		p16 = 1.631749165091437e-7,	p17 = -5.931857989915256e-9,
		p18 = -4.693392029005252e-9,	p19 = 2.571854839274148e-10,
		p20 = 4.198786822861038e-12,
		q0 = 5.540896868127855e-5,	q1 = 2.015419291097848e-1,
		q2 = -1.445310045430192e-5,	q3 = -1.567047628411722e-2,
		q4 = 2.464756294660119e-4,	q5 = -2.575458304732166e-7,
		q6 = 5.071449842454419e-3,	q7 = 9.081985795339206e-5,
		q8 = -3.635420818812898e-6,	q9 = 2.249490528450555e-8,
		q10 = -1.143810377431888e-3,	q11 = 2.066112484281530e-5,
		q12 = 7.482907137737503e-7,	q13 = 4.019321577844724e-8,
		q14 = -5.755568141370501e-10,	q15 = 1.120748754429459e-4,
		q16 = -2.420274029674485e-6,	q17 = -4.774829347564670e-8,
		q18 = -4.279037686797859e-9,	q19 = -2.045829202713288e-10,
		q20 = 5.025109163112005e-12,
		s0 = 3.432285006604888e-3,	s1 = 1.672940491817403e-1,
		s2 = 2.640304401023995e-5,	s3 = 1.082267090441036e-1,
		s4 = -6.296778883666940e-5,	s5 = -4.542775152303671e-7,
		s6 = -1.859711038699727e-1,	s7 = 7.659006320303959e-4,
		s8 = -4.794661268817618e-7,	s9 = 8.093368602891911e-9,
		s10 = 1.001140606840692e-1,	s11 = -1.038712945546608e-3,
		s12 = -6.227915160991074e-6,	s13 = 2.798564479737090e-8,
		s14 = -1.343623657549961e-10,	s15 = 1.024345179842964e-2,
		s16 = 4.981135430579384e-4,	s17 = 4.466087528793912e-6,
		s18 = 1.960872795577774e-8,	s19 = -2.723159418888634e-10,
		s20 = 1.122200786423241e-12,
		u0 = 5.180529787390576e-3,	u1 = 1.052097167201052e-3,
		u2 = 3.666193708310848e-5,	u3 = 7.112223828976632e0,
		u4 = -3.631366777096209e-4,	u5 = -7.336295318742821e-7,
		u6 = -1.576886793288888e+2,	u7 = -1.840239113483083e-3,
		u8 = 8.624279120240952e-6,	u9 = 1.233529799729501e-8,
		u10 = 1.826482800939545e+3,	u11 = 1.633903983457674e-1,
		u12 = -9.201096427222349e-5,	u13 = -9.187900959754842e-8,
		u14 = -1.442010369809705e-10,	u15 = -8.542357182595853e+3,
		u16 = -1.408635241899082e0,	u17 = 1.660164829963661e-4,
		u18 = 6.797409608973845e-7,	u19 = 3.345074990451475e-10,
		u20 = 8.285687652694768e-13;

	double	t68, ft68, x, rtx=0.0, dsp_drtx, sqrty,
		part1, part2, hill_ratio, sp_est,
		rtx_old, rt, aa, bb, cc, dd, ee, ra,r, rt_lc, rtxm,
		sp_hill_raw;

	t68	= t*1.00024e0;
	ft68	= (t68 - 15e0)/(1e0 + k*(t68 - 15e0));

	x	= sqrt(sp);

    /*
     |--------------------------------------------------------------------------
     ! Finding the starting value of Rtx, the square root of Rt, using four
     ! different polynomials of SP and t68.
     !--------------------------------------------------------------------------
    */

	if (sp >= 9.0) {
	    rtx	= p0 + x*(p1 + p4*t68 + x*(p3 + p7*t68 + x*(p6 
		  + p11*t68 + x*(p10 + p16*t68 + x*p15))))
		  + t68*(p2+ t68*(p5 + x*x*(p12 + x*p17) + p8*x
		  + t68*(p9 + x*(p13 + x*p18)+ t68*(p14 + p19*x + p20*t68))));
	} else if (sp >= 0.25 && sp < 9.0) {
	    rtx	= q0 + x*(q1 + q4*t68 + x*(q3 + q7*t68 + x*(q6
		  + q11*t68 + x*(q10 + q16*t68 + x*q15)))) 
		  + t68*(q2+ t68*(q5 + x*x*(q12 + x*q17) + q8*x 
		  + t68*(q9 + x*(q13 + x*q18)+ t68*(q14 + q19*x + q20*t68))));
	} else if (sp >= 0.003 && sp < 0.25) {
	    rtx	=  s0 + x*(s1 + s4*t68 + x*(s3 + s7*t68 + x*(s6
		  + s11*t68 + x*(s10 + s16*t68 + x*s15)))) 
		  + t68*(s2+ t68*(s5 + x*x*(s12 + x*s17) + s8*x 
		  + t68*(s9 + x*(s13 + x*s18)+ t68*(s14 + s19*x + s20*t68))));
	} else if (sp < 0.003) {
	    rtx	=  u0 + x*(u1 + u4*t68 + x*(u3 + u7*t68 + x*(u6
		  + u11*t68 + x*(u10 + u16*t68 + x*u15)))) 
		  + t68*(u2+ t68*(u5 + x*x*(u12 + x*u17) + u8*x 
		  + t68*(u9 + x*(u13 + x*u18)+ t68*(u14 + u19*x + u20*t68))));
	}

    /*
     !--------------------------------------------------------------------------
     ! Finding the starting value of dSP_dRtx, the derivative of SP with respect
     ! to Rtx.
     !--------------------------------------------------------------------------
    */
	dsp_drtx	=  a1 + (2e0*a2 + (3e0*a3 +
				(4e0*a4 + 5e0*a5*rtx)*rtx)*rtx)*rtx 
    			  + ft68*(b1 + (2e0*b2 + (3e0*b3 + (4e0*b4 +
				5e0*b5*rtx)*rtx)*rtx)*rtx);

	if (sp < 2.0) {
	    x		= 400e0*(rtx*rtx);
	    sqrty 	= 10.0*rtx;
	    part1	= 1e0 + x*(1.5e0 + x);
	    part2	= 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
	    hill_ratio	= gsw_hill_ratio_at_sp2(t);
	    dsp_drtx	= dsp_drtx
			  + a0*800e0*rtx*(1.5e0 + 2e0*x)/(part1*part1)
			  + b0*ft68*(10e0 + sqrty*(20e0 + 30e0*sqrty))/
				(part2*part2);
	    dsp_drtx	= hill_ratio*dsp_drtx;
	}

    /*
     !--------------------------------------------------------------------------
     ! One iteration through the modified Newton-Raphson method (McDougall and
     ! Wotherspoon, 2012) achieves an error in Practical Salinity of about
     ! 10^-12 for all combinations of the inputs.  One and a half iterations of
     ! the modified Newton-Raphson method achevies a maximum error in terms of
     ! Practical Salinity of better than 2x10^-14 everywhere.
     !
     ! We recommend one and a half iterations of the modified Newton-Raphson
     ! method.
     !
     ! Begin the modified Newton-Raphson method.
     !--------------------------------------------------------------------------
    */
	sp_est	= a0 + (a1 + (a2 + (a3 + (a4 + a5*rtx)*rtx)*rtx)*rtx)*rtx
		+ ft68*(b0 + (b1 + (b2+ (b3 + (b4 + b5*rtx)*rtx)*rtx)*rtx)*rtx);
	if (sp_est <  2.0) {
	    x		= 400e0*(rtx*rtx);
	    sqrty	= 10e0*rtx;
	    part1	= 1e0 + x*(1.5e0 + x);
	    part2	= 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
	    sp_hill_raw	= sp_est - a0/part1 - b0*ft68/part2;
	    hill_ratio	= gsw_hill_ratio_at_sp2(t);
	    sp_est	= hill_ratio*sp_hill_raw;
	}

	rtx_old	= rtx;
	rtx	= rtx_old - (sp_est - sp)/dsp_drtx;

	rtxm	= 0.5e0*(rtx + rtx_old); /*This mean value of Rtx, Rtxm, is the
		  value of Rtx at which the derivative dSP_dRtx is evaluated.*/

	dsp_drtx=  a1 + (2e0*a2 + (3e0*a3 + (4e0*a4 +
				5e0*a5*rtxm)*rtxm)*rtxm)*rtxm
		   + ft68*(b1 + (2e0*b2 + (3e0*b3 + (4e0*b4 +
				5e0*b5*rtxm)*rtxm)*rtxm)*rtxm);
	if (sp_est <  2.0) {
	    x	= 400e0*(rtxm*rtxm);
	    sqrty	= 10e0*rtxm;
	    part1	= 1e0 + x*(1.5e0 + x);
	    part2	= 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
	    dsp_drtx	= dsp_drtx 
			  + a0*800e0*rtxm*(1.5e0 + 2e0*x)/(part1*part1)
			  + b0*ft68*(10e0 + sqrty*(20e0 + 30e0*sqrty))/
				(part2*part2);
	    hill_ratio	= gsw_hill_ratio_at_sp2(t);
	    dsp_drtx	= hill_ratio*dsp_drtx;
	}

    /*
     !--------------------------------------------------------------------------
     ! The line below is where Rtx is updated at the end of the one full
     ! iteration of the modified Newton-Raphson technique.
     !--------------------------------------------------------------------------
    */
	rtx	= rtx_old - (sp_est - sp)/dsp_drtx;
    /*
     !--------------------------------------------------------------------------
     ! Now we do another half iteration of the modified Newton-Raphson
     ! technique, making a total of one and a half modified N-R iterations.
     !--------------------------------------------------------------------------
    */
	sp_est	= a0 + (a1 + (a2 + (a3 + (a4 + a5*rtx)*rtx)*rtx)*rtx)*rtx 
		+ ft68*(b0 + (b1 + (b2+ (b3 + (b4 + b5*rtx)*rtx)*rtx)*rtx)*rtx);
	if (sp_est <  2.0) {
	    x		= 400e0*(rtx*rtx);
	    sqrty	= 10e0*rtx;
	    part1	= 1e0 + x*(1.5e0 + x);
	    part2	= 1e0 + sqrty*(1e0 + sqrty*(1e0 + sqrty));
	    sp_hill_raw	= sp_est - a0/part1 - b0*ft68/part2;
	    hill_ratio	= gsw_hill_ratio_at_sp2(t);
	    sp_est	= hill_ratio*sp_hill_raw;
	}
	rtx	= rtx - (sp_est - sp)/dsp_drtx;

    /*
     !--------------------------------------------------------------------------
     ! Now go from Rtx to Rt and then to the conductivity ratio R at pressure p.
     !--------------------------------------------------------------------------
    */
	rt	= rtx*rtx;

	aa	= d3 + d4*t68;
	bb	= 1e0 + t68*(d1 + d2*t68);
	cc	= p*(e1 + p*(e2 + e3*p));
    /* rt_lc (i.e. rt_lower_case) corresponds to rt as defined in
       the UNESCO 44 (1983) routines. */
	rt_lc	= c0 + (c1 + (c2 + (c3 + c4*t68)*t68)*t68)*t68;

	dd	= bb - aa*rt_lc*rt;
	ee	= rt_lc*rt*aa*(bb + cc);
	ra	= sqrt(dd*dd + 4e0*ee) - dd;
	r	= 0.5e0*ra/aa;

    /*
     ! The dimensionless conductivity ratio, R, is the conductivity input, C,
     ! divided by the present estimate of C(SP=35, t_68=15, p=0) which is
     ! 42.9140 mS/cm (=4.29140 S/m^).
    */
	return (gsw_c3515*r);
}
