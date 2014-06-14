/*
**  $Id: gsw_oceanographic_toolbox.c,v d7a5468a0b8c 2014/06/14 02:35:01 fmd $
**  $Version: 3.03.0 $
**
**  This is a translation of the original f90 source code into C
**  by the Shipboard Technical Support Computing Resources group
**  at Scripps Institution of Oceanography -- sts-cr@sio.ucsd.edu.
**  The original notices follow.
**

!==========================================================================
! Gibbs SeaWater (GSW) Oceanographic Toolbox of TEOS-10 version 3.03 (Fortran)
!==========================================================================
!
! This is a subset of functions contained in the Gibbs SeaWater (GSW) 
! Oceanographic Toolbox of TEOS-10 (version 3.03).
! 
! Practical Salinity (SP), PSS-78
! gsw_sp_from_c           - Practical Salinity from conductivity
! gsw_c_from_sp           - conductivity from Practical Salinity
! gsw_sp_from_sk          - Practical Salinity from Knudsen Salinity
! 
! salinity and temperature conversions
! gsw_sa_from_sp          - Absolute Salinity from Practical Salinity
! gsw_sstar_from_sp       - Preformed Salinity from Practical Salinity
! gsw_ct_from_t           - Conservative Temperature from in-situ temperature
!
! gsw_deltasa_from_sp     - Absolute Salinity Anomaly from Practical Salinity
! gsw_sr_from_sp          - Reference Salinity from Practical Salinity
! gsw_sp_from_sr          - Practical Salinity from Reference Salinity
! gsw_sp_from_sa          - Practical Salinity from Absolute Salinity
! gsw_sstar_from_sa       - Preformed Salinity from Absolute Salinity
! gsw_sp_from_sstar       - Practical Salinity from Preformed Salinity
! gsw_sa_from_sstar       - Absolute Salinity from Preformed Salinity
! gsw_ct_from_pt          - Conservative Temperature from potential temperature
! gsw_pt_from_ct          - potential temperature from Conservative Temperature
! gsw_pt0_from_t          - potential temperature with reference pressure of 0 dbar
! gsw_pt_from_t           - potential temperature 
! gsw_z_from_p            - height from pressure
! gsw_entropy_from_t      - entropy from in-situ temperature
! gsw_adiabatic_lapse_rate_from_ct - adiabatic lapse rate from CT
!
! density and enthalpy, based on the 48-term expression for density
! gsw_rho                 - in-situ density and potential density
! gsw_alpha               - thermal expansion coefficient with respect to CT
! gsw_beta                - saline contraction coefficient at constant CT
! gsw_alpha_on_beta       - alpha divided by beta
! gsw_rho_first_derivatives  - first derivatives of density
! gsw_specvol             - specific volume
! gsw_specvol_anom        - specific volume anomaly
! gsw_sigma0              - sigma0 with reference pressure of 0 dbar
! gsw_sigma1              - sigma1 with reference pressure of 1000 dbar
! gsw_sigma2              - sigma2 with reference pressure of 2000 dbar
! gsw_sigma3              - sigma3 with reference pressure of 3000 dbar
! gsw_sigma4              - sigma4 with reference pressure of 4000 dbar
! gsw_sound_speed         - sound speed
! gsw_kappa               - isentropic compressibility
! gsw_cabbeling           - cabbeling coefficient
! gsw_thermobaric         - thermobaric coefficient
! gsw_internal_energy     - internal energy
! gsw_enthalpy            - enthalpy
! gsw_dynamic_enthalpy    - dynamic enthalpy
! gsw_sa_from_rho         - Absolute Salinity from density
!
! water column properties, based on the 48-term expression for density
! gsw_nsquared            - buoyancy (Brunt-Vaisala) frequency squared (N^2)
! gsw_turner_rsubrho      - Turner angle & Rsubrho
! gsw_ipv_vs_fnsquared_ratio  - ratio of the vertical gradient of potential density
!                               (with reference pressure, p_ref), to the vertical gradient
!                               of locally-referenced potential density
!
! freezing temperatures
! gsw_ct_freezing         - Conservative Temperature freezing temperature of seawater
! gsw_t_freezing          - in-situ temperature freezing temperature of seawater
!
! isobaric melting enthalpy and isobaric evaporation enthalpy
! gsw_latentheat_melting  - latent heat of melting
! gsw_latentheat_evap_ct  - latent heat of evaporation
! gsw_latentheat_evap_t   - latent heat of evaporation
!
! planet Earth properties
! gsw_grav                - gravitational acceleration
!
! basic thermodynamic properties in terms of in-situ t, based on the exact Gibbs function
! gsw_rho_t_exact         - in-situ density
! gsw_pot_rho_t_exact     - potential density
! gsw_alpha_wrt_t_exact   - thermal expansion coefficient with respect to in-situ temperature
! gsw_beta_const_t_exact  - saline contraction coefficient at constant in-situ temperature
! gsw_specvol_t_exact     - specific volume
! gsw_sound_speed_t_exact - sound speed
! gsw_kappa_t_exact       - isentropic compressibility
! gsw_enthalpy_t_exact    - enthalpy
! gsw_cp_t_exact          - isobaric heat capacity
!
! Library functions of the GSW toolbox
! gsw_gibbs               - the TEOS-10 Gibbs function and its derivatives
! gsw_saar                - Absolute Salinity Anomaly Ratio (excluding the Baltic Sea)
! gsw_deltasa_atlas       - Absolute Salinity Anomaly atlas value (excluding the Baltic Sea)
! gsw_fdelta              - ratio of Absolute to Preformed Salinity, minus 1
! gsw_sa_from_sp_baltic   - Absolute Salinity Anomaly from Practical Salinity in the Baltic Sea
! gsw_sp_from_sa_baltic   - Practical Salinity from Absolute Salinity in the Baltic Sea
! gsw_entropy_part        - entropy minus the terms that are a function of only SA
! gsw_entropy_part_zerop  - entropy_part evaluated at 0 dbar
! gsw_gibbs_pt0_pt0       - gibbs(0,2,0,SA,t,0)
! gsw_specvol_sso_0_p     - specvol at (SSO,CT=0,p)
! gsw_enthalpy_sso_0_p    - enthalpy at (SSO,CT=0,p)
! gsw_hill_ratio_at_sp2   - Hill ratio at Practical Salinity of 2
!
!
! Version 1.0 written by David Jackett
! Modified by Paul Barker (version 3.03)
!
! For help with this Oceanographic Toolbox email:- help_gsw@teos-10.org
!
! This software is available from http://www.teos-10.org
!
!==========================================================================
**
*/
#include <gswteos-10.h>

/*
!--------------------------------------------------------------------------
! Practical Salinity (SP), PSS-78
!--------------------------------------------------------------------------
*/
/*
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
	double	a0 = 0.0080e0, a1 = -0.1692e0, a2 = 25.3851e0,
		a3 = 14.0941e0, a4 = -7.0261e0, a5 = 2.7081e0,
		b0 = 0.0005e0, b1 = -0.0056e0, b2 = -0.0066e0,
		b3 = -0.0375e0, b4 = 0.0636e0, b5 = -0.0144e0,
		c0 = 0.6766097e0, c1 = 2.00564e-2,
		c2 = 1.104259e-4, c3 = -6.9698e-7, c4 = 1.0031e-9,
		d1 = 3.426e-2, d2 = 4.464e-4, d3 =  4.215e-1,
		d4 = -3.107e-3, e1 = 2.070e-5, e2 = -6.370e-10,
		e3 = 3.989e-15, k  = 0.0162;

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

	r = 0.023302418791070513e0*c;	/* 0.023302418791070513 = 1./42.9140 */

	/*rt_lc corresponds to rt as defined in the UNESCO 44 (1983) routines.*/
	rt_lc	= c0 + (c1 + (c2 + (c3 + c4*t68)*t68)*t68)*t68;
	rp	= 1e0 + (p*(e1 + e2*p + e3*p*p))/(1e0 + d1*t68 + d2*t68*t68 +
		  (d3 + d4*t68)*r);
	rt	= r/(rp*rt_lc);

	if (rt < 0.0) {
	    rt = GSW_INVALID_VALUE;
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

/*
!--------------------------------------------------------------------------

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
	double	a0 = 0.0080e0, a1 = -0.1692e0, a2 = 25.3851e0,
		a3 = 14.0941e0, a4 = -7.0261e0, a5 = 2.7081e0,
		b0 = 0.0005e0, b1 = -0.0056e0, b2 = -0.0066e0,
		b3 = -0.0375e0, b4 = 0.0636e0, b5 = -0.0144e0,
		c0 = 0.6766097e0, c1 = 2.00564e-2, c2 = 1.104259e-4,
		c3 = -6.9698e-7, c4 = 1.0031e-9, d1 = 3.426e-2,
		d2 = 4.464e-4, d3 = 4.215e-1, d4 = -3.107e-3,
		e1 = 2.070e-5, e2 = -6.370e-10, e3 = 3.989e-15,
		p0 = 4.577801212923119e-3, p1 = 1.924049429136640e-1,
		p2 = 2.183871685127932e-5, p3 = -7.292156330457999e-3,
		p4 = 1.568129536470258e-4, p5 = -1.478995271680869e-6,
		p6 = 9.086442524716395e-4, p7 = -1.949560839540487e-5,
		p8 = -3.223058111118377e-6, p9 = 1.175871639741131e-7,
		p10 = -7.522895856600089e-5,
		p11 = -2.254458513439107e-6,
		p12 = 6.179992190192848e-7,
		p13 = 1.005054226996868e-8,
		p14 = -1.923745566122602e-9,
		p15 = 2.259550611212616e-6,
		p16 = 1.631749165091437e-7,
		p17 = -5.931857989915256e-9,
		p18 = -4.693392029005252e-9,
		p19 = 2.571854839274148e-10,
		p20 = 4.198786822861038e-12,
		q0 = 5.540896868127855e-5, q1 = 2.015419291097848e-1,
		q2 = -1.445310045430192e-5,
		q3 = -1.567047628411722e-2,
		q4 = 2.464756294660119e-4,
		q5 = -2.575458304732166e-7,
		q6 = 5.071449842454419e-3,
		q7 = 9.081985795339206e-5,
		q8 = -3.635420818812898e-6,
		q9 = 2.249490528450555e-8,
		q10 = -1.143810377431888e-3,
		q11 = 2.066112484281530e-5,
		q12 = 7.482907137737503e-7,
		q13 = 4.019321577844724e-8,
		q14 = -5.755568141370501e-10,
		q15 = 1.120748754429459e-4,
		q16 = -2.420274029674485e-6,
		q17 = -4.774829347564670e-8,
		q18 = -4.279037686797859e-9,
		q19 = -2.045829202713288e-10,
		q20 = 5.025109163112005e-12,
		s0 = 3.432285006604888e-3, s1 = 1.672940491817403e-1,
		s2 = 2.640304401023995e-5, s3 = 1.082267090441036e-1,
		s4 = -6.296778883666940e-5, s5 = -4.542775152303671e-7,
		s6 = -1.859711038699727e-1, s7 = 7.659006320303959e-4,
		s8 = -4.794661268817618e-7, s9 = 8.093368602891911e-9,
		s10 = 1.001140606840692e-1,
		s11 = -1.038712945546608e-3,
		s12 = -6.227915160991074e-6,
		s13 = 2.798564479737090e-8,
		s14 = -1.343623657549961e-10,
		s15 = 1.024345179842964e-2,
		s16 = 4.981135430579384e-4,
		s17 = 4.466087528793912e-6,
		s18 = 1.960872795577774e-8,
		s19 = -2.723159418888634e-10,
		s20 = 1.122200786423241e-12,
		u0 = 5.180529787390576e-3, u1 = 1.052097167201052e-3,
		u2 = 3.666193708310848e-5, u3 = 7.112223828976632e0,
		u4 = -3.631366777096209e-4, u5 = -7.336295318742821e-7,
		u6 = -1.576886793288888e+2, u7 = -1.840239113483083e-3,
		u8 = 8.624279120240952e-6, u9 = 1.233529799729501e-8,
		u10 = 1.826482800939545e+3,
		u11 = 1.633903983457674e-1,
		u12 = -9.201096427222349e-5,
		u13 = -9.187900959754842e-8,
		u14 = -1.442010369809705e-10,
		u15 = -8.542357182595853e+3,
		u16 = -1.408635241899082e0,
		u17 = 1.660164829963661e-4,
		u18 = 6.797409608973845e-7,
		u19 = 3.345074990451475e-10,
		u20 = 8.285687652694768e-13, k = 0.0162e0;

	double	t68, ft68, x, rtx, dsp_drtx, sqrty,
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

	if (sp > 9.0) {
	    rtx	= p0 + x*(p1 + p4*t68 + x*(p3 + p7*t68 + x*(p6 
		  + p11*t68 + x*(p10 + p16*t68 + x*p15))))
		  + t68*(p2+ t68*(p5 + x*x*(p12 + x*p17) + p8*x
		  + t68*(p9 + x*(p13 + x*p18)+ t68*(p14 + p19*x + p20*t68))));
	}

	if (sp >= 0.25 && sp < 9.0) {
	    rtx	= q0 + x*(q1 + q4*t68 + x*(q3 + q7*t68 + x*(q6
		  + q11*t68 + x*(q10 + q16*t68 + x*q15)))) 
		  + t68*(q2+ t68*(q5 + x*x*(q12 + x*q17) + q8*x 
		  + t68*(q9 + x*(q13 + x*q18)+ t68*(q14 + q19*x + q20*t68))));
	}

	if (sp >= 0.003 && sp < 0.25) {
	    rtx	=  s0 + x*(s1 + s4*t68 + x*(s3 + s7*t68 + x*(s6
		  + s11*t68 + x*(s10 + s16*t68 + x*s15)))) 
		  + t68*(s2+ t68*(s5 + x*x*(s12 + x*s17) + s8*x 
		  + t68*(s9 + x*(s13 + x*s18)+ t68*(s14 + s19*x + s20*t68))));
	}

	if (sp < 0.003) {
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
	return (42.9140e0*r);
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_sp_from_sk(sk)
!==========================================================================

! Calculates Practical Salinity, SP, from SK
!
!  SK    : Knudsen Salinity                        [parts per thousand, ppt]
!
! gsw_sp_from_sk  : Practical Salinity                              [unitless]
*/
double
gsw_sp_from_sk(double sk)
{
	double	gsw_sp_from_sk_value;


	gsw_sp_from_sk_value = (sk - 0.03e0)*(1.80655e0/1.805e0);

	/*! This line ensures that SP is non-negative.*/
	if (gsw_sp_from_sk_value < 0e0)
            gsw_sp_from_sk_value = GSW_INVALID_VALUE;

	return (gsw_sp_from_sk_value);
}

/*
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! salinity and temperature conversions
!--------------------------------------------------------------------------
*/
/*
!==========================================================================
function gsw_sa_from_sp(sp,p,lon,lat)       
!==========================================================================

! Calculates Absolute Salinity, SA, from Practical Salinity, SP
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]
!
! gsw_sa_from_sp   : Absolute Salinity                     [g/kg]
*/
double
gsw_sa_from_sp(double sp, double p, double lon, double lat)
{
	double	saar, gsw_sa_baltic;

	saar		= gsw_saar(p,lon,lat);
	if (saar == GSW_INVALID_VALUE)
	    return (saar);
	gsw_sa_baltic	= gsw_sa_from_sp_baltic(sp,lon,lat);
	if (gsw_sa_baltic < 1e10)
	    return (gsw_sa_baltic);
	return ((35.16504e0/35.e0)*sp*(1.e0 + saar));
}

/*
!==========================================================================
function gsw_sstar_from_sp(sp,p,lon,lat)  
!==========================================================================

! Calculates Preformed Salinity, Sstar, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sstar_from_sp  : Preformed Salinity                  [g/kg]
*/
double
gsw_sstar_from_sp(double sp, double p, double lon, double lat)
{
	double	saar, sstar_baltic;

	saar		= gsw_saar(p,lon,lat);
	if (saar == GSW_INVALID_VALUE)
	    return (saar);
    /*
	!In the Baltic Sea, Sstar = SA.
    */
	sstar_baltic	= gsw_sa_from_sp_baltic(sp,lon,lat);
	if (sstar_baltic < 1e10)
	    return (sstar_baltic);
	return ((35.16504e0/35.e0)*sp*(1 - 0.35e0*saar));
}

/*
!==========================================================================
function gsw_ct_from_t(sa,t,p)  
!==========================================================================
   
! Calculates Conservative Temperature from in-situ temperature
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_ct_from_t : Conservative Temperature                 [deg C]
*/
double
gsw_ct_from_t(double sa, double t, double p)
{
	double	pt0;

	pt0	= gsw_pt0_from_t(sa,t,p);
	return (gsw_ct_from_pt(sa,pt0));
}

/*
!==========================================================================
function gsw_deltasa_from_sp(sp,p,lon,lat)  
!==========================================================================

! Calculates Absolute Salinity Anomaly, deltaSA, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_deltasa_from_sp : Absolute Salinty Anomaly           [g/kg]
*/
double
gsw_deltasa_from_sp(double sp, double p, double lon, double lat)
{
	double	res;

	res	= gsw_sa_from_sp(sp,p,lon,lat) - gsw_sr_from_sp(sp);
	if (res > 1e10)
	    res	= GSW_INVALID_VALUE;
	return (res);
}

/*
!==========================================================================
function gsw_sr_from_sp(sp)  
!==========================================================================

! Calculates Reference Salinity, SR, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
!
! gsw_sr_from_sp : Reference Salinity                      [g/kg]
*/
double
gsw_sr_from_sp(double sp)
{
	double	res;

	res	= 1.004715428571429e0*sp;
	if (res >= 1.e10)
	    res	= 9.e15;
	return (res);
}

/*
!==========================================================================
function gsw_sp_from_sr(sr)  
!==========================================================================

! Calculates Practical Salinity, sp, from Reference Salinity, sr. 
!
! sr     : Reference Salinity                              [g/kg]
!
! gsw_sp_from_sr  : Practical Salinity                     [unitless]
*/
double
gsw_sp_from_sr(double sr)
{
	double	res;

	res	= 0.995306702338459e0*sr;
	if (res > 1e10)
	    res	= GSW_INVALID_VALUE;
	return (res);
}

/*
!==========================================================================
function gsw_sp_from_sa(sa,p,lon,lat)  
!==========================================================================

! Calculates Practical salinity, sp, from Absolute salinity, sa  
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]
!
! gsw_sp_from_sa      : Practical Salinity                 [unitless]
*/
double
gsw_sp_from_sa(double sa, double p, double lon, double lat)
{
	double	saar, gsw_sp_baltic;

	saar	= gsw_saar(p,lon,lat);
	if (saar == GSW_INVALID_VALUE)
	    return (saar);
	gsw_sp_baltic	= gsw_sp_from_sa_baltic(sa,lon,lat);
	if (gsw_sp_baltic < 1e10)
	    return (gsw_sp_baltic);
	return ((35.e0/35.16504e0)*sa/(1e0 + saar));
}

/*
!==========================================================================
function gsw_sstar_from_sa(sa,p,lon,lat)  
!==========================================================================

! Calculates Preformed Salinity, Sstar, from Absolute Salinity, SA. 
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sstar_from_sa : Preformed Salinity                   [g/kg]
*/
double
gsw_sstar_from_sa(double sa, double p, double lon, double lat)
{
	double	saar;

	saar	= gsw_saar(p,lon,lat);
    /*
	! In the Baltic Sea, Sstar = sa, and note that gsw_saar returns zero
	! for saar in the Baltic.
    */
	if (saar == GSW_INVALID_VALUE)
	    return (saar);
	return (sa*(1e0 - 0.35e0*saar)/(1e0 + saar));
}

/*
!==========================================================================
function gsw_sa_from_sstar(sstar,p,lon,lat)  
!==========================================================================

! Calculates Absolute Salinity, SA, from Preformed Salinity, Sstar.
!
! Sstar  : Preformed Salinity                              [g/kg]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sa_from_sstar   : Absolute Salinity                  [g/kg]
*/
double
gsw_sa_from_sstar(double sstar, double p, double lon, double lat)
{
	double	saar;

	saar	= gsw_saar(p,lon,lat);
	if (saar == GSW_INVALID_VALUE)
	    return (saar);
    /*
    **! In the Baltic Sea, Sstar = SA, and note that gsw_saar returns zero
    **! for SAAR in the Baltic.
    */
	return (sstar*(1e0 + saar)/(1e0 - 0.35e0*saar));
}

/*
!==========================================================================
function gsw_sp_from_sstar(sstar,p,lon,lat)  
!==========================================================================

! Calculates Practical Salinity, SP, from Preformed Salinity, Sstar. 
!
! sstar  : Preformed Salinity                              [g/kg]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sp_from_Sstar : Preformed Salinity                   [g/kg]
*/
double
gsw_sp_from_sstar(double sstar, double p, double lon, double lat)
{
	double	saar, sp_baltic;

	saar	= gsw_saar(p,lon,lat);
	if (saar == GSW_INVALID_VALUE)
	    return (saar);
    /*
    **! In the Baltic Sea, SA = Sstar.
    */
	sp_baltic	= gsw_sp_from_sa_baltic(sstar,lon,lat);
	if (sp_baltic < 1810)
	    return (sp_baltic);
	return ((35.e0/35.16504e0)*sstar/(1 - 0.35e0*saar));
}

/*
!==========================================================================
function gsw_t_from_ct(sa,ct,p)  
!==========================================================================

! Calculates in-situ temperature from Conservative Temperature of seawater  
!
! sa      : Absolute Salinity                              [g/kg]
! ct      : Conservative Temperature                       [deg C]
!
! gsw_t_from_ct : in-situ temperature                      [deg C]
*/
double
gsw_t_from_ct(double sa, double ct, double p)
{
	double	pt0, p0;

	p0	= 0e0;
	pt0	= gsw_pt_from_ct(sa,ct);
	return (gsw_pt_from_t(sa,pt0,p0,p));
}

/*
!==========================================================================
function gsw_ct_from_pt(sa,pt)  
!==========================================================================

! Calculates Conservative Temperature from potential temperature of seawater  
!
! sa      : Absolute Salinity                              [g/kg]
! pt      : potential temperature with                     [deg C]
!           reference pressure of 0 dbar
!
! gsw_ct_from_pt : Conservative Temperature                [deg C]
*/
double
gsw_ct_from_pt(double sa, double pt)
{
	double	sfac, x2, x, y, cp0, pot_enthalpy;

	sfac		= 0.0248826675584615e0;
	x2		= sfac*sa;
	x		= sqrt(x2);
	y		= pt*0.025e0;	/*! normalize for F03 and F08 */
	pot_enthalpy	=  61.01362420681071e0 + y*(168776.46138048015e0 +
             y*(-2735.2785605119625e0 + y*(2574.2164453821433e0 +
             y*(-1536.6644434977543e0 + y*(545.7340497931629e0 +
             (-50.91091728474331e0 - 18.30489878927802e0*y)*y))))) +
             x2*(268.5520265845071e0 + y*(-12019.028203559312e0 +
             y*(3734.858026725145e0 + y*(-2046.7671145057618e0 +
             y*(465.28655623826234e0 + (-0.6370820302376359e0 -
             10.650848542359153e0*y)*y)))) +
             x*(937.2099110620707e0 + y*(588.1802812170108e0+
             y*(248.39476522971285e0 + (-3.871557904936333e0-
             2.6268019854268356e0*y)*y)) +
             x*(-1687.914374187449e0 + x*(246.9598888781377e0 +
             x*(123.59576582457964e0 - 48.5891069025409e0*x)) +
             y*(936.3206544460336e0 +
             y*(-942.7827304544439e0 + y*(369.4389437509002e0 +
             (-33.83664947895248e0 - 9.987880382780322e0*y)*y))))));

	cp0	= 3991.86795711963e0;
	return (pot_enthalpy/cp0);
}

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
	int	 n0, n2, no_iter;
	double	cp0, sso, s1, pt, ptm, pt_old, dentropy, dentropy_dt,
		true_entropy_part;

	n0	= 0;
	n2	= 2;
	cp0	= 3991.86795711963e0;
	sso	= 35.16504e0;
	s1	= sa*35e0/sso;
	pt	= t+(p-p_ref)*( 8.65483913395442e-6  -
		  s1 *  1.41636299744881e-6  -
		  (p+p_ref) *  7.38286467135737e-9  +
		  t  *(-8.38241357039698e-6  +
		  s1 *  2.83933368585534e-8  +
		  t  *  1.77803965218656e-8  +
		  (p+p_ref) *  1.71155619208233e-10));

	dentropy_dt	= cp0/((273.15e0 + pt)*(1e0-0.05e0*(1e0 - sa/sso)));
	true_entropy_part	= gsw_entropy_part(sa,t,p);
	for (no_iter=1; no_iter <= 2; no_iter++) {
	    pt_old	= pt;
	    dentropy	= gsw_entropy_part(sa,pt_old,p_ref) - true_entropy_part;
	    pt		= pt_old - dentropy/dentropy_dt;
	    ptm		= 0.5e0*(pt + pt_old);
	    dentropy_dt	= -gsw_gibbs(n0,n2,n0,sa,ptm,p_ref);
	    pt		= pt_old - dentropy/dentropy_dt;
	}
	return (pt);
}

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
	int	n0, n2, no_iter;
	double	pt0, pt0_old, dentropy, dentropy_dt, sso, cp0;
	double	s1, true_entropy_part, pt0m;

	n0	= 0;
	n2	= 2;

	cp0	= 3991.86795711963e0;
	sso	= 35.16504e0;

	s1	= sa*35e0/sso;

	pt0	= t+p*( 8.65483913395442e-6  -
        	  s1 *  1.41636299744881e-6  -
		   p *  7.38286467135737e-9  +
		   t *(-8.38241357039698e-6  +
		  s1 *  2.83933368585534e-8  +
		   t *  1.77803965218656e-8  +
		   p *  1.71155619208233e-10));

	dentropy_dt	= cp0/((273.15e0 + pt0)*(1e0-0.05e0*(1e0 - sa/sso)));

	true_entropy_part	= gsw_entropy_part(sa,t,p);

	for (no_iter=1; no_iter <= 2; no_iter++) {
	    pt0_old	= pt0;
	    dentropy	= gsw_entropy_part_zerop(sa,pt0_old) -
			  true_entropy_part;
	    pt0		= pt0_old - dentropy/dentropy_dt;
	    pt0m	= 0.5e0*(pt0 + pt0_old);
	    dentropy_dt	= -gsw_gibbs_pt0_pt0(sa,pt0m);
	    pt0		= pt0_old - dentropy/dentropy_dt;
	}
	return (pt0);
}

/*
!==========================================================================
function gsw_pt_from_ct(sa,ct)  
!==========================================================================

! potential temperature of seawater from conservative temperature
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_pt_from_ct : potential temperature with              [deg C]
!                  reference pressure of  0 dbar
*/
double
gsw_pt_from_ct(double sa, double ct)
{
	int	n0, n2;
	double	s1, p0, cp0 ;
	double	a0, a1, a2, a3, a4, a5, b0, b1, b2, b3;
	double	a5ct, b3ct, ct_factor, pt_num, pt_den, ct_diff;
	double	pt, pt_old, ptm, dct_dpt;

	cp0	= 3991.86795711963e0;

	n0	= 0;
	n2	= 2;

	s1	= sa*35.e0/35.16504e0;
	p0	= 0.e0;

	a0	= -1.446013646344788e-2;    
	a1	= -3.305308995852924e-3;    
	a2	=  1.062415929128982e-4;     
	a3	=  9.477566673794488e-1;     
	a4	=  2.166591947736613e-3;
	a5	=  3.828842955039902e-3;

	b0	=  1.000000000000000e0;
	b1	=  6.506097115635800e-4;
	b2	=  3.830289486850898e-3;
	b3	=  1.247811760368034e-6;

	a5ct	= a5*ct;
	b3ct	= b3*ct;
	;
	ct_factor	= (a3 + a4*s1 + a5ct);
	pt_num	= a0 + s1*(a1 + a2*s1) + ct*ct_factor;
	pt_den	= b0 + b1*s1 + ct*(b2 + b3ct);
	pt	= (pt_num)/(pt_den);

	dct_dpt	= (pt_den)/(ct_factor + a5ct - (b2 + b3ct + b3ct)*pt);

    /*
    **  Start the 1.5 iterations through the modified Newton-Raphson
    **  iterative method.
    */

	ct_diff	= gsw_ct_from_pt(sa,pt) - ct;
	pt_old	= pt;
	pt	= pt_old - (ct_diff)/dct_dpt;
	ptm	= 0.5e0*(pt + pt_old);

	dct_dpt	= -(ptm + 273.15e0)*gsw_gibbs_pt0_pt0(sa,ptm)/cp0;

	pt	= pt_old - (ct_diff)/dct_dpt;
	ct_diff	= gsw_ct_from_pt(sa,pt) - ct;
	pt_old	= pt;
	return (pt_old - (ct_diff)/dct_dpt);
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_z_from_p(p,lat)
!==========================================================================

! Calculates the height z from pressure p
!
! p      : sea pressure                                    [dbar]
! lat    : latitude                                        [deg]
!
! gsw_z_from_p : height                                    [m]
*/
double
gsw_z_from_p(double p, double lat)
{
	double	pi = 3.141592653589793e0;

	double	gamma, deg2rad, x, sin2, b, c, a;

	gamma	= 2.26e-7;
	deg2rad	= pi/180e0;
	x	= sin(lat*deg2rad);
	sin2	= x*x;
	b	= 9.780327e0*(1e0 + (5.2792e-3 + (2.32e-5*sin2))*sin2);
	a	= -0.5e0*gamma*b;
	c	= gsw_enthalpy_sso_0_p(p);

	return (-2e0*c/(b + sqrt(b*b - 4e0*a*c)));
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_entropy_from_t(sa,t,p)
!==========================================================================

! Calculates the specific entropy of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_entropy_from_t : specific entropy                    [J/(kg K)]
*/
double
gsw_entropy_from_t(double sa, double t, double p)
{
	int	n0, n1;

	n0	= 0;
	n1	= 1;

	return (-gsw_gibbs(n0,n1,n0,sa,t,p));

}

/*
!--------------------------------------------------------------------------
!==========================================================================
function gsw_adiabatic_lapse_rate_from_ct(sa,ct,p)
!==========================================================================

! Calculates the adiabatic lapse rate from Conservative Temperature
!
! sa     : Absolute Salinity                                 [g/kg]
! ct     : Conservative Temperature                          [deg C]
! p      : sea pressure                                      [dbar]
!
! gsw_adiabatic_lapse_rate_from_ct : adiabatic lapse rate    [K/Pa]
*/
double
gsw_adiabatic_lapse_rate_from_ct(double sa, double ct, double p)
{
	int	n0, n1, n2;

	double	pt0, pr0, t;

	n0	= 0;
	n1	= 1;
	n2	= 2;

	pr0	= 0e0;
	pt0	= gsw_pt_from_ct(sa,ct);
	t	= gsw_pt_from_t(sa,pt0,pr0,p);

	return (-gsw_gibbs(n0,n1,n1,sa,t,p)/(gsw_gibbs(n0,n2,n0,sa,t,p)));

}

/*
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! density and enthalpy, based on the 48-term expression for density
!--------------------------------------------------------------------------

!==========================================================================
function gsw_rho(sa,ct,p)  
!==========================================================================

!  Calculates in-situ density from Absolute Salinity and Conservative 
!  Temperature, using the computationally-efficient 48-term expression for
!  density in terms of SA, CT and p (McDougall et al., 2011).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_rho  : in-situ density (48 term equation)
*/
double
gsw_rho(double sa, double ct, double p)
{
	double	v01 =  9.998420897506056e+2, v02 =  2.839940833161907e0,
		v03 = -3.147759265588511e-2, v04 =  1.181805545074306e-3,
		v05 = -6.698001071123802e0,  v06 = -2.986498947203215e-2,
		v07 =  2.327859407479162e-4, v08 = -3.988822378968490e-2,
		v09 =  5.095422573880500e-4, v10 = -1.426984671633621e-5,
		v11 =  1.645039373682922e-7, v12 = -2.233269627352527e-2,
		v13 = -3.436090079851880e-4, v14 =  3.726050720345733e-6,
		v15 = -1.806789763745328e-4, v16 =  6.876837219536232e-7,
		v17 = -3.087032500374211e-7, v18 = -1.988366587925593e-8,
		v19 = -1.061519070296458e-11,v20 =  1.550932729220080e-10,
		v21 =  1.0e0,
		v22 =  2.775927747785646e-3, v23 = -2.349607444135925e-5,
		v24 =  1.119513357486743e-6, v25 =  6.743689325042773e-10,
		v26 = -7.521448093615448e-3, v27 = -2.764306979894411e-5,
		v28 =  1.262937315098546e-7, v29 =  9.527875081696435e-10,
		v30 = -1.811147201949891e-11, v31 = -3.303308871386421e-5,
		v32 =  3.801564588876298e-7, v33 = -7.672876869259043e-9,
		v34 = -4.634182341116144e-11, v35 =  2.681097235569143e-12,
		v36 =  5.419326551148740e-6, v37 = -2.742185394906099e-5,
		v38 = -3.212746477974189e-7, v39 =  3.191413910561627e-9,
		v40 = -1.931012931541776e-12, v41 = -1.105097577149576e-7,
		v42 =  6.211426728363857e-10, v43 = -1.119011592875110e-10,
		v44 = -1.941660213148725e-11, v45 = -1.864826425365600e-14,
		v46 =  1.119522344879478e-14, v47 = -1.200507748551599e-15,
		v48 =  6.057902487546866e-17;
	double	sqrtsa, v_hat_denominator, v_hat_numerator;

	sqrtsa			= sqrt(sa);

	v_hat_denominator	=
			v01 + ct*(v02 + ct*(v03 + v04*ct))  
			+ sa*(v05 + ct*(v06 + v07*ct) 
			+ sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct)))) 
			+ p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct) 
			+ p*(v17 + ct*(v18 + v19*ct) + v20*sa));

	v_hat_numerator		=
			v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct))) 
			+ sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct)))
			+ v36*sa 
			+ sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34+v35*ct)))))
			+ p*(v37 + ct*(v38 + ct*(v39 + v40*ct))  
			+ sa*(v41 + v42*ct) 
			+ p*(v43 + ct*(v44 + v45*ct + v46*sa) 
			+ p*(v47 + v48*ct)));

	return (v_hat_denominator/v_hat_numerator);
}

/*
!==========================================================================
function gsw_alpha(sa,ct,p)  
!==========================================================================

!  Calculates the thermal expansion coefficient of seawater with respect to 
!  Conservative Temperature using the computationally-efficient 48-term 
!  expression for density in terms of SA, CT and p (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_alpha : thermal expansion coefficient of seawater (48 term equation)
*/
double
gsw_alpha(double sa, double ct, double p)
{
	double	v01 =  9.998420897506056e+2, v02 =  2.839940833161907e0,
		v03 = -3.147759265588511e-2, v04 =  1.181805545074306e-3,
		v05 = -6.698001071123802e0,  v06 = -2.986498947203215e-2,
		v07 =  2.327859407479162e-4, v08 = -3.988822378968490e-2,
		v09 =  5.095422573880500e-4, v10 = -1.426984671633621e-5,
		v11 =  1.645039373682922e-7, v12 = -2.233269627352527e-2,
		v13 = -3.436090079851880e-4, v14 =  3.726050720345733e-6,
		v15 = -1.806789763745328e-4, v16 =  6.876837219536232e-7,
		v17 = -3.087032500374211e-7, v18 = -1.988366587925593e-8,
		v19 = -1.061519070296458e-11, v20 =  1.550932729220080e-10,
		v21 =  1.0e0,
		v22 =  2.775927747785646e-3, v23 = -2.349607444135925e-5,
		v24 =  1.119513357486743e-6, v25 =  6.743689325042773e-10,
		v26 = -7.521448093615448e-3, v27 = -2.764306979894411e-5,
		v28 =  1.262937315098546e-7, v29 =  9.527875081696435e-10,
		v30 = -1.811147201949891e-11,v31 = -3.303308871386421e-5,
		v32 =  3.801564588876298e-7, v33 = -7.672876869259043e-9,
		v34 = -4.634182341116144e-11,v35 =  2.681097235569143e-12,
		v36 =  5.419326551148740e-6, v37 = -2.742185394906099e-5,
		v38 = -3.212746477974189e-7, v39 =  3.191413910561627e-9,
		v40 = -1.931012931541776e-12,v41 = -1.105097577149576e-7,
		v42 =  6.211426728363857e-10,v43 = -1.119011592875110e-10,
		v44 = -1.941660213148725e-11,v45 = -1.864826425365600e-14,
		v46 =  1.119522344879478e-14,v47 = -1.200507748551599e-15,
		v48 =  6.057902487546866e-17,
		a01 =  2.839940833161907e0,  a02 = -6.295518531177023e-2,
		a03 =  3.545416635222918e-3, a04 = -2.986498947203215e-2,
		a05 =  4.655718814958324e-4, a06 =  5.095422573880500e-4,
		a07 = -2.853969343267241e-5, a08 =  4.935118121048767e-7,
		a09 = -3.436090079851880e-4, a10 =  7.452101440691467e-6,
		a11 =  6.876837219536232e-7, a12 = -1.988366587925593e-8,
		a13 = -2.123038140592916e-11,a14 =  2.775927747785646e-3,
		a15 = -4.699214888271850e-5, a16 =  3.358540072460230e-6,
		a17 =  2.697475730017109e-9, a18 = -2.764306979894411e-5,
		a19 =  2.525874630197091e-7, a20 =  2.858362524508931e-9,
		a21 = -7.244588807799565e-11,a22 =  3.801564588876298e-7,
		a23 = -1.534575373851809e-8, a24 = -1.390254702334843e-10,
		a25 =  1.072438894227657e-11,a26 = -3.212746477974189e-7,
		a27 =  6.382827821123254e-9, a28 = -5.793038794625329e-12,
		a29 =  6.211426728363857e-10,a30 = -1.941660213148725e-11,
		a31 = -3.729652850731201e-14,a32 =  1.119522344879478e-14,
		a33 =  6.057902487546866e-17;

	double	sqrtsa, v_hat_denominator, v_hat_numerator;
	double	dvhatden_dct, dvhatnum_dct;

	sqrtsa			= sqrt(sa);

	v_hat_denominator	=
		v01 + ct*(v02 + ct*(v03 + v04*ct))
		+ sa*(v05 + ct*(v06 + v07*ct)
		+ sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
		+ p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct)
		+ p*(v17 + ct*(v18 + v19*ct) + v20*sa));

	v_hat_numerator		=
		v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct)))
		+ sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa
		+ sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))
		+ p*(v37 + ct*(v38 + ct*(v39 + v40*ct))
		+ sa*(v41 + v42*ct)
		+ p*(v43 + ct*(v44 + v45*ct + v46*sa)
		+ p*(v47 + v48*ct)));
       
	dvhatden_dct		= a01 + ct*(a02 + a03*ct)
				+ sa*(a04 + a05*ct
				+ sqrtsa*(a06 + ct*(a07 + a08*ct)))
				+ p*(a09 + a10*ct + a11*sa
				+ p*(a12 + a13*ct));

	dvhatnum_dct		= a14 + ct*(a15 + ct*(a16 + a17*ct))
				+ sa*(a18 + ct*(a19 + ct*(a20 + a21*ct))
				+ sqrtsa*(a22 + ct*(a23 + ct*(a24 + a25*ct))))
				+ p*(a26 + ct*(a27 + a28*ct) + a29*sa
				+ p*(a30 + a31*ct + a32*sa + a33*p));
 
	return ((v_hat_denominator*dvhatnum_dct-v_hat_numerator*dvhatden_dct)/
		(v_hat_numerator*v_hat_denominator));
}

/*
!==========================================================================
function gsw_beta(sa,ct,p)  
!==========================================================================

!  Calculates the saline (i.e. haline) contraction coefficient of seawater  
!  at constant Conservative Temperature using the computationally-efficient
!  48-term expression for density in terms of SA, CT and p 
!  (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_beta : saline contraction coefficient of seawater (48 term equation)
*/
double
gsw_beta(double sa, double ct, double p)
{
	double	v01 =  9.998420897506056e+2, v02 =  2.839940833161907e0,
		v03 = -3.147759265588511e-2, v04 =  1.181805545074306e-3,
		v05 = -6.698001071123802e0,  v06 = -2.986498947203215e-2,
		v07 =  2.327859407479162e-4, v08 = -3.988822378968490e-2,
		v09 =  5.095422573880500e-4, v10 = -1.426984671633621e-5,
		v11 =  1.645039373682922e-7, v12 = -2.233269627352527e-2,
		v13 = -3.436090079851880e-4, v14 =  3.726050720345733e-6,
		v15 = -1.806789763745328e-4, v16 =  6.876837219536232e-7,
		v17 = -3.087032500374211e-7, v18 = -1.988366587925593e-8,
		v19 = -1.061519070296458e-11,v20 =  1.550932729220080e-10,
		v21 =  1.0e0,
		v22 =  2.775927747785646e-3, v23 = -2.349607444135925e-5,
		v24 =  1.119513357486743e-6, v25 =  6.743689325042773e-10,
		v26 = -7.521448093615448e-3, v27 = -2.764306979894411e-5,
		v28 =  1.262937315098546e-7, v29 =  9.527875081696435e-10,
		v30 = -1.811147201949891e-11,v31 = -3.303308871386421e-5,
		v32 =  3.801564588876298e-7, v33 = -7.672876869259043e-9,
		v34 = -4.634182341116144e-11,v35 =  2.681097235569143e-12,
		v36 =  5.419326551148740e-6, v37 = -2.742185394906099e-5,
		v38 = -3.212746477974189e-7, v39 =  3.191413910561627e-9,
		v40 = -1.931012931541776e-12,v41 = -1.105097577149576e-7,
		v42 =  6.211426728363857e-10,v43 = -1.119011592875110e-10,
		v44 = -1.941660213148725e-11,v45 = -1.864826425365600e-14,
		v46 =  1.119522344879478e-14,v47 = -1.200507748551599e-15,
		v48 =  6.057902487546866e-17,
		b01 = -6.698001071123802e0,  b02 = -2.986498947203215e-2,
		b03 =  2.327859407479162e-4, b04 = -5.983233568452735e-2,
		b05 =  7.643133860820750e-4, b06 = -2.140477007450431e-5,
		b07 =  2.467559060524383e-7, b08 = -1.806789763745328e-4,
		b09 =  6.876837219536232e-7, b10 =  1.550932729220080e-10,
		b11 = -7.521448093615448e-3, b12 = -2.764306979894411e-5,
		b13 =  1.262937315098546e-7, b14 =  9.527875081696435e-10,
		b15 = -1.811147201949891e-11,b16 = -4.954963307079632e-5,
		b17 =  5.702346883314446e-7, b18 = -1.150931530388857e-8,
		b19 = -6.951273511674217e-11,b20 =  4.021645853353715e-12,
		b21 =  1.083865310229748e-5, b22 = -1.105097577149576e-7,
		b23 =  6.211426728363857e-10,b24 =  1.119522344879478e-14;

	double	sqrtsa, v_hat_denominator, v_hat_numerator;
	double	dvhatden_dsa, dvhatnum_dsa;

	sqrtsa			= sqrt(sa);

	v_hat_denominator	=
		v01 + ct*(v02 + ct*(v03 + v04*ct))
		+ sa*(v05 + ct*(v06 + v07*ct)
		+ sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
		+ p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct)
		+ p*(v17 + ct*(v18 + v19*ct) + v20*sa));

	v_hat_numerator		=
		v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct)))
		+ sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa
		+ sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))
		+ p*(v37 + ct*(v38 + ct*(v39 + v40*ct))
		+ sa*(v41 + v42*ct)
		+ p*(v43 + ct*(v44 + v45*ct + v46*sa)
		+ p*(v47 + v48*ct)));

	dvhatden_dsa		=
		b01 + ct*(b02 + b03*ct)
		+ sqrtsa*(b04 + ct*(b05 + ct*(b06 + b07*ct)))
		+ p*(b08 + b09*ct + b10*p); 

	dvhatnum_dsa		=
		b11 + ct*(b12 + ct*(b13 + ct*(b14 + b15*ct)))
		+ sqrtsa*(b16 + ct*(b17 + ct*(b18 + ct*(b19 + b20*ct))))
		+ b21*sa
		+ p*(b22 + ct*(b23 + b24*p));

	return ((v_hat_numerator*dvhatden_dsa-v_hat_denominator*dvhatnum_dsa)/
		(v_hat_numerator*v_hat_denominator));
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_alpha_on_beta(sa,ct,p)  
!==========================================================================

!  Calculates alpha divided by beta, where alpha is the thermal expansion
!  coefficient and beta is the saline contraction coefficient of seawater 
!  from Absolute Salinity and Conservative Temperature.  This function uses
!  the computationally-efficient 48-term expression for density in terms of 
!  SA, CT and p (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_alpha_on_beta : thermal expansion coefficient of seawater (48 term equation)
*/
double
gsw_alpha_on_beta(double sa, double ct, double p)
{
	double	v01 =  9.998420897506056e2,  v02 =  2.839940833161907e0,
		v03 = -3.147759265588511e-2, v04 =  1.181805545074306e-3,
		v05 = -6.698001071123802e0,  v06 = -2.986498947203215e-2,
		v07 =  2.327859407479162e-4, v08 = -3.988822378968490e-2,
		v09 =  5.095422573880500e-4, v10 = -1.426984671633621e-5,
		v11 =  1.645039373682922e-7, v12 = -2.233269627352527e-2,
		v13 = -3.436090079851880e-4, v14 =  3.726050720345733e-6,
		v15 = -1.806789763745328e-4, v16 =  6.876837219536232e-7,
		v17 = -3.087032500374211e-7, v18 = -1.988366587925593e-8,
		v19 = -1.061519070296458e-11,v20 =  1.550932729220080e-10,
		v21 =  1.0e0,
		v22 = 2.775927747785646e-3,  v23 = -2.349607444135925e-5,
		v24 =  1.119513357486743e-6, v25 =  6.743689325042773e-10,
		v26 = -7.521448093615448e-3, v27 = -2.764306979894411e-5,
		v28 =  1.262937315098546e-7, v29 =  9.527875081696435e-10,
		v30 = -1.811147201949891e-11,v31 = -3.303308871386421e-5,
		v32 =  3.801564588876298e-7, v33 = -7.672876869259043e-9,
		v34 = -4.634182341116144e-11,v35 =  2.681097235569143e-12,
		v36 =  5.419326551148740e-6, v37 = -2.742185394906099e-5,
		v38 = -3.212746477974189e-7, v39 =  3.191413910561627e-9,
		v40 = -1.931012931541776e-12,v41 = -1.105097577149576e-7,
		v42 =  6.211426728363857e-10,v43 = -1.119011592875110e-10,
		v44 = -1.941660213148725e-11,v45 = -1.864826425365600e-14,
		v46 =  1.119522344879478e-14,v47 = -1.200507748551599e-15,
		v48 =  6.057902487546866e-17,
		a01 =  2.839940833161907e0,  a02 = -6.295518531177023e-2,
		a03 =  3.545416635222918e-3, a04 = -2.986498947203215e-2,
		a05 =  4.655718814958324e-4, a06 =  5.095422573880500e-4,
		a07 = -2.853969343267241e-5, a08 =  4.935118121048767e-7,
		a09 = -3.436090079851880e-4, a10 =  7.452101440691467e-6,
		a11 =  6.876837219536232e-7, a12 = -1.988366587925593e-8,
		a13 = -2.123038140592916e-11,a14 =  2.775927747785646e-3,
		a15 = -4.699214888271850e-5, a16 =  3.358540072460230e-6,
		a17 =  2.697475730017109e-9, a18 = -2.764306979894411e-5,
		a19 =  2.525874630197091e-7, a20 =  2.858362524508931e-9,
		a21 = -7.244588807799565e-11,a22 =  3.801564588876298e-7,
		a23 = -1.534575373851809e-8, a24 = -1.390254702334843e-10,
		a25 =  1.072438894227657e-11,a26 = -3.212746477974189e-7,
		a27 =  6.382827821123254e-9, a28 = -5.793038794625329e-12,
		a29 =  6.211426728363857e-10,a30 = -1.941660213148725e-11,
		a31 = -3.729652850731201e-14,a32 =  1.119522344879478e-14,
		a33 =  6.057902487546866e-17,
		b01 = -6.698001071123802e0,  b02 = -2.986498947203215e-2,
		b03 =  2.327859407479162e-4, b04 = -5.983233568452735e-2,
		b05 =  7.643133860820750e-4, b06 = -2.140477007450431e-5,
		b07 =  2.467559060524383e-7, b08 = -1.806789763745328e-4,
		b09 =  6.876837219536232e-7, b10 =  1.550932729220080e-10,
		b11 = -7.521448093615448e-3, b12 = -2.764306979894411e-5,
		b13 =  1.262937315098546e-7, b14 =  9.527875081696435e-10,
		b15 = -1.811147201949891e-11,b16 = -4.954963307079632e-5,
		b17 =  5.702346883314446e-7, b18 = -1.150931530388857e-8,
		b19 = -6.951273511674217e-11,b20 =  4.021645853353715e-12,
		b21 =  1.083865310229748e-5, b22 = -1.105097577149576e-7,
		b23 =  6.211426728363857e-10,b24 =  1.119522344879478e-14;

	double	sqrtsa, v_hat_denominator, v_hat_numerator;
	double	dvhatden_dct, dvhatnum_dct;
	double	dvhatden_dsa, dvhatnum_dsa;

	sqrtsa = sqrt(sa);

	v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))
			    + sa*(v05 + ct*(v06 + v07*ct)
			    + sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
			    + p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct)
			    + p*(v17 + ct*(v18 + v19*ct) + v20*sa));

	v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct)))
			    + sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct)))
			    + v36*sa + sqrtsa*
				(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))
			    + p*(v37 + ct*(v38 + ct*(v39 + v40*ct))
			    + sa*(v41 + v42*ct)
			    + p*(v43 + ct*(v44 + v45*ct + v46*sa)
			    + p*(v47 + v48*ct)));
       
	dvhatden_dct = a01 + ct*(a02 + a03*ct)
			+ sa*(a04 + a05*ct
			+ sqrtsa*(a06 + ct*(a07 + a08*ct)))
			+ p*(a09 + a10*ct + a11*sa
			+ p*(a12 + a13*ct));

	dvhatnum_dct = a14 + ct*(a15 + ct*(a16 + a17*ct))
			+ sa*(a18 + ct*(a19 + ct*(a20 + a21*ct))
			+ sqrtsa*(a22 + ct*(a23 + ct*(a24 + a25*ct))))
			+ p*(a26 + ct*(a27 + a28*ct) + a29*sa
			+ p*(a30 + a31*ct + a32*sa + a33*p));

	dvhatden_dsa = b01 + ct*(b02 + b03*ct)
			+ sqrtsa*(b04 + ct*(b05 + ct*(b06 + b07*ct)))
			+ p*(b08 + b09*ct + b10*p);

	dvhatnum_dsa = b11 + ct*(b12 + ct*(b13 + ct*(b14 + b15*ct)))
			+ sqrtsa*(b16 + ct*(b17 + ct*(b18 + ct*(b19 + b20*ct))))
			+ b21*sa + p*(b22 + ct*(b23 + b24*p));

	return ((dvhatnum_dct*v_hat_denominator - dvhatden_dct*v_hat_numerator)/
                (dvhatden_dsa*v_hat_numerator - dvhatnum_dsa*v_hat_denominator));
}

/*
!--------------------------------------------------------------------------

!==========================================================================
subroutine gsw_rho_first_derivatives(sa, ct, p, drho_dsa, drho_dct, drho_dp)
!==========================================================================

!  Calculates the three (3) partial derivatives of in situ density with 
!  respect to Absolute Salinity, Conservative Temperature and pressure.  
!  Note that the pressure derivative is done with respect to pressure in 
!  Pa, not dbar.  This function uses the computationally-efficient 48-term 
!  expression for density in terms of SA, CT and p.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! drho_dsa  : partial derivatives of density             [ kg^2/(g m^3) ]
!              with respect to Absolute Salinity
! drho_dct  : partial derivatives of density               [ kg/(K m^3) ]
!              with respect to Conservative Temperature
! drho_dp   : partial derivatives of density              [ kg/(Pa m^3) ]
!              with respect to pressure in Pa
*/
void
gsw_rho_first_derivatives(double sa, double ct, double p,
	double *drho_dsa, double *drho_dct, double *drho_dp)
{
	double	v01 =  9.998420897506056e2,  v02 =  2.839940833161907e0,
		v03 = -3.147759265588511e-2, v04 =  1.181805545074306e-3,
		v05 = -6.698001071123802e0,  v06 = -2.986498947203215e-2,
		v07 =  2.327859407479162e-4, v08 = -3.988822378968490e-2,
		v09 =  5.095422573880500e-4, v10 = -1.426984671633621e-5,
		v11 =  1.645039373682922e-7, v12 = -2.233269627352527e-2,
		v13 = -3.436090079851880e-4, v14 =  3.726050720345733e-6,
		v15 = -1.806789763745328e-4, v16 =  6.876837219536232e-7,
		v17 = -3.087032500374211e-7, v18 = -1.988366587925593e-8,
		v19 = -1.061519070296458e-11,v20 =  1.550932729220080e-10,
		v21 =  1.0e0,
		v22 = 2.775927747785646e-3,  v23 = -2.349607444135925e-5,
		v24 =  1.119513357486743e-6, v25 =  6.743689325042773e-10,
		v26 = -7.521448093615448e-3, v27 = -2.764306979894411e-5,
		v28 =  1.262937315098546e-7, v29 =  9.527875081696435e-10,
		v30 = -1.811147201949891e-11,v31 = -3.303308871386421e-5,
		v32 =  3.801564588876298e-7, v33 = -7.672876869259043e-9,
		v34 = -4.634182341116144e-11,v35 =  2.681097235569143e-12,
		v36 =  5.419326551148740e-6, v37 = -2.742185394906099e-5,
		v38 = -3.212746477974189e-7, v39 =  3.191413910561627e-9,
		v40 = -1.931012931541776e-12,v41 = -1.105097577149576e-7,
		v42 =  6.211426728363857e-10,v43 = -1.119011592875110e-10,
		v44 = -1.941660213148725e-11,v45 = -1.864826425365600e-14,
		v46 =  1.119522344879478e-14,v47 = -1.200507748551599e-15,
		v48 =  6.057902487546866e-17,
		a01 =  2.839940833161907e0,  a02 = -6.295518531177023e-2,
		a03 =  3.545416635222918e-3, a04 = -2.986498947203215e-2,
		a05 =  4.655718814958324e-4, a06 =  5.095422573880500e-4,
		a07 = -2.853969343267241e-5, a08 =  4.935118121048767e-7,
		a09 = -3.436090079851880e-4, a10 =  7.452101440691467e-6,
		a11 =  6.876837219536232e-7, a12 = -1.988366587925593e-8,
		a13 = -2.123038140592916e-11,a14 =  2.775927747785646e-3,
		a15 = -4.699214888271850e-5, a16 =  3.358540072460230e-6,
		a17 =  2.697475730017109e-9, a18 = -2.764306979894411e-5,
		a19 =  2.525874630197091e-7, a20 =  2.858362524508931e-9,
		a21 = -7.244588807799565e-11,a22 =  3.801564588876298e-7,
		a23 = -1.534575373851809e-8, a24 = -1.390254702334843e-10,
		a25 =  1.072438894227657e-11,a26 = -3.212746477974189e-7,
		a27 =  6.382827821123254e-9, a28 = -5.793038794625329e-12,
		a29 =  6.211426728363857e-10,a30 = -1.941660213148725e-11,
		a31 = -3.729652850731201e-14,a32 =  1.119522344879478e-14,
		a33 =  6.057902487546866e-17,
		b01 = -6.698001071123802e0,  b02 = -2.986498947203215e-2,
		b03 =  2.327859407479162e-4, b04 = -5.983233568452735e-2,
		b05 =  7.643133860820750e-4, b06 = -2.140477007450431e-5,
		b07 =  2.467559060524383e-7, b08 = -1.806789763745328e-4,
		b09 =  6.876837219536232e-7, b10 =  1.550932729220080e-10,
		b11 = -7.521448093615448e-3, b12 = -2.764306979894411e-5,
		b13 =  1.262937315098546e-7, b14 =  9.527875081696435e-10,
		b15 = -1.811147201949891e-11,b16 = -4.954963307079632e-5,
		b17 =  5.702346883314446e-7, b18 = -1.150931530388857e-8,
		b19 = -6.951273511674217e-11,b20 =  4.021645853353715e-12,
		b21 =  1.083865310229748e-5, b22 = -1.105097577149576e-7,
		b23 =  6.211426728363857e-10,b24 =  1.119522344879478e-14,
		c01 = -2.233269627352527e-2, c02 = -3.436090079851880e-4,
		c03 =  3.726050720345733e-6, c04 = -1.806789763745328e-4,
		c05 =  6.876837219536232e-7, c06 = -6.174065000748422e-7,
		c07 = -3.976733175851186e-8, c08 = -2.123038140592916e-11,
		c09 =  3.101865458440160e-10,c10 = -2.742185394906099e-5,
		c11 = -3.212746477974189e-7, c12 =  3.191413910561627e-9,
		c13 = -1.931012931541776e-12,c14 = -1.105097577149576e-7,
		c15 =  6.211426728363857e-10,c16 = -2.238023185750219e-10,
		c17 = -3.883320426297450e-11,c18 = -3.729652850731201e-14,
		c19 =  2.239044689758956e-14,c20 = -3.601523245654798e-15,
		c21 =  1.817370746264060e-16,pa2db = 1e-4;

	double	sqrtsa, v_hat_denominator, v_hat_numerator;
	double	dvhatden_dct, dvhatnum_dct, dvhatden_dsa, dvhatnum_dsa;
	double	dvhatden_dp, dvhatnum_dp, rho, rec_num;

	sqrtsa = sqrt(sa);

	v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))
			+ sa*(v05 + ct*(v06 + v07*ct)
			+ sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
			+ p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct)
			+ p*(v17 + ct*(v18 + v19*ct) + v20*sa));

	v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct)))
			+ sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct)))
			+ v36*sa
			+ sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))
			+ p*(v37 + ct*(v38 + ct*(v39 + v40*ct))
			+ sa*(v41 + v42*ct)
			+ p*(v43 + ct*(v44 + v45*ct + v46*sa)
			+ p*(v47 + v48*ct)));
       
	dvhatden_dct = a01 + ct*(a02 + a03*ct)
			+ sa*(a04 + a05*ct
			+ sqrtsa*(a06 + ct*(a07 + a08*ct)))
			+ p*(a09 + a10*ct + a11*sa
			+ p*(a12 + a13*ct));

	dvhatnum_dct = a14 + ct*(a15 + ct*(a16 + a17*ct))
			+ sa*(a18 + ct*(a19 + ct*(a20 + a21*ct))
			+ sqrtsa*(a22 + ct*(a23 + ct*(a24 + a25*ct))))
			+ p*(a26 + ct*(a27 + a28*ct) + a29*sa
			+ p*(a30 + a31*ct + a32*sa + a33*p));

	dvhatden_dsa = b01 + ct*(b02 + b03*ct)
			+ sqrtsa*(b04 + ct*(b05 + ct*(b06 + b07*ct)))
			+ p*(b08 + b09*ct + b10*p);

	dvhatnum_dsa = b11 + ct*(b12 + ct*(b13 + ct*(b14 + b15*ct)))
			+ sqrtsa*(b16 + ct*(b17 + ct*(b18 + ct*(b19 + b20*ct))))
			+ b21*sa
			+ p*(b22 + ct*(b23 + b24*p));
		
	dvhatden_dp = c01 + ct*(c02 + c03*ct)
			+ sa*(c04 + c05*ct)
			+ p*(c06 + ct*(c07 + c08*ct) + c09*sa);

	dvhatnum_dp = c10 + ct*(c11 + ct*(c12 + c13*ct))
			+ sa*(c14 + c15*ct)
			+ p*(c16 + ct*(c17 + c18*ct + c19*sa)
			+ p*(c20 + c21*ct));

	rec_num = 1e0/v_hat_numerator;
       
	rho = rec_num*v_hat_denominator;

	*drho_dsa = (dvhatden_dsa - dvhatnum_dsa*rho)*rec_num;

	*drho_dct = (dvhatden_dct - dvhatnum_dct*rho)*rec_num;

	*drho_dp = pa2db*(dvhatden_dp - dvhatnum_dp*rho)*rec_num;

	return;
}

/*
!==========================================================================
function gsw_specvol(sa,ct,p)
!==========================================================================

!  Calculates specific volume of seawater using the computationally-
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol  :  specific volume of seawater (48 term equation)
*/
double
gsw_specvol(double sa, double ct, double p)
{
	return (1e0/gsw_rho(sa,ct,p));
}

/*
!==========================================================================
function gsw_specvol_anom(sa,ct,p)
!==========================================================================

!  Calculates specific volume anomaly of seawater using the computationally-
!  efficient 48-term expression for density in terms of SA, CT and p
!  (McDougall et al., 2011)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol_anom  :  specific volume anomaly of seawater (48 term equation)
*/
double
gsw_specvol_anom(double sa, double ct, double p)
{
	return (gsw_specvol(sa,ct,p) - gsw_specvol_sso_0_p(p));
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_sigma0(sa,ct) 
!==========================================================================

!  Calculates potential density anomaly with reference pressure of 0 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  function has inputs of Absolute Salinity and Conservative Temperature.
!  This function uses the computationally-efficient 48-term expression for 
!  density in terms of SA, CT and p (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! gsw_sigma0  : potential density anomaly with reference pressure of 0
!                                                      (48 term equation)
*/
double
gsw_sigma0(double sa, double ct)
{
	double	v01 =  9.998420897506056e2,  v02 =  2.839940833161907e0,
		v03 = -3.147759265588511e-2, v04 =  1.181805545074306e-3,
		v05 = -6.698001071123802e0,  v06 = -2.986498947203215e-2,
		v07 =  2.327859407479162e-4, v08 = -3.988822378968490e-2,
		v09 =  5.095422573880500e-4, v10 = -1.426984671633621e-5,
		v11 =  1.645039373682922e-7,
		v21 =  1.0e0,
		v22 =  2.775927747785646e-3, v23 = -2.349607444135925e-5,
		v24 =  1.119513357486743e-6, v25 =  6.743689325042773e-10,
		v26 = -7.521448093615448e-3, v27 = -2.764306979894411e-5,
		v28 =  1.262937315098546e-7, v29 =  9.527875081696435e-10,
		v30 = -1.811147201949891e-11,v31 = -3.303308871386421e-5,
		v32 =  3.801564588876298e-7, v33 = -7.672876869259043e-9,
		v34 = -4.634182341116144e-11,v35 =  2.681097235569143e-12,
		v36 =  5.419326551148740e-6;

	double	v_hat_denominator, v_hat_numerator, sqrtsa;

	sqrtsa = sqrt(sa);

	v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))
		+ sa*(v05 + ct*(v06 + v07*ct)
		+ sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))));

	v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct)))
		+ sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa
		+ sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))));

	return (v_hat_denominator/v_hat_numerator  - 1000e0);
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_sigma1(sa,ct) 
!==========================================================================

!  Calculates potential density anomaly with reference pressure of 1000 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  function has inputs of Absolute Salinity and Conservative Temperature.
!  This function uses the computationally-efficient 48-term expression for 
!  density in terms of SA, CT and p (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! gsw_sigma1  : potential density anomaly with reference pressure of 1000
!                                                      (48 term equation)
*/
double
gsw_sigma1(double sa, double ct)
{
	return (gsw_rho(sa,ct,1000e0) - 1000e0);
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_sigma2(sa,ct) 
!==========================================================================

!  Calculates potential density anomaly with reference pressure of 2000 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  function has inputs of Absolute Salinity and Conservative Temperature.
!  This function uses the computationally-efficient 48-term expression for 
!  density in terms of SA, CT and p (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! gsw_sigma2  : potential density anomaly with reference pressure of 2000
!                                                      (48 term equation)
*/
double
gsw_sigma2(double sa, double ct)
{
	return (gsw_rho(sa,ct,2000e0) - 1000e0);
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_sigma3(sa,ct) 
!==========================================================================

!  Calculates potential density anomaly with reference pressure of 3000 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  function has inputs of Absolute Salinity and Conservative Temperature.
!  This function uses the computationally-efficient 48-term expression for 
!  density in terms of SA, CT and p (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! gsw_sigma3  : potential density anomaly with reference pressure of 3000
!                                                      (48 term equation)
*/
double
gsw_sigma3(double sa, double ct)
{
	return (gsw_rho(sa,ct,3000e0) - 1000e0);
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_sigma4(sa,ct) 
!==========================================================================

!  Calculates potential density anomaly with reference pressure of 4000 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  function has inputs of Absolute Salinity and Conservative Temperature.
!  This function uses the computationally-efficient 48-term expression for 
!  density in terms of SA, CT and p (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! gsw_sigma4  : potential density anomaly with reference pressure of 4000
!                                                      (48 term equation)
*/
double
gsw_sigma4(double sa, double ct)
{
	return (gsw_rho(sa,ct,4000e0) - 1000e0);
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_sound_speed(sa,ct,p)  
!==========================================================================

!  Calculates sound speed of seawater using the computationally-
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_sound_speed  :  sound speed of seawater (48 term equation)
*/
double
gsw_sound_speed(double sa, double ct, double p)
{
	double	v01 =  9.998420897506056e+2, v02 =  2.839940833161907e0,
		v03 = -3.147759265588511e-2, v04 =  1.181805545074306e-3,
		v05 = -6.698001071123802e0,  v06 = -2.986498947203215e-2,
		v07 =  2.327859407479162e-4, v08 = -3.988822378968490e-2,
		v09 =  5.095422573880500e-4, v10 = -1.426984671633621e-5,
		v11 =  1.645039373682922e-7, v12 = -2.233269627352527e-2,
		v13 = -3.436090079851880e-4, v14 =  3.726050720345733e-6,
		v15 = -1.806789763745328e-4, v16 =  6.876837219536232e-7,
		v17 = -3.087032500374211e-7, v18 = -1.988366587925593e-8,
		v19 = -1.061519070296458e-11,v20 =  1.550932729220080e-10,
		v21 =  1.0e0,
		v22 =  2.775927747785646e-3, v23 = -2.349607444135925e-5,
		v24 =  1.119513357486743e-6, v25 =  6.743689325042773e-10,
		v26 = -7.521448093615448e-3, v27 = -2.764306979894411e-5,
		v28 =  1.262937315098546e-7, v29 =  9.527875081696435e-10,
		v30 = -1.811147201949891e-11,v31 = -3.303308871386421e-5,
		v32 =  3.801564588876298e-7, v33 = -7.672876869259043e-9,
		v34 = -4.634182341116144e-11,v35 =  2.681097235569143e-12,
		v36 =  5.419326551148740e-6, v37 = -2.742185394906099e-5,
		v38 = -3.212746477974189e-7, v39 =  3.191413910561627e-9,
		v40 = -1.931012931541776e-12,v41 = -1.105097577149576e-7,
		v42 =  6.211426728363857e-10,v43 = -1.119011592875110e-10,
		v44 = -1.941660213148725e-11,v45 = -1.864826425365600e-14,
		v46 =  1.119522344879478e-14,v47 = -1.200507748551599e-15,
		v48 =  6.057902487546866e-17 ,
		c01 = -2.233269627352527e-2, c02 = -3.436090079851880e-4,
		c03 =  3.726050720345733e-6, c04 = -1.806789763745328e-4,
		c05 =  6.876837219536232e-7, c06 = -6.174065000748422e-7,
		c07 = -3.976733175851186e-8, c08 = -2.123038140592916e-11,
		c09 =  3.101865458440160e-10,c10 = -2.742185394906099e-5,
		c11 = -3.212746477974189e-7, c12 =  3.191413910561627e-9,
		c13 = -1.931012931541776e-12,c14 = -1.105097577149576e-7,
		c15 =  6.211426728363857e-10,c16 = -2.238023185750219e-10,
		c17 = -3.883320426297450e-11,c18 = -3.729652850731201e-14,
		c19 =  2.239044689758956e-14,c20 = -3.601523245654798e-15,
		c21 =  1.817370746264060e-16;

	double	v_hat_denominator, v_hat_numerator;
	double	sqrtsa, dvden_dp, dvnum_dp, dp_drho;

	sqrtsa			= sqrt(sa);

	v_hat_denominator	=
		v01 + ct*(v02 + ct*(v03 + v04*ct))
		+ sa*(v05 + ct*(v06 + v07*ct)
		+ sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
		+ p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct)
		+ p*(v17 + ct*(v18 + v19*ct) + v20*sa));

	v_hat_numerator		=
		v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct)))
		+ sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa
		+ sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))
		+ p*(v37 + ct*(v38 + ct*(v39 + v40*ct))
		+ sa*(v41 + v42*ct)
		+ p*(v43 + ct*(v44 + v45*ct + v46*sa)
		+ p*(v47 + v48*ct)));

	dvden_dp		=
		c01 + ct*(c02 + c03*ct)
		+ sa*(c04 + c05*ct)
		+ p*(c06 + ct*(c07 + c08*ct) + c09*sa);

	dvnum_dp		=
		c10 + ct*(c11 + ct*(c12 + c13*ct))
		+ sa*(c14 + c15*ct)
		+ p*(c16 + ct*(c17 + c18*ct + c19*sa)
		+ p*(c20 + c21*ct));

	dp_drho			=
		(v_hat_numerator*v_hat_numerator)/
		(dvden_dp*v_hat_numerator - dvnum_dp*v_hat_denominator);
    
	return (100.0*sqrt(dp_drho));
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_kappa(sa,ct,p)
!==========================================================================

!  Calculates isentropic compressibility of seawater using the computationally
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_kappa  :  isentropic compressibility (48 term equation)
*/
double
gsw_kappa(double sa, double ct, double p)
{
	double	v01 =  9.998420897506056e2,   v02 =  2.839940833161907e0,
		v03 = -3.147759265588511e-2,  v04 =  1.181805545074306e-3,
		v05 = -6.698001071123802e0,   v06 = -2.986498947203215e-2,
		v07 =  2.327859407479162e-4,  v08 = -3.988822378968490e-2,
		v09 =  5.095422573880500e-4,  v10 = -1.426984671633621e-5,
		v11 =  1.645039373682922e-7,  v12 = -2.233269627352527e-2,
		v13 = -3.436090079851880e-4,  v14 =  3.726050720345733e-6,
		v15 = -1.806789763745328e-4,  v16 =  6.876837219536232e-7,
		v17 = -3.087032500374211e-7,  v18 = -1.988366587925593e-8,
		v19 = -1.061519070296458e-11, v20 =  1.550932729220080e-10,
		v21 =  1.0e0,
		v22 = 2.775927747785646e-3,   v23 = -2.349607444135925e-5,
		v24 =  1.119513357486743e-6,  v25 =  6.743689325042773e-10,
		v26 = -7.521448093615448e-3,  v27 = -2.764306979894411e-5,
		v28 =  1.262937315098546e-7,  v29 =  9.527875081696435e-10,
		v30 = -1.811147201949891e-11, v31 = -3.303308871386421e-5,
		v32 =  3.801564588876298e-7,  v33 = -7.672876869259043e-9,
		v34 = -4.634182341116144e-11, v35 =  2.681097235569143e-12,
		v36 =  5.419326551148740e-6,  v37 = -2.742185394906099e-5,
		v38 = -3.212746477974189e-7,  v39 =  3.191413910561627e-9,
		v40 = -1.931012931541776e-12, v41 = -1.105097577149576e-7,
		v42 =  6.211426728363857e-10, v43 = -1.119011592875110e-10,
		v44 = -1.941660213148725e-11, v45 = -1.864826425365600e-14,
		v46 =  1.119522344879478e-14, v47 = -1.200507748551599e-15,
		v48 =  6.057902487546866e-17,
		c01 = -2.233269627352527e-2,  c02 = -3.436090079851880e-4,
		c03 =  3.726050720345733e-6,  c04 = -1.806789763745328e-4,
		c05 =  6.876837219536232e-7,  c06 = -6.174065000748422e-7,
		c07 = -3.976733175851186e-8,  c08 = -2.123038140592916e-11,
		c09 =  3.101865458440160e-10, c10 = -2.742185394906099e-5,
		c11 = -3.212746477974189e-7,  c12 =  3.191413910561627e-9,
		c13 = -1.931012931541776e-12, c14 = -1.105097577149576e-7,
		c15 =  6.211426728363857e-10, c16 = -2.238023185750219e-10,
		c17 = -3.883320426297450e-11, c18 = -3.729652850731201e-14,
		c19 =  2.239044689758956e-14, c20 = -3.601523245654798e-15,
		c21 =  1.817370746264060e-16;

	double	v_hat_denominator, v_hat_numerator, sqrtsa, dvden_dp, dvnum_dp;

	sqrtsa	= sqrt(sa);

	v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))
				+ sa*(v05 + ct*(v06 + v07*ct)
				+ sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
				+ p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct)
				+ p*(v17 + ct*(v18 + v19*ct) + v20*sa));

	v_hat_numerator	= v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct)))
		+ sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa
		+ sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct))))) 
		+ p*(v37 + ct*(v38 + ct*(v39 + v40*ct)) 
		+ sa*(v41 + v42*ct)
		+ p*(v43 + ct*(v44 + v45*ct + v46*sa)
		+ p*(v47 + v48*ct)));

	dvden_dp	= c01 + ct*(c02 + c03*ct)
				+ sa*(c04 + c05*ct)
				+ p*(c06 + ct*(c07 + c08*ct) + c09*sa);

	dvnum_dp	= c10 + ct*(c11 + ct*(c12 + c13*ct))
				+ sa*(c14 + c15*ct)
				+ p*(c16 + ct*(c17 + c18*ct + c19*sa)
				+ p*(c20 + c21*ct));


	return (1e-4*(dvden_dp/v_hat_denominator - dvnum_dp/v_hat_numerator));

}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_cabbeling(sa,ct,p)  
!==========================================================================

!  Calculates the cabbeling coefficient of seawater with respect to  
!  Conservative Temperature.  This function uses the computationally-
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_cabbeling  : cabbeling coefficient with respect to            [ 1/K^2 ]
!                  Conservative Temperature. (48 term equation)
*/
double
gsw_cabbeling(double sa, double ct, double p)
{
	double	v01 =  9.998420897506056e2,   v02 =  2.839940833161907e0,
		v03 = -3.147759265588511e-2,  v04 =  1.181805545074306e-3,
		v05 = -6.698001071123802e0,   v06 = -2.986498947203215e-2,
		v07 =  2.327859407479162e-4,  v08 = -3.988822378968490e-2,
		v09 =  5.095422573880500e-4,  v10 = -1.426984671633621e-5,
		v11 =  1.645039373682922e-7,  v12 = -2.233269627352527e-2,
		v13 = -3.436090079851880e-4,  v14 =  3.726050720345733e-6,
		v15 = -1.806789763745328e-4,  v16 =  6.876837219536232e-7,
		v17 = -3.087032500374211e-7,  v18 = -1.988366587925593e-8,
		v19 = -1.061519070296458e-11, v20 =  1.550932729220080e-10,
		v21 =  1.0e0,
		v22 =  2.775927747785646e-3,  v23 = -2.349607444135925e-5,
		v24 =  1.119513357486743e-6,  v25 =  6.743689325042773e-10,
		v26 = -7.521448093615448e-3,  v27 = -2.764306979894411e-5,
		v28 =  1.262937315098546e-7,  v29 =  9.527875081696435e-10,
		v30 = -1.811147201949891e-11, v31 = -3.303308871386421e-5,
		v32 =  3.801564588876298e-7,  v33 = -7.672876869259043e-9,
		v34 = -4.634182341116144e-11, v35 =  2.681097235569143e-12,
		v36 =  5.419326551148740e-6,  v37 = -2.742185394906099e-5,
		v38 = -3.212746477974189e-7,  v39 =  3.191413910561627e-9,
		v40 = -1.931012931541776e-12, v41 = -1.105097577149576e-7,
		v42 =  6.211426728363857e-10, v43 = -1.119011592875110e-10,
		v44 = -1.941660213148725e-11, v45 = -1.864826425365600e-14,
		v46 =  1.119522344879478e-14, v47 = -1.200507748551599e-15,
		v48 =  6.057902487546866e-17 ,
		a01 =  2.839940833161907e0,   a02 = -6.295518531177023e-2,
		a03 =  3.545416635222918e-3,  a04 = -2.986498947203215e-2,
		a05 =  4.655718814958324e-4,  a06 =  5.095422573880500e-4,
		a07 = -2.853969343267241e-5,  a08 =  4.935118121048767e-7,
		a09 = -3.436090079851880e-4,  a10 =  7.452101440691467e-6,
		a11 =  6.876837219536232e-7,  a12 = -1.988366587925593e-8,
		a13 = -2.123038140592916e-11, a14 =  2.775927747785646e-3,
		a15 = -4.699214888271850e-5,  a16 =  3.358540072460230e-6,
		a17 =  2.697475730017109e-9,  a18 = -2.764306979894411e-5,
		a19 =  2.525874630197091e-7,  a20 =  2.858362524508931e-9,
		a21 = -7.244588807799565e-11, a22 =  3.801564588876298e-7,
		a23 = -1.534575373851809e-8,  a24 = -1.390254702334843e-10,
		a25 =  1.072438894227657e-11, a26 = -3.212746477974189e-7,
		a27 =  6.382827821123254e-9,  a28 = -5.793038794625329e-12,
		a29 =  6.211426728363857e-10, a30 = -1.941660213148725e-11,
		a31 = -3.729652850731201e-14, a32 =  1.119522344879478e-14,
		a33 =  6.057902487546866e-17, 
		b01 = -6.698001071123802e0,   b02 = -2.986498947203215e-2,
		b03 =  2.327859407479162e-4,  b04 = -5.983233568452735e-2,
		b05 =  7.643133860820750e-4,  b06 = -2.140477007450431e-5,
		b07 =  2.467559060524383e-7,  b08 = -1.806789763745328e-4,
		b09 =  6.876837219536232e-7,  b10 =  1.550932729220080e-10,
		b11 = -7.521448093615448e-3,  b12 = -2.764306979894411e-5,
		b13 =  1.262937315098546e-7,  b14 =  9.527875081696435e-10,
		b15 = -1.811147201949891e-11, b16 = -4.954963307079632e-5,
		b17 =  5.702346883314446e-7,  b18 = -1.150931530388857e-8,
		b19 = -6.951273511674217e-11, b20 =  4.021645853353715e-12,
		b21 =  1.083865310229748e-5,  b22 = -1.105097577149576e-7,
		b23 =  6.211426728363857e-10, b24 =  1.119522344879478e-14;

	double	sqrtsa, v_hat_denominator, v_hat_numerator;
	double	dvhatden_dct, dvhatnum_dct, dvhatden_dctdct, dvhatden_dctdsa;
	double	dvhatden_dsa, dvhatnum_dsa, dvhatnum_dsadsa, dvhatden_dsadsa;
	double	dvhatnum_dctdct, dvhatnum_dsadct, dvhatnum_dctdsa;
	double	p1a, p1b, p1c, p1d, part1, factor2a, factor2b;
	double	p2a, p2b, p2c, p2d, p2e, part2;
	double	factor3a, factor3b, p3a, p3b, p3c, p3d, part3;

	sqrtsa = sqrt(sa);

	v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))
			+ sa*(v05 + ct*(v06 + v07*ct)
			+ sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
			+ p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct)
			+ p*(v17 + ct*(v18 + v19*ct) + v20*sa));

	v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct)))
			+ sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct)))
			+ v36*sa
			+ sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))
			+ p*(v37 + ct*(v38 + ct*(v39 + v40*ct))
			+ sa*(v41 + v42*ct)
			+ p*(v43 + ct*(v44 + v45*ct + v46*sa)
			+ p*(v47 + v48*ct)));

	dvhatden_dct = a01 + ct*(a02 + a03*ct)
			+ sa*(a04 + a05*ct
			+ sqrtsa*(a06 + ct*(a07 + a08*ct)))
			+ p*(a09 + a10*ct + a11*sa
			+ p*(a12 + a13*ct));

	dvhatnum_dct = a14 + ct*(a15 + ct*(a16 + a17*ct))
			+ sa*(a18 + ct*(a19 + ct*(a20 + a21*ct))
			+ sqrtsa*(a22 + ct*(a23 + ct*(a24 + a25*ct))))
			+ p*(a26 + ct*(a27 + a28*ct) + a29*sa
			+ p*(a30 + a31*ct + a32*sa + a33*p));

	dvhatden_dctdct = a02 + 2e0*a03*ct
			+ sa*(a05 + sqrtsa*(a07 + 2e0*a08*ct))
			+ p*(a10 + a13*p);
     
	dvhatden_dctdsa = a04 + a05*ct
			+ (3e0/2e0)*sqrtsa*(a06 + ct*(a07 + a08*ct))
			+ a11*p;

	dvhatden_dsa = b01 + ct*(b02 + b03*ct)
			+ sqrtsa*(b04 + ct*(b05 + ct*(b06 + b07*ct)))
			+ p*(b08 + b09*ct + b10*p);

	dvhatnum_dsa = b11 + ct*(b12 + ct*(b13 + ct*(b14 + b15*ct)))
			+ sqrtsa*(b16 + ct*(b17 + ct*(b18 + ct*(b19 + b20*ct))))
			+ b21*sa
			+ p*(b22 + ct*(b23 + b24*p));

	dvhatnum_dctdct = a15 + ct*(2*a16 + 3*a17*ct)
			+ sa*(a19 + ct*(2*a20 + 3*a21*ct)
			+sqrtsa*(a23 + ct*(2*a24 + 3*a25*ct)))
			+ p*(a27 + 2*a28*ct + a31*p);

	dvhatnum_dctdsa = a18 + ct*(a19 + ct*(a20 + a21*ct))
			+ (3e0/2e0)*sqrtsa*(a22 + ct*(a23 + ct*(a24 + a25*ct)))
			+ p*(a29 + p*a32);

	dvhatnum_dsadsa = (1e0/(2e0*sqrtsa))*(b16 +
			   ct*(b17 + ct*(b18 + ct*(b19 + b20*ct)))) + b21;

	dvhatden_dsadsa = (1e0/(2e0*sqrtsa))*(b04 +
			   ct*(b05 + ct*(b06 + b07*ct)));

	p1a = dvhatnum_dctdct/v_hat_numerator;
	p1b = (dvhatnum_dct*dvhatden_dct)/(v_hat_numerator*v_hat_denominator);
	p1c = dvhatden_dctdct/v_hat_denominator;
	p1d = dvhatden_dct/v_hat_denominator;

	part1 = p1a - 2e0*p1b - p1c + 2e0*p1d*p1d;

	factor2a = (v_hat_denominator*dvhatnum_dct
		   - v_hat_numerator*dvhatden_dct);
	factor2b = (v_hat_denominator*dvhatnum_dsa
		   - v_hat_numerator*dvhatden_dsa);

	p2a = dvhatnum_dctdsa/v_hat_numerator;
	p2b = (dvhatnum_dct*dvhatden_dsa)/(v_hat_numerator*v_hat_denominator);
	p2c = (dvhatnum_dsa*dvhatden_dct)/(v_hat_numerator*v_hat_denominator);
	p2d = dvhatden_dctdsa/v_hat_denominator;
	p2e = (dvhatden_dct*dvhatden_dsa)/(v_hat_denominator*v_hat_denominator);

	part2 = p2a - p2b - p2c - p2d + 2e0*p2e;

	factor3a = (v_hat_denominator*dvhatnum_dct
		   - v_hat_numerator*dvhatden_dct);
	factor3b = (v_hat_denominator*dvhatnum_dsa
		   - v_hat_numerator*dvhatden_dsa);

	p3a = dvhatnum_dsadsa/v_hat_numerator;
	p3b = (dvhatnum_dsa*dvhatden_dsa)/(v_hat_numerator*v_hat_denominator);
	p3c = (v_hat_numerator*dvhatden_dsadsa)/
		(v_hat_numerator*v_hat_denominator);
	p3d = dvhatden_dsa/v_hat_denominator;

	part3 = p3a - 2e0*p3b - p3c + 2e0*p3d*p3d;

	return (part1 - 2e0*(factor2a/factor2b)*part2 +
		(factor3a/factor3b)*(factor3a/factor3b)*part3);
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_thermobaric(sa,ct,p)  
!==========================================================================

!  Calculates the thermobaric coefficient of seawater with respect to
!  Conservative Temperature.  This routine calculates rho from the 
!  computationally-efficient 48-term expression for density in terms of
!  SA, CT and p (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_thermobaric  : thermobaric coefficient with          [ 1/(K Pa) ] 
!                    respect to Conservative Temperature (48 term equation)
*/
double
gsw_thermobaric(double sa, double ct, double p)
{
	double	v01 =  9.998420897506056e2,   v02 =  2.839940833161907e0,
		v03 = -3.147759265588511e-2,  v04 =  1.181805545074306e-3,
		v05 = -6.698001071123802e0,   v06 = -2.986498947203215e-2,
		v07 =  2.327859407479162e-4,  v08 = -3.988822378968490e-2,
		v09 =  5.095422573880500e-4,  v10 = -1.426984671633621e-5,
		v11 =  1.645039373682922e-7,  v12 = -2.233269627352527e-2,
		v13 = -3.436090079851880e-4,  v14 =  3.726050720345733e-6,
		v15 = -1.806789763745328e-4,  v16 =  6.876837219536232e-7,
		v17 = -3.087032500374211e-7,  v18 = -1.988366587925593e-8,
		v19 = -1.061519070296458e-11, v20 =  1.550932729220080e-10,
		v21 =  1.0e0,
		v22 =  2.775927747785646e-3,  v23 = -2.349607444135925e-5,
		v24 =  1.119513357486743e-6,  v25 =  6.743689325042773e-10,
		v26 = -7.521448093615448e-3,  v27 = -2.764306979894411e-5,
		v28 =  1.262937315098546e-7,  v29 =  9.527875081696435e-10,
		v30 = -1.811147201949891e-11, v31 = -3.303308871386421e-5,
		v32 =  3.801564588876298e-7,  v33 = -7.672876869259043e-9,
		v34 = -4.634182341116144e-11, v35 =  2.681097235569143e-12,
		v36 =  5.419326551148740e-6,  v37 = -2.742185394906099e-5,
		v38 = -3.212746477974189e-7,  v39 =  3.191413910561627e-9,
		v40 = -1.931012931541776e-12, v41 = -1.105097577149576e-7,
		v42 =  6.211426728363857e-10, v43 = -1.119011592875110e-10,
		v44 = -1.941660213148725e-11, v45 = -1.864826425365600e-14,
		v46 =  1.119522344879478e-14, v47 = -1.200507748551599e-15,
		v48 =  6.057902487546866e-17 ,
		a01 =  2.839940833161907e0,   a02 = -6.295518531177023e-2,
		a03 =  3.545416635222918e-3,  a04 = -2.986498947203215e-2,
		a05 =  4.655718814958324e-4,  a06 =  5.095422573880500e-4,
		a07 = -2.853969343267241e-5,  a08 =  4.935118121048767e-7,
		a09 = -3.436090079851880e-4,  a10 =  7.452101440691467e-6,
		a11 =  6.876837219536232e-7,  a12 = -1.988366587925593e-8,
		a13 = -2.123038140592916e-11, a14 =  2.775927747785646e-3,
		a15 = -4.699214888271850e-5,  a16 =  3.358540072460230e-6,
		a17 =  2.697475730017109e-9,  a18 = -2.764306979894411e-5,
		a19 =  2.525874630197091e-7,  a20 =  2.858362524508931e-9,
		a21 = -7.244588807799565e-11, a22 =  3.801564588876298e-7,
		a23 = -1.534575373851809e-8,  a24 = -1.390254702334843e-10,
		a25 =  1.072438894227657e-11, a26 = -3.212746477974189e-7,
		a27 =  6.382827821123254e-9,  a28 = -5.793038794625329e-12,
		a29 =  6.211426728363857e-10, a30 = -1.941660213148725e-11,
		a31 = -3.729652850731201e-14, a32 =  1.119522344879478e-14,
		a33 =  6.057902487546866e-17, 
		b01 = -6.698001071123802e0,   b02 = -2.986498947203215e-2,
		b03 =  2.327859407479162e-4,  b04 = -5.983233568452735e-2,
		b05 =  7.643133860820750e-4,  b06 = -2.140477007450431e-5,
		b07 =  2.467559060524383e-7,  b08 = -1.806789763745328e-4,
		b09 =  6.876837219536232e-7,  b10 =  1.550932729220080e-10,
		b11 = -7.521448093615448e-3,  b12 = -2.764306979894411e-5,
		b13 =  1.262937315098546e-7,  b14 =  9.527875081696435e-10,
		b15 = -1.811147201949891e-11, b16 = -4.954963307079632e-5,
		b17 =  5.702346883314446e-7,  b18 = -1.150931530388857e-8,
		b19 = -6.951273511674217e-11, b20 =  4.021645853353715e-12,
		b21 =  1.083865310229748e-5,  b22 = -1.105097577149576e-7,
		b23 =  6.211426728363857e-10, b24 =  1.119522344879478e-14,
		c01 = -2.233269627352527e-2,  c02 = -3.436090079851880e-4,
		c03 =  3.726050720345733e-6,  c04 = -1.806789763745328e-4,
		c05 =  6.876837219536232e-7,  c06 = -6.174065000748422e-7,
		c07 = -3.976733175851186e-8,  c08 = -2.123038140592916e-11,
		c09 =  3.101865458440160e-10, c10 = -2.742185394906099e-5,
		c11 = -3.212746477974189e-7,  c12 =  3.191413910561627e-9,
		c13 = -1.931012931541776e-12, c14 = -1.105097577149576e-7,
		c15 =  6.211426728363857e-10, c16 = -2.238023185750219e-10,
		c17 = -3.883320426297450e-11, c18 = -3.729652850731201e-14,
		c19 =  2.239044689758956e-14, c20 = -3.601523245654798e-15,
		c21 =  1.817370746264060e-16,rec_db2pa = 1e-4;

	double	sqrtsa, v_hat_denominator, v_hat_numerator;
	double	dvhatden_dct, dvhatnum_dct, dvhatden_dp, dvhatnum_dp;
	double	dvhatden_dsa, dvhatnum_dsa, dvhatden_dpdct, dvhatnum_dpdct;
	double	dvhatden_dpdsa, dvhatnum_dpdsa;
	double	p1a, p1b, p1c, p1d, p1e, part1, factor2;
	double	p2a, p2b, p2c, p2d, p2e, part2;

	sqrtsa = sqrt(sa);

	v_hat_denominator = v01 + ct*(v02 + ct*(v03 + v04*ct))
		+ sa*(v05 + ct*(v06 + v07*ct)
		+ sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))))
		+ p*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct)
		+ p*(v17 + ct*(v18 + v19*ct) + v20*sa));

	v_hat_numerator = v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct)))
		+ sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa
		+ sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))))
		+ p*(v37 + ct*(v38 + ct*(v39 + v40*ct))
		+ sa*(v41 + v42*ct)
		+ p*(v43 + ct*(v44 + v45*ct + v46*sa)
		+ p*(v47 + v48*ct)));

	dvhatden_dct = a01 + ct*(a02 + a03*ct)
		+ sa*(a04 + a05*ct
		+ sqrtsa*(a06 + ct*(a07 + a08*ct)))
		+ p*(a09 + a10*ct + a11*sa
		+ p*(a12 + a13*ct));

	dvhatnum_dct = a14 + ct*(a15 + ct*(a16 + a17*ct))
		+ sa*(a18 + ct*(a19 + ct*(a20 + a21*ct))
		+ sqrtsa*(a22 + ct*(a23 + ct*(a24 + a25*ct))))
		+ p*(a26 + ct*(a27 + a28*ct) + a29*sa
		+ p*(a30 + a31*ct + a32*sa + a33*p));

	dvhatden_dsa = b01 + ct*(b02 + b03*ct)
		+ sqrtsa*(b04 + ct*(b05 + ct*(b06 + b07*ct)))
		+ p*(b08 + b09*ct + b10*p);

	dvhatnum_dsa = b11 + ct*(b12 + ct*(b13 + ct*(b14 + b15*ct)))
		+ sqrtsa*(b16 + ct*(b17 + ct*(b18 + ct*(b19 + b20*ct))))
		+ b21*sa + p*(b22 + ct*(b23 + b24*p));

	dvhatden_dp = c01 + ct*(c02 + c03*ct)
		+ sa*(c04 + c05*ct)
		+ p*(c06 + ct*(c07 + c08*ct) + c09*sa);

	dvhatnum_dp = c10 + ct*(c11 + ct*(c12 + c13*ct))
		+ sa*(c14 + c15*ct)
		+ p*(c16 + ct*(c17 + c18*ct + c19*sa)
		+ p*(c20 + c21*ct));

	dvhatden_dpdct = c02 + 2e0*c03*ct + c05*sa
		+ p*(c07 + 2e0*c08*ct);
       
	dvhatnum_dpdct = c11 + ct*(2e0*c12 + 3e0*c13*ct) + c15*sa
		+ p*(c17 + ct*2e0*c18 + c19*sa + c21*p);

	dvhatden_dpdsa = c04 + c05*ct + c09*p;

	dvhatnum_dpdsa = c14 + c15*ct + c19*ct*p;

	p1a = dvhatnum_dpdct/v_hat_numerator;
	p1b = (dvhatnum_dct*dvhatden_dp)/(v_hat_numerator*v_hat_denominator);
	p1c = (dvhatnum_dp*dvhatden_dct)/(v_hat_numerator*v_hat_denominator);
	p1d = (dvhatden_dp*dvhatden_dct)/(v_hat_denominator*v_hat_denominator);
	p1e = dvhatden_dpdct/v_hat_denominator;

	part1 =  p1a - p1b - p1c + 2e0*p1d - p1e;

	factor2 = (v_hat_denominator*dvhatnum_dct
		   - v_hat_numerator*dvhatden_dct)/
           (v_hat_denominator*dvhatnum_dsa - v_hat_numerator*dvhatden_dsa);

	p2a = dvhatnum_dpdsa/v_hat_numerator;
	p2b = (dvhatnum_dsa*dvhatden_dp)/(v_hat_numerator*v_hat_denominator);
	p2c = (dvhatnum_dp*dvhatden_dsa)/(v_hat_numerator*v_hat_denominator);
	p2d = (dvhatden_dp*dvhatden_dsa)/(v_hat_denominator*v_hat_denominator);
	p2e = dvhatden_dpdsa/v_hat_denominator;

	part2 =  p2a - p2b - p2c + 2e0*p2d - p2e;

	return ((part1 - factor2*part2)*rec_db2pa);
}
/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_internal_energy(sa,ct,p)  
!==========================================================================

!  Calculates internal energy of seawater using the computationally
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_internal_energy  :  internal_energy of seawater (48 term equation)
*/
double
gsw_internal_energy(double sa, double ct, double p)
{

	double	p0 = 101325e0, db2pa = 1e4;

	return (gsw_enthalpy(sa,ct,p) - (p0 + db2pa*p)*gsw_specvol(sa,ct,p));
}

/*
!==========================================================================
function gsw_enthalpy(sa,ct,p)  
!==========================================================================

!  Calculates specific enthalpy of seawater using the computationally-
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_enthalpy  :  specific enthalpy of seawater (48 term equation)
*/
double
gsw_enthalpy(double sa, double ct, double p)
{
	double	v01 =  9.998420897506056e+2, v02 =  2.839940833161907e0,
		v03 = -3.147759265588511e-2, v04 =  1.181805545074306e-3,
		v05 = -6.698001071123802e0,  v06 = -2.986498947203215e-2,
		v07 =  2.327859407479162e-4, v08 = -3.988822378968490e-2,
		v09 =  5.095422573880500e-4, v10 = -1.426984671633621e-5,
		v11 =  1.645039373682922e-7, v12 = -2.233269627352527e-2,
		v13 = -3.436090079851880e-4, v14 =  3.726050720345733e-6,
		v15 = -1.806789763745328e-4, v16 =  6.876837219536232e-7,
		v17 = -3.087032500374211e-7, v18 = -1.988366587925593e-8,
		v19 = -1.061519070296458e-11,v20 =  1.550932729220080e-10,
		v21 =  1.0e0,
		v22 =  2.775927747785646e-3, v23 = -2.349607444135925e-5,
		v24 =  1.119513357486743e-6, v25 =  6.743689325042773e-10,
		v26 = -7.521448093615448e-3, v27 = -2.764306979894411e-5,
		v28 =  1.262937315098546e-7, v29 =  9.527875081696435e-10,
		v30 = -1.811147201949891e-11,v31 = -3.303308871386421e-5,
		v32 =  3.801564588876298e-7, v33 = -7.672876869259043e-9,
		v34 = -4.634182341116144e-11,v35 =  2.681097235569143e-12,
		v36 =  5.419326551148740e-6, v37 = -2.742185394906099e-5,
		v38 = -3.212746477974189e-7, v39 =  3.191413910561627e-9,
		v40 = -1.931012931541776e-12,v41 = -1.105097577149576e-7,
		v42 =  6.211426728363857e-10,v43 = -1.119011592875110e-10,
		v44 = -1.941660213148725e-11,v45 = -1.864826425365600e-14,
		v46 =  1.119522344879478e-14,v47 = -1.200507748551599e-15,
		v48 =  6.057902487546866e-17;

	double	db2pa, cp0, sqrtsa, a0, a1, a2, a3, b0, b1, b2, b1sq;
	double	sqrt_disc, ca, cb, cn, cm, part;

	db2pa	= 1e4;			/* factor to convert from dbar to Pa */
	cp0	= 3991.86795711963e0;   /* from Eqn. (3.3.3) of IOC et al.
					   (2010) */

	sqrtsa	= sqrt(sa);

	a0	= v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct)))
		+ sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa
		+ sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))));
 
	a1	= v37 + ct*(v38 + ct*(v39 + v40*ct)) + sa*(v41 + v42*ct);

	a2	= v43 + ct*(v44 + v45*ct + v46*sa);

	a3	= v47 + v48*ct;

	b0	= v01 + ct*(v02 + ct*(v03 + v04*ct))
		+ sa*(v05 + ct*(v06 + v07*ct)
		+ sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))));
 
	b1	= 0.5e0*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct));

	b2	= v17 + ct*(v18 + v19*ct) + v20*sa;

	b1sq	= b1*b1;
	sqrt_disc	=
		  sqrt(b1sq - b0*b2);

	cn	= a0 + (2.0*a3*b0*b1/b2 - a2*b0)/b2;

	cm	= a1 + (4.0*a3*b1sq/b2 - a3*b0 - 2*a2*b1)/b2;

	ca	= b1 - sqrt_disc;
	cb	= b1 + sqrt_disc;

	part	= (cn*b2 - cm*b1)/(b2*(cb - ca));

	return (cp0*ct +
		db2pa*(p*(a2 - 2.0*a3*b1/b2 + 0.5*a3*p)/b2 +
		(cm/(2.0*b2))*log(1 + p*(2.0*b1 + b2*p)/b0) +
		part*log(1.0 + (b2*p*(cb - ca))/(ca*(cb + b2*p)))));
}

/*
!==========================================================================
function gsw_dynamic_enthalpy(sa,ct,p)  
!==========================================================================

!  Calculates dynamic enthalpy of seawater using the computationally-
!  efficient 48-term expression for density in terms of SA, CT and p
!  (IOC et al., 2010)
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_dynamic_enthalpy  :  dynamic enthalpy of seawater (48 term equation)
*/
double
gsw_dynamic_enthalpy(double sa, double ct, double p)
{
	double	v01 =  9.998420897506056e+2, v02 =  2.839940833161907e0,
		v03 = -3.147759265588511e-2, v04 =  1.181805545074306e-3,
		v05 = -6.698001071123802e0,  v06 = -2.986498947203215e-2,
		v07 =  2.327859407479162e-4, v08 = -3.988822378968490e-2,
		v09 =  5.095422573880500e-4, v10 = -1.426984671633621e-5,
		v11 =  1.645039373682922e-7, v12 = -2.233269627352527e-2,
		v13 = -3.436090079851880e-4, v14 =  3.726050720345733e-6,
		v15 = -1.806789763745328e-4, v16 =  6.876837219536232e-7,
		v17 = -3.087032500374211e-7, v18 = -1.988366587925593e-8,
		v19 = -1.061519070296458e-11,v20 =  1.550932729220080e-10,
		v21 =  1.0e0,
		v22 =  2.775927747785646e-3, v23 = -2.349607444135925e-5,
		v24 =  1.119513357486743e-6, v25 =  6.743689325042773e-10,
		v26 = -7.521448093615448e-3, v27 = -2.764306979894411e-5,
		v28 =  1.262937315098546e-7, v29 =  9.527875081696435e-10,
		v30 = -1.811147201949891e-11,v31 = -3.303308871386421e-5,
		v32 =  3.801564588876298e-7, v33 = -7.672876869259043e-9,
		v34 = -4.634182341116144e-11,v35 =  2.681097235569143e-12,
		v36 =  5.419326551148740e-6, v37 = -2.742185394906099e-5,
		v38 = -3.212746477974189e-7, v39 =  3.191413910561627e-9,
		v40 = -1.931012931541776e-12,v41 = -1.105097577149576e-7,
		v42 =  6.211426728363857e-10,v43 = -1.119011592875110e-10,
		v44 = -1.941660213148725e-11,v45 = -1.864826425365600e-14,
		v46 =  1.119522344879478e-14,v47 = -1.200507748551599e-15,
		v48 =  6.057902487546866e-17;

	double	db2pa, sqrtsa, a0, a1, a2, a3, b0, b1, b2, b1sq;
	double	sqrt_disc, ca, cb, cn, cm, part;

	db2pa	= 1e4;		/* factor to convert from dbar to Pa */

	sqrtsa	= sqrt(sa);

	a0	= v21 + ct*(v22 + ct*(v23 + ct*(v24 + v25*ct)))
		+ sa*(v26 + ct*(v27 + ct*(v28 + ct*(v29 + v30*ct))) + v36*sa 
		+ sqrtsa*(v31 + ct*(v32 + ct*(v33 + ct*(v34 + v35*ct)))));
 
	a1	= v37 + ct*(v38 + ct*(v39 + v40*ct)) + sa*(v41 + v42*ct);

	a2	= v43 + ct*(v44 + v45*ct + v46*sa);

	a3	= v47 + v48*ct;

	b0	= v01 + ct*(v02 + ct*(v03 + v04*ct))
		+ sa*(v05 + ct*(v06 + v07*ct)
		+ sqrtsa*(v08 + ct*(v09 + ct*(v10 + v11*ct))));
 
	b1	= 0.5e0*(v12 + ct*(v13 + v14*ct) + sa*(v15 + v16*ct));

	b2	= v17 + ct*(v18 + v19*ct) + v20*sa;

	b1sq	= b1*b1;
	sqrt_disc	=
		  sqrt(b1sq - b0*b2);

	cn	= a0 + (2*a3*b0*b1/b2 - a2*b0)/b2;

	cm	= a1 + (4*a3*b1sq/b2 - a3*b0 - 2*a2*b1)/b2;

	ca	= b1 - sqrt_disc;
	cb	= b1 + sqrt_disc;

	part	= (cn*b2 - cm*b1)/(b2*(cb - ca));

	return (db2pa*(p*(a2 - 2.0*a3*b1/b2 + 0.5*a3*p)/b2
		+ (cm/(2.0*b2))*log(1.0 + p*(2.0*b1 + b2*p)/b0)
		+ part*log(1.0 + (b2*p*(cb - ca))/(ca*(cb + b2*p)))));
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_sa_from_rho(rho,ct,p)
!==========================================================================

!  Calculates the Absolute Salinity of a seawater sample, for given values
!  of its density, Conservative Temperature and sea pressure (in dbar).
!  This function uses the computationally-efficient 48-term expression for
!  density in terms of SA, CT and p (IOC et al., 2010).

!  rho =  density of a seawater sample (e.g. 1026 kg/m^3).       [ kg/m^3 ]
!   Note. This input has not had 1000 kg/m^3 subtracted from it.
!     That is, it is 'density', not 'density anomaly'.
!  ct  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!
!  sa  =  Absolute Salinity                                          [g/kg]
*/
double
gsw_sa_from_rho(double rho, double ct, double p)
{
	int	no_iter;

	double	sa, v_lab, v_0, v_50, v_sa,
		sa_old, delta_v, sa_mean, alpha, beta;

	v_lab	= 1e0/rho;
	v_0	= gsw_specvol(0e0,ct,p);
	v_50	= gsw_specvol(50e0,ct,p);

	sa	= 50e0*(v_lab - v_0)/(v_50 - v_0);
	if (sa < 0e0 || sa > 50e0)
	    sa	= GSW_INVALID_VALUE;

	v_sa	= (v_50 - v_0)/50e0;

	for (no_iter=1; no_iter <= 2; no_iter++) {
	    sa_old	= sa;
	    delta_v	= gsw_specvol(sa_old,ct,p) - v_lab;
	    sa		= sa_old - delta_v/v_sa;
	    sa_mean	= 0.5e0*(sa + sa_old);
	    alpha	= gsw_alpha(sa_mean,ct,p);
	    beta	= gsw_beta(sa_mean,ct,p);
	    v_sa	= - beta/rho;
	    sa		= sa_old - delta_v/v_sa;
	    if (sa < 0e0 || sa > 50e0)
	 	sa	= GSW_INVALID_VALUE;
	}

	return (sa);

}

/*
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! water column properties, based on the 48-term expression for density
!--------------------------------------------------------------------------

!==========================================================================
subroutine gsw_nsquared(sa,ct,p,lat,nz,n2,p_mid)
!==========================================================================

!  Calculates the buoyancy frequency squared (N^2)(i.e. the Brunt-Vaisala 
!  frequency squared) at the mid pressure from the equation,
!
!
!           2      2             beta.d(SA) - alpha.d(CT)
!         N   =  g  .rho_local. -------------------------
!                                          dP
!
!  The pressure increment, dP, in the above formula is in Pa, so that it is
!  10^4 times the pressure increment dp in dbar. 
!
!  Note. This routine uses rho from "gsw_rho", which is the computationally
!  efficient 48-term expression for density in terms of SA, CT and p.  The    
!  48-term equation has been fitted in a restricted range of parameter 
!  space, and is most accurate inside the "oceanographic funnel" described 
!  in IOC et al. (2010).  The GSW library function "gsw_infunnel(SA,CT,p)" 
!  is avaialble to be used if one wants to test if some of one's data lies
!  outside this "funnel".  

! sa     : Absolute Salinity         (a profile (length nz))     [g/kg]
! ct     : Conservative Temperature  (a profile (length nz))     [deg C]
! p      : sea pressure              (a profile (length nz))     [dbar]
! lat    : latitude                  (a profile (length nz))     [deg N]                            
! nz     : number of bottles                             
! n2     : Brunt-Vaisala Frequency squared  (length nz-1)        [s^-2]
! p_mid  : Mid pressure between p grid      (length nz-1)        [dbar]
*/
void
gsw_nsquared(double *sa, double *ct, double *p, double *lat, int nz,
	double *n2, double *p_mid)
{
	int	k;

	double	db2pa = 1e4;
	double	grav_local, dsa, sa_mid, dct, ct_mid, dp, rho_mid;
	double	alpha_mid, beta_mid;

	for (k = 0; k < nz-1; k++) {
	    grav_local = 0.5*(gsw_grav(lat[k],p[k]) + gsw_grav(lat[k+1],p[k+1]));

	    dsa = (sa[k+1] - sa[k]);
	    sa_mid = 0.5*(sa[k] + sa[k+1]);
	    dct = (ct[k+1] - ct[k]);
	    ct_mid = 0.5*(ct[k] + ct[k+1]);
	    dp = (p[k+1] - p[k]);
	    p_mid[k] = 0.5*(p[k] + p[k+1]);

	    rho_mid = gsw_rho(sa_mid,ct_mid,p_mid[k]);
	    alpha_mid = gsw_alpha(sa_mid,ct_mid,p_mid[k]);
	    beta_mid = gsw_beta(sa_mid,ct_mid,p_mid[k]);

	    n2[k] = (grav_local*grav_local)*(rho_mid/(db2pa*dp))*
			(beta_mid*dsa - alpha_mid*dct);
	}
}

/*
!--------------------------------------------------------------------------

!==========================================================================
subroutine gsw_turner_rsubrho(sa,ct,p,nz,tu,rsubrho,p_mid)
!==========================================================================

!  Calculates the Turner angle and the Rsubrho as a function of pressure 
!  down a vertical water column.  These quantities express the relative 
!  contributions of the vertical gradients of Conservative Temperature 
!  and Absolute Salinity to the vertical stability (the square of the 
!  Brunt-Vaisala Frequency squared, N^2).  Tu and Rsubrho are evaluated at 
!  the mid pressure between the individual data points in the vertical.  
!  This function uses computationally-efficient 48-term expression for 
!  density in terms of SA, CT and p (IOC et al., 2010).  Note that 
!  in the double-diffusive literature, papers concerned with the 
!  "diffusive" form of double-diffusive convection often define the 
!  stability ratio as the reciprocal of what is defined here as the 
!  stability ratio.  
!
!  Note. The 48-term equation has been fitted in a restricted range of parameter 
!  space, and is most accurate inside the "oceanographic funnel" described 
!  in IOC et al. (2010).  The GSW library function "gsw_infunnel(SA,CT,p)" 
!  is avaialble to be used if one wants to test if some of one's data lies
!  outside this "funnel".  

! sa      : Absolute Salinity         (a profile (length nz))     [g/kg]
! ct      : Conservative Temperature  (a profile (length nz))     [deg C]
! p       : sea pressure              (a profile (length nz))     [dbar]                            
! nz      : number of bottles                             
! tu      : Turner angle, on the same (nz-1) grid as p_mid.
!           Turner angle has units of:           [ degrees of rotation ]
! rsubrho : Stability Ratio, on the same (nz-1) grid as p_mid.
!           Rsubrho is dimensionless.                       [ unitless ]
! p_mid   : Mid pressure between p grid  (length nz-1)           [dbar]
*/
void
gsw_turner_rsubrho(double *sa, double *ct, double *p, int nz,
	double *tu, double *rsubrho, double *p_mid)
{
	int	k;

	double	pi = 3.141592653589793e0;
	double	dsa, sa_mid, dct, ct_mid, dp;
	double	alpha_mid, beta_mid;

	for (k = 0; k < nz-1; k++) {
	    dsa = (sa[k] - sa[k+1]);
	    sa_mid = 0.5e0*(sa[k] + sa[k+1]);
	    dct = (ct[k] - ct[k+1]);
	    ct_mid = 0.5e0*(ct[k] + ct[k+1]);
	    dp = (p[k] - p[k+1]);
	    p_mid[k] = 0.5e0*(p[k] + p[k+1]);

	    alpha_mid = gsw_alpha(sa_mid,ct_mid,p_mid[k]);
	    beta_mid = gsw_beta(sa_mid,ct_mid,p_mid[k]);

	    tu[k] = (180e0/pi)*atan2((alpha_mid*dct + beta_mid*dsa),
				(alpha_mid*dct - beta_mid*dsa));

	    if (dsa == 0e0)
		rsubrho[k] = GSW_INVALID_VALUE;
	    else 
		rsubrho[k] = (alpha_mid*dct)/(beta_mid*dsa);
	}
}

/*
!--------------------------------------------------------------------------

!==========================================================================
subroutine gsw_ipv_vs_fnsquared_ratio(sa,ct,p,nz,ipv_vs_fnsquared_ratio,p_mid)
!==========================================================================

!  Calculates the ratio of the vertical gradient of potential density to 
!  the vertical gradient of locally-referenced potential density.  This 
!  ratio is also the ratio of the planetary Isopycnal Potential Vorticity
!  (IPV) to f times N^2, hence the name for this variable,
!  IPV_vs_fNsquared_ratio (see Eqn. (3.20.5) of IOC et al. (2010)). 
!  The reference sea pressure, p_ref, of the potential density surface must
!  have a constant value.
!
!  IPV_vs_fNsquared_ratio is evaluated at the mid pressure between the 
!  individual data points in the vertical.  This function uses the 
!  computationally-efficient 48-term expression for density in terms of 
!  SA, CT and p (IOC et al., 2010). 
!  Note. The 48-term equation has been fitted in a restricted range of parameter 
!  space, and is most accurate inside the "oceanographic funnel" described 
!  in IOC et al. (2010).  The GSW library function "gsw_infunnel(SA,CT,p)" 
!  is avaialble to be used if one wants to test if some of one's data lies
!  outside this "funnel".  

! sa      : Absolute Salinity         (a profile (length nz))     [g/kg]
! ct      : Conservative Temperature  (a profile (length nz))     [deg C]
! p       : sea pressure              (a profile (length nz))     [dbar]                            
! nz      : number of bottles                             
! IPV_vs_fNsquared_ratio
!         : The ratio of the vertical gradient of potential density
!           referenced to p_ref, to the vertical gradient of locally-
!           referenced potential density.  It is ouput on the same
!           vertical (M-1)xN grid as p_mid. 
!           IPV_vs_fNsquared_ratio is dimensionless.          [ unitless ]
! p_mid   : Mid pressure between p grid  (length nz-1)           [dbar]
*/
void
gsw_ipv_vs_fnsquared_ratio(double *sa, double *ct, double *p, int nz,
	double *ipv_vs_fnsquared_ratio, double *p_mid)
{
	int	k;

	double	dsa, sa_mid, dct, ct_mid, dp, p_ref;
	double	alpha_mid, beta_mid;
	double	alpha_pref, beta_pref, numerator, denominator;

	for (k = 0; k < nz-1; k++) {
	    dsa = (sa[k+1] - sa[k]);
	    sa_mid = 0.5*(sa[k] + sa[k+1]);
	    dct = (ct[k+1] - ct[k]);
	    ct_mid = 0.5*(ct[k] + ct[k+1]);
	    dp = (p[k+1] - p[k]);
	    p_mid[k] = 0.5*(p[k] + p[k+1]);

	    alpha_mid = gsw_alpha(sa_mid,ct_mid,p_mid[k]);
	    beta_mid = gsw_beta(sa_mid,ct_mid,p_mid[k]);
	    alpha_pref = gsw_alpha(sa_mid,ct_mid,p_ref);
	    beta_pref = gsw_beta(sa_mid,ct_mid,p_ref);

	    numerator = dct*alpha_pref - dsa*beta_pref;
	    denominator = dct*alpha_mid - dsa*beta_mid;

	    if (denominator == 0)
		ipv_vs_fnsquared_ratio[k] = GSW_INVALID_VALUE;
	    else
		ipv_vs_fnsquared_ratio[k] = numerator/denominator;
	}
}
/*
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! freezing temperatures
!--------------------------------------------------------------------------

!==========================================================================
function gsw_ct_freezing(sa,p,saturation_fraction)  
!==========================================================================

! Calculates the Conservative Temperature at which of seawater freezes 
! from Absolute Salinity and pressure.
!
! sa     : Absolute Salinity                                 [g/kg]
! p      : sea pressure                                      [dbar]
! saturation_fraction : saturation fraction
!
! gsw_ct_freezing : Conservative Temperature freezing point  [deg C]
*/
double
gsw_ct_freezing(double sa, double p, double saturation_fraction)
{
	double	c0  = 0.017947064327968736e0,c1  = -6.076099099929818e0,
		c2  = 4.883198653547851e0,   c3  = -11.88081601230542e0,
		c4  = 13.34658511480257e0,   c5  = -8.722761043208607e0,
		c6  = 2.082038908808201e0,   c7  = -7.389420998107497e0,
		c8  = -2.110913185058476e0,  c9  = 0.2295491578006229e0,
		c10 = -0.9891538123307282e0, c11 = -0.08987150128406496e0,
		c12 = 0.3831132432071728e0,  c13 = 1.054318231187074e0,
		c14 = 1.065556599652796e0,   c15 = -0.7997496801694032e0,
		c16 = 0.3850133554097069e0,  c17 = -2.078616693017569e0,
		c18 = 0.8756340772729538e0,  c19 = -2.079022768390933e0,
		c20 = 1.596435439942262e0,   c21 = 0.1338002171109174e0,
		c22 = 1.242891021876471e0;

	double	ct_freezing;
	double	sa_r, x, p_r, a, b;

	sa_r	= sa*1e-2;
	x	= sqrt(sa_r);
	p_r	= p*1e-4;

	ct_freezing	=
		c0 + sa_r*(c1 + x*(c2 + x*(c3 + x*(c4 + x*(c5 + c6*x)))))
		+ p_r*(c7 + p_r*(c8 + c9*p_r))
		+ sa_r*p_r*(c10 + p_r*(c12 + p_r*(c15 + c21*sa_r))
				+ sa_r*(c13 + c17*p_r + c19*sa_r)
 				+ x*(c11 + p_r*(c14 + c18*p_r)
				+ sa_r*(c16 + c20*p_r + c22*sa_r)));
    /*
    ** Adjust for the effects of dissolved air 
    */
	a	= 0.014289763856964e0;	/* Note that a =
						0.502500117621/35.16504. */
	b	= 0.057000649899720e0;
	ct_freezing	=
		ct_freezing - saturation_fraction*(1e-3)
		*(2.4e0 - a*sa)*(1e0 + b*(1e0 - sa/35.16504e0));

	if (p > 10000e0  ||  sa > 120e0  ||
	    (p+sa*71.428571428571402e0) > 13571.42857142857e0)
	    ct_freezing	= GSW_INVALID_VALUE;

	return (ct_freezing);
}

/*
!==========================================================================
function gsw_t_freezing(sa,p,saturation_fraction)  
!==========================================================================

! Calculates the in-situ temperature at which of seawater freezes 
! from Absolute Salinity and pressure.
!
! sa     : Absolute Salinity                                 [g/kg]
! p      : sea pressure                                      [dbar]
! saturation_fraction : saturation fraction
!
! gsw_t_freezing : in-situ temperature freezing point  [deg C]
*/
double
gsw_t_freezing(double sa, double p, double saturation_fraction)
{
	double	ct_freezing;
	double	t_freezing;

	ct_freezing	= gsw_ct_freezing(sa,p,saturation_fraction);
	t_freezing	= gsw_t_from_ct(sa,ct_freezing,p);

	if (ct_freezing > 9e10)
	    t_freezing	= GSW_INVALID_VALUE;

	return (t_freezing);
}

/*
!--------------------------------------------------------------------------
! isobaric melting enthalpy and isobaric evaporation enthalpy
!--------------------------------------------------------------------------

!==========================================================================
function gsw_latentheat_melting(sa,p)  
!==========================================================================

! Calculates latent heat, or enthalpy, of melting.
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! 
! gsw_latentheat_melting : latent heat of melting          [kg/m^3]
*/
double
gsw_latentheat_melting(double sa, double p)
{
	double	c0  =  3.334265169240710e5, c1  = -2.789444646733159e0,
		c2  = -1.822150156453350e4, c3  = -4.984585692734338e3,
		c4  = -7.371966528571920e1, c5  = -7.605802553358546e3,
		c6  =  1.195857305019339e3, c7  =  1.233720336206392e3,
		c8  =  2.294798676591890e2, c9  =  9.655751370889338e2,
		c10 = -5.792068522727968e2, c11 = -1.649446955902331e3,
		c12 = -1.029021448430547e3, c13 = -3.171558017172501e2,
		c14 = -1.751401389905041e2, c15 =  6.836527214265952e2,
		c16 =  1.078283734113611e3, c17 =  5.613896351265648e2,
		c18 =  6.968934948667265e2, c19 =  1.793032021946783e2,
		c20 =  8.692558481134256e1, c21 = -2.371103254714944e2,
		c22 = -5.775033277201674e2, c23 = -3.019749254648732e2,
		c24 = -6.420420579160927e2, c25 = -2.657570848596042e2,
		c26 = -1.646738151143109e1, c27 =  4.618228988300871e0;

	double	s_u, x, y;

	s_u	= 40e0*(35.16504e0/35e0);
	x	= sqrt(sa/s_u);
	y	= p*1e-4;

	return (c0 + x*(c1 + c4*y + x*(c3
		+ y*(c7 + c12*y) + x*(c6 + y*(c11 + y*(c17 + c24*y))
		+ x*(c10  + y*(c16 + c23*y) + x*(c15 + c22*y + c21*x)))))
		+ y*(c2 + y*(c5 + c8*x + y*(c9 + x*(c13 + c18*x)
		+ y*(c14 + x*(c19 + c25*x) + y*(c20 + c26*x + c27*y))))));
}

/*
!==========================================================================
function gsw_latentheat_evap_ct(sa,ct)  
!==========================================================================

! Calculates latent heat, or enthalpy, of evaporation.
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! gsw_latentheat_evaporation : latent heat of evaporation  [J/kg]
*/
double
gsw_latentheat_evap_ct(double sa, double ct)
{
	double	c0  =  2.499065844825125e6, c1  = -1.544590633515099e-1,
		c2  = -9.096800915831875e4, c3  =  1.665513670736000e2,
		c4  =  4.589984751248335e1, c5  =  1.894281502222415e1,
		c6  =  1.192559661490269e3, c7  = -6.631757848479068e3,
		c8  = -1.104989199195898e2, c9  = -1.207006482532330e3,
		c10 = -3.148710097513822e3, c11 =  7.437431482069087e2,
		c12 =  2.519335841663499e3, c13 =  1.186568375570869e1,
		c14 =  5.731307337366114e2, c15 =  1.213387273240204e3,
		c16 =  1.062383995581363e3, c17 = -6.399956483223386e2,
		c18 = -1.541083032068263e3, c19 =  8.460780175632090e1,
		c20 = -3.233571307223379e2, c21 = -2.031538422351553e2,
		c22 =  4.351585544019463e1, c23 = -8.062279018001309e2,
		c24 =  7.510134932437941e2, c25 =  1.797443329095446e2,
		c26 = -2.389853928747630e1, c27 =  1.021046205356775e2;

	double	s_u, x, y;

	s_u	= 40e0*(35.16504e0/35e0);
	x	= sqrt(sa/s_u);
	y	= ct/40;

	return (c0 + x*(c1 + c4*y + x*(c3
		+ y*(c7 + c12*y) + x*(c6 + y*(c11 + y*(c17 + c24*y))
		+ x*(c10 + y*(c16 + c23*y) + x*(c15 + c22*y + c21*x)))))
		+ y*(c2 + y*(c5 + c8*x + y*(c9 + x*(c13 + c18*x)
		+ y*(c14 + x*(c19 + c25*x) + y*(c20 + c26*x + c27*y))))));
}

/*
!==========================================================================
function gsw_latentheat_evap_t(sa,t)  
!==========================================================================

! Calculates latent heat, or enthalpy, of evaporation.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! 
! gsw_latentheat_evap_t : latent heat of evaporation       [J/kg]
*/
double
gsw_latentheat_evap_t(double sa, double t)
{

	double	ct = gsw_ct_from_pt(sa,t);

	return (gsw_latentheat_evap_ct(sa,ct));
}

/*
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! planet Earth properties
!--------------------------------------------------------------------------

!==========================================================================
function gsw_grav(lat,p)
!==========================================================================

! Calculates acceleration due to gravity as a function of latitude and as
!  a function of pressure in the ocean.
!
! lat  =  latitude in decimal degress north                [ -90 ... +90 ]
!  p  =  sea pressure                                              [ dbar ]
!
! gsw_grav : grav  =  gravitational acceleration               [ m s^-2 ]
*/
double
gsw_grav(double lat, double p)
{
	double	pi = 3.141592653589793e0;

	double	gamma, deg2rad, x, sin2, gs, z;

	gamma	= 2.26e-7;
	deg2rad	= pi/180e0;
	x	= sin(lat*deg2rad);  /* convert to radians */
	sin2	= x*x;
	gs	= 9.780327e0*(1.0e0 + (5.2792e-3 + (2.32e-5*sin2))*sin2);

	z	= gsw_z_from_p(p,lat);

	return (gs*(1 - gamma*z));	/* z is the height corresponding to p.
					   Note. In the ocean z is negative. */

}

/*
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! basic thermodynamic properties in terms of in-situ t, based on the exact Gibbs function
!--------------------------------------------------------------------------

!==========================================================================
function gsw_rho_t_exact(sa,t,p)  
!==========================================================================

! Calculates in-situ density of seawater from Absolute Salinity and 
! in-situ temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_rho_t_exact : in-situ density                        [kg/m^3]
*/
double
gsw_rho_t_exact(double sa, double t, double p)
{
	int	n0, n1;

	n0	= 0;
	n1	= 1;

	return (1.e0/gsw_gibbs(n0,n0,n1,sa,t,p));
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_pot_rho_t_exact(sa,t,p,p_ref)  
!==========================================================================

! Calculates the potential density of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! p_ref  : reference sea pressure                          [dbar]
! 
! gsw_pot_rho_t_exact : potential density                  [kg/m^3]
*/
double
gsw_pot_rho_t_exact(double sa, double t, double p, double p_ref)
{
	double	pt = gsw_pt_from_t(sa,t,p,p_ref);

	return (gsw_rho_t_exact(sa,pt,p_ref));
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_alpha_wrt_t_exact(sa,t,p)  
!==========================================================================

! Calculates thermal expansion coefficient of seawater with respect to 
! in-situ temperature
!
! sa     : Absolute Salinity                               [g/kg]
! t      : insitu temperature                              [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_alpha_wrt_t_exact : thermal expansion coefficient    [1/K]
!                         wrt (in-situ) temperature
*/
double
gsw_alpha_wrt_t_exact(double sa, double t, double p)
{
	int	n0, n1;

	n0	= 0;
	n1	= 1;

	return (gsw_gibbs(n0,n1,n1,sa,t,p)/gsw_gibbs(n0,n0,n1,sa,t,p));
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_beta_const_t_exact(sa,t,p)  
!==========================================================================

! Calculates saline (haline) contraction coefficient of seawater at 
! constant in-situ temperature.
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_beta_const_t_exact : haline contraction coefficient  [kg/g]
*/
double
gsw_beta_const_t_exact(double sa, double t, double p)
{
	int	n0, n1;

	n0	= 0;
	n1	= 1;

	return (-gsw_gibbs(n1,n0,n1,sa,t,p)/gsw_gibbs(n0,n0,n1,sa,t,p));
}

/*
!==========================================================================
function gsw_specvol_t_exact(sa,t,p)  
!==========================================================================

! Calculates the specific volume of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_specvol_t_exact : specific volume                    [kg/m^3]
*/
double
gsw_specvol_t_exact(double sa, double t, double p)
{
	int	n0, n1;

	n0	= 0;
	n1	= 1;

	return (gsw_gibbs(n0,n0,n1,sa,t,p));
}

/*
!==========================================================================
function gsw_sound_speed_t_exact(sa,t,p)  
!==========================================================================

! Calculates the sound speed of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_sound_speed_t_exact : sound speed                    [m/s]
*/
double
gsw_sound_speed_t_exact(double sa, double t, double p)
{
	int	n0, n1, n2;
	double	g_tt, g_tp;

	n0	= 0;
	n1	= 1;
	n2	= 2;

	g_tt	= gsw_gibbs(n0,n2,n0,sa,t,p);
	g_tp	= gsw_gibbs(n0,n1,n1,sa,t,p);

	return (gsw_gibbs(n0,n0,n1,sa,t,p) *
		sqrt(g_tt/(g_tp*g_tp - g_tt*gsw_gibbs(n0,n0,n2,sa,t,p))));
}

/*
!==========================================================================
function gsw_kappa_t_exact(sa,t,p)  
!==========================================================================

! isentropic compressibility of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_kappa_t_exact : isentropic compressibility           [1/Pa]
*/
double
gsw_kappa_t_exact(double sa, double t, double p)
{
	int	n0, n1, n2;
	double	g_tt, g_tp;

	n0	= 0;
	n1	= 1;
	n2	= 2;

	g_tt	= gsw_gibbs(n0,n2,n0,sa,t,p);
	g_tp	= gsw_gibbs(n0,n1,n1,sa,t,p);

	return ((g_tp*g_tp - g_tt*gsw_gibbs(n0,n0,n2,sa,t,p)) /
		(gsw_gibbs(n0,n0,n1,sa,t,p)*g_tt));
}

/*
!==========================================================================
function gsw_enthalpy_t_exact(sa,t,p)  
!==========================================================================

! Calculates the specific enthalpy of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_enthalpy_t_exact : specific enthalpy                 [J/kg]
*/
double
gsw_enthalpy_t_exact(double sa, double t, double p)
{
	int	n0, n1;

	n0	= 0;
	n1	= 1;

	return (gsw_gibbs(n0,n0,n0,sa,t,p) -
		(t+273.15e0)*gsw_gibbs(n0,n1,n0,sa,t,p));
}

/*
!==========================================================================
function gsw_cp_t_exact(sa,t,p)
!==========================================================================

! Calculates isobaric heat capacity of seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_cp_t_exact : heat capacity                           [J/(kg K)]
*/
double
gsw_cp_t_exact(double sa, double t, double p)
{
	int	n0, n2;

	n0 = 0;
	n2 = 2;

	return (-(t+273.15e0)*gsw_gibbs(n0,n2,n0,sa,t,p));
}

/*
!--------------------------------------------------------------------------
! Library functions of the GSW toolbox
!--------------------------------------------------------------------------

!==========================================================================
function gsw_gibbs(ns,nt,np,sa,t,p)
!==========================================================================

! seawater specific Gibbs free energy and derivatives up to order 2
!
! ns     : order of s derivative
! nt     : order of t derivative
! np     : order of p derivative
! sa     : Absolute Salinity                               [g/kg]
! t      : temperature                                     [deg C]
! p      : sea pressure                                    [dbar]
! 								-1
! gsw_gibbs  : specific Gibbs energy or its derivative	   [J kg  ]
*/
double
gsw_gibbs(int ns, int nt, int np, double sa, double t, double p)
{
	double	sfac, x2, x, y, z, g03, g08, return_value = 0.0;

	sfac	= 0.0248826675584615e0;

	x2	= sfac*sa;
	x	= sqrt(x2);
	y	= t*0.025e0;
	z	= p*1e-4;

	if (ns == 0  && nt == 0  && np == 0) {
	    g03	= 101.342743139674e0 + z*(100015.695367145e0 +
		z*(-2544.5765420363e0 + z*(284.517778446287e0 +
		z*(-33.3146754253611e0 + (4.20263108803084e0 -
		   0.546428511471039e0*z)*z)))) +
		y*(5.90578347909402e0 + z*(-270.983805184062e0 +
		z*(776.153611613101e0 + z*(-196.51255088122e0 +
		   (28.9796526294175e0 - 2.13290083518327e0*z)*z))) +
		y*(-12357.785933039e0 + z*(1455.0364540468e0 +
		z*(-756.558385769359e0 + z*(273.479662323528e0 +
		   z*(-55.5604063817218e0 + 4.34420671917197e0*z)))) +
		y*(736.741204151612e0 + z*(-672.50778314507e0 +
		z*(499.360390819152e0 + z*(-239.545330654412e0 +
		   (48.8012518593872e0 - 1.66307106208905e0*z)*z))) +
		y*(-148.185936433658e0 + z*(397.968445406972e0 +
		z*(-301.815380621876e0 + (152.196371733841e0 -
		   26.3748377232802e0*z)*z)) +
		y*(58.0259125842571e0 + z*(-194.618310617595e0 +
		z*(120.520654902025e0 + z*(-55.2723052340152e0 +
		   6.48190668077221e0*z))) +
		y*(-18.9843846514172e0 + y*(3.05081646487967e0 -
		   9.63108119393062e0*z) +
		z*(63.5113936641785e0 + z*(-22.2897317140459e0 +
		   8.17060541818112e0*z))))))));
          
	    g08	= x2*(1416.27648484197e0 + z*(-3310.49154044839e0 +
		z*(384.794152978599e0 + z*(-96.5324320107458e0 +
		   (15.8408172766824e0 - 2.62480156590992e0*z)*z))) +
		x*(-2432.14662381794e0 + x*(2025.80115603697e0 +
		y*(543.835333000098e0 + y*(-68.5572509204491e0 +
		y*(49.3667694856254e0 + y*(-17.1397577419788e0 +
		   2.49697009569508e0*y))) - 22.6683558512829e0*z) +
		x*(-1091.66841042967e0 - 196.028306689776e0*y +
		x*(374.60123787784e0 - 48.5891069025409e0*x +
		   36.7571622995805e0*y) + 36.0284195611086e0*z) +
		z*(-54.7919133532887e0 + (-4.08193978912261e0 -
		   30.1755111971161e0*z)*z)) +
		z*(199.459603073901e0 + z*(-52.2940909281335e0 +
		   (68.0444942726459e0 - 3.41251932441282e0*z)*z)) +
		y*(-493.407510141682e0 + z*(-175.292041186547e0 +
		   (83.1923927801819e0 - 29.483064349429e0*z)*z) +
		y*(-43.0664675978042e0 + z*(383.058066002476e0 +
		   z*(-54.1917262517112e0 + 25.6398487389914e0*z)) +
		y*(-10.0227370861875e0 - 460.319931801257e0*z +
		   y*(0.875600661808945e0 + 234.565187611355e0*z))))) +
		y*(168.072408311545e0 + z*(729.116529735046e0 +
		z*(-343.956902961561e0 + z*(124.687671116248e0 +
		   z*(-31.656964386073e0 + 7.04658803315449e0*z)))) +
		y*(880.031352997204e0 + y*(-225.267649263401e0 +
		y*(91.4260447751259e0 + y*(-21.6603240875311e0 +
		   2.13016970847183e0*y) +
		z*(-297.728741987187e0 + (74.726141138756e0 -
		   36.4872919001588e0*z)*z)) +
		z*(694.244814133268e0 + z*(-204.889641964903e0 +
		   (113.561697840594e0 - 11.1282734326413e0*z)*z))) +
		z*(-860.764303783977e0 + z*(337.409530269367e0 +
		z*(-178.314556207638e0 + (44.2040358308e0 -
		   7.92001547211682e0*z)*z))))));
        
	    if (sa > 0.e0)
		g08	= g08 + x2*(5812.81456626732e0 +
			  851.226734946706e0*y)*log(x);

	    return_value	= g03 + g08;
  
	} else if (ns == 1  && nt == 0  && np == 0) {
        
	    g08	= 8645.36753595126e0 + z*(-6620.98308089678e0 +
		z*(769.588305957198e0 + z*(-193.0648640214916e0 +
		   (31.6816345533648e0 - 5.24960313181984e0*z)*z))) +
		x*(-7296.43987145382e0 + x*(8103.20462414788e0 +
		y*(2175.341332000392e0 + y*(-274.2290036817964e0 +
		y*(197.4670779425016e0 + y*(-68.5590309679152e0 +
		   9.98788038278032e0*y))) - 90.6734234051316e0*z) +
		x*(-5458.34205214835e0 - 980.14153344888e0*y +
		x*(2247.60742726704e0 - 340.1237483177863e0*x +
		   220.542973797483e0*y) + 180.142097805543e0*z) +
		z*(-219.1676534131548e0 + (-16.32775915649044e0 -
		   120.7020447884644e0*z)*z)) +
		z*(598.378809221703e0 + z*(-156.8822727844005e0 +
		   (204.1334828179377e0 - 10.23755797323846e0*z)*z)) +
		y*(-1480.222530425046e0 + z*(-525.876123559641e0 +
		   (249.57717834054571e0 - 88.449193048287e0*z)*z) +
		y*(-129.1994027934126e0 + z*(1149.174198007428e0 +
		   z*(-162.5751787551336e0 + 76.9195462169742e0*z)) +
		y*(-30.0682112585625e0 - 1380.9597954037708e0*z +
		   y*(2.626801985426835e0 + 703.695562834065e0*z))))) +
		y*(1187.3715515697959e0 + z*(1458.233059470092e0 +
		z*(-687.913805923122e0 + z*(249.375342232496e0 +
		   z*(-63.313928772146e0 + 14.09317606630898e0*z)))) +
		y*(1760.062705994408e0 + y*(-450.535298526802e0 +
		y*(182.8520895502518e0 + y*(-43.3206481750622e0 +
		   4.26033941694366e0*y) +
		z*(-595.457483974374e0 + (149.452282277512e0 -
		   72.9745838003176e0*z)*z)) +
		z*(1388.489628266536e0 + z*(-409.779283929806e0 +
		   (227.123395681188e0 - 22.2565468652826e0*z)*z))) +
		z*(-1721.528607567954e0 + z*(674.819060538734e0 +
		z*(-356.629112415276e0 + (88.4080716616e0 -
		   15.84003094423364e0*z)*z)))));
  
	    if (sa > 0.e0)
		g08	= g08 + (11625.62913253464e0 + 1702.453469893412e0*y)*
			  log(x);
	    else
		g08 = 0.e0;
  
	    return_value	= 0.5*sfac*g08;

	} else if (ns == 0  && nt == 1  && np == 0) {
               
	    g03	= 5.90578347909402e0 + z*(-270.983805184062e0 +
		z*(776.153611613101e0 + z*(-196.51255088122e0 +
		   (28.9796526294175e0 - 2.13290083518327e0*z)*z))) +
		y*(-24715.571866078e0 + z*(2910.0729080936e0 +
		z*(-1513.116771538718e0 + z*(546.959324647056e0 +
		   z*(-111.1208127634436e0 + 8.68841343834394e0*z)))) +
		y*(2210.2236124548363e0 + z*(-2017.52334943521e0 +
		z*(1498.081172457456e0 + z*(-718.6359919632359e0 +
		   (146.4037555781616e0 - 4.9892131862671505e0*z)*z))) +
		y*(-592.743745734632e0 + z*(1591.873781627888e0 +
		z*(-1207.261522487504e0 + (608.785486935364e0 -
		   105.4993508931208e0*z)*z)) +
		y*(290.12956292128547e0 + z*(-973.091553087975e0 +
		z*(602.603274510125e0 + z*(-276.361526170076e0 +
		   32.40953340386105e0*z))) +
		y*(-113.90630790850321e0 + y*(21.35571525415769e0 -
		   67.41756835751434e0*z) +
		z*(381.06836198507096e0 + z*(-133.7383902842754e0 +
		   49.023632509086724e0*z)))))));
              
	    g08	= x2*(168.072408311545e0 + z*(729.116529735046e0 +
		z*(-343.956902961561e0 + z*(124.687671116248e0 +
		   z*(-31.656964386073e0 + 7.04658803315449e0*z)))) +
		x*(-493.407510141682e0 + x*(543.835333000098e0 +
		   x*(-196.028306689776e0 + 36.7571622995805e0*x) +
		y*(-137.1145018408982e0 + y*(148.10030845687618e0 +
		   y*(-68.5590309679152e0 + 12.4848504784754e0*y))) -
		   22.6683558512829e0*z) + z*(-175.292041186547e0 +
		   (83.1923927801819e0 - 29.483064349429e0*z)*z) +
		y*(-86.1329351956084e0 + z*(766.116132004952e0 +
		   z*(-108.3834525034224e0 + 51.2796974779828e0*z)) +
		y*(-30.0682112585625e0 - 1380.9597954037708e0*z +
		   y*(3.50240264723578e0 + 938.26075044542e0*z)))) +
		y*(1760.062705994408e0 + y*(-675.802947790203e0 +
		y*(365.7041791005036e0 + y*(-108.30162043765552e0 +
		   12.78101825083098e0*y) +
		z*(-1190.914967948748e0 + (298.904564555024e0 -
		   145.9491676006352e0*z)*z)) +
		z*(2082.7344423998043e0 + z*(-614.668925894709e0 +
		   (340.685093521782e0 - 33.3848202979239e0*z)*z))) +
		z*(-1721.528607567954e0 + z*(674.819060538734e0 +
		z*(-356.629112415276e0 + (88.4080716616e0 -
		   15.84003094423364e0*z)*z)))));
      
	    if (sa > 0.00)
		g08	= g08 + 851.226734946706e0*x2*log(x);
  
	    return_value	= (g03 + g08)*0.025e0;

	} else if (ns == 0  && nt == 0  && np == 1) {
    
	    g03	= 100015.695367145e0 + z*(-5089.1530840726e0 +
		z*(853.5533353388611e0 + z*(-133.2587017014444e0 +
		   (21.0131554401542e0 - 3.278571068826234e0*z)*z))) +
		y*(-270.983805184062e0 + z*(1552.307223226202e0 +
		z*(-589.53765264366e0 + (115.91861051767e0 -
		   10.664504175916349e0*z)*z)) +
		y*(1455.0364540468e0 + z*(-1513.116771538718e0 +
		z*(820.438986970584e0 + z*(-222.2416255268872e0 +
		   21.72103359585985e0*z))) +
		y*(-672.50778314507e0 + z*(998.720781638304e0 +
		z*(-718.6359919632359e0 + (195.2050074375488e0 -
		   8.31535531044525e0*z)*z)) +
		y*(397.968445406972e0 + z*(-603.630761243752e0 +
		   (456.589115201523e0 - 105.4993508931208e0*z)*z) +
		y*(-194.618310617595e0 + y*(63.5113936641785e0 -
		   9.63108119393062e0*y +
		z*(-44.5794634280918e0 + 24.511816254543362e0*z)) +
		z*(241.04130980405e0 + z*(-165.8169157020456e0 +
		25.92762672308884e0*z)))))));
  
	    g08	= x2*(-3310.49154044839e0 + z*(769.588305957198e0 +
		z*(-289.5972960322374e0 + (63.3632691067296e0 -
		   13.1240078295496e0*z)*z)) +
		x*(199.459603073901e0 + x*(-54.7919133532887e0 +
		   36.0284195611086e0*x - 22.6683558512829e0*y +
		(-8.16387957824522e0 - 90.52653359134831e0*z)*z) +
		z*(-104.588181856267e0 + (204.1334828179377e0 -
		   13.65007729765128e0*z)*z) +
		y*(-175.292041186547e0 + (166.3847855603638e0 -
		   88.449193048287e0*z)*z +
		y*(383.058066002476e0 + y*(-460.319931801257e0 +
		   234.565187611355e0*y) +
		z*(-108.3834525034224e0 + 76.9195462169742e0*z)))) +
		y*(729.116529735046e0 + z*(-687.913805923122e0 +
		z*(374.063013348744e0 + z*(-126.627857544292e0 +
		   35.23294016577245e0*z))) +
		y*(-860.764303783977e0 + y*(694.244814133268e0 +
		y*(-297.728741987187e0 + (149.452282277512e0 -
		   109.46187570047641e0*z)*z) +
		z*(-409.779283929806e0 + (340.685093521782e0 -
		   44.5130937305652e0*z)*z)) +
		z*(674.819060538734e0 + z*(-534.943668622914e0 +
		   (176.8161433232e0 - 39.600077360584095e0*z)*z)))));
     
	    return_value	= (g03 + g08)*1e-8;

	} else if (ns == 0  && nt == 2  && np == 0) {

	    g03	= -24715.571866078e0 + z*(2910.0729080936e0 + z*
		(-1513.116771538718e0 + z*(546.959324647056e0 +
		 z*(-111.1208127634436e0 + 8.68841343834394e0*z)))) +
		y*(4420.4472249096725e0 + z*(-4035.04669887042e0 +
		z*(2996.162344914912e0 + z*(-1437.2719839264719e0 +
		   (292.8075111563232e0 - 9.978426372534301e0*z)*z))) +
		y*(-1778.231237203896e0 + z*(4775.621344883664e0 +
		z*(-3621.784567462512e0 + (1826.356460806092e0 -
		   316.49805267936244e0*z)*z)) +
		y*(1160.5182516851419e0 + z*(-3892.3662123519e0 +
		z*(2410.4130980405e0 + z*(-1105.446104680304e0 +
		   129.6381336154442e0*z))) +
		y*(-569.531539542516e0 + y*(128.13429152494615e0 -
		   404.50541014508605e0*z) +
		z*(1905.341809925355e0 + z*(-668.691951421377e0 +
		   245.11816254543362e0*z))))));

	    g08	= x2*(1760.062705994408e0 + x*(-86.1329351956084e0 +
		x*(-137.1145018408982e0 + y*(296.20061691375236e0 +
		   y*(-205.67709290374563e0 + 49.9394019139016e0*y))) +
		z*(766.116132004952e0 + z*(-108.3834525034224e0 +
		   51.2796974779828e0*z)) +
		y*(-60.136422517125e0 - 2761.9195908075417e0*z +
		   y*(10.50720794170734e0 + 2814.78225133626e0*z))) +
		y*(-1351.605895580406e0 + y*(1097.1125373015109e0 +
		   y*(-433.20648175062206e0 + 63.905091254154904e0*y) +
		z*(-3572.7449038462437e0 + (896.713693665072e0 -
		   437.84750280190565e0*z)*z)) +
		z*(4165.4688847996085e0 + z*(-1229.337851789418e0 +
		   (681.370187043564e0 - 66.7696405958478e0*z)*z))) +
		z*(-1721.528607567954e0 + z*(674.819060538734e0 +
		z*(-356.629112415276e0 + (88.4080716616e0 -
		   15.84003094423364e0*z)*z))));
     
	    return_value	= (g03 + g08)*0.000625e0  ;

	} else if (ns == 1  && nt == 0  && np == 1) {

	    g08	=     -6620.98308089678e0 + z*(1539.176611914396e0 +
		z*(-579.1945920644748e0 + (126.7265382134592e0 -
		   26.2480156590992e0*z)*z)) +
		x*(598.378809221703e0 + x*(-219.1676534131548e0 +
		   180.142097805543e0*x - 90.6734234051316e0*y +
		(-32.65551831298088e0 - 362.10613436539325e0*z)*z) +
		z*(-313.764545568801e0 + (612.4004484538132e0 -
		   40.95023189295384e0*z)*z) +
		y*(-525.876123559641e0 + (499.15435668109143e0 -
		   265.347579144861e0*z)*z +
		y*(1149.174198007428e0 + y*(-1380.9597954037708e0 +
		   703.695562834065e0*y) +
		z*(-325.1503575102672e0 + 230.7586386509226e0*z)))) +
		y*(1458.233059470092e0 + z*(-1375.827611846244e0 +
		z*(748.126026697488e0 + z*(-253.255715088584e0 +
		   70.4658803315449e0*z))) +
		y*(-1721.528607567954e0 + y*(1388.489628266536e0 +
		y*(-595.457483974374e0 + (298.904564555024e0 -
		   218.92375140095282e0*z)*z) +
		z*(-819.558567859612e0 + (681.370187043564e0 -
		   89.0261874611304e0*z)*z)) +
		z*(1349.638121077468e0 + z*(-1069.887337245828e0 +
		   (353.6322866464e0 - 79.20015472116819e0*z)*z))));    

	    return_value	= g08*sfac*0.5e-8;

	} else if (ns == 0  && nt == 1  && np == 1) {

	    g03	= -270.983805184062e0 + z*(1552.307223226202e0 +
		z*(-589.53765264366e0 + (115.91861051767e0 -
		   10.664504175916349e0*z)*z)) +
		y*(2910.0729080936e0 + z*(-3026.233543077436e0 +
		z*(1640.877973941168e0 + z*(-444.4832510537744e0 +
		   43.4420671917197e0*z))) +
		y*(-2017.52334943521e0 + z*(2996.162344914912e0 +
		z*(-2155.907975889708e0 + (585.6150223126464e0 -
		   24.946065931335752e0*z)*z)) +
		y*(1591.873781627888e0 + z*(-2414.523044975008e0 +
		   (1826.356460806092e0 - 421.9974035724832e0*z)*z) +
		y*(-973.091553087975e0 + z*(1205.20654902025e0 +
		   z*(-829.084578510228e0 + 129.6381336154442e0*z)) +
		y*(381.06836198507096e0 - 67.41756835751434e0*y +
		   z*(-267.4767805685508e0 + 147.07089752726017e0*z))))));
    
	    g08	= x2*(729.116529735046e0 + z*(-687.913805923122e0 +
		z*(374.063013348744e0 + z*(-126.627857544292e0 +
		   35.23294016577245e0*z))) +
		x*(-175.292041186547e0 - 22.6683558512829e0*x +
		   (166.3847855603638e0 - 88.449193048287e0*z)*z +
		y*(766.116132004952e0 + y*(-1380.9597954037708e0 +
		   938.26075044542e0*y) +
		z*(-216.7669050068448e0 + 153.8390924339484e0*z))) +
		y*(-1721.528607567954e0 + y*(2082.7344423998043e0 +
		y*(-1190.914967948748e0 + (597.809129110048e0 -
		   437.84750280190565e0*z)*z) +
		z*(-1229.337851789418e0 + (1022.055280565346e0 -
		   133.5392811916956e0*z)*z)) +
		z*(1349.638121077468e0 + z*(-1069.887337245828e0 +
		   (353.6322866464e0 - 79.20015472116819e0*z)*z))));
    
	    return_value	= (g03 + g08)*2.5e-10;

	} else if (ns == 0  && nt == 0  && np == 2) {
           
	    g03	= -5089.1530840726e0 + z*(1707.1066706777221e0 +
		z*(-399.7761051043332e0 + (84.0526217606168e0 -
		   16.39285534413117e0*z)*z)) +
		y*(1552.307223226202e0 + z*(-1179.07530528732e0 +
		   (347.75583155301e0 - 42.658016703665396e0*z)*z) +
		y*(-1513.116771538718e0 + z*(1640.877973941168e0 +
		   z*(-666.7248765806615e0 + 86.8841343834394e0*z)) +
		y*(998.720781638304e0 + z*(-1437.2719839264719e0 +
		   (585.6150223126464e0 - 33.261421241781e0*z)*z) +
		y*(-603.630761243752e0 + (913.178230403046e0 -
		   316.49805267936244e0*z)*z +
		y*(241.04130980405e0 + y*(-44.5794634280918e0 +
		   49.023632509086724e0*z) +
		z*(-331.6338314040912e0 + 77.78288016926652e0*z))))));
            
	    g08	= x2*(769.588305957198e0 + z*(-579.1945920644748e0 +
		     (190.08980732018878e0 - 52.4960313181984e0*z)*z) +
		x*(-104.588181856267e0 + x*(-8.16387957824522e0 -
		   181.05306718269662e0*z) +
		(408.2669656358754e0 - 40.95023189295384e0*z)*z +
		y*(166.3847855603638e0 - 176.898386096574e0*z +
		   y*(-108.3834525034224e0 + 153.8390924339484e0*z))) +
		y*(-687.913805923122e0 + z*(748.126026697488e0 +
		   z*(-379.883572632876e0 + 140.9317606630898e0*z)) +
		y*(674.819060538734e0 + z*(-1069.887337245828e0 +
		   (530.4484299696e0 - 158.40030944233638e0*z)*z) +
		y*(-409.779283929806e0 + y*(149.452282277512e0 -
		   218.92375140095282e0*z) +
		(681.370187043564e0 - 133.5392811916956e0*z)*z))));
    
	    return_value	= (g03 + g08)*1e-16 ;

	}

	return (return_value);
}

/*
!==========================================================================
subroutine gsw_add_barrier(input_data,lon,lat,long_grid,lat_grid,dlong_grid,dlat_grid,output_data)
!==========================================================================

!  Adds a barrier through Central America (Panama) and then averages
!  over the appropriate side of the barrier
! 
!  data_in      :  data                                         [unitless]
!  lon          :  Longitudes of data degrees east              [0 ... +360]
!  lat          :  Latitudes of data degrees north              [-90 ... +90]
!  longs_grid   :  Longitudes of regular grid degrees east      [0 ... +360]
!  lats_grid    :  Latitudes of regular grid degrees north      [-90 ... +90]
!  dlongs_grid  :  Longitude difference of regular grid degrees [deg longitude]
!  dlats_grid   :  Latitude difference of regular grid degrees  [deg latitude]
!
! gsw_add_barrier  : average of data depending on which side of the 
!                    Panama cannal it is on                                 [unitless]
*/
void
gsw_add_barrier(double *input_data, double lon, double lat,
		double long_grid, double lat_grid, double dlong_grid,
		double dlat_grid, double *output_data)
{
	int	above_line[4];
	int	k, nmean, above_line0, kk;
	double	longs_pan[6] = {260.0, 272.59, 276.50, 278.65, 280.73, 292.0},
		lats_pan[6] = {19.55, 13.97, 9.6, 8.1, 9.33, 3.4};
	double	r, lats_line, data_mean;

	k	= gsw_indx(longs_pan,6,lon);	/*   the lon/lat point */
	r	= (lon-longs_pan[k])/(longs_pan[k+1]-longs_pan[k]);
	lats_line	= lats_pan[k] + r*(lats_pan[k+1]-lats_pan[k]);

	above_line0	= (lats_line <= lat);

	k	= gsw_indx(longs_pan,6,long_grid);/*the 1 & 4 lon/lat points*/ 
	r	= (long_grid-longs_pan[k])/(longs_pan[k+1]-longs_pan[k]);
	lats_line	= lats_pan[k] + r*(lats_pan[k+1]-lats_pan[k]);

	above_line[0]	= (lats_line <= lat_grid);
	above_line[3]	= (lats_line <= lat_grid+dlat_grid);

	k		= gsw_indx(longs_pan,6,long_grid+dlong_grid);
			/*the 2 & 3 lon/lat points */
	r		= (long_grid+dlong_grid-longs_pan[k])/
			(longs_pan[k+1]-longs_pan[k]);
	lats_line	= lats_pan[k] + r*(lats_pan[k+1]-lats_pan[k]);

	above_line[1]	= (lats_line <= lat_grid);
	above_line[2]	= (lats_line <= lat_grid+dlat_grid);

	nmean		= 0;
	data_mean	= 0.0;

	for (kk=0; kk<4; kk++) {
	    if ((fabs(input_data[kk]) <= 100.0) &&
		above_line0 == above_line[kk]) {
		nmean	= nmean+1;
		data_mean	= data_mean+input_data[kk];
	    }
	}
	if (nmean == 0)
	    data_mean	= 0.0;	/*errorreturn*/
	else
	    data_mean	= data_mean/nmean;

	for (kk=0; kk<4; kk++) {
	    if ((fabs(input_data[kk]) >= 1e10) ||
		above_line0 != above_line[kk])
		output_data[kk]	= data_mean;
	    else
		output_data[kk]	= input_data[kk];
	}

	return;
}

/*
!==========================================================================
subroutine gsw_add_mean(data_in,lon,lat,data_out)
!==========================================================================

! Replaces NaN's with non-nan mean of the 4 adjacent neighbours
!
! data_in   : data set of the 4 adjacent neighbours   
! lon      : longitude
! lat       : latitude
!
! data_out : non-nan mean of the 4 adjacent neighbours     [unitless]
*/
void
gsw_add_mean(double *data_in, double lon, double lat, double *data_out)
{
	int	k, nmean;
	double	data_mean;

	nmean		= 0;
	data_mean	= 0.0;

	for (k=0; k<4; k++) {
	    if (fabs(data_in[k]) <= 100.0) {
		nmean++;
		data_mean	= data_mean+data_in[k];
	    }
	}

	if (nmean == 0)
	    data_mean	= 0.0;    /*errorreturn*/
	else
	    data_mean	= data_mean/nmean;

	for (k=0; k<4; k++) {
	    if (fabs(data_in[k]) >= 100.0)
		data_out[k]	= data_mean;
	    else
		data_out[k]	= data_in[k];
	}
	return;
}

/*
!==========================================================================
function gsw_xinterp1(x,y,n,x0)
!==========================================================================

! Linearly interpolate a real array   
!
! x      : y array (Must be monotonic)               
! y      : y array     
! n      : length of X and Y arrays
! x0     : value to be interpolated
!
! gsw_xinterp1 : Linearly interpolated value
*/
double
gsw_xinterp1(double *x, double *y, int n, double x0)
{
	int	k;
	double	r;

	k	= gsw_indx(x,n,x0);
	r	= (x0-x[k])/(x[k+1]-x[k]);
	return (y[k] + r*(y[k+1]-y[k]));
}

/*
!==========================================================================
subroutine gsw_indx(x,n,z,k)
!==========================================================================

!  Finds the index of the value in a monotonically increasing array
!
!  x	 :  array of monotonically increasing values
!  n     :  length of the array
!  z     :  value to be indexed
!
!  K      : index K - if X(K) <= Z < X(K+1), or
!  N-1     		    - if Z = X(N)
!
*/
int
gsw_indx(double *x, int n, double z)
{
	int	k, ku, kl, km;

	if (z > x[0] && z < x[n-1]) {
	    kl	= 0;
	    ku	= n-1;
	    while (ku-kl > 1) {
		km	= (ku+kl)>>1;
		if (z > x[km])
		    kl	= km;
		else
		    ku	= km;
	    }
	    k	= kl;
	    if (z == x[k+1])
		k++;
	} else if (z <= x[0])
	    k	= 0;
	else if (z >= x[n-1])
	    k	= n-2;
	else {
	    fprintf(stderr, "ERROR in function gsw_indx : out of range\n");
	    fprintf(stderr,"z = %g, n = %d, x:\n", z, n);
	    for (kl=0; kl<n; kl++)
		fprintf(stderr,"x[%d] = %g\n", kl, x[kl]);
	    k	= 0;
	}

	return (k);
}

/*
!==========================================================================
function gsw_fdelta(p,lon,lat)
!==========================================================================

! Calculates fdelta. 
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_fdelta : Absolute Salinty Anomaly                    [unitless]
*/
double
gsw_fdelta(double p, double lon, double lat)
{
	double	sa, saar;

	saar	= gsw_saar(p,lon,lat);
	sa	= ((1.0 + 0.35)*saar)/(1.0 - 0.35*saar);
	if (saar > 1e10)
	    sa	= GSW_INVALID_VALUE;
	else
	    sa	= ((1.0 + 0.35)*saar)/(1.0 - 0.35*saar);
	return (sa);
}

/*
!==========================================================================
function gsw_sa_from_sp_baltic(sp,lon,lat)
!==========================================================================

! For the Baltic Sea, calculates Absolute Salinity with a value
! computed analytically from Practical Salinity
!
! sp     : Practical Salinity                              [unitless]
! lon    : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
! p      : sea pressure                                    [dbar]
!
! gsw_sa_from_sp_baltic : Absolute Salinity                [g/kg]
*/
double
gsw_sa_from_sp_baltic(double sp, double lon, double lat)
{
	double	xx_left, xx_right, return_value;
	double	xb_left[3]={12.6, 7.0, 26.0}, yb_left[3]={50.0, 59.0, 69.0},
		xb_right[2]={45.0, 26.0}, yb_right[2]={50.0, 69.0};

	if (xb_left[1] < lon  && lon < xb_right[0]  && yb_left[0] < lat  &&
	    lat < yb_left[2]) {
  
	    xx_left	= gsw_xinterp1(yb_left, xb_left, 3, lat);
    
	    xx_right	= gsw_xinterp1(yb_right, xb_right, 2, lat);
    
	    if (xx_left <= lon  && lon <= xx_right)
		return_value	=((35.16504 - 0.087)/35.0)*sp + 0.087;
	    else
		return_value	= GSW_INVALID_VALUE;
	} else
	    return_value	= GSW_INVALID_VALUE;

	return (return_value);
}

/*
!==========================================================================
function gsw_sp_from_sa_baltic(sa,lon,lat)
!==========================================================================

! For the Baltic Sea, calculates Practical Salinity with a value
! computed analytically from Absolute Salinity
!
! sa     : Absolute Salinity                               [g/kg]
! lon   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
! p      : sea pressure                                    [dbar]
!
! gsw_sp_from_sa_baltic  : Practical Salinity              [unitless]
*/
double
gsw_sp_from_sa_baltic(double sa, double lon, double lat)
{
	double	xx_left, xx_right, return_value;
	double	xb_left[3]={12.6, 7.0, 26.0}, yb_left[3]={50.0, 59.0, 69.0},
		xb_right[2]={45.0, 26.0}, yb_right[2]={50.0, 69.0};

	if (xb_left[1] < lon  && lon < xb_right[0]  && yb_left[0] < lat  &&
	    lat < yb_left[2]) {
  
	    xx_left	= gsw_xinterp1(yb_left, xb_left, 3, lat);
    
	    xx_right	= gsw_xinterp1(yb_right, xb_right, 2, lat);
    
	    if (xx_left <= lon  && lon <= xx_right)
		return_value	= (35.0/(35.16504 - 0.087))*(sa - 0.087);
	    else
		return_value	= GSW_INVALID_VALUE;
	} else
	    return_value	= GSW_INVALID_VALUE;

	return (return_value);
}
     
/*
!==========================================================================
function gsw_entropy_part(sa,t,p)
!==========================================================================

! entropy minus the terms that are a function of only SA
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_entropy_part : entropy part
*/
double
gsw_entropy_part(double sa, double t, double p)
{
	double	sfac, x2, x, y, z, g03, g08, return_value;

	sfac	= 0.0248826675584615;

	x2	= sfac*sa;
	x	= sqrt(x2);
	y	= t*0.025;
	z	= p*1e-4;

	g03	= z*(-270.983805184062e0 +
		z*(776.153611613101e0 + z*(-196.51255088122e0 +
		   (28.9796526294175e0 - 2.13290083518327e0*z)*z))) +
		y*(-24715.571866078e0 + z*(2910.0729080936e0 +
		z*(-1513.116771538718e0 + z*(546.959324647056e0 +
		   z*(-111.1208127634436e0 + 8.68841343834394e0*z)))) +
		y*(2210.2236124548363e0 + z*(-2017.52334943521e0 +
		z*(1498.081172457456e0 + z*(-718.6359919632359e0 +
		   (146.4037555781616e0 - 4.9892131862671505e0*z)*z))) +
		y*(-592.743745734632e0 + z*(1591.873781627888e0 +
		z*(-1207.261522487504e0 + (608.785486935364e0 -
		   105.4993508931208e0*z)*z)) +
		y*(290.12956292128547e0 + z*(-973.091553087975e0 +
		z*(602.603274510125e0 + z*(-276.361526170076e0 +
		   32.40953340386105e0*z))) +
		y*(-113.90630790850321e0 + y*(21.35571525415769e0 -
		   67.41756835751434e0*z) +
		z*(381.06836198507096e0 + z*(-133.7383902842754e0 +
		   49.023632509086724e0*z)))))));

	g08	= x2*(z*(729.116529735046e0 +
		z*(-343.956902961561e0 + z*(124.687671116248e0 +
		   z*(-31.656964386073e0 + 7.04658803315449e0*z)))) +
		x*( x*(y*(-137.1145018408982e0 + y*(148.10030845687618e0 +
		   y*(-68.5590309679152e0 + 12.4848504784754e0*y))) -
		22.6683558512829e0*z) + z*(-175.292041186547e0 +
		   (83.1923927801819e0 - 29.483064349429e0*z)*z) +
		y*(-86.1329351956084e0 + z*(766.116132004952e0 +
		   z*(-108.3834525034224e0 + 51.2796974779828e0*z)) +
		y*(-30.0682112585625e0 - 1380.9597954037708e0*z +
		   y*(3.50240264723578e0 + 938.26075044542e0*z)))) +
		y*(1760.062705994408e0 + y*(-675.802947790203e0 +
		y*(365.7041791005036e0 + y*(-108.30162043765552e0 +
		   12.78101825083098e0*y) +
		z*(-1190.914967948748e0 + (298.904564555024e0 -
		   145.9491676006352e0*z)*z)) +
		z*(2082.7344423998043e0 + z*(-614.668925894709e0 +
		   (340.685093521782e0 - 33.3848202979239e0*z)*z))) +
		z*(-1721.528607567954e0 + z*(674.819060538734e0 +
		z*(-356.629112415276e0 + (88.4080716616e0 -
		   15.84003094423364e0*z)*z)))));

	return_value	= -(g03 + g08)*0.025;

	return (return_value);
}

/*
!==========================================================================
function gsw_entropy_part_zerop(sa,pt0)
!==========================================================================

! entropy part evaluated at the sea surface
!
! sa     : Absolute Salinity                               [g/kg]
! pt0    : insitu temperature                              [deg C]
! 
! gsw_entropy_part_zerop : entropy part at the sea surface
*/
double
gsw_entropy_part_zerop(double sa, double pt0)
{
	double	sfac, x2, x, y, g03, g08, return_value;

	sfac	= 0.0248826675584615;

	x2	= sfac*sa;
	x	= sqrt(x2);
	y	= pt0*0.025;

	g03	= y*(-24715.571866078e0 + y*(2210.2236124548363e0 +
		y*(-592.743745734632e0 + y*(290.12956292128547e0 +
		y*(-113.90630790850321e0 + y*21.35571525415769e0)))));

	g08	= x2*(x*(x*(y*(-137.1145018408982e0 + y*(148.10030845687618e0 +
		y*(-68.5590309679152e0 + 12.4848504784754e0*y)))) +
		y*(-86.1329351956084e0 + y*(-30.0682112585625e0 +
		   y*3.50240264723578e0))) +
		y*(1760.062705994408e0 + y*(-675.802947790203e0 +
		y*(365.7041791005036e0 + y*(-108.30162043765552e0 +
		   12.78101825083098e0*y)))));

	return_value	= -(g03 + g08)*0.025;

	return (return_value);
}

/*
!==========================================================================
function gsw_gibbs_pt0_pt0(sa,pt0)
!==========================================================================

! gibbs_tt at (sa,pt,0)
!
! sa     : Absolute Salinity                            [g/kg]
! pt0    : potential temperature                        [deg C]
! 
! gsw_gibbs_pt0_pt0 : gibbs_tt at (sa,pt,0)
*/
double
gsw_gibbs_pt0_pt0(double sa, double pt0)
{
	double	sfac, x2, x, y, g03, g08, return_value;

	sfac	= 0.0248826675584615;

	x2	= sfac*sa;
	x	= sqrt(x2);
	y	= pt0*0.025;

	g03	= -24715.571866078e0 +
		y*(4420.4472249096725e0 +
		y*(-1778.231237203896e0 +
		y*(1160.5182516851419e0 +
		y*(-569.531539542516e0 + y*128.13429152494615e0))));

	g08	= x2*(1760.062705994408e0 + x*(-86.1329351956084e0 +
		x*(-137.1145018408982e0 + y*(296.20061691375236e0 +
		y*(-205.67709290374563e0 + 49.9394019139016e0*y))) +
		y*(-60.136422517125e0 + y*10.50720794170734e0)) +
		y*(-1351.605895580406e0 + y*(1097.1125373015109e0 +
		y*(-433.20648175062206e0 + 63.905091254154904e0*y))));

	return_value	= (g03 + g08)*0.000625;

	return (return_value);
}

/*
!==========================================================================
function gsw_specvol_sso_0_p(p) 
!==========================================================================

!  This function calculates specifc volume at the Standard Ocean Salinty,
!  SSO, and at a Conservative Temperature of zero degrees C, as a function 
!  of pressure, p, in dbar, using a streamlined version of the 48-term CT
!  version of specific volume, that is, a streamlined version of the code
!  "gsw_specvol(SA,CT,p)".
!
! p      : sea pressure                                    [dbar]
! 							     3   -1
! gsw_specvol_sso_0_p : specvol(sso,0,p)		   [m  kg  ]
*/
double
gsw_specvol_sso_0_p(double p)
{
	double	v01 =  9.998420897506056e+2, v05 = -6.698001071123802e0,
		v08 = -3.988822378968490e-2, v12 = -2.233269627352527e-2,
		v15 = -1.806789763745328e-4, v17 = -3.087032500374211e-7,
		v20 =  1.550932729220080e-10,v21 =  1.0e0,
		v26 = -7.521448093615448e-3, v31 = -3.303308871386421e-5,
		v36 =  5.419326551148740e-6, v37 = -2.742185394906099e-5,
		v41 = -1.105097577149576e-7, v43 = -1.119011592875110e-10,
		v47 = -1.200507748551599e-15;

	double	sso, sqrtsso, return_value;

	sso	= 35.16504;
	sqrtsso	= 5.930011804372737;      /*sqrt(SSO) = 5.930011804372737*/

	return_value	= (v21 + sso*(v26 + v36*sso + v31*sqrtsso) 
			+ p*(v37 + v41*sso + p*(v43 + v47*p )))/
			(v01 + sso*(v05 + v08*sqrtsso)
			+ p*(v12 + v15*sso + p*(v17 + v20*sso)));

	return (return_value);
}

/*
!--------------------------------------------------------------------------

!==========================================================================
function gsw_enthalpy_sso_0_p(p)
!==========================================================================

!  This function calculates enthalpy at the Standard Ocean Salinity, SSO,
!  and at a Conservative Temperature of zero degrees C, as a function of
!  pressure, p, in dbar, using a streamlined version of the 48-term CT
!  version of the Gibbs function, that is, a streamlined version of the
!  code "gsw_enthalpy(SA,CT,p)".
!
! p      : sea pressure                                    [dbar]
!
! gsw_enthalpy_sso_0_p : enthalpy(sso,0,p)
*/
double
gsw_enthalpy_sso_0_p(double p)
{
	double	v01 =  9.998420897506056e2, v05 = -6.698001071123802e0,
		v08 = -3.988822378968490e-2, v12 = -2.233269627352527e-2,
		v15 = -1.806789763745328e-4, v17 = -3.087032500374211e-7,
		v20 =  1.550932729220080e-10, v21 =  1.0e0,
		v26 = -7.521448093615448e-3, v31 = -3.303308871386421e-5,
		v36 =  5.419326551148740e-6, v37 = -2.742185394906099e-5,
		v41 = -1.105097577149576e-7, v43 = -1.119011592875110e-10,
		v47 = -1.200507748551599e-15, db2pa = 1e4,
		sso = 35.16504e0, sqrtsso = 5.930011804372737e0;

	double	a0, a1, a2, a3, b0, b1, b2, b1sq, sqrt_disc, n, m, a, b, part;

	a0	= v21 + sso*(v26 + v36*sso + v31*sqrtsso);

	a1	= v37 + v41*sso;

	a2	= v43;

	a3	= v47;

	b0	= v01 + sso*(v05 + v08*sqrtsso);

	b1	= 0.5*(v12 + v15*sso);

	b2	= v17 + v20*sso;

	b1sq	= b1*b1;
	sqrt_disc	= sqrt(b1sq - b0*b2);

	n	= a0 + (2e0*a3*b0*b1/b2 - a2*b0)/b2;

	m	= a1 + (4e0*a3*b1sq/b2 - a3*b0 - 2*a2*b1)/b2;

	a	= b1 - sqrt_disc;
	b	= b1 + sqrt_disc;

	part	= (n*b2 - m*b1)/(b2*(b - a));

	return (db2pa*(p*(a2 - 2e0*a3*b1/b2 + 0.5e0*a3*p)/b2 + 
		(m/(2e0*b2))*log(1e0 + p*(2e0*b1 + b2*p)/b0) + 
		part*log(1e0 + (b2*p*(b - a))/(a*(b + b2*p)))));

}

/*
!--------------------------------------------------------------------------

!==========================================================================
function  gsw_hill_ratio_at_sp2(t)
!==========================================================================

!  Calculates the Hill ratio, which is the adjustment needed to apply for
!  Practical Salinities smaller than 2.  This ratio is defined at a 
!  Practical Salinity = 2 and in-situ temperature, t using PSS-78. The Hill
!  ratio is the ratio of 2 to the output of the Hill et al. (1986) formula
!  for Practical Salinity at the conductivity ratio, Rt, at which Practical
!  Salinity on the PSS-78 scale is exactly 2.
*/
double
gsw_hill_ratio_at_sp2(double t)
{
	double	
		a0 =  0.0080e0, a1 = -0.1692e0, a2 = 25.3851e0,
		a3 = 14.0941e0, a4 = -7.0261e0, a5 =  2.7081e0,
		b0 =  0.0005e0, b1 = -0.0056e0, b2 = -0.0066e0,
		b3 = -0.0375e0, b4 =  0.0636e0, b5 = -0.0144e0,
		g0 = 2.641463563366498e-1, g1 = 2.007883247811176e-4,
		g2 = -4.107694432853053e-6, g3 = 8.401670882091225e-8,
		g4 = -1.711392021989210e-9, g5 = 3.374193893377380e-11,
		g6 = -5.923731174730784e-13, g7 = 8.057771569962299e-15,
		g8 = -7.054313817447962e-17, g9 = 2.859992717347235e-19,
		rk  =  0.0162e0, sp2 = 2e0;

	double	t68, ft68, rtx0, dsp_drtx, sp_est, rtx, rtxm, x, part1, part2;
	double	sqrty, sp_hill_raw_at_sp2;

	t68	= t*1.00024;
	ft68	= (t68 - 15.0)/(1.0 + rk*(t68 - 15.0));

    /*!------------------------------------------------------------------------
    **! Find the initial estimates of Rtx (Rtx0) and of the derivative dSP_dRtx
    **! at SP = 2. 
    **!------------------------------------------------------------------------
    */
	rtx0	= g0 + t68*(g1 + t68*(g2 + t68*(g3 + t68*(g4 + t68*(g5
		+ t68*(g6 + t68*(g7 + t68*(g8 + t68*g9))))))));
     
	dsp_drtx	=
		a1 + (2*a2 + (3*a3 + (4*a4 + 5*a5*rtx0)*rtx0)*rtx0)*rtx0 +
		ft68*(b1 + (2*b2 + (3*b3 + (4*b4 + 5*b5*rtx0)*rtx0)*rtx0)*rtx0);

    /*!-------------------------------------------------------------------------
    **! Begin a single modified Newton-Raphson iteration to find Rt at SP = 2.
    **!-------------------------------------------------------------------------
    */
	sp_est	= a0 + (a1 + (a2 + (a3 + (a4 + a5*rtx0)*rtx0)*rtx0)*rtx0)*rtx0
		+ ft68*(b0 + (b1 + (b2+ (b3 + (b4 + b5*rtx0)*rtx0)*rtx0)*
		  rtx0)*rtx0);
	rtx	= rtx0 - (sp_est - sp2)/dsp_drtx;
	rtxm	= 0.5*(rtx + rtx0);
	dsp_drtx= a1 + (2*a2 + (3*a3 + (4*a4 + 5*a5*rtxm)*rtxm)*rtxm)*rtxm
		+ ft68*(b1 + (2*b2 + (3*b3 + (4*b4 + 5*b5*rtxm)*
						rtxm)*rtxm)*rtxm);
	rtx	= rtx0 - (sp_est - sp2)/dsp_drtx;
    /*
    **! This is the end of one full iteration of the modified Newton-Raphson 
    **! iterative equation solver. The error in Rtx at this point is equivalent 
    **! to an error in SP of 9e-16 psu.
    */
                                
	x	= 400.0*rtx*rtx;
	sqrty	= 10.0*rtx;
	part1	= 1.0 + x*(1.5 + x);
	part2	= 1.0 + sqrty*(1.0 + sqrty*(1.0 + sqrty));
	sp_hill_raw_at_sp2	= sp2 - a0/part1 - b0*ft68/part2;

	return (2.0/sp_hill_raw_at_sp2);
}
/*
**  The End
**!==========================================================================
*/
