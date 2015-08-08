/*
!==========================================================================
function  gsw_hill_ratio_at_sp2(t)
!==========================================================================

!  Calculates the Hill ratio, which is the adjustment needed to apply for
!  Practical Salinities smaller than 2.  This ratio is defined at a 
!  Practical Salinity = 2 and in-situ temperature, t using PSS-78. The Hill
!  ratio is the ratio of 2 to the output of the Hill et al. (1986) formula
!  for Practical Salinity at the conductivity ratio, Rt, at which Practical
!  Salinity on the PSS-78 scale is exactly 2.
!
!  t                 : in-situ temperature (ITS-90)              [deg C]
!  hill_ratio_at_sp2 : Hill ratio                                [dimensionless]
*/
double
gsw_hill_ratio_at_sp2(double t)
{
	GSW_SP_COEFFICIENTS;
	double	g0 = 2.641463563366498e-1, g1 = 2.007883247811176e-4,
		g2 = -4.107694432853053e-6, g3 = 8.401670882091225e-8,
		g4 = -1.711392021989210e-9, g5 = 3.374193893377380e-11,
		g6 = -5.923731174730784e-13, g7 = 8.057771569962299e-15,
		g8 = -7.054313817447962e-17, g9 = 2.859992717347235e-19,
		sp2 = 2.0;
	double	t68, ft68, rtx0, dsp_drtx, sp_est, rtx, rtxm, x, part1, part2;
	double	sqrty, sp_hill_raw_at_sp2;

	t68	= t*1.00024;
	ft68	= (t68 - 15.0)/(1.0 + k*(t68 - 15.0));

    /*!------------------------------------------------------------------------
    **! Find the initial estimates of Rtx (Rtx0) and of the derivative dSP_dRtx
    **! at SP = 2. 
    **!------------------------------------------------------------------------
    */
	rtx0	= g0 + t68*(g1 + t68*(g2 + t68*(g3 + t68*(g4 + t68*(g5
		+ t68*(g6 + t68*(g7 + t68*(g8 + t68*g9))))))));
     
	dsp_drtx= a1 + (2*a2 + (3*a3 + (4*a4 + 5*a5*rtx0)*rtx0)*rtx0)*rtx0 +
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
	sp_hill_raw_at_sp2 = sp2 - a0/part1 - b0*ft68/part2;

	return (2.0/sp_hill_raw_at_sp2);
}
