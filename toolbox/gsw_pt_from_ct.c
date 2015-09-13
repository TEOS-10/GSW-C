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
	GSW_TEOS10_CONSTANTS;
	double	a5ct, b3ct, ct_factor, pt_num, pt_recden, ct_diff;
	double	pt, pt_old, ptm, dpt_dct, s1;
	double	a0	= -1.446013646344788e-2,    
		a1	= -3.305308995852924e-3,    
		a2	=  1.062415929128982e-4,     
		a3	=  9.477566673794488e-1,     
		a4	=  2.166591947736613e-3,
		a5	=  3.828842955039902e-3,
		b0	=  1.000000000000000e0,
		b1	=  6.506097115635800e-4,
		b2	=  3.830289486850898e-3,
		b3	=  1.247811760368034e-6;

	s1	= sa/gsw_ups;

	a5ct	= a5*ct;
	b3ct	= b3*ct;

	ct_factor	= (a3 + a4*s1 + a5ct);
	pt_num		= a0 + s1*(a1 + a2*s1) + ct*ct_factor;
	pt_recden	= 1.0/(b0 + b1*s1 + ct*(b2 + b3ct));
	pt		= pt_num*pt_recden;

	dpt_dct	= pt_recden*(ct_factor + a5ct - (b2 + b3ct + b3ct)*pt);

    /*
    **  Start the 1.5 iterations through the modified Newton-Raphson
    **  iterative method.
    */

	ct_diff	= gsw_ct_from_pt(sa,pt) - ct;
	pt_old	= pt;
	pt	= pt_old - ct_diff*dpt_dct;
	ptm	= 0.5*(pt + pt_old);

	dpt_dct	= -gsw_cp0/((ptm + gsw_t0)*gsw_gibbs_pt0_pt0(sa,ptm));

	pt	= pt_old - ct_diff*dpt_dct;
	ct_diff	= gsw_ct_from_pt(sa,pt) - ct;
	pt_old	= pt;
	return (pt_old - ct_diff*dpt_dct);
}
