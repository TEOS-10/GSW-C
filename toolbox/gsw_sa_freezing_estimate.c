/*
!==========================================================================
elemental function gsw_sa_freezing_estimate (p, saturation_fraction, ct, t)
!==========================================================================
!
! Form an estimate of SA from a polynomial in CT and p 
!
!--------------------------------------------------------------------------
*/
double
gsw_sa_freezing_estimate(double p, double saturation_fraction, double *ct,
	double *t)
{
	GSW_TEOS10_CONSTANTS;
	double	ctx, ctsat, sa,
		/*note that aa = 0.502500117621d0/35.16504*/
		aa = 0.014289763856964,
		bb = 0.057000649899720,

		p0  =  2.570124672768757e-1,
		p1  = -1.917742353032266e1,
		p2  = -1.413382858617969e-2,
		p3  = -5.427484830917552e-1,
		p4  = -4.126621135193472e-4,
		p5  = -4.176407833276121e-7,
		p6  =  4.688217641883641e-5,
		p7  = -3.039808885885726e-8,
		p8  = -4.990118091261456e-11,
		p9  = -9.733920711119464e-9,
		p10 = -7.723324202726337e-12,
		p11 =  7.121854166249257e-16,
		p12 =  1.256474634100811e-12,
		p13 =  2.105103897918125e-15,
		p14 =  8.663811778227171e-19;

	/*A very rough estimate of sa to get the saturated ct*/
	if (ct != NULL) {
	    sa = max(-(*ct + 9e-4*p)/0.06, 0.0);
	    ctx = *ct;
	} else if (t != NULL) {
	    sa = max(-(*t + 9e-4*p)/0.06, 0.0);
	    ctx = gsw_ct_from_t(sa,*t,p);
	} else {
	    return (0.0);
	}
	/*
	! CTsat is the estimated value of CT if the seawater were saturated with
	! dissolved air, recognizing that it actually has the air fraction
	! saturation_fraction; see McDougall, Barker and Feistel, 2014).  
	*/
	ctsat = ctx - (1.0-saturation_fraction)*
	        (1e-3)*(2.4-aa*sa)*(1.0+bb*(1.0-sa/gsw_sso));

	return (p0 + p*(p2 + p4*ctsat + p*(p5 + ctsat*(p7 + p9*ctsat)
	    + p*(p8  + ctsat*(p10 + p12*ctsat) + p*(p11 + p13*ctsat + p14*p))))
	    + ctsat*(p1 + ctsat*(p3 + p6*p)));
}
