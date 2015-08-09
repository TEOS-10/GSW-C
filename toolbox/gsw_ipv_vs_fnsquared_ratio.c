/*
!==========================================================================
subroutine gsw_ipv_vs_fnsquared_ratio(sa,ct,p,pref,nz,ipv_vs_fnsquared_ratio,p_mid)
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
!  individual data points in the vertical.

! sa      : Absolute Salinity         (a profile (length nz))     [g/kg]
! ct      : Conservative Temperature  (a profile (length nz))     [deg C]
! p       : sea pressure              (a profile (length nz))     [dbar]
! p_ref   : reference sea pressure of the potential density surface
!        ( i.e. absolute reference pressure - 10.1325 dbar )      [dbar]
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
gsw_ipv_vs_fnsquared_ratio(double *sa, double *ct, double *p, double p_ref,
	int nz, double *ipv_vs_fnsquared_ratio, double *p_mid)
{
	int	k;
	double	dsa, sa_mid, dct, ct_mid;
	double	alpha_mid, beta_mid;
	double	alpha_pref, beta_pref, numerator, denominator;

        if (nz < 2) {
	    *p_mid = *ipv_vs_fnsquared_ratio = GSW_INVALID_VALUE;
	    return;
	}
	for (k = 0; k < nz-1; k++) {
	    dsa		= (sa[k] - sa[k+1]);
	    dct		= (ct[k] - ct[k+1]);
	    sa_mid	= 0.5*(sa[k] + sa[k+1]);
	    ct_mid	= 0.5*(ct[k] + ct[k+1]);
	    p_mid[k]	= 0.5*(p[k] + p[k+1]);

	    alpha_mid = gsw_alpha(sa_mid,ct_mid,p_mid[k]);
	    beta_mid = gsw_beta(sa_mid,ct_mid,p_mid[k]);
	    alpha_pref = gsw_alpha(sa_mid,ct_mid,p_ref);
	    beta_pref = gsw_beta(sa_mid,ct_mid,p_ref);

	    numerator = dct*alpha_pref - dsa*beta_pref;
	    denominator = dct*alpha_mid - dsa*beta_mid;

	    if (denominator == 0.0)
		ipv_vs_fnsquared_ratio[k] = GSW_INVALID_VALUE;
	    else
		ipv_vs_fnsquared_ratio[k] = numerator/denominator;
	}
}
