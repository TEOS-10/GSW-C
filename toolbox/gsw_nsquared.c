/*
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
! sa     : Absolute Salinity         (a profile (length nz))     [g/kg]
! ct     : Conservative Temperature  (a profile (length nz))     [deg C]
! p      : sea pressure              (a profile (length nz))     [dbar]
! lat    : latitude                  (a profile (length nz))     [deg N] 
! nz     : number of levels in the profile               
! n2     : Brunt-Vaisala Frequency squared  (length nz-1)        [s^-2]
! p_mid  : Mid pressure between p grid      (length nz-1)        [dbar]
*/
void
gsw_nsquared(double *sa, double *ct, double *p, double *lat, int nz,
	double *n2, double *p_mid)
{
	GSW_TEOS10_CONSTANTS;
	int	k;
	double	p_grav, n_grav, grav_local, dsa, sa_mid, dct, ct_mid,
		dp, rho_mid, alpha_mid, beta_mid;

	if (nz < 2)
	    return;
	p_grav	= gsw_grav(lat[0],p[0]);
	for (k = 0; k < nz-1; k++) {
	    n_grav	= gsw_grav(lat[k+1],p[k+1]);
	    grav_local	= 0.5*(p_grav + n_grav);
	    dsa		= (sa[k+1] - sa[k]);
	    sa_mid	= 0.5*(sa[k] + sa[k+1]);
	    dct		= (ct[k+1] - ct[k]);
	    ct_mid	= 0.5*(ct[k] + ct[k+1]);
	    dp		= (p[k+1] - p[k]);
	    p_mid[k]	= 0.5*(p[k] + p[k+1]);
	    rho_mid	= gsw_rho(sa_mid,ct_mid,p_mid[k]);
	    alpha_mid	= gsw_alpha(sa_mid,ct_mid,p_mid[k]);
	    beta_mid	= gsw_beta(sa_mid,ct_mid,p_mid[k]);

	    n2[k]	= (grav_local*grav_local)*(rho_mid/(db2pa*dp))*
			  (beta_mid*dsa - alpha_mid*dct);
	    p_grav	= n_grav;
	}
}
