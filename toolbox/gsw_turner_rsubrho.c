/*
!==========================================================================
subroutine gsw_turner_rsubrho(sa,ct,p,nz,tu,rsubrho,p_mid)
!==========================================================================

!  Calculates the Turner angle and the Rsubrho as a function of pressure 
!  down a vertical water column.  These quantities express the relative 
!  contributions of the vertical gradients of Conservative Temperature 
!  and Absolute Salinity to the vertical stability (the square of the 
!  Brunt-Vaisala Frequency squared, N^2).  Tu and Rsubrho are evaluated at 
!  the mid pressure between the individual data points in the vertical.  
!
!  Note that in the double-diffusive literature, papers concerned with
!  the "diffusive" form of double-diffusive convection often define the 
!  stability ratio as the reciprocal of what is defined here as the 
!  stability ratio.  

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
	GSW_TEOS10_CONSTANTS;
	int	k;
	double	dsa, sa_mid, dct, ct_mid, alpha_mid, beta_mid;

	if (nz < 2)
	    return;

	for (k = 0; k < nz-1; k++) {
	    dsa		= (sa[k] - sa[k+1]);
	    sa_mid	= 0.5e0*(sa[k] + sa[k+1]);
	    dct		= (ct[k] - ct[k+1]);
	    ct_mid	= 0.5e0*(ct[k] + ct[k+1]);
	    p_mid[k]	= 0.5e0*(p[k] + p[k+1]);
	    gsw_specvol_alpha_beta(sa_mid,ct_mid,p_mid[k],NULL,&alpha_mid,
					&beta_mid);
	    tu[k] = rad2deg*atan2((alpha_mid*dct + beta_mid*dsa),
				(alpha_mid*dct - beta_mid*dsa));
	    if (dsa == 0.0)
		rsubrho[k] = GSW_INVALID_VALUE;
	    else 
		rsubrho[k] = (alpha_mid*dct)/(beta_mid*dsa);
	}
}
