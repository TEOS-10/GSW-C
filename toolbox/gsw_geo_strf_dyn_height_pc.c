/*
!==========================================================================
pure subroutine gsw_geo_strf_dyn_height_pc (sa, ct, delta_p, &
                                            geo_strf_dyn_height_pc, p_mid)
!==========================================================================
!
!  Calculates dynamic height anomaly as the integral of specific volume 
!  anomaly from the the sea surface pressure (0 Pa) to the pressure p.
!  This function, gsw_geo_strf_dyn_height_pc, is to used when the 
!  Absolute Salinity and Conservative Temperature are piecewise constant in 
!  the vertical over sucessive pressure intervals of delta_p (such as in
!  a forward "z-coordinate" ocean model).  "geo_strf_dyn_height_pc" is
!  the dynamic height anomaly with respect to the sea surface.  That is, 
!  "geo_strf_dyn_height_pc" is the geostrophic streamfunction for the 
!  difference between the horizontal velocity at the pressure concerned, p,
!  and the horizontal velocity at the sea surface.  Dynamic height anomaly 
!  is the geostrophic streamfunction in an isobaric surface.  The reference
!  values used for the specific volume anomaly are SA = SSO = 35.16504 g/kg
!  and CT = 0 deg C.  The output values of geo_strf_dyn_height_pc are 
!  given at the mid-point pressures, p_mid, of each layer in which SA and 
!  CT are vertically piecewice constant (pc).  This function calculates 
!  enthalpy using the computationally-efficient 75-term expression for 
!  specific volume of Roquet et al., (2015). 
!
!  SA       =  Absolute Salinity                                   [ g/kg ]
!  CT       =  Conservative Temperature (ITS-90)                  [ deg C ]
!  delta_p  =  difference in sea pressure between the deep and     [ dbar ]
!              shallow extents of each layer in which SA and CT
!              are vertically constant. delta_p must be positive.
!              
!  Note. sea pressure is absolute pressure minus 10.1325 dbar.
!
!  geo_strf_dyn_height_pc =  dynamic height anomaly             [ m^2/s^2 ]
!  p_mid                  =  mid-point pressure in each layer      [ dbar ]
!--------------------------------------------------------------------------
*/
double *
gsw_geo_strf_dyn_height_pc(double *sa, double *ct, double *delta_p, int n_levels,
	double *geo_strf_dyn_height_pc, double *p_mid)
{
	int	i, np;
	double	*delta_h, delta_h_half, dyn_height_deep=0.0,
		*p_deep, *p_shallow;

	for (i=0; i<n_levels; i++)
	    if (delta_p[i] < 0.0)
		return (NULL);

	np = n_levels;
	delta_h = malloc(3*np*sizeof (double));
	p_deep = delta_h+np; p_shallow = p_deep+np;

	for (i=0; i<np; i++) {
	    p_deep[i] = (i==0)? delta_p[0] : p_deep[i-1] + delta_p[i];
	    p_shallow[i] = p_deep[i] - delta_p[i];
	    delta_h[i] = gsw_enthalpy_diff(sa[i],ct[i],p_shallow[i],p_deep[i]);
	}

	for (i=0; i<np; i++) {
	    dyn_height_deep = dyn_height_deep - delta_h[i];
		/* This is Phi minus Phi_0 of Eqn. (3.32.2) of IOC et al. (2010).*/
	    p_mid[i] = 0.5*(p_shallow[i]  + p_deep[i]);
	    delta_h_half = gsw_enthalpy_diff(sa[i],ct[i],p_mid[i],p_deep[i]);

	    geo_strf_dyn_height_pc[i] = gsw_enthalpy_sso_0(p_mid[i]) + 
	                           dyn_height_deep + delta_h_half;
	}
	free(delta_h);
	return (geo_strf_dyn_height_pc);
}
