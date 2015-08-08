/*
!==========================================================================
elemental function gsw_melting_ice_sa_ct_ratio (sa, ct, p, t_ih)
!==========================================================================
!
!  Calculates the ratio of SA to CT changes when ice melts into seawater.
!  It is assumed that a small mass of ice melts into an infinite mass of
!  seawater.  Because of the infinite mass of seawater, the ice will always
!  melt.   
!
!  The output, melting_seaice_SA_CT_ratio, is dSA/dCT rather than dCT/dSA. 
!  This is done so that when SA = 0, the output, dSA/dCT is zero whereas 
!  dCT/dSA would be infinite. 
!
!  SA   =  Absolute Salinity of seawater                           [ g/kg ]
!  CT   =  Conservative Temperature of seawater (ITS-90)          [ deg C ]
!  p    =  sea pressure at which the melting occurs                [ dbar ]
!         ( i.e. absolute pressure - 10.1325d0 dbar ) 
!  t_Ih =  the in-situ temperature of the ice (ITS-90)            [ deg C ]
!
!  melting_ice_SA_CT_ratio = the ratio of SA to CT changes when ice melts
!                            into a large mass of seawater 
!                                                          [ g kg^-1 K^-1 ] 
!--------------------------------------------------------------------------
*/
double
gsw_melting_ice_sa_ct_ratio(double sa, double ct, double p, double t_ih)
{
	double	ctf, h, h_ih, tf, h_hat_sa, h_hat_ct;
	double	saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	if (ct < ctf) {
	    /*the seawater ct input is below the freezing temperature*/
	    return (GSW_INVALID_VALUE);
	}

	tf = gsw_t_freezing(0.0,p,saturation_fraction);
	if (t_ih > tf) {
	    /*t_ih exceeds the freezing temperature at sa = 0*/
	    return (GSW_INVALID_VALUE);
	}

	h = gsw_enthalpy_ct_exact(sa,ct,p);
	h_ih = gsw_enthalpy_ice(t_ih,p);
	gsw_enthalpy_first_derivatives_ct_exact(sa,ct,p,&h_hat_sa,&h_hat_ct);
	    /*Note that h_hat_CT is equal to cp0*(273.15 + t)/(273.15 + pt0)*/

	return (sa*h_hat_ct/(h - h_ih - sa*h_hat_sa));
}
