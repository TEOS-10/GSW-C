/*
!==========================================================================
elemental function gsw_ct_freezing (sa, p, saturation_fraction)
!==========================================================================
!
!  Calculates the Conservative Temperature at which seawater freezes.  The
!  Conservative Temperature freezing point is calculated from the exact
!  in-situ freezing temperature which is found by a modified Newton-Raphson
!  iteration (McDougall and Wotherspoon, 2013) of the equality of the
!  chemical potentials of water in seawater and in ice.
!
!  An alternative GSW function, gsw_CT_freezing_poly, it is based on a
!  computationally-efficient polynomial, and is accurate to within -5e-4 K
!  and 6e-4 K, when compared with this function.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  saturation_fraction = the saturation fraction of dissolved air in
!                        seawater
!
!  CT_freezing = Conservative Temperature at freezing of seawater [ deg C ]
!--------------------------------------------------------------------------
*/
double
gsw_ct_freezing(double sa, double p, double saturation_fraction)
{
	double	t_freezing;

	t_freezing = gsw_t_freezing(sa,p,saturation_fraction);
	return (gsw_ct_from_t(sa,t_freezing,p));
}
