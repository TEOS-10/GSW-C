/*
!==========================================================================
elemental function gsw_t_freezing_poly (sa, p, saturation_fraction)
!==========================================================================
!
!  Calculates the in-situ temperature at which seawater freezes from a
!  computationally efficient polynomial.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  saturation_fraction = the saturation fraction of dissolved air in
!                        seawater
!
!  t_freezing = in-situ temperature at which seawater freezes.    [ deg C ]
!               (ITS-90)
!--------------------------------------------------------------------------
*/
double
gsw_t_freezing_poly(double sa, double p, double saturation_fraction)
{
	GSW_TEOS10_CONSTANTS;
	double ctf, return_value;

    ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
    return_value = gsw_t_from_ct(sa,ctf,p);
	return (return_value);
}
