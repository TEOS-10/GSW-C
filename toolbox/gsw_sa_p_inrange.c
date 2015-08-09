/*
!==========================================================================
elemental function gsw_sa_p_inrange (sa, p)
!==========================================================================
!
!  Check for any values that are out of the TEOS-10 range ...
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!---------------------------------------------------------------------------
*/
int
gsw_sa_p_inrange(double sa, double p)
{
	if (p > 10000.0 || sa > 120.0 ||
	    (p + sa*71.428571428571402) > 13571.42857142857)
	    return (0);
	return (1);
}
