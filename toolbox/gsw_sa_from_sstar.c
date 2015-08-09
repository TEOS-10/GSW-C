/*
!==========================================================================
function gsw_sa_from_sstar(sstar,p,lon,lat)  
!==========================================================================

! Calculates Absolute Salinity, SA, from Preformed Salinity, Sstar.
!
! Sstar  : Preformed Salinity                              [g/kg]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sa_from_sstar   : Absolute Salinity                  [g/kg]
*/
double
gsw_sa_from_sstar(double sstar, double p, double lon, double lat)
{
	double	saar;

	saar	= gsw_saar(p,lon,lat);
	if (saar == GSW_INVALID_VALUE)
	    return (saar);
    /*
    **! In the Baltic Sea, Sstar = SA, and note that gsw_saar returns zero
    **! for SAAR in the Baltic.
    */
	return (sstar*(1e0 + saar)/(1e0 - 0.35e0*saar));
}
