/*
!==========================================================================
function gsw_sstar_from_sa(sa,p,lon,lat)  
!==========================================================================

! Calculates Preformed Salinity, Sstar, from Absolute Salinity, SA. 
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sstar_from_sa : Preformed Salinity                   [g/kg]
*/
double
gsw_sstar_from_sa(double sa, double p, double lon, double lat)
{
	double	saar;

	saar	= gsw_saar(p,lon,lat);
    /*
	! In the Baltic Sea, Sstar = sa, and note that gsw_saar returns zero
	! for saar in the Baltic.
    */
	if (saar == GSW_INVALID_VALUE)
	    return (saar);
	return (sa*(1e0 - 0.35e0*saar)/(1e0 + saar));
}
