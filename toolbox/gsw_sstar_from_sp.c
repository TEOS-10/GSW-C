/*
!==========================================================================
function gsw_sstar_from_sp(sp,p,lon,lat)  
!==========================================================================

! Calculates Preformed Salinity, Sstar, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! lon    : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sstar_from_sp  : Preformed Salinity                  [g/kg]
*/
double
gsw_sstar_from_sp(double sp, double p, double lon, double lat)
{
	GSW_TEOS10_CONSTANTS;
	double	saar, sstar_baltic;

    /*
	!In the Baltic Sea, Sstar = SA.
    */
	sstar_baltic	= gsw_sa_from_sp_baltic(sp,lon,lat);
	if (sstar_baltic < GSW_ERROR_LIMIT)
	    return (sstar_baltic);
	saar		= gsw_saar(p,lon,lat);
	if (saar == GSW_INVALID_VALUE)
	    return (saar);
	return (gsw_ups*sp*(1 - 0.35e0*saar));
}
