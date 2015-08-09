/*
!==========================================================================
function gsw_sa_from_sp(sp,p,lon,lat)       
!==========================================================================

! Calculates Absolute Salinity, SA, from Practical Salinity, SP
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! lon	 : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]
!
! gsw_sa_from_sp   : Absolute Salinity                     [g/kg]
*/
double
gsw_sa_from_sp(double sp, double p, double lon, double lat)
{
	GSW_TEOS10_CONSTANTS;
	double	saar, gsw_sa_baltic;

	gsw_sa_baltic	= gsw_sa_from_sp_baltic(sp,lon,lat);
	if (gsw_sa_baltic < GSW_ERROR_LIMIT)
	    return (gsw_sa_baltic);
	saar		= gsw_saar(p,lon,lat);
	if (saar == GSW_INVALID_VALUE)
	    return (saar);
	return (gsw_ups*sp*(1.e0 + saar));
}
