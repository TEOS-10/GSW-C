/*
!==========================================================================
function gsw_sp_from_sa(sa,p,lon,lat)  
!==========================================================================

! Calculates Practical salinity, sp, from Absolute salinity, sa  
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! lon    : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]
!
! gsw_sp_from_sa      : Practical Salinity                 [unitless]
*/
double
gsw_sp_from_sa(double sa, double p, double lon, double lat)
{
	GSW_TEOS10_CONSTANTS;
	double	saar, gsw_sp_baltic;

	gsw_sp_baltic	= gsw_sp_from_sa_baltic(sa,lon,lat);
	if (gsw_sp_baltic < GSW_ERROR_LIMIT)
	    return (gsw_sp_baltic);
	saar	= gsw_saar(p,lon,lat);
	if (saar == GSW_INVALID_VALUE)
	    return (saar);
	return ((sa/gsw_ups)/(1e0 + saar));
}
