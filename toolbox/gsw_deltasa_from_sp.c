/*
!==========================================================================
function gsw_deltasa_from_sp(sp,p,lon,lat)  
!==========================================================================

! Calculates Absolute Salinity Anomaly, deltaSA, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! lon    : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_deltasa_from_sp : Absolute Salinty Anomaly           [g/kg]
*/
double
gsw_deltasa_from_sp(double sp, double p, double lon, double lat)
{
	double	res;

	res	= gsw_sa_from_sp(sp,p,lon,lat) - gsw_sr_from_sp(sp);
	if (res > GSW_ERROR_LIMIT)
	    res	= GSW_INVALID_VALUE;
	return (res);
}
