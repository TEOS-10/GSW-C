/*
!==========================================================================
function gsw_sa_from_sp_baltic(sp,lon,lat)
!==========================================================================

! For the Baltic Sea, calculates Absolute Salinity with a value
! computed analytically from Practical Salinity
!
! sp     : Practical Salinity                              [unitless]
! lon    : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! sa_from_sp_baltic : Absolute Salinity                    [g/kg]
*/
double
gsw_sa_from_sp_baltic(double sp, double lon, double lat)
{
	GSW_TEOS10_CONSTANTS;
	GSW_BALTIC_DATA;
	double	xx_left, xx_right, return_value;

	if (xb_left[1] < lon  && lon < xb_right[0]  && yb_left[0] < lat  &&
	    lat < yb_left[2]) {
  
	    xx_left	= gsw_util_xinterp1(yb_left, xb_left, 3, lat);
    
	    xx_right	= gsw_util_xinterp1(yb_right, xb_right, 2, lat);
    
	    if (xx_left <= lon  && lon <= xx_right)
		return_value	=((gsw_sso - 0.087)/35.0)*sp + 0.087;
	    else
		return_value	= GSW_INVALID_VALUE;
	} else
	    return_value	= GSW_INVALID_VALUE;

	return (return_value);
}
