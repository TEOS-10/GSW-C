/*
!==========================================================================
function gsw_sp_from_sa_baltic(sa,lon,lat)
!==========================================================================

! For the Baltic Sea, calculates Practical Salinity with a value
! computed analytically from Absolute Salinity
!
! sa     : Absolute Salinity                               [g/kg]
! lon    : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sp_from_sa_baltic  : Practical Salinity              [unitless]
*/
double
gsw_sp_from_sa_baltic(double sa, double lon, double lat)
{
	GSW_TEOS10_CONSTANTS;
	GSW_BALTIC_DATA;
	double	xx_left, xx_right, return_value;

	if (xb_left[1] < lon  && lon < xb_right[0]  && yb_left[0] < lat  &&
	    lat < yb_left[2]) {
  
	    xx_left	= gsw_util_xinterp1(yb_left, xb_left, 3, lat);
    
	    xx_right	= gsw_util_xinterp1(yb_right, xb_right, 2, lat);
    
	    if (xx_left <= lon  && lon <= xx_right)
		return_value	= (35.0/(gsw_sso - 0.087))*(sa - 0.087);
	    else
		return_value	= GSW_INVALID_VALUE;
	} else
	    return_value	= GSW_INVALID_VALUE;

	return (return_value);
}
