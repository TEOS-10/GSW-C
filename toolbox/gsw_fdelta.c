/*
!==========================================================================
function gsw_fdelta(p,lon,lat)
!==========================================================================

! Calculates fdelta. 
!
! p      : sea pressure                                    [dbar]
! lon    : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_fdelta : Absolute Salinty Anomaly                    [unitless]
*/
double
gsw_fdelta(double p, double lon, double lat)
{
	double	sa, saar;

	saar	= gsw_saar(p,lon,lat);
	if (saar >= GSW_ERROR_LIMIT)
	    sa	= GSW_INVALID_VALUE;
	else
	    sa	= ((1.0 + 0.35)*saar)/(1.0 - 0.35*saar);
	return (sa);
}
