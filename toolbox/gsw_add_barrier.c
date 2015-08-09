/*
!==========================================================================
subroutine gsw_add_barrier(input_data,lon,lat,long_grid,lat_grid,dlong_grid,dlat_grid,output_data)
!==========================================================================

!  Adds a barrier through Central America (Panama) and then averages
!  over the appropriate side of the barrier
! 
!  data_in      :  data                                         [unitless]
!  lon          :  Longitudes of data degrees east              [0 ... +360]
!  lat          :  Latitudes of data degrees north              [-90 ... +90]
!  longs_grid   :  Longitudes of regular grid degrees east      [0 ... +360]
!  lats_grid    :  Latitudes of regular grid degrees north      [-90 ... +90]
!  dlongs_grid  :  Longitude difference of regular grid degrees [deg longitude]
!  dlats_grid   :  Latitude difference of regular grid degrees  [deg latitude]
!
!  output_data  : average of data depending on which side of the 
!                 Panama canal it is on                         [unitless]
*/
void
gsw_add_barrier(double *input_data, double lon, double lat,
		double long_grid, double lat_grid, double dlong_grid,
		double dlat_grid, double *output_data)
{
	GSW_SAAR_DATA;
	int	above_line[4];
	int	k, nmean, above_line0, kk;
	double	r, lats_line, data_mean;

	k		= gsw_util_indx(longs_pan,npan,lon);
			/*   the lon/lat point */
	r		= (lon-longs_pan[k])/(longs_pan[k+1]-longs_pan[k]);
	lats_line	= lats_pan[k] + r*(lats_pan[k+1]-lats_pan[k]);

	above_line0	= (lats_line <= lat);

	k		= gsw_util_indx(longs_pan,npan,long_grid);
			/*the 1 & 4 lon/lat points*/ 
	r		= (long_grid-longs_pan[k])/
				(longs_pan[k+1]-longs_pan[k]);
	lats_line	= lats_pan[k] + r*(lats_pan[k+1]-lats_pan[k]);

	above_line[0]	= (lats_line <= lat_grid);
	above_line[3]	= (lats_line <= lat_grid+dlat_grid);

	k		= gsw_util_indx(longs_pan,6,long_grid+dlong_grid);
			/*the 2 & 3 lon/lat points */
	r		= (long_grid+dlong_grid-longs_pan[k])/
			(longs_pan[k+1]-longs_pan[k]);
	lats_line	= lats_pan[k] + r*(lats_pan[k+1]-lats_pan[k]);

	above_line[1]	= (lats_line <= lat_grid);
	above_line[2]	= (lats_line <= lat_grid+dlat_grid);

	nmean		= 0;
	data_mean	= 0.0;

	for (kk=0; kk<4; kk++) {
	    if ((fabs(input_data[kk]) <= 100.0) &&
		above_line0 == above_line[kk]) {
		nmean	= nmean+1;
		data_mean	= data_mean+input_data[kk];
	    }
	}
	if (nmean == 0)
	    data_mean	= 0.0;	/*errorreturn*/
	else
	    data_mean	= data_mean/nmean;

	for (kk=0; kk<4; kk++) {
	    if ((fabs(input_data[kk]) >= 1.0e10) ||
		above_line0 != above_line[kk])
		output_data[kk]	= data_mean;
	    else
		output_data[kk]	= input_data[kk];
	}

	return;
}
