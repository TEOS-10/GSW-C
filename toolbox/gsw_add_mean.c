/*
!==========================================================================
subroutine gsw_add_mean(data_in,data_out)
!==========================================================================

! Replaces NaN's with non-nan mean of the 4 adjacent neighbours
!
! data_in   : data set of the 4 adjacent neighbours   
!
! data_out : non-nan mean of the 4 adjacent neighbours     [unitless]
*/
void
gsw_add_mean(double *data_in, double *data_out)
{
	int	k, nmean;
	double	data_mean;

	nmean		= 0;
	data_mean	= 0.0;

	for (k=0; k<4; k++) {
	    if (fabs(data_in[k]) <= 100.0) {
		nmean++;
		data_mean	= data_mean+data_in[k];
	    }
	}

	if (nmean == 0.0)
	    data_mean	= 0.0;    /*errorreturn*/
	else
	    data_mean	= data_mean/nmean;

	for (k=0; k<4; k++) {
	    if (fabs(data_in[k]) >= 100.0)
		data_out[k]	= data_mean;
	    else
		data_out[k]	= data_in[k];
	}
	return;
}
