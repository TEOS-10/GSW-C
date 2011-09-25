/*
**  $Id: gsw_saar.c,v 0fa6ed68e79e 2011/09/25 18:18:19 fdelahoyde $
**
**  TEOS-10 V3.0
*/
#include <gswteos-10.h>
#include <gsw_saar_data.c>

static double
gsw_sum(double *x, int n)
{
	int	i;
	double	val;

	for (val=0.0, i=0; i<n; val += x[i], i++);
	return (val);
}

#define max(a,b)	(((a)>(b))?(a):(b))
#define sum(x)		gsw_sum(x, sizeof (x)/sizeof (double))

/*
!==========================================================================
function gsw_saar(p,long,lat)
!==========================================================================

! Calculates the Absolute Salinity Anomaly Ratio, SAAR.
!
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_saar : Absolute Salinity Anomaly Ratio               [unitless]
*/
double
gsw_saar(double p, double lon, double lat)
{
	int	nx=gsw_nx, ny=gsw_ny, nz=gsw_nz;
	int	indx0, indy0, indz0, i, j, k;
	int	nmean, flag_saar, ndepth_index;
	int	deli[4] = {0,1,1,0}, delj[4] = {0,0,1,1};
	double	saar[4], saar_old[4];

	double	dlong, dlat;
	double	lon0_in, sa_upper, sa_lower;
	double	r1, s1, t1, saar_mean, ndepth_max, return_value;


	return_value	 = GSW_INVALID_VALUE;

	if (lat  <  -86e0  ||  lat  >  90e0)
	    return (return_value);

	if (lon  <  0.0)
	    lon	= lon + 360.0;

	dlong	= longs_ref[1]-longs_ref[0];
	dlat	= lats_ref[1]-lats_ref[0];

	indx0	= floor(0 + (nx-1)*(lon-longs_ref[0])/(longs_ref[nx-1]-
		    longs_ref[0]));
	if (indx0 == nx-1)
	    indx0	= nx-2;

	indy0 = floor(0 + (ny-1)*(lat-lats_ref[0])/(lats_ref[ny-1]-
		    lats_ref[0]));
	if(indy0 == ny-1)
	    indy0	= ny-2;

	ndepth_max	= -1;
	for (k=0; k < 4; k++) {
	    ndepth_index	= indy0+delj[k]+(indx0+deli[k])*ny;
	    if (ndepth_ref[ndepth_index] > 0.0)
		ndepth_max = max(ndepth_max, ndepth_ref[ndepth_index]);
	}

	if (ndepth_max == -1.e0)
	    return (0.0);

	if (p > p_ref[(int)(ndepth_max)-1])
	    p	= p_ref[(int)(ndepth_max)-1];
	indz0	= gsw_indx(p_ref,nz,p);
    
	r1	= (lon-longs_ref[indx0])/(longs_ref[indx0+1]-longs_ref[indx0]);
	s1	= (lat-lats_ref[indy0])/(lats_ref[indy0+1]-lats_ref[indy0]);
	t1	= (p-p_ref[indz0])/(p_ref[indz0+1]-p_ref[indz0]);

	for (k=0; k<4; k++)
	    saar[k]	= saar_ref[indz0+nz*(indy0+delj[k]+(indx0+deli[k])*ny)];

	if (260.0 <= lon && lon <= 291.999 && 3.4 <= lat && lat <= 19.55) {
	    memmove(saar_old,saar,4*sizeof (double));
	    gsw_add_barrier(saar_old,lon,lat,longs_ref[indx0],
			lats_ref[indy0],dlong,dlat,saar);
		/* FIXME v FIXME */
	} else if (fabs(sum(saar))  >=  1e10) {
	    memmove(saar_old,saar,4*sizeof (double));
	    gsw_add_mean(saar_old,lon,lat,saar);
	}

	sa_upper	= (1.0-s1)*(saar[0] + r1*(saar[1]-saar[0])) +
			s1*(saar[3] + r1*(saar[2]-saar[3]));

	for (k=0; k<4; k++)
	    saar[k]	= saar_ref[indz0+1+nz*(indy0+delj[k]+
				(indx0+deli[k])*ny)];

	if (260.0 <= lon && lon <= 291.999 && 3.4 <= lat && lat <= 19.55) {
	    memmove(saar_old,saar,4*sizeof (double));
	    gsw_add_barrier(saar_old,lon,lat,longs_ref[indx0],
				lats_ref[indy0],dlong,dlat,saar);
		/* FIXME v FIXME */
	} else if (fabs(sum(saar))  >=  1e10) {
	    memmove(saar_old,saar,4*sizeof (double));
	    gsw_add_mean(saar_old,lon,lat,saar);
	}

	sa_lower	= (1.0-s1)*(saar[0] + r1*(saar[1]-saar[0])) +
				s1*(saar[3] + r1*(saar[2]-saar[3]));
	if (fabs(sa_lower)  >=  1e10)
	    sa_lower	= sa_upper;
	return_value	= sa_upper + t1*(sa_lower-sa_upper);

	if (fabs(return_value) >= 1e10)
	    return_value	= GSW_INVALID_VALUE;

	return (return_value);
}

/*
!==========================================================================
function gsw_delta_sa_ref(p,lon,lat)
!==========================================================================

! Calculates the Absolute Salinity Anomaly reference value, delta_SA_ref.
!
! p      : sea pressure                                    [dbar]
! lon   : longiture                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_delta_sa_ref : Absolute Salinity Anomaly reference value    [g/kg]
*/
double
gsw_delta_sa_ref(double p, double lon, double lat)
{
	int	nx=gsw_nx, ny=gsw_ny, nz=gsw_nz;

	int	indx0, indy0, indz0, i, j, k, ndepth_index;
	int	nmean, flag_dsar;
	int	deli[4]={0,1,1,0}, delj[4]={0,0,1,1};

	double	dsar[4], dsar_old[4];
	double	dlong, dlat;
	double	return_value, lon0_in, sa_upper, sa_lower;
	double	r1, s1, t1, dsar_mean, ndepth_max;

	return_value	= GSW_INVALID_VALUE;

	if (lat < -86.0  ||  lat  >  90.0)
	    return (return_value);

	if (lon < 0.0)
	    lon	= lon + 360.0;

	dlong	= longs_ref[1]-longs_ref[0];
	dlat	= lats_ref[1]-lats_ref[0];

	indx0	= floor(0 + (nx-1)*(lon-longs_ref[0])/
			(longs_ref[nx-1]-longs_ref[0]));
	if (indx0 == nx-1)
	    indx0	= nx-2;

	indy0	= floor(0 + (ny-1)*(lat-lats_ref[0])/
			(lats_ref[ny-1]-lats_ref[0]));
	if (indy0 == ny-1)
	    indy0	= ny-2;

	ndepth_max	= -1;
	for (k=0; k<4; k++) {
	    ndepth_index	= indy0+delj[k]+(indx0+deli[k])*ny;
	    if (ndepth_ref[ndepth_index] > 0.0)
		ndepth_max	= max(ndepth_max, ndepth_ref[ndepth_index]);
	}

	if (ndepth_max == -1.0)
	    return (0.0);

	if (p > p_ref[(int)(ndepth_max)-1])
	    p	= p_ref[(int)(ndepth_max)-1];
	indz0	= gsw_indx(p_ref,nz,p);
    
	r1	= (lon-longs_ref[indx0])/
			(longs_ref[indx0+1]-longs_ref[indx0]);
	s1	= (lat-lats_ref[indy0])/
			(lats_ref[indy0+1]-lats_ref[indy0]);
	t1	= (p-p_ref[indz0])/
			(p_ref[indz0+1]-p_ref[indz0]);

	for (k=0; k < 4; k++)
	    dsar[k]	= delta_sa_ref[indz0+nz*(indy0+delj[k]+
				(indx0+deli[k])*ny)];

	if (260.0 <= lon && lon <= 291.999 && 3.4 <= lat && lat <= 19.55) {
	    memmove(dsar_old,dsar,4*sizeof (double));
	    gsw_add_barrier(dsar_old,lon,lat,longs_ref[indx0],
				lats_ref[indy0],dlong,dlat,dsar);
	} else if (fabs(sum(dsar)) >= 1e10) {
	    memmove(dsar_old,dsar,4*sizeof (double));
	    gsw_add_mean(dsar_old,lon,lat,dsar);
	}

	sa_upper	= (1.0-s1)*(dsar[0] + r1*(dsar[1]-dsar[0])) +
				s1*(dsar[3] + r1*(dsar[2]-dsar[3]));

	for (k=0; k<4; k++)
	    dsar[k]	= delta_sa_ref[indz0+1+nz*(indy0+delj[k]+
				(indx0+deli[k])*ny)];

	if (260.0 <= lon && lon <= 291.999 && 3.4 <= lat && lat <= 19.55) {
	    memmove(dsar_old,dsar,4*sizeof (double));
	    gsw_add_barrier(dsar_old,lon,lat,longs_ref[indx0],
				lats_ref[indy0],dlong,dlat,dsar);
	} else if (fabs(sum(dsar)) >= 1e10) {
	    memmove(dsar_old,dsar,4*sizeof (double));
	    gsw_add_mean(dsar_old,lon,lat,dsar);
	}

	sa_lower	= (1.0-s1)*(dsar[0] + r1*(dsar[1]-dsar[0])) +
				s1*(dsar[3] + r1*(dsar[2]-dsar[3]));
	if (fabs(sa_lower) >= 1e10)
	    sa_lower	= sa_upper;
	return_value	= sa_upper + t1*(sa_lower-sa_upper);

	if (fabs(return_value) >= 1e10)
	    return (GSW_INVALID_VALUE);

	return (return_value);
}
/*
**  The End
*/
