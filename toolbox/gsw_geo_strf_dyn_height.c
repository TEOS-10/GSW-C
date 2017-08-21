/*
!==========================================================================
pure function gsw_geo_strf_dyn_height (sa, ct, p, p_ref)
!==========================================================================
!
!  Calculates dynamic height anomaly as the integral of specific volume
!  anomaly from the pressure p of the bottle to the reference pressure
!  p_ref.
!
!  Hence, geo_strf_dyn_height is the dynamic height anomaly with respect
!  to a given reference pressure.  This is the geostrophic streamfunction
!  for the difference between the horizontal velocity at the pressure
!  concerned, p, and the horizontal velocity at p_ref.  Dynamic height
!  anomaly is the geostrophic streamfunction in an isobaric surface.  The
!  reference values used for the specific volume anomaly are
!  SSO = 35.16504 g/kg and CT = 0 deg C.  This function calculates
!  specific volume anomaly using the computationally efficient
!  expression for specific volume of Roquet et al. (2015).
!
!  This function evaluates the pressure integral of specific volume using
!  SA and CT interpolated with respect to pressure using the method of
!  Reiniger and Ross (1968).  It uses a weighted mean of (i) values
!  obtained from linear interpolation of the two nearest data points, and
!  (ii) a linear extrapolation of the pairs of data above and below.  This
!  "curve fitting" method resembles the use of cubic splines.
!
!  SA    =  Absolute Salinity                                      [ g/kg ]
!  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
!  p     =  sea pressure                                           [ dbar ]
!           ( i.e. absolute pressure - 10.1325 dbar )
!  p_ref =  reference pressure                                     [ dbar ]
!           ( i.e. reference absolute pressure - 10.1325 dbar )
!
!  geo_strf_dyn_height  =  dynamic height anomaly               [ m^2/s^2 ]
!   Note. If p_ref exceeds the pressure of the deepest bottle on a
!     vertical profile, the dynamic height anomaly for each bottle
!     on the whole vertical profile is returned as NaN.
!--------------------------------------------------------------------------
*/
static void p_sequence(double p1,double p2,double max_dp_i,double *pseq,
	int *nps);	/* forward reference */

double	*	/* Returns NULL on error, dyn_height if okay */
gsw_geo_strf_dyn_height(double *sa, double *ct, double *p, double p_ref,
	int n_levels, double *dyn_height)
{
	GSW_TEOS10_CONSTANTS;
	int	m_levels = (n_levels <= 0) ? 1 : n_levels,
		p_cnt, top_pad, i, nz, ibottle, ipref, np_max, np, ibpr=0,
		*iidata;
	double	dp_min, dp_max, p_min, p_max, max_dp_i,
		*b, *b_av, *dp, *dp_i, *sa_i=NULL, *ct_i, *p_i=NULL,
		*geo_strf_dyn_height0;

/*
!--------------------------------------------------------------------------
!  This max_dp_i is the limit we choose for the evaluation of specific
!  volume in the pressure integration.  That is, the vertical integration
!  of specific volume with respect to pressure is perfomed with the pressure
!  increment being no more than max_dp_i (the default value being 1 dbar).
!--------------------------------------------------------------------------
*/
	max_dp_i = 1.0;

	if ((nz = m_levels) <= 1)
	    return (NULL);

	dp = malloc(nz*sizeof (double));
	dp_min = 11000.0;
	dp_max = -11000.0;
	for (i=0; i<nz-1; i++) {
	    if ((dp[i] = p[i+1] - p[i]) < dp_min)
		dp_min = dp[i];
	    if (dp[i] > dp_max)
		dp_max = dp[i];
	}

	if (dp_min <= 0.0) {
	    /* pressure must be monotonic */
	    free(dp);
	    return (NULL);
	}
	p_min = p[0];
	p_max = p[nz-1];

	if (p_ref > p_max) {
	    /*the reference pressure p_ref is deeper than all bottles*/
	    free(dp);
	    return (NULL);
	}

	/* Determine if there is a "bottle" at exactly p_ref */
	ipref = -1;
	for (ibottle = 0; ibottle < nz; ibottle++) {
	    if (p[ibottle] == p_ref) {
	        ipref = ibottle;
		break;
	    }
	}
	if ((dp_max <= max_dp_i) && (p[0] == 0.0) && (ipref >= 0)) {
	    /*
	    !vertical resolution is good (bottle gap is no larger than max_dp_i)
	    ! & the vertical profile begins at the surface (i.e. at p = 0 dbar)
	    ! & the profile contains a "bottle" at exactly p_ref.
 	    */
	    b = malloc(3*nz*sizeof (double));
	    b_av = b+nz; geo_strf_dyn_height0 = b_av+nz;
	    for (i=0; i<nz; i++) {
		b[i] = gsw_specvol_anom_standard(sa[i],ct[i],p[i]);
		if (i > 0)
		    b_av[i-1] = 0.5*(b[i] + b[i-1]);
	    }
	    /*
	    ! "geo_strf_dyn_height0" is the dynamic height anomaly with respect
	    ! to p_ref = 0 (the surface).
	    */
	    geo_strf_dyn_height0[0] = 0.0;
	    for (i=1; i<nz; i++)
		geo_strf_dyn_height0[i] = b_av[i]*dp[i]*db2pa;
	    for (i=1; i<nz; i++) /* cumulative sum */
		geo_strf_dyn_height0[i] = geo_strf_dyn_height0[i-1]
		                          - geo_strf_dyn_height0[i];
	    for (i=0; i<nz; i++)
		dyn_height[i] = geo_strf_dyn_height0[i]
				- geo_strf_dyn_height0[ipref];
	    free(b);
	} else {
	/*
	! Test if there are vertical gaps between adjacent "bottles" which are
	! greater than max_dp_i, and that there is a "bottle" exactly at the
	! reference pressure.
	*/
	    iidata = malloc((nz+1)*sizeof (int));

	    if ((dp_max <= max_dp_i) && (ipref >= 0)) {
	    /*
	    ! Vertical resolution is already good (no larger than max_dp_i), and
	    ! there is a "bottle" at exactly p_ref.
	    */
		sa_i = malloc(2*(nz+1)*sizeof (double));
		ct_i = sa_i+nz+1;
		p_i = malloc((nz+1)*sizeof (double));;

	        if (p_min > 0.0) {
		/*
	        ! resolution is fine and there is a bottle at p_ref, but
	        ! there is not a bottle at p = 0. So add an extra bottle.
		*/
		    for (i=0; i<nz; i++) {
			sa_i[i+1]	= sa[i];
			ct_i[i+1]	= ct[i];
			p_i[i+1]	= p[i];
		    }
		    sa_i[0] = sa[0];
		    ct_i[0] = ct[0];
		    p_i[0] = 0.0;
		    ibpr = ipref+1;
		    p_cnt = nz+1;
		    for (i=0; i<p_cnt; i++)
			iidata[i] = i;
	        } else {
		/*
	        ! resolution is fine, there is a bottle at p_ref, and
	        ! there is a bottle at p = 0
		*/
		    memmove(sa_i, sa, nz*sizeof (double));
		    memmove(ct_i, ct, nz*sizeof (double));
		    memmove(p_i, p, nz*sizeof (double));
		    ibpr = ipref;
		    for (i=0; i<nz; i++)
			iidata[i] = i;
		    p_cnt = nz;
	        }

	    } else {
	    /*
	    ! interpolation is needed.
	    */
		np_max = 2*rint(p[nz-1]/max_dp_i+0.5);
		p_i = malloc(np_max*sizeof (double));
		/* sa_i is allocated below, when its size is known */

	        if (p_min > 0.0) {
		/*
	        ! there is not a bottle at p = 0.
		*/
	            if (p_ref < p_min) {
		    /*
	            ! p_ref is shallower than the minimum bottle pressure.
		    */
			p_i[0] = 0.0;
			p_sequence(p_i[0],p_ref,max_dp_i, p_i+1,&np);
			ibpr = p_cnt = np;
			p_cnt++;
			p_sequence(p_ref,p_min,max_dp_i, p_i+p_cnt,&np);
	                p_cnt += np;
	                top_pad = p_cnt;
	            } else {
		    /*
	            ! p_ref is deeper than the minimum bottle pressure.
		    */
			p_i[0] = 0.0;
			p_i[1] = p_min;
	                top_pad = 2;
	                p_cnt = 2;
	            }
	        } else {
		/*
	        ! there is a bottle at p = 0.
		*/
	            p_i[0] = p_min;
	            top_pad = 1;
	            p_cnt = 1;
	        }

		for (ibottle=0; ibottle < nz-1; ibottle++) {

		    iidata[ibottle] = p_cnt-1;
	            if (p[ibottle] == p_ref) ibpr = p_cnt-1;

	            if (p[ibottle] < p_ref && p[ibottle+1] > p_ref) {
		    /*
	            ! ... reference pressure is spanned by bottle pairs -
	            ! need to include p_ref as an interpolated pressure.
		    */
	                p_sequence(p[ibottle],p_ref,max_dp_i, p_i+p_cnt,&np);
	                p_cnt += np;
			ibpr = p_cnt-1;
	                p_sequence(p_ref,p[ibottle+1],max_dp_i,p_i+p_cnt,&np);
	                p_cnt += np;
	            } else {
		    /*
	            ! ... reference pressure is not spanned by bottle pairs.
		    */
	                p_sequence(p[ibottle],p[ibottle+1],max_dp_i,
				p_i+p_cnt,&np);
	                p_cnt += np;
	            }

	        }

		iidata[nz-1] = p_cnt-1;
	        if (p[nz-1] == p_ref) ibpr = p_cnt-1;

		sa_i = malloc(2*p_cnt*sizeof (double));
		ct_i = sa_i+p_cnt;

		if (top_pad > 1) {
	            gsw_linear_interp_sa_ct(sa,ct,p,nz,
			p_i,top_pad-1,sa_i,ct_i);
		}
	        gsw_rr68_interp_sa_ct(sa,ct,p,nz,p_i+top_pad-1,p_cnt-top_pad+1,
		                      sa_i+top_pad-1,ct_i+top_pad-1);
	    }

	    b = malloc(4*p_cnt*sizeof (double));
	    b_av = b+p_cnt; dp_i = b_av+p_cnt;
	    geo_strf_dyn_height0 = dp_i+p_cnt;
	    for (i=0; i<p_cnt; i++) {
		b[i] = gsw_specvol_anom_standard(sa_i[i],ct_i[i],p_i[i]);
		if (i > 0) {
		    dp_i[i-1] = p_i[i]-p_i[i-1];
		    b_av[i-1] = 0.5*(b[i] + b[i-1]);
		}
	    }
	    /*
	    ! "geo_strf_dyn_height0" is the dynamic height anomaly with respect
	    ! to p_ref = 0 (the surface).
	    */
	    geo_strf_dyn_height0[0] = 0.0;
	    for (i=1; i<p_cnt; i++)
		geo_strf_dyn_height0[i] = b_av[i-1]*dp_i[i-1];
	    for (i=1; i<p_cnt; i++) /* cumulative sum */
		geo_strf_dyn_height0[i] = geo_strf_dyn_height0[i-1]
		                          - geo_strf_dyn_height0[i];
	    for (i=0; i<nz; i++)
		dyn_height[i] = (geo_strf_dyn_height0[iidata[i]]
				- geo_strf_dyn_height0[ibpr])*db2pa;

	    free(b);
	    free(iidata);
	    if (sa_i != NULL)
		free(sa_i);
	    if (p_i != NULL)
		free(p_i);

	}
	free(dp);
	return (dyn_height);
}

static void
p_sequence(double p1, double p2, double max_dp_i, double *pseq, int *nps)
{
	double	dp, pstep;
	int		n, i;

	dp = p2 - p1;
	n = ceil(dp/max_dp_i);
	pstep = dp/n;

	if (nps != NULL) *nps = n;

	/*
	! Generate the sequence ensuring that the value of p2 is exact to
	! avoid round-off issues, ie. don't do "pseq = p1+pstep*(i+1)".
	*/
	for (i=0; i<n; i++)
	    pseq[i] = p2-pstep*(n-1-i);

}
