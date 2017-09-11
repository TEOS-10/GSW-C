/*
This is a replacement for gsw_geo_strf_dyn_height, with a different
signature and interpolation algorithms.
!==========================================================================
int   (returns nonzero on error, 0 if OK)
gsw_geo_strf_dyn_height_1(double *sa, double *ct, double *p, double p_ref,
    int nz, double *dyn_height, double max_dp_i, int interp_method)
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
!  SA and CT interpolated with respect to pressure. The interpolation method
!  may be chosen as linear or "PCHIP", piecewise cubic Hermite using a shape-
!  preserving algorithm for setting the derivatives.
!
!  SA    =  Absolute Salinity                                      [ g/kg ]
!  CT    =  Conservative Temperature (ITS-90)                     [ deg C ]
!  p     =  sea pressure  (increasing with index)                  [ dbar ]
!           ( i.e. absolute pressure - 10.1325 dbar )
!  nz    =  number of points in each array
!  p_ref =  reference pressure                                     [ dbar ]
!           ( i.e. reference absolute pressure - 10.1325 dbar )
!  geo_strf_dyn_height  =  dynamic height anomaly               [ m^2/s^2 ]
!  max_dp_i = maximum pressure difference between points for triggering
!              interpolation.
!  interp_method = 1 for linear, 2 for PCHIP
!
!   Note. If p_ref falls outside the range of a
!     vertical profile, the dynamic height anomaly for each bottle
!     on the whole vertical profile is returned as NaN.
!--------------------------------------------------------------------------
*/

/*
    Make a new grid based on an original monotonic array, typically P, such that
    on the new monotonic grid, p_i:

        1) The first element is p[0], the last is p[nz-1].
        2) Approximate integer multiples of dp are included.
        3) All original p points are included.
        4) The value p_ref is included.

    This function allocates memory; it is
    the responsibility of subsequent code to free this memory.

    In addition to the new p_i grid, the function returns the array
    of indices (p_indices)
    of the original p grid points in p_i, and a pointer to
    the index of p_ref within p_i.

    Its return argument is the number of elements in p_i.
*/
static int refine_grid_for_dh(double *p, double p_ref, int nz,
    double dp,
    double *p_i, int ni_max,  /* size of p_i array; larger than needed */
    int *p_indices, int *p_ref_ind_ptr)
{
    int i, iuniform, iorig;
    double p_next;
    /* Don't add a new point if it is within p_tol of an original. */
    double p_tol = 0.001 * dp;

    p_i[0] = p[0];
    p_indices[0] = 0;
    *p_ref_ind_ptr = -1;  /* initialize to a flag value */
    if (p_ref <= p[0] + p_tol)
    {
        *p_ref_ind_ptr = 0;
    }
    for (i=1, iuniform=1, iorig=1; i<ni_max && iorig<nz; i++)
    {
        /* Candidate insertion based on uniform grid: */
        p_next = p[0] + dp * iuniform;

        /* See if we need to insert p_ref: */
        if (*p_ref_ind_ptr == -1 && p_ref <= p_next && p_ref <= p[iorig])
        {
            p_i[i] = p_ref;
            *p_ref_ind_ptr = i;
            if (p_ref == p[iorig])
            {
                p_indices[iorig] = i;
                iorig++;
            }
            if (p_ref > p_next - p_tol)
            {
                iuniform++;
            }
            continue;
        }

        /* We did not insert p_ref, so insert either p_next or p[iorig]. */
        if (p_next < p[iorig] - p_tol)
        {
            p_i[i] = p_next;
            iuniform++;
        }
        else
        {
            p_i[i] = p[iorig];
            p_indices[iorig] = i;
            /* Skip this p_next if it is close to the point we just added. */
            if (p_next < p[iorig] + p_tol)
            {
                iuniform++;
            }
            iorig++;
        }
    }

    if (i == ni_max)
    {
        return (-1);  /* error! */
    }
    return (i);  /* number of elements in p_i */
}

/*  Linearly interpolate to the grid made by define_grid_for_dh.
    We take advantage of what we know about the grids: they match
    at the end points, and both are monotonic.
*/

static int linear_interp_SA_CT_for_dh(double *sa, double *ct, double *p, int nz,
    double *p_i, int n_i,
    double *sa_i, double *ct_i)
{
    int i, ii;
    double pfac;

    sa_i[0] = sa[0];
    sa_i[n_i-1] = sa[nz-1];
    ct_i[0] = ct[0];
    ct_i[n_i-1] = ct[nz-1];
    i = 1;
    for (ii=1; ii<n_i-1; ii++)
    {
        /* Find the second point of the pair in the original grid that
           bracket the target.
        */
        while (p[i] < p_i[ii])
        {
            i++;
            if (i == nz)
            {
                return -1;  /* error! */
            }
        }
        pfac = (p_i[ii] - p[i-1]) / (p[i] - p[i-1]);
        sa_i[ii] = sa[i-1] + pfac * (sa[i] - sa[i-1]);
        ct_i[ii] = ct[i-1] + pfac * (ct[i] - ct[i-1]);
    }
    return 0;
}


int  /* returns nonzero on error, 0 if OK */
gsw_geo_strf_dyn_height_1(double *sa, double *ct, double *p, double p_ref,
    int nz, double *dyn_height, double max_dp_i, int interp_method)
{
    GSW_TEOS10_CONSTANTS;
    int i, ipref,
        *p_indices, n_i, ni_max, err;
    double    dp_min, dp_max, p_min, p_max,
        *b, *b_av, *dp, *sa_i, *ct_i, *p_i,
        *dh_i;
    double dh_ref;

    if (nz <= 2)
        return (1);

    dp = malloc((nz-1) * sizeof(double));
    dp_min = 11000.0;
    dp_max = -11000.0;
    for (i=0; i<nz-1; i++)
    {
        dp[i] = p[i+1] - p[i];
        if (dp[i] < dp_min)
        {
             dp_min = dp[i];
         }
        if (dp[i] > dp_max)
        {
            dp_max = dp[i];
        }
    }

    if (dp_min <= 0.0) {
        /* pressure must be monotonic */
        free(dp);
        return (2);
    }
    p_min = p[0];
    p_max = p[nz-1];

    if (p_ref > p_max || p_ref < p_min) {
        /* Reference pressure must be within the data range. */
        free(dp);
        return (3);
    }

    /* Determine if there is a sample at exactly p_ref */
    ipref = -1;
    for (i = 0; i < nz; i++)
    {
        if (p[i] == p_ref)
        {
            ipref = i;
            break;
        }
    }

    if ((dp_max <= max_dp_i) && (ipref >= 0)) {
        /*
        !vertical resolution is good (bottle gap is no larger than max_dp_i)
        ! & the profile contains a "bottle" at exactly p_ref.
         */
        b = malloc(nz*sizeof (double));
        b_av = malloc((nz-1) * sizeof(double));
        for (i=0; i<nz; i++) {
            b[i] = gsw_specvol_anom_standard(sa[i],ct[i],p[i]);
        }
        for (i=0; i<(nz-1); i++) {
            b_av[i] = 0.5*(b[i+1] + b[i]);
        }
        /* First calculate dynamic height relative to the first (shallowest)
           depth. */
        dyn_height[0] = 0.0;
        for (i=1; i<nz; i++)
        {
            dyn_height[i] = dyn_height[i-1] - b_av[i-1]*dp[i-1]*db2pa;
        }
        /* Then subtract out the value at the reference pressure. */
        dh_ref = dyn_height[ipref];
        for (i=0; i<nz; i++)
        {
            dyn_height[i] -= dh_ref;
        }
        free(b);
        free(b_av);
        free(dp);
        return (0);
    }

    /*
    If we got this far, then we need to interpolate: either or both of
    inserting a point for p_ref and subdividing the intervals to keep the max
    interval less than max_dp_i.
    */

    free(dp);  /* Need to recalculate, so free here and malloc when needed. */

    ni_max = nz + (int) ceil((p[nz-1] - p[0]) / max_dp_i) + 2;
    /* Maximum possible size of new grid: Original grid size plus
       the number of dp intervals plus 1 for the p_ref,
       plus 1 so that we can know we exited the loop before we hit it.
    */

    p_i = (double *) malloc(ni_max * sizeof(double));
    p_indices = (int *) malloc(nz * sizeof(int));

    n_i = refine_grid_for_dh(p, p_ref, nz, max_dp_i,
                             p_i, ni_max,
                             p_indices, &ipref);
    /* Reminder: if successful, this allocated p_i and p_indices. */
    if (n_i == -1)
    {
        free(p_i);
        free(p_indices);
        return (4);
    }

    ct_i = malloc(n_i * sizeof(double));
    sa_i = malloc(n_i * sizeof(double));

    if (interp_method == INTERP_METHOD_LINEAR)
    {
        err = linear_interp_SA_CT_for_dh(sa, ct, p, nz,
                                         p_i, n_i,
                                         sa_i, ct_i);
        if (err) err = 5;
    }
    else if (interp_method == INTERP_METHOD_PCHIP)
    {
        err = gsw_util_pchip_interp(p, sa, nz, p_i, sa_i, n_i);
        err = err || gsw_util_pchip_interp(p, ct, nz, p_i, ct_i, n_i);
        if (err) err = 6;
    }
    else
    {
        err = 7;
    }
    if (err)
    {
        free(p_i);
        free(p_indices);
        free(ct_i);
        free(sa_i);
        return (err);
    }

    dh_i = malloc(n_i * sizeof(double));
    dp = malloc((n_i-1) * sizeof(double));
    b = malloc(n_i*sizeof (double));
    b_av = malloc((n_i-1) * sizeof(double));

    for (i=0; i<n_i; i++)
    {
        b[i] = gsw_specvol_anom_standard(sa_i[i], ct_i[i], p_i[i]);
    }
    free(ct_i);
    free(sa_i);

    for (i=0; i<(n_i-1); i++)
    {
        b_av[i] = 0.5*(b[i+1] + b[i]);
        dp[i] = p_i[i+1] - p_i[i];
    }
    free(p_i);
    /* First calculate dynamic height relative to the first (shallowest)
       depth. */
    dh_i[0] = 0.0;
    for (i=1; i<n_i; i++)
    {
        dh_i[i] = dh_i[i-1] - b_av[i-1]*dp[i-1]*db2pa;
    }
    free(b);
    free(b_av);
    free(dp);

    dh_ref = dh_i[ipref];

    for (i=0; i<nz; i++)
    {
        dyn_height[i] = dh_i[p_indices[i]] - dh_ref;
    }
    free(p_indices);
    free(dh_i);

    return 0;
}
