#include "gswteos-10.h"
#include "gsw_internal_const.h"

/*
!==========================================================================
subroutine gsw_tracer_CT_interp (sa,ct,p,p_i)
!==========================================================================
!
! SA and CT interpolation to p_i on a cast
!
!  tracer   =  tracer                                                 [ ? ]
!  CT   =  Conservative Temperature (ITS-90)                      [ deg C ]
!  p    =  sea pressure                                            [ dbar ]
!           ( i.e. absolute pressure - 10.1325 dbar )
!  p_i  =  specific query points at which the interpolated tracer_i and CT_i
!            are required                                          [ dbar ]

!
!  tracer_i  =  interpolated SA values at pressures p_i                [ g/kg ]
!  CT_i  =  interpolated CT values at pressures p_i               [ deg C ]
!
!--------------------------------------------------------------------------
*/
void
gsw_tracer_ct_interp(double *sa, double *ct, double *p, int m,
        double *p_i, int m_i, double factor, double *tracer_i, double *ct_i)
{
        double  rec_factor  = 1./factor,

                sin_kpi_on_16[] = {
                    1.950903220161283e-1,
                    3.826834323650898e-1,
                    5.555702330196022e-1,
                    7.071067811865475e-1,
                    8.314696123025452e-1,
                    9.238795325112867e-1,
                    9.807852804032304e-1
                },
                cos_kpi_on_16[] = {
                    9.807852804032304e-1,
                    9.238795325112867e-1,
                    8.314696123025452e-1,
                    7.071067811865476e-1,
                    5.555702330196023e-1,
                    3.826834323650898e-1,
                    1.950903220161283e-1
                };

        int i, j, k, prof_len, 
            not_monotonic, unique_count, new_len, p_all_len,
            i_min_p_obs, i_obs_plus_interp_len, i_surf_and_obs_plus_interp_len,
            i_out_len, i_2_len, i_frozen, i_shallower, i_above, i_above_i, i_below_i;
        int *p_idx, *p_all_idx, *i_obs_plus_interp, *i_surf_and_obs_plus_interp,
            *i_out, *i_1, *i_2, *i_3;
        double d, ct_f, unique_p, ct_sum, tracer_sum, min_p_obs, max_p_obs;
        double *tracer_obs, *ct_obs, *p_obs, *p_i_tmp,
               *p_sort, *tracer_sort, *ct_sort, *p_all, *p_all_sort,
               *p_obs_plus_interp, *p_surf_and_obs_plus_interp,
               *independent_variable, *independent_variable_obs_plus_interp,
               *scaled_tracer_obs, *v_tmp, *q_tmp, *v_i, *q_i,
               *tracer_i_obs_plus_interp, *ct_i_obs_plus_interp,
               *tracer_i_tooutput, *ct_i_tooutput;

        tracer_obs = (double *) malloc(m*sizeof (double));
        ct_obs = (double *) malloc(m*sizeof (double));
        p_obs = (double *) malloc(m*sizeof (double));
        p_idx = (int *) malloc(m*sizeof (int));
        p_sort = (double *) malloc(m*sizeof (double));
        tracer_sort = (double *) malloc(m*sizeof (double));
        ct_sort = (double *) malloc(m*sizeof (double));

        if (m < 4)
            return; // There must be at least 4 bottles'

        // Check if interpolating pressure is monotonic
        for (i=0; i<m_i-1; ++i) {
            d = p_i[i + 1] - p_i[i];
            if (d < 0.) {
                return;
            }
        }

        for (i=0; i<m_i; i++) {
            tracer_i[i] = NAN;
            ct_i[i] = NAN;
        }
        
        // Find NaNs in profile
        prof_len = 0;
        for (i=0; i<m; ++i) {
            d = sa[i] + ct[i] + p[i];
            if (!isnan(d)) {
                tracer_obs[prof_len] = sa[i];
                ct_obs[prof_len] = ct[i];

                p_obs[prof_len] = 1e-3 * round(1e3 * p[i]);

                ++prof_len;
            }
        }

        if (prof_len < 2) {
            return;
        }

        p_i_tmp = (double *) malloc(m_i*sizeof (double));
        for (i=0; i<m_i; ++i) {
            p_i_tmp[i] = 1e-3 * round(1e3 * p_i[i]);
        }

        // Check if profile pressure values are monotonic
        // If they are not, then sort pressure values along with SA and CT
        not_monotonic = 0;
        for (i=0; i<prof_len-1; ++i) {
            d = p_obs[i+1] - p_obs[i];
            if (d < 0.) {
                ++not_monotonic;
            }
        }

        if (not_monotonic > 0) {
            gsw_util_sort_real(p_obs, prof_len, p_idx);
            for (i=0; i<prof_len; ++i) {
                p_sort[i] = p_obs[p_idx[i]];
                tracer_sort[i] = tracer_obs[p_idx[i]];
                ct_sort[i] = ct_obs[p_idx[i]];
            }

            // Once sorted, only save unique pressure values.
            // SA and CT observations with the same pressure
            // will be averaged.
            new_len = 0;
            unique_count = 1;
            unique_p = p_sort[0];
            tracer_sum = tracer_sort[0];
            ct_sum = ct_sort[0];
            for (i=1; i<prof_len; ++i) {
                if (unique_p == p_sort[i]) {
                    tracer_sum += tracer_sort[i];
                    ct_sum += ct_sort[i];
                    ++unique_count;
                } else {
                    p_obs[new_len] = unique_p;
                    tracer_obs[new_len] = tracer_sum/(double)unique_count;
                    ct_obs[new_len] = ct_sum/(double)unique_count;
                    ++new_len;
                    unique_p = p_sort[i];
                    tracer_sum = tracer_sort[i];
                    ct_sum = ct_sort[i];
                    unique_count = 1;
                }
            }
            p_obs[new_len] = unique_p;
            tracer_obs[new_len] = tracer_sum/(double)unique_count;
            ct_obs[new_len] = ct_sum/(double)unique_count;
            ++new_len;
            prof_len = new_len;
        }

        // Combine pressure values of observed and interpolated SA and CT
        p_all_len = prof_len + m_i;

        p_all = (double *) malloc(p_all_len*sizeof (double));
        p_all_sort = (double *) malloc(p_all_len*sizeof (double));
        p_all_idx = (int *) malloc(p_all_len*sizeof (int));

        memcpy(p_all, p_obs, prof_len*sizeof (double));
        memcpy(&p_all[prof_len], p_i_tmp, m_i*sizeof (double));
        gsw_util_sort_real(p_all, p_all_len, p_all_idx);
        for (i=0; i<p_all_len; ++i) {
            p_all_sort[i] = p_all[p_all_idx[i]];
        }
        new_len = 0;
        unique_p = p_all_sort[0];
        for (i=1; i<p_all_len; ++i) {
            if (unique_p != p_all_sort[i]) {
                p_all[new_len] = unique_p;
                ++new_len;
                unique_p = p_all_sort[i];
            }
        }
        p_all[new_len] = unique_p;
        ++new_len;
        p_all_len = new_len;

        i_min_p_obs = 0;
        min_p_obs = p_obs[0];
        max_p_obs = p_obs[0];
        for (i=1; i<prof_len; ++i) {
            if (p_obs[i] < min_p_obs) {
                i_min_p_obs = i;
                min_p_obs = p_obs[i];
            }
            max_p_obs = max(max_p_obs, p_obs[i]);
        }

        i_obs_plus_interp = (int *) malloc(p_all_len*sizeof (int));
        i_surf_and_obs_plus_interp = (int *) malloc(p_all_len*sizeof (int));
        p_obs_plus_interp = (double *) malloc(p_all_len*sizeof (double));
        p_surf_and_obs_plus_interp = (double *) malloc(p_all_len*sizeof (double));
        i_out = (int *) malloc(p_all_len*sizeof (int));
        i_1 = (int *) malloc(p_all_len*sizeof (int));
        i_2 = (int *) malloc(p_all_len*sizeof (int));
        i_3 = (int *) malloc(p_all_len*sizeof (int));

        i_obs_plus_interp_len = 0;
        i_surf_and_obs_plus_interp_len = 0;
        for (i=0; i<p_all_len; ++i) {
            if (p_all[i] >= min_p_obs && p_all[i] <= max_p_obs) {
                i_obs_plus_interp[i_obs_plus_interp_len] = i;
                p_obs_plus_interp[i_obs_plus_interp_len] = p_all[i];
                i_obs_plus_interp_len++;
            }
            if (p_all[i] <= max_p_obs) {
                i_surf_and_obs_plus_interp[i_surf_and_obs_plus_interp_len] = i;
                p_surf_and_obs_plus_interp[i_surf_and_obs_plus_interp_len] = p_all[i];
                i_surf_and_obs_plus_interp_len++;
            }
        }
        i_out_len = gsw_util_intersect(p_i_tmp, m_i, p_surf_and_obs_plus_interp, i_surf_and_obs_plus_interp_len, i_out, i_1);
        i_2_len = gsw_util_intersect(p_obs, prof_len, p_obs_plus_interp, i_obs_plus_interp_len, i_2, i_3);

        independent_variable = (double *) malloc(prof_len*sizeof (double));
        independent_variable_obs_plus_interp = (double *) malloc(i_obs_plus_interp_len*sizeof (double));

        for(i=0; i<prof_len; ++i) {
            independent_variable[i] = (float)i;
        }
        gsw_util_pchip_interp(p_obs, independent_variable, prof_len,
                                p_obs_plus_interp, independent_variable_obs_plus_interp, i_obs_plus_interp_len);
        
        scaled_tracer_obs = (double *) malloc(prof_len*sizeof (double));
        v_tmp = (double *) malloc(prof_len*sizeof (double));
        q_tmp = (double *) malloc(prof_len*sizeof (double));
        v_i = (double *) malloc(i_obs_plus_interp_len*sizeof (double));
        q_i = (double *) malloc(i_obs_plus_interp_len*sizeof (double));
        tracer_i_obs_plus_interp = (double *)malloc(i_obs_plus_interp_len*sizeof (double));
        ct_i_obs_plus_interp = (double *)malloc(i_obs_plus_interp_len*sizeof (double));

        for(i=0; i<prof_len; ++i) {
            scaled_tracer_obs[i] = factor * tracer_obs[i];
        }

        gsw_util_pchip_interp(independent_variable, ct_obs, prof_len,
                                independent_variable_obs_plus_interp, ct_i_obs_plus_interp, i_obs_plus_interp_len);
        gsw_util_pchip_interp(independent_variable, scaled_tracer_obs, prof_len,
                                independent_variable_obs_plus_interp, q_i, i_obs_plus_interp_len);

        for(i=0; i<i_obs_plus_interp_len; ++i) {
            tracer_i_obs_plus_interp[i] = rec_factor * q_i[i];
        }

        for(k=0; k<7; ++k) {
            for(i=0; i<prof_len; ++i) {
                v_tmp[i] = scaled_tracer_obs[i] * sin_kpi_on_16[k] + ct_obs[i] * cos_kpi_on_16[k];
                q_tmp[i] = scaled_tracer_obs[i] * cos_kpi_on_16[k] - ct_obs[i] * sin_kpi_on_16[k];
            }
            gsw_util_pchip_interp(independent_variable, v_tmp, prof_len,
                                    independent_variable_obs_plus_interp, v_i, i_obs_plus_interp_len);
            gsw_util_pchip_interp(independent_variable, q_tmp, prof_len,
                                    independent_variable_obs_plus_interp, q_i, i_obs_plus_interp_len);
            for(i=0; i<i_obs_plus_interp_len; ++i) {
                ct_i_obs_plus_interp[i] += -q_i[i] * sin_kpi_on_16[k] + v_i[i] * cos_kpi_on_16[k];
                tracer_i_obs_plus_interp[i] += rec_factor * (q_i[i] * cos_kpi_on_16[k] + v_i[i] * sin_kpi_on_16[k]);
            }
        }

        for(i=0; i<i_obs_plus_interp_len; ++i) {
            tracer_i_obs_plus_interp[i] *= 0.125;
            ct_i_obs_plus_interp[i] *= 0.125;
        }

        tracer_i_tooutput = (double *)malloc(i_surf_and_obs_plus_interp_len*sizeof (double));
        ct_i_tooutput = (double *)malloc(i_surf_and_obs_plus_interp_len*sizeof (double));

        if (min_p_obs != 0.) {
            for (i=0; i<i_surf_and_obs_plus_interp_len; ++i) {
                tracer_i_tooutput[i] = NAN;
                ct_i_tooutput[i] = NAN;
                if (p_i_tmp[i] < min_p_obs) {
                    tracer_i_tooutput[i] = tracer_i_obs_plus_interp[i_3[i_min_p_obs]];
                    ct_i_tooutput[i] = ct_i_obs_plus_interp[i_3[i_min_p_obs]];
                }
            }
            for (i=0; i<i_obs_plus_interp_len; ++i) {
                tracer_i_tooutput[i_obs_plus_interp[i]] = tracer_i_obs_plus_interp[i];
                ct_i_tooutput[i_obs_plus_interp[i]] = ct_i_obs_plus_interp[i];
            }
        } else {
            for (i=0; i<i_obs_plus_interp_len; ++i) {
                tracer_i_tooutput[i] = tracer_i_obs_plus_interp[i];
                ct_i_tooutput[i] = ct_i_obs_plus_interp[i];
            }

        }

        for(i=0; i<i_out_len; ++i) {
            tracer_i[i_out[i]] = tracer_i_tooutput[i_1[i]];
            ct_i[i_out[i]] = ct_i_tooutput[i_1[i]];
        }

        free(independent_variable);
        free(independent_variable_obs_plus_interp);
        free(scaled_tracer_obs);
        free(v_tmp);
        free(q_tmp);
        free(v_i);
        free(q_i);
        free(tracer_i_obs_plus_interp);
        free(ct_i_obs_plus_interp);
        free(tracer_i_tooutput);
        free(ct_i_tooutput);
        free(tracer_obs);
        free(ct_obs);
        free(p_obs);
        free(p_i_tmp);
        free(p_idx);
        free(p_sort);
        free(tracer_sort);
        free(ct_sort);
        free(p_all);
        free(p_all_sort);
        free(p_all_idx);
        free(i_obs_plus_interp);
        free(i_surf_and_obs_plus_interp);
        free(p_obs_plus_interp);
        free(p_surf_and_obs_plus_interp);
        free(i_out);
        free(i_1);
        free(i_2);
        free(i_3);

        return;
}