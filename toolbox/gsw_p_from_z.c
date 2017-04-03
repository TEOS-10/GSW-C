/*
!==========================================================================
function gsw_p_from_z(z,lat)
!==========================================================================

! Calculates the pressure p from height z
!
! z      : height                                          [m]
! lat    : latitude                                        [deg]
!
! gsw_p_from_z : pressure                                  [dbar]
*/
double
gsw_p_from_z(double z, double lat)
{
    GSW_TEOS10_CONSTANTS;
    double sinlat, sin2, gs, c1, p, df_dp, f, p_old, p_mid;

    if (z > 5) return GSW_INVALID_VALUE;

    sinlat = sin(lat*deg2rad);
    sin2 = sinlat*sinlat;
    gs = 9.780327*(1.0 + (5.2792e-3 + (2.32e-5*sin2))*sin2);

    /* get the first estimate of p from Saunders (1981) */
    c1 =  5.25e-3*sin2 + 5.92e-3;
    p  = -2.0*z/((1-c1) + sqrt((1-c1)*(1-c1) + 8.84e-6*z)) ;
    /* end of the first estimate from Saunders (1981) */

    df_dp = db2pa*gsw_specvol_sso_0(p); /* initial value of the derivative of f */

    f = gsw_enthalpy_sso_0(p) + gs*(z - 0.5*gamma*(z*z));
             /*   - (geo_strf_dyn_height + sea_surface_geopotental); */
    p_old = p;
    p = p_old - f/df_dp;
    p_mid = 0.5*(p + p_old);
    df_dp = db2pa*gsw_specvol_sso_0(p_mid);
    p = p_old - f/df_dp;

    return p;
}
