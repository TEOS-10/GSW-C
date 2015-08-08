/*
!==========================================================================
function gsw_sp_from_sstar(sstar,p,lon,lat)  
!==========================================================================

! Calculates Practical Salinity, SP, from Preformed Salinity, Sstar. 
!
! sstar  : Preformed Salinity                              [g/kg]
! p      : sea pressure                                    [dbar]
! lon   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sp_from_Sstar : Preformed Salinity                   [g/kg]
*/
double
gsw_sp_from_sstar(double sstar, double p, double lon, double lat)
{
	GSW_TEOS10_CONSTANTS;
	double	saar, sp_baltic;

    /*
    **! In the Baltic Sea, SA = Sstar.
    */
	sp_baltic	= gsw_sp_from_sa_baltic(sstar,lon,lat);
	if (sp_baltic < GSW_ERROR_LIMIT)
	    return (sp_baltic);
	saar	= gsw_saar(p,lon,lat);
	if (saar == GSW_INVALID_VALUE)
	    return (saar);
	return ((sstar/gsw_ups)/(1.0 - 0.35e0*saar));
}
