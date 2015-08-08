/*
!==========================================================================
function gsw_sr_from_sp(sp)  
!==========================================================================

! Calculates Reference Salinity, SR, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
!
! gsw_sr_from_sp : Reference Salinity                      [g/kg]
*/
double
gsw_sr_from_sp(double sp)
{
	GSW_TEOS10_CONSTANTS;

	return (sp*gsw_ups);
}
