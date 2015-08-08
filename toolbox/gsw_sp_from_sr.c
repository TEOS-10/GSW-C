/*
!==========================================================================
function gsw_sp_from_sr(sr)  
!==========================================================================

! Calculates Practical Salinity, sp, from Reference Salinity, sr. 
!
! sr     : Reference Salinity                              [g/kg]
!
! gsw_sp_from_sr  : Practical Salinity                     [unitless]
*/
double
gsw_sp_from_sr(double sr)
{
	GSW_TEOS10_CONSTANTS;

	return(sr/gsw_ups);
}
