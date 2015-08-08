/*
!==========================================================================
function gsw_sp_from_sk(sk)
!==========================================================================

! Calculates Practical Salinity, SP, from SK
!
!  SK    : Knudsen Salinity                        [parts per thousand, ppt]
!
! gsw_sp_from_sk  : Practical Salinity                              [unitless]
*/
double
gsw_sp_from_sk(double sk)
{
	GSW_TEOS10_CONSTANTS;
	double	gsw_sp_from_sk_value;


	gsw_sp_from_sk_value = (sk - 0.03e0)*(gsw_soncl/1.805e0);

	if (gsw_sp_from_sk_value < 0e0)
            gsw_sp_from_sk_value = GSW_INVALID_VALUE;

	return (gsw_sp_from_sk_value);
}
