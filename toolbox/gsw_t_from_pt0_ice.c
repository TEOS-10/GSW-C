/*
! =========================================================================
elemental function gsw_t_from_pt0_ice (pt0_ice, p)
! =========================================================================
!
!  Calculates in-situ temperature from the potential temperature of ice Ih 
!  with reference pressure, p_ref, of 0 dbar (the surface), and the 
!  in-situ pressure.
!
!  pt0_ice  =  potential temperature of ice Ih with reference pressure of 
!              zero dbar (ITS-90)                                 [ deg C ]
!  p        =  sea pressure                                        [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!--------------------------------------------------------------------------
*/
double
gsw_t_from_pt0_ice(double pt0_ice, double p)
{
	double	p0 = 0.0;

	return (gsw_pt_from_t_ice(pt0_ice,p0,p));
}
