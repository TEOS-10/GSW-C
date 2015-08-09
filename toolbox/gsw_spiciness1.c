/*
!==========================================================================
function gsw_spiciness1(sa,ct)  
!==========================================================================
!
!  Calculates spiciness from Absolute Salinity and Conservative 
!  Temperature at a pressure of 1000 dbar, as described by McDougall and 
!  Krzysik (2015).  This routine is based on the computationally-efficient 
!  expression for specific volume in terms of SA, CT and p (Roquet et al., 
!  2015).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  spiciness1  =  spiciness referenced to a pressure of 1000 dbar [ kg/m^3 ]
!
*/
double
gsw_spiciness1(double sa, double ct)
{
	GSW_TEOS10_CONSTANTS;
	double	s01 = -9.19874584868912e1,	s02 = -1.33517268529408e1,
		s03 =  2.18352211648107e1,	s04 = -2.01491744114173e1,
		s05 =  3.70004204355132e1,	s06 = -3.78831543226261e1,
		s07 =  1.76337834294554e1,	s08 =  2.87838842773396e2,
		s09 =  2.14531420554522e1,	s10 =  3.14679705198796e1,
		s11 = -4.04398864750692e1,	s12 = -7.70796428950487e1,
		s13 =  1.36783833820955e2,	s14 = -7.36834317044850e1,
		s15 = -6.41753415180701e2,	s16 =  1.33701981685590,
		s17 = -1.75289327948412e2,	s18 =  2.42666160657536e2,
		s19 =  3.17062400799114e1,	s20 = -2.28131490440865e2,
		s21 =  1.39564245068468e2,	s22 =  8.27747934506435e2,
		s23 = -3.50901590694775e1,	s24 =  2.87473907262029e2,
		s25 = -4.00227341144928e2,	s26 =  6.48307189919433e1,
		s27 =  2.16433334701578e2,	s28 = -1.48273032774305e2,
		s29 = -5.74545648799754e2,	s30 =  4.50446431127421e1,
		s31 = -2.30714981343772e2,	s32 =  3.15958389253065e2,
		s33 = -8.60635313930106e1,	s34 = -1.22978455069097e2,
		s35 =  9.18287282626261e1,	s36 =  2.12120473062203e2,
		s37 = -2.21528216973820e1,	s38 =  9.19013417923270e1,
		s39 = -1.24400776026014e2,	s40 =  4.08512871163839e1,
		s41 =  3.91127352213516e1,	s42 = -3.10508021853093e1,
		s43 = -3.24790035899152e1,	s44 =  3.91029016556786,
		s45 = -1.45362719385412e1,	s46 =  1.96136194246355e1,
		s47 = -7.06035474689088,	s48 = -5.36884688614009,
		s49 =  4.43247303092448;
	double 	xs, ys, spiciness1;
	
	xs 	= sqrt(gsw_sfac*sa + offset);
	ys 	= ct*0.025;
	spiciness1= s01+ys*(s02+ys*(s03+ys*(s04+ys*(s05+ys*(s06+s07*ys)))))
		+xs*(s08+ys*(s09+ys*(s10+ys*(s11+ys*(s12+ys*(s13+s14*ys)))))
		+xs*(s15+ys*(s16+ys*(s17+ys*(s18+ys*(s19+ys*(s20+s21*ys)))))
		+xs*(s22+ys*(s23+ys*(s24+ys*(s25+ys*(s26+ys*(s27+s28*ys)))))
		+xs*(s29+ys*(s30+ys*(s31+ys*(s32+ys*(s33+ys*(s34+s35*ys)))))
		+xs*(s36+ys*(s37+ys*(s38+ys*(s39+ys*(s40+ys*(s41+s42*ys)))))
		+xs*(s43+ys*(s44+ys*(s45+ys*(s46+ys*(s47+ys*(s48+s49*ys)))))
		))))));
	return (spiciness1);
}
