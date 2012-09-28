# This file was automatically created by FeynRules $Revision: 364 $
# Mathematica version: 7.0 for Mac OS X x86 (64-bit) (November 11, 2008)
# Date: Wed 10 Nov 2010 10:19:46


from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec



GC_1 = Coupling(name = 'GC_1',
                value = '-(ee*complex(0,1))/3.',
                order = {'QED':1})

GC_2 = Coupling(name = 'GC_2',
                value = '(2*ee*complex(0,1))/3.',
                order = {'QED':1})

GC_3 = Coupling(name = 'GC_3',
                value = '-(ee*complex(0,1))',
                order = {'QED':1})

GC_4 = Coupling(name = 'GC_4',
                value = '-G',
                order = {'QCD':1})

GC_5 = Coupling(name = 'GC_5',
                value = 'complex(0,1)*G',
                order = {'QCD':1})

GC_6 = Coupling(name = 'GC_6',
                value = 'complex(0,1)*G**2',
                order = {'QCD':2})

GC_7 = Coupling(name = 'GC_7',
                value = 'ca*cw*complex(0,1)*gw',
                order = {'QED':1})

GC_8 = Coupling(name = 'GC_8',
                value = '-(complex(0,1)*gw**2)',
                order = {'QED':2})

GC_9 = Coupling(name = 'GC_9',
                value = 'ca**2*cw**2*complex(0,1)*gw**2',
                order = {'QED':2})

GC_10 = Coupling(name = 'GC_10',
                 value = '-(cw*complex(0,1)*gw*sa)',
                 order = {'QED':1})

GC_11 = Coupling(name = 'GC_11',
                 value = '-(ca*cw**2*complex(0,1)*gw**2*sa)',
                 order = {'QED':2})

GC_12 = Coupling(name = 'GC_12',
                 value = 'cw**2*complex(0,1)*gw**2*sa**2',
                 order = {'QED':2})

GC_13 = Coupling(name = 'GC_13',
                 value = '-6*ch**3*complex(0,1)*kap*sh + 24*ch**3*complex(0,1)*rho*sh + 6*ch*complex(0,1)*kap*sh**3 - 6*ch*complex(0,1)*l*sh**3',
                 order = {'QED':2})

GC_14 = Coupling(name = 'GC_14',
                 value = '6*ch**3*complex(0,1)*kap*sh - 6*ch**3*complex(0,1)*l*sh - 6*ch*complex(0,1)*kap*sh**3 + 24*ch*complex(0,1)*rho*sh**3',
                 order = {'QED':2})

GC_15 = Coupling(name = 'GC_15',
                 value = '-2*ch**4*complex(0,1)*kap + 8*ch**2*complex(0,1)*kap*sh**2 - 6*ch**2*complex(0,1)*l*sh**2 - 24*ch**2*complex(0,1)*rho*sh**2 - 2*complex(0,1)*kap*sh**4',
                 order = {'QED':2})

GC_16 = Coupling(name = 'GC_16',
                 value = '-24*ch**4*complex(0,1)*rho - 12*ch**2*complex(0,1)*kap*sh**2 - 6*complex(0,1)*l*sh**4',
                 order = {'QED':2})

GC_17 = Coupling(name = 'GC_17',
                 value = '-6*ch**4*complex(0,1)*l - 12*ch**2*complex(0,1)*kap*sh**2 - 24*complex(0,1)*rho*sh**4',
                 order = {'QED':2})

GC_18 = Coupling(name = 'GC_18',
                 value = '(ch**2*ee**2*complex(0,1))/(2.*sw**2)',
                 order = {'QED':2})

GC_19 = Coupling(name = 'GC_19',
                 value = '(ch*ee**2*complex(0,1)*sh)/(2.*sw**2)',
                 order = {'QED':2})

GC_20 = Coupling(name = 'GC_20',
                 value = '(ee**2*complex(0,1)*sh**2)/(2.*sw**2)',
                 order = {'QED':2})

GC_21 = Coupling(name = 'GC_21',
                 value = '(ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_22 = Coupling(name = 'GC_22',
                 value = '(CKM11*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_23 = Coupling(name = 'GC_23',
                 value = '(CKM12*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_24 = Coupling(name = 'GC_24',
                 value = '(CKM21*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_25 = Coupling(name = 'GC_25',
                 value = '(CKM22*ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_26 = Coupling(name = 'GC_26',
                 value = '-(ca*cw*ee*complex(0,1))/(2.*sw)',
                 order = {'QED':1})

GC_27 = Coupling(name = 'GC_27',
                 value = '(ca*cw*ee*complex(0,1))/(2.*sw)',
                 order = {'QED':1})

GC_28 = Coupling(name = 'GC_28',
                 value = '-(cw*ee*complex(0,1)*sa)/(2.*sw)',
                 order = {'QED':1})

GC_29 = Coupling(name = 'GC_29',
                 value = '(cw*ee*complex(0,1)*sa)/(2.*sw)',
                 order = {'QED':1})

GC_30 = Coupling(name = 'GC_30',
                 value = 'complex(0,1)*gw*sw',
                 order = {'QED':1})

GC_31 = Coupling(name = 'GC_31',
                 value = '-2*ca*cw*complex(0,1)*gw**2*sw',
                 order = {'QED':2})

GC_32 = Coupling(name = 'GC_32',
                 value = '2*cw*complex(0,1)*gw**2*sa*sw',
                 order = {'QED':2})

GC_33 = Coupling(name = 'GC_33',
                 value = 'complex(0,1)*gw**2*sw**2',
                 order = {'QED':2})

GC_34 = Coupling(name = 'GC_34',
                 value = '(ee*eta*complex(0,1)*sa)/(6.*cw) - (ca*ee*complex(0,1)*sw)/(6.*cw)',
                 order = {'QED':1})

GC_35 = Coupling(name = 'GC_35',
                 value = '-(ee*eta*complex(0,1)*sa)/(2.*cw) + (ca*ee*complex(0,1)*sw)/(2.*cw)',
                 order = {'QED':1})

GC_36 = Coupling(name = 'GC_36',
                 value = '-(ee*eta*complex(0,1)*sa)/(2.*cw) + (ca*cw*ee*complex(0,1))/(2.*sw) + (ca*ee*complex(0,1)*sw)/(2.*cw)',
                 order = {'QED':1})

GC_37 = Coupling(name = 'GC_37',
                 value = '(ca*ee*eta*complex(0,1))/(6.*cw) + (ee*complex(0,1)*sa*sw)/(6.*cw)',
                 order = {'QED':1})

GC_38 = Coupling(name = 'GC_38',
                 value = '-(ca*ee*eta*complex(0,1))/(2.*cw) - (ee*complex(0,1)*sa*sw)/(2.*cw)',
                 order = {'QED':1})

GC_39 = Coupling(name = 'GC_39',
                 value = '-(ca*ee*eta*complex(0,1))/(2.*cw) - (cw*ee*complex(0,1)*sa)/(2.*sw) - (ee*complex(0,1)*sa*sw)/(2.*cw)',
                 order = {'QED':1})

GC_40 = Coupling(name = 'GC_40',
                 value = 'ca**2*ch**2*ee**2*complex(0,1) + (ch**2*ee**2*eta**2*complex(0,1)*sa**2)/(2.*cw**2) + 4*chi**2*eta**2*complex(0,1)*gX**2*sa**2*sh**2 + (ca**2*ch**2*cw**2*ee**2*complex(0,1))/(2.*sw**2) - (ca*ch**2*ee**2*eta*complex(0,1)*sa)/sw - (ca*ch**2*ee**2*eta*complex(0,1)*sa*sw)/cw**2 + (ca**2*ch**2*ee**2*complex(0,1)*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_41 = Coupling(name = 'GC_41',
                 value = '-(ca*ch**2*ee**2*complex(0,1)*sa) + (ca*ch**2*ee**2*eta**2*complex(0,1)*sa)/(2.*cw**2) + 4*ca*chi**2*eta**2*complex(0,1)*gX**2*sa*sh**2 - (ca*ch**2*cw**2*ee**2*complex(0,1)*sa)/(2.*sw**2) - (ca**2*ch**2*ee**2*eta*complex(0,1))/(2.*sw) + (ch**2*ee**2*eta*complex(0,1)*sa**2)/(2.*sw) - (ca**2*ch**2*ee**2*eta*complex(0,1)*sw)/(2.*cw**2) + (ch**2*ee**2*eta*complex(0,1)*sa**2*sw)/(2.*cw**2) - (ca*ch**2*ee**2*complex(0,1)*sa*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_42 = Coupling(name = 'GC_42',
                 value = '(ca**2*ch**2*ee**2*eta**2*complex(0,1))/(2.*cw**2) + ch**2*ee**2*complex(0,1)*sa**2 + 4*ca**2*chi**2*eta**2*complex(0,1)*gX**2*sh**2 + (ch**2*cw**2*ee**2*complex(0,1)*sa**2)/(2.*sw**2) + (ca*ch**2*ee**2*eta*complex(0,1)*sa)/sw + (ca*ch**2*ee**2*eta*complex(0,1)*sa*sw)/cw**2 + (ch**2*ee**2*complex(0,1)*sa**2*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_43 = Coupling(name = 'GC_43',
                 value = 'ca**2*ch*ee**2*complex(0,1)*sh + (ch*ee**2*eta**2*complex(0,1)*sa**2*sh)/(2.*cw**2) - 4*ch*chi**2*eta**2*complex(0,1)*gX**2*sa**2*sh + (ca**2*ch*cw**2*ee**2*complex(0,1)*sh)/(2.*sw**2) - (ca*ch*ee**2*eta*complex(0,1)*sa*sh)/sw - (ca*ch*ee**2*eta*complex(0,1)*sa*sh*sw)/cw**2 + (ca**2*ch*ee**2*complex(0,1)*sh*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_44 = Coupling(name = 'GC_44',
                 value = '-(ca*ch*ee**2*complex(0,1)*sa*sh) + (ca*ch*ee**2*eta**2*complex(0,1)*sa*sh)/(2.*cw**2) - 4*ca*ch*chi**2*eta**2*complex(0,1)*gX**2*sa*sh - (ca*ch*cw**2*ee**2*complex(0,1)*sa*sh)/(2.*sw**2) - (ca**2*ch*ee**2*eta*complex(0,1)*sh)/(2.*sw) + (ch*ee**2*eta*complex(0,1)*sa**2*sh)/(2.*sw) - (ca**2*ch*ee**2*eta*complex(0,1)*sh*sw)/(2.*cw**2) + (ch*ee**2*eta*complex(0,1)*sa**2*sh*sw)/(2.*cw**2) - (ca*ch*ee**2*complex(0,1)*sa*sh*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_45 = Coupling(name = 'GC_45',
                 value = '(ca**2*ch*ee**2*eta**2*complex(0,1)*sh)/(2.*cw**2) - 4*ca**2*ch*chi**2*eta**2*complex(0,1)*gX**2*sh + ch*ee**2*complex(0,1)*sa**2*sh + (ch*cw**2*ee**2*complex(0,1)*sa**2*sh)/(2.*sw**2) + (ca*ch*ee**2*eta*complex(0,1)*sa*sh)/sw + (ca*ch*ee**2*eta*complex(0,1)*sa*sh*sw)/cw**2 + (ch*ee**2*complex(0,1)*sa**2*sh*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_46 = Coupling(name = 'GC_46',
                 value = '4*ch**2*chi**2*eta**2*complex(0,1)*gX**2*sa**2 + ca**2*ee**2*complex(0,1)*sh**2 + (ee**2*eta**2*complex(0,1)*sa**2*sh**2)/(2.*cw**2) + (ca**2*cw**2*ee**2*complex(0,1)*sh**2)/(2.*sw**2) - (ca*ee**2*eta*complex(0,1)*sa*sh**2)/sw - (ca*ee**2*eta*complex(0,1)*sa*sh**2*sw)/cw**2 + (ca**2*ee**2*complex(0,1)*sh**2*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_47 = Coupling(name = 'GC_47',
                 value = '4*ca*ch**2*chi**2*eta**2*complex(0,1)*gX**2*sa - ca*ee**2*complex(0,1)*sa*sh**2 + (ca*ee**2*eta**2*complex(0,1)*sa*sh**2)/(2.*cw**2) - (ca*cw**2*ee**2*complex(0,1)*sa*sh**2)/(2.*sw**2) - (ca**2*ee**2*eta*complex(0,1)*sh**2)/(2.*sw) + (ee**2*eta*complex(0,1)*sa**2*sh**2)/(2.*sw) - (ca**2*ee**2*eta*complex(0,1)*sh**2*sw)/(2.*cw**2) + (ee**2*eta*complex(0,1)*sa**2*sh**2*sw)/(2.*cw**2) - (ca*ee**2*complex(0,1)*sa*sh**2*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_48 = Coupling(name = 'GC_48',
                 value = '4*ca**2*ch**2*chi**2*eta**2*complex(0,1)*gX**2 + (ca**2*ee**2*eta**2*complex(0,1)*sh**2)/(2.*cw**2) + ee**2*complex(0,1)*sa**2*sh**2 + (cw**2*ee**2*complex(0,1)*sa**2*sh**2)/(2.*sw**2) + (ca*ee**2*eta*complex(0,1)*sa*sh**2)/sw + (ca*ee**2*eta*complex(0,1)*sa*sh**2*sw)/cw**2 + (ee**2*complex(0,1)*sa**2*sh**2*sw**2)/(2.*cw**2)',
                 order = {'QED':2})

GC_49 = Coupling(name = 'GC_49',
                 value = '(ch*ee**2*complex(0,1)*v)/(2.*sw**2)',
                 order = {'QED':1})

GC_50 = Coupling(name = 'GC_50',
                 value = '(ee**2*complex(0,1)*sh*v)/(2.*sw**2)',
                 order = {'QED':1})

GC_51 = Coupling(name = 'GC_51',
                 value = '(ca**2*ee**2*eta**2*complex(0,1)*sh*v)/(2.*cw**2) + ee**2*complex(0,1)*sa**2*sh*v + (cw**2*ee**2*complex(0,1)*sa**2*sh*v)/(2.*sw**2) + (ca*ee**2*eta*complex(0,1)*sa*sh*v)/sw + (ca*ee**2*eta*complex(0,1)*sa*sh*sw*v)/cw**2 + (ee**2*complex(0,1)*sa**2*sh*sw**2*v)/(2.*cw**2) + 2*ca**2*ch*chi**2*eta**2*complex(0,1)*gX**2*xi*cmath.sqrt(2)',
                 order = {'QED':1})

GC_52 = Coupling(name = 'GC_52',
                 value = '-(ca*ee**2*complex(0,1)*sa*sh*v) + (ca*ee**2*eta**2*complex(0,1)*sa*sh*v)/(2.*cw**2) - (ca*cw**2*ee**2*complex(0,1)*sa*sh*v)/(2.*sw**2) - (ca**2*ee**2*eta*complex(0,1)*sh*v)/(2.*sw) + (ee**2*eta*complex(0,1)*sa**2*sh*v)/(2.*sw) - (ca**2*ee**2*eta*complex(0,1)*sh*sw*v)/(2.*cw**2) + (ee**2*eta*complex(0,1)*sa**2*sh*sw*v)/(2.*cw**2) - (ca*ee**2*complex(0,1)*sa*sh*sw**2*v)/(2.*cw**2) + 2*ca*ch*chi**2*eta**2*complex(0,1)*gX**2*sa*xi*cmath.sqrt(2)',
                 order = {'QED':1})

GC_53 = Coupling(name = 'GC_53',
                 value = 'ca**2*ee**2*complex(0,1)*sh*v + (ee**2*eta**2*complex(0,1)*sa**2*sh*v)/(2.*cw**2) + (ca**2*cw**2*ee**2*complex(0,1)*sh*v)/(2.*sw**2) - (ca*ee**2*eta*complex(0,1)*sa*sh*v)/sw - (ca*ee**2*eta*complex(0,1)*sa*sh*sw*v)/cw**2 + (ca**2*ee**2*complex(0,1)*sh*sw**2*v)/(2.*cw**2) + 2*ch*chi**2*eta**2*complex(0,1)*gX**2*sa**2*xi*cmath.sqrt(2)',
                 order = {'QED':1})

GC_54 = Coupling(name = 'GC_54',
                 value = '(ca**2*ch*ee**2*eta**2*complex(0,1)*v)/(2.*cw**2) + ch*ee**2*complex(0,1)*sa**2*v + (ch*cw**2*ee**2*complex(0,1)*sa**2*v)/(2.*sw**2) + (ca*ch*ee**2*eta*complex(0,1)*sa*v)/sw + (ca*ch*ee**2*eta*complex(0,1)*sa*sw*v)/cw**2 + (ch*ee**2*complex(0,1)*sa**2*sw**2*v)/(2.*cw**2) - 2*ca**2*chi**2*eta**2*complex(0,1)*gX**2*sh*xi*cmath.sqrt(2)',
                 order = {'QED':1})

GC_55 = Coupling(name = 'GC_55',
                 value = '-(ca*ch*ee**2*complex(0,1)*sa*v) + (ca*ch*ee**2*eta**2*complex(0,1)*sa*v)/(2.*cw**2) - (ca*ch*cw**2*ee**2*complex(0,1)*sa*v)/(2.*sw**2) - (ca**2*ch*ee**2*eta*complex(0,1)*v)/(2.*sw) + (ch*ee**2*eta*complex(0,1)*sa**2*v)/(2.*sw) - (ca**2*ch*ee**2*eta*complex(0,1)*sw*v)/(2.*cw**2) + (ch*ee**2*eta*complex(0,1)*sa**2*sw*v)/(2.*cw**2) - (ca*ch*ee**2*complex(0,1)*sa*sw**2*v)/(2.*cw**2) - 2*ca*chi**2*eta**2*complex(0,1)*gX**2*sa*sh*xi*cmath.sqrt(2)',
                 order = {'QED':1})

GC_56 = Coupling(name = 'GC_56',
                 value = 'ca**2*ch*ee**2*complex(0,1)*v + (ch*ee**2*eta**2*complex(0,1)*sa**2*v)/(2.*cw**2) + (ca**2*ch*cw**2*ee**2*complex(0,1)*v)/(2.*sw**2) - (ca*ch*ee**2*eta*complex(0,1)*sa*v)/sw - (ca*ch*ee**2*eta*complex(0,1)*sa*sw*v)/cw**2 + (ca**2*ch*ee**2*complex(0,1)*sw**2*v)/(2.*cw**2) - 2*chi**2*eta**2*complex(0,1)*gX**2*sa**2*sh*xi*cmath.sqrt(2)',
                 order = {'QED':1})

GC_57 = Coupling(name = 'GC_57',
                 value = '-6*ch**2*complex(0,1)*kap*sh*v - 6*complex(0,1)*l*sh**3*v - 12*ch**3*complex(0,1)*rho*xi*cmath.sqrt(2) - 3*ch*complex(0,1)*kap*sh**2*xi*cmath.sqrt(2)',
                 order = {'QED':1})

GC_58 = Coupling(name = 'GC_58',
                 value = '4*ch**2*complex(0,1)*kap*sh*v - 6*ch**2*complex(0,1)*l*sh*v - 2*complex(0,1)*kap*sh**3*v - ch**3*complex(0,1)*kap*xi*cmath.sqrt(2) + 2*ch*complex(0,1)*kap*sh**2*xi*cmath.sqrt(2) - 12*ch*complex(0,1)*rho*sh**2*xi*cmath.sqrt(2)',
                 order = {'QED':1})

GC_59 = Coupling(name = 'GC_59',
                 value = '-2*ch**3*complex(0,1)*kap*v + 4*ch*complex(0,1)*kap*sh**2*v - 6*ch*complex(0,1)*l*sh**2*v - 2*ch**2*complex(0,1)*kap*sh*xi*cmath.sqrt(2) + 12*ch**2*complex(0,1)*rho*sh*xi*cmath.sqrt(2) + complex(0,1)*kap*sh**3*xi*cmath.sqrt(2)',
                 order = {'QED':1})

GC_60 = Coupling(name = 'GC_60',
                 value = '-6*ch**3*complex(0,1)*l*v - 6*ch*complex(0,1)*kap*sh**2*v + 3*ch**2*complex(0,1)*kap*sh*xi*cmath.sqrt(2) + 12*complex(0,1)*rho*sh**3*xi*cmath.sqrt(2)',
                 order = {'QED':1})

GC_61 = Coupling(name = 'GC_61',
                 value = '-((ch*complex(0,1)*yb)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_62 = Coupling(name = 'GC_62',
                 value = '-((complex(0,1)*sh*yb)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_63 = Coupling(name = 'GC_63',
                 value = '-((ch*complex(0,1)*yc)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_64 = Coupling(name = 'GC_64',
                 value = '-((complex(0,1)*sh*yc)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_65 = Coupling(name = 'GC_65',
                 value = '-((ch*complex(0,1)*yt)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_66 = Coupling(name = 'GC_66',
                 value = '-((complex(0,1)*sh*yt)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_67 = Coupling(name = 'GC_67',
                 value = '-((ch*complex(0,1)*ytau)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_68 = Coupling(name = 'GC_68',
                 value = '-((complex(0,1)*sh*ytau)/cmath.sqrt(2))',
                 order = {'QED':1})

GC_69 = Coupling(name = 'GC_69',
                 value = '(ee*complex(0,1)*complexconjugate(CKM11))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_70 = Coupling(name = 'GC_70',
                 value = '(ee*complex(0,1)*complexconjugate(CKM12))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_71 = Coupling(name = 'GC_71',
                 value = '(ee*complex(0,1)*complexconjugate(CKM21))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_72 = Coupling(name = 'GC_72',
                 value = '(ee*complex(0,1)*complexconjugate(CKM22))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

