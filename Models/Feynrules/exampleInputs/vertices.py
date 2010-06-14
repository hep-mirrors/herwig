# This file was automatically created by FeynRules $Revision: 189 $
# Mathematica version: 7.0 for Linux x86 (32-bit) (November 11, 2008)
# Date: Thu 10 Jun 2010 10:33:50


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L


V_1 = Vertex(name = 'V_1',
             particles = [ P.G, P.G, P.G ],
             color = [ 'f(1,2,3)' ],
             lorentz = [ L.L_9 ],
             couplings = {(0,0):C.GC_5})

V_2 = Vertex(name = 'V_2',
             particles = [ P.G, P.G, P.G, P.G ],
             color = [ 'f(2,3,a1)*f(a1,1,4)', 'f(2,4,a1)*f(a1,1,3)', 'f(3,4,a1)*f(a1,1,2)' ],
             lorentz = [ L.L_12, L.L_14, L.L_15 ],
             couplings = {(1,1):C.GC_6,(2,0):C.GC_6,(0,2):C.GC_6})

V_3 = Vertex(name = 'V_3',
             particles = [ P.A, P.W__plus__, P.W__minus__ ],
             color = [ '1' ],
             lorentz = [ L.L_9 ],
             couplings = {(0,0):C.GC_30})

V_4 = Vertex(name = 'V_4',
             particles = [ P.A, P.A, P.W__plus__, P.W__minus__ ],
             color = [ '1' ],
             lorentz = [ L.L_13 ],
             couplings = {(0,0):C.GC_33})

V_5 = Vertex(name = 'V_5',
             particles = [ P.W__plus__, P.W__minus__, P.Z ],
             color = [ '1' ],
             lorentz = [ L.L_9 ],
             couplings = {(0,0):C.GC_7})

V_6 = Vertex(name = 'V_6',
             particles = [ P.W__plus__, P.W__minus__, P.Zp ],
             color = [ '1' ],
             lorentz = [ L.L_9 ],
             couplings = {(0,0):C.GC_10})

V_7 = Vertex(name = 'V_7',
             particles = [ P.W__plus__, P.W__plus__, P.W__minus__, P.W__minus__ ],
             color = [ '1' ],
             lorentz = [ L.L_13 ],
             couplings = {(0,0):C.GC_8})

V_8 = Vertex(name = 'V_8',
             particles = [ P.A, P.W__plus__, P.W__minus__, P.Z ],
             color = [ '1' ],
             lorentz = [ L.L_16 ],
             couplings = {(0,0):C.GC_31})

V_9 = Vertex(name = 'V_9',
             particles = [ P.W__plus__, P.W__minus__, P.Z, P.Z ],
             color = [ '1' ],
             lorentz = [ L.L_13 ],
             couplings = {(0,0):C.GC_9})

V_10 = Vertex(name = 'V_10',
              particles = [ P.A, P.W__plus__, P.W__minus__, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_16 ],
              couplings = {(0,0):C.GC_32})

V_11 = Vertex(name = 'V_11',
              particles = [ P.W__plus__, P.W__minus__, P.Z, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_13 ],
              couplings = {(0,0):C.GC_11})

V_12 = Vertex(name = 'V_12',
              particles = [ P.W__plus__, P.W__minus__, P.Zp, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_13 ],
              couplings = {(0,0):C.GC_12})

V_13 = Vertex(name = 'V_13',
              particles = [ P.h1, P.h1, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_10 ],
              couplings = {(0,0):C.GC_15})

V_14 = Vertex(name = 'V_14',
              particles = [ P.h1, P.h1, P.h1, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_10 ],
              couplings = {(0,0):C.GC_14})

V_15 = Vertex(name = 'V_15',
              particles = [ P.h1, P.h2, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_10 ],
              couplings = {(0,0):C.GC_13})

V_16 = Vertex(name = 'V_16',
              particles = [ P.h1, P.h1, P.h1, P.h1 ],
              color = [ '1' ],
              lorentz = [ L.L_10 ],
              couplings = {(0,0):C.GC_17})

V_17 = Vertex(name = 'V_17',
              particles = [ P.h2, P.h2, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_10 ],
              couplings = {(0,0):C.GC_16})

V_18 = Vertex(name = 'V_18',
              particles = [ P.h1, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_1 ],
              couplings = {(0,0):C.GC_59})

V_19 = Vertex(name = 'V_19',
              particles = [ P.h1, P.h1, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_1 ],
              couplings = {(0,0):C.GC_58})

V_20 = Vertex(name = 'V_20',
              particles = [ P.h2, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_1 ],
              couplings = {(0,0):C.GC_57})

V_21 = Vertex(name = 'V_21',
              particles = [ P.h1, P.h1, P.h1 ],
              color = [ '1' ],
              lorentz = [ L.L_1 ],
              couplings = {(0,0):C.GC_60})

V_22 = Vertex(name = 'V_22',
              particles = [ P.h1, P.h1, P.W__plus__, P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_18})

V_23 = Vertex(name = 'V_23',
              particles = [ P.h1, P.h2, P.W__plus__, P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_19})

V_24 = Vertex(name = 'V_24',
              particles = [ P.h2, P.h2, P.W__plus__, P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_20})

V_25 = Vertex(name = 'V_25',
              particles = [ P.h1, P.W__plus__, P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_49})

V_26 = Vertex(name = 'V_26',
              particles = [ P.h2, P.W__plus__, P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_50})

V_27 = Vertex(name = 'V_27',
              particles = [ P.h1, P.h1, P.Z, P.Z ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_40})

V_28 = Vertex(name = 'V_28',
              particles = [ P.h1, P.h2, P.Z, P.Z ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_43})

V_29 = Vertex(name = 'V_29',
              particles = [ P.h2, P.h2, P.Z, P.Z ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_46})

V_30 = Vertex(name = 'V_30',
              particles = [ P.h1, P.Z, P.Z ],
              color = [ '1' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_56})

V_31 = Vertex(name = 'V_31',
              particles = [ P.h2, P.Z, P.Z ],
              color = [ '1' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_53})

V_32 = Vertex(name = 'V_32',
              particles = [ P.h1, P.h1, P.Z, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_41})

V_33 = Vertex(name = 'V_33',
              particles = [ P.h1, P.h2, P.Z, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_44})

V_34 = Vertex(name = 'V_34',
              particles = [ P.h2, P.h2, P.Z, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_47})

V_35 = Vertex(name = 'V_35',
              particles = [ P.h1, P.Z, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_55})

V_36 = Vertex(name = 'V_36',
              particles = [ P.h2, P.Z, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_52})

V_37 = Vertex(name = 'V_37',
              particles = [ P.h1, P.h1, P.Zp, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_42})

V_38 = Vertex(name = 'V_38',
              particles = [ P.h1, P.h2, P.Zp, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_45})

V_39 = Vertex(name = 'V_39',
              particles = [ P.h2, P.h2, P.Zp, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_48})

V_40 = Vertex(name = 'V_40',
              particles = [ P.h1, P.Zp, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_54})

V_41 = Vertex(name = 'V_41',
              particles = [ P.h2, P.Zp, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_51})

V_42 = Vertex(name = 'V_42',
              particles = [ P.A, P.d__tilde__, P.d ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_1})

V_43 = Vertex(name = 'V_43',
              particles = [ P.A, P.s__tilde__, P.s ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_1})

V_44 = Vertex(name = 'V_44',
              particles = [ P.A, P.b__tilde__, P.b ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_1})

V_45 = Vertex(name = 'V_45',
              particles = [ P.A, P.e__plus__, P.e__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_3})

V_46 = Vertex(name = 'V_46',
              particles = [ P.A, P.m__plus__, P.m__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_3})

V_47 = Vertex(name = 'V_47',
              particles = [ P.A, P.tt__plus__, P.tt__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_3})

V_48 = Vertex(name = 'V_48',
              particles = [ P.A, P.u__tilde__, P.u ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_2})

V_49 = Vertex(name = 'V_49',
              particles = [ P.A, P.c__tilde__, P.c ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_2})

V_50 = Vertex(name = 'V_50',
              particles = [ P.A, P.t__tilde__, P.t ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_2})

V_51 = Vertex(name = 'V_51',
              particles = [ P.G, P.d__tilde__, P.d ],
              color = [ 'T(1,2,3)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_4})

V_52 = Vertex(name = 'V_52',
              particles = [ P.G, P.s__tilde__, P.s ],
              color = [ 'T(1,2,3)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_4})

V_53 = Vertex(name = 'V_53',
              particles = [ P.G, P.b__tilde__, P.b ],
              color = [ 'T(1,2,3)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_4})

V_54 = Vertex(name = 'V_54',
              particles = [ P.h1, P.b__tilde__, P.b ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_61})

V_55 = Vertex(name = 'V_55',
              particles = [ P.h2, P.b__tilde__, P.b ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_62})

V_56 = Vertex(name = 'V_56',
              particles = [ P.Z, P.d__tilde__, P.d ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5, L.L_6 ],
              couplings = {(0,1):C.GC_34,(0,0):C.GC_26})

V_57 = Vertex(name = 'V_57',
              particles = [ P.Z, P.s__tilde__, P.s ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5, L.L_6 ],
              couplings = {(0,1):C.GC_34,(0,0):C.GC_26})

V_58 = Vertex(name = 'V_58',
              particles = [ P.Z, P.b__tilde__, P.b ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5, L.L_6 ],
              couplings = {(0,1):C.GC_34,(0,0):C.GC_26})

V_59 = Vertex(name = 'V_59',
              particles = [ P.Zp, P.d__tilde__, P.d ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5, L.L_6 ],
              couplings = {(0,1):C.GC_37,(0,0):C.GC_29})

V_60 = Vertex(name = 'V_60',
              particles = [ P.Zp, P.s__tilde__, P.s ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5, L.L_6 ],
              couplings = {(0,1):C.GC_37,(0,0):C.GC_29})

V_61 = Vertex(name = 'V_61',
              particles = [ P.Zp, P.b__tilde__, P.b ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5, L.L_6 ],
              couplings = {(0,1):C.GC_37,(0,0):C.GC_29})

V_62 = Vertex(name = 'V_62',
              particles = [ P.W__plus__, P.u__tilde__, P.d ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_22})

V_63 = Vertex(name = 'V_63',
              particles = [ P.W__plus__, P.c__tilde__, P.d ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_24})

V_64 = Vertex(name = 'V_64',
              particles = [ P.W__plus__, P.u__tilde__, P.s ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_23})

V_65 = Vertex(name = 'V_65',
              particles = [ P.W__plus__, P.c__tilde__, P.s ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_25})

V_66 = Vertex(name = 'V_66',
              particles = [ P.W__plus__, P.t__tilde__, P.b ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_21})

V_67 = Vertex(name = 'V_67',
              particles = [ P.W__minus__, P.d__tilde__, P.u ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_69})

V_68 = Vertex(name = 'V_68',
              particles = [ P.W__minus__, P.d__tilde__, P.c ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_71})

V_69 = Vertex(name = 'V_69',
              particles = [ P.W__minus__, P.s__tilde__, P.u ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_70})

V_70 = Vertex(name = 'V_70',
              particles = [ P.W__minus__, P.s__tilde__, P.c ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_72})

V_71 = Vertex(name = 'V_71',
              particles = [ P.W__minus__, P.b__tilde__, P.t ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_21})

V_72 = Vertex(name = 'V_72',
              particles = [ P.G, P.u__tilde__, P.u ],
              color = [ 'T(1,2,3)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_4})

V_73 = Vertex(name = 'V_73',
              particles = [ P.G, P.c__tilde__, P.c ],
              color = [ 'T(1,2,3)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_4})

V_74 = Vertex(name = 'V_74',
              particles = [ P.G, P.t__tilde__, P.t ],
              color = [ 'T(1,2,3)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_4})

V_75 = Vertex(name = 'V_75',
              particles = [ P.h1, P.tt__plus__, P.tt__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_67})

V_76 = Vertex(name = 'V_76',
              particles = [ P.h1, P.c__tilde__, P.c ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_63})

V_77 = Vertex(name = 'V_77',
              particles = [ P.h1, P.t__tilde__, P.t ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_65})

V_78 = Vertex(name = 'V_78',
              particles = [ P.h2, P.tt__plus__, P.tt__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_68})

V_79 = Vertex(name = 'V_79',
              particles = [ P.h2, P.c__tilde__, P.c ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_64})

V_80 = Vertex(name = 'V_80',
              particles = [ P.h2, P.t__tilde__, P.t ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_66})

V_81 = Vertex(name = 'V_81',
              particles = [ P.Z, P.e__plus__, P.e__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_5, L.L_7 ],
              couplings = {(0,1):C.GC_35,(0,0):C.GC_26})

V_82 = Vertex(name = 'V_82',
              particles = [ P.Z, P.m__plus__, P.m__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_5, L.L_7 ],
              couplings = {(0,1):C.GC_35,(0,0):C.GC_26})

V_83 = Vertex(name = 'V_83',
              particles = [ P.Z, P.tt__plus__, P.tt__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_5, L.L_7 ],
              couplings = {(0,1):C.GC_35,(0,0):C.GC_26})

V_84 = Vertex(name = 'V_84',
              particles = [ P.Zp, P.e__plus__, P.e__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_5, L.L_7 ],
              couplings = {(0,1):C.GC_38,(0,0):C.GC_29})

V_85 = Vertex(name = 'V_85',
              particles = [ P.Zp, P.m__plus__, P.m__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_5, L.L_7 ],
              couplings = {(0,1):C.GC_38,(0,0):C.GC_29})

V_86 = Vertex(name = 'V_86',
              particles = [ P.Zp, P.tt__plus__, P.tt__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_5, L.L_7 ],
              couplings = {(0,1):C.GC_38,(0,0):C.GC_29})

V_87 = Vertex(name = 'V_87',
              particles = [ P.W__plus__, P.ve__tilde__, P.e__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_21})

V_88 = Vertex(name = 'V_88',
              particles = [ P.W__plus__, P.vm__tilde__, P.m__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_21})

V_89 = Vertex(name = 'V_89',
              particles = [ P.W__plus__, P.vt__tilde__, P.tt__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_21})

V_90 = Vertex(name = 'V_90',
              particles = [ P.W__minus__, P.e__plus__, P.ve ],
              color = [ '1' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_21})

V_91 = Vertex(name = 'V_91',
              particles = [ P.W__minus__, P.m__plus__, P.vm ],
              color = [ '1' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_21})

V_92 = Vertex(name = 'V_92',
              particles = [ P.W__minus__, P.tt__plus__, P.vt ],
              color = [ '1' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_21})

V_93 = Vertex(name = 'V_93',
              particles = [ P.Z, P.u__tilde__, P.u ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5, L.L_8 ],
              couplings = {(0,1):C.GC_34,(0,0):C.GC_27})

V_94 = Vertex(name = 'V_94',
              particles = [ P.Z, P.c__tilde__, P.c ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5, L.L_8 ],
              couplings = {(0,1):C.GC_34,(0,0):C.GC_27})

V_95 = Vertex(name = 'V_95',
              particles = [ P.Z, P.t__tilde__, P.t ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5, L.L_8 ],
              couplings = {(0,1):C.GC_34,(0,0):C.GC_27})

V_96 = Vertex(name = 'V_96',
              particles = [ P.Zp, P.u__tilde__, P.u ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5, L.L_8 ],
              couplings = {(0,1):C.GC_37,(0,0):C.GC_28})

V_97 = Vertex(name = 'V_97',
              particles = [ P.Zp, P.c__tilde__, P.c ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5, L.L_8 ],
              couplings = {(0,1):C.GC_37,(0,0):C.GC_28})

V_98 = Vertex(name = 'V_98',
              particles = [ P.Zp, P.t__tilde__, P.t ],
              color = [ 'Identity(2,3)' ],
              lorentz = [ L.L_5, L.L_8 ],
              couplings = {(0,1):C.GC_37,(0,0):C.GC_28})

V_99 = Vertex(name = 'V_99',
              particles = [ P.Z, P.ve__tilde__, P.ve ],
              color = [ '1' ],
              lorentz = [ L.L_5 ],
              couplings = {(0,0):C.GC_36})

V_100 = Vertex(name = 'V_100',
               particles = [ P.Z, P.vm__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.L_5 ],
               couplings = {(0,0):C.GC_36})

V_101 = Vertex(name = 'V_101',
               particles = [ P.Z, P.vt__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.L_5 ],
               couplings = {(0,0):C.GC_36})

V_102 = Vertex(name = 'V_102',
               particles = [ P.Zp, P.ve__tilde__, P.ve ],
               color = [ '1' ],
               lorentz = [ L.L_5 ],
               couplings = {(0,0):C.GC_39})

V_103 = Vertex(name = 'V_103',
               particles = [ P.Zp, P.vm__tilde__, P.vm ],
               color = [ '1' ],
               lorentz = [ L.L_5 ],
               couplings = {(0,0):C.GC_39})

V_104 = Vertex(name = 'V_104',
               particles = [ P.Zp, P.vt__tilde__, P.vt ],
               color = [ '1' ],
               lorentz = [ L.L_5 ],
               couplings = {(0,0):C.GC_39})

