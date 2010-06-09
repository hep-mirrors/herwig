# This file was automatically created by FeynRules $Revision: 173 $
# Mathematica version: 7.0 for Mac OS X x86 (64-bit) (November 11, 2008)
# Date: Sun 6 Jun 2010 11:38:46


from object_library import all_vertices, Vertex
from particles import Particle as P
from couplings import Coupling as C
from lorentz import Lorentz as L


V_1 = Vertex(particles = [ P.G, P.G, P.G ],
             color = [ 'f(1,2,3)' ],
             lorentz = [ L.L_9 ],
             couplings = {(0,0):C.GC_5})

V_2 = Vertex(particles = [ P.G, P.G, P.G, P.G ],
             color = [ 'f(2,3,a1)*f(a1,1,4)', 'f(2,4,a1)*f(a1,1,3)', 'f(3,4,a1)*f(a1,1,2)' ],
             lorentz = [ L.L_12, L.L_14, L.L_15 ],
             couplings = {(1,1):C.GC_6,(2,0):C.GC_6,(0,2):C.GC_6})

V_3 = Vertex(particles = [ P.A, P.W__plus__, P.W__minus__ ],
             color = [ '1' ],
             lorentz = [ L.L_9 ],
             couplings = {(0,0):C.GC_30})

V_4 = Vertex(particles = [ P.A, P.A, P.W__plus__, P.W__minus__ ],
             color = [ '1' ],
             lorentz = [ L.L_13 ],
             couplings = {(0,0):C.GC_33})

V_5 = Vertex(particles = [ P.W__plus__, P.W__minus__, P.Z ],
             color = [ '1' ],
             lorentz = [ L.L_9 ],
             couplings = {(0,0):C.GC_7})

V_6 = Vertex(particles = [ P.W__plus__, P.W__minus__, P.Zp ],
             color = [ '1' ],
             lorentz = [ L.L_9 ],
             couplings = {(0,0):C.GC_10})

V_7 = Vertex(particles = [ P.W__plus__, P.W__plus__, P.W__minus__, P.W__minus__ ],
             color = [ '1' ],
             lorentz = [ L.L_13 ],
             couplings = {(0,0):C.GC_8})

V_8 = Vertex(particles = [ P.A, P.W__plus__, P.W__minus__, P.Z ],
             color = [ '1' ],
             lorentz = [ L.L_16 ],
             couplings = {(0,0):C.GC_31})

V_9 = Vertex(particles = [ P.W__plus__, P.W__minus__, P.Z, P.Z ],
             color = [ '1' ],
             lorentz = [ L.L_13 ],
             couplings = {(0,0):C.GC_9})

V_10 = Vertex(particles = [ P.A, P.W__plus__, P.W__minus__, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_16 ],
              couplings = {(0,0):C.GC_32})

V_11 = Vertex(particles = [ P.W__plus__, P.W__minus__, P.Z, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_13 ],
              couplings = {(0,0):C.GC_11})

V_12 = Vertex(particles = [ P.W__plus__, P.W__minus__, P.Zp, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_13 ],
              couplings = {(0,0):C.GC_12})

V_13 = Vertex(particles = [ P.h1, P.h1, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_10 ],
              couplings = {(0,0):C.GC_15})

V_14 = Vertex(particles = [ P.h1, P.h1, P.h1, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_10 ],
              couplings = {(0,0):C.GC_14})

V_15 = Vertex(particles = [ P.h1, P.h2, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_10 ],
              couplings = {(0,0):C.GC_13})

V_16 = Vertex(particles = [ P.h1, P.h1, P.h1, P.h1 ],
              color = [ '1' ],
              lorentz = [ L.L_10 ],
              couplings = {(0,0):C.GC_17})

V_17 = Vertex(particles = [ P.h2, P.h2, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_10 ],
              couplings = {(0,0):C.GC_16})

V_18 = Vertex(particles = [ P.h1, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_1 ],
              couplings = {(0,0):C.GC_59})

V_19 = Vertex(particles = [ P.h1, P.h1, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_1 ],
              couplings = {(0,0):C.GC_58})

V_20 = Vertex(particles = [ P.h2, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_1 ],
              couplings = {(0,0):C.GC_57})

V_21 = Vertex(particles = [ P.h1, P.h1, P.h1 ],
              color = [ '1' ],
              lorentz = [ L.L_1 ],
              couplings = {(0,0):C.GC_60})

V_22 = Vertex(particles = [ P.W__plus__, P.W__minus__, P.h1, P.h1 ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_18})

V_23 = Vertex(particles = [ P.W__plus__, P.W__minus__, P.h1, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_19})

V_24 = Vertex(particles = [ P.W__plus__, P.W__minus__, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_20})

V_25 = Vertex(particles = [ P.W__plus__, P.W__minus__, P.h1 ],
              color = [ '1' ],
              lorentz = [ L.L_8 ],
              couplings = {(0,0):C.GC_49})

V_26 = Vertex(particles = [ P.W__plus__, P.W__minus__, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_8 ],
              couplings = {(0,0):C.GC_50})

V_27 = Vertex(particles = [ P.Z, P.Z, P.h1, P.h1 ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_40})

V_28 = Vertex(particles = [ P.Z, P.Z, P.h1, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_43})

V_29 = Vertex(particles = [ P.Z, P.Z, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_46})

V_30 = Vertex(particles = [ P.Z, P.Z, P.h1 ],
              color = [ '1' ],
              lorentz = [ L.L_8 ],
              couplings = {(0,0):C.GC_56})

V_31 = Vertex(particles = [ P.Z, P.Z, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_8 ],
              couplings = {(0,0):C.GC_53})

V_32 = Vertex(particles = [ P.Z, P.Zp, P.h1, P.h1 ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_41})

V_33 = Vertex(particles = [ P.Z, P.Zp, P.h1, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_44})

V_34 = Vertex(particles = [ P.Z, P.Zp, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_47})

V_35 = Vertex(particles = [ P.Z, P.Zp, P.h1 ],
              color = [ '1' ],
              lorentz = [ L.L_8 ],
              couplings = {(0,0):C.GC_55})

V_36 = Vertex(particles = [ P.Z, P.Zp, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_8 ],
              couplings = {(0,0):C.GC_52})

V_37 = Vertex(particles = [ P.Zp, P.Zp, P.h1, P.h1 ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_42})

V_38 = Vertex(particles = [ P.Zp, P.Zp, P.h1, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_45})

V_39 = Vertex(particles = [ P.Zp, P.Zp, P.h2, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_11 ],
              couplings = {(0,0):C.GC_48})

V_40 = Vertex(particles = [ P.Zp, P.Zp, P.h1 ],
              color = [ '1' ],
              lorentz = [ L.L_8 ],
              couplings = {(0,0):C.GC_54})

V_41 = Vertex(particles = [ P.Zp, P.Zp, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_8 ],
              couplings = {(0,0):C.GC_51})

V_42 = Vertex(particles = [ P.d, P.d, P.A ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_1})

V_43 = Vertex(particles = [ P.s, P.s, P.A ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_1})

V_44 = Vertex(particles = [ P.b, P.b, P.A ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_1})

V_45 = Vertex(particles = [ P.e__minus__, P.e__minus__, P.A ],
              color = [ '1' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_3})

V_46 = Vertex(particles = [ P.m__minus__, P.m__minus__, P.A ],
              color = [ '1' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_3})

V_47 = Vertex(particles = [ P.tt__minus__, P.tt__minus__, P.A ],
              color = [ '1' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_3})

V_48 = Vertex(particles = [ P.u, P.u, P.A ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_2})

V_49 = Vertex(particles = [ P.c, P.c, P.A ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_2})

V_50 = Vertex(particles = [ P.t, P.t, P.A ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_2})

V_51 = Vertex(particles = [ P.d, P.d, P.G ],
              color = [ 'T(3,2,1)' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_4})

V_52 = Vertex(particles = [ P.s, P.s, P.G ],
              color = [ 'T(3,2,1)' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_4})

V_53 = Vertex(particles = [ P.b, P.b, P.G ],
              color = [ 'T(3,2,1)' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_4})

V_54 = Vertex(particles = [ P.b, P.b, P.h1 ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_61})

V_55 = Vertex(particles = [ P.b, P.b, P.h2 ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_62})

V_56 = Vertex(particles = [ P.d, P.d, P.Z ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4, L.L_5 ],
              couplings = {(0,1):C.GC_34,(0,0):C.GC_26})

V_57 = Vertex(particles = [ P.s, P.s, P.Z ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4, L.L_5 ],
              couplings = {(0,1):C.GC_34,(0,0):C.GC_26})

V_58 = Vertex(particles = [ P.b, P.b, P.Z ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4, L.L_5 ],
              couplings = {(0,1):C.GC_34,(0,0):C.GC_26})

V_59 = Vertex(particles = [ P.d, P.d, P.Zp ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4, L.L_5 ],
              couplings = {(0,1):C.GC_37,(0,0):C.GC_29})

V_60 = Vertex(particles = [ P.s, P.s, P.Zp ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4, L.L_5 ],
              couplings = {(0,1):C.GC_37,(0,0):C.GC_29})

V_61 = Vertex(particles = [ P.b, P.b, P.Zp ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4, L.L_5 ],
              couplings = {(0,1):C.GC_37,(0,0):C.GC_29})

V_62 = Vertex(particles = [ P.d, P.u, P.W__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_22})

V_63 = Vertex(particles = [ P.d, P.c, P.W__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_24})

V_64 = Vertex(particles = [ P.s, P.u, P.W__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_23})

V_65 = Vertex(particles = [ P.s, P.c, P.W__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_25})

V_66 = Vertex(particles = [ P.b, P.t, P.W__plus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_21})

V_67 = Vertex(particles = [ P.u, P.d, P.W__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_69})

V_68 = Vertex(particles = [ P.c, P.d, P.W__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_71})

V_69 = Vertex(particles = [ P.u, P.s, P.W__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_70})

V_70 = Vertex(particles = [ P.c, P.s, P.W__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_72})

V_71 = Vertex(particles = [ P.t, P.b, P.W__minus__ ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_21})

V_72 = Vertex(particles = [ P.u, P.u, P.G ],
              color = [ 'T(3,2,1)' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_4})

V_73 = Vertex(particles = [ P.c, P.c, P.G ],
              color = [ 'T(3,2,1)' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_4})

V_74 = Vertex(particles = [ P.t, P.t, P.G ],
              color = [ 'T(3,2,1)' ],
              lorentz = [ L.L_3 ],
              couplings = {(0,0):C.GC_4})

V_75 = Vertex(particles = [ P.tt__minus__, P.tt__minus__, P.h1 ],
              color = [ '1' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_67})

V_76 = Vertex(particles = [ P.c, P.c, P.h1 ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_63})

V_77 = Vertex(particles = [ P.t, P.t, P.h1 ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_65})

V_78 = Vertex(particles = [ P.tt__minus__, P.tt__minus__, P.h2 ],
              color = [ '1' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_68})

V_79 = Vertex(particles = [ P.c, P.c, P.h2 ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_64})

V_80 = Vertex(particles = [ P.t, P.t, P.h2 ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_2 ],
              couplings = {(0,0):C.GC_66})

V_81 = Vertex(particles = [ P.e__minus__, P.e__minus__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.L_4, L.L_6 ],
              couplings = {(0,1):C.GC_35,(0,0):C.GC_26})

V_82 = Vertex(particles = [ P.m__minus__, P.m__minus__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.L_4, L.L_6 ],
              couplings = {(0,1):C.GC_35,(0,0):C.GC_26})

V_83 = Vertex(particles = [ P.tt__minus__, P.tt__minus__, P.Z ],
              color = [ '1' ],
              lorentz = [ L.L_4, L.L_6 ],
              couplings = {(0,1):C.GC_35,(0,0):C.GC_26})

V_84 = Vertex(particles = [ P.e__minus__, P.e__minus__, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_4, L.L_6 ],
              couplings = {(0,1):C.GC_38,(0,0):C.GC_29})

V_85 = Vertex(particles = [ P.m__minus__, P.m__minus__, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_4, L.L_6 ],
              couplings = {(0,1):C.GC_38,(0,0):C.GC_29})

V_86 = Vertex(particles = [ P.tt__minus__, P.tt__minus__, P.Zp ],
              color = [ '1' ],
              lorentz = [ L.L_4, L.L_6 ],
              couplings = {(0,1):C.GC_38,(0,0):C.GC_29})

V_87 = Vertex(particles = [ P.e__minus__, P.ve, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_21})

V_88 = Vertex(particles = [ P.m__minus__, P.vm, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_21})

V_89 = Vertex(particles = [ P.tt__minus__, P.vt, P.W__plus__ ],
              color = [ '1' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_21})

V_90 = Vertex(particles = [ P.ve, P.e__minus__, P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_21})

V_91 = Vertex(particles = [ P.vm, P.m__minus__, P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_21})

V_92 = Vertex(particles = [ P.vt, P.tt__minus__, P.W__minus__ ],
              color = [ '1' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_21})

V_93 = Vertex(particles = [ P.u, P.u, P.Z ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4, L.L_7 ],
              couplings = {(0,1):C.GC_34,(0,0):C.GC_27})

V_94 = Vertex(particles = [ P.c, P.c, P.Z ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4, L.L_7 ],
              couplings = {(0,1):C.GC_34,(0,0):C.GC_27})

V_95 = Vertex(particles = [ P.t, P.t, P.Z ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4, L.L_7 ],
              couplings = {(0,1):C.GC_34,(0,0):C.GC_27})

V_96 = Vertex(particles = [ P.u, P.u, P.Zp ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4, L.L_7 ],
              couplings = {(0,1):C.GC_37,(0,0):C.GC_28})

V_97 = Vertex(particles = [ P.c, P.c, P.Zp ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4, L.L_7 ],
              couplings = {(0,1):C.GC_37,(0,0):C.GC_28})

V_98 = Vertex(particles = [ P.t, P.t, P.Zp ],
              color = [ 'Identity(1,2)' ],
              lorentz = [ L.L_4, L.L_7 ],
              couplings = {(0,1):C.GC_37,(0,0):C.GC_28})

V_99 = Vertex(particles = [ P.ve, P.ve, P.Z ],
              color = [ '1' ],
              lorentz = [ L.L_4 ],
              couplings = {(0,0):C.GC_36})

V_100 = Vertex(particles = [ P.vm, P.vm, P.Z ],
               color = [ '1' ],
               lorentz = [ L.L_4 ],
               couplings = {(0,0):C.GC_36})

V_101 = Vertex(particles = [ P.vt, P.vt, P.Z ],
               color = [ '1' ],
               lorentz = [ L.L_4 ],
               couplings = {(0,0):C.GC_36})

V_102 = Vertex(particles = [ P.ve, P.ve, P.Zp ],
               color = [ '1' ],
               lorentz = [ L.L_4 ],
               couplings = {(0,0):C.GC_39})

V_103 = Vertex(particles = [ P.vm, P.vm, P.Zp ],
               color = [ '1' ],
               lorentz = [ L.L_4 ],
               couplings = {(0,0):C.GC_39})

V_104 = Vertex(particles = [ P.vt, P.vt, P.Zp ],
               color = [ '1' ],
               lorentz = [ L.L_4 ],
               couplings = {(0,0):C.GC_39})

