# This file was automatically created by FeynRules $Revision: 364 $
# Mathematica version: 7.0 for Mac OS X x86 (64-bit) (November 11, 2008)
# Date: Wed 10 Nov 2010 10:19:46



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec

cabi = Parameter(name = 'cabi',
                 nature = 'external',
                 type = 'real',
                 value = 0.488,
                 texname = '\\theta _c',
                 lhablock = 'CKMBLOCK',
                 lhacode = [ 1 ])

aXM1 = Parameter(name = 'aXM1',
                 nature = 'external',
                 type = 'real',
                 value = 127.9,
                 texname = '\\text{aXM1}',
                 lhablock = 'HIDDEN',
                 lhacode = [ 1 ])

eta = Parameter(name = 'eta',
                nature = 'external',
                type = 'real',
                value = 0.01,
                texname = '\\eta ',
                lhablock = 'HIDDEN',
                lhacode = [ 2 ])

rho = Parameter(name = 'rho',
                nature = 'external',
                type = 'real',
                value = 0.010142,
                texname = '\\rho ',
                lhablock = 'HIDDEN',
                lhacode = [ 3 ])

kap = Parameter(name = 'kap',
                nature = 'external',
                type = 'real',
                value = 0.0977392,
                texname = '\\text{kap}',
                lhablock = 'HIDDEN',
                lhacode = [ 4 ])

l = Parameter(name = 'l',
              nature = 'external',
              type = 'real',
              value = 0.42568,
              texname = 'l',
              lhablock = 'HIGGS',
              lhacode = [ 1 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.9,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

Gf = Parameter(name = 'Gf',
               nature = 'external',
               type = 'real',
               value = 0.000011663900000000002,
               texname = '\\text{Gf}',
               lhablock = 'SMINPUTS',
               lhacode = [ 2 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.118,
               texname = '\\text{aS}',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

ymc = Parameter(name = 'ymc',
                nature = 'external',
                type = 'real',
                value = 1.42,
                texname = '\\text{ymc}',
                lhablock = 'YUKAWA',
                lhacode = [ 4 ])

ymb = Parameter(name = 'ymb',
                nature = 'external',
                type = 'real',
                value = 4.7,
                texname = '\\text{ymb}',
                lhablock = 'YUKAWA',
                lhacode = [ 5 ])

ymt = Parameter(name = 'ymt',
                nature = 'external',
                type = 'real',
                value = 174.3,
                texname = '\\text{ymt}',
                lhablock = 'YUKAWA',
                lhacode = [ 6 ])

ymtau = Parameter(name = 'ymtau',
                  nature = 'external',
                  type = 'real',
                  value = 1.777,
                  texname = '\\text{ymtau}',
                  lhablock = 'YUKAWA',
                  lhacode = [ 15 ])

MTA = Parameter(name = 'MTA',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{MTA}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MC = Parameter(name = 'MC',
               nature = 'external',
               type = 'real',
               value = 1.42,
               texname = '\\text{MC}',
               lhablock = 'MASS',
               lhacode = [ 4 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 174.3,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MB = Parameter(name = 'MB',
               nature = 'external',
               type = 'real',
               value = 4.7,
               texname = '\\text{MB}',
               lhablock = 'MASS',
               lhacode = [ 5 ])

MZ = Parameter(name = 'MZ',
               nature = 'external',
               type = 'real',
               value = 91.188,
               texname = '\\text{MZ}',
               lhablock = 'MASS',
               lhacode = [ 23 ])

MZp = Parameter(name = 'MZp',
                nature = 'external',
                type = 'real',
                value = 500,
                texname = '\\text{MZp}',
                lhablock = 'MASS',
                lhacode = [ 1023 ])

MW = Parameter(name = 'MW',
               nature = 'external',
               type = 'real',
               value = 80.419,
               texname = '\\text{MW}',
               lhablock = 'MASS',
               lhacode = [ 24 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.50833649,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

WZ = Parameter(name = 'WZ',
               nature = 'external',
               type = 'real',
               value = 2.44140351,
               texname = '\\text{WZ}',
               lhablock = 'DECAY',
               lhacode = [ 23 ])

WZp = Parameter(name = 'WZp',
                nature = 'external',
                type = 'real',
                value = 0.0008252,
                texname = '\\text{WZp}',
                lhablock = 'DECAY',
                lhacode = [ 1023 ])

WW = Parameter(name = 'WW',
               nature = 'external',
               type = 'real',
               value = 2.04759951,
               texname = '\\text{WW}',
               lhablock = 'DECAY',
               lhacode = [ 24 ])

WH1 = Parameter(name = 'WH1',
                nature = 'external',
                type = 'real',
                value = 0.00282299,
                texname = '\\text{WH1}',
                lhablock = 'DECAY',
                lhacode = [ 25 ])

WH2 = Parameter(name = 'WH2',
                nature = 'external',
                type = 'real',
                value = 5.23795,
                texname = '\\text{WH2}',
                lhablock = 'DECAY',
                lhacode = [ 35 ])

cw = Parameter(name = 'cw',
               nature = 'internal',
               type = 'real',
               value = 'MW/MZ',
               texname = 'c_w')

aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1/aEWM1',
                texname = '\\text{aEW}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

aX = Parameter(name = 'aX',
               nature = 'internal',
               type = 'real',
               value = '1/aXM1',
               texname = '\\text{aX}')

MZ0 = Parameter(name = 'MZ0',
                nature = 'internal',
                type = 'real',
                value = 'MZ',
                texname = '\\text{MZ0}')

MX = Parameter(name = 'MX',
               nature = 'internal',
               type = 'real',
               value = 'MZp',
               texname = '\\text{MX}')

v = Parameter(name = 'v',
              nature = 'internal',
              type = 'real',
              value = '1/(2**0.25*cmath.sqrt(Gf))',
              texname = 'v')

chi = Parameter(name = 'chi',
                nature = 'internal',
                type = 'real',
                value = '(-1 + cmath.sqrt(1 + 4*eta**2))/(2.*eta)',
                texname = '\\chi ')

CKM11 = Parameter(name = 'CKM11',
                  nature = 'internal',
                  type = 'complex',
                  value = 'cmath.cos(cabi)',
                  texname = '\\text{CKM11}')

CKM12 = Parameter(name = 'CKM12',
                  nature = 'internal',
                  type = 'complex',
                  value = 'cmath.sin(cabi)',
                  texname = '\\text{CKM12}')

CKM21 = Parameter(name = 'CKM21',
                  nature = 'internal',
                  type = 'complex',
                  value = '-cmath.sin(cabi)',
                  texname = '\\text{CKM21}')

CKM22 = Parameter(name = 'CKM22',
                  nature = 'internal',
                  type = 'complex',
                  value = 'cmath.cos(cabi)',
                  texname = '\\text{CKM22}')

DZ = Parameter(name = 'DZ',
               nature = 'internal',
               type = 'real',
               value = 'MX**2/MZ0**2',
               texname = '\\text{DZ}')

sw = Parameter(name = 'sw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - cw**2)',
               texname = 's_w')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)',
               texname = 'e')

gX = Parameter(name = 'gX',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aX)*cmath.sqrt(cmath.pi)',
               texname = 'g_X')

yb = Parameter(name = 'yb',
               nature = 'internal',
               type = 'real',
               value = '(ymb*cmath.sqrt(2))/v',
               texname = '\\text{yb}')

yc = Parameter(name = 'yc',
               nature = 'internal',
               type = 'real',
               value = '(ymc*cmath.sqrt(2))/v',
               texname = '\\text{yc}')

yt = Parameter(name = 'yt',
               nature = 'internal',
               type = 'real',
               value = '(ymt*cmath.sqrt(2))/v',
               texname = '\\text{yt}')

ytau = Parameter(name = 'ytau',
                 nature = 'internal',
                 type = 'real',
                 value = '(ymtau*cmath.sqrt(2))/v',
                 texname = '\\text{ytau}')

alp = Parameter(name = 'alp',
                nature = 'internal',
                type = 'real',
                value = '-cmath.atan((2*eta*sw)/(1 - DZ - eta**2*sw**2))/2.',
                texname = '\\text{alp}')

g1 = Parameter(name = 'g1',
               nature = 'internal',
               type = 'real',
               value = 'ee/cw',
               texname = 'g_1')

gw = Parameter(name = 'gw',
               nature = 'internal',
               type = 'real',
               value = 'ee/sw',
               texname = 'g_w')

xi = Parameter(name = 'xi',
               nature = 'internal',
               type = 'real',
               value = 'MX/gX',
               texname = '\\xi ')

ca = Parameter(name = 'ca',
               nature = 'internal',
               type = 'real',
               value = 'cmath.cos(alp)',
               texname = 'c_{\\alpha }')

MH1 = Parameter(name = 'MH1',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(l*v**2 + rho*xi**2 - cmath.sqrt(kap**2*v**2*xi**2 + (l*v**2 - rho*xi**2)**2))',
                texname = '\\text{MH1}')

MH2 = Parameter(name = 'MH2',
                nature = 'internal',
                type = 'real',
                value = 'cmath.sqrt(l*v**2 + rho*xi**2 + cmath.sqrt(kap**2*v**2*xi**2 + (l*v**2 - rho*xi**2)**2))',
                texname = '\\text{MH2}')

muH2 = Parameter(name = 'muH2',
                 nature = 'internal',
                 type = 'real',
                 value = '(kap*v**2 + l*xi**2)/2.',
                 texname = '\\text{muH2}')

muSM2 = Parameter(name = 'muSM2',
                  nature = 'internal',
                  type = 'real',
                  value = '(rho*v**2 + kap*xi**2)/2.',
                  texname = '\\text{muSM2}')

sa = Parameter(name = 'sa',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sin(alp)',
               texname = 's_{\\alpha }')

th = Parameter(name = 'th',
               nature = 'internal',
               type = 'real',
               value = 'cmath.atan((kap*v*xi)/(-(l*v**2) + rho*xi**2))/2.',
               texname = '\\text{th}')

ch = Parameter(name = 'ch',
               nature = 'internal',
               type = 'real',
               value = 'cmath.cos(th)',
               texname = 'c_h')

sh = Parameter(name = 'sh',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sin(th)',
               texname = 's_h')

