from __future__ import print_function
import copy
from .helpers import SkipThisVertex,def_from_model
from .converter import py2cpp
import string,re
from string import Template

epsValue=[[[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
           [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]],
          [[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
           [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]],
          [[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
           [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]],
          [[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],
           [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]],[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]]]

epsValue[0][1][2][3] = -1.
epsValue[0][1][3][2] =  1.
epsValue[0][2][1][3] =  1.
epsValue[0][2][3][1] = -1.
epsValue[0][3][1][2] = -1.
epsValue[0][3][2][1] =  1.
epsValue[1][0][2][3] =  1.
epsValue[1][0][3][2] = -1.
epsValue[1][2][0][3] = -1.
epsValue[1][2][3][0] =  1.
epsValue[1][3][0][2] =  1.
epsValue[1][3][2][0] = -1.
epsValue[2][0][1][3] = -1.
epsValue[2][0][3][1] =  1.
epsValue[2][1][0][3] =  1.
epsValue[2][1][3][0] = -1.
epsValue[2][3][0][1] = -1.
epsValue[2][3][1][0] =  1.
epsValue[3][0][1][2] =  1.
epsValue[3][0][2][1] = -1.
epsValue[3][1][0][2] = -1.
epsValue[3][1][2][0] =  1.
epsValue[3][2][0][1] =  1.
epsValue[3][2][1][0] = -1.

# self contracted tensor propagator
tPropA=[[],[],[],[]]
tPropA[0].append(Template("-2. / 3. * (M${iloc}2 + 2 * P${iloc}t ** 2) * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[0].append(Template("-4. / 3. * P${iloc}t * P${iloc}x * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[0].append(Template("-4. / 3. * P${iloc}t * P${iloc}y * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[0].append(Template("-4. / 3. * P${iloc}t * P${iloc}z * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[1].append(Template("-4. / 3. * P${iloc}t * P${iloc}x * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[1].append(Template(" 2. / 3. * (M${iloc}2 - 2 * P${iloc}x ** 2) * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[1].append(Template("-4. / 3. * P${iloc}x * P${iloc}y * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[1].append(Template("-4. / 3. * P${iloc}x * P${iloc}z * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[2].append(Template("-4. / 3. * P${iloc}t * P${iloc}y * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[2].append(Template("-4. / 3. * P${iloc}x * P${iloc}y * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[2].append(Template(" 2. / 3. * (M${iloc}2 - 2 * P${iloc}y ** 2) * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[2].append(Template("-4. / 3. * P${iloc}y * P${iloc}z * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[3].append(Template("-4. / 3. * P${iloc}t * P${iloc}z * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[3].append(Template("-4. / 3. * P${iloc}x * P${iloc}z * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[3].append(Template("-4. / 3. * P${iloc}y * P${iloc}z * (M${iloc}2 -p2)*OM${iloc}**2"))
tPropA[3].append(Template(" 2. / 3. * (M${iloc}2 - 2 * P${iloc}z ** 2) * (M${iloc}2 -p2)*OM${iloc}**2"))

# tensor propagator 1 index contracted
tPropB=[[[],[],[],[]],[[],[],[],[]],[[],[],[],[]],[[],[],[],[]]]
tPropB[0][0].append(Template("4. / 3. * (${V}t - ${dot}*OM${iloc} * P${iloc}t) * (1. - OM${iloc} * P${iloc}t ** 2)"))
tPropB[0][0].append(Template("-2 * (${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}t * P${iloc}x - 2. / 3. * (1. - OM${iloc} * P${iloc}t ** 2) * (${V}x - ${dot}*OM${iloc} * P${iloc}x)"))
tPropB[0][0].append(Template(" -2 * (${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}t * P${iloc}y - 2. / 3. * (1. - OM${iloc} * P${iloc}t ** 2) * (${V}y - ${dot}*OM${iloc} * P${iloc}y)"))
tPropB[0][0].append(Template(" -2 * (${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}t * P${iloc}z - 2. / 3. * (1. - OM${iloc} * P${iloc}t ** 2) * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[0][1].append(Template(" -(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}t * P${iloc}x / 3. + (1. - OM${iloc} * P${iloc}t ** 2) * (${V}x - ${dot}*OM${iloc} * P${iloc}x)"))
tPropB[0][1].append(Template(" (${V}t - ${dot}*OM${iloc} * P${iloc}t) * (-1. - OM${iloc} * P${iloc}x ** 2) - OM${iloc} * P${iloc}t * P${iloc}x * (${V}x - ${dot}*OM${iloc} * P${iloc}x) / 3."))
tPropB[0][1].append(Template("-(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}x * P${iloc}y - OM${iloc} * P${iloc}t * P${iloc}y * (${V}x - ${dot}*OM${iloc} * P${iloc}x) + 2. / 3.*OM${iloc} * P${iloc}t * P${iloc}x * (${V}y - ${dot}*OM${iloc} * P${iloc}y)"))
tPropB[0][1].append(Template("-(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}x * P${iloc}z - OM${iloc} * P${iloc}t * P${iloc}z * (${V}x - ${dot}*OM${iloc} * P${iloc}x) + 2. / 3.*OM${iloc} * P${iloc}t * P${iloc}x * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[0][2].append(Template(" -(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}t * P${iloc}y / 3. + (1. - OM${iloc} * P${iloc}t ** 2) * (${V}y - ${dot}*OM${iloc} * P${iloc}y)"))
tPropB[0][2].append(Template("-(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}x * P${iloc}y - OM${iloc} * P${iloc}t * P${iloc}x * (${V}y - ${dot}*OM${iloc} * P${iloc}y) + 2. / 3.*OM${iloc} * P${iloc}t * P${iloc}y * (${V}x - ${dot}*OM${iloc} * P${iloc}x)"))
tPropB[0][2].append(Template(" (${V}t - ${dot}*OM${iloc} * P${iloc}t) * (-1. - OM${iloc} * P${iloc}y ** 2) - OM${iloc} * P${iloc}t * P${iloc}y * (${V}y - ${dot}*OM${iloc} * P${iloc}y) / 3."))
tPropB[0][2].append(Template("-(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}y * P${iloc}z - OM${iloc} * P${iloc}t * P${iloc}z * (${V}y - ${dot}*OM${iloc} * P${iloc}y) + 2. / 3.*OM${iloc} * P${iloc}t * P${iloc}y * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[0][3].append(Template(" -(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}t * P${iloc}z / 3. + (1. - OM${iloc} * P${iloc}t ** 2) * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[0][3].append(Template("-(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}x * P${iloc}z - OM${iloc} * P${iloc}t * P${iloc}x * (${V}z - ${dot}*OM${iloc} * P${iloc}z) + 2. / 3.*OM${iloc} * P${iloc}t * P${iloc}z * (${V}x - ${dot}*OM${iloc} * P${iloc}x)"))
tPropB[0][3].append(Template(" -(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}y * P${iloc}z - OM${iloc} * P${iloc}t * P${iloc}y * (${V}z - ${dot}*OM${iloc} * P${iloc}z) + 2. / 3.*OM${iloc} * P${iloc}t * P${iloc}z * (${V}y - ${dot}*OM${iloc} * P${iloc}y)"))
tPropB[0][3].append(Template("(${V}t - ${dot}*OM${iloc} * P${iloc}t) * (-1. - OM${iloc} * P${iloc}z ** 2) - OM${iloc} * P${iloc}t * P${iloc}z * (${V}z - ${dot}*OM${iloc} * P${iloc}z) / 3."))

tPropB[1][0].append(Template("-(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}t * P${iloc}x / 3. + (1 - OM${iloc} * P${iloc}t ** 2) * (${V}x - ${dot}*OM${iloc} * P${iloc}x)"))
tPropB[1][0].append(Template(" (${V}t - ${dot}*OM${iloc} * P${iloc}t) * (-1 - OM${iloc} * P${iloc}x ** 2) - OM${iloc} * P${iloc}t * P${iloc}x * (${V}x - ${dot}*OM${iloc} * P${iloc}x) / 3."))
tPropB[1][0].append(Template(" -(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}x * P${iloc}y - OM${iloc} * P${iloc}t * P${iloc}y * (${V}x - ${dot}*OM${iloc} * P${iloc}x) + 2. / 3.*OM${iloc} * P${iloc}t * P${iloc}x * (${V}y - ${dot}*OM${iloc} * P${iloc}y)"))
tPropB[1][0].append(Template(" -(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}x * P${iloc}z - OM${iloc} * P${iloc}t * P${iloc}z * (${V}x - ${dot}*OM${iloc} * P${iloc}x) + 2. / 3.*OM${iloc} * P${iloc}t * P${iloc}x * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[1][1].append(Template(" -2*OM${iloc} * P${iloc}t * P${iloc}x * (${V}x - ${dot}*OM${iloc} * P${iloc}x) - 2. / 3. * (${V}t - ${dot}*OM${iloc} * P${iloc}t) * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropB[1][1].append(Template(" 4. / 3. * (${V}x - ${dot}*OM${iloc} * P${iloc}x) * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropB[1][1].append(Template(" -2 * (${V}x - ${dot}*OM${iloc} * P${iloc}x)*OM${iloc} * P${iloc}x * P${iloc}y - 2. / 3. * (-1 - OM${iloc} * P${iloc}x ** 2) * (${V}y - ${dot}*OM${iloc} * P${iloc}y)"))
tPropB[1][1].append(Template(" -2 * (${V}x - ${dot}*OM${iloc} * P${iloc}x)*OM${iloc} * P${iloc}x * P${iloc}z - 2. / 3. * (-1 - OM${iloc} * P${iloc}x ** 2) * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[1][2].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}y * (${V}x - ${dot}*OM${iloc} * P${iloc}x) - OM${iloc} * P${iloc}t * P${iloc}x * (${V}y - ${dot}*OM${iloc} * P${iloc}y) + 2. / 3. * (${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}x * P${iloc}y"))
tPropB[1][2].append(Template(" -(${V}x - ${dot}*OM${iloc} * P${iloc}x)*OM${iloc} * P${iloc}x * P${iloc}y / 3. + (-1 - OM${iloc} * P${iloc}x ** 2) * (${V}y - ${dot}*OM${iloc} * P${iloc}y)"))
tPropB[1][2].append(Template(" (${V}x - ${dot}*OM${iloc} * P${iloc}x) * (-1 - OM${iloc} * P${iloc}y ** 2) - OM${iloc} * P${iloc}x * P${iloc}y * (${V}y - ${dot}*OM${iloc} * P${iloc}y) / 3."))
tPropB[1][2].append(Template("-(${V}x - ${dot}*OM${iloc} * P${iloc}x)*OM${iloc} * P${iloc}y * P${iloc}z - OM${iloc} * P${iloc}x * P${iloc}z * (${V}y - ${dot}*OM${iloc} * P${iloc}y) + 2. / 3.*OM${iloc} * P${iloc}x * P${iloc}y * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[1][3].append(Template("-OM${iloc} * P${iloc}t * P${iloc}z * (${V}x - ${dot}*OM${iloc} * P${iloc}x) - OM${iloc} * P${iloc}t * P${iloc}x * (${V}z - ${dot}*OM${iloc} * P${iloc}z) + 2. / 3. * (${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}x * P${iloc}z"))
tPropB[1][3].append(Template(" -(${V}x - ${dot}*OM${iloc} * P${iloc}x)*OM${iloc} * P${iloc}x * P${iloc}z / 3. + (-1 - OM${iloc} * P${iloc}x ** 2) * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[1][3].append(Template(" -(${V}x - ${dot}*OM${iloc} * P${iloc}x)*OM${iloc} * P${iloc}y * P${iloc}z - OM${iloc} * P${iloc}x * P${iloc}y * (${V}z - ${dot}*OM${iloc} * P${iloc}z) + 2. / 3.*OM${iloc} * P${iloc}x * P${iloc}z * (${V}y - ${dot}*OM${iloc} * P${iloc}y)"))
tPropB[1][3].append(Template("(${V}x - ${dot}*OM${iloc} * P${iloc}x) * (-1 - OM${iloc} * P${iloc}z ** 2) - OM${iloc} * P${iloc}x * P${iloc}z * (${V}z - ${dot}*OM${iloc} * P${iloc}z) / 3."))

tPropB[2][0].append(Template(" -(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}t * P${iloc}y / 3. + (1 - OM${iloc} * P${iloc}t ** 2) * (${V}y - ${dot}*OM${iloc} * P${iloc}y)"))
tPropB[2][0].append(Template(" -(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}x * P${iloc}y - OM${iloc} * P${iloc}t * P${iloc}x * (${V}y - ${dot}*OM${iloc} * P${iloc}y) + 2. / 3.*OM${iloc} * P${iloc}t * P${iloc}y * (${V}x - ${dot}*OM${iloc} * P${iloc}x)"))
tPropB[2][0].append(Template(" (${V}t - ${dot}*OM${iloc} * P${iloc}t) * (-1 - OM${iloc} * P${iloc}y ** 2) - OM${iloc} * P${iloc}t * P${iloc}y * (${V}y - ${dot}*OM${iloc} * P${iloc}y) / 3."))
tPropB[2][0].append(Template(" -(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}y * P${iloc}z - OM${iloc} * P${iloc}t * P${iloc}z * (${V}y - ${dot}*OM${iloc} * P${iloc}y) + 2. / 3.*OM${iloc} * P${iloc}t * P${iloc}y * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[2][1].append(Template("-OM${iloc} * P${iloc}t * P${iloc}y * (${V}x - ${dot}*OM${iloc} * P${iloc}x) - OM${iloc} * P${iloc}t * P${iloc}x * (${V}y - ${dot}*OM${iloc} * P${iloc}y) + 2. / 3. * (${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}x * P${iloc}y"))
tPropB[2][1].append(Template("-(${V}x - ${dot}*OM${iloc} * P${iloc}x)*OM${iloc} * P${iloc}x * P${iloc}y / 3. + (-1 - OM${iloc} * P${iloc}x ** 2) * (${V}y - ${dot}*OM${iloc} * P${iloc}y)"))
tPropB[2][1].append(Template(" (${V}x - ${dot}*OM${iloc} * P${iloc}x) * (-1 - OM${iloc} * P${iloc}y ** 2) - OM${iloc} * P${iloc}x * P${iloc}y * (${V}y - ${dot}*OM${iloc} * P${iloc}y) / 3."))
tPropB[2][1].append(Template("-(${V}x - ${dot}*OM${iloc} * P${iloc}x)*OM${iloc} * P${iloc}y * P${iloc}z - OM${iloc} * P${iloc}x * P${iloc}z * (${V}y - ${dot}*OM${iloc} * P${iloc}y) + 2. / 3.*OM${iloc} * P${iloc}x * P${iloc}y * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[2][2].append(Template(" -2*OM${iloc} * P${iloc}t * P${iloc}y * (${V}y - ${dot}*OM${iloc} * P${iloc}y) - 2. / 3. * (${V}t - ${dot}*OM${iloc} * P${iloc}t) * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropB[2][2].append(Template(" -2*OM${iloc} * P${iloc}x * P${iloc}y * (${V}y - ${dot}*OM${iloc} * P${iloc}y) - 2. / 3. * (${V}x - ${dot}*OM${iloc} * P${iloc}x) * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropB[2][2].append(Template("4. / 3. * (${V}y - ${dot}*OM${iloc} * P${iloc}y) * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropB[2][2].append(Template(" -2 * (${V}y - ${dot}*OM${iloc} * P${iloc}y)*OM${iloc} * P${iloc}y * P${iloc}z - 2. / 3. * (-1 - OM${iloc} * P${iloc}y ** 2) * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[2][3].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}z * (${V}y - ${dot}*OM${iloc} * P${iloc}y) - OM${iloc} * P${iloc}t * P${iloc}y * (${V}z - ${dot}*OM${iloc} * P${iloc}z) + 2. / 3. * (${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}y * P${iloc}z"))
tPropB[2][3].append(Template(" -OM${iloc} * P${iloc}x * P${iloc}z * (${V}y - ${dot}*OM${iloc} * P${iloc}y) - OM${iloc} * P${iloc}x * P${iloc}y * (${V}z - ${dot}*OM${iloc} * P${iloc}z) + 2. / 3. * (${V}x - ${dot}*OM${iloc} * P${iloc}x)*OM${iloc} * P${iloc}y * P${iloc}z"))
tPropB[2][3].append(Template(" -(${V}y - ${dot}*OM${iloc} * P${iloc}y)*OM${iloc} * P${iloc}y * P${iloc}z / 3. + (-1 - OM${iloc} * P${iloc}y ** 2) * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[2][3].append(Template(" (${V}y - ${dot}*OM${iloc} * P${iloc}y) * (-1 - OM${iloc} * P${iloc}z ** 2) - OM${iloc} * P${iloc}y * P${iloc}z * (${V}z - ${dot}*OM${iloc} * P${iloc}z) / 3."))

tPropB[3][0].append(Template(" -(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}t * P${iloc}z / 3. + (1 - OM${iloc} * P${iloc}t ** 2) * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[3][0].append(Template(" -(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}x * P${iloc}z - OM${iloc} * P${iloc}t * P${iloc}x * (${V}z - ${dot}*OM${iloc} * P${iloc}z) + 2. / 3.*OM${iloc} * P${iloc}t * P${iloc}z * (${V}x - ${dot}*OM${iloc} * P${iloc}x)"))
tPropB[3][0].append(Template("-(${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}y * P${iloc}z - OM${iloc} * P${iloc}t * P${iloc}y * (${V}z - ${dot}*OM${iloc} * P${iloc}z) + 2. / 3.*OM${iloc} * P${iloc}t * P${iloc}z * (${V}y - ${dot}*OM${iloc} * P${iloc}y)"))
tPropB[3][0].append(Template(" (${V}t - ${dot}*OM${iloc} * P${iloc}t) * (-1 - OM${iloc} * P${iloc}z ** 2) - OM${iloc} * P${iloc}t * P${iloc}z * (${V}z - ${dot}*OM${iloc} * P${iloc}z) / 3."))
tPropB[3][1].append(Template("-OM${iloc} * P${iloc}t * P${iloc}z * (${V}x - ${dot}*OM${iloc} * P${iloc}x) - OM${iloc} * P${iloc}t * P${iloc}x * (${V}z - ${dot}*OM${iloc} * P${iloc}z) + 2. / 3. * (${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}x * P${iloc}z"))
tPropB[3][1].append(Template(" -(${V}x - ${dot}*OM${iloc} * P${iloc}x)*OM${iloc} * P${iloc}x * P${iloc}z / 3. + (-1 - OM${iloc} * P${iloc}x ** 2) * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[3][1].append(Template(" -(${V}x - ${dot}*OM${iloc} * P${iloc}x)*OM${iloc} * P${iloc}y * P${iloc}z - OM${iloc} * P${iloc}x * P${iloc}y * (${V}z - ${dot}*OM${iloc} * P${iloc}z) + 2. / 3.*OM${iloc} * P${iloc}x * P${iloc}z * (${V}y - ${dot}*OM${iloc} * P${iloc}y)"))
tPropB[3][1].append(Template(" (${V}x - ${dot}*OM${iloc} * P${iloc}x) * (-1 - OM${iloc} * P${iloc}z ** 2) - OM${iloc} * P${iloc}x * P${iloc}z * (${V}z - ${dot}*OM${iloc} * P${iloc}z) / 3."))
tPropB[3][2].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}z * (${V}y - ${dot}*OM${iloc} * P${iloc}y) - OM${iloc} * P${iloc}t * P${iloc}y * (${V}z - ${dot}*OM${iloc} * P${iloc}z) + 2. / 3. * (${V}t - ${dot}*OM${iloc} * P${iloc}t)*OM${iloc} * P${iloc}y * P${iloc}z"))
tPropB[3][2].append(Template("-OM${iloc} * P${iloc}x * P${iloc}z * (${V}y - ${dot}*OM${iloc} * P${iloc}y) - OM${iloc} * P${iloc}x * P${iloc}y * (${V}z - ${dot}*OM${iloc} * P${iloc}z) + 2. / 3. * (${V}x - ${dot}*OM${iloc} * P${iloc}x)*OM${iloc} * P${iloc}y * P${iloc}z"))
tPropB[3][2].append(Template(" -(${V}y - ${dot}*OM${iloc} * P${iloc}y)*OM${iloc} * P${iloc}y * P${iloc}z / 3. + (-1 - OM${iloc} * P${iloc}y ** 2) * (${V}z - ${dot}*OM${iloc} * P${iloc}z)"))
tPropB[3][2].append(Template(" (${V}y - ${dot}*OM${iloc} * P${iloc}y) * (-1 - OM${iloc} * P${iloc}z ** 2) - OM${iloc} * P${iloc}y * P${iloc}z * (${V}z - ${dot}*OM${iloc} * P${iloc}z) / 3."))
tPropB[3][3].append(Template(" -2*OM${iloc} * P${iloc}t * P${iloc}z * (${V}z - ${dot}*OM${iloc} * P${iloc}z) - 2. / 3. * (${V}t - ${dot}*OM${iloc} * P${iloc}t) * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropB[3][3].append(Template("-2*OM${iloc} * P${iloc}x * P${iloc}z * (${V}z - ${dot}*OM${iloc} * P${iloc}z) - 2. / 3. * (${V}x - ${dot}*OM${iloc} * P${iloc}x) * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropB[3][3].append(Template("-2*OM${iloc} * P${iloc}y * P${iloc}z * (${V}z - ${dot}*OM${iloc} * P${iloc}z) - 2. / 3. * (${V}y - ${dot}*OM${iloc} * P${iloc}y) * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropB[3][3].append(Template("4. / 3. * (${V}z - ${dot}*OM${iloc} * P${iloc}z) * (-1 - OM${iloc} * P${iloc}z ** 2)"))

# tensor propagator, no contracted indices
tPropC=[[[[],[],[],[]],[[],[],[],[]],[[],[],[],[]],[[],[],[],[]]],
        [[[],[],[],[]],[[],[],[],[]],[[],[],[],[]],[[],[],[],[]]],
        [[[],[],[],[]],[[],[],[],[]],[[],[],[],[]],[[],[],[],[]]],
        [[[],[],[],[]],[[],[],[],[]],[[],[],[],[]],[[],[],[],[]]]]
tPropC[0][0][0].append(Template("4./3. * (1 - OM${iloc} * P${iloc}t ** 2) ** 2"))
tPropC[0][0][0].append(Template("-4./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}t * P${iloc}x"))
tPropC[0][0][0].append(Template("-4./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}t * P${iloc}y"))
tPropC[0][0][0].append(Template("-4./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}t * P${iloc}z"))
tPropC[0][0][1].append(Template("-4./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}t * P${iloc}x"))
tPropC[0][0][1].append(Template("2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x ** 2 - 2./3. * (1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[0][0][1].append(Template("2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}y + 2./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}y"))
tPropC[0][0][1].append(Template("2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}z + 2./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}z"))
tPropC[0][0][2].append(Template("-4./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}t * P${iloc}y"))
tPropC[0][0][2].append(Template("2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}y + 2./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}y"))
tPropC[0][0][2].append(Template("2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y ** 2 - 2./3. * (1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[0][0][2].append(Template("2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y * P${iloc}z + 2./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}y * P${iloc}z"))
tPropC[0][0][3].append(Template("-4./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}t * P${iloc}z"))
tPropC[0][0][3].append(Template("2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}z + 2./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}z"))
tPropC[0][0][3].append(Template("2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y * P${iloc}z + 2./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}y * P${iloc}z"))
tPropC[0][0][3].append(Template("2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}z ** 2 - 2./3. * (1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2)"))

tPropC[0][1][0].append(Template("-4./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}t * P${iloc}x"))
tPropC[0][1][0].append(Template("(1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}x ** 2) + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x ** 2 /3."))
tPropC[0][1][0].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}y + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}y /3."))
tPropC[0][1][0].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}z /3."))
tPropC[0][1][1].append(Template("(1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}x ** 2) + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x ** 2 /3."))
tPropC[0][1][1].append(Template("-4./3.*OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[0][1][1].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}y /3. - OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[0][1][1].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[0][1][2].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}y + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}y /3."))
tPropC[0][1][2].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}y /3. - OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[0][1][2].append(Template("2*OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}x + 2./3.*OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[0][1][2].append(Template("4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z"))
tPropC[0][1][3].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}z /3."))
tPropC[0][1][3].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[0][1][3].append(Template("4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z"))
tPropC[0][1][3].append(Template("2*OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}x + 2./3.*OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}z ** 2)"))

tPropC[0][2][0].append(Template("-4./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}t * P${iloc}y"))
tPropC[0][2][0].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}y + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}y /3."))
tPropC[0][2][0].append(Template("(1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y ** 2 /3."))
tPropC[0][2][0].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y * P${iloc}z /3."))
tPropC[0][2][1].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}y + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}y /3."))
tPropC[0][2][1].append(Template("2*OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}y + 2./3.*OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[0][2][1].append(Template("-OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}x /3."))
tPropC[0][2][1].append(Template("4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z"))
tPropC[0][2][2].append(Template("(1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y ** 2 /3."))
tPropC[0][2][2].append(Template("-OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}x /3."))
tPropC[0][2][2].append(Template("-4./3.*OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[0][2][2].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[0][2][3].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y * P${iloc}z /3."))
tPropC[0][2][3].append(Template("4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z"))
tPropC[0][2][3].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[0][2][3].append(Template("2*OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}y + 2./3.*OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2)"))

tPropC[0][3][0].append(Template("-4./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}t * P${iloc}z"))
tPropC[0][3][0].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}z /3."))
tPropC[0][3][0].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y * P${iloc}z /3."))
tPropC[0][3][0].append(Template("(1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}z ** 2 /3."))
tPropC[0][3][1].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}z /3."))
tPropC[0][3][1].append(Template("2*OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}z + 2./3.*OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[0][3][1].append(Template("4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z"))
tPropC[0][3][1].append(Template("-OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}x /3."))
tPropC[0][3][2].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y * P${iloc}z /3."))
tPropC[0][3][2].append(Template("4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z"))
tPropC[0][3][2].append(Template("2*OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}z + 2./3.*OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[0][3][2].append(Template("-OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[0][3][3].append(Template("(1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}z ** 2 /3."))
tPropC[0][3][3].append(Template("-OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}x /3."))
tPropC[0][3][3].append(Template("-OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[0][3][3].append(Template("-4./3.*OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}z ** 2)"))

tPropC[1][0][0].append(Template(" -4./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}t * P${iloc}x"))
tPropC[1][0][0].append(Template(" (1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}x ** 2) + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x ** 2 /3."))
tPropC[1][0][0].append(Template(" -(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}y + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}y /3."))
tPropC[1][0][0].append(Template(" -(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}z /3."))
tPropC[1][0][1].append(Template(" (1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}x ** 2) + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x ** 2 /3."))
tPropC[1][0][1].append(Template(" -4./3.*OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][0][1].append(Template(" OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}y /3. - OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][0][1].append(Template(" OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][0][2].append(Template(" -(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}y + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}y /3."))
tPropC[1][0][2].append(Template(" OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}y /3. - OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][0][2].append(Template(" 2*OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}x + 2./3.*OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[1][0][2].append(Template(" 4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z"))
tPropC[1][0][3].append(Template(" -(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}z /3."))
tPropC[1][0][3].append(Template(" OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][0][3].append(Template(" 4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z"))
tPropC[1][0][3].append(Template(" 2*OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}x + 2./3.*OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[1][1][0].append(Template(" 2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x ** 2 - 2./3. * (1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][1][0].append(Template(" -4./3.*OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][1][0].append(Template(" 2*OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}y + 2./3.*OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][1][0].append(Template(" 2*OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}z + 2./3.*OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][1][1].append(Template(" -4./3.*OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][1][1].append(Template(" 4./3. * (-1 - OM${iloc} * P${iloc}x ** 2) ** 2"))
tPropC[1][1][1].append(Template(" -4./3. * (-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}x * P${iloc}y"))
tPropC[1][1][1].append(Template(" -4./3. * (-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}x * P${iloc}z"))
tPropC[1][1][2].append(Template(" 2*OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}y + 2./3.*OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][1][2].append(Template(" -4./3. * (-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}x * P${iloc}y"))
tPropC[1][1][2].append(Template(" 2*OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y ** 2 - 2./3. * (-1 - OM${iloc} * P${iloc}x ** 2) * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[1][1][2].append(Template(" 2*OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y * P${iloc}z + 2./3. * (-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}y * P${iloc}z"))
tPropC[1][1][3].append(Template(" 2*OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}z + 2./3.*OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][1][3].append(Template(" -4./3. * (-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}x * P${iloc}z"))
tPropC[1][1][3].append(Template(" 2*OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y * P${iloc}z + 2./3. * (-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}y * P${iloc}z"))
tPropC[1][1][3].append(Template(" 2*OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}z ** 2 - 2./3. * (-1 - OM${iloc} * P${iloc}x ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2)"))

tPropC[1][2][0].append(Template("2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}y + 2./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}y"))
tPropC[1][2][0].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}y /3. - OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][2][0].append(Template("-OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}x /3."))
tPropC[1][2][0].append(Template("4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[1][2][1].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}y /3. - OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][2][1].append(Template(" -4./3. * (-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}x * P${iloc}y							"))
tPropC[1][2][1].append(Template(" (-1 - OM${iloc} * P${iloc}x ** 2) * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y ** 2 /3.	"))
tPropC[1][2][1].append(Template(" -(-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y * P${iloc}z /3."))
tPropC[1][2][2].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}x /3."))
tPropC[1][2][2].append(Template(" (-1 - OM${iloc} * P${iloc}x ** 2) * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y ** 2 /3.	"))
tPropC[1][2][2].append(Template(" -4./3.*OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}y ** 2)							"))
tPropC[1][2][2].append(Template(" OM${iloc}**2 * P${iloc}x * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[1][2][3].append(Template(" 4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[1][2][3].append(Template(" -(-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y * P${iloc}z /3."))
tPropC[1][2][3].append(Template(" OM${iloc}**2 * P${iloc}x * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[1][2][3].append(Template(" 2*OM${iloc}**2 * P${iloc}x * P${iloc}z ** 2 * P${iloc}y + 2./3.*OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2)"))

tPropC[1][3][0].append(Template("2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}z + 2./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}z"))
tPropC[1][3][0].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][3][0].append(Template("4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[1][3][0].append(Template("-OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}x /3."))
tPropC[1][3][1].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[1][3][1].append(Template("-4./3. * (-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}x * P${iloc}z							"))
tPropC[1][3][1].append(Template("-(-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y * P${iloc}z /3."))
tPropC[1][3][1].append(Template("(-1 - OM${iloc} * P${iloc}x ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}z ** 2 /3.	"))
tPropC[1][3][2].append(Template("4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[1][3][2].append(Template("-(-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y * P${iloc}z /3."))
tPropC[1][3][2].append(Template("2*OM${iloc}**2 * P${iloc}x * P${iloc}y ** 2 * P${iloc}z + 2./3.*OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[1][3][2].append(Template("-OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}x * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[1][3][3].append(Template("-OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}x /3."))
tPropC[1][3][3].append(Template("(-1 - OM${iloc} * P${iloc}x ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}z ** 2 /3.	"))
tPropC[1][3][3].append(Template("-OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}x * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[1][3][3].append(Template("-4./3.*OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}z ** 2)"))

tPropC[2][0][0].append(Template("-4./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}t * P${iloc}y							"))
tPropC[2][0][0].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}y + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}y /3.	"))
tPropC[2][0][0].append(Template("(1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y ** 2 /3.	"))
tPropC[2][0][0].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y * P${iloc}z /3.	"))
tPropC[2][0][1].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}y + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}y /3.	"))
tPropC[2][0][1].append(Template("2*OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}y + 2./3.*OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[2][0][1].append(Template("-OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}x /3."))
tPropC[2][0][1].append(Template("4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[2][0][2].append(Template("(1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y ** 2 /3.	"))
tPropC[2][0][2].append(Template("-OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}x /3."))
tPropC[2][0][2].append(Template("-4./3.*OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}y ** 2)							"))
tPropC[2][0][2].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[2][0][3].append(Template("-(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y * P${iloc}z /3.	"))
tPropC[2][0][3].append(Template("4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[2][0][3].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[2][0][3].append(Template("2*OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}y + 2./3.*OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2)"))

tPropC[2][1][0].append(Template("2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}y + 2./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}y"))
tPropC[2][1][0].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}y /3. - OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[2][1][0].append(Template("-OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}x /3."))
tPropC[2][1][0].append(Template("4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[2][1][1].append(Template("OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}y /3. - OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[2][1][1].append(Template("-4./3. * (-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}x * P${iloc}y							"))
tPropC[2][1][1].append(Template("(-1 - OM${iloc} * P${iloc}x ** 2) * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y ** 2 /3.	"))
tPropC[2][1][1].append(Template("-(-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y * P${iloc}z /3."))
tPropC[2][1][2].append(Template("-OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}x /3."))
tPropC[2][1][2].append(Template("(-1 - OM${iloc} * P${iloc}x ** 2) * (-1 - OM${iloc} * P${iloc}y ** 2) + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y ** 2 /3.	"))
tPropC[2][1][2].append(Template("-4./3.*OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}y ** 2)							"))
tPropC[2][1][2].append(Template("OM${iloc}**2 * P${iloc}x * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[2][1][3].append(Template("4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[2][1][3].append(Template("-(-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y * P${iloc}z /3."))
tPropC[2][1][3].append(Template("OM${iloc}**2 * P${iloc}x * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[2][1][3].append(Template("2*OM${iloc}**2 * P${iloc}x * P${iloc}z ** 2 * P${iloc}y + 2./3.*OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2)"))

tPropC[2][2][0].append(Template("2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y ** 2 - 2./3. * (1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}y ** 2)		"))
tPropC[2][2][0].append(Template("2*OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}x + 2./3.*OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}y ** 2)	"))
tPropC[2][2][0].append(Template("-4./3.*OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}y ** 2)								"))
tPropC[2][2][0].append(Template("2*OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}z + 2./3.*OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)	"))
tPropC[2][2][1].append(Template("2*OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}x + 2./3.*OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}y ** 2)	"))
tPropC[2][2][1].append(Template("2*OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y ** 2 - 2./3. * (-1 - OM${iloc} * P${iloc}x ** 2) * (-1 - OM${iloc} * P${iloc}y ** 2)	"))
tPropC[2][2][1].append(Template("-4./3.*OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}y ** 2)								"))
tPropC[2][2][1].append(Template("2*OM${iloc}**2 * P${iloc}x * P${iloc}y ** 2 * P${iloc}z + 2./3.*OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)	"))
tPropC[2][2][2].append(Template("-4./3.*OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}y ** 2)								"))
tPropC[2][2][2].append(Template("-4./3.*OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}y ** 2)								"))
tPropC[2][2][2].append(Template("4./3. * (-1 - OM${iloc} * P${iloc}y ** 2) ** 2												"))
tPropC[2][2][2].append(Template("-4./3. * (-1 - OM${iloc} * P${iloc}y ** 2)*OM${iloc} * P${iloc}y * P${iloc}z								"))
tPropC[2][2][3].append(Template("2*OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}z + 2./3.*OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)	"))
tPropC[2][2][3].append(Template("2*OM${iloc}**2 * P${iloc}x * P${iloc}y ** 2 * P${iloc}z + 2./3.*OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)	"))
tPropC[2][2][3].append(Template("-4./3. * (-1 - OM${iloc} * P${iloc}y ** 2)*OM${iloc} * P${iloc}y * P${iloc}z								"))
tPropC[2][2][3].append(Template("2*OM${iloc}**2 * P${iloc}y ** 2 * P${iloc}z ** 2 - 2./3. * (-1 - OM${iloc} * P${iloc}y ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2)"))

tPropC[2][3][0].append(Template(" 2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y * P${iloc}z + 2./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}y * P${iloc}z"))
tPropC[2][3][0].append(Template(" 4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[2][3][0].append(Template(" OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[2][3][0].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[2][3][1].append(Template(" 4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[2][3][1].append(Template(" 2*OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y * P${iloc}z + 2./3. * (-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}y * P${iloc}z"))
tPropC[2][3][1].append(Template(" OM${iloc}**2 * P${iloc}x * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[2][3][1].append(Template(" -OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}x * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[2][3][2].append(Template(" OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[2][3][2].append(Template(" OM${iloc}**2 * P${iloc}x * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[2][3][2].append(Template(" -4./3. * (-1 - OM${iloc} * P${iloc}y ** 2)*OM${iloc} * P${iloc}y * P${iloc}z							"))
tPropC[2][3][2].append(Template(" (-1 - OM${iloc} * P${iloc}y ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}y ** 2 * P${iloc}z ** 2 /3.	"))
tPropC[2][3][3].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[2][3][3].append(Template(" -OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}x * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[2][3][3].append(Template(" (-1 - OM${iloc} * P${iloc}y ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}y ** 2 * P${iloc}z ** 2 /3.	"))
tPropC[2][3][3].append(Template(" -4./3.*OM${iloc} * P${iloc}y * P${iloc}z * (-1 - OM${iloc} * P${iloc}z ** 2)							"))

tPropC[3][0][0].append(Template(" -4./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}t * P${iloc}z							"))
tPropC[3][0][0].append(Template(" -(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}z /3.	"))
tPropC[3][0][0].append(Template(" -(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y * P${iloc}z /3.	"))
tPropC[3][0][0].append(Template(" (1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}z ** 2 /3.	"))
tPropC[3][0][1].append(Template(" -(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}z /3.	"))
tPropC[3][0][1].append(Template(" 2*OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}z + 2./3.*OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[3][0][1].append(Template(" 4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[3][0][1].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}x /3."))
tPropC[3][0][2].append(Template(" -(1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y * P${iloc}z /3.	"))
tPropC[3][0][2].append(Template(" 4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[3][0][2].append(Template(" 2*OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}z + 2./3.*OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[3][0][2].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[3][0][3].append(Template(" (1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}z ** 2 /3.	"))
tPropC[3][0][3].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}x /3."))
tPropC[3][0][3].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[3][0][3].append(Template(" -4./3.*OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}z ** 2)"))

tPropC[3][1][0].append(Template(" 2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}x * P${iloc}z + 2./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}x * P${iloc}z"))
tPropC[3][1][0].append(Template(" OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[3][1][0].append(Template(" 4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[3][1][0].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}x /3."))
tPropC[3][1][1].append(Template(" OM${iloc}**2 * P${iloc}t * P${iloc}x ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}x ** 2)"))
tPropC[3][1][1].append(Template(" -4./3. * (-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}x * P${iloc}z							"))
tPropC[3][1][1].append(Template(" -(-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y * P${iloc}z /3."))
tPropC[3][1][1].append(Template(" (-1 - OM${iloc} * P${iloc}x ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}z ** 2 /3.	"))
tPropC[3][1][2].append(Template(" 4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z								"))
tPropC[3][1][2].append(Template(" -(-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}y * P${iloc}z + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y * P${iloc}z /3."))
tPropC[3][1][2].append(Template(" 2*OM${iloc}**2 * P${iloc}x * P${iloc}y ** 2 * P${iloc}z + 2./3.*OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[3][1][2].append(Template(" -OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}x * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[3][1][3].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}x /3."))
tPropC[3][1][3].append(Template(" (-1 - OM${iloc} * P${iloc}x ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}z ** 2 /3.	"))
tPropC[3][1][3].append(Template(" -OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}x * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[3][1][3].append(Template(" -4./3.*OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}z ** 2)							"))

tPropC[3][2][0].append(Template(" 2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}y * P${iloc}z + 2./3. * (1 - OM${iloc} * P${iloc}t ** 2)*OM${iloc} * P${iloc}y * P${iloc}z"))
tPropC[3][2][0].append(Template(" 4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z"))
tPropC[3][2][0].append(Template(" OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[3][2][0].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[3][2][1].append(Template(" 4./3.*OM${iloc}**2 * P${iloc}t * P${iloc}y * P${iloc}x * P${iloc}z"))
tPropC[3][2][1].append(Template(" 2*OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}y * P${iloc}z + 2./3. * (-1 - OM${iloc} * P${iloc}x ** 2)*OM${iloc} * P${iloc}y * P${iloc}z"))
tPropC[3][2][1].append(Template(" OM${iloc}**2 * P${iloc}x * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[3][2][1].append(Template(" -OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}x * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[3][2][2].append(Template(" OM${iloc}**2 * P${iloc}t * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[3][2][2].append(Template(" OM${iloc}**2 * P${iloc}x * P${iloc}y ** 2 * P${iloc}z /3. - OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}y ** 2)"))
tPropC[3][2][2].append(Template(" -4./3. * (-1 - OM${iloc} * P${iloc}y ** 2)*OM${iloc} * P${iloc}y * P${iloc}z"))
tPropC[3][2][2].append(Template(" (-1 - OM${iloc} * P${iloc}y ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}y ** 2 * P${iloc}z ** 2 /3.	"))
tPropC[3][2][3].append(Template(" -OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[3][2][3].append(Template(" -OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}x * P${iloc}z ** 2 * P${iloc}y /3."))
tPropC[3][2][3].append(Template(" (-1 - OM${iloc} * P${iloc}y ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2) + OM${iloc}**2 * P${iloc}y ** 2 * P${iloc}z ** 2 /3.	"))
tPropC[3][2][3].append(Template(" -4./3.*OM${iloc} * P${iloc}y * P${iloc}z * (-1 - OM${iloc} * P${iloc}z ** 2)"))

tPropC[3][3][0].append(Template(" 2*OM${iloc}**2 * P${iloc}t ** 2 * P${iloc}z ** 2 - 2./3. * (1 - OM${iloc} * P${iloc}t ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[3][3][0].append(Template(" 2*OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}x + 2./3.*OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[3][3][0].append(Template(" 2*OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}y + 2./3.*OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[3][3][0].append(Template(" -4./3.*OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}z ** 2)							"))
tPropC[3][3][1].append(Template(" 2*OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}x + 2./3.*OM${iloc} * P${iloc}t * P${iloc}x * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[3][3][1].append(Template(" 2*OM${iloc}**2 * P${iloc}x ** 2 * P${iloc}z ** 2 - 2./3. * (-1 - OM${iloc} * P${iloc}x ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[3][3][1].append(Template(" 2*OM${iloc}**2 * P${iloc}x * P${iloc}z ** 2 * P${iloc}y + 2./3.*OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[3][3][1].append(Template(" -4./3.*OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}z ** 2)							"))
tPropC[3][3][2].append(Template(" 2*OM${iloc}**2 * P${iloc}t * P${iloc}z ** 2 * P${iloc}y + 2./3.*OM${iloc} * P${iloc}t * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[3][3][2].append(Template(" 2*OM${iloc}**2 * P${iloc}x * P${iloc}z ** 2 * P${iloc}y + 2./3.*OM${iloc} * P${iloc}x * P${iloc}y * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[3][3][2].append(Template(" 2*OM${iloc}**2 * P${iloc}y ** 2 * P${iloc}z ** 2 - 2./3. * (-1 - OM${iloc} * P${iloc}y ** 2) * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[3][3][2].append(Template(" -4./3.*OM${iloc} * P${iloc}y * P${iloc}z * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[3][3][3].append(Template(" -4./3.*OM${iloc} * P${iloc}t * P${iloc}z * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[3][3][3].append(Template(" -4./3.*OM${iloc} * P${iloc}x * P${iloc}z * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[3][3][3].append(Template(" -4./3.*OM${iloc} * P${iloc}y * P${iloc}z * (-1 - OM${iloc} * P${iloc}z ** 2)"))
tPropC[3][3][3].append(Template(" 4./3. * (-1 - OM${iloc} * P${iloc}z ** 2) ** 2"))

imap=["t","x","y","z"]

RSDotProduct = Template("${s}ts${si}*${v}t-${s}xs${si}*${v}x-${s}ys${si}*${v}y-${s}zs${si}*${v}z")

vTemplateT="""\
{header} {{
         if({type}W{iloc}.id()=={id}) {{
            return {normal};
         }}
         else {{
            return {transpose};
         }}
     }};
"""

vTemplate4="""\
{header} {{
{swap}
    if(id{iloc1}=={id1}) {{
        if(id{iloc2}=={id2}) {{
          return {res1};
        }}
        else {{
          return {res2};
        }}
    }}
    else {{
        if(id{iloc2}=={id2}) {{
          return {res3};
        }}
        else {{
          return {res4};
        }}
    }}
}};

"""

vecTemplate="""\
    Energy2 p2 = P{iloc}.m2();
    LorentzPolarizationVector vtemp = {res};
    Complex fact = -Complex(0.,1.)*({cf})*propagator(iopt,p2,out,mass,width);
    if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
    complex<Energy2> mass2 = sqr(mass);
    if(mass.real()==ZERO) {{
        vtemp =fact*vtemp;
    }}
    else {{
        complex<Energy> dot = P{iloc}*vtemp;
        vtemp = fact*(vtemp-dot/mass2*P{iloc});
    }}
    return VectorWaveFunction(P{iloc},out,vtemp.x(),vtemp.y(),vtemp.z(),vtemp.t());
"""


sTemplate="""\
    Energy2 p2 = P{iloc}.m2();
    if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
    Complex fact = Complex(0.,1.)*({cf})*propagator(iopt,p2,out,mass,width);
    Lorentz{offTypeA}<double> newSpin = fact*({res});
    return {offTypeB}(P{iloc},out,newSpin);
"""

RSTemplate="""\
    if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
    Energy2 p2 = P{iloc}.m2();
    Complex fact = Complex(0.,-1.)*({cf})*propagator(iopt,p2,out,mass,width);
    complex<InvEnergy> Omass = mass.real()==ZERO ? InvEnergy(ZERO) : 1./mass;
    Lorentz{offTypeA}<double> newSpin = fact*({res});
    return {offTypeB}(P{iloc},out,newSpin);
"""

scaTemplate="""\
    if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
    Energy2 p2 = P{iloc}.m2();
    Complex fact = Complex(0.,1.)*({cf})*propagator(iopt,p2,out,mass,width);
    complex<double> output = fact*({res});
    return ScalarWaveFunction(P{iloc},out,output);
"""

tenTemplate="""\
    if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
    InvEnergy2 OM{iloc} = mass.real()==ZERO ? InvEnergy2(ZERO) : 1./sqr(mass.real());
    Energy2 M{iloc}2 = sqr(mass.real());
    Energy2 p2 = P{iloc}.m2();
    Complex fact = Complex(0.,1.)*({cf})*propagator(iopt,p2,out,mass,width);
    LorentzTensor<double> output = fact*({res});
    return TensorWaveFunction(P{iloc},out,output);
"""

# various strings for matrixes
I4 = "Matrix([[1.,0,0,0],[0,1.,0,0],[0,0,1.,0],[0,0,0,1.]])"
G5 = "Matrix([[-1.,0,0,0],[0,-1.,0,0],[0,0,1.,0],[0,0,0,1.]])"
PM = "Matrix([[1.,0,0,0],[0,1.,0,0],[0,0,0,0],[0,0,0,0]])"
PP = "Matrix([[0,0,0,0],[0,0,0,0],[0,0,1.,0],[0,0,0,1.]])"


vslash  = Template("Matrix([[0,0,${v}TMZ,-${v}XMY],[0,0,-${v}XPY,${v}TPZ],[${v}TPZ,${v}XMY,0,0],[${v}XPY,${v}TMZ,0,0]])")
vslashS = Template("${v}TPZ=Symbol(\"${v}TPZ\")\n${v}TMZ=Symbol(\"${v}TMZ\")\n${v}XPY=Symbol(\"${v}XPY\")\n${v}XMY=Symbol(\"${v}XMY\")\n")
momCom  = Template("${v}t = Symbol(\"${v}t\")\n${v}x = Symbol(\"${v}x\")\n${v}y = Symbol(\"${v}y\")\n${v}z = Symbol(\"${v}z\")\n")
vslashD = Template("complex<${var}> ${v}TPZ = ${v}.t()+${v}.z();\n    complex<${var}> ${v}TMZ = ${v}.t()-${v}.z();\n    complex<${var}> ${v}XPY = ${v}.x()+Complex(0.,1.)*${v}.y();\n    complex<${var}> ${v}XMY = ${v}.x()-Complex(0.,1.)*${v}.y();")
vslashM  = Template("Matrix([[$m,0,${v}TMZ,-${v}XMY],[0,$m,-${v}XPY,${v}TPZ],[${v}TPZ,${v}XMY,$m,0],[${v}XPY,${v}TMZ,0,$m]])")
vslashM2 = Template("Matrix([[$m,0,-${v}TMZ,${v}XMY],[0,$m,${v}XPY,-${v}TPZ],[-${v}TPZ,-${v}XMY,$m,0],[-${v}XPY,-${v}TMZ,0,$m]])")
vslashMS = Template("${v}TPZ=Symbol(\"${v}TPZ\")\n${v}TMZ=Symbol(\"${v}TMZ\")\n${v}XPY=Symbol(\"${v}XPY\")\n${v}XMY=Symbol(\"${v}XMY\")\n${m}=Symbol(\"${m}\")\nO${m}=Symbol(\"O${m}\")\n")

rslash   = Template("Matrix([[$m,0,${v}TMZ,-${v}XMY],[0,$m,-${v}XPY,${v}TPZ],[${v}TPZ,${v}XMY,$m,0],[${v}XPY,${v}TMZ,0,$m]])*( ($${eta}-2./3.*O${m}**2*${v}$${A}*${v}$${B})*Matrix([[1.,0,0,0],[0,1.,0,0],[0,0,1.,0],[0,0,0,1.]]) -1./3.*$${DA}*$${DB} -1./3.*O${m}*(${v}$${B}*$${DA}-${v}$${A}*$${DB}))")
rslashB  = Template("Matrix([[$m,0,${v}TMZ,-${v}XMY],[0,$m,-${v}XPY,${v}TPZ],[${v}TPZ,${v}XMY,$m,0],[${v}XPY,${v}TMZ,0,$m]])*( (${v2}$${A}-2./3.*O${m}**2*${v}$${A}*${dot})*Matrix([[1.,0,0,0],[0,1.,0,0],[0,0,1.,0],[0,0,0,1.]]) -1./3.*$${DA}*${DB} -1./3.*O${m}*(${dot}*$${DA}-${v}$${A}*${DB}))")



rslash2  = Template("Matrix([[$m,0,-${v}TMZ,${v}XMY],[0,$m,${v}XPY,-${v}TPZ],[-${v}TPZ,-${v}XMY,$m,0],[-${v}XPY,-${v}TMZ,0,$m]])*( ($${eta}-2./3.*O${m}**2*${v}$${A}*${v}$${B})*Matrix([[1.,0,0,0],[0,1.,0,0],[0,0,1.,0],[0,0,0,1.]]) -1./3.*$${DA}*$${DB} +1./3.*O${m}*(${v}$${B}*$${DA}-${v}$${A}*$${DB}))")
rslash2B  = Template("Matrix([[$m,0,-${v}TMZ,${v}XMY],[0,$m,${v}XPY,-${v}TPZ],[-${v}TPZ,-${v}XMY,$m,0],[-${v}XPY,-${v}TMZ,0,$m]])*( (${v2}$${B}-2./3.*O${m}**2*${dot}*${v}$${B})*Matrix([[1.,0,0,0],[0,1.,0,0],[0,0,1.,0],[0,0,0,1.]]) -1./3.*${DA}*$${DB} +1./3.*O${m}*(${v}$${B}*${DA}-${dot}*$${DB}))")

dirac=["Matrix([[0,0,1.,0],[0,0,0,1.],[1.,0,0,0],[0,1.,0,0]])","Matrix([[0,0,0,1.],[0,0,1.,0],[0,-1.,0,0],[-1.,0,0,0]])",
       "Matrix([[0,0,0,complex(0, -1.)],[0,0,complex(0, 1.),0],[0,complex(0, 1.),0,0],[complex(0, -1.),0,0,0]])",
       "Matrix([[0,0,1.,0],[0,0,0,-1.],[-1.,0,0,0],[0,1.,0,0]])"]
CC = "Matrix([[0,1.,0,0],[-1.,0,0,0],[0,0,0,-1.],[0,0,1.,0]])"
CD = "Matrix([[0,-1.,0,0],[1.,0,0,0],[0,0,0,1.],[0,0,-1.,0]])"

evaluateTemplate = """\
{decl} {{
    {momenta}
    {waves}
{swap}
    {symbols}
    {couplings}
{defns}
    {result}
}}
"""
spinor   = Template("Matrix([[${s}s1],[${s}s2],[${s}s3],[${s}s4]])")
sbar     = Template("Matrix([[${s}s1,${s}s2,${s}s3,${s}s4]])")
sline    = Template("${s}s1=Symbol(\"${s}s1\")\n${s}s2=Symbol(\"${s}s2\")\n${s}s3=Symbol(\"${s}s3\")\n${s}s4=Symbol(\"${s}s4\")\n")

RSSpinorTemplate = Template("${type}<double>(${outxs1},\n${outxs2},\n${outxs3},\n${outxs4},\n${outys1},\n${outys2},\n${outys3},\n${outys4},\n${outzs1},\n${outzs2},\n${outzs3},\n${outzs4},\n${outts1},\n${outts2},\n${outts3},\n${outts4})")
SpinorTemplate = Template("${type}<double>(${outs1},\n${outs2},\n${outs3},\n${outs4})")

class LorentzIndex :
    """ A simple classs to store a Lorentz index """
    type=""
    value=0
    dimension=0
    def __repr__(self):
        if(self.type=="V" and not isinstance(self.value,int)) :
            return self.value
        else :
            return "%s%s" % (self.type,self.value)

    def __init__(self,val) :
        if(isinstance(val,int)) :
            self.dimension=0
            if(val<0) :
                self.type="D"
                self.value = val
            elif(val>0 and val/1000==0) :
                self.type="E"
                self.value = val
            elif(val>0 and val/1000==1) :
                self.type="T1"
                self.value = val%1000
            elif(val>0 and val/1000==2) :
                self.type="T2"
                self.value = val%1000
            else :
                print("Unknown value in Lorentz index:",val)
                raise SkipThisVertex()
        else :
            print("Unknown value in Lorentz index:",val)
            raise SkipThisVertex()
            
    def __eq__(self,other):
        if(not isinstance(other,LorentzIndex)) :
            return False
        return ( (self.type, self.value) 
                 == (other.type, other.value) )

    def __hash__(self) :
        return hash((self.type,self.value))

class DiracMatrix:
    """A simple class to store Dirac matrices"""
    name =""
    value=""
    index=0
    def __init(self) :
        self.name=""
        self.value=""
        self.index=0

    def __repr__(self) :
        if(self.value==0) :
            return "%s" % self.index
        else :
            return "%s" % self.value
    
class LorentzStructure:
    """A simple class to store a Lorentz structures"""
    name=""
    value=0
    lorentz=[]
    spin=[]

    def __init(self) :
        self.name=""
        self.value=0
        self.lorentz=[]
        self.spin=[]
    
    def __repr__(self):
        output = self.name
        if((self.name=="P" or self.name=="Tensor") and self.value!=0) :
            output += "%s" % self.value
        if(self.name=="int" or self.name=="sign") :
            output += "=%s" % self.value
        elif(len(self.spin)==0) :
            output += "("
            for val in self.lorentz :
                output += "%s," % val
            output=output.rstrip(",")
            output+=")"
        elif(len(self.lorentz)==0) :
            output += "("
            for val in self.spin :
                output += "%s," % val
            output=output.rstrip(",")
            output+=")"
        else :
            output += "("
            for val in self.lorentz :
                output += "%s," % val
            for val in self.spin :
                output += "%s," % val
            output=output.rstrip(",")
            output+=")"
        return output

def LorentzCompare(a,b) :
    if(a.name=="int" and b.name=="int") :
        return int(abs(b.value)-abs(a.value))
    elif(a.name=="int") :
        return -1
    elif(b.name=="int") :
        return 1
    elif(len(a.spin)==0) :
        if(len(b.spin)==0) :
            return len(b.lorentz)-len(a.lorentz)
        else :
            return -1
    elif(len(b.spin)==0) :
         return 1
    else :
        if(len(a.spin)==0 or len(b.spin)==0) :
            print('Index problem in lorentz compare',
                  a.name,b.name,a.spin,b.spin)
            raise SkipThisVertex()
        if(a.spin[0]>0 or b.spin[1]>0 ) : return -1
        if(a.spin[1]>0 or b.spin[0]>0 ) : return  1
        if(a.spin[1]==b.spin[0]) : return -1
        if(b.spin[1]==a.spin[0]) : return 1
    return 0

def extractIndices(struct) :
    if(struct.find("(")<0) : return []
    temp=struct.split("(")[1].split(")")[0].split(",")
    output=[]
    for val in temp :
        output.append(int(val))
    return output

def parse_structure(structure,spins) :
    output=[]
    found = True
    while(found) :
        found = False
        # signs between terms
        if(structure=="+" or structure=="-") :
            output.append(LorentzStructure())
            output[0].name="sign"
            output[0].value=structure[0]+"1."
            output[0].value=float(output[0].value)
            output[0].lorentz=[]
            output[0].spin=[]
            return output
        # simple numeric pre/post factors
        elif((structure[0]=="-" or structure[0]=="+") and
             structure[-1]==")" and structure[1]=="(") :
            output.append(LorentzStructure())
            output[-1].name="int"
            output[-1].value=structure[0]+"1."
            output[-1].value=float(output[-1].value)
            output[-1].lorentz=[]
            output[-1].spin=[]
            structure=structure[2:-1]
            found=True
        elif(structure[0]=="(") :
            temp=structure.rsplit(")",1)
            structure=temp[0][1:]
            output.append(LorentzStructure())
            output[-1].name="int"
            output[-1].value="1."+temp[1]
            output[-1].value=float(eval(output[-1].value))
            output[-1].lorentz=[]
            output[-1].spin=[]
            found=True
        elif(structure[0:2]=="-(") :
            temp=structure.rsplit(")",1)
            structure=temp[0][2:]
            output.append(LorentzStructure())
            output[-1].name="int"
            output[-1].value="-1."+temp[1]
            output[-1].value=float(eval(output[-1].value))
            output[-1].lorentz=[]
            output[-1].spin=[]
            found=True
    # special handling for powers
    power = False
    if("**" in structure ) :
        power = True
        structure = structure.replace("**","^")
    structures = structure.split("*")
    if(power) :
        for j in range(0,len(structures)):
            if(structures[j].find("^")>=0) :
                temp = structures[j].split("^")
                structures[j] = temp[0]
                for i in range(0,int(temp[1])-1) :
                    structures.append(temp[0])
    # split up the structure
    for struct in structures:
        ind = extractIndices(struct)
        # different types of object
        # object with only spin indices
        if(struct.find("Identity")==0 or
           struct.find("Proj")==0 or
           struct.find("Gamma5")==0) :
            output.append(LorentzStructure())
            output[-1].spin=ind
            output[-1].lorentz=[]
            output[-1].name=struct.split("(")[0]
            output[-1].value=0
            if(len(struct.replace("%s(%s,%s)" % (output[-1].name,ind[0],ind[1]),""))!=0) :
                print("Problem handling %s structure " % output[-1].name)
                raise SkipThisVertex()
        # objects with 2 lorentz indices
        elif(struct.find("Metric")==0) :
            output.append(LorentzStructure())
            output[-1].lorentz=[LorentzIndex(ind[0]),LorentzIndex(ind[1])]
            output[-1].name=struct.split("(")[0]
            output[-1].value=0
            output[-1].spin=[]
            if(len(struct.replace("%s(%s,%s)" % (output[-1].name,ind[0],ind[1]),""))!=0) :
                print("Problem handling %s structure " % output[-1].name)
                raise SkipThisVertex()
        elif(struct.find("P(")==0) :
            output.append(LorentzStructure())
            output[-1].lorentz=[LorentzIndex(ind[0])]
            output[-1].name=struct.split("(")[0]
            output[-1].value=ind[1]
            output[-1].spin=[]
            if(len(struct.replace("%s(%s,%s)" % (output[-1].name,ind[0],ind[1]),""))!=0) :
                print("Problem handling %s structure " % output[-1].name)
                raise SkipThisVertex()
        # 1 lorentz and 1 spin index
        elif(struct.find("Gamma")==0) :
            output.append(LorentzStructure())
            output[-1].lorentz=[LorentzIndex(ind[0])]
            output[-1].spin=[ind[1],ind[2]]
            output[-1].name=struct.split("(")[0]
            output[-1].value=1
            if(len(struct.replace("%s(%s,%s,%s)" % (output[-1].name,ind[0],ind[1],ind[2]),""))!=0) :
                print("problem parsing gamma matrix",struct)
                raise SkipThisVertex()
        # objects with 4 lorentz indices
        elif(struct.find("Epsilon")==0) :
            output.append(LorentzStructure())
            output[-1].lorentz=[]
            for i in range(0,len(ind)) :
                output[-1].lorentz.append(LorentzIndex(ind[i]))
            output[-1].spin=[]
            output[-1].name=struct.split("(")[0]
            output[-1].value=1
            if(len(struct.replace("%s(%s,%s,%s,%s)" % (output[-1].name,ind[0],ind[1],ind[2],ind[3]),""))!=0) :
                print('Problem parsing epsilon',struct)
                raise SkipThisVertex()
        # scalars
        else :
            try :
                output.append(LorentzStructure())
                output[-1].value=float(struct)
                output[-1].name="int"
                output[-1].lorentz=[]
                output[-1].spin=[]
            except :
                if(struct.find("complex")==0) :
                    vals = struct[0:-1].replace("complex(","").split(",")
                    output[-1].value=complex(float(vals[0]),float(vals[1]))
                    output[-1].name="int"
                    output[-1].lorentz=[]
                    output[-1].spin=[]
                else :
                    print('Problem parsing scalar',struct)
                    raise SkipThisVertex()
    # now do the sorting
    if(len(output)==1) : return output
    output = sorted(output,cmp=LorentzCompare)
    # fix indices in the RS case
    if(4 in spins) :
        for i in range(0,len(output)) :
            for ll in range(0,len(output[i].lorentz)) :
                if(spins[output[i].lorentz[ll].value-1]==4 and
                   output[i].lorentz[ll].type=="E") :
                    output[i].lorentz[ll].type="R"
    # return the answer
    return output

def constructDotProduct(ind1,ind2,defns) :
    (ind1,ind2) = sorted((ind1,ind2),cmp=indSort)
    dimension=ind1.dimension+ind2.dimension
    # this product already dealt with ?
    if((ind1,ind2) in defns) :
        name = defns[(ind1,ind2)][0]
    # handle the product
    else :
        name = "dot%s" % (len(defns)+1)
        unit = computeUnit(dimension)
        defns[(ind1,ind2)] = [name,"complex<%s> %s = %s*%s;" % (unit,name,ind1,ind2)]
    return (name,dimension)

def contract(parsed) :
    for j in range(0,len(parsed)) :
        if(parsed[j]=="") : continue
        if(parsed[j].name=="P") :
            # simplest case
            if(parsed[j].lorentz[0].type=="E" or
               parsed[j].lorentz[0].type=="P") :
                newIndex = LorentzIndex(parsed[j].value)
                newIndex.type="P"
                newIndex.dimension=1
                parsed[j].name="Metric"
                parsed[j].lorentz.append(newIndex)
                parsed[j].lorentz = sorted(parsed[j].lorentz,cmp=indSort)
                continue
            ll=1
            found=False
            for k in range(0,len(parsed)) :
                if(j==k or parsed[k]=="" ) : continue
                for i in range(0,len(parsed[k].lorentz)) :
                    if(parsed[k].lorentz[i] == parsed[j].lorentz[0]) :
                        parsed[k].lorentz[i].type="P"
                        parsed[k].lorentz[i].value = parsed[j].value
                        parsed[k].lorentz[i].dimension=1
                        if(parsed[k].name=="P") :
                            parsed[k].lorentz.append(LorentzIndex(parsed[k].value))
                            parsed[k].lorentz[1].type="P"
                            parsed[k].lorentz[1].dimension=1
                            parsed[k].name="Metric"
                            parsed[k].value = 0
                        found=True
                        break
                if(found) :
                    parsed[j]=""
                    break
    return [x for x in parsed if x != ""]

def computeUnit(dimension) :
    if(isinstance(dimension,int)) :
        dtemp = dimension
    else :
        dtemp=dimension[1]+dimension[2]
    if(dtemp==0) :
        unit="double"
    elif(dtemp==1) :
        unit="Energy"
    elif(dtemp==-1) :
        unit="InvEnergy"
    elif(dtemp>0) :
        unit="Energy%s" % (dtemp)
    elif(dtemp<0) :
        unit="InvEnergy%s" % (dtemp)
    return unit

def computeUnit2(dimension,vDim) :
    # first correct for any coupling power in vertex
    totalDim = int(dimension[0])+dimension[2]+vDim-4
    output=""
    if(totalDim!=0) :
        if(totalDim>0) :
            if(totalDim==1) :
                output = "1./GeV"
            elif(totalDim==2) :
                output = "1./GeV2"
            else :
                output="1."
                for i in range(0,totalDim) :
                    output +="/GeV"
        else :
            if(totalDim==-1) :
                output = "GeV"
            elif(totalDim==-2) :
                output = "GeV2"
            else :
                output="1."
                for i in range(0,-totalDim) :
                    output +="*GeV"
    expr=""
    # now remove the remaining dimensionality
    removal=dimension[1]-int(dimension[0])-vDim+4
    if(removal!=0) :
        if(removal>0) :
            if(removal==1) :
                expr = "UnitRemovalInvE"
            else :
                expr = "UnitRemovalInvE%s" % removal
        else :
            if(removal==-1) :
                expr = "UnitRemovalE"
            else :
                expr = "UnitRemovalE%s" % (-removal)
    if(output=="") : return expr
    elif(expr=="") : return output
    else           : return "%s*%s" %(output,expr)

# order the indices of a dot product
def indSort(a,b) :
    if(not isinstance(a,LorentzIndex) or
       not isinstance(b,LorentzIndex)) :
       print("Trying to sort something that's not a Lorentz index",a,b)
       raise SkipThisVertex()
    if(a.type==b.type) :
        i1=a.value
        i2=b.value
        if(i1>i2) :
            return 1
        elif(i1<i2) :
            return -1
        else :
            return 0
    else :
        if(a.type=="E") :
            return 1
        else :
            return -1

def finishParsing(parsed,dimension,lorentztag,iloc,defns,eps) :
    output=1.
    # replace signs
    if(len(parsed)==1 and parsed[0].name=="sign") :
        if(parsed[0].value>0) :
            output="+"
        else :
            output="-"
        parsed=[]
        return (output,parsed,dimension,eps)
    # replace integers (really lorentz scalars)
    for j in range(0,len(parsed)) :
        if(parsed[j]!="" and parsed[j].name=="int") :
            output *= parsed[j].value
            parsed[j]=""
    # bracket this for safety
    if(output!="") : output = "(%s)" % output
    # special for tensor indices
    if("T" in lorentztag) :
        for j in range(0,len(parsed)) :
            if(parsed[j]=="") :continue
            # check for tensor index
            found=False
            for li in parsed[j].lorentz :
                if(li.type[0]=="T") : 
                    index = li
                    found=True
                    break
            if(not found) : continue
            # workout the other index for the tensor
            index2 = LorentzIndex(li.value)
            if(index.type=="T1") :
                index2.type="T2"
            else :
                index2.type="T1"
            # special is tensor contracted with itself
            if(parsed[j].name=="Metric" and index2 == parsed[j].lorentz[1]) :
                parsed[j].name = "Tensor"
                parsed[j].value = index.value
                parsed[j].lorentz = []
                if(iloc!=index.value) : 
                    name= "traceT%s" % parsed[j].value
                    if( name  in defns ) :
                        output += "*(%s)" % defns[name][0]
                    else :
                        defns[name] = [name,"Complex %s = T%s.trace();" % (name,parsed[j].value)]
                        output += "*(%s)" % defns[name][0]
                    parsed[j]=""
                continue
            # otherwise search for the match
            for k in range(j+1,len(parsed)) :
                found = False
                for li in parsed[k].lorentz :
                    if(li == index2) : 
                        found=True
                        break
                if(not found) : continue
                if(parsed[j].name=="P") :
                    newIndex1 = LorentzIndex(parsed[j].value)
                    newIndex1.type="P"
                    newIndex1.dimension=1
                elif(parsed[j].name=="Metric") :
                    for li in parsed[j].lorentz :
                        if(li != index) :
                            newIndex1=li
                            break
                else :
                    print('Unknown object with tensor index, first object',parsed[j])
                    raise SkipThisVertex()
                if(parsed[k].name=="P") :
                    newIndex2 = LorentzIndex(parsed[k].value)
                    newIndex2.type="P"
                    newIndex2.dimension=1
                elif(parsed[k].name=="Metric") :
                    for li in parsed[k].lorentz :
                        if(li != index) :
                            newIndex2=li
                            break
                elif(parsed[k].name=="Gamma") :
                    # if we can't contract
                    if(index.value==iloc or (newIndex1.type=="E" and newIndex1.value==iloc)) :
                        newIndex2=index2
                    # otherwise contract
                    else : 
                        unit=computeUnit(newIndex1.dimension)
                        if(index.type=="T1") :
                            name="T%s%sF" % (index.value,newIndex1)
                            defns[name] = [name,"LorentzVector<complex<%s> > %s = T%s.preDot(%s);" % (unit,name,index.value,newIndex1)]
                        else :
                            name="T%s%sS" % (index.value,newIndex1)
                            defns[name] = [name,"LorentzVector<complex<%s> > %s = T%s.postDot(%s);" % (unit,name,index.value,newIndex1)]
                        parsed[j]=""
                        gIndex=LorentzIndex(-1)
                        gIndex.type="V"
                        gIndex.value=name
                        gIndex.dimension=newIndex1.dimension
                        parsed[k].lorentz[0] = gIndex
                        break
                else :
                    print('Unknown object with tensor index, second object',parsed[j],parsed[k])
                    raise SkipThisVertex()
                if(index2.type=="T1") :
                    newIndex1,newIndex2=newIndex2,newIndex1
                parsed[j].name = "Tensor"
                parsed[j].value= int(index.value)
                parsed[j].lorentz= [newIndex1,newIndex2]
                if(parsed[k].name!="Gamma") : parsed[k]=""
                break
    # main handling of lorentz structures
    for j in range(0,len(parsed)) :
        if(parsed[j]=="") : continue
        if(parsed[j].name=="Metric") :
            # check whether or not we can contract
            canContract=True
            for ll in parsed[j].lorentz :
                if((ll.type=="E" and ll.value==iloc) or ll.type=="R") :
                    canContract = False
                    break
            if(not canContract) : continue
            # if we can do it
            (name,dTemp) = constructDotProduct(parsed[j].lorentz[0],parsed[j].lorentz[1],defns)
            output += "*(%s)" % name
            dimension[2] += dTemp
            parsed[j]=""
        elif(parsed[j].name=="Epsilon") :
            if(not eps) : eps = True
            # work out which, if any of the indices can be summed over
            summable=[]
            for ix in range(0,len(parsed[j].lorentz)) :
                if(parsed[j].lorentz[ix].type=="P" or
                   (parsed[j].lorentz[ix].type=="E" and iloc !=parsed[j].lorentz[ix].value)) :
                    summable.append(True)
                else :
                    summable.append(False)
            sc = summable.count(True)
            # less than 3 contractable indices, leave for later
            if(sc<3) :
                continue
            # can contract to a vector
            elif(sc==3) :
                offLoc = -1
                for i in range(0,len(summable)):
                    if(not summable[i]) :
                        offLoc = i
                        break
            else :
                offLoc = 0
            indices=[]
            dTemp=0
            for ix in range(0,len(parsed[j].lorentz)) :
                dTemp += parsed[j].lorentz[ix].dimension
                if(ix!=offLoc) : indices.append(parsed[j].lorentz[ix])
            # contract all the indices
            if(sc==4) :
                dimension[2] += dTemp
                iTemp = (parsed[j].lorentz[0],parsed[j].lorentz[1],
                         parsed[j].lorentz[2],parsed[j].lorentz[3])
                if(iTemp in defns) :
                    output += "*(%s)" % defns[iTemp][0]
                    parsed[j]=""
                else :
                    name = "dot%s" % (len(defns)+1)
                    unit = computeUnit(dTemp)
                    defns[iTemp] = [name,"complex<%s> %s =-%s*epsilon(%s,%s,%s);" % (unit,name,parsed[j].lorentz[0],
                                                                                     indices[0],indices[1],indices[2]) ]
                    output += "*(%s)" % name
                parsed[j]=""
            # contract 3 indices leaving a vector
            else :
                iTemp = (indices[0],indices[1],indices[2])
                sign = "1"
                if(offLoc%2!=0) : sign="-1"
                if(iTemp in defns) :
                    name = defns[iTemp][0]
                else :
                    name = "V%s" % (len(defns)+1)
                    unit = computeUnit(dTemp)
                    defns[iTemp] = [name,"LorentzVector<complex<%s> > %s =-epsilon(%s,%s,%s);" % (unit,name,
                                                                                                  indices[0],indices[1],indices[2]) ]
                newIndex = LorentzIndex(int(name[1:]))
                newIndex.type="V"
                newIndex.dimension=dTemp
                output += "*(%s)" % (sign)
                oi = parsed[j].lorentz[offLoc]
                if(oi.type!="D") :
                    parsed[j].name="Metric"
                    parsed[j].spins=[]
                    parsed[j].value=0
                    parsed[j].lorentz=[newIndex,oi]
                else :
                    found=False
                    for k in range(0,len(parsed)):
                        if(k==j or parsed[k]=="") : continue
                        for ll in range(0,len(parsed[k].lorentz)) :
                            if(parsed[k].lorentz[ll]==oi) :
                                found=True
                                parsed[k].lorentz[ll]=newIndex
                                break
                        if(found) : break
                    if(not found) :
                        print("Problem contracting indices of Epsilon tensor")
                        raise SkipThisVertex()
                    parsed[j]=""
        elif(parsed[j].name=="Tensor") :
            # not an external tensor
            if(parsed[j].value!=iloc) :
                # now check the lorentz indices
                con=[]
                uncon=[]
                dtemp=0
                for li in parsed[j].lorentz :
                    if(li.type=="P" or li.type=="V") :
                        con.append(li)
                        dtemp+=li.dimension
                    elif(li.type=="E") :
                        if(li.value!=iloc) :
                            con.append(li)
                        else :
                            uncon.append(li)
                    else :
                        print("Can't handle index ",li,"in tensor",parsed[j])
                        raise SkipThisVertex()
                if(len(con)==2) :
                    iTemp = ("T%s%s%s"% (parsed[j].value,con[0],con[1]))
                    dimension[2]+=dtemp
                    if(iTemp in defns) :
                        output += "*(%s)" % defns[iTemp][0]
                    else :
                        unit=computeUnit(dtemp)
                        name = "dot%s" % (len(defns)+1)
                        defns[iTemp] = [name,"complex<%s> %s = T%s.preDot(%s)*%s;" % (unit,name,parsed[j].value,con[0],con[1])]
                        output += "*(%s)" % name
                    parsed[j]=""
                # handled in final stage
                else :
                    continue
        elif(parsed[j].name.find("Proj")>=0 or
             parsed[j].name.find("Gamma")>=0 or
             parsed[j].name.find("Identity")>=0) :
            continue
        elif(parsed[j].name=="P" and parsed[j].lorentz[0].type=="R") :
            continue
        else :
            print('Lorentz structure',parsed[j],'not handled')
            raise SkipThisVertex()
    # remove leading *
    if(output!="" and output[0]=="*") : output = output[1:]
    # remove any (now) empty elements
    parsed = [x for x in parsed if x != ""]
    return (output,parsed,dimension,eps)

def finalContractions(output,parsed,dimension,lorentztag,iloc,defns) :
    if(len(parsed)==0) :
       return (output,dimension)
    elif(len(parsed)!=1) :
        print("Summation can't be handled",parsed)
        raise SkipThisVertex()
    if(parsed[0].name=="Tensor") :
        # contracted with off-shell vector
        if(parsed[0].value!=iloc) :
            found = False
            for ll in parsed[0].lorentz :
                if(ll.type=="E" and ll.value==iloc) :
                    found = True
                else :
                    lo=ll
            if(found) :
                dimension[2]+= lo.dimension
                unit=computeUnit(lo.dimension)
                if(lo==parsed[0].lorentz[0]) :
                    name="T%s%sF" % (parsed[0].value,lo)
                    defns[name] = [name,"LorentzVector<complex<%s> > %s = T%s.preDot(%s);" % (unit,name,parsed[0].value,lo)]
                else :
                    name="T%s%sS" % (parsed[0].value,lo)
                    defns[name] = [name,"LorentzVector<complex<%s> > %s = T%s.postDot(%s);" % (unit,name,parsed[0].value,lo)]
                parsed[0]=""
                if(output=="") : output="1."
                output = "(%s)*(%s)" %(output,name)
            else :
                print("Can\'t contract tensor",lo,iloc)
                raise SkipThisVertex()
        # off-shell tensor
        else :
            if(len(parsed[0].lorentz)!=0) :
                dimension[2]+=parsed[0].lorentz[0].dimension+parsed[0].lorentz[1].dimension
            tensor = tensorPropagator(parsed[0],defns)
            if(output=="") : output="1."
            output = [output,tensor,()]
    elif(parsed[0].name=="Metric") :
        found = False
        for ll in parsed[0].lorentz :
            if(ll.type=="E" and ll.value==iloc) :
                found = True
            else :
                lo=ll
        if(found) :
            parsed[0]=""
            dimension[2] += lo.dimension
            if(output=="") : output="1."
            output = "(%s)*(%s)" %(output,lo)
    else :
        print("Structure can't be handled",parsed,iloc)
        raise SkipThisVertex()
    return (output,dimension)
    
def tensorPropagator(struct,defns) :
    # dummy index
    i0 = LorentzIndex(-1000)
    # index for momentum of propagator
    ip = LorentzIndex(struct.value)
    ip.type="P"
    ip.dimension=1
    # the metric tensor
    terms=[]
    if(len(struct.lorentz)==0) :
        (dp,dTemp) = constructDotProduct(ip,ip,defns)
        pre = "-1./3.*(1.-%s*OM%s)" % (dp,struct.value)
        terms.append((pre,i0,i0))
        pre = "-2./3.*(1.-%s*OM%s)" % (dp,struct.value)
        terms.append(("%s*OM%s" %(pre,struct.value),ip,ip))
    else :
        # indices of the tensor
        ind1 = struct.lorentz[0]
        ind2 = struct.lorentz[1]
        # the dot products we need
        (d1,dtemp) = constructDotProduct(ind1,ip,defns)
        (d2,dtemp) = constructDotProduct(ind2,ip,defns)
        (d3,dtemp) = constructDotProduct(ind1,ind2,defns)
        # various terms in the propagator
        terms.append(("0.5",ind1,ind2))
        terms.append(("-0.5*OM%s*%s"%(struct.value,d1),ip,ind2))
        terms.append(("-0.5*OM%s*%s"%(struct.value,d2),ind1,ip))
        terms.append(("0.5",ind2,ind1))
        terms.append(("-0.5*OM%s*%s"%(struct.value,d2),ip,ind1))
        terms.append(("-0.5*OM%s*%s"%(struct.value,d1),ind2,ip))
        terms.append(("-1./3.*"+d3,i0,i0))
        terms.append(("1./3.*OM%s*%s*%s"%(struct.value,d1,d2),i0,i0))
        terms.append(("1./3.*OM%s*%s"%(struct.value,d3),ip,ip))
        terms.append(("2./3.*OM%s*OM%s*%s*%s"%(struct.value,struct.value,d1,d2),ip,ip))
    # compute the output as a dict
    output={}
    for i1 in imap:
        for i2 in imap :
            val=""
            for term in terms:
                if(term[0][0]!="-") :
                    pre = "+"+term[0]
                else :
                    pre = term[0]
                if(term[1]==i0) :
                    if(i1==i2) :
                        if(i1=="t") :
                            val += pre
                        else :
                            if(pre[0]=="+") :
                                val +="-"+pre[1:]
                            else :
                                val +="+"+pre[1:]
                                    
                            
                else :
                    val += "%s*%s%s*%s%s" % (pre,term[1],i1,term[2],i2)
            output["%s%s" % (i1,i2) ] = val.replace("+1*","+").replace("-1*","-")
    return output
        
def generateVertex(iloc,L,parsed,lorentztag,vertex,defns) :
    # try to import sympy and exit if required
    try :
        import sympy
        from sympy import Matrix,Symbol
    except :
        print('ufo2herwig uses the python sympy module to translate general lorentz structures.')
        print('This must be installed if you wish to use this option.')
        print('EXITTING')
        quit()
    eps=False
    # parse the lorentz structures
    output    = [1.]*len(parsed)
    dimension=[]
    for i in range(0,len(parsed)) :
        dimension.append([0,0,0])
    for i in range (0,len(parsed)) :
        (output[i],parsed[i],dimension[i],eps) = finishParsing(parsed[i],dimension[i],lorentztag,iloc,defns,eps)
    # still need to process gamma matrix strings for fermions
    if(lorentztag[0] in ["F","R"] ) :
        return convertDirac(output,dimension,eps,iloc,L,parsed,lorentztag,vertex,defns)
    # return the answer
    else :
        handled=True
        for i in range (0,len(parsed)) :
            if(len(parsed[i])!=0) :
                handled = False
                break
        if(not handled) :
            for i in range (0,len(parsed)) :
                (output[i],dimension[i]) = finalContractions(output[i],parsed[i],dimension[i],lorentztag,iloc,defns)
        return (output,dimension,eps)

def convertDirac(output,dimension,eps,iloc,L,parsed,lorentztag,vertex,defns) :
    for i in range(0,len(parsed)):
        # skip empty elements
        if(len(parsed[i])==0 or (len(parsed[i])==1 and parsed[i][0]=="")) : continue
        # parse the string
        (output[i],dimension[i],defns) = convertDiracStructure(parsed[i],output[i],dimension[i],
                                                               defns,iloc,L,lorentztag,vertex)
    return (output,dimension,eps)

# parse the gamma matrices
def convertMatrix(structure,spins,unContracted,Symbols,dtemp,defns,iloc) :
    i1 = structure.spin[0]
    i2 = structure.spin[1]
    if(structure.name=="Identity") :
        output = DiracMatrix()
        output.value = I4
        output.index=0
        output.name="M"
        structure=""
    elif(structure.name=="Gamma5") :
        output = DiracMatrix()
        output.value = G5
        output.index=0
        output.name="M"
        structure=""
    elif(structure.name=="ProjM") :
        output = DiracMatrix()
        output.value = PM
        output.index=0
        output.name="M"
        structure=""
    elif(structure.name=="ProjP") :
        output = DiracMatrix()
        output.value = PP
        output.index=0
        output.name="M"
        structure=""
    elif(structure.name=="Gamma") :
        # gamma(mu) lorentz matrix contracted with dummy index
        if(structure.lorentz[0].type=="D" or structure.lorentz[0].type=="R") :
            if(structure.lorentz[0] not in unContracted) :
                unContracted[structure.lorentz[0]] = 0
            output = DiracMatrix()
            output.value=0
            output.index=structure.lorentz[0]
            output.name="GMU"
            structure=""
        elif(structure.lorentz[0].type  == "E" and
             structure.lorentz[0].value == iloc ) :
            if(structure.lorentz[0] not in unContracted) :
                unContracted[structure.lorentz[0]] = 0
            output = DiracMatrix()
            output.value=0
            output.index=structure.lorentz[0]
            output.name="GMU"
            structure=""
        elif(structure.lorentz[0].type  == "T1" or
              structure.lorentz[0].type  == "T2") :
            if(structure.lorentz[0] not in unContracted) :
                unContracted[structure.lorentz[0]] = 0
            output = DiracMatrix()
            output.value=0
            output.index=structure.lorentz[0]
            output.name="GMU"
            structure=""
        else :
            output=DiracMatrix()
            output.name="M"
            output.value = vslash.substitute({ "v" : structure.lorentz[0]})
            Symbols += vslashS.substitute({ "v" : structure.lorentz[0]})
            variable = computeUnit(structure.lorentz[0].dimension)
            #if(structure.lorentz[0].type!="V" or
            #   structure.lorentz[0].type=="V") :
            dtemp[2] += structure.lorentz[0].dimension
            defns["vv%s" % structure.lorentz[0] ] = \
                                                    ["vv%s" % structure.lorentz[0],
                                                     vslashD.substitute({ "var" : variable,
                                                                          "v" : structure.lorentz[0]})]
            structure=""
    else :
        print('Unknown Gamma matrix structure',structure)
        raise SkipThisVertex()
    return (i1,i2,output,structure,Symbols)


def checkRSContract(parsed,loc,dtemp) :
    rindex=LorentzIndex(loc)
    rindex.type="R"
    contract=""
    for i in range(0,len(parsed)) :
        if(parsed[i]=="") : continue
        found = False
        for ll in range(0,len(parsed[i].lorentz)) :
            if(parsed[i].lorentz[ll]==rindex) :
                found=True
                break
        if(not found) :
            continue
        if(parsed[i].name=="P") :
            dtemp[2]+=1
            contract = LorentzIndex(parsed[i].value)
            contract.type="P"
            contract.dimension=1
            parsed[i]=""
            break
        elif(parsed[i].name=="Metric") :
            for ll in parsed[i].lorentz :
                if(ll==rindex) :
                    continue
                else :
                    break
            if(ll.type=="P") :
                dtemp[2]+=1
            contract=ll
            parsed[i]=""
            break
        elif(parsed[i].name=="Epsilon") :
            continue
        else :
            print("Unkonwn type contracted with RS spinor",parsed[i])
            raise SkipThisVertex()
    return contract

def processChain(dtemp,parsed,spins,Symbols,unContracted,defns,iloc) :
    # piece of dimension which is common (0.5 for sbar and spinor)
    dtemp[0]+=1
    # set up the spin indices
    sind = 0
    lind = 0
    expr=[]
    # now find the next thing in the matrix string
    ii = 0
    index=0
    while True :
        # already handled
        if(parsed[ii]=="") :
            ii+=1
            continue
        # start of the chain
        elif(sind==0 and len(parsed[ii].spin)==2 and parsed[ii].spin[0]>0 ) :
            (sind,index,matrix,parsed[ii],Symbols) \
                = convertMatrix(parsed[ii],spins,unContracted,Symbols,dtemp,defns,iloc)
            expr.append(matrix)
        # next element in the chain
        elif(index!=0 and len(parsed[ii].spin)==2 and parsed[ii].spin[0]==index) :
            (crap,index,matrix,parsed[ii],Symbols) \
                = convertMatrix(parsed[ii],spins,unContracted,Symbols,dtemp,defns,iloc)
            expr.append(matrix)
        # check the index to see if we're at the end
        if(index>0) :
            lind=index
            break
        ii+=1
        if(ii>=len(parsed)) :
            print("Can't parsed the chain of dirac matrices")
            print(parsed)
            raise SkipThisVertex()
    # start and end of the spin chains
    # first particle spin 1/2
    if(spins[sind-1]==2) :
        start = DiracMatrix()
        endT  = DiracMatrix()
        start.index=0
        endT .index=0
        # start of chain and end of transpose
        # off shell
        if(sind==iloc) :
            start.name="M"
            endT .name="M"
            start.value = vslashM .substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
            Symbols    += vslashMS.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
            endT.value  = vslashM2.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
            defns["vvP%s" % sind ] = ["vvP%s" % sind ,
                                      vslashD.substitute({ "var" : "Energy",
                                                    "v" :  "P%s" % sind })]
            dtemp[1]+=1
        # onshell
        else :
            start.name="S"
            endT .name="S"
            subs = {'s' : ("sbar%s" % sind)}
            start.value = sbar .substitute(subs)
            Symbols += sline.substitute(subs)
            subs = {'s' : ("s%s" % sind)}
            endT.value = spinor.substitute(subs)
            Symbols += sline.substitute(subs)
    # spin 3/2 fermion
    elif spins[sind-1]==4 :
        # check if we can easily contract
        contract=checkRSContract(parsed,sind,dtemp)
        # off-shell
        if(sind==iloc) :
            oindex = LorentzIndex(sind)
            oindex.type="O"
            unContracted[oindex]=0
            Symbols += vslashMS.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
            Symbols += momCom.substitute({"v" : "P%s" %sind })
            defns["vvP%s" % sind ] = ["vvP%s" % sind ,
                                      vslashD.substitute({ "var" : "Energy",
                                                           "v" :  "P%s" % sind })]
            dtemp[1] += 1
            if(contract=="") :
                rindex = LorentzIndex(sind)
                rindex.type="R"
                start=DiracMatrix()
                start.value = Template(rslash.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind, "loc" : sind }))
                start.name  = "RP"
                start.index = (oindex,rindex)
                endT=DiracMatrix()
                endT.value = Template(rslash2.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind, "loc" : sind }))
                endT.name = "RP"
                endT.index = (rindex,oindex)
            else :
                # construct dot product
                pi = LorentzIndex(sind)
                pi.type="P"
                pi.dimension=1
                (name,dummy) = constructDotProduct(pi,contract,defns)
                Symbols += momCom.substitute({"v" : contract })
                RB = vslash.substitute({ "v" : contract})
                Symbols += vslashS.substitute({ "v" : contract })
                Symbols += "%s = Symbol('%s')\n" % (name,name)
                defns["vv%s" % contract ] = ["vv%s" % contract,
                                             vslashD.substitute({ "var" : computeUnit(contract.dimension),
                                                                  "v" :  "%s" % contract })]
                start=DiracMatrix()
                start.name="RQ"
                start.index = oindex
                start.value =Template( rslashB.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind, "loc" : sind,
                                                            "DB" : RB, "dot" : name , "v2" : contract}))
                
                endT=DiracMatrix()
                endT.name = "RQ"
                endT.index = oindex
                endT.value = Template(rslash2B.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind, "loc" : sind ,
                                                            "DA" : RB, "dot" : name , "v2" : contract}))
        # on-shell
        else :
            start = DiracMatrix()
            endT  = DiracMatrix()
            # no contraction
            if(contract=="" or (contract.type=="E" and contract.value==iloc) ) :
                if contract == "" :
                    contract = LorentzIndex(sind)
                    contract.type="R"
                # start of matrix string
                start.value = Template(sbar  .substitute({'s' : ("Rsbar%s${L}" % sind)}))
                start.name="RS"
                start.index = contract
                # end of transpose string
                endT.value=Template(spinor.substitute({'s' : ("Rs%s${L}" % sind)}))
                endT.name="RS"
                endT.index = contract
                unContracted[contract]=0
                # variables for sympy
                for LI in imap :
                    Symbols += sline.substitute({'s' : ("Rs%s%s"    % (sind,LI))})
                    Symbols += sline.substitute({'s' : ("Rsbar%s%s" % (sind,LI))})
            else :
                # start of matrix string
                start.name="S"
                start.value = "Matrix([[%s,%s,%s,%s]])" % (RSDotProduct.substitute({'s' : ("Rsbar%s" % sind), 'v':contract, 'si' : 1}),
                                                           RSDotProduct.substitute({'s' : ("Rsbar%s" % sind), 'v':contract, 'si' : 2}),
                                                           RSDotProduct.substitute({'s' : ("Rsbar%s" % sind), 'v':contract, 'si' : 3}),
                                                           RSDotProduct.substitute({'s' : ("Rsbar%s" % sind), 'v':contract, 'si' : 4}))
                endT.name="S"
                endT.value = "Matrix([[%s],[%s],[%s],[%s]])" % (RSDotProduct.substitute({'s' : ("Rs%s" % sind), 'v':contract, 'si' : 1}),
                                                                RSDotProduct.substitute({'s' : ("Rs%s" % sind), 'v':contract, 'si' : 2}),
                                                                RSDotProduct.substitute({'s' : ("Rs%s" % sind), 'v':contract, 'si' : 3}),
                                                                RSDotProduct.substitute({'s' : ("Rs%s" % sind), 'v':contract, 'si' : 4}))
                Symbols += momCom.substitute({"v" : contract })
                for LI in ["x","y","z","t"] :
                    Symbols += sline.substitute({'s' : ("Rs%s%s"    % (sind,LI))})
                    Symbols += sline.substitute({'s' : ("Rsbar%s%s" % (sind,LI))})
    # last particle spin 1/2
    if( spins[lind-1]==2 ) :
        end    = DiracMatrix()
        startT = DiracMatrix()
        end   .index=0
        startT.index=0
        # end of chain
        if(lind==iloc) :
            end.name   ="M"
            startT.name="M"
            end.value    = vslashM2.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
            startT.value =  vslashM.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
            Symbols += vslashMS.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
            defns["vvP%s" % lind ] = ["vvP%s" % lind ,
                                      vslashD.substitute({ "var" : "Energy",
                                                           "v" :  "P%s" % lind })]
            dtemp[1] += 1
        else :
            startT.name="S"
            end   .name="S"
            subs = {'s' : ("s%s" % lind)}
            end.value = spinor.substitute(subs)
            Symbols += sline.substitute(subs)
            subs = {'s' : ("sbar%s" % lind)}
            startT.value = sbar .substitute(subs)
            Symbols += sline.substitute(subs)
    # last particle spin 3/2
    elif spins[lind-1]==4 :
        # check if we can easily contract
        contract=checkRSContract(parsed,lind,dtemp)
        # off-shell
        if(lind==iloc) :
            oindex = LorentzIndex(lind)
            oindex.type="O"
            unContracted[oindex]=0
            Symbols += vslashMS.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
            Symbols += momCom.substitute({"v" : "P%s" %lind })
            defns["vvP%s" % lind ] = ["vvP%s" % lind ,
                                      vslashD.substitute({ "var" : "Energy",
                                                           "v" :  "P%s" % lind })]
            dtemp[1] += 1
            if(contract=="") :
                rindex = LorentzIndex(lind)
                rindex.type="R"
                end=DiracMatrix()
                end.value = Template(rslash2.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind, "loc" : lind }))
                end.name  = "RP"
                end.index =   (rindex,oindex)
                startT=DiracMatrix()
                startT.value=Template(rslash.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind, "loc" : lind }))
                startT.name = "RP"
                startT.index = (oindex,rindex)
            else :
                # construct dot product
                pi = LorentzIndex(lind)
                pi.type="P"
                pi.dimension=1
                (name,unit) = constructDotProduct(pi,contract,defns)
                Symbols += momCom.substitute({"v" : contract })
                RB = vslash.substitute({ "v" : contract})
                Symbols += vslashS.substitute({ "v" : contract })
                Symbols += "%s = Symbol('%s')\n" % (name,name)
                defns["vv%s" % contract ] = ["vv%s" % contract,
                                             vslashD.substitute({ "var" : computeUnit(contract.dimension),
                                                                  "v" :  "%s" % contract })]
                end=DiracMatrix()
                end.name="RQ"
                end.index = oindex
                end.value =Template( rslash2B.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind, "loc" : lind,
                                                           "DA" : RB, "dot" : name , "v2" : contract}))
                
                startT=DiracMatrix()
                startT.name = "RQ"
                startT.index = oindex
                startT.value = Template(rslashB.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind, "loc" : lind ,
                                                             "DB" : RB, "dot" : name , "v2" : contract}))
        # on-shell
        else :
            end    = DiracMatrix()
            startT = DiracMatrix()
            # no contraction
            if(contract=="" or (contract.type=="E" and contract.value==iloc) ) :
                if contract == "" :
                    contract = LorentzIndex(lind)
                    contract.type="R"
                # end of matrix string
                end.value = Template(spinor.substitute({'s' : ("Rs%s${L}" % lind)}))
                end.name  = "RS"
                end.index = contract
                # start of matrix string
                startT.value = Template(sbar  .substitute({'s' : ("Rsbar%s${L}" % lind)}))
                startT.name  = "RS"
                startT.index = contract
                unContracted[contract]=0
                # variables for sympy
                for LI in imap :
                    Symbols += sline.substitute({'s' : ("Rs%s%s"    % (lind,LI))})
                    Symbols += sline.substitute({'s' : ("Rsbar%s%s" % (lind,LI))})
            # contraction
            else :
                # end of the matrix string
                end.name = "S"
                end.value = "Matrix([[%s],[%s],[%s],[%s]])" % (RSDotProduct.substitute({'s' : ("Rs%s" % lind), 'v':contract, 'si' : 1}),
                                                               RSDotProduct.substitute({'s' : ("Rs%s" % lind), 'v':contract, 'si' : 2}),
                                                               RSDotProduct.substitute({'s' : ("Rs%s" % lind), 'v':contract, 'si' : 3}),
                                                               RSDotProduct.substitute({'s' : ("Rs%s" % lind), 'v':contract, 'si' : 4}))
                startT.name = "S"
                startT.value = "Matrix([[%s,%s,%s,%s]])" % (RSDotProduct.substitute({'s' : ("Rsbar%s" % lind), 'v':contract, 'si' : 1}),
                                                            RSDotProduct.substitute({'s' : ("Rsbar%s" % lind), 'v':contract, 'si' : 2}),
                                                            RSDotProduct.substitute({'s' : ("Rsbar%s" % lind), 'v':contract, 'si' : 3}),
                                                            RSDotProduct.substitute({'s' : ("Rsbar%s" % lind), 'v':contract, 'si' : 4}))
                Symbols += momCom.substitute({"v" : contract })
                for LI in ["x","y","z","t"] :
                    Symbols += sline.substitute({'s' : ("Rs%s%s"    % (lind,LI))})
                    Symbols += sline.substitute({'s' : ("Rsbar%s%s" % (lind,LI))})
    return(sind,lind,expr,start,startT,end,endT,Symbols)

def calculateDirac(expr,start,end,startT,endT,sind,lind,Symbols,iloc) :
    res=[]
    for ichain in range(0,len(start)) :
        # calculate the matrix string
        etemp="*".join(str(x) for x in expr[ichain])
        temp={}
        exec("import sympy\nfrom sympy import Symbol,Matrix\n"+Symbols+"result="+
             ( "%s*%s*%s" %(start[ichain],etemp,end[ichain]))) in temp
        res.append(temp["result"])
        tempT={}
        exec("import sympy\nfrom sympy import Symbol,Matrix,Transpose\n"+Symbols+"result="+
             ( "%s*%s*Transpose(%s)*%s*%s" %(startT[ichain],CC,etemp,CD,endT[ichain]))) in tempT
        res.append(tempT["result"])
    if(len(start)==1) :
        if(iloc==0 or (iloc!=sind[ichain] and iloc!=lind[ichain])) :
            sVal = {'s' : temp ["result"][0,0],'sT' : tempT["result"][0,0]}
        else :
            sVal={}
            for jj in range(1,5) :
                sVal["s%s"  % jj] = temp ["result"][jj-1]
                sVal["sT%s" % jj] = tempT["result"][jj-1]
    else :
        sVal={}
        sVal["s"   ] = res[0][0,0]*res[2][0,0]
        sVal["sT2" ] = res[0][0,0]*res[3][0,0]
        sVal["sT1" ] = res[1][0,0]*res[2][0,0]
        sVal["sT12"] = res[1][0,0]*res[3][0,0]
    return sVal

def addToOutput(res,nchain,sign,rTemp) :
    # 1 spin chain
    if(nchain==1) :
        for ii in range(0,2) :
            if(rTemp[ii][0].shape[0]==1) :
                # result is scalar
                if(rTemp[ii][0].shape[1]==1) :
                    if(len(res[ii])==0) :
                        res[ii].append(sign*rTemp [ii][0][0,0])
                    else :
                        res[ii][0]  += sign*rTemp [ii][0][0,0]
                # result is a spinor
                elif(rTemp[ii][0].shape[1]==4) :
                    if(len(res[ii])==0) :
                        for j in range(0,4) :
                            res[ii].append(sign*rTemp[ii][0][0,j])
                    else :
                        for j in range(0,4) :
                            res[ii][j] += sign*rTemp[ii][0][0,j]
                else :
                    print("Size problem adding results A",sign,rTemp[ii].shape)
                    raise SkipThisVertex()
            # spinor
            elif(rTemp[ii][0].shape[0]==4 and rTemp[ii][0].shape[1]==1 ) :
                if(len(res[ii])==0) :
                    for j in range(0,4) :
                        res[ii].append(sign*rTemp[ii][0][j,0])
                else :
                    for j in range(0,4) :
                        res[ii][j] += sign*rTemp[ii][0][j,0]
            else :
                print("Size problem adding results A",sign,rTemp[ii][0].shape)
                raise SkipThisVertex()
    # 2 spin chains, should only be for a vertex
    else :
        for j1 in range(0,2) :
            for j2 in range (0,2) :
                val = sign*rTemp[j1][0]*rTemp[j2][1]
                if(len(res[3])==0) :
                    res[2*j1+j2].append(val[0,0])
                else :
                    res[2*j1+j2][0] += val[0,0]

def calculateDirac2(expr,start,end,startT,endT,sind,lind,Symbols,defns,
                    iloc,unContracted,spins,lorentz) :
    tDot=""
    # output
    sVal={}
    # no of chains
    nchain=len(expr)
    # now deal with the uncontracted cases
    contracted={}
    # sort out contracted and uncontracted indices
    keys = unContracted.keys()
    for key in keys:
        # summed dummy index
        if key.type=="D" :
            contracted[key]=0
            del unContracted[key]
        # RS index
        elif key.type =="R" :
            contracted[key]=0
            del unContracted[key]
        # tensor index
        elif key.type == "T1" or key.type=="T2" :
            contracted[key]=0
            del unContracted[key]
        # external index
        elif key.type == "O" :
            continue
        # uncontracted vector index
        elif key.type=="E" or key.type=="Q":
            continue
        else :
            print('Unknown type of uncontracted index',key)
            raise SkipThisVertex()
    # check the lorentz structures
    for lstruct in lorentz :
        if(lstruct.name=="Epsilon" or
           lstruct.name=="Vector") :
            for index in lstruct.lorentz :
                if(index.type=="E" and index.value==iloc) :
                    unContracted[index]=0
                elif(index.type=="P" or index.type=="E"
                     or index.type=="R" or index.type=="D") :
                    contracted[index]=0
                else :
                    print('Unknown index',index, 'in ',lstruct)
                    raise SkipThisVertex()
        elif(lstruct.name=="Tensor") :
            if(iloc==lstruct.value) :
                Symbols += momCom.substitute({"v": "P%s"%lstruct.value})
                Symbols += "OM%s = Symbol(\"OM%s\")\n" % (lstruct.value,lstruct.value) 
                Symbols += "M%s2 = Symbol(\"M%s2\")\n" % (lstruct.value,lstruct.value) 
                Symbols += "p2 = Symbol(\"p2\")\n"
                for ival in range(1,3) :
                    newIndex=LorentzIndex(ival)
                    newIndex.type="O"
                    newIndex.dimension=0
                    lstruct.lorentz.append(newIndex)
                    unContracted[newIndex]=0
                # contracted with self
                if(len(lstruct.lorentz)==0) :
                    pass
                # both indices uncontracted, deal with later
                elif lstruct.lorentz[0].type=="T1" and lstruct.lorentz[1].type=="T2":
                    pass
                elif lstruct.lorentz[0].type=="T1":
                    pIndex = LorentzIndex(lstruct.value)
                    pIndex.dimension=1
                    pIndex.type="P"
                    (tDot,dtemp) = constructDotProduct(pIndex,lstruct.lorentz[1],defns)
                    Symbols+="%s = Symbol(\"%s\")\n" %(tDot,tDot)
                    Symbols += momCom.substitute({"v": lstruct.lorentz[1]})
                elif lstruct.lorentz[1].type=="T2" :
                    pIndex = LorentzIndex(lstruct.value)
                    pIndex.dimension=1
                    pIndex.type="P"
                    (tDot,dtemp) = constructDotProduct(pIndex,lstruct.lorentz[0],defns)
                    Symbols+="%s = Symbol(\"%s\")\n" %(tDot,tDot)
                    Symbols += momCom.substitute({"v": lstruct.lorentz[0]})
            # both indices still to be contracted
            else :
                contracted[lstruct.lorentz[0].type]=0
                contracted[lstruct.lorentz[1].type]=0
        else :
            print('Unknown lorentz object in calculateDirac2',lstruct,iloc)
            raise SkipThisVertex()
    # iterate over the uncontracted indices
    while True :
        # loop over the unContracted indices
        res = []
        for i in range(0,nchain) :
            res.append([])
            res.append([])
        # loop over the contracted indices
        while True :
            # sign from metric tensor in contractions
            sign = 1
            for key,val in contracted.items() :
                if(val>0) : sign *=-1
            eTemp =[]
            sTemp =[]
            fTemp =[]
            sTTemp=[]
            fTTemp=[]
            # make the necessary replacements for remaining indices
            for ichain in range(0,nchain) :
                # compute the main expression
                eTemp.append([])
                for val in expr[ichain] :
                    # already a matrix
                    if(val.name=="M") :
                        eTemp[ichain].append(val)
                    # gamma(mu), replace with correct dirac matrix
                    elif(val.name=="GMU") :
                        if(val.index in contracted) :
                            eTemp[ichain].append(dirac[contracted[val.index]])
                        elif(val.index in unContracted) :
                            eTemp[ichain].append(dirac[unContracted[val.index]])
                        else :
                            print('Unknown index for gamma matrix',val)
                            raise SkipThisVertex()
                    # unknown to be sorted out
                    else :
                        print('Unknown type in expr',val)
                        raise SkipThisVertex()
                # start and end
                # start
                if(start[ichain].name=="S" or start[ichain].name=="M" ) :
                    sTemp.append(start[ichain].value)
                elif(start[ichain].name=="RS") :
                    if(start[ichain].index in contracted) :
                        sTemp.append(start[ichain].value.substitute({"L" : imap[  contracted[start[ichain].index]] }))
                    else :
                        sTemp.append(start[ichain].value.substitute({"L" : imap[unContracted[start[ichain].index]] }))
                elif(start[ichain].name=="RQ") :
                    i1 = unContracted[start[ichain].index]
                    sTemp.append(start[ichain].value.substitute({"A" : imap[i1], "DA" : dirac[i1] }))
                elif(start[ichain].name=="RP") :
                    i1 = unContracted[start[ichain].index[0]]
                    i2 = contracted[start[ichain].index[1]]
                    eta=0
                    if(i1==i2) :
                        if(i1==0) : eta =  1
                        else      : eta = -1
                    sTemp.append(start[ichain].value.substitute({"eta" : eta, "A":imap[i1] , "B":imap[i2] ,
                                                                 "DA": dirac[i1], "DB": dirac[i2]}))
                else :
                    print('Barred spinor not a spinor',start[ichain])
                    raise SkipThisVertex()
                if(startT[ichain].name=="S" or startT[ichain].name=="M" ) :
                    sTTemp.append(startT[ichain].value)
                elif(startT[ichain].name=="RS") :
                    if(startT[ichain].index in contracted) :
                        sTTemp.append(startT[ichain].value.substitute({"L" : imap[  contracted[startT[ichain].index]] }))
                    else :
                        sTTemp.append(startT[ichain].value.substitute({"L" : imap[unContracted[startT[ichain].index]] }))
                elif(startT[ichain].name=="RQ") :
                    i1 = unContracted[startT[ichain].index]
                    sTTemp.append(startT[ichain].value.substitute({"A" : imap[i1], "DA" : dirac[i1] }))
                elif(startT[ichain].name=="RP") :
                    i1 = unContracted[startT[ichain].index[0]]
                    i2 = contracted[startT[ichain].index[1]]
                    eta=0
                    if(i1==i2) :
                        if(i1==0) : eta =  1
                        else      : eta = -1
                    sTTemp.append(startT[ichain].value.substitute({"eta" : eta, "A":imap[i1] , "B":imap[i2] ,
                                                                   "DA": dirac[i1], "DB": dirac[i2]}))
                else :
                    print('barred spinorT not a spinor',startT[ichain])
                    raise SkipThisVertex()
                # end
                if(end[ichain].name=="S" or end[ichain].name=="M" ) :
                    fTemp.append(end[ichain].value)
                elif(end[ichain].name=="RS") :
                    if(end[ichain].index in contracted) :
                        fTemp.append(end[ichain].value.substitute({"L" : imap[  contracted[end[ichain].index]] }))
                    else :
                        fTemp.append(end[ichain].value.substitute({"L" : imap[unContracted[end[ichain].index]] }))
                elif(end[ichain].name=="RQ") :
                    i1 = unContracted[end[ichain].index]
                    fTemp.append(end[ichain].value.substitute({"B" : imap[i1], "DB": dirac[i1] }))
                elif(end[ichain].name=="RP") :
                    i1 =   contracted[end[ichain].index[0]]
                    i2 = unContracted[end[ichain].index[1]]
                    eta=0
                    if(i1==i2) :
                        if(i1==0) : eta =  1
                        else      : eta = -1
                    fTemp.append(end[ichain].value.substitute({"eta" : eta, "A":imap[i1] , "B":imap[i2] ,
                                                               "DA": dirac[i1], "DB": dirac[i2]}))
                else :
                    print('spinor not a spinor',end[ichain])
                    raise SkipThisVertex()
                if(endT[ichain].name=="S" or endT[ichain].name=="M" ) :
                    fTTemp.append(endT[ichain].value)
                elif(endT[ichain].name=="RS") :
                    if(endT[ichain].index in contracted) :
                        fTTemp.append(endT[ichain].value.substitute({"L" : imap[  contracted[endT[ichain].index]] }))
                    else :
                        fTTemp.append(endT[ichain].value.substitute({"L" : imap[unContracted[endT[ichain].index]] }))
                elif(endT[ichain].name=="RQ") :
                    i1 = unContracted[endT[ichain].index]
                    fTTemp.append(endT[ichain].value.substitute({"B" : imap[i1], "DB": dirac[i1] }))
                elif(endT[ichain].name=="RP") :
                    i1 =   contracted[endT[ichain].index[0]]
                    i2 = unContracted[endT[ichain].index[1]]
                    eta=0
                    if(i1==i2) :
                        if(i1==0) : eta =  1
                        else      : eta = -1
                    fTTemp.append(endT[ichain].value.substitute({"eta" : eta, "A":imap[i1] , "B":imap[i2] ,
                                                                 "DA": dirac[i1], "DB": dirac[i2]}))
                else :
                    print('spinorT not a spinor',endT[ichain])
                    raise SkipThisVertex()
            # and none dirac lorentz structures
            isZero = False
            for li in lorentz:
                # uncontracted vector
                if(li.name=="Vector") :
                    index = unContracted[li.lorentz[0]]
                    Symbols += momCom.substitute({"v":li.value})
                    for ichain in range(0,nchain) :
                        eTemp[ichain].append("%s%s"% (li.value,imap[index]) )
                elif(li.name=="Epsilon") :
                    value=""
                    ival=[]
                    for index in li.lorentz :
                        if(index in contracted) :
                            if(index.type=="P" or index.type=="E") :
                                value += "*%s%s" % (index,imap[contracted[index]])
                                ival.append(contracted[index])
                            elif(index.type=="R" or index.type=="D") :
                                ival.append(contracted[index])
                            else :
                                print('Unknown index in Epsilon Tensor',index)
                                raise SkipThisVertex()
                        elif(index in unContracted) :
                            ival.append(unContracted[index])
                    if(len(value)!=0 and value[0]=="*") :
                        value = value[1:]
                    eVal = epsValue[ival[0]][ival[1]][ival[2]][ival[3]]
                    if(eVal==0) :
                        isZero = True
                    else :
                        for ichain in range(0,nchain) :
                            eTemp[ichain].append("(%s*%s)"% (eVal,value) )
                elif(li.name=="Tensor") :
                    if(li.lorentz[0] in unContracted and li.lorentz[1] in unContracted) :
                        value='0.5*(%s)'%tPropA[unContracted[li.lorentz[0]]][unContracted[li.lorentz[0]]].substitute({"iloc" : li.value})
                        for ichain in range(0,nchain) :
                            eTemp[ichain].append("(%s)"% (value) )
                    elif(len(li.lorentz)==4) :
                        if li.lorentz[0].type=="T1" and li.lorentz[1].type=="T2" :
                            value='0.5*(%s)'%tPropC[unContracted[li.lorentz[2]]][unContracted[li.lorentz[3]]][contracted[li.lorentz[0]]][contracted[li.lorentz[1]]].substitute({"iloc" : li.value})
                        elif li.lorentz[0].type=="T1":
                            value='0.5*(%s)'%tPropB[unContracted[li.lorentz[2]]][unContracted[li.lorentz[3]]][contracted[li.lorentz[0]]].substitute({"iloc" : li.value,
                                                                                                                                          "V" : li.lorentz[1],
                                                                                                                                          "dot" : tDot})
                        elif li.lorentz[1].type=="T2" :
                            value= '0.5*(%s)'%tPropB[unContracted[li.lorentz[2]]][unContracted[li.lorentz[3]]][contracted[li.lorentz[1]]].substitute({"iloc" : li.value,
                                                                                                                                          "V" : li.lorentz[0],
                                                                                                                                          "dot" : tDot})
                        else :
                            print('Both contracted tensor indices contracted')
                            raise SkipThisVertex()
                        for ichain in range(0,nchain) :
                            eTemp[ichain].append("(%s)"% (value) )
                    else :
                        print('Uncontracted on-shell tensor')
                        raise SkipThisVertex()
                # unknown
                else :
                    print('Unknown expression in lorentz loop',li)
                    raise SkipThisVertex()
            # now evaluate the result
            if(not isZero) :
                rTemp =[[],[]]
                for ichain in range(0,nchain) :
                    core = "*".join(str(x) for x in eTemp[ichain])
                    temp={}
                    exec("import sympy\nfrom sympy import Symbol,Matrix\n"+Symbols+"result="+
                         ( "(%s)*(%s)*(%s)" %(sTemp[ichain],core,fTemp[ichain]))) in temp
                    rTemp[0].append(temp["result"])
                    temp={}
                    exec("import sympy\nfrom sympy import Symbol,Matrix,Transpose\n"+Symbols+"result="+
                         ( "(%s)*(%s)*(Transpose(%s))*(%s)*(%s)" %(sTTemp[ichain],CC,core,CD,fTTemp[ichain]))) in temp
                    rTemp[1].append(temp["result"])
                # and add it to the output
                addToOutput(res,nchain,sign,rTemp)
            #### END OF THE CONTRACTED LOOP #####
            # increment the indices being summed over
            keys=contracted.keys()
            ii = len(keys)-1
            while ii >=0 :
                if(contracted[keys[ii]]<3) :
                    contracted[keys[ii]]+=1
                    break
                else :
                    contracted[keys[ii]]=0
                    ii-=1
            nZero=0
            for (key,val) in contracted.items() :
                if(val==0) : nZero+=1
            if(nZero==len(contracted)) : break
        ###### END OF THE UNCONTRACTED LOOP ######
        # no uncontracted indices
        if(len(unContracted)==0) :
            if(len(res[0])==1) :
                # scalar
               if(len(res)==2) :
                   sVal["s" ] = res[0]
                   sVal["sT"] = res[1]
                # 4 fermion
               else :
                   sVal["s"   ] = res[0]
                   sVal["sT2" ] = res[1]
                   sVal["sT1" ] = res[2]
                   sVal["sT12"] = res[3]
            # spinor
            elif(len(res[0])==4) :
                for k in range(0,4) :
                    sVal[ "s%s"  % (k+1) ] = res[0][k]
                    sVal[ "sT%s" % (k+1) ] = res[1][k]
            else :
                print('Sum problem',len(res),len(res[0]))
                raise SkipThisVertex()
            break
        # uncontracted indices
        else :
            istring = ""
            for (key,val) in unContracted.items() :
                istring +=imap[val]
            if(len(istring)>2) :
                print('Index problem',istring)
                raise SkipThisVertex()
            sVal[istring]     = res[0]
            sVal[istring+"T"] = res[1]
            # increment the unsummed indices
            keys=unContracted.keys()
            ii = len(keys)-1
            while ii >=0 :
                if(unContracted[keys[ii]]<3) :
                    unContracted[keys[ii]]+=1
                    break
                else :
                    unContracted[keys[ii]]=0
                    ii-=1
            nZero=0
            for (key,val) in unContracted.items() :
                if(val==0) : nZero+=1
            if(nZero==len(unContracted)) : break
    # handle the vector case
    if( "tt" in sVal) :
        if(len(sVal)==32 and "tt" in sVal and len(sVal["tt"])==1) :
            for key in sVal:
                sVal[key] = sVal[key][0]
        else :
            print('Tensor sum problem',len(sVal))
            raise SkipThisVertex()
    elif( "t" in sVal ) :
        # deal with pure vectors
        if(len(sVal)==8 and "t" in sVal and len(sVal["t"])==1) :
            pass
        # RS spinors
        elif(len(sVal)==8 and "t" in sVal and len(sVal["t"])==4) :
            pass
        else :
            print('Value problem',len(sVal))
            raise SkipThisVertex()
    else :
        if("s" in sVal) :
            for key in sVal:
                sVal[key] = sVal[key][0]
    return sVal
            
def convertDiracStructure(parsed,output,dimension,defns,iloc,L,lorentztag,vertex) :
    # get the spins of the particles
    spins = vertex.lorentz[0].spins
    # check if we have one or two spin chains
    nchain = (lorentztag.count("F")+lorentztag.count("R"))/2
    # storage of the intermediate results
    expr  =[]
    start =[]
    end   =[]
    startT=[]
    endT  =[]
    sind=[0]*nchain
    lind=[0]*nchain
    unContracted={}
    Symbols=""
    dtemp=[0,0,0]
    # parse the dirac matrix strings
    for ichain in range(0,nchain) :
        expr  .append([])
        start .append("")
        startT.append("")
        end   .append("")
        endT  .append("")
        sind[ichain],lind[ichain],expr[ichain],start[ichain],startT[ichain],end[ichain],endT[ichain],Symbols =\
            processChain(dtemp,parsed,spins,Symbols,unContracted,defns,iloc)
    # standard ordering of the chains
    for ichain in range(0,nchain) :
        for jchain in range(ichain+1,nchain) :
            if(sind[jchain]<sind[ichain]) :
                sind[ichain]  ,sind[jchain]   = sind[jchain]  ,sind[ichain]
                lind[ichain]  ,lind[jchain]   = lind[jchain]  ,lind[ichain]
                expr[ichain]  ,expr[jchain]   = expr[jchain]  ,expr[ichain]
                start[ichain] ,start[jchain]  = start[jchain] ,start[ichain]
                startT[ichain],startT[jchain] = startT[jchain],startT[ichain]
                end[ichain]   ,end[jchain]    = end[jchain]   ,end[ichain]
                endT[ichain]  ,endT[jchain]   = endT[jchain]  ,endT[ichain]
    # clean up parsed
    # check we've dealt with everything
    parsed = [x for x in parsed if x != ""]
    lorentz=[]
    if(len(parsed)!=0) :
        for i in range(0,len(parsed)) :
            if(parsed[i].name=="Metric") :
                found = False
                for ll in parsed[i].lorentz :
                    if(ll.type=="E" and ll.value==iloc) :
                        found = True
                        un=ll
                    else :
                        lo=ll
                if(found) :
                    lstruct = LorentzStructure()
                    lstruct.name="Vector"
                    lstruct.value=lo
                    lstruct.lorentz=[un]
                    lstruct.spin=[]
                    lorentz.append(lstruct)
                    parsed[i]=""
                    unContracted[un]=0
                    dimension[2] += lo.dimension
            elif(parsed[i].name=="Epsilon") :
                lstruct = LorentzStructure()
                lstruct.name="Epsilon"
                lstruct.lorentz=parsed[i].lorentz
                lstruct.value=0
                lstruct.spin=[]
                for index in lstruct.lorentz:
                    if(index.type=="P") :
                        dimension[2]+=1
                        if( not index in defns) :
                            defns[(index,)]=["",""]
                    if(index.type=="P" or
                       (index.type=="E" and index.value!=iloc)) :
                        Symbols += momCom.substitute( {"v": index} )
                lorentz.append(lstruct)
                parsed[i]=""
            elif(parsed[i].name=="Tensor") :
                lstruct = LorentzStructure()
                lstruct.name="Tensor"
                lstruct.value=parsed[i].value
                lstruct.lorentz=parsed[i].lorentz
                lstruct.spin=[]
                for index in parsed[i].lorentz:
                    dimension[2] += index.dimension
                parsed[i]=""
                lorentz.append(lstruct)
            else :
                print('Unknown lorentz structure',parsed[i])
                raise SkipThisVertex()
                
        parsed = [x for x in parsed if x != ""]
        if(len(parsed)!=0) :
            print("Can't parse ",parsed,iloc)
            raise SkipThisVertex()
    sVal ={}
    dimension = list(map(lambda x, y: x + y, dtemp, dimension))
    # deal with the simplest case first
    if len(unContracted) == 0 and len(lorentz) == 0:
        sVal = calculateDirac(expr,start,end,startT,endT,sind,lind,Symbols,iloc)
    else :
        sVal = calculateDirac2(expr,start,end,startT,endT,sind,lind,Symbols,defns,
                               iloc,unContracted,spins,lorentz)
    # set up the output
    old = output
    if(nchain==1) :
        output = [old,sVal,(lind[0],sind[0])]
    else :
        output = [old,sVal,(lind[0],sind[0],lind[1],sind[1])]
    return (output,dimension,defns)
            
def convertLorentz(Lstruct,lorentztag,order,vertex,iloc,defns,evalVertex) :
    eps = False
    # split the structure into individual terms
    structures=Lstruct.structure.split()
    parsed=[]
    # parse structures and convert lorentz contractions to dot products
    for struct in structures :
        ptemp = parse_structure(struct,Lstruct.spins)
        parsed.append(contract(ptemp))
    # now in a position to generate the code
    vals=generateVertex(iloc,Lstruct,parsed,lorentztag,vertex,defns)
    evalVertex.append((vals[0],vals[1]))
    if(vals[2]) : eps=True
    return eps

def swapOrder(vertex,iloc,momenta,fIndex) :
    names=['','sca','sp','v']
    waves=['','sca',''  ,'E']
    output=""
    for i in range(1,4) :
        ns = vertex.lorentz[0].spins.count(i)
        if((ns<=1 and i!=2) or (ns<=2 and i==2)) : continue
        if(i!=3 and i!=1) :
            print('Swap problem',i)
            raise SkipThisVertex()
        sloc=[]
        for j in range(0,len(vertex.lorentz[0].spins)) :
            if(vertex.lorentz[0].spins[j]==i) : sloc.append(j+1)
        if iloc in sloc : sloc.remove(iloc)
        if(len(sloc)==1) : continue
        for j in range(0,len(sloc)) :
            output += "    long id%s = %sW%s.id();\n" % (sloc[j],names[i],sloc[j])
        for j in range(0,len(sloc)) :
            for k in range(j+1,len(sloc)) :
                code = vertex.particles[sloc[j]-1].pdg_code
                output += "    if(id%s!=%s) {\n" % (sloc[j],code)
                output += "        swap(id%s,id%s);\n" % (sloc[j],sloc[k])
                output += "        swap(%s%s,%s%s);\n" % (waves[i],sloc[j],waves[i],sloc[k])
                if(momenta[sloc[j]-1][0] or momenta[sloc[k]-1][0]) :
                    momenta[sloc[j]-1][0] = True
                    momenta[sloc[k]-1][0] = True
                    output += "        swap(P%s,P%s);\n" % (sloc[j],sloc[k])
                output += "    };\n"
    return output

def swapOrderFFFF(vertex,iloc,fIndex) :
    output=""
    for j in range(0,len(fIndex)) :
        if(j%2==0) :
            output += "    SpinorWaveFunction s%s = sW%s;\n" % (fIndex[j],fIndex[j])
        else :
            output += "    SpinorBarWaveFunction sbar%s = sbarW%s;\n" % (fIndex[j],fIndex[j])

    for j in range(0,len(fIndex)) :
        if(j%2==0) :
            output += "    long id%s = sW%s.id();\n" % (fIndex[j],fIndex[j])
        else :
            output += "    long id%s = sbarW%s.id();\n" % (fIndex[j],fIndex[j])

    for j in range(0,2) :
        code = vertex.particles[fIndex[j]-1].pdg_code
        output += "    if(id%s!=%s) {\n" % (fIndex[j],code)
        output += "        swap(id%s,id%s);\n" % (fIndex[j],fIndex[j+2])
        wave="s"
        if(j==1) : wave = "sbar"
        output += "        swap(%s%s,%s%s);\n" % (wave,fIndex[j],wave,fIndex[j+2])
        output += "    };\n"
    return output

def constructSignature(vertex,order,iloc,decls,momenta,waves,fIndex) :
    nf=0
    poff=""
    offType="Complex"
    for i in order : 
        spin = vertex.lorentz[0].spins[i-1]
        if(i==iloc) :
            if(spin==1) :
                offType="ScalarWaveFunction"
            elif(spin==2) :
                if(i%2==1) :
                    offType="SpinorBarWaveFunction"
                else :
                    offType="SpinorWaveFunction"
                nf+=1
            elif(spin==4) :
                if(i%2==1) :
                    offType="RSSpinorBarWaveFunction"
                else :
                    offType="RSSpinorWaveFunction"
                nf+=1
            elif(spin==3) :
                offType="VectorWaveFunction"
            elif(spin==5) :
                offType="TensorWaveFunction"
            else :
                print('Unknown spin',spin)
                raise SkipThisVertex()
            momenta.append([False,""])
        else :
            if(spin==1) :
                decls.append("ScalarWaveFunction & scaW%s" % (i))
                momenta.append([False,"Lorentz5Momentum P%s =-scaW%s.momentum();" % (i,i)])
                waves.append("Complex sca%s = scaW%s.wave();" % (i,i))
            elif(spin==2) :
                if(i%2==1) :
                    decls.append("SpinorWaveFunction & sW%s" % (fIndex[i-1]))
                    momenta.append([False,"Lorentz5Momentum P%s =-sW%s.momentum();" % (fIndex[i-1],fIndex[i-1])])
                    waves.append("LorentzSpinor<double> s%s = sW%s.wave();" % (fIndex[i-1],fIndex[i-1]))
                    nf+=1
                else :
                    decls.append("SpinorBarWaveFunction & sbarW%s" % (fIndex[i-1]))
                    momenta.append([False,"Lorentz5Momentum P%s =-sbarW%s.momentum();" % (fIndex[i-1],fIndex[i-1])])
                    waves.append("LorentzSpinorBar<double> sbar%s = sbarW%s.wave();" % (fIndex[i-1],fIndex[i-1]))
                    nf+=1
            elif(spin==3) :
                decls.append("VectorWaveFunction & vW%s" % (i))
                momenta.append([False,"Lorentz5Momentum P%s =-vW%s.momentum();" % (i,i)])
                waves.append("LorentzPolarizationVector E%s = vW%s.wave();" % (i,i))
            elif(spin==4) :
                if(i%2==1) :
                    decls.append("RSSpinorWaveFunction & RsW%s" % (i))
                    momenta.append([False,"Lorentz5Momentum P%s =-RsW%s.momentum();" % (i,i)])
                    waves.append("LorentzRSSpinor<double> Rs%s = RsW%s.wave();" % (i,i))
                    nf+=1
                else :
                    decls.append("RSSpinorBarWaveFunction & RsbarW%s" % (i))
                    momenta.append([False,"Lorentz5Momentum P%s =-RsbarW%s.momentum();" % (i,i)])
                    waves.append("LorentzRSSpinorBar<double> Rsbar%s = RsbarW%s.wave();" % (i,i))
                    nf+=1
            elif(spin==5) :
                decls.append("TensorWaveFunction & tW%s" % (i))
                momenta.append([False,"Lorentz5Momentum P%s =-tW%s.momentum();" % (i,i)])
                waves.append("LorentzTensor<double> T%s = tW%s.wave();" % (i,i))
            else :
                print('Unknown spin',spin)
                raise SkipThisVertex()
            poff += "-P%s" % (i)
    # ensure unbarred spinor first
    ibar=-1
    isp=-1
    for i in range(0,len(decls)) :
        if(decls[i].find("Bar")>0 and ibar==-1) :
            ibar=i
        elif(decls[i].find("Spinor")>=0 and isp==-1) :
            isp=i
    if(isp!=-1 and ibar!=-1 and isp>ibar) :
        decls[ibar],decls[isp] = decls[isp],decls[ibar]
    # constrct the signature
    poff = ("Lorentz5Momentum P%s = " % iloc ) + poff
    sig=""
    if(iloc==0) :
        sig="%s evaluate(Energy2, const %s)" % (offType,", const ".join(decls))
        # special for VVT vertex
        if(len(vertex.lorentz[0].spins)==3 and vertex.lorentz[0].spins.count(3)==2 and
           vertex.lorentz[0].spins.count(5)==1) :
            sig = sig[0:-1] + ", Energy vmass=-GeV)"
    else :
        sig="%s evaluate(Energy2, int iopt, tcPDPtr out, const %s, complex<Energy> mass=-GeV, complex<Energy> width=-GeV)" % (offType,", const ".join(decls))
        momenta.append([True,poff+";"])
        # special for VVT vertex
        if(len(vertex.lorentz[0].spins)==3 and vertex.lorentz[0].spins.count(3)==2 and
           vertex.lorentz[0].spins.count(5)==1 and vertex.lorentz[0].spins[iloc-1]==5) :
            sig=sig.replace("complex<Energy> mass=-GeV","Energy vmass=-GeV, complex<Energy> mass=-GeV")
        for i in range(0,len(momenta)) : momenta[i][0]=True
    return offType,nf,poff,sig

def combineResult(res,nf,ispin,vertex) :
    # extract the vals and dimensions
    (vals,dim) = res
    # construct the output structure
    # vertex and off-shell scalars
    if(ispin<=1) :
        otype={'res':""}
    # spins
    elif(ispin==2) :
        otype={'s1':"",'s2':"",'s3':"",'s4':""}
    # vectors
    elif(ispin==3) :
        if( "t" in vals[0][1] ) :
            otype={'t':"",'x':"",'y':"",'z':""}
        else :
            otype={"res":""}
    # RS spinors
    elif(ispin==4) :
        otype={}
        for i1 in imap :
            for i in range(1,5) :
                otype["%ss%s"% (i1,i)]=""
    # off-shell tensors
    elif(ispin==5) :
        otype={}
        for i1 in imap :
            for i2 in imap :
                otype["%s%s"%(i1,i2)]=""
    else :      
        print('Unknown spin',ispin)
        raise SkipThisVertex()
    expr=[otype]
    for i in range(0,nf-1) :
        expr.append(copy.copy(otype))
    # dimension for checking
    dimCheck=dim[0]
    for i in range(0,len(vals)) :
        # simple signs
        if(vals[i]=="+" or vals[i]=="-") :
            for ii in range(0,len(expr)) :
                for(key,val) in expr[ii].items() :
                    expr[ii][key] = expr[ii][key]+vals[i]
            continue
        # check the dimensions
        if(dimCheck[0]!=dim[i][0] or dimCheck[1]!=dim[i][1] or
           dimCheck[2]!=dim[i][2]) :
            print("Dimension problem in result",i,dimCheck,dim,vertex)
            print(vertex.lorentz)
            for j in range(0,len(vals)) :
                print(j,dim[j],vals[j])
            raise SkipThisVertex()
        # simplest case 
        if(isinstance(vals[i], basestring)) :
            for ii in range(0,len(expr)) :
                for(key,val) in expr[ii].items() :
                    expr[ii][key] = expr[ii][key]+vals[i]
            continue
        # more complex structures
        pre = vals[i][0]
        if(pre=="(1.0)") : pre=""
        if(not isinstance(vals[i][1],dict)) :
            print('must be a dict here')
            raise SkipThisVertex()
        # tensors
        if("tt" in vals[i][1]) :
            for i1 in imap :
                for i2 in imap :
                    key="%s%s"%(i1,i2)
                    if(pre=="") :
                        expr[0][key] += "(%s)" % vals[i][1][key]
                    else :
                        expr[0][key] += "%s*(%s)" % (pre,vals[i][1][key])
                    if(len(expr)==2) :
                        if(pre=="") :
                            expr[1][key] +="(%s)" % vals[i][1][key+"T"]
                        else :
                            expr[1][key] +="%s*(%s)" % (pre,vals[i][1][key+"T"])
        # standard fermion vertex case
        elif(len(vals[i][1])==2 and "s" in vals[i][1] and "sT" in vals[i][1]) :
            if(pre=="") :
                expr[0]["res"] += "(%s)" % vals[i][1]["s"]
                expr[1]["res"] += "(%s)" % vals[i][1]["sT"]
            else :
                expr[0]["res"] += "%s*(%s)" % (pre,vals[i][1]["s"])
                expr[1]["res"] += "%s*(%s)" % (pre,vals[i][1]["sT"])
        # spinor case
        elif(len(vals[i][1])==8 and "s1" in vals[i][1]) :
            for jj in range(1,5) :
                if(pre=="") :
                    expr[0]["s%s" % jj] += "(%s)" % vals[i][1]["s%s"  % jj]
                    expr[1]["s%s" % jj] += "(%s)" % vals[i][1]["sT%s" % jj]
                else :
                    expr[0]["s%s" % jj] += "%s*(%s)" % (pre,vals[i][1]["s%s"  % jj])
                    expr[1]["s%s" % jj] += "%s*(%s)" % (pre,vals[i][1]["sT%s" % jj])
        # vector
        elif(len(vals[i][1])%4==0 and "t" in vals[i][1] and len(vals[i][1]["t"])==1 ) :
            for i1 in imap :
                if(pre=="") :
                    expr[0][i1] += "(%s)" % vals[i][1][i1][0]
                else :
                    expr[0][i1] += "%s*(%s)" % (pre,vals[i][1][i1][0])
                if(len(expr)==2) :
                    if(pre=="") :
                        expr[1][i1] +="(%s)" % vals[i][1][i1+"T"][0]
                    else :
                        expr[1][i1] +="%s*(%s)" % (pre,vals[i][1][i1+"T"][0])
        # 4 fermion vertex case
        elif(len(vals[i][1])==4 and "sT12" in vals[i][1]) :
            if(pre=="") :
                expr[0]["res"] += "(%s)" % vals[i][1]["s"]
                expr[1]["res"] += "(%s)" % vals[i][1]["sT2"]
                expr[2]["res"] += "(%s)" % vals[i][1]["sT1"]
                expr[3]["res"] += "(%s)" % vals[i][1]["sT12"]
            else :
                expr[0]["res"] += "%s*(%s)" % (pre,vals[i][1]["s"])
                expr[1]["res"] += "%s*(%s)" % (pre,vals[i][1]["sT2"])
                expr[2]["res"] += "(%s)" % vals[i][1]["sT1"]
                expr[3]["res"] += "(%s)" % vals[i][1]["sT12"]
        # RS spinor
        elif(len(vals[i][1])%4==0 and "t" in vals[i][1] and len(vals[i][1]["t"])==4 ) :
            for i1 in imap :
                for k in range(1,5) :
                    key = "%ss%s" % (i1,k)
                    if(pre=="") :
                        expr[0][key] += "(%s)" % vals[i][1][i1][k-1]
                        expr[1][key] += "(%s)" % vals[i][1][i1+"T"][k-1]
                    else :
                        expr[0][key] += "%s*(%s)" % (pre,vals[i][1][i1][k-1])
                        expr[1][key] += "%s*(%s)" % (pre,vals[i][1][i1+"T"][k-1])
        else :
            print('problem with type',vals[i])
            raise SkipThisVertex()
    # no of particles in the vertex
    vDim = len(vertex.lorentz[0].spins)
    # compute the unit and apply it
    unit = computeUnit2(dimCheck,vDim)
    if(unit!="") :
        for ii in range(0,len(expr)) :
            for (key,val) in expr[ii].items() :
                expr[ii][key] = "(%s)*(%s)" % (val,unit)
    return expr

def headerReplace(inval) :
    return inval.replace("virtual","").replace("ScalarWaveFunction","").replace("SpinorWaveFunction","") \
                .replace("SpinorBarWaveFunction","").replace("VectorWaveFunction","").replace("TensorWaveFunction","") \
                .replace("Energy2","q2").replace("int","").replace("complex<Energy>","").replace("Energy","").replace("=-GeV","") \
                .replace("const  &","").replace("tcPDPtr","").replace("  "," ").replace("Complex","")

def combineComponents(result,offType,RS) :
    for i in range(0,len(result)) :
        for (key,value) in result[i].items() :
            output=py2cpp(value.strip(),True)
            result[i][key]=output[0]
    # simplest case, just a value
    if(len(result[0])==1 and "res" in result[0]) :
        for i in range(0,len(result)) :
            result[i] = result[i]["res"]
            result[i]=result[i].replace("1j","ii")
        return
    # calculate the substitutions
    if(not isinstance(result[0],basestring)) :
        subs=[]
        for ii in range(0,len(result)) :
            subs.append({})
            for (key,val) in result[ii].items() :
                subs[ii]["out%s" % key]= val
    # spinors
    if("s1" in result[0]) :
        stype  = "LorentzSpinor"
        sbtype = "LorentzSpinorBar"
        if(offType.find("Bar")>0) : (stype,sbtype)=(sbtype,stype)
        subs[0]["type"] = stype
        result[0]  = SpinorTemplate.substitute(subs[0])
        subs[1]["type"] = sbtype
        result[1]  = SpinorTemplate.substitute(subs[1])
    # tensors
    elif("tt" in result[0]) :
        for ii in range(0,len(result)) :
            result[ii] = Template("LorentzTensor<double>(${outxx},\n${outxy},\n${outxz},\n${outxt},\n${outyx},\n${outyy},\n${outyz},\n${outyt},\n${outzx},\n${outzy},\n${outzz},\n${outzt},\n${outtx},\n${outty},\n${outtz},\n${outtt})").substitute(subs[ii])
            result[ii]=result[ii].replace("(+","(")
    # vectors
    elif("t" in result[0]) :
        for ii in range(0,len(result)) :
            result[ii] = Template("LorentzVector<Complex>(${outx},\n${outy},\n${outz},\n${outt})").substitute(subs[ii])
            result[ii]=result[ii].replace("(+","(")
    # RS spinors
    elif("ts1" in result[0]) :
        stype  = "LorentzRSSpinor"
        sbtype = "LorentzRSSpinorBar"
        if(offType.find("Bar")>0) : (stype,sbtype)=(sbtype,stype)
        subs[0]["type"] = stype      
        result[0] = RSSpinorTemplate.substitute(subs[0])
        subs[1]["type"] = sbtype
        result[1] = RSSpinorTemplate.substitute(subs[1])
    else :
        print('Type not implemented',result)
        raise SkipThisVertex()
    for i in range(0,len(result)) :
        result[i]=result[i].replace("1j","ii")

def generateEvaluateFunction(model,vertex,iloc,values,defns,vertexEval,cf,order) :
    RS = "R" in vertex.lorentz[0].name
    FM = "F" in vertex.lorentz[0].name
    # extract the start and end of the spin chains
    if( RS or FM ) :
        fIndex = vertexEval[0][0][0][2]
    else :
        fIndex=0
    # first construct the signature of the function
    decls=[]
    momenta=[]
    waves=[]
    offType,nf,poff,sig = constructSignature(vertex,order,iloc,decls,momenta,waves,fIndex)
    # combine the different terms in the result
    symbols=set()
    localCouplings=[]
    result=[]
    ispin = 0
    if(iloc!=0) :
        ispin = vertex.lorentz[0].spins[iloc-1]
    # put the lorentz structures and couplings together
    for j in range(0,len(vertexEval)) :
        # get the lorentz structure piece
        expr = combineResult(vertexEval[j],nf,ispin,vertex)
        # get the coupling for this bit
        val, sym = py2cpp(values[j])
        localCouplings.append("Complex local_C%s = %s;\n" % (j,val))
        symbols |=sym
        # put them together
        vtype="Complex"
        if("res" in expr[0] and offType=="VectorWaveFunction") :
            vtype="LorentzPolarizationVector"
        if(len(result)==0) :
            for ii in range(0,len(expr)) :
                result.append({})
                for (key,val) in expr[ii].items() :
                    result[ii][key] = " %s((local_C%s)*(%s)) " % (vtype,j,val)
        else :
            for ii in range(0,len(expr)) :
                for (key,val) in expr[ii].items(): 
                    result[ii][key] += " + %s((local_C%s)*(%s)) " % (vtype,j,val)
    # for more complex types merge the spin/lorentz components
    combineComponents(result,offType,RS)
    # multiple by scalar wavefunctions
    scalars=""
    for i in range (0,len(vertex.lorentz[0].spins)) :
        if(vertex.lorentz[0].spins[i]==1 and i+1!=iloc) :
            scalars += "sca%s*" % (i+1)
    if(scalars!="") :
        for ii in range(0,len(result)) :
            result[ii]  = "(%s)*(%s)" % (result[ii],scalars[0:-1])
    # vertex, just return the answer
    if(iloc==0) :
        result[0] = "return (%s)*(%s);\n" % (result[0],py2cpp(cf)[0])
        if(FM or RS) :
            for i in range(1,len(result)) :
                result[i] = "return (%s)*(%s);\n" % (result[i],py2cpp(cf)[0])
    # off-shell particle
    else :
        # off-shell scalar
        if(vertex.lorentz[0].spins[iloc-1] == 1 ) :
            result[0] = scaTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],res=result[0])
            if(FM or RS) :
                result[1] = scaTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],res=result[1])
        # off-shell fermion
        elif(vertex.lorentz[0].spins[iloc-1] == 2 ) :
            result[0] = sTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offType.replace("WaveFunction",""),
                                         res=result[0].replace( "M%s" % iloc, "mass" ),offTypeB=offType)
            if(FM or RS) :
                if(offType.find("Bar")>0) :
                    offTypeT=offType.replace("Bar","")
                else :
                    offTypeT=offType.replace("Spinor","SpinorBar")
                result[1] = sTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offTypeT.replace("WaveFunction",""),
                                             res=result[1].replace( "M%s" % iloc, "mass" ),offTypeB=offTypeT)
        # off-shell vector
        elif(vertex.lorentz[0].spins[iloc-1] == 3 ) :
            result[0] = vecTemplate.format(iloc=iloc,res=result[0],cf=py2cpp(cf)[0])
            if(FM or RS) :
                result[1] = vecTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],res=result[1])
        elif(vertex.lorentz[0].spins[iloc-1] == 4 ) :
            if(offType.find("Bar")>0) :
                offTypeT=offType.replace("Bar","")
            else :
                offTypeT=offType.replace("Spinor","SpinorBar")
            result[1] = RSTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offTypeT.replace("WaveFunction",""),
                                          res=result[1].replace( "M%s" % iloc, "mass" ),offTypeB=offTypeT)
            result[0] = RSTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offType.replace("WaveFunction",""),
                                          res=result[0].replace( "M%s" % iloc, "mass" ),offTypeB=offType)
        # tensors
        elif(vertex.lorentz[0].spins[iloc-1]) :
            if(RS) :
                print("RS spinors and tensors not handled")
                raise SkipThisVertex()
            result[0] = tenTemplate.format(iloc=iloc,cf=py2cpp(cf)[0],res=result[0])
            if(FM or RS) :
                result[1] = tenTemplate.format(iloc=iloc,cf=py2cpp(cf)[0],res=result[1])
        else :
            print('Unknown spin for off-shell particle',vertex.lorentz[0].spins[iloc-1])
            raise SkipThisVertex()
    # check if momenta defns needed to clean up compile of code
    for (key,val) in defns.items() :
        if( isinstance(key, basestring)) :
            if(key.find("vvP")==0) :
                momenta[int(key[3])-1][0] = True
        else :
            for vals in key :
                if(vals.type=="P") :
                    momenta[vals.value-1][0] = True
    # cat the definitions
    defString=""
    for (key,value) in defns.items() :
        if(value[0]=="") : continue
        if(value[0][0]=="V" or value[0][0]=="T") :
            defString+="    %s\n" %value[1]
    for (key,value) in defns.items() :
        if(value[0]=="") : continue
        if(value[0][0]!="V" and  value[0][0]!="T") :
            defString+="    %s\n" %value[1]
    if(len(result)<=2) :
        sorder=swapOrder(vertex,iloc,momenta,fIndex)
    else :
        sorder=""
    momentastring=""
    for i in range(0,len(momenta)) :
        if(momenta[i][0] and momenta[i][1]!="")  :
            momentastring+=momenta[i][1]+"\n    "
    # special for 4-point VVVV
    if(vertex.lorentz[0].spins.count(3)==4 and iloc==0) :
        sig=sig.replace("Energy2","Energy2,int")
            
    header="virtual %s" % sig
    sig=sig.replace("=-GeV","")
    symboldefs = [ def_from_model(model,s) for s in symbols ]
    function = evaluateTemplate.format(decl=sig,momenta=momentastring,defns=defString,
                                       waves="\n    ".join(waves),symbols='\n    '.join(symboldefs),
                                       couplings="\n    ".join(localCouplings),
                                       result=result[0],swap=sorder)

    # special for transpose
    if(FM or RS) :
        h2=header
        if(not RS) :
            h2=header.replace("evaluate","evaluateN")
            function=function.replace("evaluate(","evaluateN(")
        headers=[]
        newHeader=""
        for ifunction in range(1,len(result)) :
            waveNew=[]
            momentastring=""
            htemp = header.split(",")
            irs=-1
            isp=-1
            # RS case
            if(RS) :
                # sort out the wavefunctions
                for i in range(0,len(waves)) :
                    if(waves[i].find("Spinor")<0) :
                        waveNew.append(waves[i])
                        continue
                    if(waves[i].find("Bar")>0) :
                        waveNew.append(waves[i].replace("Bar","").replace("bar",""))
                    else :
                        waveNew.append(waves[i].replace("Spinor","SpinorBar").replace(" s"," sbar").replace("Rs","Rsbar"))
                # sort out the momenta definitions
                for i in range(0,len(momenta)) :
                    if(momenta[i][0] and momenta[i][1]!="")  :
                        if(momenta[i][1].find("barW")>0) :
                            momentastring+=momenta[i][1].replace("barW","W")+"\n   "
                        elif(momenta[i][1].find("sW")>0) :
                            momentastring+=momenta[i][1].replace("sW","sbarW")+"\n   "
                        else :
                            momentastring+=momenta[i][1]+"\n    "
                # header string
                for i in range(0,len(htemp)) :
                    if(htemp[i].find("RS")>0) :
                        if(htemp[i].find("Bar")>0) :
                            htemp[i]=htemp[i].replace("Bar","").replace("RsbarW","RsW")
                        else :
                            htemp[i]=htemp[i].replace("Spinor","SpinorBar").replace("RsW","RsbarW")
                        if(i>0) : irs=i
                    elif(htemp[i].find("Spinor")>0) :
                        if(htemp[i].find("Bar")>0) :
                            htemp[i]=htemp[i].replace("Bar","").replace("sbarW","sW")
                        else :
                            htemp[i]=htemp[i].replace("Spinor","SpinorBar").replace("sW","sbarW")
                        if(i>0) : isp=i
                if(irs>0 and isp >0) :
                    htemp[irs],htemp[isp] = htemp[isp],htemp[irs]
            # fermion case
            else :
                htemp2 = header.split(",")
                # which fermions to exchange
                if(len(fIndex)==2) :
                    isp  = (fIndex[0],)
                    ibar = (fIndex[1],)
                else :
                    if(ifunction==1) :
                        isp  = (fIndex[2],)
                        ibar = (fIndex[3],)
                    elif(ifunction==2) :
                        isp  = (fIndex[0],)
                        ibar = (fIndex[1],)
                    elif(ifunction==3) :
                        isp  = (fIndex[0],fIndex[2])
                        ibar = (fIndex[1],fIndex[3])
                # wavefunctions
                for i in range(0,len(waves)) :
                    if(waves[i].find("Spinor")<0) :
                        waveNew.append(waves[i])
                        continue
                    if(waves[i].find("Bar")>0) :
                        found=False
                        for itest in range(0,len(ibar)) :
                            if(waves[i].find("sbarW%s"%ibar[itest])>=0) :
                                waveNew.append(waves[i].replace("Bar","").replace("bar",""))
                                found=True
                                break
                        if(not found) : waveNew.append(waves[i])
                    else :
                        found=False
                        for itest in range(0,len(isp)) :
                            if(waves[i].find("sW%s"%isp[itest])>=0) :
                                waveNew.append(waves[i].replace("Spinor","SpinorBar").replace(" s"," sbar"))
                                found=True
                                break
                        if(not found) : waveNew.append(waves[i])
                # momenta definitions
                for i in range(0,len(momenta)) :
                    if(momenta[i][0] and momenta[i][1]!="")  :
                        if(momenta[i][1].find("barW")>0) :
                            found = False
                            for itest in range(0,len(ibar)) :
                                if(momenta[i][1].find("barW%s"%ibar[itest])>=0) :
                                    momentastring+=momenta[i][1].replace("barW","W")+"\n   "
                                    found=True
                                    break
                            if(not found) :
                                momentastring+=momenta[i][1]+"\n    "
                        elif(momenta[i][1].find("sW")>0) :
                            found=False
                            for itest in range(0,len(isp)) :
                                if(momenta[i][1].find("sW%s"%isp[itest])>=0) :
                                    momentastring+=momenta[i][1].replace("sW","sbarW")+"\n   "
                                    found=True
                                    break
                            if(not found) :
                                momentastring+=momenta[i][1]+"\n    "
                        else :
                            momentastring+=momenta[i][1]+"\n    "
                # header
                for i in range(0,len(htemp)) :
                    if(htemp[i].find("Spinor")<0) : continue
                    if(htemp[i].find("Bar")>0) :
                        if(i==0) :
                            htemp[i] =htemp [i].replace("Bar","")
                            htemp2[i]=htemp2[i].replace("Bar","")
                            continue
                        for itest in range(0,len(ibar)) :
                            if(htemp[i].find("sbarW%s"%ibar[itest])>=0) :
                                htemp[i] =htemp [i].replace("Bar","").replace("sbarW","sW")
                                htemp2[i]=htemp2[i].replace("Bar","").replace("sbarW%s"%ibar[itest],"sW%s"%isp[itest])
                                break
                    else :
                        if(i==0) :
                            htemp [i]=htemp [i].replace("Spinor","SpinorBar")
                            htemp2[i]=htemp2[i].replace("Spinor","SpinorBar")
                            continue
                        for itest in range(0,len(isp)) :
                            if(htemp[i].find("sW%s"%isp[itest])>=0) :
                                htemp [i]=htemp [i].replace("Spinor","SpinorBar").replace("sW","sbarW")
                                htemp2[i]=htemp2[i].replace("Spinor","SpinorBar").replace("sW%s"%isp[itest],"sbarW%s"%ibar[itest])
                                break
            # header for transposed function
            hnew = ','.join(htemp)
            hnew = hnew.replace("virtual ","").replace("=-GeV","")
            if(not RS) :
                theader = ','.join(htemp2)
                theader = theader.replace("virtual ","").replace("=-GeV","")
                if(len(result)==2) :
                    hnew   =hnew   .replace("evaluate","evaluateT")
                    theader=theader.replace("evaluate","evaluateT")
                else :
                    hnew   =hnew   .replace("evaluate","evaluateT%s" % ifunction)
                    theader=theader.replace("evaluate","evaluateT%s" % ifunction)
                if(iloc not in fIndex) :
                    theader = headerReplace(theader)
                else :
                    theader = headerReplace(h2).replace("evaluateN","evaluateT")
                headers.append(theader)
                newHeader += hnew +";\n"
            else :
                newHeader += hnew
            fnew = evaluateTemplate.format(decl=hnew,momenta=momentastring,defns=defString,
                                           waves="\n    ".join(waveNew),symbols='\n    '.join(symboldefs),
                                           couplings="\n    ".join(localCouplings),
                                           result=result[ifunction],swap=sorder)
            function +="\n" + fnew


        if(FM and not RS) :
            if(len(result)==2) :
                if(iloc!=fIndex[1]) :
                    fi=1
                    stype="sbar"
                else :
                    fi=0
                    stype="s"
                header = vTemplateT.format(header=header.replace("Energy2,","Energy2 q2,"),
                                           normal=headerReplace(h2),
                                           transpose=theader,type=stype,
                                           iloc=fIndex[fi],id=vertex.particles[fIndex[fi]-1].pdg_code) \
                                           +newHeader+h2.replace("virtual","")
            else :
                sorder = swapOrderFFFF(vertex,iloc,fIndex)
                header = vTemplate4.format(header=header.replace("Energy2,","Energy2 q2,"),
                                           iloc1=fIndex[1],iloc2=fIndex[3],swap=sorder,
                                           id1=vertex.particles[fIndex[1]-1].pdg_code,
                                           id2=vertex.particles[fIndex[3]-1].pdg_code,
                                           cf=py2cpp(cf)[0],
                                           res1=headerReplace(h2).replace("W",""),
                                           res2=headers[0].replace("W",""),
                                           res3=headers[1].replace("W",""),
                                           res4=headers[2].replace("W",""))\
                                           +newHeader+h2.replace("virtual","")
        else :
            header=header + ";\n" + newHeader
    return (header,function)


evaluateMultiple = """\
{decl} {{
{code}
}}
"""

def multipleEvaluate(vertex,spin,defns) :
    if(spin==1) :
        name="scaW"
    elif(spin==3) :
        name="vW"
    else :
        print('Evaluate multiple problem',spin)
        raise SkipThisVertex()
    if(len(defns)==0) : return ("","")
    header = defns[0]
    ccdefn = header.replace("=-GeV","").replace("virtual ","").replace("Energy2","Energy2 q2")
    code=""
    spins=vertex.lorentz[0].spins
    iloc=1
    waves=[]
    for i in range(0,len(spins)) :
        if(spins[i]==spin) :
            waves.append("%s%s" %(name,i+1))
    for i in range(0,len(spins)) :
        if(spins[i]==spin) :
            if(iloc==1) : el=""
            else        : el="else "
            call = headerReplace(defns[iloc])
            if(iloc!=1) :
                call = call.replace(waves[0],waves[iloc-1])
            pdgid = vertex.particles[i].pdg_code
            code += "   %sif(out->id()==%s) return %s;\n" % (el,pdgid,call)
            iloc+=1
    code+="   else assert(false);\n"
    return (header,evaluateMultiple.format(decl=ccdefn,code=code))

            
