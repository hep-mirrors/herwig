from __future__ import print_function
from string import Template
from os import path
import sys,cmath,glob
import re

"""
Helper functions for the Herwig Feynrules converter
"""

class CheckUnique:
    """Uniqueness checker.

    An object of this class remembers the value it was called with first.
    Any subsequent call to it will only succeed if the same value is passed again.
    For example,

    >>> f = CheckUnique()
    >>> f(5)
    >>> f(5)
    >>> f(4)
    Traceback (most recent call last):
        ...
    AssertionError
    """
    def __init__(self):
        self.val = None

    def __call__(self,val):
        """Store value on first call, then compare."""
        if self.val is None:
            self.val = val
        else:
            assert( val == self.val )



def is_number(s):
    """Check if a value is a number."""
    try:
        float(s)
    except ValueError:
        return False
    return True



def getTemplate(name):
    """Create a python string template from a file."""
    templatename = '{name}.template'.format(name=name)
    # assumes the template files sit next to this script
    moduledir = path.dirname(path.abspath(__file__))
    templatepath = path.join(moduledir,templatename)
    with open(templatepath, 'r') as f:
        templateText = f.read()
    return Template( templateText )

def writeFile(filename, text):
    """Write text to a filename."""
    with open(filename,'w') as f:
        f.write(text)

def coupling_orders(vertex, coupling, defns):
    # if more  than one take QCD over QED and then lowest order in QED
    if type(coupling) is list:
        print('not sure this happens')
        quit()
        qed = 0
        qcd = 0
        for coup in coupling :
            qed1 = coup.order.get('QED',0)
            qcd1 = coup.order.get('QCD',0)
            if qcd1 != 0:
                if qcd == 0 or (qcd1 != 0 and qcd1 < qcd):
                    qcd=qcd1
                    qed=qed1
            else:
                if qed == 0 or (qed1 != 0 and qed1 < qed):
                    qed=qed1
    else:
        output={}
        for ctype in defns :
            output[ctype]=coupling.order.get(ctype,0)
    return output

def def_from_model(FR,s):
    """Return a C++ line that defines parameter s as coming from the model file."""
    if("hw_kine" in s) :return ""
    stype = typemap(getattr(FR.parameters,s).type)
    return '{t} {s} = model_->{s}();'.format(t=stype,s=s)

_typemap = {'complex':'Complex',
            'real':'double'}

def typemap(s):
    return _typemap[s]

def add_brackets(expr, syms):
    result = expr
    for s in syms:
        pattern = r'({symb})(\W|$)'.format(symb=s)
        result = re.sub(pattern, r'\1()\2', result)
    return result




def banner():
    return """\
===============================================================================================================
______                  ______        _                 __ _   _                        _
|  ___|                 | ___ \      | |               / /| | | |                      (_)          _      _
| |_  ___  _   _  _ __  | |_/ /_   _ | |  ___  ___    / / | |_| |  ___  _ __ __      __ _   __ _  _| |_  _| |_
|  _|/ _ \| | | || \_ \ |    /| | | || | / _ \/ __|  / /  |  _  | / _ \| \__|\ \ /\ / /| | / _` ||_   _||_   _|
| | |  __/| |_| || | | || |\ \| |_| || ||  __/\__ \ / /   | | | ||  __/| |    \ V  V / | || (_| |  |_|    |_|
\_|  \___| \__, ||_| |_|\_| \_|\__,_||_| \___||___//_/    \_| |_/ \___||_|     \_/\_/  |_| \__, |
            __/ |                                                                           __/ |
           |___/                                                                           |___/
===============================================================================================================
generating model/vertex/.model/.in files
please be patient!
===============================================================================================================
"""




#################### ??? #######################


# function that replaces alphaS (aS)-dependent variables
# with their explicit form which also contains strongCoupling
def aStoStrongCoup(stringin, paramstoreplace, paramstoreplace_expressions):
    #print stringin
    for xx in range(0,len(paramstoreplace)):
        #print paramstoreplace[xx], paramstoreplace_expressions[xx]
        stringout = stringin.replace(paramstoreplace[xx], '(' +  PyMathToThePEGMath(paramstoreplace_expressions[xx],allparams) + ')')
    stringout = stringout.replace('aS', '(sqr(strongCoupling(q2))/(4.0*Constants::pi))')
    #print 'resulting string', stringout
    return stringout


# function that replaces alphaEW (aEW)-dependent variables
# with their explicit form which also contains weakCoupling
def aEWtoWeakCoup(stringin, paramstoreplace, paramstoreplace_expressions):
    #print stringin
    for xx in range(0,len(paramstoreplace)):
        #print paramstoreplace[xx], paramstoreplace_expressions[xx]
        stringout = stringin.replace(paramstoreplace[xx], '(' +  PyMathToThePEGMath(paramstoreplace_expressions[xx],allparams) + ')')
    stringout = stringout.replace('aEWM1', '(1/(sqr(electroMagneticCoupling(q2))/(4.0*Constants::pi)))')
    #print 'resulting string', stringout
    return stringout




if __name__ == "__main__":
    import doctest
    doctest.testmod()


if False:

# Check if the Vertex is self-conjugate or not
    pdgcode = [0,0,0,0]
    notsmvertex = False
    vhasw = 0
    vhasz = 0
    vhasf = 0
    vhasg = 0
    vhash = 0
    vhasp = 0
#   print 'printing particles in vertex'
    for i in range(len(v.particles)):
#       print v.particles[i].pdg_code
        pdgcode[i] = v.particles[i].pdg_code
        if pdgcode[i] == 23:
            vhasz += 1
        if pdgcode[i] == 22:
            vhasp += 1
        if pdgcode[i] == 25:
            vhash += 1
        if pdgcode[i] == 21:
            vhasg += 1
        if pdgcode[i] == 24:
            vhasw += 1
        if abs(pdgcode[i]) < 7 or (abs(pdgcode[i]) > 10 and abs(pdgcode[i]) < 17):
            vhasf += 1
        if pdgcode[i] not in SMPARTICLES:
            notsmvertex = True


#  treat replacement of SM vertices with BSM vertices?
    if notsmvertex == False:
        if( (vhasf == 2 and vhasz == 1) or (vhasf == 2 and vhasw == 1) or (vhasf == 2 and vhash == 1) or (vhasf == 2 and vhasg == 1) or (vhasf == 2 and vhasp == 0) or (vhasg == 3) or (vhasg == 4) or (vhasw == 2 and vhash == 1) or (vhasw == 3) or (vhasw == 4) or (vhash == 1 and vhasg == 2) or (vhash == 1 and vhasp == 2)):
            #print 'VERTEX INCLUDED IN STANDARD MODEL!'
            v.include = 0
            notincluded += 1
            #continue


    selfconjugate = 0
    for j in range(len(pdgcode)):
        for k in range(len(pdgcode)):
               if  j != k and j != 0 and abs(pdgcode[j]) == abs(pdgcode[k]):
                   selfconjugate = 1
                   #print 'self-conjugate vertex'
#        print pdgcode[j]

# if the Vertex is not self-conjugate, then add the conjugate vertex
# automatically
    scfac = [1,1,1,1]
    if selfconjugate == 0:
        #first find the self-conjugate particles
        for u in range(len(v.particles)):
              if v.particles[u].selfconjugate == 0:
                  scfac[u] = -1
#                  print 'particle ', v.particles[u].pdg_code, ' found not to be self-conjugate'

    if selfconjugate == 0:
        plistarray[1] += str(scfac[1] * v.particles[1].pdg_code) + ',' + str(scfac[0] * v.particles[0].pdg_code) + ',' + str(scfac[2] * v.particles[2].pdg_code)
        if len(v.particles) == 4:
            plistarray[1] += ',' + str(scfac[3] * v.particles[3].pdg_code)
        #print 'Conjugate vertex:', plistarray[1]


class SkipThisVertex(Exception):
    pass

def extractAntiSymmetricIndices(instring,funct) :
    terms = instring.strip(funct).strip(")").split(",")
    sign=1.
    for iy in range(0,len(terms)) :
        for ix in range(-1,-len(terms)+iy,-1) :
            swap = False
            if(len(terms[ix])==1 and len(terms[ix-1])==1) :
                swap = int(terms[ix])<int(terms[ix-1])
            elif(len(terms[ix])==2 and len(terms[ix-1])==2) :
                swap = int(terms[ix][1])<int(terms[ix-1][1])
            elif(len(terms[ix])==1 and len(terms[ix-1])==2) :
                swap = True
            if(swap) :
                sign *=-1.
                terms[ix],terms[ix-1] = terms[ix-1],terms[ix]
    return (terms,sign)

def isGoldstone(p) :
    """check if particle is a Goldstone"""
    def gstest(name):
        try:
            return getattr(p,name)
        except AttributeError:
            return False
    # names of goldstone bosons
    gsnames = ['goldstone','goldstoneboson','GoldstoneBoson']
    if any(map(gstest, gsnames)):
        return True
    return False

def isGhost(p) :
    """Check if particle is a ghost"""
    try:
        getattr(p,'GhostNumber')
    except AttributeError:
        return False
    return p.GhostNumber != 0

def convertToPython3(ufodir) :
    # find all the python files
    fNames=glob.glob(path.abspath(ufodir)+"/*.py")
    mNames=[]
    for i in range(0,len(fNames)) :
        mNames.append(path.split(fNames[i])[-1].replace(".py",""))
    # convert them
    for fName in fNames :
        convertFileToPython3(fName,mNames)

def copy(path1, path2):
    path1 = format_path(path1)
    path2 = format_path(path2)
    shutil.copy(path1, path2)

def prepForConversion(ufodir) :
    import os
    path=os.path.abspath(ufodir)
    if not os.path.isdir(path):
        raise Exception( 'path to the UFO directory seems to be wrong!')
    model_dir = path
    text = open(os.path.join(model_dir, 'object_library.py')).read()
    text = text.replace('.iteritems()', '.items()')
    text = re.sub('raise (\w+)\s*,\s*["\']([^"]+)["\']','raise \g<1>("\g<2>")', text)
    text = open(os.path.join(model_dir, 'object_library.py'),'w').write(text)
    text = open(os.path.join(model_dir, '__init__.py')).read()
    mod = False
    to_check =  ['object_library', 'function_library']
    for lib in to_check:
        if 'import %s' % lib in text:
            continue
        mod = True
        text = "import %s \n" % lib + text
    if mod:
        open(os.path.join(model_dir, '__init__.py'),'w').write(text)

def convertFileToPython3(fName,names) :
    output=""
    inFile=open(fName)
    line=inFile.readline()
    isInit = "__init__.py" in fName
    tNames=[]
    for val in names : tNames.append(val)
    while line :
        # iteritems -> items
        line=line.replace("iteritems","items")
        # fix imports
        if("import" in line) :
            for val in names :
                if("import %s" %val in line) :
                    line=line.replace("import %s"  %val,
                                      "from . import %s" %val)
                    if(val in tNames) : tNames.remove(val)
                if("from %s" %val in line) :
                    line=line.replace("from %s"  %val,
                                      "from .%s" %val)
        # add brackets to print statements
        if("print" in line) :
            if line.strip()[0:5] == "print" :
                line=line.strip("\n").replace("print","print(")+")\n"
        output += line
        line=inFile.readline()
    inFile.close()
    if(isInit) :
        if "__init__" in tNames : tNames.remove("__init__")
        temp=""
        for val in tNames :
            temp+= "from . import %s\n" % val
        output = temp+output
    with open(fName,'w') as dest:
        dest.write(output)
