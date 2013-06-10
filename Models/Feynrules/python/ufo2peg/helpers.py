from string import Template
from os import path
#import re

"""
Helper functions for the Herwig++ Feynrules converter
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
    """Create a template from a file."""
    templatename = '{}.template'.format(name)
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



def get_lorentztag(spin):
    """Produce a ThePEG spin tag for the given numeric FR spins."""
    spins = { 1 : 'S', 2 : 'F', 3 : 'V', -1 : 'U', 5 : 'T' }
    result = [ spins[s] for s in spin ]

    def spinsort(a,b):
        """Helper function for ThePEG's FVST spin tag ordering."""
        if a == b: return 0
        for letter in 'FVST':
            if a == letter: return -1
            if b == letter: return  1

    result = sorted(result, cmp=spinsort)
    return ''.join(result)

def unique_lorentztag(vertex):
    """Check and return the Lorentz tag of the vertex."""
    unique = CheckUnique()
    for l in vertex.lorentz:
        lorentztag = get_lorentztag(l.spins)
        unique( lorentztag )
        assert( lorentztag == l.name[:len(lorentztag)] )
    return lorentztag


def spindirectory(lt):
    """Return the spin directory name for a given Lorentz tag."""
    if 'T' in lt: 
        spin_directory = 'Tensor'
    elif 'S' in lt: 
        spin_directory = 'Scalar'
    elif 'V' in lt: 
        spin_directory = 'Vector'
    else:
        raise Exception("Unknown Lorentz tag {}.".format(lt))
    return spin_directory




def colorpositions(struct):
    positions = { 
        1 : [],
        3 : [],
        -3 : [],
        8 : [],
    }
    for i,s in enumerate(struct,1):
        positions[s].append(i)
    return positions
        
def colorfactor(vertex):
    try:
        unique = CheckUnique()
        for pl in vertex.particles_list:
            struct = [ p.color for p in pl ]
            unique(struct)
    except:
        struct = [ p.color for p in vertex.particles ]
    pos = colorpositions(struct)

    def match(patterns):
        result = [ p == t
                   for p,t in zip(patterns,vertex.color) ]
        return all(result)

    label = None
    L = len(struct)
    l = lambda c: len(pos[c])
    if l(1) == L:
        label = ('1',)
        if match(label): return ('1',)

    elif l(3) == l(-3) == 1 and l(1) == L-2:
        nums = [pos[3][0], pos[-3][0]]
        label = ('Identity({},{})'.format(*sorted(nums)),)
        if match(label): return ('1',)

    elif l(8) == l(3) == l(-3) == 1 and l(1) == L-3:
        label = ('T({},{},{})'.format(pos[8][0],pos[3][0],pos[-3][0]),)
        if match(label): return ('1',)

    elif l(8) == L == 3:
        # if lorentz is FFV get extra minus sign
        lorentztag = unique_lorentztag(vertex)
        factor = '*(-1)' if lorentztag in ['FFV'] else ''
        label = ('f(1,2,3)',)
        if match(label): return ('-complex(0,1)%s'%factor,)
        label = ('f(3,2,1)',)
        if match(label): return ('complex(0,1)%s'%factor,)

    elif l(8) == L == 4:
        label = ('f(-1,1,2)*f(3,4,-1)',
                 'f(-1,1,3)*f(2,4,-1)',
                 'f(-1,1,4)*f(2,3,-1)',
             )
        if match(label): return ('1','1','1')

    elif l(8) == 2 and l(3) == l(-3) == 1:
        subs = {
            'g1' : pos[8][0],
            'g2' : pos[8][1],
            'qq' : pos[3][0],
            'qb' : pos[-3][0] 
        }
        label = ('T({g1},-1,{qb})*T({g2},{qq},-1)'.format(**subs),
                 'T({g1},{qq},-1)*T({g2},-1,{qb})'.format(**subs))
        if match(label): return ('0.5','0.5')

    raise Exception("Unknown colour tag {}.".format(vertex.color))


def def_from_model(FR,s):
    """Return a C++ line that defines parameter s as coming from the model file."""
    stype = typemap(getattr(FR.parameters,s).type)
    return '{t} {s} = model_->{s}();'.format(t=stype,s=s)

_typemap = {'complex':'Complex',
            'real':'double'}

def typemap(s):
    return _typemap[s]

def add_brackets(expr, syms):
    result = expr
    for s in syms:
        result = result.replace(s,s+'()')
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
        if len(v.particles) is 4:                                                                                                                      
            plistarray[1] += ',' + str(scfac[3] * v.particles[3].pdg_code)
        #print 'Conjugate vertex:', plistarray[1]
