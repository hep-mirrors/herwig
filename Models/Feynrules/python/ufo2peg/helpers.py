from string import Template
from os import path
import sys,itertools,cmath
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





def qcd_qed_orders(vertex, coupling):
    # if more  than one take QCD over QED and then lowest order in QED
    if type(coupling) is list:
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
        qed = coupling.order.get('QED',0)
        qcd = coupling.order.get('QCD',0)
    # WARNING: FIX FOR CASES WHEN BOTH ARE ZERO
    # Is there a better way to treat this?
    if qed + qcd + 2 != len(vertex.particles):
        qed = len(vertex.particles) - qcd - 2

    return qcd, qed

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
        if len(v.particles) is 4:                                                                                                                      
            plistarray[1] += ',' + str(scfac[3] * v.particles[3].pdg_code)
        #print 'Conjugate vertex:', plistarray[1]

# ordering for EW VVV vertices
def VVVordering(vertex) :
    pattern = "if((p1->id()==%s&&p2->id()==%s&&p3->id()==%s)"+\
        "||(p1->id()==%s&&p2->id()==%s&&p3->id()==%s)||"+\
        "(p1->id()==%s&&p2->id()==%s&&p3->id()==%s)) {norm(-norm());}"
    ordering = pattern%(vertex.particles[1].pdg_code,
                        vertex.particles[0].pdg_code,
                        vertex.particles[2].pdg_code,
                        vertex.particles[0].pdg_code,
                        vertex.particles[2].pdg_code,
                        vertex.particles[1].pdg_code,
                        vertex.particles[2].pdg_code,
                        vertex.particles[1].pdg_code,
                        vertex.particles[0].pdg_code)
    return ordering

def tensorCouplings(vertex,coupling,prefactors,L,lorentztag,pos) :
    # split the structure into its different terms for analysis
    ordering=""
    structure1 = L.structure.split()
    structures =[]
    sign=''
    for struct in structure1 :
        if(struct=='+') :
            continue
        elif(struct=='-') :
            sign='-'
        else :
            structures.append(sign+struct.strip())
            sign=''
    lterms=[]
    rterms=[]
    if(lorentztag == 'SST') :
        terms=[['P(1003,2)','P(2003,1)'],
               ['P(1003,1)','P(2003,2)'],
               ['P(-1,1)','P(-1,2)','Metric(1003,2003)']]
        signs=[1.,1.,-1.]
    elif(lorentztag == 'FFT' ) :
        terms=[['P(2003,1)','Gamma(1003,2,1)'],
               ['P(2003,2)','Gamma(1003,2,1)'],
               ['P(1003,1)','Gamma(2003,2,1)'],
               ['P(1003,2)','Gamma(2003,2,1)'],
               ['P(-1,1)','Gamma(-1,2,1)','Metric(1003,2003)'],
               ['P(-1,2)','Gamma(-1,2,1)','Metric(1003,2003)']]
        signs=[1.,-1.,1.,-1.,-0.5,0.5]
    elif(lorentztag == 'VVT' ) :
        terms=[['P(-1,1)','P(-1,2)','Metric(1,2003)','Metric(2,1003)'],
               ['P(-1,1)','P(-1,2)','Metric(1,1003)','Metric(2,2003)'],
               ['P(-1,1)','P(-1,2)','Metric(1,2)','Metric(1003,2003)'],
               ['P(1,2)','P(2,1)','Metric(1003,2003)'],
               ['P(1,2)','P(2003,1)','Metric(2,1003)'],
               ['P(1,2)','P(1003,1)','Metric(2,2003)'],
               ['P(2,1)','P(2003,2)','Metric(1,1003)'],
               ['P(2,1)','P(1003,2)','Metric(1,2003)'],
               ['P(1003,2)','P(2003,1)','Metric(1,2)'],
               ['P(1003,1)','P(2003,2)','Metric(1,2)']]
        signs=[1.,1.,-1.,1.,-1.,-1.,-1.,-1.,1.,1.]
    elif(lorentztag == 'FFVT' ) :
        terms = [['Gamma(2004,2,1)','Metric(3,1004)'],
                 ['Gamma(1004,2,1)','Metric(3,2004)'],
                 ['Gamma(3,2,1)','Metric(1004,2004)']]
        lterms=[['Gamma(2004,2,-1)','Metric(3,1004)','ProjM(-1,1)'],
                ['Gamma(1004,2,-1)','Metric(3,2004)','ProjM(-1,1)'],
                ['Gamma(3,2,-1)','Metric(1004,2004)','ProjM(-1,1)']]
        rterms=[['Gamma(2004,2,-1)','Metric(3,1004)','ProjP(-1,1)'],
                ['Gamma(1004,2,-1)','Metric(3,2004)','ProjP(-1,1)'],
                ['Gamma(3,2,-1)','Metric(1004,2004)','ProjP(-1,1)']]
        signs=[1.,1.,-0.5]
    elif(lorentztag == 'VVVT' ) :
        # the F(mu nu,rho sigma lambda) terms first
        terms = [['P(2004,2)','Metric(1,1004)','Metric(2,3)'],['P(2004,3)','Metric(1,1004)','Metric(2,3)'],
                 ['P(1004,2)','Metric(1,2004)','Metric(2,3)'],['P(1004,3)','Metric(1,2004)','Metric(2,3)'],
                 ['P(2004,3)','Metric(1,3)','Metric(2,1004)'],['P(2004,1)','Metric(1,3)','Metric(2,1004)'],
                 ['P(1004,3)','Metric(1,3)','Metric(2,2004)'],['P(1004,1)','Metric(1,3)','Metric(2,2004)'],
                 ['P(2004,1)','Metric(1,2)','Metric(3,1004)'],['P(2004,2)','Metric(1,2)','Metric(3,1004)'],
                 ['P(1004,1)','Metric(1,2)','Metric(3,2004)'],['P(1004,2)','Metric(1,2)','Metric(3,2004)'],
                 ['P(3,1)','Metric(1,2004)','Metric(2,1004)'],['P(3,2)','Metric(1,2004)','Metric(2,1004)'], 
                 ['P(3,1)','Metric(1,1004)','Metric(2,2004)'],['P(3,2)','Metric(1,1004)','Metric(2,2004)'],
                 ['P(3,1)','Metric(1,2)','Metric(1004,2004)'],['P(3,2)','Metric(1,2)','Metric(1004,2004)'],
                 ['P(2,3)','Metric(1,2004)','Metric(3,1004)'],['P(2,1)','Metric(1,2004)','Metric(3,1004)'],
                 ['P(2,3)','Metric(1,1004)','Metric(3,2004)'],['P(2,1)','Metric(1,1004)','Metric(3,2004)'],
                 ['P(2,3)','Metric(1,3)','Metric(1004,2004)'],['P(2,1)','Metric(1,3)','Metric(1004,2004)'],
                 ['P(1,2)','Metric(2,2004)','Metric(3,1004)'],['P(1,3)','Metric(2,2004)','Metric(3,1004)'],
                 ['P(1,2)','Metric(2,1004)','Metric(3,2004)'],['P(1,3)','Metric(2,1004)','Metric(3,2004)'],
                 ['P(1,2)','Metric(2,3)','Metric(1004,2004)'],['P(1,3)','Metric(2,3)','Metric(1004,2004)']]
        signs = [1.,-1.,1.,-1.,1.,-1.,1.,-1.,1.,-1.,1.,-1.,
                 1.,-1.,1.,-1.,-1.,1.,1.,-1.,1.,-1.,-1.,1.,1.,-1.,1.,-1.,-1.,1.]
        signs = [1.,-1.,1.,-1.,1.,-1.,1.,-1.,1.,-1.,1.,-1.,
                 1.,-1.,1.,-1.,-1.,1.,1.,-1.,1.,-1.,-1.,1.,1.,-1.,1.,-1.,-1.,1.]
        l = lambda c: len(pos[c])
        if l(8)!=3 :
            ordering = VVVordering(vertex)
    # unknown
    else :
        raise Exception('Unknown data type "%s".' % lorentztag)
    sum   = [0.,0.,0.]
    itype = 0
    for types in (terms,lterms,rterms) :
        i=0
        for term in types:
            for perm in itertools.permutations(term):
                label = '*'.join(perm)
                for struct in structures :
                    if label in struct :
                        reminder = struct.replace(label,'1.',1)
                        sum[itype] += eval(reminder, {'cmath':cmath} )*signs[i]
            i+=1
        itype += 1
    all_coup   = []
    left_coup  = []
    right_coup = []
    if(len(lterms)==0) :
        all_coup.append('(%s) *(%s) * (%s)' % (sum[0]/float(len(signs)), prefactors,coupling.value))
    else :
        sum[1] += sum[0]
        sum[2] += sum[0]
        left_coup .append('(%s) * (%s) * (%s)' % (prefactors,sum[1]/float(len(signs)),coupling.value))
        right_coup.append('(%s) * (%s) * (%s)' % (prefactors,sum[2]/float(len(signs)),coupling.value))
    return (all_coup,left_coup,right_coup,ordering)

def EWVVVVCouplings(vertex,L) :
    terms=['Metric(1,2)*Metric(3,4)',
           'Metric(1,3)*Metric(2,4)',
           'Metric(1,4)*Metric(2,3)']

    structure1 = L.structure.split()

    structures =[]
    sign=''
    for struct in structure1 :
        if(struct=='+') :
            continue
        elif(struct=='-') :
            sign='-'
        else :
            structures.append(sign+struct.strip())
            sign=''
    factors=[]
    for term in terms:
        for struct in structures :
            if term in struct :
                reminder = struct.replace(term,'1.',1)
                try:
                    factors.append(eval(reminder, {'cmath':cmath} ))
                except NameError:
                    name_error = True
                else:
                    name_error = False

    if len(factors) != 3 or name_error:
        sys.stderr.write(
            'Warning: unsupported {tag} ( {ps} ) Lorentz structure in {name}:\n{lstr}\n'
            .format(tag=unique_lorentztag(vertex), name=vertex.name, 
                    lstr=L.structure, ps=' '.join(map(str,vertex.particles)))
        )
        raise SkipThisVertex()

    factor=0.
    order=[]
    if(factors[0]==-2.*factors[1] and factors[0]==-2.*factors[2] ) :
        order=[0,1,2,3]
        factor = factors[0]/2.
    elif(factors[1]==-2.*factors[0] and factors[1]==-2.*factors[2] ) :
        order=[0,2,1,3]
        factor = factors[1]/2.
    elif(factors[2]==-2.*factors[0] and factors[2]==-2.*factors[1] ) :
        order=[0,3,1,2]
        factor = factors[2]/2.
    else:
        sys.stderr.write(
            'Warning: unsupported {tag} ( {ps} ) Lorentz structure in {name}:\n{lstr}\n'
            .format(tag=unique_lorentztag(vertex), name=vertex.name, 
                    lstr=L.structure, ps=' '.join(map(str,vertex.particles)))
        )
        raise SkipThisVertex()


    pattern = \
        "bool done[4]={false,false,false,false};\n" + \
        "    tcPDPtr part[4]={p1,p2,p3,p4};\n" + \
        "    unsigned int iorder[4]={0,0,0,0};\n" + \
        "    for(unsigned int ix=0;ix<4;++ix) {\n" + \
        "       if(!done[0] && part[ix]->id()==%s) {done[0]=true; iorder[%s] = ix; continue;}\n" + \
        "       if(!done[1] && part[ix]->id()==%s) {done[1]=true; iorder[%s] = ix; continue;}\n" + \
        "       if(!done[2] && part[ix]->id()==%s) {done[2]=true; iorder[%s] = ix; continue;}\n" + \
        "       if(!done[3] && part[ix]->id()==%s) {done[3]=true; iorder[%s] = ix; continue;}\n" + \
        "    }\n" + \
        "    setType(2);\n" + \
        "    setOrder(iorder[0],iorder[1],iorder[2],iorder[3]);"
    ordering=pattern % ( vertex.particles[0].pdg_code,order[0],
                         vertex.particles[1].pdg_code,order[1],
                         vertex.particles[2].pdg_code,order[2],
                         vertex.particles[3].pdg_code,order[3] )
    return (ordering,factor)

def changeSign(sign1,sign2) :
    if((sign1=="+" and sign2=="+") or\
       (sign1=="-" and sign2=="-")) :
        return "+"
    else :
        return "-"

def epsilonOrder(eps) :
    terms = eps.strip("Epsilon(").strip(")").split(",")
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
    return (sign,"Epsilon(%s,%s,%s,%s)" % (terms[0],terms[1],terms[2],terms[3]))

    
def VVSEpsilon(couplings,struct) :
    if(struct.find("Epsilon")<0) :
        return
    fact=""
    sign="+"
    if(struct[-1]==")") :
        fact=struct.split("(")[0]
        if(fact.find("Epsilon")>=0) :
            fact=""
        else :
            struct=struct.split("(",1)[1][:-1]
            if(fact[0]=="-") :
                sign="-"
                fact=fact[1:]
    split = struct.split("*")
    # find the epsilon
    eps=""
    for piece in split :
        if(piece.find("Epsilon")>=0) :
            eps=piece
            split.remove(piece)
            break
    # and any prefactors
    for piece in split :
        if(piece.find("P(")<0) :
            split.remove(piece)
            if(piece[0]=="+" or piece[0]=="-") :
                sign=changeSign(sign,piece[0])
                piece=piece[1:]
            if(fact=="") :
                fact=piece
            else :
                fact = "( %s ) * ( %s )" % (fact , piece) 
    # now sort out the momenta
    for piece in split :
        terms=piece.split(",")
        terms[0]=terms[0].strip("P(")
        terms[1]=terms[1].strip(")")
        eps=eps.replace(terms[0],"P%s"%terms[1])
    (nsign,eps)=epsilonOrder(eps)
    if(nsign<0) : sign=changeSign(sign,"-")
    if(fact=="") : fact="1."
    if(eps!="Epsilon(1,2,P1,P2)") :
        return
    if(couplings[6]==0.) :
        couplings[6] = "( %s%s )" % (sign,fact)
    else :
        couplings[6] = "( %s ) + ( %s%s )" % (couplings[6],sign,fact)

def VVSCouplings(vertex,coupling,prefactors,L,lorentztag) :
    # split the structure into its different terms for analysis
    ordering=""
    structure1 = L.structure.split()
    structures =[]
    sign=''
    # extract the lorentz structures
    for struct in structure1 :
        if(struct=='+') :
            continue
        elif(struct=='-') :
            sign='-'
        else :
            structures.append(sign+struct.strip())
            sign=''
    # handle the scalar couplings
    couplings=[0.,0.,0.,0.,0.,0.,0.]
    terms=[['P(-1,1)','P(-1,2)','Metric(1,2)'],
           ['P(1,1)','P(2,1)'],
           ['P(1,1)','P(2,2)'],
           ['P(1,2)','P(2,1)'],
           ['P(1,2)','P(2,2)'],['Metric(1,2)']]
    itype=-1
    for term in terms:
        itype+=1
        for perm in itertools.permutations(term):
            label = '*'.join(perm)
            for istruct in range(0,len(structures)) :
                if label in structures[istruct] :
                    reminder = structures[istruct].replace(label,'1.',1)
                    couplings[itype]+=eval(reminder, {'cmath':cmath} )
                    structures[istruct]='Done'
    # handle the pseudoscalar couplings
    for struct in structures :
        if(struct != "Done" ) :
            VVSEpsilon(couplings,struct)
    # evaluate the prefactor
    if type(coupling) is not list:
        value = coupling.value
    else:
        value = "("
        for coup in coupling :
            value += '+(%s)' % coup.value
        value +=")"
    # put it all together
    for ic in range(0,len(couplings)) :
        if(couplings[ic]!=0.) :
            couplings[ic] = '(%s) * (%s) * (%s)' % (prefactors,value,couplings[ic])
    return couplings
