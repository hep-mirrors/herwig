import itertools,cmath,re,sys,copy
from .helpers import SkipThisVertex,extractAntiSymmetricIndices,def_from_model
from .converter import py2cpp
from .lorentzparser import parse_lorentz
import sympy,string
from string import Template
from sympy import Matrix,Symbol

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

imap=["t","x","y","z"]

RSDotProduct = Template("${s}ts${si}*${v}t-${s}xs${si}*${v}x-${s}ys${si}*${v}y-${s}zs${si}*${v}z")

vTemplateT="""\
    if(sbarW{iloc}.id()=={id}) {{
        return ({res})*({cf});
    }}
    else {{
        return ({resT})*({cf});
    }}
"""

vTemplate4="""\
    if(id{iloc1}=={id1}) {{
        if(id{iloc2}=={id2}) {{
          return ({res1})*({cf});
        }}
        else {{
          return ({res2})*({cf});
        }}
    }}
    else {{
        if(id{iloc2}=={id2}) {{
          return ({res3})*({cf});
        }}
        else {{
          return ({res4})*({cf});
        }}
    }}
"""

vecTemplate="""\
    LorentzPolarizationVector vtemp = {res};
    Energy2 p2 = P{iloc}.m2();
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

vecTTemplate="""\
    LorentzPolarizationVector vtemp;
    if(sbarW{isp}.id()=={id}) {{
       vtemp = {res};
    }}
    else {{
       vtemp  = {resT};
    }}
    Energy2 p2 = P{iloc}.m2();
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
    if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
    Energy2 p2 = P{iloc}.m2();
    Complex fact = Complex(0.,1.)*({cf})*propagator(iopt,p2,out,mass,width);
    Lorentz{offTypeA}<double> newSpin = fact*({res});
    return {offTypeB}(P{iloc},out,newSpin);
"""

RSTemplate="""\
    if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
    Energy2 p2 = P{iloc}.m2();
    Complex fact = Complex(0.,1.)*({cf})*propagator(iopt,p2,out,mass,width);
    complex<InvEnergy> Omass = mass.real()==ZERO ? InvEnergy(ZERO) : 1./mass;
    Lorentz{offTypeA}<double> newSpin = fact*({res});
    return {offTypeB}(P{iloc},out,newSpin);
"""

sTTemplate="""\
    if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
    Energy2 p2 = P{iloc}.m2();
    Complex fact = Complex(0.,1.)*({cf})*propagator(iopt,p2,out,mass,width);
    if(out->id()=={id}) {{
        Lorentz{offTypeA}<double> newSpin = fact*({res});
        return {offTypeB}(P{iloc},out,newSpin);
    }}
    else {{
        Lorentz{offTypeA}<double> newSpin = fact*({resT});
        return {offTypeB}(P{iloc},out,newSpin);
    }}
"""

scaTemplate="""\
    if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
    Energy2 p2 = P{iloc}.m2();
    Complex fact = Complex(0.,1.)*({cf})*propagator(iopt,p2,out,mass,width);
    complex<double> output = fact*({res});
    return ScalarWaveFunction(P{iloc},out,output);
"""

scaTTemplate="""\
    if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
    Energy2 p2 = P{iloc}.m2();
    Complex fact = Complex(0.,1.)*({cf})*propagator(iopt,p2,out,mass,width);
    complex<double> output = sbarW{isp}.id()=={id} ? fact*({res}) : fact*({resT});
    return ScalarWaveFunction(P{iloc},out,output);
"""

tenTemplate="""\
    if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
    InvEnergy2 OM{iloc} = mass.real()==ZERO ? InvEnergy2(ZERO) : 1./sqr(mass.real()); 
    Energy2 p2 = P{iloc}.m2();
    Complex fact = Complex(0.,1.)*({cf})*propagator(iopt,p2,out,mass,width);
    LorentzTensor<double> output = fact*({res});
    return TensorWaveFunction(P{iloc},out,output);
"""

tenTTemplate="""\
    if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();
    InvEnergy2 OM{iloc} = mass.real()==ZERO ? InvEnergy2(ZERO) : 1./sqr(mass.real()); 
    Energy2 p2 = P{iloc}.m2();
    Complex fact = Complex(0.,1.)*({cf})*propagator(iopt,p2,out,mass,width);
    LorentzTensor<double> output = sbarW{isp}.id()=={id} ? fact*({res}) : fact*({resT});
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

def compare(a,b) :
    num=abs(a-b)
    den=abs(a+b)
    if(den == 0. and 1e-10) :
        return True
    return num/den<1e-10

def evaluate(x,model,parmsubs):
    import cmath
    return eval(x, 
                {'cmath':cmath,
                 'complexconjugate':model.function_library.complexconjugate}, 
                parmsubs)

def pdgCC(particle) :
    pdgid = particle.pdg_code
    if(particle.name!=particle.antiname) :
        pdgid *= -1
    return pdgid
# ordering for EW VVV vertices (ordering not an issue as all same spin)
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

def tensorCouplings(vertex,value,prefactors,L,lorentztag,pos,all_couplings,order) :
    # split the structure into its different terms for analysis
    ordering=""
    structures = extractStructures(L)
    if(lorentztag == 'SST') :
        terms=[['P(1003,2)','P(2003,1)'],
               ['P(1003,1)','P(2003,2)'],
               ['P(-1,1)','P(-1,2)','Metric(1003,2003)'],
               ['Metric(1003,2003)']]
        signs=[1.,1.,-1.,-1.]
        new_couplings=[False]*len(terms)
    elif(lorentztag == 'FFT' ) :
        terms=[['P(2003,1)','Gamma(1003,2,1)'],
               ['P(2003,2)','Gamma(1003,2,1)'],
               ['P(1003,1)','Gamma(2003,2,1)'],
               ['P(1003,2)','Gamma(2003,2,1)'],
               ['P(-1,1)','Gamma(-1,2,1)','Metric(1003,2003)'],
               ['P(-1,2)','Gamma(-1,2,1)','Metric(1003,2003)'],
               ['Metric(1003,2003)']]
        signs=[1.,-1.,1.,-1.,-0.5,0.5,1.]
        new_couplings=[False]*3*len(terms)
    elif(lorentztag == 'VVT' ) :
        terms=[['P(-1,1)','P(-1,2)','Metric(1,2003)','Metric(2,1003)'], # from C term
               ['P(-1,1)','P(-1,2)','Metric(1,1003)','Metric(2,2003)'], # from C term
               ['P(-1,1)','P(-1,2)','Metric(1,2)','Metric(1003,2003)'], # from C term
               ['P(1,2)','P(2,1)','Metric(1003,2003)'], # from D term (sym)
               ['P(1,2)','P(2003,1)','Metric(2,1003)'], # 1st term
               ['P(1,2)','P(1003,1)','Metric(2,2003)'], # 1st swap
               ['P(2,1)','P(2003,2)','Metric(1,1003)'], # 2nd term
               ['P(2,1)','P(1003,2)','Metric(1,2003)'], # 2nd swap
               ['P(1003,2)','P(2003,1)','Metric(1,2)'], # 3rd term 
               ['P(1003,1)','P(2003,2)','Metric(1,2)'], # 3rd swap
               ['Metric(1,2003)','Metric(2,1003)'], # from mass term
               ['Metric(1,1003)','Metric(2,2003)'], # from mass term
               ['Metric(1,2)','Metric(1003,2003)'], # from mass term
               ['P(1,1)','P(2,1)','Metric(1003,2003)'], # gauge terms
               ['P(1,2)','P(2,2)','Metric(1003,2003)'], # gauge terms
               ['P(1,1)','P(2,2)','Metric(1003,2003)'], # gauge terms
               ['P(1003,1)','P(1,1)','Metric(2,2003)'], # gauge terms
               ['P(1003,2)','P(2,2)','Metric(1,2003)'], # gauge terms
               ['P(2003,1)','P(1,1)','Metric(2,1003)'], # gauge terms
               ['P(2003,2)','P(2,2)','Metric(1,1003)'], # gauge terms
        ]
        signs=[1.,1.,-1.,1.,-1.,-1.,-1.,-1.,1.,1.,1.,1.,-1.,1.,1.,0.25,-1.,-1.,-1.,-1.]
        new_couplings=[False]*len(terms)
    elif(lorentztag == 'FFVT' ) :
        terms = [['Gamma(2004,2,1)','Metric(3,1004)'],
                 ['Gamma(1004,2,1)','Metric(3,2004)'],
                 ['Gamma(3,2,1)','Metric(1004,2004)'],
                 ['Gamma(2004,2,-1)','Metric(3,1004)'],
                 ['Gamma(1004,2,-1)','Metric(3,2004)'],
                 ['Gamma(3,2,-1)','Metric(1004,2004)']]
        signs=[1.,1.,-0.5,1.,1.,-0.5]
        new_couplings=[False]*3*len(terms)
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
        new_couplings=[False]*len(terms)
        l = lambda c: len(pos[c])
        if l(8)!=3 :
            ordering = VVVordering(vertex)
    # unknown
    else :
        raise Exception('Unknown data type "%s".' % lorentztag)
    iterm=0
    try :
        for term in terms:
            for perm in itertools.permutations(term):
                label = '*'.join(perm)
                for istruct in range(0,len(structures)) :
                    if label in structures[istruct] :
                        reminder = structures[istruct].replace(label,'1.',1)
                        loc=iterm
                        if(reminder.find("ProjM")>=0) :
                            reminder=re.sub("\*ProjM\(.*,.\)","",reminder)
                            loc+=len(terms)
                        elif(reminder.find("ProjP")>=0) :
                            reminder=re.sub("\*ProjP\(.*,.\)","",reminder)
                            loc+=2*len(terms)
                        structures[istruct] = "Done"
                        val = eval(reminder, {'cmath':cmath} )*signs[iterm]
                        if(new_couplings[loc]) :
                            new_couplings[loc] += val
                        else :
                            new_couplings[loc] = val
            iterm+=1
    except :
        SkipThisVertex()
    # check we've handled all the terms
    for val in structures:
        if(val!="Done") :
            raise SkipThisVertex()
    # special for FFVT
    if(lorentztag=="FFVT") :
        t_couplings=new_couplings
        new_couplings=[False]*9
        for i in range(0,9) :
            j = i+3*(i/3)
            k = i+3+3*(i/3)
            if( not t_couplings[j]) :
                new_couplings[i] = t_couplings[k]
            else :
                new_couplings[i] = t_couplings[j]
    # set the couplings
    for icoup in range(0,len(new_couplings)) :
        if(new_couplings[icoup]) :
            new_couplings[icoup] = '(%s) * (%s) *(%s)' % (new_couplings[icoup],prefactors,value)
    if(len(all_couplings)==0) :
        all_couplings=new_couplings
    else :
        for icoup in range(0,len(new_couplings)) :
            if(new_couplings[icoup] and all_couplings[icoup]) :
                all_couplings[icoup] = '(%s) + (%s) ' % (new_couplings[icoup],all_couplings[icoup])
            elif(new_couplings[icoup]) :
                all_couplings[icoup] = new_couplings[icoup]
    # return the results
    return (ordering,all_couplings)

def processTensorCouplings(lorentztag,vertex,model,parmsubs,all_couplings,order) :
    # check for fermion vertices (i.e. has L/R couplings)
    fermions = "FF" in lorentztag
    # test and process the values of the couplings
    tval  = [False]*3
    value = [False]*3
    # loop over the colours
    for icolor in range(0,len(all_couplings)) :
        lmax = len(all_couplings[icolor])
        if(fermions) : lmax /=3
        # loop over the different terms
        for ix in range(0,lmax) :
            test = [False]*3
            # normal case
            if( not fermions ) :
                test[0] = all_couplings[icolor][ix]
            else :
                # first case vector but no L/R couplings
                if( not all_couplings[icolor][lmax+ix]  and
                    not all_couplings[icolor][2*lmax+ix] ) :
                    test[0] = all_couplings[icolor][ix]
                    # special for mass terms and massless particles
                    if(not all_couplings[icolor][ix]) :
                        code = abs(vertex.particles[order[0]-1].pdg_code)
                        if(ix==6 and code ==12 or code ==14 or code==16) :
                            continue
                        else :
                            raise SkipThisVertex()
                # second case L/R couplings
                elif( not all_couplings[icolor][ix] ) :
                    # L=R, replace with vector
                    if(all_couplings[icolor][lmax+ix] ==
                       all_couplings[icolor][2*lmax+ix]) :
                        test[0]  = all_couplings[icolor][lmax+ix]
                    else :
                        test[1] = all_couplings[icolor][lmax+ix]
                        test[2] = all_couplings[icolor][2*lmax+ix]
                else :
                    raise SkipThisVertex()
            # special handling of mass terms
            # scalar divide by m**2
            if((ix==3 and lorentztag=="SST") or
               ( ix>=10 and ix<=12 and lorentztag=="VVT" )) :
                for i in range(0,len(test)) :
                    if(test[i]) :
                        test[i] = '(%s)/%s**2' % (test[i],vertex.particles[order[0]-1].mass.value)
            # fermions divide by 4*m
            elif(ix==6 and lorentztag=="FFT" and
                 float(vertex.particles[order[0]-1].mass.value) != 0. ) :
                for i in range(0,len(test)) :
                    if(test[i]) :
                        test[i] = '-(%s)/%s/4' % (test[i],vertex.particles[order[0]-1].mass.value)
            # set values on first pass
            if(not tval[0] and not tval[1] and not tval[2]) :
                value = test
                for i in range(0,len(test)) :
                    if(test[i]) : tval[i] = evaluate(test[i],model,parmsubs)
            else :
                for i in range(0,len(test)) :
                    if(not test[i] and not tval[i]) :
                        continue
                    if(not test[i] or not tval[i]) :
                        # special for mass terms and vectors
                        if(lorentztag=="VVT" and ix >=10 and ix <=12 and
                           float(vertex.particles[order[0]-1].mass.value) == 0. ) :
                            continue
                        # special for vector gauge terms
                        if(lorentztag=="VVT" and ix>=13) :
                            continue
                        raise SkipThisVertex()
                    tval2 = evaluate(test[i],model,parmsubs)
                    if(abs(tval[i]-tval2)>1e-6) :
                        # special for fermion mass term if fermions massless
                        if(lorentztag=="FFT" and ix ==6 and tval2 == 0. and
                           float(vertex.particles[order[0]-1].mass.value) == 0. ) :
                            continue
                        raise SkipThisVertex()
    # simple clean up
    for i in range(0,len(value)):
        if(value[i]) :
            value[i] = value[i].replace("(1.0) * ","").replace(" * (1)","")
    # put everything together
    coup_left  = 0.
    coup_right = 0.
    coup_norm  = 0.
    if(lorentztag ==  "SST" or lorentztag == "VVT" or
       lorentztag == "VVVT" or lorentztag == "FFT" ) :
        coup_norm = value[0]
        if(value[1] or value[2]) :
            raise SkipThisVertex()
    elif(lorentztag=="FFVT") :
        if(not value[1] and not value[2]) :
            coup_norm  = value[0]
            coup_left  = "1."
            coup_right = "1."
        elif(not value[0]) :
            coup_norm = "1."
            if(value[1] and value[2]) :
                coup_left  = value[1]
                coup_right = value[2]
            elif(value[1]) :
                coup_left  = value[1]
                coup_right = "0."
            elif(value[2]) :
                coup_left  = "0."
                coup_right = value[2]
            else :
                raise SkipThisVertex()
        else :
            raise SkipThisVertex()
    else :
        raise SkipThisVertex()
    # return the answer
    return (coup_left,coup_right,coup_norm)

def extractStructures(L) :
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
    return structures



def changeSign(sign1,sign2) :
    if((sign1=="+" and sign2=="+") or\
       (sign1=="-" and sign2=="-")) :
        return "+"
    else :
        return "-"
    
def epsilonOrder(eps) :
    terms,sign = extractAntiSymmetricIndices(eps,"Epsilon(")
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
    if(nsign>0) : sign=changeSign(sign,"-")
    if(fact=="") : fact="1."
    if(eps!="Epsilon(1,2,P1,P2)") :
        return
    if(couplings[6]==0.) :
        couplings[6] = "( %s%s )" % (sign,fact)
    else :
        couplings[6] = "( %s ) + ( %s%s )" % (couplings[6],sign,fact)


def scalarVectorCouplings(value,prefactors,L,lorentztag,all_couplings,order) :
    # set up the types of term we are looking for
    if(lorentztag=="VVS") :
        couplings=[0.,0.,0.,0.,0.,0.,0.]
        terms=[['P(-1,%s)' % order[0],
                'P(-1,%s)' % order[1],
                'Metric(%s,%s)' %(order[0],order[1])],
               ['P(1,%s)' % order[0],
                'P(2,%s)' % order[0]],
               ['P(1,%s)' % order[0],
                'P(2,%s)' % order[1]],
               ['P(1,%s)' % order[1],
                'P(2,%s)' % order[0]],
               ['P(1,%s)' % order[1],
                'P(2,%s)' % order[1]],
               ['Metric(%s,%s)'%(order[0],order[1])]]
    elif(lorentztag=="VVSS") :
        couplings=[0.]
        terms=[['Metric(%s,%s)' % (order[0],order[1])]]
    elif(lorentztag=="VSS"):
         couplings=[0.,0.]
         terms=[['P(%s,%s)' % (order[0],order[2])],
                ['P(%s,%s)' % (order[0],order[1])]]
    # extract the lorentz structures
    structures = extractStructures(L)
    # handle the scalar couplings
    itype=-1
    try :
        for term in terms:
            itype+=1
            for perm in itertools.permutations(term):
                label = '*'.join(perm)
                for istruct in range(0,len(structures)) :
                    if label in structures[istruct] :
                        reminder = structures[istruct].replace(label,'1.',1)
                        couplings[itype]+=eval(reminder, {'cmath':cmath} )
                        structures[istruct]='Done'
    except :
        raise SkipThisVertex()
    # special for VVS and epsilon
    # handle the pseudoscalar couplings
    for struct in structures :
        if(struct != "Done" ) :
            if(lorentztag=="VVS") :
                VVSEpsilon(couplings,struct)
            else :
                raise SkipThisVertex()
    # put it all together
    if(len(all_couplings)==0) :
        for ic in range(0,len(couplings)) :
            if(couplings[ic]!=0.) :
                all_couplings.append('(%s) * (%s) * (%s)' % (prefactors,value,couplings[ic]))
            else :
                all_couplings.append(False)
    else :
        for ic in range(0,len(couplings)) :
            if(couplings[ic]!=0. and all_couplings[ic]) :
                all_couplings[ic] = '(%s) * (%s) * (%s) + (%s) ' % (prefactors,value,
                                                                    couplings[ic],all_couplings[ic])
            elif(couplings[ic]!=0) :
                all_couplings[ic] = '(%s) * (%s) * (%s) ' % (prefactors,value,couplings[ic])
    return all_couplings

def processScalarVectorCouplings(lorentztag,vertex,model,parmsubs,all_couplings,header,order) :
    # check the values
    tval = [False]*len(all_couplings[0])
    value =[False]*len(all_couplings[0])
    for icolor in range(0,len(all_couplings)) :
        for ix in range(0,len(all_couplings[icolor])) :
            if(not value[ix]) :
                value[ix] = all_couplings[icolor][ix]
            if(value[ix] and not tval[ix]) :
                tval[ix] = evaluate(value[ix],model,parmsubs)
            elif(value[ix]) :
                tval2 = evaluate(all_couplings[icolor][0],model,parmsubs)
                if(abs(tval[ix]-tval2)>1e-6) :
                    raise SkipThisVertex()

    append = ""
    symbols = set()
    coup_norm=0.
    if(lorentztag=="VVS") :
        if(not value[0] and not value[1] and not value[2] and \
           not value[3] and not value[4] and not value[6] and value[5]) :
            coup_norm=value[5]
        else :
            for ix in range(0,len(value)) :
                if(value[ix]) :
                    value[ix], sym = py2cpp(value[ix])
                    symbols |= sym
                else :
                    value[ix]=0.
            lorentztag = 'GeneralVVS'
            header="kinematics(true);"
            # g_mu,nv piece of coupling 
            if(value[5]!=0.) :
                append +='a00( %s + Complex(( %s )* GeV2/invariant(1,2)));\n' % ( value[0],value[5])
            else :
                append +='a00( %s );\n' % value[0]
            # other couplings
            append += 'a11( %s );\n    a12( %s );\n    a21( %s );\n    a22( %s );\n aEp( %s );\n' % \
                      ( value[1],value[2],value[3],value[4],value[6] )
            coup_norm="1."
    elif(lorentztag=="VVSS") :
        coup_norm = value[0]
    elif(lorentztag=="VSS") :
        if(abs(tval[0]+tval[1])>1e-6) :
            for ix in range(0,len(value)) :
                if(value[ix]) :
                    value[ix], sym = py2cpp(value[ix])
                    symbols |= sym
                else :
                    value[ix]=0.
            coup_norm = "1."
            append = 'if(p2->id()==%s) { a( %s ); b( %s);}\n else { a( %s ); b( %s);}' \
                     % (vertex.particles[order[1]-1].pdg_code,
                        value[0],value[1],value[1],value[0])
        else :
            coup_norm = value[1]
            append = 'if(p2->id()!=%s){norm(-norm());}' \
                     % vertex.particles[order[1]-1].pdg_code
    # return the answer
    return (coup_norm,append,lorentztag,header,symbols)

def getIndices(term) :
    if(term[0:2]=="P(") :
        indices = term.strip(")").strip("P(").split(",")
        mom   = int(indices[1])
        index = int(indices[0])
        return (True,mom,index)
    else :
        return (False,0,0)
    

def lorentzScalar(vertex,L) :
    dotProduct = """(invariant( i[{i1}], i[{i2}] )/GeV2)"""
    structures=L.structure.split()
    output="("
    for struct in structures:
        if(struct=="+" or struct=="-") :
            output+=struct
            continue
        structure = struct.split("*")
        worked = False
        mom=-1
        newTerm=""
        while True :
            term = structure[-1]
            structure.pop()
            (momentum,mom,index) = getIndices(term)
            if( not momentum) : break
            # look for the matching momenta
            for term in structure :
                (momentum,mom2,index2) = getIndices(term)
                if(index2==index) :
                    structure.remove(term)
                    dot = dotProduct.format(i1=mom-1,i2=mom2-1)
                    if(newTerm=="") :
                        newTerm = dot
                    else :
                        newTerm = " ( %s) * ( %s ) " % (newTerm,dot)
            if(len(structure)==0) :
                worked = True
                break
        if(not worked) :
            return False
        else :
            output+=newTerm
    output+=")"
    return output

kinematicsline = """\
long id [3]={{{id1},{id2},{id3}}};
    long id2[3]={{p1->id(),p2->id(),p3->id()}};
    unsigned int i[3];
    for(unsigned int ix=0;ix<3;++ix) {{
      for(unsigned int iy=0;iy<3;++iy) {{
	if(id[ix]==id2[iy]) {{
	  i[ix] = iy;
	  id2[iy]=0;
	  break;
	}}
      }}
    }}
    double hw_kine1 = {kine};
"""

kinematicsline2 = """\
long id [4]={{{id1},{id2},{id3},{id4}}};
    long id2[4]={{p1->id(),p2->id(),p3->id(),p4->id()}};
    unsigned int i[4];
    for(unsigned int ix=0;ix<4;++ix) {{
      for(unsigned int iy=0;iy<4;++iy) {{
	if(id[ix]==id2[iy]) {{
	  i[ix] = iy;
	  id2[iy]=0;
	  break;
	}}
      }}
    }}
    double hw_kine1 = {kine};
"""

kinematicsline3 ="""\
    double hw_kine{i} = {kine};
"""

def scalarCouplings(vertex,value,prefactors,L,lorentztag,
                    all_couplings,prepend,header) :
    try :
        val = int(L.structure)
    except :
        output = lorentzScalar(vertex,L)
        if( not output ) :
            raise SkipThisVertex()
        else :
            if(prepend=="") :
                if(lorentztag=="SSS") :
                    # order doesn't matter here, all same spin
                    prepend = kinematicsline.format(id1=vertex.particles[0].pdg_code,
                                                    id2=vertex.particles[1].pdg_code,
                                                    id3=vertex.particles[2].pdg_code,
                                                    kine=output)
                else :
                    # order doesn't matter here, all same spin
                    prepend = kinematicsline2.format(id1=vertex.particles[0].pdg_code,
                                                     id2=vertex.particles[1].pdg_code,
                                                     id3=vertex.particles[2].pdg_code,
                                                     id4=vertex.particles[3].pdg_code,
                                                     kine=output)
                value = "(%s) *(hw_kine1)" % value
            else :
                osplit=prepend.split("\n")
                i=-1
                while osplit[i]=="":
                    i=i-1
                ikin=int(osplit[i].split("=")[0].replace("double hw_kine",""))+1
                prepend +=kinematicsline3.format(kine=output,i=ikin)
                value = "(%s) *(hw_kine%s)" % (value,ikin)
            header="kinematics(true);"
    if(len(all_couplings)==0) :
        all_couplings.append('(%s) * (%s)' % (prefactors,value))
    else :
        all_couplings[0] = '(%s) * (%s) + (%s)' % (prefactors,value,all_couplings[0])
    return (prepend, header,all_couplings)

def processScalarCouplings(model,parmsubs,all_couplings) :
    tval = False
    value = False
    for icolor in range(0,len(all_couplings)) :
        if(len(all_couplings[icolor])!=1) :
            raise SkipThisVertex()
        if(not value) :
            value = all_couplings[icolor][0]
        m = re.findall('hw_kine[0-9]*',  all_couplings[icolor][0])
        if m:
            for kine in m:
                # bizarre number for checks, must be a better option
                parmsubs[kine] = 987654321.
        if(not tval) :
            tval = evaluate(value,model,parmsubs)
        else :
            tval2 = evaluate(all_couplings[icolor][0],model,parmsubs)
            if(abs(tval[i]-tval2)>1e-6) :
                raise SkipThisVertex()
    # cleanup and return the answer
    return value.replace("(1.0) * ","").replace(" * (1)","")

def vectorCouplings(vertex,value,prefactors,L,lorentztag,pos,
                    all_couplings,append,qcd,order) :
    structures=extractStructures(L)
    terms=[]
    signs=[]
    if(lorentztag=="VVV") :
        terms=[['P(%s,%s)' % (order[2],order[0]),'Metric(%s,%s)' % (order[0],order[1])],
               ['P(%s,%s)' % (order[2],order[1]),'Metric(%s,%s)' % (order[0],order[1])],
               ['P(%s,%s)' % (order[1],order[0]),'Metric(%s,%s)' % (order[0],order[2])],
               ['P(%s,%s)' % (order[1],order[2]),'Metric(%s,%s)' % (order[0],order[2])],
               ['P(%s,%s)' % (order[0],order[1]),'Metric(%s,%s)' % (order[1],order[2])],
               ['P(%s,%s)' % (order[0],order[2]),'Metric(%s,%s)' % (order[1],order[2])]]
        signs=[1.,-1.,-1.,1.,1.,-1.]
    elif(lorentztag=="VVVV") :
        terms=[['Metric(%s,%s)' % (order[0],order[3]),'Metric(%s,%s)' % (order[1],order[2])],
               ['Metric(%s,%s)' % (order[0],order[2]),'Metric(%s,%s)' % (order[1],order[3])],
               ['Metric(%s,%s)' % (order[0],order[1]),'Metric(%s,%s)' % (order[2],order[3])]]
        signs=[1.,1.,1.]
    elif(lorentztag=="VVVS") :
        terms=[['P(%s,%s)' % (order[2],order[0]),'Metric(%s,%s)' % (order[0],order[1])],
               ['P(%s,%s)' % (order[2],order[1]),'Metric(%s,%s)' % (order[0],order[1])],
               ['P(%s,%s)' % (order[1],order[0]),'Metric(%s,%s)' % (order[0],order[2])],
               ['P(%s,%s)' % (order[1],order[2]),'Metric(%s,%s)' % (order[0],order[2])],
               ['P(%s,%s)' % (order[0],order[1]),'Metric(%s,%s)' % (order[1],order[2])],
               ['P(%s,%s)' % (order[0],order[2]),'Metric(%s,%s)' % (order[1],order[2])],
               ['Epsilon(1,2,3,-1)','P(-1,1)'],['Epsilon(1,2,3,-1)','P(-1,2)'],
               ['Epsilon(1,2,3,-1)','P(-1,3)']]
        signs=[1.,-1.,-1.,1.,1.,-1.,1.,1.,1.]

    # extract the couplings
    new_couplings  = [False]*len(terms)
    iterm=0
    try :
        for term in terms:
            for perm in itertools.permutations(term):
                label = '*'.join(perm)
                for istruct in range(0,len(structures)) :
                    if label in structures[istruct] :
                        reminder = structures[istruct].replace(label,'1.',1)
                        structures[istruct] = "Done"
                        val = eval(reminder, {'cmath':cmath} )*signs[iterm]
                        if(new_couplings[iterm]) :
                            new_couplings[iterm] += val
                        else :
                            new_couplings[iterm] = val
            iterm += 1
    except :
        raise SkipThisVertex()
    # check we've handled all the terms
    for val in structures:
        if(val!="Done") :
            raise SkipThisVertex()
    # set the couplings
    for icoup in range(0,len(new_couplings)) :
        if(new_couplings[icoup]) :
            new_couplings[icoup] = '(%s) * (%s) *(%s)' % (new_couplings[icoup],prefactors,value)
    if(len(all_couplings)==0) :
        all_couplings=new_couplings
    else :
        for icoup in range(0,len(new_couplings)) :
            if(new_couplings[icoup] and all_couplings[icoup]) :
                all_couplings[icoup] = '(%s) * (%s) *(%s) + (%s) ' % (new_couplings[icoup],prefactors,value,all_couplings[icoup])
            elif(new_couplings[icoup]) :
                all_couplings[icoup] = new_couplings[icoup]
    # ordering for VVV type vertices
    if(len(pos[8]) != 3 and (lorentztag=="VVV" or lorentztag=="VVVS")) :
        append = VVVordering(vertex)
    return all_couplings,append

def processVectorCouplings(lorentztag,vertex,model,parmsubs,all_couplings,append,header) :
    value = False
    tval  = False
    if(lorentztag=="VVV") :
        for icolor in range(0,len(all_couplings)) :
            # loop over the different terms
            for ix in range(0,len(all_couplings[icolor])) :
                if(not value) :
                    value = all_couplings[icolor][ix]
                    tval = evaluate(value,model,parmsubs)
                else :
                    tval2 = evaluate(all_couplings[icolor][ix],model,parmsubs)
                    if(abs(tval-tval2)>1e-6) :
                        raise SkipThisVertex()
    elif(lorentztag=="VVVV") :
        order=[]
        colours = vertex.color
        if(len(colours)==1) :
            tval=[]
            for i in range(0,3) :
                tval.append(evaluate(all_couplings[0][i],model,parmsubs))
            if(compare(tval[2],-2.*tval[1]) and
               compare(tval[2],-2.*tval[0]) ) :
                order=[0,1,2,3]
                value = "0.5*(%s)" % all_couplings[0][2]
            elif(compare(tval[1],-2.*tval[2]) and
                 compare(tval[1],-2.*tval[0]) ) :
                order=[0,2,1,3]
                value = "0.5*(%s)" % all_couplings[0][1]
            elif(compare(tval[0],-2.*tval[2]) and
                 compare(tval[0],-2.*tval[1]) ) :
                order=[0,3,1,2]
                value = "0.5*(%s)" % all_couplings[0][0]
            else:
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
            # order doesn't matter here same spin
            append = pattern % ( vertex.particles[0].pdg_code,order[0],
                                 vertex.particles[1].pdg_code,order[1],
                                 vertex.particles[2].pdg_code,order[2],
                                 vertex.particles[3].pdg_code,order[3] )
        else :
            for icolor in range(0,len(all_couplings)) :
                col=colours[icolor].split("*")
                if(len(col)==2 and "f(" in col[0] and "f(" in col[1]) :
                    sign = 1
                    for i in range(0,2) :
                        col[i],stemp = extractAntiSymmetricIndices(col[i],"f(")
                        for ix in range(0,len(col[i])): col[i][ix]=int(col[i][ix])
                        sign *=stemp
                    if(col[0][0]>col[1][0]) : col[0],col[1] = col[1],col[0]
                    # first flow
                    if(col[0][0]==1 and col[0][1]==2 and col[1][0] ==3 and col[1][1] == 4) :
                        if(all_couplings[icolor][2] or not all_couplings[icolor][0] or
                           not all_couplings[icolor][1]) :
                            raise SkipThisVertex()
                        if(not value) :
                            value = all_couplings[icolor][1]
                            tval  = evaluate(value,model,parmsubs)
                        tval2 = -evaluate(all_couplings[icolor][0],model,parmsubs)
                        tval3 =  evaluate(all_couplings[icolor][1],model,parmsubs)
                    elif(col[0][0]==1 and col[0][1]==3 and col[1][0] ==2 and col[1][1] == 4) : 
                        if(all_couplings[icolor][1] or not all_couplings[icolor][0] or
                           not all_couplings[icolor][2]) :
                            raise SkipThisVertex()
                        if(not value) :
                            value = all_couplings[icolor][2]
                            tval  = evaluate(value,model,parmsubs)
                        tval2 = -evaluate(all_couplings[icolor][0],model,parmsubs)
                        tval3 =  evaluate(all_couplings[icolor][2],model,parmsubs)
                    elif(col[0][0]==1 and col[0][1]==4 and col[1][0] ==2 and col[1][1] == 3) : 
                        if(all_couplings[icolor][0] or not all_couplings[icolor][1] or
                           not all_couplings[icolor][2]) :
                            raise SkipThisVertex()
                        if(not value) :
                            value = all_couplings[icolor][2]
                            tval  = evaluate(value,model,parmsubs)
                        tval2 = -evaluate(all_couplings[icolor][1],model,parmsubs)
                        tval3 =  evaluate(all_couplings[icolor][2],model,parmsubs)
                    else :
                        raise SkipThisVertex()
                    if(abs(tval-tval2)>1e-6 or abs(tval-tval3)>1e-6 ) :
                        raise SkipThisVertex()
                    append = 'setType(1);\nsetOrder(0,1,2,3);'
                else :
                    print 'unknown colour structure for VVVV vertex'
                    raise SkipThisVertex()
    elif(lorentztag=="VVVS") :
        try :
            # two distinct cases 0-5 = , 6-8=
            if(all_couplings[0][0]) :
                imin=0
                imax=6
                header="scalar(true);"
            else :
                imin=6
                imax=9
                header="scalar(false);"
            for icolor in range(0,len(all_couplings)) :
                # loop over the different terms
                for ix in range(imin,imax) :
                    if(not value) :
                        value = all_couplings[icolor][ix]
                        tval = evaluate(value,model,parmsubs)
                    else :
                        tval2 = evaluate(value,model,parmsubs)
                        if(abs(tval-tval2)>1e-6) :
                            raise SkipThisVertex()
        except :
            SkipThisVertex()
    # cleanup and return the answer
    value = value.replace("(1.0) * ","").replace(" * (1)","")
    return (value,append,header)

def fermionCouplings(value,prefactors,L,all_couplings,order) :
    new_couplings=[False,False]
    try :
        new_couplings[0],new_couplings[1] = parse_lorentz(L.structure)
    except :
        raise SkipThisVertex()
    for i in range(0,2) :
        if new_couplings[i]:
            new_couplings[i] = '(%s) * (%s) * (%s)' % (prefactors,new_couplings[i],value)
    if(len(all_couplings)==0) :
        all_couplings=new_couplings
    else :
        for i in range(0,len(new_couplings)) :
            if(new_couplings[i] and all_couplings[i]) :
                all_couplings[i] = '(%s) + (%s) ' % (new_couplings[i],all_couplings[i])
            elif(new_couplings[i]) :
                all_couplings[i] = new_couplings[i]
    return all_couplings

def processFermionCouplings(lorentztag,vertex,model,parmsubs,all_couplings,order) :
    leftcontent  = all_couplings[0][0] if all_couplings[0][0] else "0."
    rightcontent = all_couplings[0][1] if all_couplings[0][1] else "0."
    tval=[evaluate( leftcontent,model,parmsubs),
          evaluate(rightcontent,model,parmsubs)]
    for icolor in range(0,len(all_couplings)) :
        # loop over the different terms
        for ix in range(0,len(all_couplings[icolor])) :
            tval2 = evaluate(all_couplings[icolor][ix],model,parmsubs) if all_couplings[icolor][ix] else 0.
            if(abs(tval[ix]-tval2)>1e-6) :
                raise SkipThisVertex()
    normcontent  = "1."
    append=""
    if lorentztag == 'FFV':
        append = ('if(p1->id()!=%s) {Complex ltemp=left(), rtemp=right(); left(-rtemp); right(-ltemp);}' 
                  % vertex.particles[order[0]-1].pdg_code)
    return normcontent,leftcontent,rightcontent,append

def RSCouplings(value,prefactors,L,all_couplings,order) :
    raise SkipThisVertex()

class LorentzIndex :
    """ A simple classs to store a Lorentz index """
    type=""
    value=0
    dimension=0
    def __repr__(self):
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
                print "INDEX",val
                quit()
        else :
            print 'unknown value in lorentz index',val
            quit()
            
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
            print 'index problem',a.name,b.name
            print a.spin,b.spin
            quit()
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
    # special handling for powers , assume only 2
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
                print "Problem handling %s structure " % output[-1].name
                raise SkipThisVertex()
        # objects with 2 lorentz indices
        elif(struct.find("Metric")==0) :
            output.append(LorentzStructure())
            output[-1].lorentz=[LorentzIndex(ind[0]),LorentzIndex(ind[1])]
            output[-1].name=struct.split("(")[0]
            output[-1].value=0
            output[-1].spin=[]
            if(len(struct.replace("%s(%s,%s)" % (output[-1].name,ind[0],ind[1]),""))!=0) :
                print "Problem handling %s structure " % output[-1].name
                raise SkipThisVertex()
        elif(struct.find("P(")==0) :
            output.append(LorentzStructure())
            output[-1].lorentz=[LorentzIndex(ind[0])]
            output[-1].name=struct.split("(")[0]
            output[-1].value=ind[1]
            output[-1].spin=[]
            if(len(struct.replace("%s(%s,%s)" % (output[-1].name,ind[0],ind[1]),""))!=0) :
                print "Problem handling %s structure " % output[-1].name
                raise SkipThisVertex()
        # 1 lorentz and 1 spin index
        elif(struct.find("Gamma")==0) :
            output.append(LorentzStructure())
            output[-1].lorentz=[LorentzIndex(ind[0])]
            output[-1].spin=[ind[1],ind[2]]
            output[-1].name=struct.split("(")[0]
            output[-1].value=1
            if(len(struct.replace("%s(%s,%s,%s)" % (output[-1].name,ind[0],ind[1],ind[2]),""))!=0) :
                print "problem parsing gamma matrix",struct
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
                print 'problem parsing epsilon',struct
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
                    print 'scalar problem',struct
                    print complex(struct)
                    quit()
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
       quit()
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
                    print 'unknown type'
                    print parsed[j]
                    quit()
                if(parsed[k].name=="P") :
                    newIndex2 = LorentzIndex(parsed[k].value)
                    newIndex2.type="P"
                    newIndex2.dimension=1
                elif(parsed[k].name=="Metric") :
                    for li in parsed[k].lorentz :
                        if(li != index) :
                            newIndex2=li
                            break
                else :
                    print 'unknown type'
                    print parsed[j]
                    quit()
                parsed[j].name = "Tensor"
                parsed[j].value= int(index.value)
                parsed[j].lorentz= [newIndex1,newIndex2]
                parsed[k]=""
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
            if(sc<3) :
                continue
            elif(sc==3) :
                offLoc = -1
                for i in summable:
                    if(not summable[i]) :
                        offLoc = i
                        break
            else :
                offLoc = 0
            indices=[]
            dTemp=0
            for ix in range(0,len(parsed[j].lorentz)) :
                if(parsed[j].lorentz[ix].type=="P") : dTemp+=1
                if(ix!=offLoc) : indices.append(parsed[j].lorentz[ix])
            dimension[2] += dTemp
            if(sc==4) :
                iTemp = (parsed[j].lorentz[0],parsed[j].lorentz[1],
                         parsed[j].lorentz[2],parsed[j].lorentz[3])
                if(iTemp in defns) :
                    output += "*(%s)" % defns[iTemp][0]
                    parsed[j]=""
                else :
                    name = "dot%s" % (len(defns)+1)
                    unit = computeUnit(dTemp)
                    print 'eps construct A',name,unit
                    defns[iTemp] = [name,"complex<%s> %s =-%s*epsilon(%s,%s,%s);" % (unit,name,parsed[j].lorentz[0],
                                                                                     indices[0],indices[1],indices[2]) ]
                    output += "*(%s)" % name
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
                newIndex = LorentzIndex(int(name[1]))
                newIndex.type="V"
                newIndex.dimension=dTemp
                output += "*(%s)" % (sign)
                oi = parsed[j].lorentz[offLoc]
                if(oi.type!="D") :
                    parsed[j].name="Metric"
                    parsed[j].spins=[]
                    parsed[j].value=0
                    parsed[j].lorentz=[newIndex,oi]
                    print 'replace with metric'
                    print parsed[j]
                    quit()
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
                        print "problem in eps"
                        quit()
                    parsed[j]=""
        elif(parsed[j].name=="Tensor") :
            # not an exteral tensor
            if(parsed[j].value!=iloc) :
                # now check the lorentz indices
                con=[]
                uncon=[]
                dtemp=0
                for li in parsed[j].lorentz :
                    if(li.type=="P") :
                        con.append(li)
                        dtemp+=1
                    elif(li.type=="E") :
                        if(li.value!=iloc) :
                            con.append(li)
                        else :
                            uncon.append(li)
                    else :
                        print 'need to handle ',li,'in tensor',parsed[j]
                        print li
                        quit()
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
                elif(len(con)==1 and len(uncon)==1) :
                    print 'uncon'
                    quit()
                else :
                    print "can't happen"
                    quit()
            else :
                dtemp=0
                for li in parsed[j].lorentz :
                    if(li.type=="P") : dtemp+=1
                dimension[2]+=dtemp
        elif(parsed[j].name.find("Proj")>=0 or
             parsed[j].name.find("Gamma")>=0) :
            continue
        elif(parsed[j].name=="P" and parsed[j].lorentz[0].type=="R") :
            continue
        else :
            print 'not handled',parsed[j],iloc
            quit()
    # remove leading *
    if(output!="" and output[0]=="*") : output = output[1:]
    # remove any (now) empty elements
    parsed = [x for x in parsed if x != ""]
    return (output,parsed,dimension,eps)

def finalContractions(output,parsed,dimension,lorentztag,iloc,defns) :
    if(len(parsed)==0) :
       return (output,dimension)
    elif(len(parsed)!=1) :
        print "summation can't be handled",parsed,iloc,output
        quit()
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
                unit="double"
                if(lo.type=="P") :
                    dimension[2]+=1
                    unit="Energy"
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
                print 'problem with tensor',iloc
                print parsed
                quit()
        # off-shell tensor
        else :
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
            if(lo.type=="P") :
                dimension[2]+=1
            if(output=="") : output="1."
            output = "(%s)*(%s)" %(output,lo)
    else :
        print "structure can't be handled",parsed,iloc
        quit()
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
        pre = "-2./3.*(1.-%s*OM%s)" % (dp,struct.value)
        terms.append((pre,i0,i0))
        pre = "-4./3.*(1.-%s*OM%s)" % (dp,struct.value)
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
        terms.append(("1",ind1,ind2))
        terms.append(("-OM%s*%s"%(struct.value,d1),ip,ind2))
        terms.append(("-OM%s*%s"%(struct.value,d2),ind1,ip))
        terms.append(("1",ind2,ind1))
        terms.append(("-OM%s*%s"%(struct.value,d2),ip,ind1))
        terms.append(("-OM%s*%s"%(struct.value,d1),ind2,ip))
        terms.append(("-2./3.*"+d3,i0,i0))
        terms.append(("2./3.*OM%s*%s*%s"%(struct.value,d1,d2),i0,i0))
        terms.append(("2./3.*OM%s*%s"%(struct.value,d3),ip,ip))
        terms.append(("4./3.*OM%s*OM%s*%s*%s"%(struct.value,struct.value,d1,d2),ip,ip))
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
                    val += "%s*%s.%s()*%s.%s()" % (pre,term[1],i1,term[2],i2)
            output["%s%s" % (i1,i2) ] = val.replace("+1*","+").replace("-1*","-")
    return output
        
def generateVertex(iloc,L,parsed,lorentztag,vertex,defns) :
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
        else :
            output=DiracMatrix()
            output.name="M"
            output.value = vslash.substitute({ "v" : structure.lorentz[0]})
            Symbols += vslashS.substitute({ "v" : structure.lorentz[0]})
            variable = computeUnit(structure.lorentz[0].dimension)
            if(structure.lorentz[0].type!="V") :
                dtemp[2] += structure.lorentz[0].dimension
            defns["vv%s" % structure.lorentz[0] ] = \
                                                    ["vv%s" % structure.lorentz[0],
                                                     vslashD.substitute({ "var" : variable,
                                                                          "v" : structure.lorentz[0]})]
            structure=""
    else :
        print 'Unknown Gamma matrix structure',structure
        quit()
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
            print " CONT TEST A",parsed
            print 'unkonwn type',parsed[i]
            quit()
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
            print "Can't parsed the chain of dirac matrices"
            quit()
            SkipThisVertex()
    # start and end of the spin chains
    # easy case, both spin 1/2
    if(spins[sind-1]==2 and spins[lind-1]==2) :
        start = DiracMatrix()
        end   = DiracMatrix()
        start.index=0
        end  .index=0
        start.name="S"
        end  .name="S"
        # start of chain
        # off shell
        if(sind==iloc) :
            start.name="M"
            start.value = vslashM.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
            Symbols += vslashMS.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
            defns["vvP%s" % sind ] = ["vvP%s" % sind ,
                                      vslashD.substitute({ "var" : "Energy",
                                                    "v" :  "P%s" % sind })]
            dtemp[1]+=1
        # onshell
        else :
            subs = {'s' : ("sbar%s" % sind)}
            start.value = sbar .substitute(subs)
            Symbols += sline.substitute(subs)
        # end of chain
        if(lind==iloc) :
            end.name="M"
            end.value = vslashM2.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
            Symbols += vslashMS.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
            defns["vvP%s" % lind ] = ["vvP%s" % lind ,
                                      vslashD.substitute({ "var" : "Energy",
                                                           "v" :  "P%s" % lind })]
            dtemp[1] += 1
        else :
            subs = {'s' : ("s%s" % lind)}
            end.value = spinor.substitute(subs)
            Symbols += sline.substitute(subs)
        startT = start
        endT   = end
    # start 1/2 and end 3/2
    elif spins[sind-1]==2 and spins[lind-1]==4 :
        start = DiracMatrix()
        endT  = DiracMatrix()
        start.index=0
        endT .index=0
        # spin 1/2 fermion
        if(sind==iloc) :
            start.value  = vslashM.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
            Symbols     += vslashMS.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
            endT.value   = vslashM2.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
            defns["vvP%s" % sind ] = ["vvP%s" % sind ,
                                      vslashD.substitute({ "var" : "Energy",
                                                           "v" :  "P%s" % sind })]
            dtemp[1]+=1
            start.name="M"
            endT .name="M"
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
    # start 3/2 and end 1/2
    elif spins[sind-1]==4 and spins[lind-1]==2 :
        # spin 1/2 fermion
        end    = DiracMatrix()
        startT = DiracMatrix()
        if(lind==iloc) :
            end.value    = vslashM2.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
            startT.value =  vslashM.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
            Symbols += vslashMS.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
            defns["vvP%s" % lind ] = ["vvP%s" % lind ,
                                      vslashD.substitute({ "var" : "Energy",
                                                           "v" :  "P%s" % lind })]
            dtemp[1] += 1
            startT.name="M"
            end   .name="M"
        else :
            subs = {'s' : ("s%s" % lind)}
            end.value  = spinor.substitute(subs)
            Symbols += sline.substitute(subs)
            subs = {'s' : ("sbar%s" % lind)}
            startT.value = sbar .substitute(subs)
            Symbols += sline.substitute(subs)
            startT.name="S"
            end   .name="S"
        # spin 3/2 fermion
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
                    contract = LorentzIndex(lind)
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
    else :
        print 'At most one R-S spinor allowed in a vertex'
        raise SkipThisVertex()
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
             ( "%s*%s*Transpose(%s)*%s*%s" %(startT[0],CC,etemp,CD,endT[0]))) in tempT
        res.append(tempT["result"])
    if(len(start)==1) :
        if(iloc==0 or (iloc!=sind[0] and iloc!=lind[0])) :
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
                    print "SIZE PROBLEM A",sign,rTemp[ii].shape
                    quit()
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
                print "SIZE PROBLEM B ",sign,rTemp[ii][0].shape
                quit()
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
        # external index
        elif key.type == "O" :
            continue
        # uncontracted vector index
        elif key.type=="E" or key.type=="Q":
            continue
        else :
            print 'uknown type of uncontracted index',key
            quit()
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
                    print 'unknown index',index
                    quit()
        else :
            print 'unknown lorentz object',lstruct,iloc
            print lstruct.value
            print lorentz
            print expr
            print start
            print end
            
            quit()
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
            for key,val in contracted.iteritems() :
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
                            print 'unknown index'
                            print val
                            print val.name
                            print val.value
                            print val.index
                            quit()
                    # unknown to be sorted out
                    else :
                        print 'unkonwn type in expr'
                        print val
                        print val.name
                        print val.value
                        print val.index
                        quit()
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
                    print 'barred spinor not a spinor'
                    print start[ichain].name
                    print start[ichain].value
                    print start[ichain].index
                    quit()
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
                    print 'barred spinorT not a spinor'
                    print startT[ichain].name
                    print startT[ichain].value
                    print startT[ichain].index
                    quit()
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
                    print 'spinor not a spinor'
                    print end[ichain].name
                    print end[ichain].value
                    print end[ichain].index
                    quit()
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
                    print 'spinorT not a spinor'
                    print endT[ichain].name
                    print endT[ichain].value
                    print endT[ichain].index
                    quit()
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
                                print 'unknown index in eps',index
                                quit()
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
                # unknown
                else :
                    print 'unknown expression in lorentz loop'
                    print li.name
                    print li.value
                    quit()
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
            for (key,val) in contracted.iteritems() :
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
                   sVal["s"   ] = res[0][0]
                   sVal["sT2" ] = res[1][0]
                   sVal["sT1" ] = res[2][0]
                   sVal["sT12"] = res[3][0]
            # spinor
            elif(len(res[0])==4) :
                for k in range(0,4) :
                    sVal[ "s%s"  % (k+1) ] = res[0][k]
                    sVal[ "sT%s" % (k+1) ] = res[1][k]
            else :
                print 'summ problem',len(res),len(res[0])
                quit()
            break
        # uncontracted indices
        else :
            istring = ""
            for (key,val) in unContracted.iteritems() :
                istring +=imap[val]
            if(len(istring)!=1) :
                print 'index problem',istring
                print unContracted
                print unI
                quit()
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
            for (key,val) in unContracted.iteritems() :
                if(val==0) : nZero+=1
            if(nZero==len(unContracted)) : break
    # handle the vector case
    if( "t" in sVal ) :
        # deal with pure vectors
        if(len(sVal)==8 and "t" in sVal and len(sVal["t"])==1) :
            pass
        # RS spinors
        elif(len(sVal)==8 and "t" in sVal and len(sVal["t"])==4) :
            pass
        else :
            print sVal
            print 'val problrm A',len(sVal)
            quit()
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
                    if(lo.type=="P") : dimension[2]+=1
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
            else :
                print 'unknown lorentz structure',parsed[i]
                print parsed
                print expr
                print start
                print end
                quit()
                
        parsed = [x for x in parsed if x != ""]
        if(len(parsed)!=0) :
            print expr
            print "Can't parse ",parsed,iloc
            quit()
            raise SkipThisVertex()
    sVal ={}
    # deal with the simplest case first
    if len(unContracted) == 0 :
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
    dimension = list(map(lambda x, y: x + y, dtemp, dimension))
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
    # special for 4 fermion case
    if(vertex.lorentz[0].spins.count(2)==4) :
        return swapOrderFFFF(vertex,iloc,momenta,fIndex)

    
    names=['','sca','sp','v']
    waves=['','sca',''  ,'E']
    output=""
    for i in range(1,4) :
        ns = vertex.lorentz[0].spins.count(i)
        if((ns<=1 and i!=2) or (ns<=2 and i==2)) : continue
        if(i!=3 and i!=1) :
            print 'swap problem',i
            quit()
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

def swapOrderFFFF(vertex,iloc,momenta,fIndex) :
    output=""
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
        if(momenta[fIndex[j]-1][0] or momenta[fIndex[j+2]-1][0]) :
            momenta[fIndex[j]-1][0] = True
            momenta[fIndex[j+2]-1][0] = True
            output += "        swap(P%s,P%s);\n" % (fIndex[j],fIndex[j+2])
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
                print 'unknown spin',spin
                quit()
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
                print 'unknown spin',spin
                quit()
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
        print vals
        print 'spin problem',ispin
        quit()
    expr=[otype]
    for i in range(0,nf-1) :
        expr.append(copy.copy(otype))
    # dimension for checking
    dimCheck=dim[0]
    for i in range(0,len(vals)) :
        # simple signs
        if(vals[i]=="+" or vals[i]=="-") :
            for ii in range(0,len(expr)) :
                for(key,val) in expr[ii].iteritems() :
                    expr[ii][key] = expr[ii][key]+vals[i]
            continue
        # check the dimensions
        if(dimCheck[0]!=dim[i][0] or dimCheck[1]!=dim[i][1] or
           dimCheck[2]!=dim[i][2]) :
            print vertex.lorentz
            for j in range(0,len(vals)) :
                print dim[j],vals[j]
            print "DIMENSION PROBLEM",i,dimCheck,dim,vertex
            quit()
        # simplest case 
        if(isinstance(vals[i], basestring)) :
            for ii in range(0,len(expr)) :
                for(key,val) in expr[ii].iteritems() :
                    expr[ii][key] = expr[ii][key]+vals[i]
            continue
        # more complex structures
        pre = vals[i][0]
        if(pre=="(1.0)") : pre=""
        if(not isinstance(vals[i][1],dict)) :
            print 'must be a dict here'
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
            print vals[i]
            print 'problem with type'
            quit()
    # no of particles in the vertex
    vDim = len(vertex.lorentz[0].spins)
    # compute the unit and apply it
    unit = computeUnit2(dimCheck,vDim)
    if(unit!="") :
        for ii in range(0,len(expr)) :
            for (key,val) in expr[ii].iteritems() :
                expr[ii][key] = "(%s)*(%s)" % (val,unit)
    return expr

def combineComponents(result,offType,RS) :
    for i in range(0,len(result)) :
        for (key,value) in result[i].iteritems() :
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
            for (key,val) in result[ii].iteritems() :
                subs[ii]["out%s" % key]= val
    # spinors
    if("s1" in result[0]) :
        stype  = "LorentzSpinor"
        sbtype = "LorentzSpinorBar"
        if(offType.find("Bar")>0) : (stype,sbtype)=(sbtype,stype)
        if(not RS) : sbtype=stype
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
        print subs[0]
        print result
        print 'type not implemented'
        quit()
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
                for (key,val) in expr[ii].iteritems() :
                    result[ii][key] = " %s((local_C%s)*(%s)) " % (vtype,j,val)
        else :
            for ii in range(0,len(expr)) :
                for (key,val) in expr[ii].iteritems(): 
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
        if(FM and not RS) :
            if(nf!=4) :
                result = [vTemplateT.format(iloc=fIndex[1],id=vertex.particles[fIndex[1]-1].pdg_code,
                                              cf=py2cpp(cf)[0],res=result[0],resT=result[1])]
            else :
                result = [vTemplate4.format(iloc1=fIndex[1],iloc2=fIndex[3],
                                            id1=vertex.particles[fIndex[1]-1].pdg_code,
                                            id2=vertex.particles[fIndex[3]-1].pdg_code,
                                            cf=py2cpp(cf)[0],res1=result[0],res2=result[1],res3=result[2],res4=result[3])]
        else :
            result[0] = "return (%s)*(%s);\n" % (result[0],py2cpp(cf)[0])
        if(RS) :
            result[1] = "return (%s)*(%s);\n" % (result[1],py2cpp(cf)[0])
    # off-shell particle
    else :
        # off-shell scalar
        if(vertex.lorentz[0].spins[iloc-1] == 1 ) :
            if(RS) :
                result[1] = scaTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],res=result[1])
            if(FM and not RS) :
                result[0] = scaTTemplate.format(iloc=iloc,res=result[0],resT=result[1],isp=fIndex[1],
                                                id=vertex.particles[fIndex[1]-1].pdg_code,cf=py2cpp(cf[0])[0])
            else :
                result[0] = scaTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],res=result[0])
        # off-shell fermion
        elif(vertex.lorentz[0].spins[iloc-1] == 2 ) :
            if(RS) :
                if(offType.find("Bar")>0) :
                    offTypeT=offType.replace("Bar","")
                else :
                    offTypeT=offType.replace("Spinor","SpinorBar")
                result[1] = sTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offTypeT.replace("WaveFunction",""),
                                             res=result[1].replace( "M%s" % iloc, "mass" ),offTypeB=offTypeT)
            if(FM and not RS) :
                result = [sTTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offType.replace("WaveFunction",""),
                                            res=result[0].replace( "M%s" % iloc, "mass" ),resT=result[1].replace( "M%s" % iloc, "mass" ),
                                            offTypeB=offType,id=vertex.particles[iloc-1].pdg_code)]
            else :
                result[0] = sTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offType.replace("WaveFunction",""),
                                             res=result[0].replace( "M%s" % iloc, "mass" ),offTypeB=offType)
        # off-shell vector
        elif(vertex.lorentz[0].spins[iloc-1] == 3 ) :
            if(RS) :
                result[1] = vecTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],res=result[1])
            if(FM and not RS) :
                result[0] = vecTTemplate.format(iloc=iloc,res=result[0],resT=result[1],isp=fIndex[1],
                                                id=vertex.particles[fIndex[1]-1].pdg_code,
                                                cf=py2cpp(cf[0])[0])
            else :
                result[0] = vecTemplate.format(iloc=iloc,res=result[0],cf=py2cpp(cf)[0])
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
                print "RS spinors and tensors not handled"
                quit()
            if(FM) :
                result = [tenTTemplate.format(iloc=iloc,res=result[0],resT=result[1],isp=fIndex[1],
                                              id=vertex.particles[fIndex[1]-1].pdg_code,cf=py2cpp(cf[0])[0])]
            else :
                result = [tenTemplate.format(iloc=iloc,cf=py2cpp(cf)[0],res=result[0])]
        else :
            print 'unknown spin for off-shell particle',vertex.lorentz[0].spins[iloc-1]
            quit()
    # check if momenta defns needed to clean up compile of code
    for (key,val) in defns.iteritems() :
        if( isinstance(key, basestring)) :
            if(key.find("vvP")==0) :
                momenta[int(key[3])-1][0] = True
        else :
            for vals in key :
                if(vals.type=="P") :
                    momenta[vals.value-1][0] = True
    # cat the definitions
    defString=""
    for (key,value) in defns.iteritems() :
        if(value[0]=="") : continue
        if(value[0][0]=="V") :
            defString+="    %s\n" %value[1]
    for (key,value) in defns.iteritems() :
        if(value[0]=="") : continue
        if(value[0][0]!="V") :
            defString+="    %s\n" %value[1]
    sorder=swapOrder(vertex,iloc,momenta,fIndex)
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

    # special for transpose in the RS case
    if(FM and RS) :
        htemp = header.split(",")
        irs=-1
        isp=-1
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
        hnew = ','.join(htemp)
        header += ";\n" + hnew
        hnew = hnew.replace("virtual ","").replace("=-GeV","")
        for i in range(0,len(waves)) :
            if(waves[i].find("Spinor")<0) : continue
            if(waves[i].find("Bar")>0) :
                waves[i] = waves[i].replace("Bar","").replace("bar","")
            else :
                waves[i] = waves[i].replace("Spinor","SpinorBar").replace(" s"," sbar").replace("Rs","Rsbar")
        momentastring=""
        for i in range(0,len(momenta)) :
            if(momenta[i][0] and momenta[i][1]!="")  :
                if(momenta[i][1].find("barW")>0) :
                    momentastring+=momenta[i][1].replace("barW","W")+"\n   "
                elif(momenta[i][1].find("sW")>0) :
                    momentastring+=momenta[i][1].replace("sW","sbarW")+"\n   "
                else :
                    momentastring+=momenta[i][1]+"\n    "
            
        fnew = evaluateTemplate.format(decl=hnew,momenta=momentastring,defns=defString,
                                       waves="\n    ".join(waves),symbols='\n    '.join(symboldefs),
                                       couplings="\n    ".join(localCouplings),
                                       result=result[1],swap=sorder)
        function +="\n" + fnew
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
        print 'testing evaluate multiple problem',spin
        quit()
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
            call = defns[iloc].replace("virtual","").replace("ScalarWaveFunction","").replace("SpinorWaveFunction","") \
                              .replace("SpinorBarWaveFunction","").replace("VectorWaveFunction","").replace("TensorWaveFunction","") \
                              .replace("Energy2","q2").replace("int","").replace("complex<Energy>","").replace("Energy","").replace("=-GeV","") \
                              .replace("const  &","").replace("tcPDPtr","").replace("  "," ")
            if(iloc!=1) :
                call = call.replace(waves[0],waves[iloc-1])
            pdgid = vertex.particles[i].pdg_code
            code += "   %sif(out->id()==%s) return %s;\n" % (el,pdgid,call)
            iloc+=1
    code+="   else assert(false);\n"
    return (header,evaluateMultiple.format(decl=ccdefn,code=code))
