import itertools,cmath,re,sys,copy
from .helpers import SkipThisVertex,extractAntiSymmetricIndices,def_from_model
from .converter import py2cpp
from .lorentzparser import parse_lorentz
import sympy,string
from string import Template
from sympy import Matrix,Symbol


imap=["t","x","y","z"]

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
    }}
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
    InvEnergy2 OM{iloc} = mass.real()==ZERO ? InvEnergy2(ZERO) : 1./sqr(mass); 
    Energy2 p2 = P{iloc}.m2();
    Complex fact = Complex(0.,1.)*({cf})*propagator(iopt,p2,out,mass,width);
    LorentzTensor<double> output = sbarW{isp}.id()=={id} ? fact*({res}) : fact*({resT});
    return TensorWaveFunction(P{iloc},out,output);
"""

# various strings for matrixes
I4 = "Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])"
G5 = "Matrix([[-1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,1]])"
PM = "Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]])"
PP = "Matrix([[0,0,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,1]])"


vslash  = Template("Matrix([[0,0,${v}tmz,-${v}xmy],[0,0,-${v}xpy,${v}tpz],[${v}tpz,${v}xmy,0,0],[${v}xpy,${v}tmz,0,0]])")
vslashS = Template("${v}tpz=Symbol(\"${v}tpz\")\n${v}tmz=Symbol(\"${v}tmz\")\n${v}xpy=Symbol(\"${v}xpy\")\n${v}xmy=Symbol(\"${v}xmy\")\n")
momCom  = Template("${v}t = Symbol(\"${v}t\")\n${v}x = Symbol(\"${v}x\")\n${v}y = Symbol(\"${v}y\")\n${v}z = Symbol(\"${v}z\")\n")
vslashD = Template("complex<${var}> ${v}tpz = ${v}.t()+${v}.z();\n    complex<${var}> ${v}tmz = ${v}.t()-${v}.z();\n    complex<${var}> ${v}xpy = ${v}.x()+Complex(0.,1.)*${v}.y();\n    complex<${var}> ${v}xmy = ${v}.x()-Complex(0.,1.)*${v}.y();")
vslashM  = Template("Matrix([[$m,0,${v}tmz,-${v}xmy],[0,$m,-${v}xpy,${v}tpz],[${v}tpz,${v}xmy,$m,0],[${v}xpy,${v}tmz,0,$m]])")
vslashM2 = Template("Matrix([[$m,0,-${v}tmz,${v}xmy],[0,$m,${v}xpy,-${v}tpz],[-${v}tpz,-${v}xmy,$m,0],[-${v}xpy,-${v}tmz,0,$m]])")
vslashMS = Template("${v}tpz=Symbol(\"${v}tpz\")\n${v}tmz=Symbol(\"${v}tmz\")\n${v}xpy=Symbol(\"${v}xpy\")\n${v}xmy=Symbol(\"${v}xmy\")\n${m}=Symbol(\"${m}\")\nO${m}=Symbol(\"O${m}\")\n")

rslash   = Template("Matrix([[$m,0,${v}tmz,-${v}xmy],[0,$m,-${v}xpy,${v}tpz],[${v}tpz,${v}xmy,$m,0],[${v}xpy,${v}tmz,0,$m]])*( (ETA(B!,A!)-2/3*O${m}**2*${v}A!*${v}B!)*Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]) -1/3*R${loc}A*R${loc}B -1/3*O${m}*(${v}B!*R${loc}A-${v}A!*R${loc}B))")
rslash2  = Template("Matrix([[$m,0,-${v}tmz,${v}xmy],[0,$m,${v}xpy,-${v}tpz],[-${v}tpz,-${v}xmy,$m,0],[-${v}xpy,-${v}tmz,0,$m]])*( (ETA(B!,A!)-2/3*O${m}**2*${v}B!*${v}A!)*Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]) -1/3*R${loc}B*R${loc}A +1/3*O${m}*(${v}A!*R${loc}B-${v}B!*R${loc}A))")

dirac=["Matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])","Matrix([[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]])",
       "Matrix([[0,0,0,complex(0, -1)],[0,0,complex(0, 1),0],[0,complex(0, 1),0,0],[complex(0, -1),0,0,0]])",
       "Matrix([[0,0,1,0],[0,0,0,-1],[-1,0,0,0],[0,1,0,0]])"]
CC = "Matrix([[0,1,0,0],[-1,0,0,0],[0,0,0,-1],[0,0,1,0]])"
CD = "Matrix([[0,-1,0,0],[1,0,0,0],[0,0,0,1],[0,0,-1,0]])"

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
    def __repr__(self):
        return "%s%s" % (self.type,self.value)

    def __init__(self,val) :
        if(isinstance(val,int)) :
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
        return ( (self.type, self.value) 
                 == (other.type, other.value) )

    def __hash__(self) :
        return hash((self.type,self.value))

# def indexEqual(i1,i2):
#     return ( (i1.type, i1.value) == (i2.type, i2.value) )

class LorentzStructure:
    """A simple class to store a Lorentz structures"""
    name=""
    value=0
    lorentz=[]
    spin=[]
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
        return int(b.value-a.value)
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

def parse_structure(structure) :
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
            return output
        # simple numeric pre/post factors
        elif((structure[0]=="-" or structure[0]=="+") and
             structure[-1]==")" and structure[1]=="(") :
            output.append(LorentzStructure())
            output[-1].name="int"
            output[-1].value=structure[0]+"1."
            output[-1].value=float(output[-1].value)
            structure=structure[2:-1]
            found=True
        elif(structure[0]=="(") :
            temp=structure.rsplit(")",1)
            structure=temp[0][1:]
            output.append(LorentzStructure())
            output[-1].name="int"
            output[-1].value="1."+temp[1]
            output[-1].value=float(eval(output[-1].value))
            found=True
        elif(structure[0:2]=="-(") :
            temp=structure.rsplit(")",1)
            structure=temp[0][2:]
            output.append(LorentzStructure())
            output[-1].name="int"
            output[-1].value="-1."+temp[1]
            output[-1].value=float(eval(output[-1].value))
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
            print 'found proj'
            output.append(LorentzStructure())
            output[-1].spin=ind
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
            if(len(struct.replace("%s(%s,%s)" % (output[-1].name,ind[0],ind[1]),""))!=0) :
                print "Problem handling %s structure " % output[-1].name
                raise SkipThisVertex()
        elif(struct.find("P(")==0) :
            output.append(LorentzStructure())
            output[-1].lorentz=[LorentzIndex(ind[0])]
            output[-1].name=struct.split("(")[0]
            output[-1].value=ind[1]
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
            print 'eps'
            quit()
            output.append(LorentzStructure())
            output[-1].lorentz=ind
            output[-1].name=struct.split("(")[0]
            output[-1].value=1.
            if(len(struct.replace("%s(%s,%s,%s,%s)" % (output[-1].name,ind[0],ind[1],ind[2],ind[3]),""))!=0) :
                print "problem D"
                quit()
        # scalars
        else :
            try :
                output.append(LorentzStructure())
                output[-1].value=float(struct)
                output[-1].name="int"
            except :
                if(struct.find("complex")==0) :
                    vals = struct[0:-1].replace("complex(","").split(",")
                    output[-1].value=complex(float(vals[0]),float(vals[1]))
                    output[-1].name="int"
                else :
                    print 'scalar problem',struct
                    print complex(struct)
                    quit()
    # now do the sorting
    if(len(output)==1) : return output
    output = sorted(output,cmp=LorentzCompare)
    print output
    return output

def contract(parsed) :
    print "!!!!!!!!!!!!!!!!!!!!!!!!! IN CONTRACT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print parsed
    for j in range(0,len(parsed)) :
        if(parsed[j]=="") : continue
        if(parsed[j].name=="P") :
            # simplest case
            if(parsed[j].lorentz[0].type=="E" or
               parsed[j].lorentz[0].type=="P") :
                newIndex = LorentzIndex(parsed[j].value)
                newIndex.type="P"
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
                        if(parsed[k].name=="P") :
                            parsed[k].lorentz.append(LorentzIndex(parsed[k].value))
                            parsed[k].lorentz[1].type="P"
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
                expr = "UnitRemoval::InvE"
            else :
                expr = "UnitRemoval::InvE%s" % removal
        else :
            if(removal==-1) :
                expr = "UnitRemoval::E"
            else :
                expr = "UnitRemoval::E%s" % (-removal)
    if(output=="") : return expr
    elif(expr=="") : return output
    else           : return "%s*%s" %(output,expr)

# order the indices of a dot product
def indSort(a,b) :
    print 'in sort',a,b
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


def finishParsing(parsed,dimension,lorentztag,iloc,defns) :
    output=1.
    # replace signs
    if(len(parsed)==1 and parsed[0].name=="sign") :
        if(parsed[0].value>0) :
            output="+"
        else :
            output="-"
        parsed=[]
        return (output,parsed,dimension)
    # replace integers (really lorentz scalars)
    print parsed
    for j in range(0,len(parsed)) :
        if(parsed[j]!="" and parsed[j].name=="int") :
            print output,parsed[j]
            output *= parsed[j].value
            parsed[j]=""
    # bracket this for safety
    if(output!="") : output = "(%s)" % output
    # special for tensor indices
    if("T" in lorentztag) :
        for j in range(0,len(parsed)) :
            if(parsed[j]=="") :continue
            print "!!!!!!!!!",parsed[j]
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
            print index,index2
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
                print parsed[k]
                found = False
                for li in parsed[k].lorentz :
                    if(li == index2) : 
                        found=True
                        break
                if(not found) : continue
                if(parsed[j].name=="P") :
                    newIndex1 = LorentzIndex(parsed[j].value)
                    newIndex1.type="P"
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
            print 'before sort',parsed[j].lorentz[0],parsed[j].lorentz[1]
            (ind1,ind2) = sorted((parsed[j].lorentz[0],parsed[j].lorentz[1]),cmp=indSort)
            print "after sort ",ind1,ind2
            # this product already dealt with ?
            if((ind1,ind2) in defns) :
                output += "*(%s)" % defns[(ind1,ind2)][0]
                parsed[j]=""
                if(ind1.type=="P") : dimension[2] +=1
                if(ind2.type=="P") : dimension[2] +=1
                continue
            # handle the product
            name = "dot%s" % (len(defns)+1)
            if(ind1.type=="P") :
                # dot product of two momenta
                if(ind2.type=="P") :
                    dimension[2] +=2
                    defns[(ind1,ind2)] = [name,"Energy2 %s = %s*%s;" % (name,ind1,ind2)]
                    output += "*(%s)" % name
                    parsed[j]=""
                elif(ind2.type=="E") :
                    if(ind2.value!=iloc) :
                        dimension[2] += 1
                        defns[(ind1,ind2)] = [name,"complex<Energy> %s = %s*%s;" % (name,ind1,ind2)]
                        output += "*(%s)" % name
                        parsed[j]=""
        #                 
        #                 if(int(ind2[1])==iloc) :
        #                     output += "*(%s)" % ind1
        #                 else :
        #                     
        #                     
            elif(ind1.type=="E") :
                if(ind2.type!="E") :
                    print "EE problem",ind1,ind2
                    quit()
                elif(ind1.value!=iloc and ind2.value!=iloc) :
                    defns[(ind1,ind2)] = [name,"complex<double> %s = %s*%s;" % (name,ind1,ind2)]
                    output += "*(%s)" % name
                    parsed[j]=""
                    
        #             if(int(ind1[1])==iloc) :
        #                 output += "*(%s)" % ind2
        #             elif(int(ind2[1])==iloc) :
        #                 output += "*(%s)" % ind1
        #             else :
        #                 
        #                 
        #     elif(parsed[j].name=="Epsilon") :
        #         if(not eps) : eps = True
        #         offLoc = -1
        #         indices=[]
        #         dTemp=0
        #         for ix in range(0,len(parsed[j].lorentz)) :
        #             if(isinstance(parsed[j].lorentz[ix],int)) :
        #                 offLoc = ix
        #                 break
        #             elif(parsed[j].lorentz[ix][0]=="E" and int(parsed[j].lorentz[ix][1])==iloc ) :
        #                 offLoc = ix
        #                 break
        #         for ix in range(0,len(parsed[j].lorentz)) :
        #             if(isinstance(parsed[j].lorentz[ix],basestring) and
        #                parsed[j].lorentz[ix][0]=="P") : dTemp+=1
        #             if((offLoc<0 and ix != 0) or
        #                (offLoc>=0 and offLoc!=ix) ) :
        #                 indices.append(parsed[j].lorentz[ix])
        #         dimension[2] += dTemp
        #         if(offLoc<0) :
        #             iTemp = (parsed[j].lorentz[0],parsed[j].lorentz[1],
        #                      parsed[j].lorentz[2],parsed[j].lorentz[3])
        #             if(iTemp in defns) :
        #                 output += "*(%s)" % defns[iTemp][0]
        #                 parsed[j]=""
        #             else :
        #                 name = "dot%s" % (len(defns)+1)
        #                 unit = computeUnit(dTemp)
        #                 defns[iTemp] = [name,"complex<%s> %s =-%s*epsilon(%s,%s,%s);" % (unit,name,parsed[j].lorentz[0],
        #                                                                                 indices[0],indices[1],indices[2]) ]
        #                 output += "*(%s)" % name
        #         else :
        #             iTemp = (indices[0],indices[1],indices[2])
        #             sign = ""
        #             if(offLoc%2!=0) : sign="-"
        #             print iTemp
        #             if(iTemp in defns) :
        #                 output += "*(%s%s)" % (sign,defns[iTemp][0])
        #                 parsed[j]=""
        #             else :
        #                 name = "vec%s" % (len(defns)+1)
        #                 output += "*(%s%s)" % (sign,name)
        #                 unit = computeUnit(dTemp)
        #                 defns[iTemp] = [name,"LorentzVector<complex<%s> > %s =-epsilon(%s,%s,%s);" % (unit,name,
        #                                                                                              indices[0],indices[1],indices[2]) ]
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
                else :
                    print "can't happen"
                    quit()
                
                print 'tensor  on-shell'
            else :
                print 'off-shell tensor'
                dtemp=0
                for li in parsed[j].lorentz :
                    if(li.type=="P") : dtemp+=1
                dimension[2]+=dtemp
        elif(parsed[j].name.find("Proj")>=0 or
             parsed[j].name.find("Gamma")>=0) :
            print 'skip gamma'
            continue
        else :
            print 'not handled',parsed[j]
            quit()
    # remove leading *
    if(output!="" and output[0]=="*") : output = output[1:]
    # remove any (now) empty elements
    parsed = [x for x in parsed if x != ""]
    return (output,parsed,dimension)

def finalContractions(output,parsed,dimension,lorentztag,iloc,defns) :
    if(len(parsed)==0) :
       return (output,dimension)
    elif(len(parsed)!=1) :
        print "summation can't be handled",parsed
        raise skipThisVertex()
    if(parsed[0].name=="Tensor") :
        tensor = tensorPropagator(parsed[0],defns)
        if(output=="") : output="1."
        output = [output,tensor,()]
    else :
        print "structure can't be handled",parsed
        raise skipThisVertex()
    return (output,dimension)
    
def tensorPropagator(struct,defns) :
    # dummy index
    i0 = LorentzIndex(-1000)
    # index for momentum of propagator
    ip = LorentzIndex(struct.value)
    ip.type="P"
    # the metric tensor
    terms=[]
    if(len(struct.lorentz)==0) :
        iTemp=(ip,ip)
        if(iTemp in defns) :
            dp = defns[iTemp][0]
        else :
            dp = "dot%s" % (len(defns)+1)
            defns[iTemp] = [dp,"Energy2 %s = %s*%s;" % (dp,ip,ip)]
        pre = "2./3.*(1.-%s*OM%s)" % (dp,struct.value)
        terms.append(("-"+pre,i0,i0))
        terms.append(("%s*OM%s" %(pre,struct.value),ip,ip))
    else :
        # indices of the tensor
        ind1 = struct.lorentz[0]
        ind2 = struct.lorentz[1]
        # the dot products we need
        iTemp = tuple(sorted((ind1,ip),cmp=indSort))
        if(iTemp in defns) :
            d1 = defns[iTemp][0]
        else :
            d1 = "dot%s" % (len(defns)+1)
            unit = "Energy"
            if(ind1.type=="P") : unit="Energy2"
            defns[iTemp] = [d1,"complex<%s> %s = %s*%s;" % (unit,d1,ind1,ip)]
        iTemp = tuple(sorted((ind2,ip),cmp=indSort))
        if(iTemp in defns) :
            d2 = defns[iTemp][0]
        else :
            d2 = "dot%s" % (len(defns)+1)
            unit = "Energy"
            if(ind2.type=="P") : unit="Energy2"
            defns[iTemp] = [d2,"complex<%s> %s = %s*%s;" % (unit,d2,ind2,ip)]
        iTemp = tuple(sorted((ind1,ind2),cmp=indSort))
        if(iTemp in defns) :
            d3 = defns[iTemp][0]
        else :
            d3 = "dot%s" % (len(defns)+1)
            dtemp=0
            if(ind1.type=="P") : dtemp+=1
            if(ind2.type=="P") : dtemp+=1
            unit=computeUnit(dtemp)
            defns[iTemp] = [d3,"complex<%s> %s = %s*%s;" % (unit,d3,ind1,ind2)]
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
        terms.append(("4./3.*OM3*OM%s*%s*%s"%(struct.value,d1,d2),ip,ip))
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
                    val += "%s*%s.%s()*%s.%s()" % (pre,term[1],i1,term[1],i2)
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
        (output[i],parsed[i],dimension[i]) = finishParsing(parsed[i],dimension[i],lorentztag,iloc,defns)
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
            print 'in NOT HANDLED',parsed,output
            for i in range (0,len(parsed)) :
                (output[i],dimension[i]) = finalContractions(output[i],parsed[i],dimension[i],lorentztag,iloc,defns)
                print output[i]
        return (output,dimension,eps)

def convertDirac(output,dimension,eps,iloc,L,parsed,lorentztag,vertex,defns) :
    for i in range(0,len(parsed)):
        # skip empty elements
        if(len(parsed[i])==0 or (len(parsed[i])==1 and parsed[i][0]=="")) : continue
        # parse the string
        (output[i],dimension[i],defns) = convertDiracStructure(parsed[i],output[i],dimension[i],
                                                               defns,iloc,L,lorentztag,vertex)
    print 'after convert loop',output
    print dimension
    return (output,dimension,eps)

# parse the gamma matrices
def convertMatrix(structure,expr,spins,unContracted,Symbols,dtemp,iloc,defns) :
    i1 = structure.spin[0]
    i2 = structure.spin[1]
    if(structure.name=="Identity") :
        expr += "*%s" % I4
        structure=""
    elif(structure.name=="Gamma5") :
        expr += "*%s" % G5
        structure=""
    elif(structure.name=="ProjM") :
        expr += "*%s" % PM
        structure=""
    elif(structure.name=="ProjP") :
        expr += "*%s" % PP
        structure=""
    elif(structure.name=="Gamma") :
        print structure
        # lorentz matrix contracted with the propagator
        if(structure.lorentz[0].type=="D") :
            if(abs(structure.lorentz[0].value) not in unContracted) :
                unContracted.append(abs(structure.lorentz[0].value))
            expr += "*U%s" % abs(structure.lorentz[0].value)
            structure=""
        elif(structure.lorentz[0].type=="E" and
             spins[structure.lorentz[0].value-1]==4) :
            expr += "*R%s" % structure.lorentz[0][1]
            unContracted.append("R%s" % structure.lorentz[0][1])
            structure=""
        elif(structure.lorentz[0].type  == "E" and
             structure.lorentz[0].value == iloc ) :
            expr += "*V%s" % iloc
            unContracted.append("V%s" % iloc)
            structure=""
        else :
            expr += "*"+vslash.substitute({ "v" : structure.lorentz[0]})
            Symbols += vslashS.substitute({ "v" : structure.lorentz[0]})
            if(structure.lorentz[0].type=="P") :
                dtemp[2] += 1
                variable="Energy"
            else :
                variable="double"
            defns["vv%s" % structure.lorentz[0] ] = \
                                                    ["vv%s" % structure.lorentz[0],
                                                     vslashD.substitute({ "var" : variable,
                                                                          "v" : structure.lorentz[0]})]
            structure=""
    else :
        print 'Unknown Gamma matrix structure',structure
        raise SkipThisVertex()
    return (i1,i2,expr,structure,unContracted,Symbols,dtemp)

def convertDiracStructure(parsed,output,dimension,defns,iloc,L,lorentztag,vertex) :
    print 'start of structure',iloc
    # templates
    spinor   = Template("Matrix([[${s}s1],[${s}s2],[${s}s3],[${s}s4]])")
    sbar     = Template("Matrix([[${s}s1,${s}s2,${s}s3,${s}s4]])")
    sline    = Template("${s}s1=Symbol(\"${s}s1\")\n${s}s2=Symbol(\"${s}s2\")\n${s}s3=Symbol(\"${s}s3\")\n${s}s4=Symbol(\"${s}s4\")\n")
    # get the spins of the particles
    spins = vertex.lorentz[0].spins
    # check if we have one or two spin chains
    nchain = (lorentztag.count("F")+lorentztag.count("R"))/2
    expr=[]
    sind=[]
    lind=[]
    start=[]
    end=[]
    startT=[]
    endT=[]
    unContracted=[]
    Symbols=""
    dtemp=[0,0,0]
    for ichain in range(0,nchain) :
        # piece of dimension which is common (0.5 for sbar and spinor)
        dtemp[0]+=1
        # set up the spin indices
        sind.append(0)
        lind.append(0)
        # and the expresion
        expr.append("")
        # now find the next thing in the string
        ii = 0
        index=0
        while True :
            # already handled
            if(parsed[ii]=="") :
                ii+=1
                continue
            # start of the chain
            elif(sind[ichain]==0 and len(parsed[ii].spin)==2 and parsed[ii].spin[0]>0 ) :
                (sind[ichain],index,expr[ichain],parsed[ii],unContracted,Symbols,dtemp) \
                = convertMatrix(parsed[ii],expr[ichain],spins,unContracted,Symbols,dtemp,iloc,defns)
            # next element in the chain
            elif(index!=0 and len(parsed[ii].spin)==2 and parsed[ii].spin[0]==index) :
                (crap,index,expr[ichain],parsed[ii],unContracted,Symbols,dtemp) \
                = convertMatrix(parsed[ii],expr[ichain],spins,unContracted,Symbols,dtemp,iloc,defns)
            # check the index to see if we're at the end
            if(index>0) :
                lind[ichain]=index
                break
            ii+=1
            if(ii>=len(parsed)) :break
        # start and end of the spin chains
        # easy case, both spin 1/2
        if(spins[sind[ichain]-1]==2 and spins[lind[ichain]-1]==2) :
            # start of chain
            # off shell
            if(sind[ichain]==iloc) :
                start.append(vslashM.substitute({ "v" : "P%s" % sind[ichain], "m" : "M%s" % sind[ichain]} ))
                Symbols += vslashMS.substitute({ "v" : "P%s" % sind[ichain], "m" : "M%s" % sind[ichain]} )
                defns["vvP%s" % sind[ichain] ] = ["vvP%s" % sind[ichain] ,
                                                  vslashD.substitute({ "var" : "Energy",
                                                                       "v" :  "P%s" % sind[ichain] })]
                dtemp[1]+=1
            # onshell
            else :
                subs = {'s' : ("sbar%s" % sind[ichain])}
                start.append(sbar .substitute(subs))
                Symbols += sline.substitute(subs)
            # end of chain
            if(lind[ichain]==iloc) :
                end.append(vslashM2.substitute({ "v" : "P%s" % lind[ichain], "m" : "M%s" % lind[ichain]} ))
                Symbols += vslashMS.substitute({ "v" : "P%s" % lind[ichain], "m" : "M%s" % lind[ichain]} )
                defns["vvP%s" % lind[ichain] ] = ["vvP%s" % lind[ichain] ,
                                                  vslashD.substitute({ "var" : "Energy",
                                                                       "v" :  "P%s" % lind[ichain] })]
                dtemp[1] += 1
            else :
                subs = {'s' : ("s%s" % lind[ichain])}
                end.append(spinor.substitute(subs))
                Symbols += sline.substitute(subs)
            startT.append(start[ichain])
            endT  .append(end[ichain])
        # start 1/2 and end 3/2
        elif spins[sind[ichain]-1]==2 and spins[lind[ichain]-1]==4 :
            # spin 1/2 fermion
            if(sind[ichain]==iloc) :
                start.append(vslashM.substitute({ "v" : "P%s" % sind[ichain], "m" : "M%s" % sind[ichain]} ))
                Symbols += vslashMS.substitute({ "v" : "P%s" % sind[ichain], "m" : "M%s" % sind[ichain]} )
                endT.append(vslashM2.substitute({ "v" : "P%s" % sind[ichain], "m" : "M%s" % sind[ichain]} ))
                defns["vvP%s" % sind[ichain] ] = ["vvP%s" % sind[ichain] ,
                                                  vslashD.substitute({ "var" : "Energy",
                                                                       "v" :  "P%s" % sind[ichain] })]
                dtemp[1]+=1
            else :
                  subs = {'s' : ("sbar%s" % sind[ichain])}
                  start.append(sbar .substitute(subs))
                  Symbols += sline.substitute(subs)
                  subs = {'s' : ("s%s" % sind[ichain])}
                  endT.append(spinor.substitute(subs))
                  Symbols += sline.substitute(subs)
            # spin 3/2 fermion
            if(lind[ichain]==iloc) :
                contract=""
                for k in range(1,len(vertex.particles)+1) :
                    if(output.find("(P%s)" %k)>=0) :
                        output = output.replace("(P%s)" %k,"1.")
                        contract="P%s" %k
                end.append(rslash2.substitute({ "v" : "P%s" % lind[ichain], "m" : "M%s" % lind[ichain],
                                                "loc" : lind[ichain] }))
                Symbols += vslashMS.substitute({ "v" : "P%s" % lind[ichain], "m" : "M%s" % lind[ichain]} )
                startT.append(rslash.substitute({ "v" : "P%s" % lind[ichain],
                                                  "m" : "M%s" % lind[ichain], "loc" : lind[ichain] }))
                defns["vvP%s" % lind[ichain] ] = ["vvP%s" % lind[ichain] ,
                                                  vslashD.substitute({ "var" : "Energy",
                                                                       "v" :  "P%s" % lind[ichain] })]
                Symbols += momCom.substitute({"v" : "P%s" %lind[ichain] })
                dtemp[1] += 1
                if(contract!="") :
                    Symbols += momCom.substitute({"v" : contract })
                    RB = vslash.substitute({ "v" : contract})
                    Symbols += vslashS.substitute({ "v" : contract })
                    startT[ichain] = startT[ichain].replace("R%sB"%lind[ichain],RB)
                    end   [ichain] = end   [ichain].replace("R%sB"%lind[ichain],RB)
                    name = "dot%s" % (len(defns)+1)
                    defns[('P%s'%lind[ichain],contract)] = [name,"complex<Energy2> %s = P%s*%s;" % (name,lind[ichain],contract) ]
                    Symbols += "%s = Symbol('%s')\n" % (name,name)
                    startT[ichain] = startT[ichain].replace("P%sB!" % lind[ichain], name).replace("ETA(B!,","ETA(%s," % contract)
                    end   [ichain] = end   [ichain].replace("P%sB!" % lind[ichain], name).replace("ETA(B!,","ETA(%s," % contract)
                    unContracted.append("RC%s"%lind[ichain])
            else :
                # check if we have a contraction issue
                contract=""
                for key,val in defns.iteritems() :
                    if(not isinstance(key,tuple)) : continue
                    if(key[0]== "E%s" %lind[ichain]) :
                        if(output.find("(%s)"%val[0])>=0) :
                            contract=key[1]
                            output = output.replace("(%s)"%val[0],"1.")
                    elif(key[1]=="E%s" %lind[ichain]) :
                        if(output.find("(%s)"%val[0])>=0) :
                            contract=key[0]
                            output = output.replace("(%s)"%val[0],"1.")
                if(contract=="") :
                    end.append(spinor.substitute({'s' : ("Rs%sL" % lind[ichain])}))
                    startT.append(sbar  .substitute({'s' : ("Rsbar%sL" % lind[ichain])}))
                    for LI in imap :
                        Symbols += sline.substitute({'s' : ("Rs%s%s"    % (lind[ichain],LI))})
                        Symbols += sline.substitute({'s' : ("Rsbar%s%s" % (lind[ichain],LI))})
                else :
                    dTemp = Template("${s}ts${si}*${v}t-${s}xs${si}*${v}x-${s}ys${si}*${v}y-${s}zs${si}*${v}z")
                    startT.append("Matrix([[%s,%s,%s,%s]])" % (dTemp.substitute({'s' : ("Rsbar%s" % lind[ichain]), 'v':contract, 'si' : 1}),
                                                               dTemp.substitute({'s' : ("Rsbar%s" % lind[ichain]), 'v':contract, 'si' : 2}),
                                                               dTemp.substitute({'s' : ("Rsbar%s" % lind[ichain]), 'v':contract, 'si' : 3}),
                                                               dTemp.substitute({'s' : ("Rsbar%s" % lind[ichain]), 'v':contract, 'si' : 4})))
                    end   .append("Matrix([[%s],[%s],[%s],[%s]])" % (dTemp.substitute({'s' : ("Rs%s" % lind[ichain]), 'v':contract, 'si' : 1}),
                                                                     dTemp.substitute({'s' : ("Rs%s" % lind[ichain]), 'v':contract, 'si' : 2}),
                                                                     dTemp.substitute({'s' : ("Rs%s" % lind[ichain]), 'v':contract, 'si' : 3}),
                                                                     dTemp.substitute({'s' : ("Rs%s" % lind[ichain]), 'v':contract, 'si' : 4})))
                    Symbols += momCom.substitute({"v" : contract })
                    for LI in ["x","y","z","t"] :
                        Symbols += sline.substitute({'s' : ("Rs%s%s"    % (lind[ichain],LI))})
                        Symbols += sline.substitute({'s' : ("Rsbar%s%s" % (lind[ichain],LI))})
        # start 1/2 and end 3/2
        elif spins[sind[ichain]-1]==4 and spins[lind[ichain]-1]==2 :
            # spin 3/2 fermion
            if(sind[ichain]==iloc) :
                contract=""
                for k in range(1,len(vertex.particles)+1) :
                    if(output.find("(P%s)" %k)>=0) :
                        output = output.replace("(P%s)" %k,"1.")
                        contract="P%s" %k
                start.append(rslash.substitute({ "v" : "P%s" % sind[ichain], "m" : "M%s" % sind[ichain], "loc" : sind[ichain] }))
                Symbols += vslashMS.substitute({ "v" : "P%s" % sind[ichain], "m" : "M%s" % sind[ichain]} )
                endT.append(rslash2.substitute({ "v" : "P%s" % sind[ichain], "m" : "M%s" % sind[ichain], "loc" : sind[ichain] }))
                defns["vvP%s" % sind[ichain] ] = ["vvP%s" % sind[ichain] ,
                                                  vslashD.substitute({ "var" : "Energy",
                                                                       "v" :  "P%s" % sind[ichain] })]
                Symbols += momCom.substitute({"v" : "P%s" %sind[ichain] })
                dtemp[1] += 1
                if(contract!="") :
                    Symbols += momCom.substitute({"v" : contract })
                    RB = vslash.substitute({ "v" : contract})
                    Symbols += vslashS.substitute({ "v" : contract })
                    start[ichain] = start[ichain].replace("R%sB"%sind[ichain],RB)
                    endT [ichain] = endT [ichain].replace("R%sB"%sind[ichain],RB)
                    name = "dot%s" % (len(defns)+1)
                    defns[('P%s'%sind[ichain],contract)] = [name,"complex<Energy2> %s = P%s*%s;" % (name,sind[ichain],contract) ]
                    Symbols += "%s = Symbol('%s')\n" % (name,name)
                    start[ichain] = start[ichain].replace("P%sB!" % sind[ichain], name).replace("ETA(B!,","ETA(%s," % contract)
                    endT [ichain] = endT [ichain].replace("P%sB!" % sind[ichain], name).replace("ETA(B!,","ETA(%s," % contract)
                    unContracted.append("RC%s"%sind[ichain])
            else :
                # check if we have a contraction issue
                contract=""
                for key,val in defns.iteritems() :
                    if(not isinstance(key,tuple)) : continue
                    if(key[0]== "E%s" %sind[ichain]) :
                        if(output.find("(%s)"%val[0])>=0) :
                            contract=key[1]
                            output = output.replace("(%s)"%val[0],"1.")
                    elif(key[1]=="E%s" %sind[ichain]) :
                        if(output.find("(%s)"%val[0])>=0) :
                            contract=key[0]
                            output = output.replace("(%s)"%val[0],"1.")
                if(contract=="") :
                    start = sbar  .substitute({'s' : ("Rsbar%sL" % sind[ichain])})
                    endT  = spinor.substitute({'s' : ("Rs%sL" % sind[ichain])})
                    for LI in imap :
                        Symbols += sline.substitute({'s' : ("Rs%s%s"    % (sind[ichain],LI))})
                        Symbols += sline.substitute({'s' : ("Rsbar%s%s" % (sind[ichain],LI))})
                else :
                    dTemp = Template("${s}ts${si}*${v}t-${s}xs${si}*${v}x-${s}ys${si}*${v}y-${s}zs${si}*${v}z")
                    start.append("Matrix([[%s,%s,%s,%s]])" % (dTemp.substitute({'s' : ("Rsbar%s" % sind[ichain]), 'v':contract, 'si' : 1}),
                                                              dTemp.substitute({'s' : ("Rsbar%s" % sind[ichain]), 'v':contract, 'si' : 2}),
                                                              dTemp.substitute({'s' : ("Rsbar%s" % sind[ichain]), 'v':contract, 'si' : 3}),
                                                              dTemp.substitute({'s' : ("Rsbar%s" % sind[ichain]), 'v':contract, 'si' : 4})))
                    endT .append("Matrix([[%s],[%s],[%s],[%s]])" % (dTemp.substitute({'s' : ("Rs%s" % sind[ichain]), 'v':contract, 'si' : 1}),
                                                                    dTemp.substitute({'s' : ("Rs%s" % sind[ichain]), 'v':contract, 'si' : 2}),
                                                                    dTemp.substitute({'s' : ("Rs%s" % sind[ichain]), 'v':contract, 'si' : 3}),
                                                                    dTemp.substitute({'s' : ("Rs%s" % sind[ichain]), 'v':contract, 'si' : 4})))
                    Symbols += momCom.substitute({"v" : contract })
                    for LI in ["x","y","z","t"] :
                        Symbols += sline.substitute({'s' : ("Rs%s%s"    % (sind[ichain],LI))})
                        Symbols += sline.substitute({'s' : ("Rsbar%s%s" % (sind[ichain],LI))})
            if(lind[ichain]==iloc) :
                end   .append(vslashM2.substitute({ "v" : "P%s" % lind[ichain], "m" : "M%s" % lind[ichain]} ))
                startT.append(vslashM.substitute({ "v" : "P%s" % lind[ichain], "m" : "M%s" % lind[ichain]} ))
                Symbols += vslashMS.substitute({ "v" : "P%s" % lind[ichain], "m" : "M%s" % lind[ichain]} )
                defns["vvP%s" % lind[ichain] ] = ["vvP%s" % lind[ichain] ,
                                                  vslashD.substitute({ "var" : "Energy",
                                                                       "v" :  "P%s" % lind[ichain] })]
                dtemp[1] += 1
            else :
                subs = {'s' : ("s%s" % lind[ichain])}
                end      .append(spinor.substitute(subs))
                Symbols += sline.substitute(subs)
                subs = {'s' : ("sbar%s" % lind[ichain])}
                startT   .append(sbar .substitute(subs))
                Symbols += sline.substitute(subs)
        else :
            print 'At most one R-S spinor allowed in a vertex'
            raise SkipThisVertex()
    # remove leading *
    for i in range(0,nchain) : expr[i]=expr[i][1:]
    # check we've dealt with everything
    parsed = [x for x in parsed if x != ""]
    if(len(parsed)!=0) :
        for i in range(0,len(parsed)) :
            if(parsed[i].name=="Metric") :
                found = False
                for ll in parsed[i].lorentz :
                    if(ll.type=="E" and ll.value==iloc) :
                        found = True
                    else :
                        lo=ll
                if(found) :
                    parsed[i]=""
                    for i in range(0,nchain) :
                        expr[i]+="*%s(V%s)" % (lo,iloc)
                    unContracted.append("V%s" % iloc)
                    if(lo.type=="P") : dimension[2]+=1
        parsed = [x for x in parsed if x != ""]
        if(len(parsed)!=0) :
            print expr
            print "Can't parse ",parsed,iloc
            quit()
            raise SkipThisVertex()
    sVal ={}
    # deal with the simplest case first
    if len(unContracted) == 0 :
        print 'testing in simple case'
        temp={}
        exec("import sympy\nfrom sympy import Symbol,Matrix\n"+Symbols+"result="+
             ( "%s*%s*%s" %(start[0],expr[0],end[0]))) in temp
        tempT={}
        exec("import sympy\nfrom sympy import Symbol,Matrix,Transpose\n"+Symbols+"result="+
             ( "%s*%s*Transpose(%s)*%s*%s" %(startT[0],CC,expr[0],CD,endT[0]))) in tempT
        if(iloc==0 or (iloc!=sind[ichain] and iloc!=lind[ichain])) :
            sVal = {'s' : temp ["result"][0,0],'sT' : tempT["result"][0,0]}
        else :
            for jj in range(1,5) :
                sVal["s%s"  % jj] = temp ["result"][jj-1]
                sVal["sT%s" % jj] = tempT["result"][jj-1]
    print 'before un',sVal
    # now deal with the uncontracted cases
    contracted=[]
    # sort out contracted and uncontracted indices
    for j in range(0,len(unContracted)) :
        # summed index
        if(isinstance(unContracted[j],int)) :
            contracted.append(unContracted[j])
            unContracted.remove(unContracted[j])
        # RS index
        elif unContracted[j][0]=="R" :
            if(unContracted[j][1]=="C") :
                unContracted[j] = unContracted[j].replace("C","")
            else :
                contracted.append(unContracted[j])
                if(int(unContracted[j][1])!=iloc):
                    unContracted.remove(unContracted[j])
        # uncontracted vector index
        elif unContracted[j][0]=="V" :
            continue
        else :
            print 'uknown type of uncontracted index',unContracted[j]
            raise SkipThisVertex()
    # iterate over the uncontracted indices
    unI=[0]*len(unContracted)
    defns["I"] = ["I","static Complex I(0.,1.);"]
    print 'before un loop',unContracted,contracted
    while True :
        if(len(contracted)==0 and len(unContracted)==0) : break
        coI=[0]*len(contracted)
        # loop over the contracted indices
        res = []
        for i in range(0,nchain) :
            res.append([])
            res.append([])
        while True :
            sign = 1
            sTemp  = copy.copy(start )
            sTTemp = copy.copy(startT)
            eTemp  = copy.copy(expr  )
            fTemp  = copy.copy(end   )
            fTTemp = copy.copy(endT  )
            # make the necessary replacements for uncontracted indices
            for ichain in range(0,nchain) :
                for j in range(0,len(unContracted)) :
                    if(unContracted[j][0]=="R") :
                        ii = int(unContracted[j][1])
                        sTemp[ichain]   = sTemp[ichain].replace(unContracted[j]+"A",dirac[unI[j]])
                        sTTemp[ichain]  = sTTemp[ichain].replace(unContracted[j]+"A",dirac[unI[j]])
                        fTemp[ichain]   = fTemp[ichain].replace(unContracted[j]+"A",dirac[unI[j]])
                        fTTemp[ichain]  = fTTemp[ichain].replace(unContracted[j]+"A",dirac[unI[j]])
                        sTemp[ichain]  =  sTemp[ichain].replace("A!",imap[unI[j]])
                        sTTemp[ichain] = sTTemp[ichain].replace("A!",imap[unI[j]])
                        fTemp[ichain]  =  fTemp[ichain].replace("A!",imap[unI[j]])
                        fTTemp[ichain] = fTTemp[ichain].replace("A!",imap[unI[j]])
                        # handle metric tensors
                        eta=re.search("ETA\(.*,.\)",sTemp[ichain])
                        if(eta) :
                            eta = eta.group(0)
                            if(eta.find("!")<0) : 
                                temp = eta.split("(")[1].split(",")
                                replace = temp[0]+temp[1][0]
                                sTemp[ichain]  =  sTemp[ichain].replace(eta,replace)
                                sTTemp[ichain] = sTTemp[ichain].replace(eta,replace)
                                fTemp[ichain]  =  fTemp[ichain].replace(eta,replace)
                                fTTemp[ichain] = fTTemp[ichain].replace(eta,replace)
                    elif(unContracted[j][0]=="V") :
                        print unContracted[j][0]
                        print eTemp[ichain]
                        loc = eTemp[ichain].find("(%s)"%unContracted[j])
                        if(loc>=0) :
                            vec=eTemp[ichain][loc-2:loc]
                            Symbols += momCom.substitute({"v":vec})
                            print "!!!!!!",vec
                            eTemp[ichain]   = eTemp[ichain].replace("(%s)"%unContracted[j],imap[unI[j]]) 
                        else :
                           eTemp[ichain]   = eTemp[ichain].replace(unContracted[j],dirac[unI[j]])
                        print eTemp[ichain]
                    else :
                        print 'uknown type of uncontracted index',unContracted[j]
                        raise SkipThisVertex()
            # make the necessary replacements for contracted indices
            for ichain in range(0,nchain) :
                for j in range(0,len(contracted)) :
                    if(isinstance(contracted[j],int)) :
                        eTemp[ichain]=eTemp[ichain].replace("U%s"%contracted[j],dirac[coI[j]])
                        if(ichain==0 and coI[j]>0) : sign *= -1
                    elif(contracted[j][0]=="R") :
                        ii = int(contracted[j][1])
                        # replace metric
                        for k in range(0,4) :
                            test = "ETA(B!,%s)" % (imap[k])
                            esign="0"
                            if(coI[j]==k) : esign="1"
                            sTemp[ichain]  =  sTemp[ichain].replace(test,esign)
                            sTTemp[ichain] = sTTemp[ichain].replace(test,esign)
                            fTemp[ichain]  =  fTemp[ichain].replace(test,esign)
                            fTTemp[ichain] = fTTemp[ichain].replace(test,esign)
                        # replace dirac matrices
                        eTemp[ichain]   = eTemp[ichain].replace(contracted[j],dirac[coI[j]])
                        sTemp[ichain]   = sTemp[ichain].replace(contracted[j]+"B",dirac[coI[j]])
                        sTTemp[ichain]  = sTTemp[ichain].replace(contracted[j]+"B",dirac[coI[j]])
                        fTemp[ichain]   = fTemp[ichain].replace(contracted[j]+"B",dirac[coI[j]])
                        fTTemp[ichain]  = fTTemp[ichain].replace(contracted[j]+"B",dirac[coI[j]])
                        # replacements for start
                        sTemp[ichain]  =  sTemp[ichain].replace("Rsbar%sL" % ii,"Rsbar%s%s" % (ii,imap[coI[j]]))
                        sTTemp[ichain] = sTTemp[ichain].replace("Rsbar%sL" % ii,"Rsbar%s%s" % (ii,imap[coI[j]]))
                        sTemp[ichain]  =  sTemp[ichain].replace("B!",imap[coI[j]])
                        sTTemp[ichain] = sTTemp[ichain].replace("B!",imap[coI[j]])
                        # replacements for end
                        fTemp[ichain]  = fTemp[ichain].replace("Rs%sL" % ii,"Rs%s%s" % (ii,imap[coI[j]]))
                        fTTemp[ichain] = fTTemp[ichain].replace("Rs%sL" % ii,"Rs%s%s" % (ii,imap[coI[j]]))
                        fTemp[ichain]  =  fTemp[ichain].replace("B!",imap[coI[j]])
                        fTTemp[ichain] = fTTemp[ichain].replace("B!",imap[coI[j]])
                        if(ichain==0 and coI[j]>0) : sign *= -1
                    else :
                        print 'need to implment',contracted[j]
                        quit()
            # now evaluate the result
            rTemp =[[],[]]
            for ichain in range(0,nchain) :
                temp={}
                exec("import sympy\nfrom sympy import Symbol,Matrix\n"+Symbols+"result="+
                     ( "(%s)*(%s)*(%s)" %(sTemp[ichain],eTemp[ichain],fTemp[ichain]))) in temp
                rTemp[0].append(temp["result"])
                temp={}
                exec("import sympy\nfrom sympy import Symbol,Matrix,Transpose\n"+Symbols+"result="+
                    ( "(%s)*(%s)*(Transpose(%s))*(%s)*(%s)" %(sTTemp[ichain],CC,eTemp[ichain],CD,fTTemp[ichain]))) in temp
                rTemp[1].append(temp["result"])
            # and add it to the output
            # 1 spin chain
            if(nchain==1) :
                for ii in range(0,2) :
                    if(rTemp[ii][0].shape[0]==1) :
                        if(rTemp[ii][0].shape[1]==1) :
                            if(len(res[ii])==0) :
                                res[ii].append(sign*rTemp [ii][0][0,0])
                            else :
                                res[ii][0]  += sign*rTemp [ii][0][0,0]
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
            #### END OF THE CONTRACTED LOOP #####
            # increment the indices being summed over
            ii = len(coI)-1
            while ii >=0 :
                if(coI[ii]<3) :
                    coI[ii]+=1
                    break
                else :
                    coI[ii]=0
                    ii-=1
            if(coI.count(0)==len(coI)) : break
        ###### END OF THE UNCONTRACTED LOOP ######
        # no uncontracted indices
        if(len(unI)==0) :
            if(len(res[0])==1) :
                if(len(res)==2) :
                    sVal["s" ] = res[0]
                    sVal["sT"] = res[1]
                else :
                    sVal["s"   ] = res[0][0]
                    sVal["sT2" ] = res[1][0]
                    sVal["sT1" ] = res[2][0]
                    sVal["sT12"] = res[3][0]
            elif(len(res)==4) :
                for k in range(0,4) :
                    sVal[ "s%s"  % (k+1) ] = res[0][k]
                    sVal[ "sT%s" % (k+1) ] = res[1][k]
            break
        # uncontracted indices
        else :
            istring = "" 
            for val in unI:
                istring +=imap[val]
            if(len(istring)!=1) :
                print 'index problem',istring
                print unContracted
                print unI
                quit()
            sVal[istring]     = res[0]
            sVal[istring+"T"] = res[1]
            ii = len(unI)-1
            while ii >=0 :
                if(unI[ii]<3) :
                    unI[ii]+=1
                    break
                else :
                    unI[ii]=0
                    ii-=1
            if(unI.count(0)==len(unI)) : break
    # handle the vector case
    if( "t" in sVal ) :
        # deal with pure vectors
        if(len(sVal)==8 and "t" in sVal and len(sVal["t"])==1) :
            for key in sVal:
                sVal[key] = sVal[key][0]

            
            # unit = computeUnit(dtemp)
            # sVal   = { "s"  : "LorentzVector<complex<%s> >(%s,%s,%s,%s)" % (unit,sVal["x" ][0],sVal["y" ][0],sVal["z" ][0],sVal["t" ][0]),
            #            "sT" : "LorentzVector<complex<%s> >(%s,%s,%s,%s)" % (unit,sVal["xT"][0],sVal["yT"][0],sVal["zT"][0],sVal["tT"][0]) }
            # print sVal
        else :
            print 'val problrm A'
            quit()
    #             # print 'testing vector problem',len(sVal),len(sVal["t"])
    #             # quit()
    #             #                     # if(expr=="") :
    #             #                     #     expr ={}
    #             #                     #     exprT={}
    #             #                     #     defns["I"] = ["I","static Complex I(0.,1.);"]
    #             #                     #     unit = computeUnit(dtemp)
    #             #                     #     if(pre=="") :
    #             #                     #         exprT["s"] = "LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
    #             #                     #                      (unit,vals[i][2]["x"],vals[i][2]["y"],vals[i][2]["z"],vals[i][2]["t"])
    #             #                     #     else :
    #             #                     #         expr ["s"] = "%s*LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
    #             #                     #                     (pre,unit,vals[i][1]["x"],vals[i][1]["y"],vals[i][1]["z"],vals[i][1]["t"])
    #             #                     #         exprT["s"] = "%s*LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
    #             #                     #                      (pre,unit,vals[i][2]["x"],vals[i][2]["y"],vals[i][2]["z"],vals[i][2]["t"])
    #             #                     # else :
    #             #                     #     unit = computeUnit(dtemp)
    #             #                     #     if(pre=="") :
    #             #                     #         expr ["s"] += "LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
    #             #                     #                     (unit,vals[i][1]["x"],vals[i][1]["y"],vals[i][1]["z"],vals[i][1]["t"])
    #             #                     #         exprT["s"] += "LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
    #             #                     #                      (unit,vals[i][2]["x"],vals[i][2]["y"],vals[i][2]["z"],vals[i][2]["t"])
    #             #                     #     else :
    #             #                     #         expr ["s"] += "%s*LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
    #             #                     #                     (pre,unit,vals[i][1]["x"],vals[i][1]["y"],vals[i][1]["z"],vals[i][1]["t"])
    #             #                     #         exprT["s"] += "%s*LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
    #             #                     #                      (pre,unit,vals[i][2]["x"],vals[i][2]["y"],vals[i][2]["z"],vals[i][2]["t"])
    #             #                     # print 'CCCCC',len(vals[i])
    old = output
    if(nchain==1) :
        print "SETTING UP THE OUTPUT ",old
        print "B",sVal
        output = [old,sVal,(lind[0],sind[0])]
    else :
        output = [old,sVal,(lind[0],sind[0],lind[1],sind[1])]
    dimension = list(map(lambda x, y: x + y, dtemp, dimension))
    # remove any dot products involving RS fermions
    if("R" in lorentztag ) :
        for key in defns.keys() :
            if(not isinstance(key,tuple)) : continue
            if(key[0]== "E%s" %sind or key[1]=="E%s" %sind or
               key[0]== "E%s" %lind[ichain] or key[1]=="E%s" %lind[ichain]) :
                del defns[key]
    return (output,dimension,defns)
            
def convertLorentz(Lstruct,lorentztag,order,vertex,iloc,defns,evalVertex) :
    eps = False
    # split the structure into individual terms
    structures=Lstruct.structure.split()
    parsed=[]
    for struct in structures :
        print struct
        parsed.append(parse_structure(struct))
    # convert lorentz contractions to dot products
    print "before loop ",  Lstruct.structure
    print parsed
    for l in range(0,len(parsed)) :
        parsed[l] = contract(parsed[l])
    print 'after loop',parsed
    # now in a position to generate the code
    vals=generateVertex(iloc,Lstruct,parsed,lorentztag,vertex,defns)
    print iloc,Lstruct,lorentztag
    print vals
    evalVertex.append((vals[0],vals[1]))
    if(vals[2]) : eps=True
    print "END OF CONVERT LORENTZ",Lstruct,iloc,len(evalVertex)
    print evalVertex
    return eps

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

def constructSignature(vertex,order,iloc,decls,momenta,waves,fermionReplace,fIndex) :
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
                if(i==1) :
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
                    fermionReplace.append("s%s"%(fIndex[i-1]))
                    fermionReplace.append("sbar%s"%(fIndex[i-1]))
                    nf+=1
                else :
                    decls.append("SpinorBarWaveFunction & sbarW%s" % (fIndex[i-1]))
                    momenta.append([False,"Lorentz5Momentum P%s =-sbarW%s.momentum();" % (fIndex[i-1],fIndex[i-1])])
                    waves.append("LorentzSpinorBar<double> sbar%s = sbarW%s.wave();" % (fIndex[i-1],fIndex[i-1]))
                    fermionReplace.append("sbar%s"%(fIndex[i-1]))
                    fermionReplace.append("s%s"%(fIndex[i-1]))
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
                    fermionReplace.append("Rs%s"%(i))
                    fermionReplace.append("Rsbar%s"%(i))
                    nf+=1
                else :
                    decls.append("RSSpinorBarWaveFunction & RsbarW%s" % (i))
                    momenta.append([False,"Lorentz5Momentum P%s =-RsbarW%s.momentum();" % (i,i)])
                    waves.append("LorentzRSSpinorBar<double> Rsbar%s = RsbarW%s.wave();" % (i,i))
                    fermionReplace.append("Rs%s"%(i))
                    fermionReplace.append("Rsbar%s"%(i))
                    nf+=1
            elif(spin==5) :
                decls.append("TensorWaveFunction & tW%s" % (i))
                momenta.append([False,"Lorentz5Momentum P%s =-tW%s.momentum();" % (i,i)])
                waves.append("LorentzTensor<double> T%s = tW%s.wave();" % (i,i))
            else :
                print 'unknown spin',spin
                quit()
            poff += "-P%s" % (i)
    poff = ("Lorentz5Momentum P%s = " % iloc ) + poff
    sig=""
    if(iloc==0) :
        sig="%s evaluate(Energy2, const %s)" % (offType,", const ".join(decls))
    else :
        sig="%s evaluate(Energy2, int iopt, tcPDPtr out, const %s, complex<Energy> mass=-GeV, complex<Energy> width=-GeV)" % (offType,", const ".join(decls))
        momenta.append([True,poff+";"])
        for i in range(0,len(momenta)) : momenta[i][0]=True
    return offType,nf,poff,sig

def combineResult(res,nf,ispin,fermionReplace,vertex) :
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
            print vals
            print 'spin problem',ispin
            quit()
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
            print "DIMENSION PROBLEM",i,dimCheck,dim,vertex,vals
            print vals
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
            print expr
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
        elif(len(vals[i][1])%4==0 and "t" in vals[i][1]) :
            for i1 in imap :
                if(pre=="") :
                    expr[0][i1] += "(%s)" % vals[i][1][i1]
                else :
                    expr[0][i1] += "%s*(%s)" % (pre,vals[i][1][i1])
                if(len(expr)==2) :
                    if(pre=="") :
                        expr[1][i1] +="(%s)" % vals[i][1][i1+"T"]
                    else :
                        expr[1][i1] +="%s*(%s)" % (pre,vals[i][1][i1+"T"])
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
        else :
            print vals[i]
            print 'problem with type'
            quit()
        #             # unknown
        #             else :
        #                 print 'unrecongnised case',vals[i]
        #                 quit()
        #                 #     elif(len(vals[i][1])==4 and "t" in vals[i][1]) :
        #                 #         # two cases, either already a 'simple' vector
        #                 #         # or its an RS spinor
        #                 #         # RS first
        #                 #         if(len(vals[i][1]["t"])==4) :
        #                 #             if(expr=="") :
        #                 #                 expr ={}
        #                 #                 exprT={}
        #                 #                 for jj in range(0,4) :
        #                 #                     for k in range(1,5) :
        #                 #                         if(pre=="") :
        #                 #                             expr ["%ss%s" % (imap[jj],k)] = "(%s)" % vals[i][1][imap[jj]][k-1]
        #                 #                             exprT["%ss%s" % (imap[jj],k)] = "(%s)" % vals[i][2][imap[jj]][k-1]
        #                 #                         else :
        #                 #                             expr ["%ss%s" % (imap[jj],k)] = "%s*(%s)" % (pre,vals[i][1][imap[jj]][k-1])
        #                 #                             exprT["%ss%s" % (imap[jj],k)] = "%s*(%s)" % (pre,vals[i][2][imap[jj]][k-1])
        #                 #             else :
        #                 #                 for jj in range(0,4) :
        #                 #                     for k in range(1,5) :
        #                 #                         if(pre=="") :
        #                 #                             expr ["%ss%s" % (imap[jj],k)] += "(%s)" % vals[i][1][imap[jj]][k-1]
        #                 #                             exprT["%ss%s" % (imap[jj],k)] += "(%s)" % vals[i][2][imap[jj]][k-1]
        #                 #                         else :
        #                 #                             expr ["%ss%s" % (imap[jj],k)] += "%s*(%s)" % (pre,vals[i][1][imap[jj]][k-1])
        #                 #                             exprT["%ss%s" % (imap[jj],k)] += "%s*(%s)" % (pre,vals[i][2][imap[jj]][k-1])
        #                 #         # simple vector
        #                 #         else :
        #                 #             print 'vector case'
        #                 #             print vals[i]
        #                 #             print expr
        #                 #             quit()

        #                 #     else :
        #                 #         print "AAAAAAA"
        #                 #         print vertex,len(vals[i][1])
        #                 #         quit()


    # tidy up fermion defns
    for rep in fermionReplace :
        for i in range(1,5) :
            kmax=1
            if(rep[0]=="R") : kmax=4
            for k in range(0,kmax) :
                if(rep[0]=="R") :
                    oldVal = "%s%ss%s" % (rep,imap[k],i)
                    newVal = "%s.%ss%s()" % (rep,imap[k],i)
                else :
                    oldVal = "%ss%s" % (rep,i)
                    newVal = "%s.s%s()" % (rep,i)
                for ii in range (0,len(expr)) :
                    for (key,val) in expr[ii].iteritems() :
                        expr[ii][key]  = val.replace(oldVal,newVal)
    # no of particles in the vertex
    vDim = len(vertex.lorentz[0].spins)
    # tidy up momentum components
    for i in range(1,vDim+1) :
        for k in range(0,4) :
            oldVal = "P%s%s*" % (i,imap[k])
            newVal = "P%s.%s()*" % (i,imap[k])
            for ii in range(0,len(expr)) :
                for (key,val) in expr[ii].iteritems() :
                    expr[ii][key]  = val.replace(oldVal,newVal)
    
    unit = computeUnit2(dimCheck,vDim)
    if(unit!="") :
        for ii in range(0,len(expr)) :
            for (key,val) in expr[ii].iteritems() :
                expr[ii][key]  = "(%s)*(%s)" % (val,unit)
    return expr

def combineComponents(result,offType,RS) :
    # simplest case, just a value
    if(len(result[0])==1 and "res" in result[0]) :
        for i in range(0,len(result)) :
            result[i] = result[i]["res"]
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
        result[0]  = Template("${type}<double>(${outs1},\n${outs2},\n${outs3},\n${outs4})").substitute(subs[0])
        subs[1]["type"] = sbtype
        result[1]  = Template("${type}<double>(${outs1},\n${outs2},\n${outs3},\n${outs4})").substitute(subs[1])
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
    else :
        
        print subs[0]
        print result
        print 'type not implemented'
        quit()
    
    # # final defns for more complicated types
                
    #     if("ts1" in result) :
    #         stype  = "LorentzRSSpinor"
    #         sbtype = "LorentzRSSpinorBar"
    #         if(offType.find("Bar")>0) : (stype,sbtype)=(sbtype,stype)
    #         subs["type"] = stype      
    #         result  = Template("${type}<double>(${outxs1},\n${outxs2},\n${outxs3},\n${outxs4},\n${outys1},\n${outys2},\n${outys3},\n${outys4},\n${outzs1},\n${outzs2},\n${outzs3},\n${outzs4},\n${outts1},\n${outts2},\n${outts3},\n${outts4})").substitute(subs)
    #         subsT["type"] = sbtype
    #         resultT = Template("${type}<double>(${outxs1},\n${outxs2},\n${outxs3},\n${outxs4},\n${outys1},\n${outys2},\n${outys3},\n${outys4},\n${outzs1},\n${outzs2},\n${outzs3},\n${outzs4},\n${outts1},\n${outts2},\n${outts3},\n${outts4})").substitute(subsT)
    #     else :
    #         print result
    #         print 'type problem'
    #         quit()



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
    fermionReplace=[]
    offType,nf,poff,sig = constructSignature(vertex,order,iloc,decls,momenta,waves,fermionReplace,fIndex)
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
        print 'calling combine',j,vertex.lorentz[j]
        expr = combineResult(vertexEval[j],nf,ispin,fermionReplace,vertex)
        # get the coupling for this bit
        val, sym = py2cpp(values[j])
        localCouplings.append("Complex local_C%s = %s;\n" % (j,val))
        symbols |=sym
        # put them together
        if(len(result)==0) :
            for ii in range(0,len(expr)) :
                result.append({})
                for (key,val) in expr[ii].iteritems() :
                    result[ii][key] = " (local_C%s)*Complex(%s) " % (j,val)
        else :
            for ii in range(0,len(expr)) :
                for (key,val) in expr[ii].iteritems(): 
                    result[ii][key] += " + (local_C%s)*Complex(%s) " % (j,val)
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
                result = vTemplateT.format(iloc=fIndex[1],id=vertex.particles[fIndex[1]-1].pdg_code,
                                           cf=py2cpp(cf)[0],res=result[0],resT=result[1])
            else :
                result = vTemplate4.format(iloc1=fIndex[1],iloc2=fIndex[3],
                                           id1=vertex.particles[fIndex[1]-1].pdg_code,
                                           id2=vertex.particles[fIndex[3]-1].pdg_code,
                                           cf=py2cpp(cf)[0],res1=result[0],res2=result[1],res3=result[2],res4=result[3])
        else :
            result = "return (%s)*(%s);\n" % (result[0],py2cpp(cf[0])[0])
        if(RS) :
            result[1] = "return (%s)*(%s);\n" % (result[1],py2cpp(cf[0])[0])
    # off-shell particle
    else :
        # off-shell scalar
        if(vertex.lorentz[0].spins[iloc-1] == 1 ) :
            if(RS) :
                resultT = scaTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],res=result[1])
            if(FM and not RS) :
                result = scaTTemplate.format(iloc=iloc,res=result[0],resT=result[1],isp=fIndex[1],
                                             id=vertex.particles[fIndex[1]-1].pdg_code,cf=py2cpp(cf[0])[0])
            else :
                result = scaTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],res=result[0])
        # off-shell fermion
        elif(vertex.lorentz[0].spins[iloc-1] == 2 ) :
            if(RS) :
                if(offType.find("Bar")>0) :
                    offTypeT=offType.replace("Bar","")
                else :
                    offTypeT=offType.replace("Spinor","SpinorBar")
                resultT = sTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offTypeT.replace("WaveFunction",""),
                                           res=result[1].replace( "M%s" % iloc, "mass" ),offTypeB=offTypeT)
            if(FM and not RS) :
                result = sTTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offType.replace("WaveFunction",""),
                                           res=result[0].replace( "M%s" % iloc, "mass" ),resT=result[1].replace( "M%s" % iloc, "mass" ),
                                           offTypeB=offType,id=vertex.particles[iloc-1].pdg_code)
            else :
                result = sTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offType.replace("WaveFunction",""),
                                          res=result[0].replace( "M%s" % iloc, "mass" ),offTypeB=offType)
        # off-shell vector
        elif(vertex.lorentz[0].spins[iloc-1] == 3 ) :
            if(FM and not RS) :
                result = vecTTemplate.format(iloc=iloc,res=result[0],resT=result[1],isp=fIndex[1],
                                             id=vertex.particles[fIndex[1]-1].pdg_code,
                                             cf=py2cpp(cf[0])[0])
            else :
                result = vecTemplate.format(iloc=iloc,res=result[0],cf=py2cpp(cf[0])[0])
        elif(vertex.lorentz[0].spins[iloc-1] == 4 ) :
            if(offType.find("Bar")>0) :
                offTypeT=offType.replace("Bar","")
            else :
                offTypeT=offType.replace("Spinor","SpinorBar")
            resultT = RSTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offTypeT.replace("WaveFunction",""),
                                        res=result[1].replace( "M%s" % iloc, "mass" ),offTypeB=offTypeT)
            result = RSTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offType.replace("WaveFunction",""),
                                       res=result[0].replace( "M%s" % iloc, "mass" ),offTypeB=offType)
        # tensors
        elif(vertex.lorentz[0].spins[iloc-1]) :
            if(RS) :
                print "RS spinors and tensors not handled"
                quit()
            if(FM) :
                result = tenTTemplate.format(iloc=iloc,res=result[0],resT=result[1],isp=fIndex[1],
                                             id=vertex.particles[fIndex[1]-1].pdg_code,cf=py2cpp(cf[0])[0])
            else :
                result = tenTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],res=result[0])
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
                                       result=result,swap=sorder)

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
                              .replace("Energy2","q2").replace("int","").replace("complex<Energy>","").replace("=-GeV","") \
                              .replace("const  &","").replace("tcPDPtr","").replace("  "," ")
            if(iloc!=1) :
                call = call.replace(waves[0],waves[iloc-1])
            pdgid = vertex.particles[i].pdg_code
            code += "   %sif(out->id()==%s) return %s;\n" % (el,pdgid,call)
            iloc+=1
    code+="   else assert(false);\n"
    return (header,evaluateMultiple.format(decl=ccdefn,code=code))
