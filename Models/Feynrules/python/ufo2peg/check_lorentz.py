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
                        tval2 = evaluate(all_couplings[icolor][0],model,parmsubs)
                        tval3 = -evaluate(all_couplings[icolor][1],model,parmsubs)
                    elif(col[0][0]==1 and col[0][1]==3 and col[1][0] ==2 and col[1][1] == 4) : 
                        if(all_couplings[icolor][1] or not all_couplings[icolor][0] or
                           not all_couplings[icolor][2]) :
                            raise SkipThisVertex()
                        if(not value) :
                            value = all_couplings[icolor][2]
                            tval  = evaluate(value,model,parmsubs)
                        tval2 = evaluate(all_couplings[icolor][0],model,parmsubs)
                        tval3 = -evaluate(all_couplings[icolor][2],model,parmsubs)
                    elif(col[0][0]==1 and col[0][1]==4 and col[1][0] ==2 and col[1][1] == 3) : 
                        if(all_couplings[icolor][0] or not all_couplings[icolor][1] or
                           not all_couplings[icolor][2]) :
                            raise SkipThisVertex()
                        if(not value) :
                            value = all_couplings[icolor][2]
                            tval  = evaluate(value,model,parmsubs)
                        tval2 = evaluate(all_couplings[icolor][1],model,parmsubs)
                        tval3 = -evaluate(all_couplings[icolor][2],model,parmsubs)
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

class LorentzStructure:
    """A simple example class to store a Lorentz structures"""
    name=""
    value=0.
    lorentz=[]
    spin=[]
    def __repr__(self):
        output = self.name
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
        ind=extractIndices(struct)
        # different types of object
        # object with only spin indices
        if(struct.find("Identity")==0 or
           struct.find("Proj")==0 or
           struct.find("Gamma5")==0) :
            output.append(LorentzStructure())
            output[-1].spin=ind
            output[-1].name=struct.split("(")[0]
            output[-1].value=1.
            if(len(struct.replace("%s(%s,%s)" % (output[-1].name,ind[0],ind[1]),""))!=0) :
                print "problem A"
                quit()
        # objects with 2 lorentz indices
        elif(struct.find("Metric")==0 or struct.find("P(")==0) :
            output.append(LorentzStructure())
            output[-1].lorentz=ind
            output[-1].name=struct.split("(")[0]
            output[-1].value=1.
            if(len(struct.replace("%s(%s,%s)" % (output[-1].name,ind[0],ind[1]),""))!=0) :
                print "problem B"
                quit()
        # 1 lorentz and 1 spin index
        elif(struct.find("Gamma")==0) :
            output.append(LorentzStructure())
            output[-1].lorentz=[ind[0]]
            output[-1].spin=[ind[1],ind[2]]
            output[-1].name=struct.split("(")[0]
            output[-1].value=1.
            if(len(struct.replace("%s(%s,%s,%s)" % (output[-1].name,ind[0],ind[1],ind[2]),""))!=0) :
                print "problem C",struct
                quit()
        # objects with 4 lorentz indices
        elif(struct.find("Epsilon")==0) :
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
    return sorted(output,cmp=LorentzCompare)

def contractLorentz(L,parsed,lorentztag,order) :
    for l in range(0,len(parsed)) :
        for j in range(0,len(parsed[l])) :
            # replace indices with polarization vectors
            ll = len(parsed[l][j].lorentz)
            if(parsed[l][j].name=="P") :
                ll=1
                found=False
                for k in range(0,len(parsed[l])) :
                    if(j==k or parsed[l][k]=="" ) : continue
                    imax = len(parsed[l][k].lorentz)
                    if(parsed[l][k].name=="P") : imax=1
                    for i in range(0,imax) :
                        if(parsed[l][k].lorentz[i]==parsed[l][j].lorentz[0]) :
                            parsed[l][k].lorentz[i] = "P%s" % parsed[l][j].lorentz[1]
                            if(parsed[l][k].name=="P") :
                                parsed[l][k].lorentz[1] = "P%s" % parsed[l][k].lorentz[1]
                                parsed[l][k].name="Metric"
                            found=True
                            break
                if(found) :
                    parsed[l][j]=''
                    ll=0
            for i in range(0,ll) :
                if(parsed[l][j]!="" and parsed[l][j].lorentz[i]>=0
                   and isinstance(parsed[l][j].lorentz[i],(int,long)) and
                   L.spins[parsed[l][j].lorentz[i]-1] in [3,4] ) :
                    parsed[l][j].lorentz[i]=  "E%s" % parsed[l][j].lorentz[i]
                    if(parsed[l][j].name=="P") :
                        parsed[l][j].lorentz[1] = "P%s" % parsed[l][j].lorentz[1]
                        parsed[l][j].name="Metric"
        parsed[l] = [x for x in parsed[l] if x != ""]

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
                expr = "UnitRemoval::E%s" % (-dimension)
    if(output=="") : return expr
    elif(expr=="") : return output
    else           : return "%s*%s" %(output,expr)

def generateVertex(iloc,L,parsed,lorentztag,vertex,defns) :
    eps=False
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
    spins = vertex.lorentz[0].spins
    # order the indices of a dot product
    def indSort(a,b) :
        if(a[0]==b[0]) :
            i1=int(a[1])
            i2=int(b[1])
            if(i1>i2) :
                return 1
            elif(i1<i2) :
                return -1
            else :
                return 0
        else :
            if(a[0]=="E") :
                return 1
            else :
                return -1
    # parse the lorentz structures
    output    = [1.]*len(parsed)
    dimension=[]
    for i in range(0,len(parsed)) :
        dimension.append([0,0,0])
    for i in range (0,len(parsed)) :
        # replace signs
        if(len(parsed[i])==1 and parsed[i][0].name=="sign") :
            if(parsed[i][0].value>0) :
                output[i]="+"
            else :
                output[i]="-"
            parsed[i][0]=''
            continue
        # replace integers
        for j in range(0,len(parsed[i])) :
            if(parsed[i][j]!="" and parsed[i][j].name=="int") :
                output[i] *= parsed[i][j].value
                parsed[i][j]=""
                continue
        output[i] = "(%s)" % output[i]
        for j in range(0,len(parsed[i])) :
            if(parsed[i][j]!="" and parsed[i][j].name=="Metric") :
                (ind1,ind2) = sorted((parsed[i][j].lorentz[0],parsed[i][j].lorentz[1]),cmp=indSort)
                # this product already dealt with ?
                if((ind1,ind2) in defns) :
                    output[i] += "*(%s)" % defns[(ind1,ind2)][0]
                    parsed[i][j]=""
                    if(ind1[0]=="P") : dimension[i][2] +=1
                    if(ind2[0]=="P") : dimension[i][2] +=1
                    continue
                # handle the product
                name = "dot%s" % (len(defns)+1)
                parsed[i][j]=""
                if(ind1[0]=="P") :
                    # dot product of two momenta
                    if(ind2[0]=="P") :
                        dimension[i][2] +=2
                        defns[(ind1,ind2)] = [name,"Energy2 %s = %s*%s;" % (name,ind1,ind2)]
                        output[i] += "*(%s)" % name
                    elif(ind2[0]=="E") :
                        dimension[i][2] +=1
                        if(int(ind2[1])==iloc) :
                            output[i] += "*(%s)" % ind1
                        else :
                            defns[(ind1,ind2)] = [name,"complex<Energy> %s = %s*%s;" % (name,ind1,ind2)]
                            output[i] += "*(%s)" % name
                elif(ind1[0]=="E") :
                    if(ind2[0]!="E") :
                        print "EE problem"
                        quit()
                    if(int(ind1[1])==iloc) :
                        output[i] += "*(%s)" % ind2
                    elif(int(ind2[1])==iloc) :
                        output[i] += "*(%s)" % ind1
                    else :
                        defns[(ind1,ind2)] = [name,"complex<double> %s = %s*%s;" % (name,ind1,ind2)]
                        output[i] += "*(%s)" % name
            elif(parsed[i][j]!="" and parsed[i][j].name=="Epsilon") :
                if(not eps) : eps = True
                offLoc = -1
                indices=[]
                dTemp=0
                for ix in range(0,len(parsed[i][j].lorentz)) :
                    if(isinstance(parsed[i][j].lorentz[ix],int)) :
                        offLoc = ix
                        break
                    elif(parsed[i][j].lorentz[ix][0]=="E" and int(parsed[i][j].lorentz[ix][1])==iloc ) :
                        offLoc = ix
                        break
                for ix in range(0,len(parsed[i][j].lorentz)) :
                    if(isinstance(parsed[i][j].lorentz[ix],basestring) and
                       parsed[i][j].lorentz[ix][0]=="P") : dTemp+=1
                    if((offLoc<0 and ix != 0) or
                       (offLoc>=0 and offLoc!=ix) ) :
                        indices.append(parsed[i][j].lorentz[ix])
                dimension[i][2] += dTemp
                if(offLoc<0) :
                    iTemp = (parsed[i][j].lorentz[0],parsed[i][j].lorentz[1],
                             parsed[i][j].lorentz[2],parsed[i][j].lorentz[3])
                    if(iTemp in defns) :
                        output[i] += "*(%s)" % defns[iTemp][0]
                        parsed[i][j]=""
                    else :
                        name = "dot%s" % (len(defns)+1)
                        unit = computeUnit(dTemp)
                        defns[iTemp] = [name,"complex<%s> %s =-%s*epsilon(%s,%s,%s);" % (unit,name,parsed[i][j].lorentz[0],
                                                                                        indices[0],indices[1],indices[2]) ]
                        output[i] += "*(%s)" % name
                else :
                    iTemp = (indices[0],indices[1],indices[2])
                    sign = ""
                    if(offLoc%2!=0) : sign="-"
                    print iTemp
                    if(iTemp in defns) :
                        output[i] += "*(%s%s)" % (sign,defns[iTemp][0])
                        parsed[i][j]=""
                    else :
                        name = "vec%s" % (len(defns)+1)
                        output[i] += "*(%s%s)" % (sign,name)
                        unit = computeUnit(dTemp)
                        defns[iTemp] = [name,"LorentzVector<complex<%s> > %s =-epsilon(%s,%s,%s);" % (unit,name,
                                                                                                     indices[0],indices[1],indices[2]) ]
    # remove any (now) empty elements
    for i in range (0,len(parsed)) :
        parsed[i] = [x for x in parsed[i] if x != ""]
    # now for gamma matrix strings
    if(lorentztag[0] in ["F","R"] ) :
        if(lorentztag[0]=="R") :
            print 'in spin bit',parsed
            print lorentztag
        result=""
        for i in range(0,len(parsed)):
            unContracted=[]
            print "PARSED A ",i
            if(len(parsed[i])==0) :
                continue
            if(len(parsed[i])==1 and parsed[i][0]=="") :
                continue
            # first and last lorentz indices
            sind=parsed[i][0].spin[0]
            lind=0
            # first piece of the expression we need to evaluate
            dtemp=[0,0,0]
            expr=""
            Symbols=""
            print "PARSED A ",i
            # parse the gamma matrices
            for j in range(0,len(parsed[i])) :
                print parsed[i][j]
                if(parsed[i][j].name=="Identity") :
                    expr += "*%s" % I4
                    lind = parsed[i][j].spin[1]
                    parsed[i][j]=""
                elif(parsed[i][j].name=="Gamma5") :
                    expr += "*%s" % G5
                    lind = parsed[i][j].spin[1]
                    parsed[i][j]=""
                elif(parsed[i][j].name=="ProjM") :
                    expr += "*%s" % PM
                    lind = parsed[i][j].spin[1]
                    parsed[i][j]=""
                elif(parsed[i][j].name=="ProjP") :
                    expr += "*%s" % PP
                    lind = parsed[i][j].spin[1]
                    parsed[i][j]=""
                elif(parsed[i][j].name=="Gamma") :
                    lind = parsed[i][j].spin[1]
                    # lorentz matrix contracted with the propagator
                    if(parsed[i][j].lorentz[0][0]=="E" and
                         spins[int(parsed[i][j].lorentz[0][1])-1]==4) :
                        expr += "*R%s" % parsed[i][j].lorentz[0][1]
                        unContracted.append("R%s" % parsed[i][j].lorentz[0][1])
                    elif(parsed[i][j].lorentz[0] == ("E%s" % iloc ) ) :
                        expr += "*V%s" % iloc
                        unContracted.append("V%s" % iloc)
                    else :
                        expr += "*"+vslash.substitute({ "v" : parsed[i][j].lorentz[0]})
                        Symbols += vslashS.substitute({ "v" : parsed[i][j].lorentz[0]})
                        if(parsed[i][j].lorentz[0][0]=="P") :
                            dtemp[2] += 1
                            variable="Energy"
                        else :
                            variable="double"
                        defns["vv%s" % parsed[i][j].lorentz[0] ] = \
                                                                   ["vv%s" % parsed[i][j].lorentz[0],
                                                                    vslashD.substitute({ "var" : variable,
                                                                                         "v" : parsed[i][j].lorentz[0]})]
                    parsed[i][j]=""
                else :
                    print 'FFFFFF'
                    print parsed[i][j]
                    quit()
            print "PARSED C ",i
            # piece of dimension which is common (0.5 for sbar and spinor)
            dtemp[0]+=1
            # start and end of the spin chains
            # easy case, both spin 1/2
            spinor   = Template("Matrix([[${s}s1],[${s}s2],[${s}s3],[${s}s4]])")
            sbar     = Template("Matrix([[${s}s1,${s}s2,${s}s3,${s}s4]])")
            sline    = Template("${s}s1=Symbol(\"${s}s1\")\n${s}s2=Symbol(\"${s}s2\")\n${s}s3=Symbol(\"${s}s3\")\n${s}s4=Symbol(\"${s}s4\")\n")
            print "PARSED D ",i
            if(spins[sind-1]==2 and spins[lind-1]==2) :
                print 'spin block A'
                if(sind==iloc) :
                    start    = vslashM.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
                    Symbols += vslashMS.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
                    defns["vvP%s" % sind ] = ["vvP%s" % sind ,
                                              vslashD.substitute({ "var" : "Energy",
                                                                   "v" :  "P%s" % sind })]
                    dtemp[1]+=1
                else :
                    subs = {'s' : ("sbar%s" % sind)}
                    start    = sbar .substitute(subs)
                    Symbols += sline.substitute(subs)
                if(lind==iloc) :
                    end = vslashM2.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
                    Symbols += vslashMS.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
                    defns["vvP%s" % lind ] = ["vvP%s" % lind ,
                                              vslashD.substitute({ "var" : "Energy",
                                                                   "v" :  "P%s" % lind })]
                    dtemp[1] += 1
                else :
                    subs = {'s' : ("s%s" % lind)}
                    end      = spinor.substitute(subs)
                    Symbols += sline.substitute(subs)
                startT = start
                endT   = end
            elif spins[sind-1]==2 and spins[lind-1]==4 :
                print 'spin block B'
                # spin 1/2 fermion
                if(sind==iloc) :
                    start    = vslashM.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
                    Symbols += vslashMS.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
                    endT = vslashM2.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
                    defns["vvP%s" % sind ] = ["vvP%s" % sind ,
                                              vslashD.substitute({ "var" : "Energy",
                                                                   "v" :  "P%s" % sind })]
                    dtemp[1]+=1
                else :
                    subs = {'s' : ("sbar%s" % sind)}
                    start    = sbar .substitute(subs)
                    Symbols += sline.substitute(subs)
                    subs = {'s' : ("s%s" % sind)}
                    endT     = spinor.substitute(subs)
                    Symbols += sline.substitute(subs)
                # spin 3/2 fermion
                if(lind==iloc) :
                    print 'testing 3/2 lind',lind
                    end      = rslash2.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind, "loc" : lind })
                    Symbols += vslashMS.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
                    startT   = rslash.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind, "loc" : lind })
                    defns["vvP%s" % lind ] = ["vvP%s" % lind ,
                                              vslashD.substitute({ "var" : "Energy",
                                                                   "v" :  "P%s" % lind })]
                    Symbols += momCom.substitute({"v" : "P%s" %lind })
                    dtemp[1] += 1
                else :
                    end    = spinor.substitute({'s' : ("Rs%sL" % lind)})
                    startT = sbar  .substitute({'s' : ("Rsbar%sL" % lind)})
                    for LI in ["x","y","z","t"] :
                        Symbols += sline.substitute({'s' : ("Rs%s%s"    % (lind,LI))})
                        Symbols += sline.substitute({'s' : ("Rsbar%s%s" % (lind,LI))})
            elif spins[sind-1]==4 and spins[lind-1]==2 :
                # spin 3/2 fermion
                if(sind==iloc) :
                    contract=""
                    for k in range(1,len(vertex.particles)+1) :
                        if(output[i].find("(P%s)" %k)>=0) :
                            output[i] = output[i].replace("(P%s)" %k,"1.")
                            contract="P%s" %k
                    start = rslash.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind, "loc" : sind })
                    Symbols += vslashMS.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
                    endT = rslash2.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind, "loc" : sind })
                    defns["vvP%s" % sind ] = ["vvP%s" % sind ,
                                              vslashD.substitute({ "var" : "Energy",
                                                                   "v" :  "P%s" % sind })]
                    Symbols += momCom.substitute({"v" : "P%s" %sind })
                    Symbols += momCom.substitute({"v" : contract })
                    dtemp[1] += 1
                    if(contract!="") :
                        RB = vslash.substitute({ "v" : contract})
                        Symbols += vslashS.substitute({ "v" : contract })
                        start = start.replace("R%sB"%sind,RB)
                        endT  = endT .replace("R%sB"%sind,RB)
                        name = "dot%s" % (len(defns)+1)
                        defns[('P%s'%sind,contract)] = [name,"complex<Energy2> %s = P%s*%s;" % (name,sind,contract) ]
                        Symbols += "%s = Symbol('%s')\n" % (name,name)
                        start = start.replace("P%sB!" % sind, name).replace("ETA(B!,","ETA(%s," % contract)
                        endT  = endT .replace("P%sB!" % sind, name).replace("ETA(B!,","ETA(%s," % contract)
                        unContracted.append("RC%s"%sind)
                else :
                    # check if we have a contraction issue
                    contract=""
                    for key,val in defns.iteritems() :
                        if(not isinstance(key,tuple)) : continue
                        if(key[0]== "E%s" %sind) :
                            if(output[i].find("(%s)"%val[0])>=0) :
                                contract=key[1]
                                output[i] = output[i].replace("(%s)"%val[0],"1.")
                        elif(key[1]=="E%s" %sind) :
                            if(output[i].find("(%s)"%val[0])>=0) :
                                contract=key[0]
                                output[i] = output[i].replace("(%s)"%val[0],"1.")
                    if(contract=="") :
                        print 'testing 3/2 sind on',sind
                        start = sbar  .substitute({'s' : ("Rsbar%sL" % sind)})
                        endT  = spinor.substitute({'s' : ("Rs%sL" % sind)})
                        for LI in ["x","y","z","t"] :
                            Symbols += sline.substitute({'s' : ("Rs%s%s"    % (sind,LI))})
                            Symbols += sline.substitute({'s' : ("Rsbar%s%s" % (sind,LI))})
                    else :
                        print "testing in special bit !!!",contract
                        dTemp = Template("${s}ts${si}*${v}t-${s}xs${si}*${v}x-${s}ys${si}*${v}y-${s}zs${si}*${v}z")
                        start = "Matrix([[%s,%s,%s,%s]])" % (dTemp.substitute({'s' : ("Rsbar%s" % sind), 'v':contract, 'si' : 1}),
                                                             dTemp.substitute({'s' : ("Rsbar%s" % sind), 'v':contract, 'si' : 2}),
                                                             dTemp.substitute({'s' : ("Rsbar%s" % sind), 'v':contract, 'si' : 3}),
                                                             dTemp.substitute({'s' : ("Rsbar%s" % sind), 'v':contract, 'si' : 4}))
                        endT = "Matrix([[%s],[%s],[%s],[%s]])" % (dTemp.substitute({'s' : ("Rs%s" % sind), 'v':contract, 'si' : 1}),
                                                                  dTemp.substitute({'s' : ("Rs%s" % sind), 'v':contract, 'si' : 2}),
                                                                  dTemp.substitute({'s' : ("Rs%s" % sind), 'v':contract, 'si' : 3}),
                                                                  dTemp.substitute({'s' : ("Rs%s" % sind), 'v':contract, 'si' : 4}))
                        Symbols += momCom.substitute({"v" : contract })
                        for LI in ["x","y","z","t"] :
                            Symbols += sline.substitute({'s' : ("Rs%s%s"    % (sind,LI))})
                            Symbols += sline.substitute({'s' : ("Rsbar%s%s" % (sind,LI))})
                if(lind==iloc) :
                    print 'testing 1/2 sind off',sind
                    end    = vslashM2.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
                    startT = vslashM.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
                    Symbols += vslashMS.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
                    defns["vvP%s" % lind ] = ["vvP%s" % lind ,
                                              vslashD.substitute({ "var" : "Energy",
                                                                   "v" :  "P%s" % lind })]
                    dtemp[1] += 1
                else :
                    subs = {'s' : ("s%s" % lind)}
                    end      = spinor.substitute(subs)
                    Symbols += sline.substitute(subs)
                    subs = {'s' : ("sbar%s" % lind)}
                    startT   = sbar .substitute(subs)
                    Symbols += sline.substitute(subs)
            else :
                print "RR"
                quit()
            print "PARSED E ",i
            # check we've dealt with everything
            parsed[i] = [x for x in parsed[i] if x != ""]
            expr=expr[1:]
            if(len(parsed[i])!=0) :
                print "ERROR"
                quit()
            # deal with the simplest case first
            print "PARSED F ",i,iloc
            if len(unContracted) == 0 :
                temp={}
                exec("import sympy\nfrom sympy import Symbol,Matrix\n"+Symbols+"result="+
                     ( "%s*%s*%s" %(start,expr,end))) in temp
                tempT={}
                exec("import sympy\nfrom sympy import Symbol,Matrix,Transpose\n"+Symbols+"result="+
                     ( "%s*%s*Transpose(%s)*%s*%s" %(startT,CC,expr,CD,endT))) in tempT
                if(iloc==0 or (iloc!=sind and iloc!=lind)) :
                    old = output[i]
                    output[i] = [("%s*(%s)" % (old,temp ["result"][0,0])),
                                 ("%s*(%s)" % (old,tempT["result"][0,0])),
                                 (sind,lind)]
                    output[i] = [old,{'s' : temp ["result"][0,0]},{'s' : tempT["result"][0,0]},(sind,lind)]
                else :
                    old = output[i]
                    sVal ={}
                    sTVal={}
                    for jj in range(1,5) :
                        sVal ["s%s" % jj] = temp ["result"][jj-1]
                        sTVal["s%s" % jj] = tempT["result"][jj-1]
                    output[i] = [old,sVal,sTVal,(sind,lind)]
                dimension[i] = list(map(lambda x, y: x + y, dtemp, dimension[i]))
                continue
            print "PARSED G ",i,iloc
            # now deal with the uncontracted cases
            contracted=[]
            for j in range(0,len(unContracted)) :
                print 'in unco loop',unContracted[j]
                # sort out contracted and uncontracted indices
                if unContracted[j][0]=="R" :
                    if(unContracted[j][1]=="C") :
                        unContracted[j] = unContracted[j].replace("C","")
                    else :
                        contracted.append(unContracted[j])
                        if(int(unContracted[j][1])!=iloc):
                            unContracted.remove(unContracted[j])
                elif unContracted[j][0]=="V" :
                    continue
                else :
                    print 'unA',unContracted[j]
                    quit()
            
            unI=[0]*len(unContracted)
            print "PARSED H ",i
            print unContracted
            print   contracted
            sVal ={}
            sTVal={}
            defns["I"] = ["I","static Complex I(0.,1.);"]
            while True :
                print "PARSED I ",i
                coI=[0]*len(contracted)
                # loop over the contracted indices
                res =[]
                resT=[]
                while True :
                    sign = 1
                    sTemp  = copy.copy(start )
                    sTTemp = copy.copy(startT)
                    eTemp  = copy.copy(expr  )
                    fTemp  = copy.copy(end   )
                    fTTemp = copy.copy(endT  )
                    # make the necessary replacements for uncontracted indices
                    for j in range(0,len(unContracted)) :
                        if(unContracted[j][0]=="R") :
                            ii = int(unContracted[j][1])
                            sTemp   = sTemp.replace(unContracted[j]+"A",dirac[unI[j]])
                            sTTemp  = sTTemp.replace(unContracted[j]+"A",dirac[unI[j]])
                            fTemp   = fTemp.replace(unContracted[j]+"A",dirac[unI[j]])
                            fTTemp  = fTTemp.replace(unContracted[j]+"A",dirac[unI[j]])
                            sTemp  =  sTemp.replace("A!",imap[unI[j]])
                            sTTemp = sTTemp.replace("A!",imap[unI[j]])
                            fTemp  =  fTemp.replace("A!",imap[unI[j]])
                            fTTemp = fTTemp.replace("A!",imap[unI[j]])
                            # handle metric tensors
                            eta=re.search("ETA\(.*,.\)",sTemp)
                            if(eta) :
                                eta = eta.group(0)
                                if(eta.find("!")<0) : 
                                    temp = eta.split("(")[1].split(",")
                                    replace = temp[0]+temp[1][0]
                                    sTemp  =  sTemp.replace(eta,replace)
                                    sTTemp = sTTemp.replace(eta,replace)
                                    fTemp  =  fTemp.replace(eta,replace)
                                    fTTemp = fTTemp.replace(eta,replace)
                        elif(unContracted[j][0]=="V") :
                            print 'did vector'
                            eTemp   = eTemp.replace(unContracted[j],dirac[unI[j]])
                        else :
                            print "uncon",unContracted[j]
                            quit()
                    # make the necessary replacements for contracted indices
                    for j in range(0,len(contracted)) :
                        if(contracted[j][0]=="R") :
                            ii = int(contracted[j][1])
                            # replace metric
                            for k in range(0,4) :
                                test = "ETA(B!,%s)" % (imap[k])
                                esign="0"
                                if(coI[j]==k) : esign="1"
                                sTemp  =  sTemp.replace(test,esign)
                                sTTemp = sTTemp.replace(test,esign)
                                fTemp  =  fTemp.replace(test,esign)
                                fTTemp = fTTemp.replace(test,esign)
                            # replace dirac matrices
                            eTemp   = eTemp.replace(contracted[j],dirac[coI[j]])
                            sTemp   = sTemp.replace(contracted[j]+"B",dirac[coI[j]])
                            sTTemp  = sTTemp.replace(contracted[j]+"B",dirac[coI[j]])
                            fTemp   = fTemp.replace(contracted[j]+"B",dirac[coI[j]])
                            fTTemp  = fTTemp.replace(contracted[j]+"B",dirac[coI[j]])
                            # replacements for start
                            sTemp  =  sTemp.replace("Rsbar%sL" % ii,"Rsbar%s%s" % (ii,imap[coI[j]]))
                            sTTemp = sTTemp.replace("Rsbar%sL" % ii,"Rsbar%s%s" % (ii,imap[coI[j]]))
                            sTemp  =  sTemp.replace("B!",imap[coI[j]])
                            sTTemp = sTTemp.replace("B!",imap[coI[j]])
                            # replacements for end
                            fTemp  = fTemp.replace("Rs%sL" % ii,"Rs%s%s" % (ii,imap[coI[j]]))
                            fTTemp = fTTemp.replace("Rs%sL" % ii,"Rs%s%s" % (ii,imap[coI[j]]))
                            fTemp  =  fTemp.replace("B!",imap[coI[j]])
                            fTTemp = fTTemp.replace("B!",imap[coI[j]])
                            if(coI[j]>0) : sign *= -1
                        else :
                            print 'need to implment'
                            quit()
                    # evaluate the results
                    temp={}
                    exec("import sympy\nfrom sympy import Symbol,Matrix\n"+Symbols+"result="+
                         ( "(%s)*(%s)*(%s)" %(sTemp,eTemp,fTemp))) in temp
                    tempT={}
                    exec("import sympy\nfrom sympy import Symbol,Matrix,Transpose\n"+Symbols+"result="+
                         ( "(%s)*(%s)*(Transpose(%s))*(%s)*(%s)" %(sTTemp,CC,eTemp,CD,fTTemp))) in tempT
                    # add the result
                    print temp["result"].shape
                    if(temp["result"].shape[0]==1) :
                        if(temp["result"].shape[1]==1) :
                            print 'doing scalar'
                            if(len(res)==0) :
                                if(sign==1) :
                                    res .append(temp ["result"][0,0])
                                    resT.append(tempT["result"][0,0])
                                else :
                                    res .append(-temp ["result"][0,0])
                                    resT.append(-tempT["result"][0,0])
                            else :
                                if(sign==1) :
                                    res [0] += temp ["result"][0,0]
                                    resT[0] += tempT["result"][0,0]
                                else :
                                    res [0] -= temp ["result"][0,0]
                                    resT[0] -= tempT["result"][0,0]
                        elif(temp["result"].shape[1]==4) :
                            print 'doing column'
                            if(len(res)==0) :
                                for j in range(0,4) :
                                    if(sign==1) :
                                        res .append(temp ["result"][0,j])
                                        resT.append(tempT["result"][j,0])
                                    else :
                                        res .append(-temp ["result"][0,j])
                                        resT.append(-tempT["result"][j,0])
                            else :
                                for j in range(0,4) :
                                    if(sign==1) :
                                        res [j] += temp ["result"][0,j]
                                        resT[j] += tempT["result"][j,0]
                                    else :
                                        res [j] -= temp ["result"][0,j]
                                        resT[j] -= tempT["result"][j,0]
                        else :
                            print "SIZE PROBLEM",sign,temp["result"].shape
                            quit()
                    elif(temp["result"].shape[0]==4 and
                         temp["result"].shape[1]==1 ) :
                        print 'doing row'
                        if(len(res)==0) :
                            for j in range(0,4) :
                                if(sign==1) :
                                    res .append(temp ["result"][j,0])
                                    resT.append(tempT["result"][0,j])
                                else :
                                    res .append(-temp ["result"][j,0])
                                    resT.append(-tempT["result"][0,j])
                        else :
                            for j in range(0,4) :
                                if(sign==1) :
                                    res [j] += temp ["result"][j,0]
                                    resT[j] += tempT["result"][0,j]
                                else :
                                    res [j] -= temp ["result"][j,0]
                                    resT[j] -= tempT["result"][0,j]
                    else :
                        print "SIZE PROBLEM",sign,temp["result"].shape
                        quit()
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
                # no uncontracted indices
                if(len(unI)==0) :
                    if(len(res)==1) :
                        sVal ["s"] = res[0]
                        sTVal["s"] = resT[0]
                    elif(len(res)==4) :
                        for k in range(0,4) :
                            sVal [ "s%s" % (k+1) ] = res [k]
                            sTVal[ "s%s" % (k+1) ] = resT[k]
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
                    sVal [istring] = res
                    sTVal[istring] = resT
                    print 'in loop indice',istring,res
                    ii = len(unI)-1
                    while ii >=0 :
                        if(unI[ii]<3) :
                            unI[ii]+=1
                            break
                        else :
                            unI[ii]=0
                            ii-=1
                    if(unI.count(0)==len(unI)) : break
            # deal with pure vectors
            if(len(sVal)==4 and "t" in sVal and len(sVal["t"])==1) :
                unit = computeUnit(dtemp)
                sVal   = { "s" : "LorentzVector<complex<%s> >(%s,%s,%s,%s)" % (unit,sVal ["x"][0],sVal ["y"][0],sVal ["z"][0],sVal ["t"][0]) }
                sTVal  = { "s" : "LorentzVector<complex<%s> >(%s,%s,%s,%s)" % (unit,sTVal["x"][0],sTVal["y"][0],sTVal["z"][0],sTVal["t"][0]) }
                # print newV
                
                # print 'testing vector problem',len(sVal),len(sVal["t"])
                # quit()
                #                     # if(expr=="") :
                #                     #     expr ={}
                #                     #     exprT={}
                #                     #     defns["I"] = ["I","static Complex I(0.,1.);"]
                #                     #     unit = computeUnit(dtemp)
                #                     #     if(pre=="") :
                #                     #         exprT["s"] = "LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
                #                     #                      (unit,vals[i][2]["x"],vals[i][2]["y"],vals[i][2]["z"],vals[i][2]["t"])
                #                     #     else :
                #                     #         expr ["s"] = "%s*LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
                #                     #                     (pre,unit,vals[i][1]["x"],vals[i][1]["y"],vals[i][1]["z"],vals[i][1]["t"])
                #                     #         exprT["s"] = "%s*LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
                #                     #                      (pre,unit,vals[i][2]["x"],vals[i][2]["y"],vals[i][2]["z"],vals[i][2]["t"])
                #                     # else :
                #                     #     unit = computeUnit(dtemp)
                #                     #     if(pre=="") :
                #                     #         expr ["s"] += "LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
                #                     #                     (unit,vals[i][1]["x"],vals[i][1]["y"],vals[i][1]["z"],vals[i][1]["t"])
                #                     #         exprT["s"] += "LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
                #                     #                      (unit,vals[i][2]["x"],vals[i][2]["y"],vals[i][2]["z"],vals[i][2]["t"])
                #                     #     else :
                #                     #         expr ["s"] += "%s*LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
                #                     #                     (pre,unit,vals[i][1]["x"],vals[i][1]["y"],vals[i][1]["z"],vals[i][1]["t"])
                #                     #         exprT["s"] += "%s*LorentzVector<complex<%s> >(%s,%s,%s,%s)" % \
                #                     #                      (pre,unit,vals[i][2]["x"],vals[i][2]["y"],vals[i][2]["z"],vals[i][2]["t"])
                #                     # print 'CCCCC',len(vals[i])
            # number of cases to deal with
            print len(output),i
            old = output[i]
            output[i] = [old,sVal,sTVal,(sind,lind)]
            dimension[i] = list(map(lambda x, y: x + y, dtemp, dimension[i]))
    # remove any dot products involving RS fermions
    if(lorentztag[0] in ["F","R"] ) :
        for key in defns.keys() :
            if(not isinstance(key,tuple)) : continue
            if(key[0]== "E%s" %sind or key[1]=="E%s" %sind or
               key[0]== "E%s" %lind or key[1]=="E%s" %lind) :
                del defns[key]
    # return the answer
    return (output,dimension,eps)
            
def convertLorentz(Lstruct,lorentztag,order,vertex,iloc,defns,evalVertex) :
    eps = False
    # split the structure into individual terms
    structures=Lstruct.structure.split()
    parsed=[]
    for struct in structures :
        print struct
        parsed.append(parse_structure(struct))
    # convert lorentz contractions to dot products
    contractLorentz(Lstruct,parsed,lorentztag,order)
    # now in a position to generate the code
    vals=generateVertex(iloc,Lstruct,parsed,lorentztag,vertex,defns)
    evalVertex.append((vals[0],vals[1]))
    if(vals[2]) : eps=True
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

def swapOrder(vertex,iloc,momenta) :
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
                if(vertex.particles[sloc[j]-1].name!=vertex.particles[sloc[j]-1].antiname) :
                    code *= -1
                output += "    if(id%s!=%s) {\n" % (sloc[j],code)
                output += "        swap(id%s,id%s);\n" % (sloc[j],sloc[k])
                output += "        swap(%s%s,%s%s);\n" % (waves[i],sloc[j],waves[i],sloc[k])
                if(momenta[sloc[j]-1][0] or momenta[sloc[k]-1][0]) :
                    momenta[sloc[j]-1][0] = True
                    momenta[sloc[k]-1][0] = True
                    output += "        swap(P%s,P%s);\n" % (sloc[j],sloc[k])
                output += "    };\n"
    return output
    
def generateEvaluateFunction(model,vertex,iloc,values,defns,vertexEval,cf,order) :
    RS = "R" in vertex.lorentz[0].name
    print 'START OF FUNCTION WRITE',RS,vertex,iloc,order
    print vertexEval
    # first construct the signature of the function
    iferm=0
    decls=[]
    offType="Complex"
    momenta=[]
    waves=[]
    poff=""
    fermionReplace=[]
    # for i in range(0,len(vertex.lorentz[0].spins)) :
    #     spin = vertex.lorentz[0].spins[i]
    nf=0
    for i in order : 
        spin = vertex.lorentz[0].spins[i-1]
        print 'in loop',i,spin
        if(i==iloc) :
            if(spin==1) :
                offType="ScalarWaveFunction"
            elif(spin==2) :
                if(i==1) :
                    offType="SpinorBarWaveFunction"
                else :
                    offType="SpinorWaveFunction"
            elif(spin==4) :
                if(i==1) :
                    offType="RSSpinorBarWaveFunction"
                else :
                    offType="RSSpinorWaveFunction"
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
                if(i==1) :
                    decls.append("SpinorWaveFunction & sW%s" % (i))
                    momenta.append([False,"Lorentz5Momentum P%s =-sW%s.momentum();" % (i,i)])
                    waves.append("LorentzSpinor<double> s%s = sW%s.wave();" % (i,i))
                    fermionReplace.append("s%s"%(i))
                    fermionReplace.append("sbar%s"%(i))
                    nf+=1
                else :
                    decls.append("SpinorBarWaveFunction & sbarW%s" % (i))
                    momenta.append([False,"Lorentz5Momentum P%s =-sbarW%s.momentum();" % (i,i)])
                    waves.append("LorentzSpinorBar<double> sbar%s = sbarW%s.wave();" % (i,i))
                    fermionReplace.append("sbar%s"%(i))
                    fermionReplace.append("s%s"%(i))
                    nf+=1
            elif(spin==3) :
                decls.append("VectorWaveFunction & vW%s" % (i))
                momenta.append([False,"Lorentz5Momentum P%s =-vW%s.momentum();" % (i,i)])
                waves.append("LorentzPolarizationVector E%s = vW%s.wave();" % (i,i))
            elif(spin==4) :
                if(i==1) :
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
                waves.append("LorentzTensor t%s = tW%s.wave()" % (i,i))
            else :
                print 'unknown spin',spin
                quit()
            poff += "-P%s" % (i)
    poff = ("Lorentz5Momentum P%s = " % iloc ) + poff
    # make sure unbar first
    if(nf==2 and decls[0].find("Bar")>0) :
        decls[0],decls[1]=decls[1],decls[0]

    
    sig=""
    if(iloc==0) :
        sig="%s evaluate(Energy2, const %s)" % (offType,", const ".join(decls))
    else :
        sig="%s evaluate(Energy2, int iopt, tcPDPtr out, const %s, complex<Energy> mass=-GeV, complex<Energy> width=-GeV)" % (offType,", const ".join(decls))
        momenta.append([True,poff+";"])
        for i in range(0,len(momenta)) : momenta[i][0]=True
    # cat the definitions
    defString=""
    for (key,value) in defns.iteritems() :
        defString+="    %s\n" %value[1]
    oval=""
    symbols=set()
    localCouplings=[]
    result =""
    resultT=""
    hasTranspose = False
    fermions=0
    for j in range(0,len(vertexEval)) :
        (vals,dim) = vertexEval[j]
        expr =""
        exprT=""
        dimCheck=dim[0]
        for i in range(0,len(vals)) :
            if(vals[i]=="+" or vals[i]=="-") :
                if(isinstance(expr,basestring)) :
                    expr  +=vals[i]
                    exprT +=vals[i]
                else :
                    for(key,val) in expr.iteritems() :
                        expr[key] = expr[key]+vals[i]
                    for(key,val) in exprT.iteritems() :
                        exprT[key] = exprT[key]+vals[i]
            else :
                # check the dimensions
                if(dimCheck[0]!=dim[i][0] or dimCheck[1]!=dim[i][1] or
                   dimCheck[2]!=dim[i][2]) :
                    print defns
                    print vertex.lorentz[j]
                    print vertex.lorentz[j].structure
                    print "DIMENSION PROBLEM",i,j,dimCheck,dim[i],vertex,vals[i]
                    print vals
                    quit()
                # simplest case 
                if(isinstance(vals[i], basestring)) :
                    expr  += "(%s)" % vals[i]
                    exprT += "(%s)" % vals[i]
                else :
                    # old default behaviour
                    if(len(vals[i])==3 and isinstance(vals[i][0],basestring) ) :
                        hasTranspose = True
                        expr  += "(%s)" % vals[i][0]
                        exprT += "(%s)" % vals[i][1]
                        fermions=vals[i][2]
                    else :
                        hasTranspose = True
                        pre = vals[i][0]
                        fermions=vals[i][3]
                        if(pre=="(1.0)") : pre=""
                        if(isinstance(vals[i][1],dict)) :
                            if(len(vals[i][1])==1 and "s" in vals[i][1]) :
                                if(pre=="") :
                                    expr  += "(%s)" % vals[i][1]["s"]
                                    exprT += "(%s)" % vals[i][2]["s"]
                                else :
                                    print "ADD TEST",expr,"\n",pre,vals[i][1]
                                    expr  += "%s*(%s)" % (pre,vals[i][1]["s"])
                                    exprT += "%s*(%s)" % (pre,vals[i][2]["s"])
                            elif(len(vals[i][1])==4 and "t" in vals[i][1]) :
                                # two cases, either already a 'simple' vector
                                # or its an RS spinor
                                # RS first
                                if(len(vals[i][1]["t"])==4) :
                                    if(expr=="") :
                                        expr ={}
                                        exprT={}
                                        for jj in range(0,4) :
                                            for k in range(1,5) :
                                                if(pre=="") :
                                                    expr ["%ss%s" % (imap[jj],k)] = "(%s)" % vals[i][1][imap[jj]][k-1]
                                                    exprT["%ss%s" % (imap[jj],k)] = "(%s)" % vals[i][2][imap[jj]][k-1]
                                                else :
                                                    expr ["%ss%s" % (imap[jj],k)] = "%s*(%s)" % (pre,vals[i][1][imap[jj]][k-1])
                                                    exprT["%ss%s" % (imap[jj],k)] = "%s*(%s)" % (pre,vals[i][2][imap[jj]][k-1])
                                    else :
                                        for jj in range(0,4) :
                                            for k in range(1,5) :
                                                if(pre=="") :
                                                    expr ["%ss%s" % (imap[jj],k)] += "(%s)" % vals[i][1][imap[jj]][k-1]
                                                    exprT["%ss%s" % (imap[jj],k)] += "(%s)" % vals[i][2][imap[jj]][k-1]
                                                else :
                                                    expr ["%ss%s" % (imap[jj],k)] += "%s*(%s)" % (pre,vals[i][1][imap[jj]][k-1])
                                                    exprT["%ss%s" % (imap[jj],k)] += "%s*(%s)" % (pre,vals[i][2][imap[jj]][k-1])
                                # simple vector
                                else :
                                    print 'vector case'
                                    print vals[i]
                                    print expr
                                    quit()

                            # spinor case
                            elif(len(vals[i][1])==4 and "s1" in vals[i][1]) :
                                if(expr=="") :
                                    expr ={}
                                    exprT={}
                                    print vals[i]
                                    for jj in range(1,5) :
                                        if(pre=="") :
                                            expr ["s%s" % jj] = "(%s)" % vals[i][1]["s%s" % jj]
                                            exprT["s%s" % jj] = "(%s)" % vals[i][2]["s%s" % jj]
                                        else :
                                            expr ["s%s" % jj] = "%s*(%s)" % (pre,vals[i][1]["s%s" % jj])
                                            exprT["s%s" % jj] = "%s*(%s)" % (pre,vals[i][2]["s%s" % jj])
                                else :
                                    for jj in range(1,5) :
                                        if(pre=="") :
                                            expr ["s%s" % jj] += "(%s)" % vals[i][1]["s%s" % jj]
                                            exprT["s%s" % jj] += "(%s)" % vals[i][2]["s%s" % jj]
                                        else :
                                            expr ["s%s" % jj] += "%s*(%s)" % (pre,vals[i][1]["s%s" % jj])
                                            exprT["s%s" % jj] += "%s*(%s)" % (pre,vals[i][2]["s%s" % jj])
                            else :
                                print "AAAAAAA"
                                print vertex,len(vals[i][1])
                                quit()
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
                    if(isinstance(expr,basestring)) :
                        expr =expr .replace(oldVal,newVal)
                        exprT=exprT.replace(oldVal,newVal)
                    else :
                        for (key,val) in expr.iteritems() :
                            expr[key]  = val.replace(oldVal,newVal)
                        for (key,val) in exprT.iteritems() :
                            exprT[key] = val.replace(oldVal,newVal)
        # # of particles in the vertex
        vDim = len(vertex.lorentz[0].spins)
        # and momentum components
        for i in range(1,vDim+1) :
            for k in range(0,4) :
                oldVal = "P%s%s*" % (i,imap[k])
                newVal = "P%s.%s()*" % (i,imap[k])
                if(isinstance(expr,basestring)) :
                    expr =expr .replace(oldVal,newVal)
                    exprT=exprT.replace(oldVal,newVal)
                else :
                    for (key,val) in expr.iteritems() :
                        expr[key]  = val.replace(oldVal,newVal)
                    for (key,val) in exprT.iteritems() :
                        exprT[key] = val.replace(oldVal,newVal)

        unit = computeUnit2(dimCheck,vDim)
        if(unit!="") :
            if(isinstance(expr,basestring)) :
                expr  = "(%s)*(%s)" % (expr ,unit)  
                exprT = "(%s)*(%s)" % (exprT,unit)
            else :
                for (key,val) in expr.iteritems() :
                    expr[key]  = "(%s)*(%s)" % (val,unit)  
                for (key,val) in exprT.iteritems() :
                    exprT[key] = "(%s)*(%s)" % (val,unit)  
        # get the coupling for this bit
        val, sym = py2cpp(values[j])
        localCouplings.append("Complex local_C%s = %s;\n" % (j,val))
        symbols |=sym
        if(result=="") :
            if(isinstance(expr,basestring)) :
                if(iloc==0 or vertex.lorentz[0].spins[iloc-1]==1) :
                    result  += " (local_C%s)*Complex(%s) " % (j,expr )
                    resultT += " (local_C%s)*Complex(%s) " % (j,exprT)
                else :
                    result  += " (local_C%s)*(%s) " % (j,expr )
                    resultT += " (local_C%s)*(%s) " % (j,exprT)
            else :
                result ={}
                resultT={}
                for (key,val) in expr.iteritems() :
                    result [key] = " (local_C%s)*Complex(%s) " % (j,val)
                for (key,val) in exprT.iteritems() :
                    resultT[key] = " (local_C%s)*Complex(%s) " % (j,val)
        else :
            if(isinstance(expr,basestring)) :
                if(iloc==0 or vertex.lorentz[0].spins[iloc-1]==1) :
                    result  += " + (local_C%s)*Complex(%s)" % (j,expr )
                    resultT += " + (local_C%s)*Complex(%s)" % (j,exprT)
                else :
                    print result,expr
                    result  += " + (local_C%s)*(%s)" % (j,expr )
                    resultT += " + (local_C%s)*(%s)" % (j,exprT)
            else :
                for (key,val) in expr.iteritems() :
                    result [key] += " + (local_C%s)*Complex(%s) " % (j,val)
                for (key,val) in exprT.iteritems() :
                    resultT[key] += " + (local_C%s)*Complex(%s) " % (j,val)
    # final defns for more complicated types
    print "TESTING RESULT !!!",result
    if(not isinstance(result,basestring)) :
        subs ={}
        for (key,val) in result.iteritems() :
            subs["out%s" % key]= val
        subsT={}
        for (key,val) in resultT.iteritems() :
            subsT["out%s" % key] = val
        
        if("ts1" in result) :
            stype  = "LorentzRSSpinor"
            sbtype = "LorentzRSSpinorBar"
            if(offType.find("Bar")>0) : (stype,sbtype)=(sbtype,stype)
            subs["type"] = stype      
            result  = Template("${type}<double>(${outxs1},\n${outxs2},\n${outxs3},\n${outxs4},\n${outys1},\n${outys2},\n${outys3},\n${outys4},\n${outzs1},\n${outzs2},\n${outzs3},\n${outzs4},\n${outts1},\n${outts2},\n${outts3},\n${outts4})").substitute(subs)
            subsT["type"] = sbtype
            resultT = Template("${type}<double>(${outxs1},\n${outxs2},\n${outxs3},\n${outxs4},\n${outys1},\n${outys2},\n${outys3},\n${outys4},\n${outzs1},\n${outzs2},\n${outzs3},\n${outzs4},\n${outts1},\n${outts2},\n${outts3},\n${outts4})").substitute(subsT)
        elif("s1" in result) :
            stype  = "LorentzSpinor"
            sbtype = "LorentzSpinorBar"
            if(offType.find("Bar")>0) : (stype,sbtype)=(sbtype,stype)
            if(not RS) : sbtype=stype
            print "TESTING !!! ",offType,stype,sbtype
            subs["type"] = stype
            result  = Template("${type}<double>(${outs1},\n${outs2},\n${outs3},\n${outs4})").substitute(subs)
            subsT["type"] = sbtype
            resultT  = Template("${type}<double>(${outs1},\n${outs2},\n${outs3},\n${outs4})").substitute(subsT)
        else :
            print result
            print 'type problem'
            quit()
    # multiple by scalar wavefunctions
    scalars=""
    for i in range (0,len(vertex.lorentz[0].spins)) :
        if(vertex.lorentz[0].spins[i]==1 and i+1!=iloc) :
            scalars += "sca%s*" % (i+1)
    if(scalars!="") :
        result  = "(%s)*(%s)" % (result ,scalars[0:-1])
        resultT = "(%s)*(%s)" % (resultT,scalars[0:-1])
    # vertex, just return the answer
    if(iloc==0) :
        if(hasTranspose and not RS) :
            result = vTemplateT.format(iloc=fermions[0],id=vertex.particles[fermions[0]-1].pdg_code,
                                      cf=py2cpp(cf[0])[0],res=result,resT=resultT)
        else :
            result  = "return (%s)*(%s);\n" % (result ,py2cpp(cf[0])[0])
        if(RS) :
            resultT = "return (%s)*(%s);\n" % (resultT,py2cpp(cf[0])[0])
    # off-shell particle
    else :
        # off-shell scalar
        if(vertex.lorentz[0].spins[iloc-1] == 1 ) :
            if(hasTranspose and not RS) :
                result = scaTTemplate.format(iloc=iloc,res=result,resT=resultT,isp=fermions[0],
                                             id=vertex.particles[fermions[0]-1].pdg_code,cf=py2cpp(cf[0])[0])
            else :
                result = scaTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],res=result)
            if(RS) :
                resultT = scaTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],res=resultT)
        # off-shell fermion
        elif(vertex.lorentz[0].spins[iloc-1] == 2 ) :
            if(hasTranspose and not RS) :
                result = sTTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offType.replace("WaveFunction",""),
                                           res=result.replace( "M%s" % iloc, "mass" ),resT=resultT.replace( "M%s" % iloc, "mass" ),
                                           offTypeB=offType,id=vertex.particles[iloc-1].pdg_code)
            else :
                result = sTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offType.replace("WaveFunction",""),
                                          res=result.replace( "M%s" % iloc, "mass" ),offTypeB=offType)
            if(RS) :
                if(offType.find("Bar")>0) :
                    offTypeT=offType.replace("Bar","")
                else :
                    offTypeT=offType.replace("Spinor","SpinorBar")
                resultT = sTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offTypeT.replace("WaveFunction",""),
                                           res=resultT.replace( "M%s" % iloc, "mass" ),offTypeB=offTypeT)
        # off-shell vector
        elif(vertex.lorentz[0].spins[iloc-1] == 3 ) :
            if(hasTranspose and not RS) :
                result = vecTTemplate.format(iloc=iloc,res=result,resT=resultT,isp=fermions[0],
                                             id=vertex.particles[fermions[0]-1].pdg_code,
                                             cf=py2cpp(cf[0])[0])
            else :
                result = vecTemplate.format(iloc=iloc,res=result,cf=py2cpp(cf[0])[0])
        elif(vertex.lorentz[0].spins[iloc-1] == 4 ) :
            result = RSTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offType.replace("WaveFunction",""),
                                       res=result.replace( "M%s" % iloc, "mass" ),offTypeB=offType)
            if(offType.find("Bar")>0) :
                offTypeT=offType.replace("Bar","")
            else :
                offTypeT=offType.replace("Spinor","SpinorBar")
            resultT = RSTemplate.format(iloc=iloc,cf=py2cpp(cf[0])[0],offTypeA=offTypeT.replace("WaveFunction",""),
                                        res=resultT.replace( "M%s" % iloc, "mass" ),offTypeB=offTypeT)
            


            
    # check if momenta defns needed to clean up compile of code
    for (key,val) in defns.iteritems() :
        if( isinstance(key, basestring)) :
            if(key.find("vvP")==0) :
                momenta[int(key[3])-1][0] = True
        else :
            for vals in key :
                if(vals[0]=="P") :
                    momenta[int(vals[1])-1][0] = True
    sorder=swapOrder(vertex,iloc,momenta)
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

    # special for tgranspose in the RS case
    if(hasTranspose and RS) :
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
                                       result=resultT,swap=sorder)
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
        print 'testing evaluate multiple porblem',spin
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
