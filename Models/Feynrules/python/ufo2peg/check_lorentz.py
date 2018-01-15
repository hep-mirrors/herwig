import itertools,cmath,re,sys
from .helpers import SkipThisVertex,extractAntiSymmetricIndices,def_from_model
from .converter import py2cpp
from .lorentzparser import parse_lorentz
import sympy,string
from string import Template
from sympy import Matrix,Symbol

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

def processTensorCouplings(lorentztag,vertex,model,parmsubs,all_couplings) :
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
                        code = abs(vertex.particles[0].pdg_code)
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
                        test[i] = '(%s)/%s**2' % (test[i],vertex.particles[0].mass.value)
            # fermions divide by 4*m
            elif(ix==6 and lorentztag=="FFT" and
                 float(vertex.particles[0].mass.value) != 0. ) :
                for i in range(0,len(test)) :
                    if(test[i]) :
                        test[i] = '-(%s)/%s/4' % (test[i],vertex.particles[0].mass.value)
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
                           float(vertex.particles[0].mass.value) == 0. ) :
                            continue
                        # special for vector gauge terms
                        if(lorentztag=="VVT" and ix>=13) :
                            continue
                        raise SkipThisVertex()
                    tval2 = evaluate(test[i],model,parmsubs)
                    if(abs(tval[i]-tval2)>1e-6) :
                        # special for fermion mass term if fermions massless
                        if(lorentztag=="FFT" and ix ==6 and tval2 == 0. and
                           float(vertex.particles[0].mass.value) == 0. ) :
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
                    prepend = kinematicsline.format(id1=vertex.particles[0].pdg_code,
                                                     id2=vertex.particles[1].pdg_code,
                                                     id3=vertex.particles[2].pdg_code,
                                                     kine=output)
                else :
                    prepend = kinematicsline2.format(id1=vertex.particles[0].pdg_code,
                                                      id2=vertex.particles[1].pdg_code,
                                                      id3=vertex.particles[2].pdg_code,
                                                      id4=vertex.particles[2].pdg_code,
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
                sys.stderr.write(
                    'Warning: unsupported {tag} ( {ps} ) Lorentz structure in {name}:\n'
                    .format(tag="VVVV", name=vertex.name, ps=' '.join(map(str,vertex.particles)))
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
                            value = all_couplings[icolor][0]
                            tval  = evaluate(value,model,parmsubs)
                        tval2 = evaluate(all_couplings[icolor][0],model,parmsubs)
                        tval3 = -evaluate(all_couplings[icolor][1],model,parmsubs)
                    elif(col[0][0]==1 and col[0][1]==3 and col[1][0] ==2 and col[1][1] == 4) : 
                        if(all_couplings[icolor][1] or not all_couplings[icolor][0] or
                           not all_couplings[icolor][2]) :
                            raise SkipThisVertex()
                        if(not value) :
                            value = all_couplings[icolor][0]
                            tval  = evaluate(value,model,parmsubs)
                        tval2 = evaluate(all_couplings[icolor][0],model,parmsubs)
                        tval3 = -evaluate(all_couplings[icolor][2],model,parmsubs)
                    elif(col[0][0]==1 and col[0][1]==4 and col[1][0] ==2 and col[1][1] == 3) : 
                        if(all_couplings[icolor][0] or not all_couplings[icolor][1] or
                           not all_couplings[icolor][2]) :
                            raise SkipThisVertex()
                        if(not value) :
                            value = all_couplings[icolor][1]
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

def processFermionCouplings(lorentztag,vertex,model,parmsubs,all_couplings) :
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
                  % vertex.particles[0].pdg_code)
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
        return b.value-a.value
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
        output[0].name="int"
        output[0].value=structure[0]+"1."
        output[0].value=float(output[0].value)
        structure=structure[2:-1]
    elif(structure[0]=="(") :
        temp=structure.rsplit(")",1)
        structure=temp[0][1:]
        output.append(LorentzStructure())
        output[0].name="int"
        output[0].value="1."+temp[1]
        output[0].value=float(eval(output[0].value))
    # split up the structure
    for struct in structure.split("*"):
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
                output[0].name="int"
            except :
                print struct
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
                    if(j==k) : continue
                    for i in range(0,len(parsed[l][k].lorentz)) :
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
                   L.spins[parsed[l][j].lorentz[i]-1]==3 ) :
                    parsed[l][j].lorentz[i]=  "E%s" % parsed[l][j].lorentz[i]
                    if(parsed[l][j].name=="P") :
                        parsed[l][j].lorentz[1] = "P%s" % parsed[l][j].lorentz[1]
                        parsed[l][j].name="Metric"
        parsed[l] = [x for x in parsed[l] if x != ""]

def computeUnit(dimension) :
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
            else :
                output = "1./GeV%s" % totalDim
        else :
            if(totalDim==-1) :
                output = "GeV"
            else :
                output = "GeV%s" % (-totalDim)
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

def generateVertex(iloc,L,parsed,lorentztag,order,defns) :
    # various strings for matrixes
    I4 = "Matrix([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])"
    G5 = "Matrix([[-1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,1]])"
    PM = "Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]])"
    PP = "Matrix([[0,0,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,1]])"
    vslash  = Template("Matrix([[0,0,${v}tmz,-${v}xmy],[0,0,-${v}xpy,${v}tpz],[${v}tpz,${v}xmy,0,0],[${v}xpy,${v}tmz,0,0]])")
    vslashS = Template("${v}tpz=Symbol(\"${v}tpz\")\n${v}tmz=Symbol(\"${v}tmz\")\n${v}xpy=Symbol(\"${v}xpy\")\n${v}xmy=Symbol(\"${v}xmy\")\n")
    vslashD = Template("complex<${var}> ${v}tpz = ${v}.t()+${v}.z();\n    complex<${var}> ${v}tmz = ${v}.t()-${v}.z();\n    complex<${var}> ${v}xpy = ${v}.x()+Complex(0.,1.)*${v}.y();\n    complex<${var}> ${v}xmy = ${v}.x()-Complex(0.,1.)*${v}.y();")
    vslashM  = Template("Matrix([[$m,0,${v}tmz,-${v}xmy],[0,$m,-${v}xpy,${v}tpz],[${v}tpz,${v}xmy,$m,0],[${v}xpy,${v}tmz,0,$m]])")
    vslashM2 = Template("Matrix([[$m,0,-${v}tmz,${v}xmy],[0,$m,${v}xpy,-${v}tpz],[-${v}tpz,-${v}xmy,$m,0],[-${v}xpy,-${v}tmz,0,$m]])")
    vslashMS = Template("${v}tpz=Symbol(\"${v}tpz\")\n${v}tmz=Symbol(\"${v}tmz\")\n${v}xpy=Symbol(\"${v}xpy\")\n${v}xmy=Symbol(\"${v}xmy\")\n${m}=Symbol(\"${m}\")\n")
    dirac=["Matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])","Matrix([[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]])",
           "Matrix([[0,0,0,complex(0, -1)],[0,0,complex(0, 1),0],[0,complex(0, 1),0,0],[complex(0, -1),0,0,0]])",
           "Matrix([[0,0,1,0],[0,0,0,-1],[-1,0,0,0],[0,1,0,0]])"]
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
                    if(ind1[1]=="P") : dimension[i][2] +=1
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
                        print "problem"
                        quit()
                    if(int(ind1[1])==iloc) :
                        output[i] += "*(%s)" % ind2
                    elif(int(ind2[1])==iloc) :
                        output[i] += "*(%s)" % ind1
                    else :
                        defns[(ind1,ind2)] = [name,"complex<double> %s = %s*%s;" % (name,ind1,ind2)]
                        output[i] += "*(%s)" % name
    # remove any (now) empty elements
    for i in range (0,len(parsed)) :
        parsed[i] = [x for x in parsed[i] if x != ""]
    # now for gamma matrix strings
    if(lorentztag[0]=="F") :
        result=""
        for i in range(0,len(parsed)): 
            if(len(parsed[i])==0) :
                continue
            if(len(parsed[i])==1 and parsed[i][0]=="") :
                result+=parsed[i][0]
                continue
            # first and last lorentz indices
            sind=parsed[i][0].spin[0]
            lind=0
            # first piece of the expression we need to evaluate
            dtemp=[0,0,0]
            if(sind==iloc) :
                expr = vslashM2.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
                Symbols = vslashMS.substitute({ "v" : "P%s" % sind, "m" : "M%s" % sind} )
                defns["vvP%s" % sind ] = ["vvP%s" % sind ,
                                          vslashD.substitute({ "var" : "Energy",
                                                               "v" :  "P%s" % sind })]
                dtemp[1]+=1
                dtemp[0]+=0.5
            else :
                Symbols=Template("sbar${s}s1=Symbol(\"sbar${s}s1\")\nsbar${s}s2=Symbol(\"sbar${s}s2\")\nsbar${s}s3=Symbol(\"sbar${s}s3\")\nsbar${s}s4=Symbol(\"sbar${s}s4\")\n").substitute({'s' :parsed[i][0].spin[0]})
                expr=Template("Matrix([[sbar${s}s1,sbar${s}s2,sbar${s}s3,sbar${s}s4]])").substitute({'s' : sind })
                dtemp[0]+=0.5
            # parse the remaining structures
            for j in range(0,len(parsed[i])) :
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
                    if(parsed[i][j].lorentz[0] == ("E%s" % iloc ) ) :
                        expr += "*DUMMY"
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
            # last piece of spin chain
            if(lind==iloc) :
                expr += "*"+vslashM.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
                Symbols += vslashMS.substitute({ "v" : "P%s" % lind, "m" : "M%s" % lind} )
                defns["vvP%s" % lind ] = ["vvP%s" % lind ,
                                          vslashD.substitute({ "var" : "Energy",
                                                               "v" :  "P%s" % lind })]
                dtemp[1] += 1
                dtemp[0] += 0.5
            else :
                expr+=Template("*Matrix([[s${s}s1],[s${s}s2],[s${s}s3],[s${s}s4]])").substitute({'s' : lind})
                Symbols+=Template("s${s}s1=Symbol(\"s${s}s1\")\ns${s}s2=Symbol(\"s${s}s2\")\ns${s}s3=Symbol(\"s${s}s3\")\ns${s}s4=Symbol(\"s${s}s4\")\n").substitute({'s' : lind})
                dtemp[0] += 0.5
            parsed[i] = [x for x in parsed[i] if x != ""]
            if(len(parsed[i])!=0) :
                print "ERROR"
                quit()
            # off-shell vector
            if(expr.find("DUMMY")>=0) :
                vtemp=[]
                defns["I"] = ["I","static Complex I(0.,1.);"]
                for matrix in dirac :
                    temp={}
                    exec("import sympy\nfrom sympy import Symbol,Matrix\n"+Symbols+"result="+expr.replace("DUMMY",matrix)) in temp
                    vtemp.append(temp["result"][0,0])
                unit = computeUnit(dtemp)
                output[i] += "*LorentzVector<complex<%s> >(%s,%s,%s,%s)" % (unit,vtemp[1],vtemp[2],vtemp[3],vtemp[0])
            else :
                temp={}
                if(iloc==0 or (iloc!=sind and iloc!=lind)) :
                    exec("import sympy\nfrom sympy import Symbol,Matrix\n"+Symbols+"result="+expr) in temp
                    output[i] += "*(%s)" % temp["result"][0,0]
                else :
                    exec("import sympy\nfrom sympy import Symbol,Matrix\n"+Symbols+"result="+expr) in temp
                    unit = computeUnit(dtemp)
                    if(iloc==sind) :
                        output[i] += "*LorentzSpinor<%s>(%s,%s,%s,%s)" % (unit,temp["result"][0],temp["result"][1],
                                                                          temp["result"][2],temp["result"][3])
                    else :
                        output[i] += "*LorentzSpinorBar<%s >(%s,%s,%s,%s)" % (unit,temp["result"][0],temp["result"][1],
                                                                             temp["result"][2],temp["result"][3])
            dimension[i] = list(map(lambda x, y: x + y, dtemp, dimension[i]))
    return (output,dimension)
            
def convertLorentz(Lstruct,lorentztag,order,iloc,defns,evalVertex) :
    # split the structure into individual terms
    structures=Lstruct.structure.split()
    parsed=[]
    for struct in structures :
        parsed.append(parse_structure(struct))
    # convert lorentz contractions to dot products
    contractLorentz(Lstruct,parsed,lorentztag,order)
    # now in a position to generate the code
    evalVertex.append(generateVertex(iloc,Lstruct,parsed,lorentztag,order,defns))

evaluateTemplate = """\
{decl} {{
    {momenta}
    {waves}
    {defns}
    {symbols}
    {couplings}
    {result}
}}
"""
    
def generateEvaluateFunction(model,vertex,iloc,values,defns,vertexEval,cf) :
    # first construct the signature of the function
    iferm=0
    decls=[]
    offType="Complex"
    momenta=[]
    waves=[]
    poff=""
    fermionReplace=[]
    for i in range(0,len(vertex.lorentz[0].spins)) :
        spin = vertex.lorentz[0].spins[i]
        if(i+1==iloc) :
            if(spin==1) :
                offType="ScalarWaveFunction"
            elif(spin==2) :
                if(iferm==0) :
                    offType="SpinorBarWaveFunction"
                else :
                    offType="SpinorWaveFunction"
                iferm+=1
            elif(spin==3) :
                offType="VectorWaveFunction"
            elif(spin==4) :
                offType="TensorWaveFunction"
            else :
                print 'unknown spin',spin
                quit()
            poff = ("Lorentz5Momentum P%s = " % (i+1) ) + poff
        else :
            if(spin==1) :
                decls.append("ScalarWaveFunction & scaW%s" % (i+1))
                momenta.append("Lorentz5Momentum P%s = scaW%s.momentum();" % (i+1,i+1))
                waves.append("Complex sca%s = scaW%s.wave();" % (i+1,i+1))
            elif(spin==2) :
                if(iferm==0) :
                    decls.append("SpinorWaveFunction & sW%s" % (i+1))
                    momenta.append("Lorentz5Momentum P%s = sW%s.momentum();" % (i+1,i+1))
                    waves.append("LorentzSpinor<double> s%s = sW%s.wave();" % (i+1,i+1))
                    fermionReplace.append("s%s"%(i+1))
                else :
                    decls.append("SpinorBarWaveFunction & sbarW%s" % (i+1))
                    momenta.append("Lorentz5Momentum P%s = sbarW%s.momentum();" % (i+1,i+1))
                    waves.append("LorentzSpinorBar<double> sbar%s = sbarW%s.wave();" % (i+1,i+1))
                    fermionReplace.append("sbar%s"%(i+1))
                iferm +=1
            elif(spin==3) :
                decls.append("VectorWaveFunction & vW%s" % (i+1))
                momenta.append("Lorentz5Momentum P%s = vW%s.momentum();" % (i+1,i+1))
                waves.append("LorentzPolarizationVector E%s = vW%s.wave();" % (i+1,i+1))
            elif(spin==4) :
                decls.append("TensorWaveFunction & tW%s" % (i+1))
                momenta.append("Lorentz5Momentum P%s = tW%s.momentum();" % (i+1,i+1))
                waves.append("LorentzTensor t%s = tW%s.wave();" % (i+1,i+1))
            else :
                print 'unknown spin',spin
                quit()
            poff += "-P%s" % (i+1)
    sig=""
    if(iloc==0) :
        sig="%s evaluate(Energy2, const %s)" % (offType,", const ".join(decls))
    else :
        sig="%s evaluate(Energy2, int iopt, tcPDPtr out, const %s, complex<Energy> mass=-GeV, complex<Energy> width=-GeV)" % (offType,", const ".join(decls))
        momenta.append(poff+";")
    # cat the definitions
    defString=""
    for (key,value) in defns.iteritems() :
        defString+="    %s\n" %value[1]

    oval=""
    symbols=set()
    localCouplings=[]
    result=""
    for j in range(0,len(vertexEval)) :
        (vals,dim) = vertexEval[j]
        expr=""
        dimCheck=dim[0]
        for i in range(0,len(vals)) :
            if(vals[i]=="+" or vals[i]=="-") :
                expr +=vals[i]
            else :
                if(dimCheck[0]!=dim[i][0] or dimCheck[1]!=dim[i][1] or
                   dimCheck[2]!=dim[i][2]) :
                    print "DIMENSION PROBLEM",dimCheck,dim[i],vertex
                    quit()
                expr += "(%s)" % vals[i]
        for rep in fermionReplace :
            for i in range(1,5) :
                oldVal = "%ss%s" % (rep,i)
                newVal = "%s.s%s()" % (rep,i)
                expr=expr.replace(oldVal,newVal)
        vDim = len(vertex.lorentz[0].spins)

        unit = computeUnit2(dimCheck,vDim)
        if(unit!="") :
            expr = "(%s)*(%s)" % (expr,unit)  
        val, sym = py2cpp(values[j])
        localCouplings.append("Complex local_C%s = %s;\n" % (j,val))
        symbols |=sym
        if(result!="") :
            result += " + (local_C%s)*(%s)" % (j,expr)
        else :
            result += " (local_C%s)*(%s) " % (j,expr)
    if(iloc==0) :
        result = "return (%s)*(%s);\n" % (result,py2cpp(cf[0])[0])
    else :
        if(vertex.lorentz[0].spins[iloc-1] == 3 ) :
            result = "LorentzPolarizationVector vtemp = %s;\n" % result
            result +="    Energy2 p2 = P%s.m2();\nComplex fact = -Complex(0.,1.)*(%s)*propagator(iopt,p2,out,mass,width);\n    if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();\n    complex<Energy2> mass2 = sqr(mass);\n    if(mass.real()==ZERO) {\n vtemp =fact*vtemp;\n}\n    else {\ncomplex<Energy> dot = P%s*vtemp;\n vtemp = fact*(vtemp-dot/mass2*P%s);\n}\nreturn VectorWaveFunction(-P%s,out,vtemp.x(),vtemp.y(),vtemp.z(),vtemp.t());\n" % (iloc,py2cpp(cf[0])[0],iloc,iloc,iloc)
        elif(vertex.lorentz[0].spins[iloc-1] == 2 ) :
            result = "if(mass.real() < ZERO) mass  = (iopt==5) ? ZERO : out->mass();\n     Energy2 p2 = P%s.m2();\n    Complex fact = Complex(0.,1.)*(%s)*propagator(iopt,p2,out,mass,width);\n Lorentz%s<double> newSpin = fact*(%s);\n    return %s(-P%s,out,newSpin.s1(),newSpin.s2(),newSpin.s3(),newSpin.s4());" % \
                                                         (iloc,py2cpp(cf[0])[0],offType.replace("WaveFunction",""),result.replace( "M%s" % iloc, "mass" ),offType,iloc)
            
    header="virtual %s" % sig
    sig=sig.replace("=-GeV","")
    symboldefs = [ def_from_model(model,s) for s in symbols ]
    function = evaluateTemplate.format(decl=sig,momenta="\n    ".join(momenta),defns=defString,
                                       waves="\n    ".join(waves),symbols='\n    '.join(symboldefs),
                                       couplings="\n    ".join(localCouplings),
                                       result=result)
    return (header,function)
