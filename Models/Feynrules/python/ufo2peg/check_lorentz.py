from __future__ import print_function
import itertools,cmath,re
from .helpers import SkipThisVertex,extractAntiSymmetricIndices
from .converter import py2cpp
from .lorentzparser import parse_lorentz
import string,re

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
            j = i+3*(i//3)
            k = i+3+3*(i//3)
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
    tval  = ["Unknown"]*3
    value = ["Unknown"]*3
    # loop over the colours
    for icolor in range(0,len(all_couplings)) :
        lmax = len(all_couplings[icolor])
        if(fermions) : lmax //=3
        # loop over the different terms
        for ix in range(0,lmax) :
            test = [False]*3
            imax=3
            # normal case
            if( not fermions ) :
                test[0] = all_couplings[icolor][ix]
                imax=1
            else :
                # first case vector but no L/R couplings
                if( not all_couplings[icolor][lmax+ix]  and
                    not all_couplings[icolor][2*lmax+ix] ) :
                    test[0] = all_couplings[icolor][ix]
                    imax=1
                    # special for mass terms and massless particles
                    if(not all_couplings[icolor][ix]) :
                        code = abs(vertex.particles[order[0]-1].pdg_code)
                        if(ix==6 and (code ==12 or code ==14 or code==16) ) :
                            continue
                        else :
                            raise SkipThisVertex()
                # second case L/R couplings
                elif( not all_couplings[icolor][ix] ) :
                    # L=R, replace with vector
                    if(all_couplings[icolor][lmax+ix] ==
                       all_couplings[icolor][2*lmax+ix]) :
                        test[0]  = all_couplings[icolor][lmax+ix]
                        imax=1
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
            if((tval[0]=="Unknown" and fermions ) or
               (not fermions and tval[0]=="Unknown" and tval[1]=="Unknown" and tval[2]=="Unknown")) :
                value = test
                for i in range(0,len(test)) :
                    if(test[i]) : tval[i] = evaluate(test[i],model,parmsubs)
            else :
                for i in range(0,imax) :
                    if(not test[i] and not tval[i]) :
                        continue
                    # special for vector gauge terms
                    if(lorentztag=="VVT" and ix>=13) :
                        continue
                    if(not test[i] or tval[i]=="Unknown") :
                        # special for mass terms and vectors
                        if(lorentztag=="VVT" and ix >=10 and ix <=12 and
                           float(vertex.particles[order[0]-1].mass.value) == 0. ) :
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
    structure1 = L.structure.replace(")-P",") - P").replace(")+P",") + P")
    structure1 = structure1.split()
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
            value = all_couplings[icolor][0]
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
                    if all_couplings[icolor][ix] == False:
                        raise SkipThisVertex()
                    else:
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
                    print('unknown colour structure for VVVV vertex')
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
