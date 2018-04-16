import sys,pprint

from .helpers import CheckUnique,getTemplate,writeFile,qcd_qed_orders,def_from_model
from .converter import py2cpp
from .collapse_vertices import collapse_vertices
from .check_lorentz import tensorCouplings,VVVordering,lorentzScalar,\
    processTensorCouplings,scalarCouplings,processScalarCouplings,scalarVectorCouplings,\
    processScalarVectorCouplings,vectorCouplings,processVectorCouplings,fermionCouplings,processFermionCouplings,\
    RSCouplings
from .helpers import SkipThisVertex,extractAntiSymmetricIndices,isGoldstone

# prefactors for vertices
lfactors = { 
    'FFV'  : '-complex(0,1)',  # ok
    'VVV'  : 'complex(0,1)',   # changed to fix ttbar
    'VVVS' : 'complex(0,1)',   # should be as VVV
    'VVVV' : 'complex(0,1)',
    'VVS'  : '-complex(0,1)',
    'VSS'  : '-complex(0,1)', # changed to minus to fix dL ->x1 W- d
    'SSS'  : '-complex(0,1)',  # ok
    'VVSS' : '-complex(0,1)',  # ok
    'VVT'  : 'complex(0,2)',
    'VVVT' : '-complex(0,2)',
    'SSSS' : '-complex(0,1)',  # ok
    'FFS'  : '-complex(0,1)',  # ok
    'SST'  : 'complex(0,2)',
    'FFT'  : '-complex(0,8)',
    'FFVT' : '-complex(0,4)',
    'RFS'  : 'complex(0,1)',
    'RFV'  : 'complex(0,1)',
}

# template for the header for a vertex
VERTEXHEADER = """\
#include "ThePEG/Helicity/Vertex/{spindirectory}/{lorentztag}Vertex.h"
"""

# template for the implmentation for a vertex
VERTEXCLASS = getTemplate('Vertex_class')

# template for the .cc file for vertices
VERTEX = getTemplate('Vertex.cc')

vertexline = """\
create Herwig::{modelname}V_{vname} /Herwig/{modelname}/V_{vname}
insert {modelname}:ExtraVertices 0 /Herwig/{modelname}/V_{vname}
"""


def get_lorentztag(spin):
    """Produce a ThePEG spin tag for the given numeric FR spins."""
    spins = { 1 : 'S', 2 : 'F', 3 : 'V', 4 : 'R', 5 : 'T', -1 : 'U' }
    result=[]
    for i in range(0,len(spin)) :
        result.append((spins[spin[i]],i+1))
    def spinsort(a,b):
        """Helper function for ThePEG's FVST spin tag ordering."""
        (a1,a2) = a
        (b1,b2) = b
        if a1 == b1: return 0
        for letter in 'URFVST':
            if a1 == letter: return -1
            if b1 == letter: return  1

    result = sorted(result, cmp=spinsort)
    order=[]
    output=""
    for i in range(0,len(result)) :
        (a,b) = result[i]
        order.append(b)
        output+=a
    return (output,order)

def unique_lorentztag(vertex):
    """Check and return the Lorentz tag of the vertex."""
    unique = CheckUnique()
    for l in vertex.lorentz:
        (lorentztag,order) = get_lorentztag(l.spins)
        unique( lorentztag )
        lname = l.name[:len(lorentztag)]
        if sorted(lorentztag) != sorted(lname):
            raise Exception("Lorentztags: %s is not %s in %s" 
                            % (lorentztag,lname,vertex))
    return (lorentztag,order)

def colors(vertex) :
    try:
        unique = CheckUnique()
        for pl in vertex.particles_list:
            struct = [ p.color for p in pl ]
            unique(struct)
    except:
        struct = [ p.color for p in vertex.particles ]
    pos = colorpositions(struct)
    L = len(struct)
    return (L,pos)

def colorfactor(vertex,L,pos,lorentztag):
    def match(patterns,color=vertex.color):
        result = [ p == t
                   for p,t in zip(patterns,color) ]
        return all(result)

    label = None
    l = lambda c: len(pos[c])
    if l(1) == L:
        label = ('1',)
        if match(label): return ('1',)

    elif l(3) == l(-3) == 1 and l(1) == L-2:
        nums = [pos[3][0], pos[-3][0]]
        label = ('Identity({0},{1})'.format(*sorted(nums)),)
        if match(label): return ('1',)

    elif l(6) == l(-6) == 1 and l(1) == L-2:
        nums = [pos[6][0], pos[-6][0]]
        label = ('Identity({0},{1})'.format(*sorted(nums)),)
        if match(label): return ('1',)

    elif l(6) == l(-6) == 2 and L==4:
        sys.stderr.write(
            'Warning: Unknown colour structure 6 6 6~ 6~ ( {ps} ) in {name}.\n'
            .format(name=vertex.name, ps=' '.join(map(str,vertex.particles)))
        )
        raise SkipThisVertex()

    elif l(8) == l(3) == l(-3) == 1 and l(1) == L-3:
        label = ('T({g},{q},{qb})'.format(g=pos[8][0],q=pos[3][0],qb=pos[-3][0]),)
        if match(label): return ('1',)

    elif l(8) == l(6) == l(-6) == 1 and l(1) == L-3:
        label = ('T6({g},{s},{sb})'.format(g=pos[8][0],s=pos[6][0],sb=pos[-6][0]),)
        if match(label): return ('1',)

    elif l(6) == 1 and l(-3) == 2 and L==3:
        label = ('K6({s},{qb1},{qb2})'.format(s=pos[6][0],qb1=pos[-3][0],qb2=pos[-3][1]),)
        if match(label): return ('1',)

    elif l(-6) == 1 and l(3) == 2 and L==3:
        label = ('K6Bar({sb},{q1},{q2})'.format(sb=pos[-6][0],q1=pos[3][0],q2=pos[3][1]),)
        if match(label): return ('1',)

    elif l(3) == L == 3:
        colors=[]
        for color in vertex.color :
            order,sign  = extractAntiSymmetricIndices(color,"Epsilon(")
            colors.append("Epsilon(%s,%s,%s)" % (order[0],order[1],order[2]))
        label = ('Epsilon(1,2,3)',)
        if match(label,colors): return ('1',) # TODO check factor!

    elif l(-3) == L == 3:
        colors=[]
        for color in vertex.color :
            order,sign  = extractAntiSymmetricIndices(color,"EpsilonBar(")
            colors.append("Epsilon(%s,%s,%s)" % (order[0],order[1],order[2]))
        label = ('EpsilonBar(1,2,3)',)
        if match(label): return ('1',) # TODO check factor!

    elif l(8) == L == 3:
        colors=[]
        for color in vertex.color :
            order,sign  = extractAntiSymmetricIndices(color,"f(")
            colors.append("f(%s,%s,%s)" % (order[0],order[1],order[2]))
        # if lorentz is FFV get extra minus sign
        if lorentztag in ['FFV'] : sign *=-1
        label = ('f(1,2,3)',)
        if match(label,colors): return ('-complex(0,1)*(%s)'%sign,)

    elif l(8) == L == 4:
        colors=[]
        for color in vertex.color :
            f = color.split("*")
            (o1,s1) = extractAntiSymmetricIndices(f[0],"f(")
            (o2,s2) = extractAntiSymmetricIndices(f[1],"f(")
            if(o2[0]<o1[0]) : o1,o2=o2,o1
            colors.append("f(%s)*f(%s)" % (",".join(o1),",".join(o2)))
        def coloursort(a,b) :
            if a == b: return 0
            i1=int(a[4])
            i2=int(b[4])
            if(i1==i2)  : return 0
            elif(i1<i2) : return -1
            else        : return 1
        colors=sorted(colors,cmp=coloursort)
        label = ('f(1,2,-1)*f(3,4,-1)',
                 'f(1,3,-1)*f(2,4,-1)',
                 'f(1,4,-1)*f(2,3,-1)')
        nmatch=0
        for c1 in label:
            for  c2 in colors :
                if(c1==c2) : nmatch+=1
        if(nmatch==2 and lorentztag=="VVSS") :
            return ('1','1')
        elif(nmatch==3 and lorentztag=="VVVV") :
            return ('-1.','-1.','-1.')

    elif l(8) == 2 and l(3) == l(-3) == 1 and L==4:
        subs = {
            'g1' : pos[8][0],
            'g2' : pos[8][1],
            'qq' : pos[3][0],
            'qb' : pos[-3][0] 
        }
        label = ('T({g1},-1,{qb})*T({g2},{qq},-1)'.format(**subs),
                 'T({g1},{qq},-1)*T({g2},-1,{qb})'.format(**subs))
        if match(label): return ('1.','1.')
        
    elif l(8) == 2 and l(6) == l(-6) == 1 and L==4:
        subs = {
            'g1' : pos[8][0],
            'g2' : pos[8][1],
            'qq' : pos[6][0],
            'qb' : pos[-6][0] 
        }
        label = ('T6({g1},-1,{qb})*T6({g2},{qq},-1)'.format(**subs),
                 'T6({g1},{qq},-1)*T6({g2},-1,{qb})'.format(**subs))
        if match(label): return ('1.','1.')

    elif l(8) == 2 and l(8)+l(1)==L :
        subs = { 'g1' : pos[8][0], 'g2' : pos[8][1] }
        label = ('Identity({g1},{g2})'.format(**subs),)
        if match(label) : return ('1.',)

    elif l(8) == 3 and l(1)==1 and L==4 :
        colors=[]
        for color in vertex.color :
            order,sign  = extractAntiSymmetricIndices(color,"f(")
            colors.append("f(%s,%s,%s)" % (order[0],order[1],order[2]))
        label = ('f(1,2,3)',)
        if match(label,colors): return ('-complex(0,1)*(%s)'%sign,)

    sys.stderr.write(
        "Warning: Unknown colour structure {color} ( {ps} ) in {name}.\n"
        .format(color = ' '.join(vertex.color), name = vertex.name,
                ps = ' '.join(map(str,vertex.particles)))
    )
    raise SkipThisVertex()

def colorpositions(struct):
    positions = { 
        1 : [],
        3 : [],
        -3 : [],
        6 : [],
        -6 : [],
        8 : [],
    }
    for i,s in enumerate(struct,1):
        positions[s].append(i)
    return positions

def spindirectory(lt):
    """Return the spin directory name for a given Lorentz tag."""
    if 'T' in lt:
        spin_directory = 'Tensor'
    elif 'S' in lt:
        spin_directory = 'Scalar'
    elif 'V' in lt:
        spin_directory = 'Vector'
    else:
        raise Exception("Unknown Lorentz tag {lt}.".format(lt=lt))
    return spin_directory

def write_vertex_file(subs):
    'Write the .cc file for some vertices'
    newname = '%s_Vertices_%03d.cc' % (subs['ModelName'],subs['vertexnumber'])
    subs['newname'] = newname
    writeFile( newname, VERTEX.substitute(subs) )
    
def checkGhostGoldstoneVertex(lorentztag,vertex) :
    'check if vertex has ghosts or goldstones'
    # remove vertices involving ghost fields
    if 'U' in lorentztag:
        return True
    # remove vertices involving goldstones
    for p in vertex.particles:
        if(isGoldstone(p)):
            return True
    return False

def calculatePrefactor(globalsign,lorentztag,lf,cf) :
    if(globalsign!=1.) :
        prefactors = '(%s) * (%s) * (%s)' \
                     % (globalsign**(len(lorentztag)-2),lf,cf)
    else :
        prefactors = '(%s) * (%s)' \
                     % (lf,cf)
    return prefactors
    # fact=[]
    # if(globalsign!=1.) :
    #     fact.append(globalsign**(len(lorentztag)-2))
    # if(lf!="1") :
    #     fact.append(lf)
    # if(cf!="1") :
    #     fact.append(cf)
    # if(len(fact)==0) : return "1"
    # prefactor = '(%s)' % fact[0]
    # for ix in range(1,len(fact)) :
    #     prefactor = '%s * (%s)' % (prefactor,fact[ix])
    # return prefactor

def couplingValue(coupling) :
    if type(coupling) is not list:
        value = coupling.value
    else:
        value = "("
        for coup in coupling :
            value += '+(%s)' % coup.value
            value +=")"
    return value

def epsilonSign(vertex,couplingptrs,append) :
    EPSSIGN = """\
    double sign = {epssign};
    if((p1->id()=={id1} && p2->id()=={id3} && p3->id()=={id2}) ||
       (p1->id()=={id2} && p2->id()=={id1} && p3->id()=={id3}) ||
       (p1->id()=={id3} && p2->id()=={id2} && p3->id()=={id1})) {{
       sign *= -1.;
    }}
    norm(norm()*sign);
"""
    if(not "p1" in couplingptrs[0]) :
        couplingptrs[0] += ' p1'
    if(not "p2" in couplingptrs[1]) :
        couplingptrs[1] += ' p2'
    if(not "p3" in couplingptrs[2]) :
        couplingptrs[2] += ' p3'
    if("Bar" not in vertex.color[0]) :
        order,sign = extractAntiSymmetricIndices(vertex.color[0],"Epsilon(")
    else :
        order,sign = extractAntiSymmetricIndices(vertex.color[0],"EpsilonBar(")
    subs = {"id1" : vertex.particles[int(order[0])-1].pdg_code,
            "id2" : vertex.particles[int(order[1])-1].pdg_code,
            "id3" : vertex.particles[int(order[2])-1].pdg_code,
            "epssign" : sign }
    append+=EPSSIGN.format(**subs)
    return couplingptrs,append

class VertexConverter:
    'Convert the vertices in a FR model to extract the information ThePEG needs.'
    def __init__(self,model,parmsubs) :
        'Initialize the parameters'
        self.ONE_EACH=True
        self.verbose=False
        self.vertex_skipped=False
        self.ignore_skipped=False
        self.model=model
        self.all_vertices= []
        self.modelname=""
        self.globalsign=self.global_sign()
        self.no_generic_loop_vertices = False
        self.parmsubs = parmsubs

    def global_sign(self):
        'Initial pass to find global sign at the moment does nothing'
        return  1.0
        # for v in self.model.all_vertices:
        #     pids = sorted([ p.pdg_code for p in v.particles ])
        #     if pids != [-11,11,22]: continue
        #     coup = v.couplings
        #     assert( len(coup) == 1 )
        #     val = coup.values()[0].value
        #     val = evaluate(val)
        #     assert( val.real == 0 )
        #     return 1 if val.imag > 0 else -1
        
    def readArgs(self,args) :
        'Extract the relevant command line arguments'
        self.ignore_skipped = args.ignore_skipped
        self.verbose        = args.verbose
        self.modelname = args.name
        self.no_generic_loop_vertices = args.no_generic_loop_vertices

    def should_print(self) :
        'Check if we should output the results'
        return not self.vertex_skipped or self.ignore_skipped

    def convert(self) :
        'Convert the vertices'
        if(self.verbose) :
            print 'verbose mode on: printing all vertices'
            print '-'*60
            labels = ('vertex', 'particles', 'Lorentz', 'C_L', 'C_R', 'norm')
            pprint.pprint(labels)
        # check if we should merge vertices
        if(self.ONE_EACH) :
            self.all_vertices = self.model.all_vertices
        else:
            self.all_vertices = collapse_vertices(self.model.all_vertices)
        # convert the vertices
        vertexclasses, vertexheaders = [], set()
        for vertexnumber,vertex in enumerate(self.all_vertices,1) :
           # process the vertex 
           (skip,vertexClass,vertexHeader) = \
           self.processVertex(vertexnumber,vertex)
           # check it can be handled
           if(skip) : continue
           # add to the list
           vertexclasses.append(vertexClass)
           vertexheaders.add(vertexHeader)
           WRAP = 25
           if vertexnumber % WRAP == 0:
               write_vertex_file({'vertexnumber' : vertexnumber//WRAP,
                                  'vertexclasses' : '\n'.join(vertexclasses),
                                  'vertexheaders' : ''.join(vertexheaders),
                                  'ModelName' : self.modelname})
               vertexclasses = []
               vertexheaders = set()
        # exit if there's vertices we can't handle
        if not self.should_print():
            sys.stderr.write(
"""
Error: The conversion was unsuccessful, some vertices could not be
generated. If you think the missing vertices are not important 
and want to go ahead anyway, use --ignore-skipped. 
Herwig may not give correct results, though.
"""
            )
            sys.exit(1)
        # if still stuff to output it do it
        if vertexclasses:
            write_vertex_file({'vertexnumber' : vertexnumber//WRAP + 1,
                               'vertexclasses' : '\n'.join(vertexclasses),
                               'vertexheaders' : ''.join(vertexheaders),
                               'ModelName' : self.modelname})
        
        print '='*60

    def setCouplingPtrs(self,lorentztag,qcd,append,prepend) :
        couplingptrs = [',tcPDPtr']*len(lorentztag)
        if lorentztag == 'VSS':
            couplingptrs[1] += ' p2'
        elif lorentztag == 'FFV':
            couplingptrs[0] += ' p1'
        elif (lorentztag == 'VVV' or lorentztag == 'VVVS' or
              lorentztag == "SSS" or lorentztag == "VVVT" ) \
             and (append  or prepend ) :
            couplingptrs[0] += ' p1'
            couplingptrs[1] += ' p2'
            couplingptrs[2] += ' p3'
        elif (lorentztag == 'VVVV' and qcd != 2) or\
             (lorentztag == "SSSS" and prepend ):
            couplingptrs[0] += ' p1'
            couplingptrs[1] += ' p2'
            couplingptrs[2] += ' p3'
            couplingptrs[3] += ' p4'
        return couplingptrs

    def processVertex(self,vertexnumber,vertex) :
        # get the Lorentz tag for the vertex
        lorentztag,order = unique_lorentztag(vertex)
        # check if we should skip the vertex
        vertex.herwig_skip_vertex = checkGhostGoldstoneVertex(lorentztag,vertex)
        if(vertex.herwig_skip_vertex) :
            return (True,"","")
        # get the factor for the vertex
        try:
            lf = lfactors[lorentztag]
        except KeyError:
            msg = 'Warning: Lorentz structure {tag} ( {ps} ) in {name} ' \
                  'is not supported.\n'.format(tag=lorentztag, name=vertex.name, 
                                               ps=' '.join(map(str,vertex.particles)))
            sys.stderr.write(msg)
            vertex.herwig_skip_vertex = True
            self.vertex_skipped=True
            return (True,"","")
        # get the ids of the particles at the vertex
        if self.ONE_EACH:
            plistarray = [ ','.join([ str(vertex.particles[o-1].pdg_code) for o in order ]) ]
        else:
            plistarray = [ ','.join([ str(p.pdg_code) for p in pl ])
                           for pl in vertex.particles_list ]
        # parse the colour structure for the vertex
        try:
            L,pos = colors(vertex)
            cf = colorfactor(vertex,L,pos,lorentztag)
        except SkipThisVertex:
            msg = 'Warning: Color structure for vertex ( {ps} ) in {name} ' \
                  'is not supported.\n'.format(tag=lorentztag, name=vertex.name, 
                                               ps=' '.join(map(str,vertex.particles)))
            sys.stderr.write(msg)
            vertex.herwig_skip_vertex = True
            self.vertex_skipped=True
            return (True,"","")
    
        ### classname
        classname = 'V_%s' % vertex.name
        # try to extract the couplings
        try:
            (all_couplings,header,qcd,qed,prepend,append) = \
            self.extractCouplings(lorentztag,pos,lf,cf,vertex,order)
        except SkipThisVertex:
            msg = 'Warning: Lorentz structure {tag} ( {ps} ) in {name} ' \
                  'is not supported, may have a non-perturbative form.\n'.format(tag=lorentztag, name=vertex.name, 
                                               ps=' '.join(map(str,vertex.particles)))
            sys.stderr.write(msg)
            vertex.herwig_skip_vertex = True
            self.vertex_skipped=True
            return (True,"","")

        # set the coupling ptrs in the setCoupling call
        couplingptrs = self.setCouplingPtrs(lorentztag,qcd,append != '',prepend != '')

        # final processing of the couplings
        try :
            symbols = set()
            if(lorentztag in ['FFS','FFV']) :
                (normcontent,leftcontent,rightcontent,append) = processFermionCouplings(lorentztag,vertex,
                                                                                        self.model,self.parmsubs,
                                                                                        all_couplings)
            elif('T' in lorentztag) :
                (leftcontent,rightcontent,normcontent) = processTensorCouplings(lorentztag,vertex,self.model,
                                                                                self.parmsubs,all_couplings)
            elif(lorentztag=="SSS" or lorentztag=="SSSS") :
                normcontent = processScalarCouplings(self.model,self.parmsubs,all_couplings)
            elif(lorentztag=="VVS" or lorentztag =="VVSS" or lorentztag=="VSS") :
                normcontent,append,lorentztag,header,sym = processScalarVectorCouplings(lorentztag,vertex,
                                                                                        self.model,self.parmsubs,
                                                                                        all_couplings,header,order)
                symbols |=sym
            elif("VVV" in lorentztag) :
                normcontent,append,header =\
                                            processVectorCouplings(lorentztag,vertex,self.model,self.parmsubs,all_couplings,append,header)
            else :
                SkipThisVertex()
        except SkipThisVertex:
            msg = 'Warning: Lorentz structure {tag} ( {ps} ) in {name} ' \
                  'is not supported, may have a non-perturbative form.\n'.format(tag=lorentztag, name=vertex.name,
                                               ps=' '.join(map(str,vertex.particles)))
            sys.stderr.write(msg)
            vertex.herwig_skip_vertex = True
            self.vertex_skipped=True
            return (True,"","")

        
        ### do we need left/right?
        if 'FF' in lorentztag and lorentztag != "FFT":
            #leftcalc = aStoStrongCoup(py2cpp(leftcontent)[0], paramstoreplace_, paramstoreplace_expressions_)
            #rightcalc = aStoStrongCoup(py2cpp(rightcontent)[0], paramstoreplace_, paramstoreplace_expressions_)
            leftcalc, sym = py2cpp(leftcontent)
            symbols |= sym
            rightcalc, sym = py2cpp(rightcontent)
            symbols |= sym
            left = 'left(' + leftcalc + ');'
            right = 'right(' + rightcalc + ');'
        else:
            left = ''
            right = ''
            leftcalc = ''
            rightcalc = ''
            #normcalc = aStoStrongCoup(py2cpp(normcontent)[0], paramstoreplace_, paramstoreplace_expressions_)
        normcalc, sym = py2cpp(normcontent)
        symbols |= sym
        # UFO is GeV by default (?)
        if lorentztag in ['VVS','SSS']:
            normcalc = 'Complex((%s) * GeV / UnitRemoval::E)' % normcalc
        elif lorentztag in ['GeneralVVS']:
            normcalc = 'Complex(-(%s) * UnitRemoval::E / GeV )' % normcalc
        elif lorentztag in ['FFT','VVT', 'SST', 'FFVT', 'VVVT' , 'VVVS' ]:
            normcalc = 'Complex((%s) / GeV * UnitRemoval::E)' % normcalc
        norm = 'norm(' + normcalc + ');'
        # finally special handling for eps tensors
        if(len(vertex.color)==1 and vertex.color[0].find("Epsilon")>=0) :
            couplingptrs, append = epsilonSign(vertex,couplingptrs,append)
        # define unkown symbols from the model
        symboldefs = [ def_from_model(self.model,s) for s in symbols ]
        ### assemble dictionary and fill template
        subs = { 'lorentztag' : lorentztag,                   # ok
                 'classname'  : classname,            # ok
                 'symbolrefs' : '\n    '.join(symboldefs),
                 'left'       : left,                 # doesn't always exist in base
                 'right'      : right,                 # doesn't always exist in base 
                 'norm'      : norm,                 # needs norm, too

                 #################### need survey which different setter methods exist in base classes

                 'addToPlist' : '\n'.join([ 'addToList(%s);'%s for s in plistarray]),
                 'parameters' : '',
                 'setCouplings' : '',
                 'qedorder'   : qed,
                 'qcdorder' : qcd,
                 'couplingptrs' : ''.join(couplingptrs),
                 'spindirectory' : spindirectory(lorentztag),
                 'ModelName' : self.modelname,
                 'prepend' : prepend,
                 'append' : append,
                 'header' : header
                 }             # ok

        # print info if required
        if self.verbose:
            print '-'*60
            pprint.pprint(( classname, plistarray, leftcalc, rightcalc, normcalc ))

        return (False,VERTEXCLASS.substitute(subs),VERTEXHEADER.format(**subs))

    def get_vertices(self,libname):
        vlist = ['library %s\n' % libname]
        for v in self.all_vertices:
            if v.herwig_skip_vertex: continue
            vlist.append( vertexline.format(modelname=self.modelname, vname=v.name) )
        if( not self.no_generic_loop_vertices) :
            vlist.append('insert {modelname}:ExtraVertices 0 /Herwig/{modelname}/V_GenericHPP\n'.format(modelname=self.modelname) )
            vlist.append('insert {modelname}:ExtraVertices 0 /Herwig/{modelname}/V_GenericHGG\n'.format(modelname=self.modelname) )
        return ''.join(vlist)


    def extractCouplings(self,lorentztag,pos,lf,cf,vertex,order) :
        coup_left  = []
        coup_right = []
        coup_norm = []
        header = ""
        qcd=0
        qed=0
        prepend=""
        append=""
        unique_qcd = CheckUnique()
        unique_qed = CheckUnique()
        maxColour=0
        for (color_idx,lorentz_idx),coupling in vertex.couplings.iteritems():
            maxColour=max(maxColour,color_idx)            
        all_couplings=[]
        for ix in range(0,maxColour+1) :
            all_couplings.append([])
        for colour in range(0,maxColour+1) :
            for (color_idx,lorentz_idx),coupling in vertex.couplings.iteritems() :
                if(color_idx!=colour) : continue
                qcd, qed = qcd_qed_orders(vertex, coupling)
                try :
                    unique_qcd( qcd )
                    unique_qed( qed )
                except :
                    msg = 'Different powers of QCD and QED couplings for the same vertex'\
                          ' is not currently supported for {ps} in {name}.\n'.format(tag=lorentztag, name=vertex.name,
                                                                                   ps=' '.join(map(str,vertex.particles)))
                    sys.stderr.write(msg)
                    raise SkipThisVertex()
                L = vertex.lorentz[lorentz_idx]
                prefactors = calculatePrefactor(self.globalsign,lorentztag,lf,cf[color_idx])
                # calculate the value of the coupling
                value = couplingValue(coupling)
                # handling of the different types of couplings
                if lorentztag in ['FFS','FFV']:
                    all_couplings[color_idx] = fermionCouplings(value,prefactors,L,all_couplings[color_idx],order)
                elif 'T' in lorentztag :
                    append, all_couplings[color_idx] = tensorCouplings(vertex,value,prefactors,L,lorentztag,pos,
                                                                       all_couplings[color_idx],order)
                elif 'R' in lorentztag :
                    all_couplings[color_idx] = RSCouplings(value,prefactors,L,all_couplings[color_idx],order)
                elif lorentztag == 'VVS' or lorentztag == "VVSS" or lorentztag == "VSS" :
                    all_couplings[color_idx] = scalarVectorCouplings(value,prefactors,L,lorentztag,
                                                                     all_couplings[color_idx],order)
                elif lorentztag == "SSS" or lorentztag == "SSSS" :
                    prepend, header, all_couplings[color_idx] = scalarCouplings(vertex,value,prefactors,L,lorentztag,
                                                                                all_couplings[color_idx],prepend,header)
                elif "VVV" in lorentztag :
                    all_couplings[color_idx],append = vectorCouplings(vertex,value,prefactors,L,lorentztag,pos,
                                                                      all_couplings[color_idx],append,qcd,order)
                else:
                    raise SkipThisVertex()

        # return the result
        return (all_couplings,header,qcd,qed,prepend,append)
