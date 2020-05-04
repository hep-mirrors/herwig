from __future__ import print_function
import sys,pprint

from .helpers import CheckUnique,getTemplate,writeFile,coupling_orders,def_from_model
from .converter import py2cpp
from .collapse_vertices import collapse_vertices
from .check_lorentz import tensorCouplings,VVVordering,lorentzScalar,\
    processTensorCouplings,scalarCouplings,processScalarCouplings,scalarVectorCouplings,\
    processScalarVectorCouplings,vectorCouplings,processVectorCouplings,fermionCouplings,processFermionCouplings,\
    RSCouplings
from .general_lorentz import convertLorentz,generateEvaluateFunction,multipleEvaluate
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

genericVertices=['FFFF','FFVV','FFSS','FFVS','VVVV','VVVT','FFVT',
                 'RFVV','RFVS','RFSS','SSST','VVST','FFST']

skipped5Point=False

# template for the header for a vertex
VERTEXHEADER = """\
#include "ThePEG/Helicity/Vertex/{spindirectory}/{lorentztag}Vertex.h"
"""
GENERALVERTEXHEADER = """\
#include "ThePEG/Helicity/Vertex/Abstract{lorentztag}Vertex.h"
"""

# template for the implmentation for a vertex
VERTEXCLASS = getTemplate('Vertex_class')
GENERALVERTEXCLASS = getTemplate('GeneralVertex_class')

# template for the .cc file for vertices
VERTEX = getTemplate('Vertex.cc')

vertexline = """\
create Herwig::FRModel{classname} /Herwig/{modelname}/{classname}
insert {modelname}:ExtraVertices 0 /Herwig/{modelname}/{classname}
"""

def get_lorentztag(spin):
    """Produce a ThePEG spin tag for the given numeric FR spins."""
    spins = { 1 : 'S', 2 : 'F', 3 : 'V', 4 : 'R', 5 : 'T', -1 : 'U' }
    result=[]
    for i in range(0,len(spin)) :
        result.append((spins[spin[i]],i+1))
    def spinsort(inVal):
        output=0
        vals=['U','R','F','V','S','T']
        for val in vals:
            if inVal==val : break
            output+=1
        return output

    result = sorted(result, key=spinsort)
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

def coloursort(a) :
    return int(a[4])
    
def colorfactor(vertex,L,pos,lorentztag):
    def match(patterns,color=vertex.color):
        result = [ p == t
                   for p,t in zip(patterns,color) ]
        return all(result)

    label = None
    l = lambda c: len(pos[c])
    if l(1) == L:
        label = ('1',)
        if match(label): return ("SINGLET",('1.',))

    elif l(3) == l(-3) == 1 and l(1) == L-2:
        nums = [pos[3][0], pos[-3][0]]
        label = ('Identity({0},{1})'.format(*sorted(nums)),)
        if match(label): return ("DELTA",('1.',))

    elif l(6) == l(-6) == 1 and l(1) == L-2:
        nums = [pos[6][0], pos[-6][0]]
        label = ('Identity({0},{1})'.format(*sorted(nums)),)
        if match(label): return ("DELTA",('1.',))

    elif l(6) == l(-6) == 2 and L==4:
        sys.stderr.write(
            'Warning: Unknown colour structure 6 6 6~ 6~ ( {ps} ) in {name}.\n'
            .format(name=vertex.name, ps=' '.join(map(str,vertex.particles)))
        )
        raise SkipThisVertex()

    elif l(8) == l(3) == l(-3) == 1 and l(1) == L-3:
        label = ('T({g},{q},{qb})'.format(g=pos[8][0],q=pos[3][0],qb=pos[-3][0]),)
        if match(label): return ("SU3TFUND",('1.',))

    elif l(8) == l(6) == l(-6) == 1 and l(1) == L-3:
        label = ('T6({g},{s},{sb})'.format(g=pos[8][0],s=pos[6][0],sb=pos[-6][0]),)
        if match(label): return ("SU3T6",('1.',))

    elif l(6) == 1 and l(-3) == 2 and L==3:
        label = ('K6({s},{qb1},{qb2})'.format(s=pos[6][0],qb1=pos[-3][0],qb2=pos[-3][1]),)
        if match(label): return ("SU3K6",('1.',))

    elif l(-6) == 1 and l(3) == 2 and L==3:
        label = ('K6Bar({sb},{q1},{q2})'.format(sb=pos[-6][0],q1=pos[3][0],q2=pos[3][1]),)
        if match(label): return ("SU3K6",('1.',))

    elif l(3) == L == 3:
        colors=[]
        for color in vertex.color :
            order,sign  = extractAntiSymmetricIndices(color,"Epsilon(")
            colors.append("Epsilon(%s,%s,%s)" % (order[0],order[1],order[2]))
        label = ('Epsilon(1,2,3)',)
        if match(label,colors): return ("EPS",('1.',)) # TODO check factor!

    elif l(-3) == L == 3:
        colors=[]
        for color in vertex.color :
            order,sign  = extractAntiSymmetricIndices(color,"EpsilonBar(")
            colors.append("Epsilon(%s,%s,%s)" % (order[0],order[1],order[2]))
        label = ('EpsilonBar(1,2,3)',)
        if match(label): return ("EPS",('1.',)) # TODO check factor!

    elif l(8) == L == 3:
        colors=[]
        for color in vertex.color :
            order,sign  = extractAntiSymmetricIndices(color,"f(")
            colors.append("f(%s,%s,%s)" % (order[0],order[1],order[2]))
        # if lorentz is FFV get extra minus sign
        if lorentztag in ['FFV'] : sign *=-1
        label = ('f(1,2,3)',)
        if match(label,colors): return ("SU3F",('-complex(0,1.)*(%s)'%sign,))

    elif l(8) == 3 and l(1)==1 and L == 4:
        colors=[]
        for color in vertex.color :
            order,sign  = extractAntiSymmetricIndices(color,"f(")
            colors.append("f(%s,%s,%s)" % (order[0],order[1],order[2]))
        if(pos[1][0]==1) :
            label = ('f(2,3,4)',)
        elif(pos[1][0]==2) :
            label = ('f(1,3,4)',)
        elif(pos[1][0]==3) :
            label = ('f(1,2,4)',)
        elif(pos[1][0]==4) :
            label = ('f(1,2,3)',)
        if match(label,colors): return ("SU3F",('-complex(0,1.)*(%s)'%sign,))
    elif l(8) == L == 4:
        colors=[]
        for color in vertex.color :
            f = color.split("*")
            (o1,s1) = extractAntiSymmetricIndices(f[0],"f(")
            (o2,s2) = extractAntiSymmetricIndices(f[1],"f(")
            if(o2[0]<o1[0]) : o1,o2=o2,o1
            colors.append("f(%s)*f(%s)" % (",".join(o1),",".join(o2)))
        colors=sorted(colors,key=coloursort)
        label = ('f(1,2,-1)*f(3,4,-1)',
                 'f(1,3,-1)*f(2,4,-1)',
                 'f(1,4,-1)*f(2,3,-1)')
        nmatch=0
        for c1 in label:
            for  c2 in colors :
                if(c1==c2) : nmatch+=1
        if(nmatch==2 and lorentztag=="VVSS") :
            return ("SU3FF",('1.','1.'))
        elif(nmatch==3 and lorentztag=="VVVV") :
            return ("SU3FF",('1.','1.','1.'))

    elif l(8) == 2 and l(3) == l(-3) == 1 and L==4:
        subs = {
            'g1' : pos[8][0],
            'g2' : pos[8][1],
            'qq' : pos[3][0],
            'qb' : pos[-3][0] 
        }
        if(vertex.lorentz[0].spins.count(1)==2) :
            label = ('T({g1},-1,{qb})*T({g2},{qq},-1)'.format(**subs),
                     'T({g1},{qq},-1)*T({g2},-1,{qb})'.format(**subs))
            if match(label): return ("SU3TTFUNDS",('1.','1.'))
        elif(vertex.lorentz[0].spins.count(2)==2) :
            label = ('f({g1},{g2},-1)*T(-1,{qq},{qb})'.format(**subs),)
            if match(label): return ("SU3TTFUNDD",('-complex(0.,1.)',))
            label = ('f(-1,{g1},{g2})*T(-1,{qq},{qb})'.format(**subs),)
            if match(label): return ("SU3TTFUNDD",('-complex(0.,1.)',))
        elif(vertex.lorentz[0].spins.count(3)==4) :
            label = ('f(-1,{g1},{g2})*T(-1,{qq},{qb})'.format(**subs),
                     'T({g1},-1,{qb})*T({g2},{qq},-1)'.format(**subs),
                     'T({g1},{qq},-1)*T({g2},-1,{qb})'.format(**subs))
            if match(label): return ("SU3TTFUNDS",(('-complex(0.,1.)','complex(0.,1.)'),'1.','1.'))
            
        
    elif l(8) == 2 and l(6) == l(-6) == 1 and L==4:
        subs = {
            'g1' : pos[8][0],
            'g2' : pos[8][1],
            'qq' : pos[6][0],
            'qb' : pos[-6][0] 
        }
        label = ('T6({g1},-1,{qb})*T6({g2},{qq},-1)'.format(**subs),
                 'T6({g1},{qq},-1)*T6({g2},-1,{qb})'.format(**subs))
        if match(label): return ("SU3TT6",('1.','1.'))

    elif l(8) == 2 and l(8)+l(1)==L :
        subs = { 'g1' : pos[8][0], 'g2' : pos[8][1] }
        label = ('Identity({g1},{g2})'.format(**subs),)
        if match(label) : return ("DELTA",('1.',))

    elif l(8) == 3 and l(1)==1 and L==4 :
        colors=[]
        for color in vertex.color :
            order,sign  = extractAntiSymmetricIndices(color,"f(")
            colors.append("f(%s,%s,%s)" % (order[0],order[1],order[2]))
        label = ('f(1,2,3)',)
        if match(label,colors): return ("SU3F",('-complex(0.,1.)*(%s)'%sign,))

    elif l(3)==2 and l(-3) == 2 and L==4  and lorentztag=="FFFF" :
        labels=["Identity(1,2)*Identity(3,4)",
                "Identity(1,4)*Identity(2,3)",
                "T(-1,2,1)*T(-1,4,3)",
                "T(-1,2,3)*T(-1,4,1)"]
        cstruct=["SU3I12I34","SU3I14I23","SU3T21T43","SU3T23T41"]
        oname=[]
        ovalue=[]
        for color in vertex.color :
            for i in range(0,len(labels)) :
                if labels[i]==color : break
            if(i<len(labels)) :
                oname.append(cstruct[i])
                ovalue.append("1.")
            else :
                sys.stderr.write(
                    "Warning: Unknown colour structure {color} ( {ps} ) in {name} for FFFF vertex.\n"
                    .format(color = ' '.join(vertex.color), name = vertex.name,
                            ps = ' '.join(map(str,vertex.particles)))
                )
                raise SkipThisVertex()
        return(oname,ovalue)
        
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

def calculatePrefactor(lf,cf) :
    prefactor = '(%s) * (%s)' % (lf,cf)
    return prefactor

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
    def __init__(self,model,parmsubs,defns) :
        'Initialize the parameters'
        self.verbose=False
        self.vertex_skipped=False
        self.ignore_skipped=False
        self.model=model
        self.all_vertices= []
        self.vertex_names = {}
        self.modelname=""
        self.no_generic_loop_vertices = False
        self.parmsubs = parmsubs
        self.couplingDefns = defns
        self.genericTensors = False
        self.hw_higgs = False
        
    def readArgs(self,args) :
        'Extract the relevant command line arguments'
        self.ignore_skipped = args.ignore_skipped
        self.verbose        = args.verbose
        self.modelname = args.name
        self.no_generic_loop_vertices = args.no_generic_loop_vertices
        self.include_generic = args.include_generic
        self.genericTensors = args.use_generic_for_tensors
        self.hw_higgs = args.use_Herwig_Higgs
        
    def should_print(self) :
        'Check if we should output the results'
        return not self.vertex_skipped or self.ignore_skipped

    def convert(self) :
        'Convert the vertices'
        if(self.verbose) :
            print('verbose mode on: printing all vertices')
            print('-'*60)
            labels = ('vertex', 'particles', 'Lorentz', 'C_L', 'C_R', 'norm')
            pprint.pprint(labels)
        # extract the vertices
        self.all_vertices = self.model.all_vertices
        # find the SM Higgs boson
        higgs=None
        if(self.hw_higgs) :
            for part in self.model.all_particles:
                if(part.pdg_code==25) :
                    higgs=part
                    break
        # convert the vertices
        vertexclasses, vertexheaders = [], set()
        ifile=1
        icount=0
        for vertexnumber,vertex in enumerate(self.all_vertices,1) :
           # process the vertex 
           (skip,vertexClass,vertexHeader) = \
           self.processVertex(vertexnumber,vertex)
           # check it can be handled
           if(skip) : continue
           # check if Higgs and skip if using Hw higgs sector
           if higgs in vertex.particles :
               nH = vertex.particles.count(higgs) 
               # skip trilinear and quartic higgs vertices
               if ( nH == len(vertex.particles) ) :
                   vertex.herwig_skip_vertex = True
                   continue
               elif (len(vertex.particles)-nH==2) :
                   p=[]
                   for part in vertex.particles :
                       if(part!=higgs) : p.append(part)
                   if(nH==1 and p[0].pdg_code==-p[1].pdg_code and
                      abs(p[0].pdg_code) in [1,2,3,4,5,6,11,13,15,24]) :
                       vertex.herwig_skip_vertex = True
                       continue
                   elif((nH==1 or nH==2) and p[0].pdg_code==p[1].pdg_code and
                        p[0].pdg_code ==23) :
                       vertex.herwig_skip_vertex = True
                       continue
                   elif(nH==2 and  p[0].pdg_code==-p[1].pdg_code and
                        abs(p[0].pdg_code)==24) :
                       vertex.herwig_skip_vertex = True
                       continue
           # add to the list
           icount +=1
           vertexclasses.append(vertexClass)
           vertexheaders.add(vertexHeader)
           WRAP = 25
           if icount % WRAP == 0 or vertexHeader.find("Abstract")>=0:
               write_vertex_file({'vertexnumber' : ifile,
                                  'vertexclasses' : '\n'.join(vertexclasses),
                                  'vertexheaders' : ''.join(vertexheaders),
                                  'ModelName' : self.modelname})
               vertexclasses = []
               vertexheaders = set()
               icount=0
               ifile+=1
        # exit if there's vertices we can't handle
        if not self.should_print():
            sys.stderr.write(
"""
Error: The conversion was unsuccessful, some vertices could not be
generated. The new --include-generic option should be able to generate
these. Otherwise, if you think the missing vertices are not important 
and want to go ahead anyway, use --ignore-skipped. 
Herwig may not give correct results, though.
"""
            )
            sys.exit(1)
        # if still stuff to output it do it
        if vertexclasses:
            write_vertex_file({'vertexnumber' : ifile + 1,
                               'vertexclasses' : '\n'.join(vertexclasses),
                               'vertexheaders' : ''.join(vertexheaders),
                               'ModelName' : self.modelname})
        
        print('='*60)

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
        elif (lorentztag == 'VVVV' and (qcd < 2 or append)) or\
             (lorentztag == "SSSS" and prepend ):
            couplingptrs[0] += ' p1'
            couplingptrs[1] += ' p2'
            couplingptrs[2] += ' p3'
            couplingptrs[3] += ' p4'
        return couplingptrs

    def processVertex(self,vertexnumber,vertex) :
        global skipped5Point
        # get the Lorentz tag for the vertex
        lorentztag,order = unique_lorentztag(vertex)
        # check if we should skip the vertex
        vertex.herwig_skip_vertex = checkGhostGoldstoneVertex(lorentztag,vertex)
        # check the order of the vertex and skip 5 points
        if(len(lorentztag)>=5) :
            vertex.herwig_skip_vertex = True
            if(not skipped5Point) :
                skipped5Point = True
                print("Skipping 5 point vertices which aren\'t used in Herwig7")
                
        if(vertex.herwig_skip_vertex) :
            return (True,"","")
        # check if we support this at all
        if( lorentztag not in lfactors and
            lorentztag not in genericVertices) :
            msg = 'Warning: Lorentz structure {tag} ( {ps} ) in {name} ' \
                  'is not supported.\n'.format(tag=lorentztag, name=vertex.name, 
                                               ps=' '.join(map(str,vertex.particles)))
            sys.stderr.write(msg)
            vertex.herwig_skip_vertex = True
            self.vertex_skipped=True
            return (True,"","")
        # get the factor for the vertex
        generic = False
        try:
            lf = lfactors[lorentztag]
            if( self.genericTensors and "T" in lorentztag) :
                raise KeyError()
        except KeyError:
            if(not self.include_generic) :
                msg = 'Warning: Lorentz structure {tag} ( {ps} ) in {name} ' \
                      'is not supported.\n'.format(tag=lorentztag, name=vertex.name, 
                                                   ps=' '.join(map(str,vertex.particles)))
                sys.stderr.write(msg)
                vertex.herwig_skip_vertex = True
                self.vertex_skipped=True
                return (True,"","")
            else :
                lf=1.
                generic=True
        # get the ids of the particles at the vertex
        plistarray = [ ','.join([ str(vertex.particles[o-1].pdg_code) for o in order ]) ]
        # parse the colour structure for the vertex
        try:
            L,pos = colors(vertex)
            cs,cf = colorfactor(vertex,L,pos,lorentztag)
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
        if(not generic) :
            try :
                return self.extractGeneric(vertex,order,lorentztag,classname,plistarray,pos,lf,cf,cs)
            except SkipThisVertex:
                if(not self.include_generic) :
                    msg = 'Warning: Lorentz structure {tag} ( {ps} ) in {name} ' \
                          'is not supported, may have a non-perturbative form.\n'.format(tag=lorentztag, name=vertex.name, 
                                                                                         ps=' '.join(map(str,vertex.particles)))
                    
                    sys.stderr.write(msg)
                    vertex.herwig_skip_vertex = True
                    self.vertex_skipped=True
                    return (True,"","")
                else :
                    try :
                        return self.extractGeneral(vertex,order,lorentztag,classname,plistarray,pos,cf,cs)
                    except SkipThisVertex:
                        msg = 'Warning: Lorentz structure {tag} ( {ps} ) in {name} ' \
                              'is not supported, may have a non-perturbative form.\n'.format(tag=lorentztag, name=vertex.name, 
                                                                                             ps=' '.join(map(str,vertex.particles)))
                    
                        sys.stderr.write(msg)
                        vertex.herwig_skip_vertex = True
                        self.vertex_skipped=True
                        return (True,"","")
        else :
            try :
                return self.extractGeneral(vertex,order,lorentztag,classname,plistarray,pos,cf,cs)
            except SkipThisVertex:
                msg = 'Warning: Lorentz structure {tag} ( {ps} ) in {name} ' \
                      'is not supported, may have a non-perturbative form.\n'.format(tag=lorentztag, name=vertex.name, 
                                                                                     ps=' '.join(map(str,vertex.particles)))
                
                sys.stderr.write(msg)
                vertex.herwig_skip_vertex = True
                self.vertex_skipped=True
                return (True,"","")
            
            
    def extractGeneric(self,vertex,order,lorentztag,classname,plistarray,pos,lf,cf,cs) :
        classes=""
        headers=""
        # identify the maximum colour flow and orders of the couplings
        maxColour=0
        couplingOrders=[]
        self.vertex_names[vertex.name] = [classname]
        for (color_idx,lorentz_idx),coupling in vertex.couplings.items():
            maxColour=max(maxColour,color_idx)
            orders = coupling_orders(vertex, coupling, self.couplingDefns)
            if(orders not in couplingOrders) : couplingOrders.append(orders)
        # loop the order of the couplings
        iorder = 0
        for corder in couplingOrders :
            iorder +=1
            cname=classname
            if(iorder!=1) :
                cname= "%s_%s" % (classname,iorder)
                self.vertex_names[vertex.name].append(cname)
            header = ""
            prepend=""
            append=""
            all_couplings=[]
            for ix in range(0,maxColour+1) :
                all_couplings.append([])
            # loop over the colour structures
            for colour in range(0,maxColour+1) :
                for (color_idx,lorentz_idx),coupling in vertex.couplings.items() :
                    # check colour structure and coupling order
                    if(color_idx!=colour) : continue
                    if(coupling_orders(vertex, coupling, self.couplingDefns)!=corder) : continue
                    # get the prefactor for the lorentz structure
                    L = vertex.lorentz[lorentz_idx]
                    prefactors = calculatePrefactor(lf,cf[color_idx])
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
                                                                          all_couplings[color_idx],append,corder["QCD"],order)
                    else:
                        raise SkipThisVertex()
            # final processing of the couplings
            symbols = set()
            if(lorentztag in ['FFS','FFV']) :
                (normcontent,leftcontent,rightcontent,append) = processFermionCouplings(lorentztag,vertex,
                                                                                        self.model,self.parmsubs,
                                                                                        all_couplings,order)
            elif('T' in lorentztag) :
                (leftcontent,rightcontent,normcontent) = processTensorCouplings(lorentztag,vertex,self.model,
                                                                                self.parmsubs,all_couplings,order)
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
            # set the coupling ptrs in the setCoupling call
            couplingptrs = self.setCouplingPtrs(lorentztag.replace("General",""),
                                                                   corder["QCD"],append != '',prepend != '')
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
            # UFO is GeV by default
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
            couplingOrder=""
            for coupName,coupVal in corder.items() :
                couplingOrder+="    orderInCoupling(CouplingType::%s,%s);\n" %(coupName,coupVal)
            ### assemble dictionary and fill template
            subs = { 'lorentztag' : lorentztag,                   # ok
                     'classname'  : cname,            # ok
                     'symbolrefs' : '\n    '.join(symboldefs),
                     'left'       : left,                 # doesn't always exist in base
                     'right'      : right,                 # doesn't always exist in base 
                     'norm'      : norm,                 # needs norm, too
                     'addToPlist' : '\n'.join([ 'addToList(%s);'%s for s in plistarray]),
                     'parameters' : '',
                     'couplingOrders' : couplingOrder,
                     'colourStructure' : cs,
                     'couplingptrs' : ''.join(couplingptrs),
                     'spindirectory' : spindirectory(lorentztag),
                     'ModelName' : self.modelname,
                     'prepend' : prepend,
                     'append' : append,
                     'header' : header
                   }             # ok
        
            # print info if required
            if self.verbose:
                print('-'*60)
                pprint.pprint(( classname, plistarray, leftcalc, rightcalc, normcalc ))
            headers+=VERTEXHEADER.format(**subs)
            classes+=VERTEXCLASS.substitute(subs)
        return (False,classes,headers)

    def extractGeneral(self,vertex,order,lorentztag,classname,plistarray,pos,cf,cs) :
        eps=False
        classes=""
        headers=""
        # check the colour flows, three cases supported either 1 flow or 3 in gggg
        # or multiple wierd ones in FFFF
        cidx=-1
        gluon4point = (len(pos[8])==4 and vertex.lorentz[0].spins.count(3)==4)
        FFFF        = (len(pos[3])==2 and len(pos[-3])==2 and vertex.lorentz[0].spins.count(2)==4)
        couplingOrders=[]
        colours={}
        
        for (color_idx,lorentz_idx),coupling in vertex.couplings.items() :
            orders = coupling_orders(vertex, coupling, self.couplingDefns)
            if(orders not in couplingOrders) : couplingOrders.append(orders)
            if(gluon4point) :
                color =  vertex.color[color_idx]
                f = color.split("*")
                (o1,s1) = extractAntiSymmetricIndices(f[0],"f(")
                (o2,s2) = extractAntiSymmetricIndices(f[1],"f(")
                if(o2[0]<o1[0]) : o1,o2=o2,o1
                color = "f(%s)*f(%s)" % (",".join(o1),",".join(o2))
                label = 'f(1,3,-1)*f(2,4,-1)'
                if(label==color) :
                    cidx=color_idx
                    colours[cidx] = (cs,cf[cidx])
            elif (FFFF) :
                colours[color_idx] = (cs[color_idx],cf[color_idx])
            else :
                cidx=color_idx
                if(color_idx!=0) :
                    vertex.herwig_skip_vertex = True
                    self.vertex_skipped=True
                    msg = 'Warning: General spin structure code currently only '\
                          'supports 1 colour structure for  {tag} ( {ps} ) in {name}\n'.format(tag=lorentztag, name=vertex.name,
                                                                                               ps=' '.join(map(str,vertex.particles)))
                    sys.stderr.write(msg)
                    return (True,"","")
                if(isinstance(cs,basestring)) :
                    colours[cidx] = (cs,cf[cidx])
                else :
                    vertex.herwig_skip_vertex = True
                    self.vertex_skipped=True
                    msg = 'Warning: General spin structure code currently only '\
                          'supports 1 colour structure for  {tag} ( {ps} ) in {name}\n'.format(tag=lorentztag, name=vertex.name,
                                                                                               ps=' '.join(map(str,vertex.particles)))
                    sys.stderr.write(msg)
                    return (True,"","")
        if(len(colours)==0) :
            msg = 'Warning: General spin structure code currently only '\
                  'supports 1 colour structure for  {tag} ( {ps} ) in {name}\n'.format(tag=lorentztag, name=vertex.name,
                                                                                       ps=' '.join(map(str,vertex.particles)))
            sys.stderr.write(msg)
            vertex.herwig_skip_vertex = True
            self.vertex_skipped=True
            return (True,"","")
        # loop over the different orders in the couplings
        # and colour structures
        iorder=0
        self.vertex_names[vertex.name]=[classname]
        for corder in couplingOrders :
            for (cidx,(cstruct,cfactor)) in colours.items() :
                iorder +=1
                cname=classname
                if(iorder!=1) :
                    cname= "%s_%s" % (classname,iorder)
                    self.vertex_names[vertex.name].append(cname)
                defns=[]
                vertexEval=[]
                values=[]
                imax = len(vertex.particles)+1
                if lorentztag in genericVertices :
                    imax=1
                for (color_idx,lorentz_idx),coupling in vertex.couplings.items() :
                    # only the colour structre and coupling order we want
                    if(color_idx != cidx) : continue
                    if(coupling_orders(vertex, coupling, self.couplingDefns)!=corder) : continue
                    # calculate the value of the coupling
                    values.append(couplingValue(coupling))
                    # now to convert the spin structures
                    for i in range(0,imax) :
                        if(len(defns)<i+1) :
                            defns.append({})
                            vertexEval.append([])
                        eps |= convertLorentz(vertex.lorentz[lorentz_idx],lorentztag,order,vertex,
                                              i,defns[i],vertexEval[i])
                # we can now generate the evaluate member functions
                header=""
                impls=""
                spins=vertex.lorentz[0].spins
                mult={}
                for i in range(1,6) :
                    if( spins.count(i)>1 and i!=2) : mult[i] = []
                for i in range(0,imax) :
                    (evalHeader,evalCC) = generateEvaluateFunction(self.model,vertex,i,values,defns[i],vertexEval[i],cfactor,order)
                    if(i!=0 and spins[i-1] in mult) :
                        if(len(mult[spins[i-1]])==0) : mult[spins[i-1]].append(evalHeader)
                        evalHeader=evalHeader.replace("evaluate(","evaluate%s(" % i)
                        evalCC    =evalCC    .replace("evaluate(","evaluate%s(" % i)
                        mult[spins[i-1]].append(evalHeader)
                    header+="    "+evalHeader+";\n"
                    impls+=evalCC
                # combine the multiple defn if needed
                for (key,val) in mult.items() :
                    (evalHeader,evalCC) = multipleEvaluate(vertex,key,val)
                    if(evalHeader!="") : header += "    "+evalHeader+";\n"
                    if(evalCC!="")     : impls   += evalCC
                impls=impls.replace("evaluate", "FRModel%s::evaluate" % cname)
                couplingOrder=""
                for coupName,coupVal in corder.items() :
                    couplingOrder+="    orderInCoupling(CouplingType::%s,%s);\n" %(coupName,coupVal)
                ### assemble dictionary and fill template
                subs = { 'lorentztag' : lorentztag,
                         'classname'  : cname,
                         'addToPlist' : '\n'.join([ 'addToList(%s);'%s for s in plistarray]),
                         'ModelName' : self.modelname,
                        'couplingOrders' : couplingOrder,
                         'colourStructure' : cstruct,
                         'evaldefs'  : header,
                         'evalimpls' : impls}
                newHeader = GENERALVERTEXHEADER.format(**subs)
                if(eps) : newHeader +="#include \"ThePEG/Helicity/epsilon.h\"\n"
                headers+=newHeader
                classes+=GENERALVERTEXCLASS.substitute(subs)
        return (False,classes,headers)

    def get_vertices(self,libname):
        vlist = ['library %s\n' % libname]
        for v in self.all_vertices:
            if v.herwig_skip_vertex: continue
            for name in self.vertex_names[v.name] :
                vlist.append( vertexline.format(modelname=self.modelname, classname=name) )
        
        if( not self.no_generic_loop_vertices and not self.hw_higgs ) :
            vlist.append('insert {modelname}:ExtraVertices 0 /Herwig/{modelname}/V_GenericHPP\n'.format(modelname=self.modelname) )
            vlist.append('insert {modelname}:ExtraVertices 0 /Herwig/{modelname}/V_GenericHGG\n'.format(modelname=self.modelname) )
        # add Hw higgs vertices if required
        if(self.hw_higgs) :
            for vertex in ["FFH","HGG","HHH","HPP","HZP","WWHH","WWH"]:
                vlist.append('insert %s:ExtraVertices 0 /Herwig/Vertices/%sVertex\n' % (self.modelname,vertex) )
        return ''.join(vlist)
