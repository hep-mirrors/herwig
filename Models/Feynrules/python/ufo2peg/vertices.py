import sys,pprint

from .helpers import CheckUnique,getTemplate,writeFile,qcd_qed_orders,def_from_model
from .lorentzparser import parse_lorentz
from .converter import py2cpp
from .collapse_vertices import collapse_vertices
from .check_lorentz import tensorCouplings,VVVordering,VVSCouplings,EWVVVVCouplings

# names of goldstone bosons
gsnames = ['goldstone','goldstoneboson','GoldstoneBoson']

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

class SkipThisVertex(Exception):
    pass

def get_lorentztag(spin):
    """Produce a ThePEG spin tag for the given numeric FR spins."""
    spins = { 1 : 'S', 2 : 'F', 3 : 'V', -1 : 'U', 5 : 'T' }
    result = [ spins[s] for s in spin ]

    def spinsort(a,b):
        """Helper function for ThePEG's FVST spin tag ordering."""
        if a == b: return 0
        for letter in 'UFVST':
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
        lname = l.name[:len(lorentztag)]
        if sorted(lorentztag) != sorted(lname):
            raise Exception("Lorentztags: %s is not %s in %s" 
                            % (lorentztag,lname,vertex))
#        if lorentztag != lname:
#            sys.stderr.write("Warning: Lorentz tag ordering: %s is not %s in %s\n"
#                             % (lorentztag,lname,vertex))

    return lorentztag

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

def colorfactor(vertex,L,pos):
    def match(patterns):
        result = [ p == t
                   for p,t in zip(patterns,vertex.color) ]
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
        label = ('Epsilon(1,2,3)',)
        if match(label): return ('1',) # TODO check factor!

    elif l(-3) == L == 3:
        label = ('EpsilonBar(1,2,3)',)
        if match(label): return ('1',) # TODO check factor!

    elif l(8) == L == 3:
        # if lorentz is FFV get extra minus sign
        lorentztag = unique_lorentztag(vertex)
        factor = '*(-1)' if lorentztag in ['FFV'] else ''
        label = ('f(1,2,3)',)
        if match(label): return ('-complex(0,1)%s'%factor,)
        label = ('f(3,2,1)',)
        if match(label): return ('complex(0,1)%s'%factor,)
        label = ('f(2,1,3)',)
        if match(label): return ('complex(0,1)%s'%factor,)

    elif l(8) == L == 4:
        label = ('f(-1,1,2)*f(3,4,-1)',
                 'f(-1,1,3)*f(2,4,-1)',
                 'f(-1,1,4)*f(2,3,-1)',
             )
        if match(label): return ('-1./3.','-1./3.','-1./3.')

    elif l(8) == 2 and l(3) == l(-3) == 1 and L==4:
        subs = {
            'g1' : pos[8][0],
            'g2' : pos[8][1],
            'qq' : pos[3][0],
            'qb' : pos[-3][0] 
        }
        label = ('T({g1},-1,{qb})*T({g2},{qq},-1)'.format(**subs),
                 'T({g1},{qq},-1)*T({g2},-1,{qb})'.format(**subs))
        if match(label): return ('0.5','0.5')
        
    elif l(8) == 2 and l(6) == l(-6) == 1 and L==4:
        subs = {
            'g1' : pos[8][0],
            'g2' : pos[8][1],
            'qq' : pos[6][0],
            'qb' : pos[-6][0] 
        }
        label = ('T6({g1},-1,{qb})*T6({g2},{qq},-1)'.format(**subs),
                 'T6({g1},{qq},-1)*T6({g2},-1,{qb})'.format(**subs))
        if match(label): return ('0.5','0.5')

    elif l(8) == 2 and l(8)+l(1)==L :
        subs = { 'g1' : pos[8][0], 'g2' : pos[8][1] }
        label = ('Identity({g1},{g2})'.format(**subs),)
        if match(label) : return ('1.',)

    elif l(8) == 3 and l(1)==1 and L==4 :
        label = ('f(1,2,3)',)
        if match(label): return ('-complex(0,1)',)
        label = ('f(3,2,1)',)
        if match(label): return ('complex(0,1)',)
        label = ('f(2,1,3)',)
        if match(label): return ('complex(0,1)',)

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
        def gstest(name):
            try:
                return getattr(p,name)
            except AttributeError:
                return False
        if any(map(gstest, gsnames)):
            return True
    return False

class VertexConverter:
    'Convert the vertices in a FR model to extract the information ThePEG needs.'
    def __init__(self,model) :
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

    def processVertex(self,vertexnumber,vertex) :
        # get the Lorentz tag for the vertex
        lorentztag = unique_lorentztag(vertex)
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
            plistarray = [ ','.join([ str(p.pdg_code) for p in vertex.particles ]) ]
        else:
            plistarray = [ ','.join([ str(p.pdg_code) for p in pl ])
                           for pl in vertex.particles_list ]
        # parse the colour structure for the vertex
        try:
            L,pos = colors(vertex)
            cf = colorfactor(vertex,L,pos)
        except SkipThisVertex:
            vertex.herwig_skip_vertex = True
            self.vertex_skipped=True
            return (True,"","")
    
        ### classname
        classname = 'V_%s' % vertex.name

        ### parse couplings
        unique_qcd = CheckUnique()
        unique_qed = CheckUnique()

        coup_left  = []
        coup_right = []
        coup_norm = []
        couplings_VVS = []
        kinematics = "false"
        if self.ONE_EACH:
            items = vertex.couplings.iteritems()
        else:
            items = vertex.couplings

        try:
            for (color_idx,lorentz_idx),coupling in items:

                qcd, qed = qcd_qed_orders(vertex, coupling)
                unique_qcd( qcd )
                unique_qed( qed )
                L = vertex.lorentz[lorentz_idx]
                prefactors = '(%s) * (%s) * (%s)' \
                                % (self.globalsign**(len(lorentztag)-2),lf,cf[color_idx])

                ordering = ''
                if lorentztag in ['FFS','FFV']:
                    left,right = parse_lorentz(L.structure)
                    if left:
                        coup_left.append('(%s) * (%s) * (%s)' % (prefactors,left,coupling.value))
                    if right:
                        coup_right.append('(%s) * (%s) * (%s)' % (prefactors,right,coupling.value))
                    if lorentztag == 'FFV':
                        ordering = ('if(p1->id()!=%s) {Complex ltemp=left(), rtemp=right(); left(-rtemp); right(-ltemp);}' 
                                    % vertex.particles[0].pdg_code)
                elif 'T' in lorentztag :
                    all_coup, left_coup, right_coup, ordering = \
                        tensorCouplings(vertex,coupling,prefactors,L,lorentztag,pos)
                    coup_norm  += all_coup
                    coup_left  += left_coup
                    coup_right += right_coup
                elif lorentztag == 'VVS' :
                    tc = VVSCouplings(vertex,coupling,prefactors,L,lorentztag)
                    if(len(couplings_VVS)==0) :
                        couplings_VVS=tc
                    else :
                        for ix in range(0,len(couplings_VVS)) :
                            if(tc[ix] == 0.) :
                                continue
                            elif(couplings_VVS[ix]==0.) :
                                couplings_VVS[ix]=tc[ix]
                            else :
                                couplings_VVS[ix] = '(( %s ) + ( %s ) )' % (couplings_VVS[ix],tc[ix])
                else:
                    if lorentztag == 'VSS':
                        if L.structure == 'P(1,3) - P(1,2)':
                            prefactors += ' * (-1)'
                        ordering = 'if(p2->id()!=%s){norm(-norm());}' \
                                       % vertex.particles[1].pdg_code
                    elif lorentztag == 'VVVV':
                        if qcd==2:
                            ordering = 'setType(1);\nsetOrder(0,1,2,3);'
                        else:
                            ordering, factor = EWVVVVCouplings(vertex,L)
                            prefactors += ' * (%s)' % factor
                    elif lorentztag == 'VVV':
                        if len(pos[8]) != 3:
                            ordering = VVVordering(vertex)
                    elif lorentztag == 'VVVS' :
                        if len(pos[8]) == 0 :
                            ordering = VVVordering(vertex)
                
                    if type(coupling) is not list:
                        value = coupling.value
                    else:
                        value = "("
                        for coup in coupling :
                            value += '+(%s)' % coup.value
                        value +=")"
                    coup_norm.append('(%s) * (%s)' % (prefactors,value))


                #print 'Colour  :',vertex.color[color_idx]
                #print 'Lorentz %s:'%L.name, L.spins, L.structure
                #print 'Coupling %s:'%C.name, C.value, '\nQED=%s'%qed, 'QCD=%s'%qcd
                #print '---------------'
        except SkipThisVertex:
            vertex.herwig_skip_vertex = True
            self.vertex_skipped=True
            return (True,"","")

        leftcontent  = ' + '.join(coup_left)  if coup_left  else '0'
        rightcontent = ' + '.join(coup_right) if coup_right else '0'
        normcontent  = ' + '.join(coup_norm)  if coup_norm  else '1'

        #print 'Left:',leftcontent
        #print 'Right:',rightcontent
        #print 'Norm:',normcontent
        #print '---------------'
        couplingptrs = [',tcPDPtr']*len(lorentztag)
        if lorentztag == 'VSS':
            couplingptrs[1] += ' p2'
        elif lorentztag == 'FFV':
            couplingptrs[0] += ' p1'
        elif (lorentztag == 'VVV' or lorentztag == 'VVVS' ) \
             and ordering != '':
            couplingptrs[0] += ' p1'
            couplingptrs[1] += ' p2'
            couplingptrs[2] += ' p3'
        elif lorentztag == 'VVVT' and ordering != '':
            couplingptrs[0] += ' p1'
            couplingptrs[1] += ' p2'
            couplingptrs[2] += ' p3'
        elif lorentztag == 'VVVV' and qcd != 2:
            couplingptrs[0] += ' p1'
            couplingptrs[1] += ' p2'
            couplingptrs[2] += ' p3'
            couplingptrs[3] += ' p4'
        
        ### do we need left/right?
        symbols = set()
        if 'FF' in lorentztag and lorentztag != "FFT":
            #leftcalc = aStoStrongCoup(py2cpp(leftcontent)[0], paramstoreplace_, paramstoreplace_expressions_)
            #rightcalc = aStoStrongCoup(py2cpp(rightcontent)[0], paramstoreplace_, paramstoreplace_expressions_)
            leftcalc, sym = py2cpp(leftcontent)
            symbols |= sym
            rightcalc, sym = py2cpp(rightcontent)
            symbols |= sym
            left = 'left(' + leftcalc + ');'
            right = 'right(' + rightcalc + ');'
        elif lorentztag == 'VVS':
            if(couplings_VVS[0]==0. and couplings_VVS[1]==0. and couplings_VVS[2]==0. and \
               couplings_VVS[3]==0. and couplings_VVS[4]==0. and couplings_VVS[5]!=0) :
                normcontent  = couplings_VVS[5]
                left=''
                right=''
            else :
                for ix in range(0,len(couplings_VVS)) :
                    if(couplings_VVS[ix] !=0.) :
                        couplings_VVS[ix], sym = py2cpp(couplings_VVS[ix])
                        symbols |= sym
                lorentztag = 'GeneralVVS'
                kinematics='true'
                # g_mu,nv piece of coupling 
                if(couplings_VVS[5]!=0.) :
                    left  = 'a00( %s + Complex(( %s )* GeV2/invariant(1,2)));' % ( couplings_VVS[0],couplings_VVS[5])
                else :
                    left  = 'a00( %s );' % couplings_VVS[0]
                # other couplings
                right = 'a11( %s );\n    a12( %s );\n    a21( %s );\n    a22( %s ); aEp( %s ); ' % \
                       ( couplings_VVS[1],couplings_VVS[2],couplings_VVS[3],couplings_VVS[4],couplings_VVS[6] )
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
                 'ordering' : ordering,
                 'kinematics' : kinematics
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
