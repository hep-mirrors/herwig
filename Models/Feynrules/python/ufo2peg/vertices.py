import sys,pprint

from .helpers import CheckUnique,getTemplate,writeFile,unique_lorentztag,colors,colorfactor,qcd_qed_orders,\
    tensorCouplings,VVVordering,VVSCouplings,EWVVVVCouplings,def_from_model,spindirectory
from .lorentzparser import parse_lorentz
from .converter import py2cpp
from .collapse_vertices import collapse_vertices

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
