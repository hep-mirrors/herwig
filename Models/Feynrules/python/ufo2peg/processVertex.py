import sys

from .helpers import CheckUnique,SkipThisVertex,colors,colorfactor,qcd_qed_orders,def_from_model,spindirectory,getTemplate,unique_lorentztag,tensorCouplings,VVVordering,VVSCouplings,EWVVVVCouplings
from .lorentzparser import parse_lorentz
from .converter import py2cpp

ONE_EACH = True

globalsign=False

gsnames = ['goldstone','goldstoneboson','GoldstoneBoson']

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

verbose=False
vertex_skipped=False
ignore_skipped=False

VERTEXHEADER = """\
#include "ThePEG/Helicity/Vertex/{spindirectory}/{lorentztag}Vertex.h"
"""

VERTEXCLASS = getTemplate('Vertex_class')


def should_print():
    return not vertex_skipped or ignore_skipped

### initial pass to find global sign
# at the moment does nothing
def global_sign(FR):
    global globalsign
    globalsign=1.
    # for v in FR.all_vertices:
    #     pids = sorted([ p.pdg_code for p in v.particles ])
    #     if pids != [-11,11,22]: continue
    #     coup = v.couplings
    #     assert( len(coup) == 1 )
    #     val = coup.values()[0].value
    #     val = evaluate(val)
    #     assert( val.real == 0 )
    #     if val.imag > 0 :
    #         globalsign= 1.
    #     else :
    #         globalsign=-1.

def checkGhostGoldstoneVertex(lorentztag,vertex) :
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

def processVertex(FR,modelname,vertexnumber,vertex) :
    global vertex_skipped
    # set the global sign if required
    if(not globalsign) : global_sign(FR)
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
        vertex_skipped=True
        return (True,"","")
    # get the ids of the particles at the vertex
    if ONE_EACH:
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
        vertex_skipped=True
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
    if ONE_EACH:
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
                            % (globalsign**(len(lorentztag)-2),lf,cf[color_idx])

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
        vertex_skipped=True
        return (True,"","")

    leftcontent  = ' + '.join(coup_left)  if coup_left  else '0'
    rightcontent = ' + '.join(coup_right) if coup_right else '0'
    normcontent  = ' + '.join(coup_norm)  if coup_norm  else '1'

#    #print 'Left:',leftcontent
#    #print 'Right:',rightcontent
#    #print 'Norm:',normcontent
#    #print '---------------'


        
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
        
#    ### do we need left/right?
    symbols = set()
    if 'FF' in lorentztag and lorentztag != "FFT":
#        leftcalc = aStoStrongCoup(py2cpp(leftcontent)[0], paramstoreplace_, paramstoreplace_expressions_)
#        rightcalc = aStoStrongCoup(py2cpp(rightcontent)[0], paramstoreplace_, paramstoreplace_expressions_)
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
#    normcalc = aStoStrongCoup(py2cpp(normcontent)[0], paramstoreplace_, paramstoreplace_expressions_)
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
    symboldefs = [ def_from_model(FR,s) for s in symbols ]
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
             'ModelName' : modelname,
             'ordering' : ordering,
             'kinematics' : kinematics
             }             # ok

    # print info if required
    if verbose:
        print '-'*60
        pprint.pprint(( classname, plistarray, leftcalc, rightcalc, normcalc ))

    return (False,VERTEXCLASS.substitute(subs),VERTEXHEADER.format(**subs))

