#! /usr/bin/env python
from __future__ import with_statement
import cmath, os, sys, fileinput, pprint
import argparse

from string import Template, strip
from helpers import CheckUnique, getTemplate, writeFile, get_lorentztag

# set up the option parser for command line input 
parser = argparse.ArgumentParser(description='Create Herwig++ model files from Feynrules UFO input.')
parser.add_argument('ufodir', metavar='UFO_directory', help='the UFO model directory')
parser.add_argument('-v', '--verbose', action="store_true", help="print verbose output")
parser.add_argument('-n','--name', default="FRModel", help="set custom nametag for the model")

args = parser.parse_args()

# check if the given "Model" path exists 
modeldir = args.ufodir.rstrip('/')
modelpath, module = os.path.split(modeldir)
if modelpath:
    sys.path.append(os.path.abspath(modelpath))

FR = __import__(module)

#print banner()
          

##################################################
##################################################
##################################################


def PyMathToThePEGMath(a,b):
    return a

def aStoStrongCoup(a,b,c):
    return a

def aEWtoWeakCoup(a,b,c):
    return a



##################################################
##################################################
##################################################

# get templates for Model header and .cc file,
# Herwig++ run input file
MODEL_H  = getTemplate('Model.h')
MODEL_CC = getTemplate('Model.cc')
MODEL_HWIN = getTemplate('LHC-FR.in')

# get the Model name from the arguments
modelname = args.name

# copy the Makefile-FR to current directory,
# replace with the modelname for compilation
with open('../Makefile-FR','r') as orig:
    with open('Makefile','w') as dest:
        dest.write(orig.read().replace("FeynrulesModel.so", modelname+".so"))


# define arrays and variables     
allplist = ""
parmdecls = []
parmgetters = []
parmconstr = []

parmfuncmap = []
paramsforev = []
paramstoreplace_ = []
paramstoreplace_expressions_ = []

# get external parameters for printing
parmsubs = dict( [ (p.name, float(p.value)) 
                   for p in FR.all_parameters 
                   if p.nature == 'external' ] ) 

# evaluate python cmath
def evaluate(x):
    return eval(x, 
                {'cmath':cmath,
                 'complexconjugate':FR.function_library.complexconjugate}, 
                parmsubs)


# get internal, external and all params into arrays
internal = [ p 
             for p in FR.all_parameters 
             if p.nature == 'internal' ]

external = [ p 
             for p in FR.all_parameters 
             if p.nature == 'external' ]

allparams =  [ p.name 
               for p in FR.all_parameters ]

paramstoreplaceEW_ = []
paramstoreplaceEW_expressions_ = []
# calculate internal parameters
for p in internal:
    if 'aS' in p.value and p.name != 'aS':
        paramstoreplace_.append(p.name)
        paramstoreplace_expressions_.append(p.value)
    if 'aEWM1' in p.value and p.name != 'aEWM1':
        paramstoreplaceEW_.append(p.name)
        paramstoreplaceEW_expressions_.append(p.value)
    newval = evaluate(p.value)
    parmsubs.update( { p.name : newval } )
        
    
# more arrays used for substitution in templates 
paramvertexcalc = []
paramsforstream = []
parmmodelconstr = []

# loop over parameters and fill in template stuff according to internal/external and complex/real
# WARNING: Complex external parameter input not tested!
if args.verbose:
    print 'verbose mode on: printing all parameters'
    print '-'*60
    paramsstuff = ('name', 'expression', 'default value', 'nature')
    pprint.pprint(paramsstuff)



interfacedecl_T = Template(
"""
static Parameter<$modelname, $type> interface$pname
  ("$pname",
   "The interface for parameter $pname",
   &$modelname::$pname, $value, 0, 0,
   false, false, Interface::nolimits);
"""
)

interfaceDecls = []

typemap = {'complex':'Complex',
           'real':'double'}

for parmnumber,p in enumerate(FR.all_parameters):
    value = parmsubs[p.name]

    if p.nature == 'external':
        interfaceDecls.append( 
            interfacedecl_T.substitute(modelname=modelname,
                                       pname=p.name,
                                       value=value,
                                       type=typemap[p.type]) 
        )

    if p.type == 'real':
        assert( value.imag < 1.0e-16 )
        value = value.real
        parmconstr.append('%s(%s)' % (p.name, value))
        if p.nature == 'external':
            parmmodelconstr.append('set %s:%s %s' % (modelname, p.name, value))
    elif p.type == 'complex':
        value = complex(value)
        parmconstr.append('%s(%s,%s)' % (p.name, value.real, value.imag))
        if p.nature == 'external':
            parmmodelconstr.append('set %s:%s (%s,%s)' % (modelname, p.name, value.real, value.imag))
    else:
        raise Exception('Unknown data type "%s".' % p.type)

    parmsubs[p.name] = value
    parmdecls.append('  %s %s;' % (typemap[p.type], p.name))
    parmgetters.append('  %s %s_() const { return %s; }' % (typemap[p.type],p.name, p.name))
    parmfuncmap.append('   case %s:  return %s_();' % (parmnumber, p.name))
    paramsforev.append('%s' % p.name)
    paramsforstream.append('%s' % p.name)

    if p.name == 'aS':
        funcvertex = '0.25 * sqr(strongCoupling(q2)) / Constants::pi'
    elif p.name == 'aEWM1':
        funcvertex = '4.0 * Constants::pi / sqr(electroMagneticCoupling(q2))'
    elif p.name == 'Gf':
        funcvertex = 'generator()->standardModel()->fermiConstant() * GeV2'
    elif p.name == 'MZ':
        funcvertex = 'getParticleData(ThePEG::ParticleID::Z0)->mass() / GeV'
    else:
        funcvertex = 'hw%s_ptr->%s_()' % (modelname, p.name)
    if p.lhablock == None:
        funcvertex = PyMathToThePEGMath(p.value, allparams)
    paramvertexcalc.append('%s = %s;' % (p.name,funcvertex))

    if args.verbose:
        pprint.pprint((p.name,p.value, value, p.nature))

parmtextsubs = { 'parmgetters' : '\n'.join(parmgetters),
                 'parmdecls' : '\n'.join(parmdecls),
                 'parmconstr' : ': ' + ',\n  '.join(parmconstr),
                 'getters' : '',
                 'decls' : '',
                 'addVertex' : '',
                 'ostream' : '\n\t<< '.join(paramsforstream),
                 'istream' : '\n\t>> '.join(paramsforstream),
                 'refs' : '',
                 'parmextinter': ''.join(interfaceDecls),
                 'num_params': len(FR.all_parameters),
                 'parmfuncmap': '\n'.join(parmfuncmap),
                 'paramsforev': ','.join(paramsforev),
                 'ModelName': modelname
                 }

print '-'*60

# write the files from templates according to the above subs
writeFile( modelname + '.h', MODEL_H.substitute(parmtextsubs) )
writeFile( modelname +'.cc', MODEL_CC.substitute(parmtextsubs) )
writeFile( 'LHC-' + modelname +'.in', MODEL_HWIN.substitute(parmtextsubs) )

##################################################
##################################################
##################################################

# ignore these, they're in Hw++ already # TODO reset Hw++ settings instead
SMPARTICLES = {

1:'d',
2:'u',
3:'s',
4:'c',
5:'b',
6:'t',

11:'e-',
12:'nu_e',
13:'mu-',
14:'nu_mu',
15:'tau-',
16:'nu_tau',

21:'g',
22:'gamma',
23:'Z0',
24:'W+',

-1:'dbar',
-2:'ubar',
-3:'sbar',
-4:'cbar',
-5:'bbar',
-6:'tbar',

-11:'e+',
-12:'nu_ebar',
-13:'mu+',
-14:'nu_mubar',
-15:'tau+',
-16:'nu_taubar',

-24:'W-',

}



particleT = Template(
"""
create ThePEG::ParticleData $name
setup $name $pdg_code $name $mass $width $wcut $ctau $charge $color $spin 0
insert /Herwig/NewPhysics/NewModel:DecayParticles 0 $name
"""
)
class ParticleConverter:
    'Convert a FR particle to extract the information ThePEG needs.'
    def __init__(self,p):
        self.name = p.name
        self.pdg_code = p.pdg_code
        self.spin = p.spin
        self.color = p.color
        self.selfconjugate = 0
        if self.color == 1:
            self.color = 0
        self.mass = parmsubs[str(p.mass)]
        self.width = parmsubs[str(p.width)]
        try:
            self.mass = self.mass.real
        except:
            pass
        hbarc = 197.3269631e-15 # GeV mm (I hope ;-) )
        if self.width != 0: self.ctau = hbarc / self.width
        else:               self.ctau = 0
        self.wcut = 10 * self.width
        self.charge = int(3 * p.charge)

    def subs(self):
        return self.__dict__

def get_all_thepeg_particles():
    plist = ''
    antis = {}
    for p in FR.all_particles:
        if p.spin == -1 or p.goldstoneboson:
            continue

        if p.pdg_code in SMPARTICLES:
            #add stuff to plist to set params
            pass
        else:
            if p.pdg_code == 25:
                plist += """
set /Herwig/Particles/h0:Mass_generator NULL
set /Herwig/Particles/h0:Width_generator NULL
rm /Herwig/Masses/HiggsMass
rm /Herwig/Widths/HiggsWidth
"""
            subs = ParticleConverter(p).subs()
            plist += particleT.substitute(subs)

            pdg, name = subs['pdg_code'],  subs['name']
            if -pdg in antis:
                plist += 'makeanti %s %s\n' % (antis[-pdg], name)
                
            else:
                antis[pdg] = name
                selfconjugate = 1
    return plist





##################################################
##################################################
##################################################

# count number of vertices not included
notincluded = 0

# get vertex template
VERTEX = getTemplate('Vertex.cc')

def produce_vertex_file(subs):
    newname = modelname + subs['classname'] + '.cc'
    writeFile( newname, VERTEX.substitute(subs) )

if args.verbose:
    #print 'vertex\tLorentz\t\t\tC_L\t\t\tC_R\t\t\t\tnorm\t'
    print 'verbose mode on: printing all vertices'
    print '-'*60
    labels = ('vertex', 'particles', 'Lorentz', 'C_L', 'C_R', 'norm')
    pprint.pprint(labels)


for v in FR.all_vertices:

    #print v.name
    #print map(str,v.particles)
    #print '---------------'
    v.include = 1

    ### Spin structure
    unique = CheckUnique()
    for l in v.lorentz:
        lt = get_lorentztag(l.spins)
        unique( lt )

    if 'T' in lt:   spin_directory = 'Tensor'
    elif 'S' in lt: spin_directory = 'Scalar'
    elif 'V' in lt: spin_directory = 'Vector'
    elif 'U' in lt: spin_directory = 'Ghost'
    
    ### Particle ids #################### sort order? ####################
    plistarray = ['','']    
    plistarray[0] = ','.join([ str(p.pdg_code) for p in v.particles ])
    plist = ','.join([ str(p.pdg_code) for p in v.particles ])
    #print plist

# Check if the Vertex is self-conjugate or not
    pdgcode = [0,0,0,0]
    notsmvertex = False
    vhasw = 0
    vhasz = 0
    vhasf = 0
    vhasg = 0
    vhash = 0
    vhasp = 0
#   print 'printing particles in vertex'
    for i in range(len(v.particles)):
#       print v.particles[i].pdg_code
        pdgcode[i] = v.particles[i].pdg_code
        if pdgcode[i] == 23:
            vhasz += 1
        if pdgcode[i] == 22:
            vhasp += 1
        if pdgcode[i] == 25:
            vhash += 1
        if pdgcode[i] == 21:
            vhasg += 1
        if pdgcode[i] == 24:
            vhasw += 1
        if abs(pdgcode[i]) < 7 or (abs(pdgcode[i]) > 10 and abs(pdgcode[i]) < 17):
            vhasf += 1
        if pdgcode[i] not in SMPARTICLES:
            notsmvertex = True
        

#  treat replacement of SM vertices with BSM vertices?               
    if notsmvertex == False:
        if( (vhasf == 2 and vhasz == 1) or (vhasf == 2 and vhasw == 1) or (vhasf == 2 and vhash == 1) or (vhasf == 2 and vhasg == 1) or (vhasf == 2 and vhasp == 0) or (vhasg == 3) or (vhasg == 4) or (vhasw == 2 and vhash == 1) or (vhasw == 3) or (vhasw == 4) or (vhash == 1 and vhasg == 2) or (vhash == 1 and vhasp == 2)):
            #print 'VERTEX INCLUDED IN STANDARD MODEL!'
            v.include = 0
            notincluded += 1
            continue
            
    
    selfconjugate = 0
    for j in range(len(pdgcode)):
        for k in range(len(pdgcode)):
               if  j != k and j != 0 and abs(pdgcode[j]) == abs(pdgcode[k]):
                   selfconjugate = 1
                   #print 'self-conjugate vertex'
#        print pdgcode[j]

# if the Vertex is not self-conjugate, then add the conjugate vertex
# automatically
    scfac = [1,1,1,1]
    if selfconjugate == 0:
        #first find the self-conjugate particles
        for u in range(len(v.particles)):
              if v.particles[u].selfconjugate == 0:
                  scfac[u] = -1
#                  print 'particle ', v.particles[u].pdg_code, ' found not to be self-conjugate'
                  
    if selfconjugate == 0:
        plistarray[1] += str(scfac[1] * v.particles[1].pdg_code) + ',' + str(scfac[0] * v.particles[0].pdg_code) + ',' + str(scfac[2] * v.particles[2].pdg_code)
        if len(v.particles) is 4:                                                                                                                      
            plistarray[1] += ',' + str(scfac[3] * v.particles[3].pdg_code)
        #print 'Conjugate vertex:', plistarray[1]
    
    ### Colour structure
    if v.color == '1': qcdord = '0'
    else:              qcdord = ''

    ### classname
    classname = 'V_%03d' % int(v.name[2:])

    ### parse couplings
    unique_qcd = CheckUnique()
    unique_qed = CheckUnique()
    
    coup_left  = []
    coup_right = []
    coup_norm = []
 
    for (ci,li),C in v.couplings.iteritems():
        qed = C.order.get('QED',0)
        qcd = C.order.get('QCD',0)
        # WARNING: FIX FOR CASES WHEN BOTH ARE ZERO
        # Is there a better way to treat this?
        if qed == 0 and qcd == 0:
            qed = 1
        unique_qcd( qed )
        unique_qed( qcd )
        L = v.lorentz[li]

        if lt in ['FFS','FFV']:
            #print 'PRINTING LORENTZ STRUCTURE'
            #print L.structure
            for lor in map(strip, L.structure.split('+')):
                breakdown = lor.split('*')
                prefactor='1'
                #print 'breakdown', breakdown, 'length', len(breakdown)
                if len(breakdown) == 3:
                    prefactor = breakdown[0]
                    breakdown = breakdown[1:]
                if len(breakdown) == 2:
                    assert(breakdown[0][:5] == 'Gamma')
                    if breakdown[1][:5] == 'ProjM':
                        #print 'LEFT HANDED'
                        coup_left.append(prefactor+' * '+C.value)
                    elif breakdown[1][:5] == 'ProjP':
                        #print 'RIGHT HANDED'
                        coup_right.append(prefactor+' * '+C.value)
                    else:
                        coup_left.append(C.value)
                        coup_right.append(C.value)
                if len(breakdown) == 1:
                    if breakdown[0][:5] == 'ProjM':
                        #print 'LEFT HANDED'
                        coup_left.append(prefactor+' * '+C.value)
                    elif breakdown[0][:5] == 'ProjP':
                        #print 'RIGHT HANDED'
                        coup_right.append(prefactor+' * '+C.value)
                    else:
                        coup_left.append(C.value)
                        coup_right.append(C.value)
        else:
            coup_norm.append(C.value)
                
            
        #print 'Colour  :',v.color[ci]
        #print 'Lorentz %s:'%L.name, L.spins, L.structure
        #print 'Coupling %s:'%C.name, C.value, '\nQED=%s'%qed, 'QCD=%s'%qcd
        #print '---------------'


    leftcontent = ' + '.join(coup_left) if len(coup_left)!=0 else '0j'
    rightcontent = ' + '.join(coup_right) if len(coup_right)!=0 else '0j'
    normcontent = ' + '.join(coup_norm) if len(coup_norm)!=0 else '1.'

    #print 'Left:',leftcontent
    #print 'Right:',rightcontent
    #print 'Norm:',normcontent
    #print '---------------'


    #leftexplicit = complex(evaluate(leftcontent))
    #rightexplicit = complex(evaluate(rightcontent))
    #    normexplicit = complex(evaluate(normcontent))
    leftdebug = ''
    rightdebug = ''
    normdebug = ''
    ### do we need left/right?
    if 'FF' in lt:
        leftcalc = aStoStrongCoup(PyMathToThePEGMath(leftcontent, allparams), paramstoreplace_, paramstoreplace_expressions_)
        rightcalc = aStoStrongCoup(PyMathToThePEGMath(rightcontent, allparams), paramstoreplace_, paramstoreplace_expressions_)
        left = 'left(' + leftcalc + ');'
        right = 'right(' + rightcalc + ');'
        #leftdebug  = 'left(Complex(%s,%s));'  % (leftexplicit.real,leftexplicit.imag)
        #rightdebug = 'right(Complex(%s,%s));' % (rightexplicit.real,rightexplicit.imag)
    else:
        left = ''
        right = ''
        leftcalc = ''
        rightcalc = ''
        leftdebug = ''
        rightdebug = ''
        

    normcalc = aStoStrongCoup(PyMathToThePEGMath(normcontent, allparams), paramstoreplace_, paramstoreplace_expressions_)
    norm = 'norm(' + normcalc + ');'
    #normdebug = 'norm(Complex(%s,%s));' % (normexplicit.real,normexplicit.imag)

    if plistarray[1] is '':
        plist2 = ''
    else:
        plist2 = 'addToList(%s);' % plistarray[1]


    # input q2 or not, depending on whether it is necessary
    if 'q2' in norm or 'q2' in left or 'q2' in right:
        q2var = ' q2'
    else:
        q2var = ''
    
        
    ### assemble dictionary and fill template
    subs = { 'lorentztag' : lt,                   # ok
             'classname'  : classname,            # ok
             'left'       : left,                 # doesn't always exist in base
             'right'      : right,                 # doesn't always exist in base 
             'norm'      : norm,                 # needs norm, too

             #################### need survey which different setter methods exist in base classes

             'addToPlist' : 'addToList(%s);' % plistarray[0], # ok
             'addToPlist2' : plist2, # ok
             'parameters' : '',
             'setCouplings' : '',
             'qedorder'   : qed,
             'qcdorder' : qcd,
             'q2'        :  q2var,
             'couplingptrs' : ',tcPDPtr'*len(v.particles),
             'spindirectory' : spin_directory,
             'ModelName' : modelname,
             'num_params' : len(FR.all_parameters),
             'leftcontent' : leftcontent,
             'rightcontent' : rightcontent,
             'normcontent' : normcontent,
             'leftcalc': leftcalc,
             'rightcalc' : rightcalc,
             'normcalc' : normcalc,
             'paramdecl': '\n'.join(parmdecls),
             'ostream' : ' << '.join(paramsforstream),
             'istream' : ' >> '.join(paramsforstream),
             'paramvertexcalc': '\n\t'.join(paramvertexcalc),
             'leftdebug': leftdebug,
             'rightdebug' : rightdebug,
             'normdebug' : normdebug
             }             # ok

    if args.verbose:
        print '-'*60
        if  selfconjugate:
            stuff = ( classname, plistarray[0], leftcalc.replace('Complex(0,1.)','i').replace('Complex(0,0)','0'), rightcalc.replace('Complex(0,1.)','i').replace('Complex(0,0)','0'), normcalc.replace('Complex(0,1.)','i').replace('Complex(0,0)','0') )
        else:
            stuff = ( classname, plistarray[0], plistarray[1], leftcalc.replace('Complex(0,1.)','i').replace('Complex(0,0)','0'), rightcalc.replace('Complex(0,1.)','i').replace('Complex(0,0)','0'), normcalc.replace('Complex(0,1.)','i').replace('Complex(0,0)','0') )
        pprint.pprint(stuff)
        
        
    if  L.spins[0] != -1 and L.spins[1] != -1 and L.spins[2] != -1 and plistarray[0] not in allplist and plistarray[1] not in allplist:
        produce_vertex_file(subs)
        allplist += plistarray[0]
        allplist += plistarray[1]
    elif  L.spins[0] != -1 and L.spins[1] != -1 and L.spins[2] != -1 and selfconjugate:
        produce_vertex_file(subs)
        allplist += plistarray[0]
    else:
        #print 'VERTEX ALREADY INCLUDED'
        v.include = 0
        
        #print '============================================================'

print '='*60

##################################################
##################################################
##################################################

vertexline = Template("""\
create $classname $name
insert ${ModelName}:ExtraVertices 0 $name
""")


def get_vertices():
    vlist = 'library ' + modelname + '.so\n'
    for v in FR.all_vertices:
        for l in v.lorentz:
            lt = get_lorentztag(l.spins)
            #print lt
        if "U" not in lt and v.include == 1:
            vlist += vertexline.substitute(
                { 'classname' : 'Herwig::%sV_%03d' % (modelname, int(v.name[2:])),
                  'name' : '/Herwig/%s/%s' % (modelname,v.name), 
                  'ModelName' : modelname } )
    return vlist


modelfilesubs = { 'plist' : get_all_thepeg_particles(),
                  'vlist' : get_vertices(),
                  'setcouplings': '\n'.join(parmmodelconstr),
                  'ModelName': modelname
                  }

#print get_all_thepeg_particles()

MODELINFILE = getTemplate('FR.model')

writeFile( modelname +'.model', MODELINFILE.substitute(modelfilesubs) )

print 'finished generating model:\t', modelname
print 'model directory:\t\t', args.ufodir
print 'generated:\t\t\t', len(FR.all_vertices)-notincluded, 'vertices'
print '='*60
print 'library:\t\t\t', modelname +'.so'
print 'input file:\t\t\t', 'LHC-' + modelname +'.in'
print 'model file:\t\t\t', modelname +'.model'
print '='*60
print """To complete the installation, compile by typing "make", 
copy the generated .so file, into the Herwig++ lib directory
and the .model and .in files in the directory that you wish to run in.
"""
print 'DONE!'
print '='*60

