#! /usr/bin/env python
from __future__ import with_statement
import cmath, string, os, sys, fileinput
from optparse import OptionParser


# set up the option parser for command line input 
parser = OptionParser(usage="%prog -d [UFO model directory] -n [custom model nametag]")

parser.add_option("-d", "--directory", dest="MODELDIR",
                  default="Model", help="UFO model directory.")

parser.add_option("-n", "--name", dest="MODELNAME",
                  default="FeynRulesModel", help="custom model nametag.")
opts, args = parser.parse_args()

# check if the given "Model" path exists 
if(os.path.exists(os.getcwd() + '/'+opts.MODELDIR) is False):
    print 'Path', os.getcwd() + '/'+ opts.MODELDIR, 'does not exist, exiting'
    sys.exit()

# if the Model path exists, then import the UFO FeynRules module
FR = __import__(opts.MODELDIR)


print 'generating Model and Vertex files/.model file/.input file'
print 'please be patient!'
print '-------------------'

# check if a function is a number
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# function that replaces ** with pow(,): used in PyMathToThePEGMath below
def powstring(stringrep,power):
    #if(stringrep[0] == '('):
    #    return 'pow(' + stringrep[0:len(stringrep)-2] + stringrep[len(stringrep)-3:len(stringrep)].replace('**','') + ',' + power + ')'
    #else:
    return 'pow(' + stringrep[0:len(stringrep)-3] + stringrep[len(stringrep)-3:len(stringrep)+1].replace('**','') + ',' + power + ')'

def BrackParams(stringin):
    # to take care of the ** powers over brackets (), append new "parameters" of type (xxxx)
    #	th = acos(1/sqrt(1 + (-pow(MM1,2) + pow(MM2,2) + sqrt(4*pow(MM12,4) + (pow(MM1,2) - pow(MM2,2))**2))**2/(4.*pow(MM12,4))));
    # start by finding the positions of the left and right brackets 
    stringbrack = stringin
    d_string_length = len(stringbrack)
    brack_pos_left = []
    brack_pos_right = []
    for ll in range(0,d_string_length):
        if(stringbrack[ll] is '('):
            brack_pos_left.append(ll)
        if(stringbrack[ll] is ')'):
            brack_pos_right.append(ll)

    paramsin_brack = []
    # loop over the left bracket positions, moving left, and count the number of right brackets
    for le in brack_pos_left:
        LRNUM = -1
        # print 'brack_pos_left', le
        # loop over the input string starting from the position of the bracket
        # count the left and right brackets to the right of the left bracket
        # left brackets are negative, right brackets are positive
        # the matched bracket is found when the left + right = LRNUM = 0
        for lm in range(le+1,d_string_length):
            if(lm in brack_pos_left):
                LRNUM = LRNUM - 1
            if(lm in brack_pos_right):
                LRNUM = LRNUM + 1
            if(LRNUM is 0):
                #print 'BRACKET', le, 'matched with', lm, stringbrack
                #print 'appending', stringbrack[le:lm+1]
                paramsin_brack.append(stringbrack[le:lm+1])
                break
            
    # sort the paramsin_brack according to length, from longest to shortest
    # insert to start of string
    paramsin_brack.sort(key=len, reverse=True)
    # for iii in range(0,len(paramsin_brack)):
    #paramsin.insert(0,paramsin_brack[iii])
    return paramsin_brack

# function that converts the vertex expressions in Python math format:
# called twice: once for the parameters, and once for terms in parentheses (...)
def PyMathToThePEGMath(stringin, paramsin):
    stringout = PyMathToThePEGMath_nb(stringin, paramsin)
    paramsbrack = BrackParams(stringout)
    #print 'paramsbrack', paramsbrack
    if(paramsbrack is not []):
        stringreturn =  PyMathToThePEGMath_nb(stringout, paramsbrack)
    return stringreturn

# function that converts the vertex expressions in Python math format to
# format that can be calculated using ThePEG 
def PyMathToThePEGMath_nb(stringin, paramsin):
    
    # define an array that contains the numbers 0-9 in string form
    numbersarray = ['0','1','2','3','4','5','6','7','8','9']

    # define counters and variables used to detect the positions of powers
    ii = 0
    pos = ''
    posnew = ''
    powpos = ''
    pow_ddg = 0
    power = ''
  

    # add '**' to the end of paramsin[ss], the array of given parameters of the model
    for ss in range(0,len(paramsin)):
        paramsin[ss] = paramsin[ss] + '**'

        #print 'in progress', stringin

    #print 'paramsin', paramsin
    powerchange = []
    # loop over the array of the model parameters and search for them in the given mathematical expression
    # each time a new position with the ** notation is found, replace with the C++ pow(,) notation
    for xx in range(0,len(paramsin)):
        # powerchange contains information on the positions
        # of necessary changes. OBSOLETE: for testing purposes only
        powerchange.append([])
        # reset counter for next variable and position variables
        ii = 0
        pos = 0
        posnew = 0
        # save the length of the string at the beginning of the loop for a parameter
        initial_string_length = len(stringin)
        # scan the string from right to left
        while (initial_string_length-ii >= 0):
            # set the new position of the found 
            posnew = stringin.find(paramsin[xx],initial_string_length-ii)
            #print 'param found', posnew, paramsin[xx]
            # if the position is new, do stuff
            if(posnew is not pos and posnew is not -1):
                # get the position of the power
                powpos = posnew + len(paramsin[xx]) 
                # check if power is single or double digit
                # i.e. -> assuming there are no powers beyond "99"
                power = stringin[powpos]
                if(powpos+1 < initial_string_length):
                    if(stringin[powpos+1] in numbersarray):
                        #print 'power is double digit ', (stringin[powpos+1])
                        pow_ddg = 1
                        power = stringin[powpos] + stringin[powpos+1]
                powerchange[xx].append([ posnew, power ])
                # do the replacement of the ** to pow(,)
                stringin = stringin[:posnew] + stringin[posnew:posnew+len(paramsin[xx])+len(power)].replace(paramsin[xx]+power,powstring(paramsin[xx],power)) + stringin[posnew+len(paramsin[xx])+len(power):]
                #print 'in progress', stringin
            # reset position variable for next point in string
            # increment the counter for the position in string
            pos = posnew
            ii += 1
    # do replacements of 'complex'
    # integers multiplying stuff (add the "."
    stringin = stringin.replace('0j','Complex(0,0)')
    stringin = stringin.replace('complex(0,1)','Complex(0,1.)')
    stringin = stringin.replace('complex','Complex')
    stringin = stringin.replace('cmath.pi', 'M_PI')
    stringin = stringin.replace('cmath.', '')
    stringin = stringin.replace('-(','(-1.)*(')
    for nn in range(0,len(numbersarray)):
        numbersnn = numbersarray[nn]
        stringin = stringin.replace(numbersnn +' *', numbersnn +'. *')
        stringin = stringin.replace(numbersnn+'*Complex',numbersnn+'.*Complex')
        
   
    # print 'final string:'
    # print 'final string', stringin, paramsin
    # reset the parameters with ** for next run of function and return
    for ss in range(0,len(paramsin)):
        paramsin[ss] = paramsin[ss].replace('**','')
    return stringin

# function that replaces alphaS (aS)-dependent variables
# with their explicit form which also contains strongCoupling
def aStoStrongCoup(stringin, paramstoreplace, paramstoreplace_expressions):
    #print stringin
    for xx in range(0,len(paramstoreplace)):
        #print paramstoreplace[xx], paramstoreplace_expressions[xx]
        stringout = stringin.replace(paramstoreplace[xx], '(' +  PyMathToThePEGMath(paramstoreplace_expressions[xx],allparams) + ')')
    stringout = stringout.replace('aS', '(sqr(strongCoupling(q2))/(4.0*Constants::pi))')
    #print 'resulting string', stringout
    return stringout


# function that replaces alphaEW (aEW)-dependent variables
# with their explicit form which also contains weakCoupling
def aEWtoWeakCoup(stringin, paramstoreplace, paramstoreplace_expressions):
    #print stringin
    for xx in range(0,len(paramstoreplace)):
        #print paramstoreplace[xx], paramstoreplace_expressions[xx]
        stringout = stringin.replace(paramstoreplace[xx], '(' +  PyMathToThePEGMath(paramstoreplace_expressions[xx],allparams) + ')')
    stringout = stringout.replace('aEWM1', '(1/(sqr(electroMagneticCoupling(q2))/(4.0*Constants::pi)))')
    #print 'resulting string', stringout
    return stringout
          
          
# function to get template
def getTemplate(basename):
    with open('../%s.template' % basename, 'r') as f:
        templateText = f.read()
    return string.Template( templateText )

# write a filename
def writeFile(filename, text):
    with open(filename,'w') as f:
        f.write(text)

class CheckUnique:
    def __init__(self):
        self.val = None

    def __call__(self,val):
        if self.val is None:
            self.val = val
        else:
            assert( val == self.val )


##################################################
##################################################
##################################################

# get templates for Model header and .cc file,
# Herwig++ run input file
MODEL_H  = getTemplate('Model.h')
MODEL_CC = getTemplate('Model.cc')
MODEL_HWIN = getTemplate('LHC-FR.in')

# get the Model name from the arguments
ModelName = opts.MODELNAME

# copy the Makefile-FR to current directory,
# replace with the ModelName for compilation
os.system('cp ../Makefile-FR Makefile')
for line in fileinput.input('Makefile', inplace = 1):
      print line.replace("FeynRulesModel.so", ModelName+".so"),


# define arrays and variables     
allplist = ""
parmdecls = []
parmgetters = []
parmconstr = []
parmextinter = []
parmfuncmap = []
paramsforev = []
paramstoreplace_ = []
paramstoreplace_expressions_ = []

# get external parameters for printing
parmsubs = dict( [ (p.name, float(p.value)) 
                   for p in FR.all_parameters 
                   if p.nature == 'external' ] ) 


#print parmsubs
#print

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

#print 'external parms:'
#print external
#print

#print 'internal parms:'
#print internal
#print
paramstoreplaceEW_ = []
paramstoreplaceEW_expressions_ = []
# calculate internal parameters
for p in internal:
    #print p.name,'=',p.value
    if('aS' in p.value and p.name is not 'aS'):
        #print 'PARAM', p.name, 'contains aS'
        #print p.value
        paramstoreplace_.append(p.name)
        paramstoreplace_expressions_.append(p.value)
    if('aEWM1' in p.value and p.name is not 'aEWM1'):
        #print 'PARAM', p.name, 'contains aEW'
        #print p.value
        paramstoreplaceEW_.append(p.name)
        paramstoreplaceEW_expressions_.append(p.value)
        #if(is_number(p.value)):
    newval = evaluate(p.value)
    parmsubs.update( { p.name : newval } )
        

#print parmsubs
#print

# put external parameters into list of parameters to be interfaced
for p in external:
    #print p.name,'=',p.value
    extinter = '%s' % (p.name)
 
#print 'NUMBER OF PARAMS', len(FR.all_parameters)
#print 'PARAMETER NAMES'
    
# more arrays used for substitution in templates 
paramvertexcalc = []
paramsforstream = []
parmmodelconstr = []
parmnumber = 0

# loop over parameters and fill in template stuff according to internal/external and complex/real
# WARNING: Complex external parameter input not tested!
for p in FR.all_parameters:
    value = parmsubs[p.name]
    extinter = ''
    #print p.name
    if (p.nature == 'external' and p.type == 'real'):
    #extinter = '%s' % (p.name)
       extinter = 'static Parameter<%s, double> interfaceg%s' % (ModelName, p.name)
       extinter += '\n'+' ("%s",' % (p.name)
       extinter += '\n'+' "The interface to the parameter %s",' % (p.name)
       extinter += '\n'+' &%s::%s, %s, -10000., 10000.,' % (ModelName, p.name, value)
       extinter += '\n'+' false, false, Interface::limited);\n'
    if (p.nature == 'external' and p.type == 'complex'):
        #extinter = '%s' % (p.name)
       extinter = 'static Parameter<%s, Complex> interfaceg%s' % (ModelName, p.name)
       extinter += '\n'+' ("%s",' % (p.name)
       extinter += '\n'+' "The interface to the parameter %s",' % (p.name)
       extinter += '\n'+' false, false, Interface::limited);\n'
    if p.type == 'real':
        try:
            assert( value.imag < 1.0e-16 )
            value = value.real
        except:
            pass
        parmsubs[p.name] = value
        decl = '  double %s;' % p.name
        constr = '%s(%s)' % (p.name, value)
        if(p.nature == 'external'):
            modelconstr = 'set ' + ModelName + ':%s %s' % (p.name, value) 
        getter = '  double %s_() const { return %s; }' % (p.name, p.name)
        funcmap = '   case %s:  return %s_();' % (parmnumber, p.name)
        forev = '%s' % p.name
        funcvertex = '%s = hw%s_ptr->%s_();' % (p.name, ModelName, p.name)
        parmnumber += 1
    elif p.type == 'complex':
        value = complex(value)
        parmsubs[p.name] = value
        decl = '  Complex %s;' % p.name
        constr = '%s(%s,%s)' % (p.name, value.real, value.imag)
        if(p.nature == 'external'):
            modelconstr = 'set ' + ModelName + ':%s (%s,%s)' % (p.name, value.real, value.imag)
        getter = '  Complex %s_() const { return %s; }' % (p.name, p.name)
        funcmap = '   case %s:  return %s_();' % (parmnumber, p.name)
        forev = '%s' % p.name
        funcvertex = '%s = hw%s_ptr->%s_();' % (p.name, ModelName, p.name)
        parmnumber += 1
    else:
        raise Exception('Unknown data type "%s".' % p.type)
    if(p.name == 'aS'):
        funcvertex = '%s = (sqr(strongCoupling(q2))/(4.0*Constants::pi));' % p.name
    if(p.name == 'aEWM1'):
        funcvertex = '%s = ((4.0*Constants::pi)/sqr(electroMagneticCoupling(q2)));' % p.name
    if(p.name == 'Gf'):
        funcvertex = '%s = generator()->standardModel()->fermiConstant()*GeV*GeV;' % p.name
    if(p.name == 'MZ'):
        funcvertex = '%s = getParticleData(ThePEG::ParticleID::Z0)->mass()/GeV;' % p.name
    if(p.lhablock == None):
        funcvertex = p.name +' = ' + PyMathToThePEGMath(p.value, allparams) + ';' 
        #print 'NO LHABLOCK:', p.name, funcvertex

    # do calc in C++, add interfaces for externals
    paramvertexcalc.append(funcvertex)    
    parmdecls.append(decl)
    parmgetters.append(getter)
    parmconstr.append(constr)
    if(p.nature == 'external'):
        parmmodelconstr.append(modelconstr)
    parmextinter.append(extinter)
    parmfuncmap.append(funcmap)
    paramsforev.append(forev)
    paramsforstream.append(forev)
    if extinter != '':
        parmextinter.append('\n')

parmtextsubs = { 'parmgetters' : '\n'.join(parmgetters),
                 'parmdecls' : '\n'.join(parmdecls),
                 'parmconstr' : ': ' + ',\n  '.join(parmconstr),
                 'getters' : '',
                 'decls' : '',
                 'addVertex' : '',
                 'ostream' : '\n\t<< '.join(paramsforstream),
                 'istream' : '\n\t>> '.join(paramsforstream),
                 'refs' : '',
                 'parmextinter': ''.join(parmextinter),
                 'num_params': len(FR.all_parameters),
                 'parmfuncmap': '\n'.join(parmfuncmap),
                 'paramsforev': ','.join(paramsforev),
                 'ModelName': ModelName
                 }



#for k,v in parmtextsubs.iteritems():
    #print k
    #print v
    #print

# write the files from templates according to the above subs
writeFile( ModelName + '.h', MODEL_H.substitute(parmtextsubs) )
writeFile( ModelName +'.cc', MODEL_CC.substitute(parmtextsubs) )
writeFile( 'LHC-' + ModelName +'.in', MODEL_HWIN.substitute(parmtextsubs) )

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



particleT = string.Template(
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


def get_lorentztag(spin):
    'Produce a ThePEG spin tag for the given numeric FR spins.'
    spins = { 1 : 'S', 2 : 'F', 3 : 'V', -1 : 'U', 5 : 'T' }
    result = [ spins[s] for s in spin ]

    def spinsort(a,b):
        "Helper function for ThePEG's FVST spin tag ordering."
        if a == b: return 0
        for letter in 'FVST':
            if a == letter: return -1
            if b == letter: return  1

    result = sorted(result, cmp=spinsort)
    return ''.join(result)




##################################################
##################################################
##################################################

# count number of vertices not included
notincluded = 0

# get vertex template
VERTEX = getTemplate('Vertex.cc')

def produce_vertex_file(subs):
    newname = ModelName + subs['classname'] + '.cc'
    writeFile( newname, VERTEX.substitute(subs) )    

# loop over all vertices
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

    if 'T' in lt:   spind = 'Tensor'
    elif 'S' in lt: spind = 'Scalar'
    elif 'V' in lt: spind = 'Vector'
    elif 'U' in lt: spind = 'Ghost'
    
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
        if(pdgcode[i] == 23):
            vhasz += 1
        if(pdgcode[i] == 22):
            vhasp += 1
        if(pdgcode[i] == 25):
            vhash += 1
        if(pdgcode[i] == 21):
            vhasg += 1
        if(pdgcode[i] == 24):
            vhasw += 1
        if(abs(pdgcode[i]) < 7 or (abs(pdgcode[i]) > 10 and abs(pdgcode[i]) < 17)):
            vhasf += 1
        if(pdgcode[i] not in SMPARTICLES):
            notsmvertex = True
        

#  treat replacement of SM vertices with BSM vertices?               
    if(notsmvertex == False):
        if( (vhasf == 2 and vhasz == 1) or (vhasf == 2 and vhasw == 1) or (vhasf == 2 and vhash == 1) or (vhasf == 2 and vhasg == 1) or (vhasf == 2 and vhasp == 0) or (vhasg == 3) or (vhasg == 4) or (vhasw == 2 and vhash == 1) or (vhasw == 3) or (vhasw == 4) or (vhash == 1 and vhasg == 2) or (vhash == 1 and vhasp == 2)):
            #print 'VERTEX INCLUDED IN STANDARD MODEL!'
            v.include = 0
            notincluded += 1
            continue
            
    
    selfconjugate = 0
    for j in range(len(pdgcode)):
        for k in range(len(pdgcode)):
               if( j != k and j != 0 and abs(pdgcode[j]) == abs(pdgcode[k])):
                   selfconjugate = 1
                   #print 'self-conjugate vertex'
#        print pdgcode[j]

# if the Vertex is not self-conjugate, then add the conjugate vertex
# automatically
    scfac = [1,1,1,1]
    if(selfconjugate == 0):
        #first find the self-conjugate particles
        for u in range(len(v.particles)):
              if(v.particles[u].selfconjugate == 0):
                  scfac[u] = -1
#                  print 'particle ', v.particles[u].pdg_code, ' found not to be self-conjugate'
                  
    if(selfconjugate == 0):
        plistarray[1] += str(scfac[1] * v.particles[1].pdg_code) + ',' + str(scfac[0] * v.particles[0].pdg_code) + ',' + str(scfac[2] * v.particles[2].pdg_code)
        if(len(v.particles) is 4):                                                                                                                      
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
        if(qed == 0 and qcd == 0):
            qed = 1
        unique_qcd( qed )
        unique_qed( qcd )
        L = v.lorentz[li]

        if lt in ['FFS','FFV']:
            #print 'PRINTING LORENTZ STRUCTURE'
            #print L.structure
            for lor in map(string.strip, L.structure.split('+')):
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

    if(plistarray[1] is ''):
        plist2 = ''
    else:
        plist2 = 'addToList(%s);' % plistarray[1]


    # input q2 or not, depending on whether it is necessary
    if('q2' in norm or 'q2' in left or 'q2' in right):
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
             'spindirectory' : spind,
             'ModelName' : ModelName,
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


    #print plistarray[0]
#    if plist in allplist:
#        print 'PLIST IN ALLPLIST'
        
        
    if( L.spins[0] != -1 and L.spins[1] != -1 and L.spins[2] != -1 and plistarray[0] not in allplist and plistarray[1] not in allplist):
        produce_vertex_file(subs)
        allplist += plistarray[0]
        allplist += plistarray[1]
    elif( L.spins[0] != -1 and L.spins[1] != -1 and L.spins[2] != -1 and selfconjugate):
        produce_vertex_file(subs)
        allplist += plistarray[0]
    else:
        #print 'VERTEX ALREADY INCLUDED'
        v.include = 0
        
        #print '============================================================'

##################################################
##################################################
##################################################

vertexline = string.Template("""\
create $classname $name
insert ${ModelName}:ExtraVertices 0 $name
""")


def get_vertices():
    vlist = 'library ' + ModelName + '.so\n'
    for v in FR.all_vertices:
        for l in v.lorentz:
            lt = get_lorentztag(l.spins)
            #print lt
        if("U" not in lt and v.include == 1):
            vlist += vertexline.substitute(
                { 'classname' : 'Herwig::' + ModelName + 'V_%03d' % int(v.name[2:]),
                'name' : '/Herwig/' + ModelName + '/%s'%v.name, 'ModelName' : ModelName } )
    return vlist


modelfilesubs = { 'plist' : get_all_thepeg_particles(),
                  'vlist' : get_vertices(),
                  'setcouplings': '\n'.join(parmmodelconstr),
                  'ModelName': ModelName
                  }

#print get_all_thepeg_particles()

MODELINFILE = getTemplate('FR.model')

writeFile( ModelName +'.model', MODELINFILE.substitute(modelfilesubs) )

print 'finished generating model', ModelName
print 'model directory:', opts.MODELDIR
print 'generated', len(FR.all_vertices)-notincluded, 'vertices'

