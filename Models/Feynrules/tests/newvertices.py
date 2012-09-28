#! /usr/bin/env python
from __future__ import with_statement
import Model as FR
import cmath, string

def getTemplate(basename):
    with open('../%s.template' % basename, 'r') as f:
        templateText = f.read()
    return string.Template( templateText )

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

MODEL_H  = getTemplate('Model.h')
MODEL_CC = getTemplate('Model.cc')

allplist = ""

parmdecls = []
parmgetters = []
parmconstr = []

parmsubs = dict( [ (p.name, float(p.value)) 
                   for p in FR.all_parameters 
                   if p.nature == 'external' ] ) 


print parmsubs
print


def evaluate(x):
    return eval(x, 
                {'cmath':cmath,
                 'complexconjugate':FR.function_library.complexconjugate}, 
                parmsubs)


internal = [ p 
             for p in FR.all_parameters 
             if p.nature == 'internal' ] 

print internal
print

for p in internal:
    print p.name,'=',p.value
    newval = evaluate(p.value)
    parmsubs.update( { p.name : newval } )

print parmsubs
print

for p in FR.all_parameters:
    value = parmsubs[p.name]
    if p.type == 'real':
        try:
            assert( value.imag < 1.0e-16 )
            value = value.real
        except:
            pass
        parmsubs[p.name] = value
        decl = '  double %s_;' % p.name
        constr = '%s_(%s)' % (p.name, value)
        getter = '  double %s() { return %s_; }' % (p.name, p.name)
    elif p.type == 'complex':
        value = complex(value)
        parmsubs[p.name] = value
        decl = '  Complex %s_;' % p.name
        constr = '%s_(%s,%s)' % (p.name, value.real, value.imag)
        getter = '  Complex %s() { return %s_; }' % (p.name, p.name)
    else:
        raise Exception('Unknown data type "%s".' % p.type)

    # do calc in C++, add interfaces for externals

    parmdecls.append(decl)
    parmgetters.append(getter)
    parmconstr.append(constr)

parmtextsubs = { 'parmgetters' : '\n'.join(parmgetters),
                 'parmdecls' : '\n'.join(parmdecls),
                 'parmconstr' : ': ' + ',\n  '.join(parmconstr),
                 'getters' : '',
                 'decls' : '',
                 'addVertex' : '',
                 'ostream' : '',
                 'istream' : '',
                 'refs' : ''
                 }
for k,v in parmtextsubs.iteritems():
    print k
    print v
    print

writeFile( 'FeynRulesModel.h', MODEL_H.substitute(parmtextsubs) )
writeFile( 'FeynRulesModel.cc', MODEL_CC.substitute(parmtextsubs) )


#exit(0)
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


VERTEX = getTemplate('Vertex.cc')

def produce_vertex_file(subs):
    newname = 'FR' + subs['classname'] + '.cc'
    writeFile( newname, VERTEX.substitute(subs) )





for v in FR.all_vertices:

    print v.name
    print map(str,v.particles)
    print '---------------'
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
    print plist

# Check if the Vertex is self-conjugate or not
    pdgcode = [0,0,0,0]
#    print 'printing particles in vertex'
    for i in range(len(v.particles)):
#       print v.particles[i].pdg_code
        pdgcode[i] = v.particles[i].pdg_code

    selfconjugate = 0
    for j in range(len(pdgcode)):
        for k in range(len(pdgcode)):
               if( j != k and j != 0 and abs(pdgcode[j]) == abs(pdgcode[k])):
                   selfconjugate = 1
                   print 'self-conjugate vertex'
#        print pdgcode[j]

# if the Vertex is not self-conjugate, then add the conjugate vertex
# WARNING:
# TO DO: What if the coupling is complex? Need to add the complex conjugate of that coupling in that case
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
        print 'Conjugate vertex:', plistarray[1]
    
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
        if(qed == 0 and qcd == 0): qed = 1
        unique_qcd( qed )
        unique_qed( qcd )
        L = v.lorentz[li]

        if lt in ['FFS','FFV']:
            print L.structure
            for lor in map(string.strip, L.structure.split('+')):
                breakdown = lor.split('*')
                prefactor='1'
                if len(breakdown) == 3:
                    prefactor = breakdown[0]
                    breakdown = breakdown[1:]

                if len(breakdown) == 2:
                    assert(breakdown[0][:5] == 'Gamma')
                    if breakdown[1][:5] == 'ProjM':
                        coup_left.append(prefactor+' * '+C.value)
                    elif breakdown[1][:5] == 'ProjP':
                        coup_right.append(prefactor+' * '+C.value)
                else:
                    coup_left.append(C.value)
                    coup_right.append(C.value)
        else:
            coup_norm.append(C.value)
                

        print 'Colour  :',v.color[ci]
        print 'Lorentz %s:'%L.name, L.spins, L.structure
        print 'Coupling %s:'%C.name, C.value, '\nQED=%s'%qed, 'QCD=%s'%qcd
        print '---------------'


    leftcontent = ' + '.join(coup_left) if len(coup_left)!=0 else '0j'
    rightcontent = ' + '.join(coup_right) if len(coup_right)!=0 else '0j'
    normcontent = ' + '.join(coup_norm) if len(coup_norm)!=0 else '1.'

    print 'Left:',leftcontent
    print 'Right:',rightcontent
    print 'Norm:',normcontent
    print '---------------'
  
    leftcontent = complex(evaluate(leftcontent))
    rightcontent = complex(evaluate(rightcontent))
    normcontent = complex(evaluate(normcontent))

    print 'Left:',leftcontent
    print 'Right:',rightcontent
    print 'Norm:',normcontent

    ### do we need left/right?
    if 'FF' in lt:
        left  = 'left(Complex(%s,%s));'  % (leftcontent.real,leftcontent.imag)
        right = 'right(Complex(%s,%s));' % (rightcontent.real,rightcontent.imag)
    else:
        left = ''
        right = ''

    if(plistarray[1] is ''):
        plist2 = ''
    else:
        plist2 = 'addToList(%s);' % plistarray[1]
        
    norm = 'norm(Complex(%s,%s));' % (normcontent.real,normcontent.imag)


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
             'q2'        :  '',
             'couplingptrs' : ',tcPDPtr'*len(v.particles),
             'spindirectory' : spind}             # ok


    print plistarray[0]
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
        print 'VERTEX ALREADY INCLUDED'
        v.include = 0
        
    print '============================================================'

##################################################
##################################################
##################################################

vertexline = string.Template("""\
create $classname $name
insert FRModel:ExtraVertices 0 $name
""")


def get_vertices():
    vlist = 'library FeynrulesModel.so\n'
    for v in FR.all_vertices:
        for l in v.lorentz:
            lt = get_lorentztag(l.spins)
            print lt
        if("U" not in lt and v.include == 1):
            vlist += vertexline.substitute(
                { 'classname' : 'Herwig::FRV_%03d' % int(v.name[2:]),
                  'name' : '/Herwig/Feynrules/%s'%v.name } )
    return vlist


modelfilesubs = { 'plist' : get_all_thepeg_particles(),
                  'vlist' : get_vertices() }

print get_all_thepeg_particles()

MODELINFILE = getTemplate('FR.model')

writeFile( 'FR.model', MODELINFILE.substitute(modelfilesubs) )

print len(FR.all_vertices)
