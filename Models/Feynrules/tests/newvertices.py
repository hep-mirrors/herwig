import Model as FR
import cmath

def getTemplate(basename):
    import string
    f = open('../%s.template' % basename, 'r')
    templateText = f.read()
    f.close()
    return string.Template( templateText )

def writeFile(filename, text):
    f = open(filename,'w')
    f.write(text)
    f.close()



##################################################
##################################################
##################################################



MODEL_H  = getTemplate('Model.h')
MODEL_CC = getTemplate('Model.cc')

parmdecls = []
parmgetters = []
parmconstr = []

parmsubs = dict( [ (p.name, float(p.value)) 
                   for p in FR.all_parameters 
                   if p.nature == 'external' ] ) 
parmsubs.update( { 'ZERO' : 0 } )

internal = [ p for p in FR.all_parameters 
             if p.nature == 'internal' ] 

for p in internal:
    newval = eval(p.value, { 'cmath' : cmath }, parmsubs)
    parmsubs.update( { p.name : newval } )

for p in FR.all_parameters:
    value = parmsubs[p.name]
    if p.type == 'real':
        try:
            assert( value.imag < 1.0e-16 )
            value = value.real
        except:
            pass
        decl = '  double %s_;' % p.name
        constr = '%s_(%s)' % (p.name, value)
        getter = '  double %s() { return %s_; }' % (p.name, p.name)
    elif p.type == 'complex':
        decl = '  Complex %s_;' % p.name
        constr = '%s_(%s,%s)' % (p.name, value.real, value.imag)
        getter = '  Complex %s() { return %s_; }' % (p.name, p.name)
    else:
        raise Exception('Unknown data type "%s".' % p.type)

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

writeFile( 'FeynRulesModel.h', MODEL_H.substitute(parmtextsubs) )
writeFile( 'FeynRulesModel.cc', MODEL_CC.substitute(parmtextsubs) )


#exit(0)
##################################################
##################################################
##################################################
import string



particleT = string.Template(
"""
create ThePEG::ParticleData $name
setup $name $pdg_code $name $mass $width $wcut $ctau $charge $color $spin 0
"""
)
class ParticleConverter:
    def __init__(self,p):
        self.name = p.name
        self.pdg_code = p.pdg_code
        self.spin = p.spin
        self.color = p.color
        self.mass = parmsubs[p.mass]
        self.width = parmsubs[p.width]
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

def get_table():
    plist = ''
    antis = {}
    for p in FR.all_particles:
        if p.spin == -1 or p.goldstoneboson:
            continue

        subs = ParticleConverter(p).subs()
        plist += particleT.substitute(subs)

        pdg, name = subs['pdg_code'],  subs['name']
        if -pdg in antis:
            plist += 'makeanti %s %s\n' % (antis[-pdg], name)
        else:
            antis[pdg] = name

    return plist

modelfilesubs = { 'plist' : get_table() }

print get_table()

MODELINFILE = getTemplate('FR.model')

writeFile( 'FR.model', MODELINFILE.substitute(modelfilesubs) )

##################################################
##################################################
##################################################


VERTEX = getTemplate('Vertex.cc')

def produce_vertex_file(subs):
    newname = 'FR' + subs['classname'] + '.cc'
    writeFile( newname, VERTEX.substitute(subs) )

def get_lorentztag(spin):
    spins = { 1 : 'S', 2 : 'F', 3 : 'V', 5 : 'T' }
    result = [ spins[s] for s in spin ]

    def spinsort(a,b):
        if a == b: return 0
        for letter in 'FVST':
            if a == letter: return -1
            if b == letter: return  1

    result = sorted(result, cmp=spinsort)
    return ''.join(result)



for v in FR.all_vertices:

    print v.name
    print map(str,v.particles)
    print '---------------'
    for (col,lor),C in v.couplings.iteritems():
        L = v.lorentz[lor]
        print 'Colour  :',v.color[col]
        print 'Lorentz :',L.name, L.spins, L.structure
        print 'Coupling:',C.name, C.value, C.order
        print '---------------'
    print '============================================================'








    ### Spin structure
    lt = None
    for l in v.lorentz:
        newLt = get_lorentztag(l.spins)
        if lt is None:
            lt = newLt
        else:
            # multiple lorentz structures must still refer to the same
            # spin content
            assert( newLt == lt )

    if 'T' in lt:   spind = 'Tensor'
    elif 'S' in lt: spind = 'Scalar'
    else:           spind = 'Vector'

    ### Particle ids #################### sort order? ####################
    plist = ','.join([ str(p.pdg_code) for p in v.particles ])
    

    ### Colour structure
    if v.color == '1': qcdord = '0'
    else:              qcdord = ''

    ### classname
    classname = 'V_%03d' % int(v.name[2:])

    ### parse couplings


    leftcontent = 1.
    rightcontent = 1.

    ### do we need left/right?
    if 'FF' in lt:
        left  = 'left(%s);'  % leftcontent
        right = 'right(%s);' % rightcontent
    else:
        left = ''
        right = ''

    ### assemble dictionary and fill template
    subs = { 'lorentztag' : lt,                   # ok
             'classname'  : classname,            # ok
             'left'       : left,                 # doesn't always exist in base
             'right'      : right,                 # doesn't always exist in base 
              'norm'      : '1.',                 # needs norm, too

             #################### need survey which different setter methods exist in base classes

             'addToPlist' : 'addToList(%s);' % plist, # ok
             'parameters' : '{}',
             'setCouplings' : '{}',
             'qedorder'   : '',
             'qcdorder' : qcdord,
             'q2'        :  '',
             'couplingptrs' : ',tcPDPtr'*len(v.particles),
             'spindirectory' : spind}             # ok
    
    produce_vertex_file(subs)

print len(FR.all_vertices)


