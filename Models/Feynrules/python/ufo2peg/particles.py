from __future__ import print_function
from string import Template
import os
import numpy as np

# ignore these, they're in Hw++ already # TODO reset Hw++ settings instead
SMPARTICLES = {

1:'d',
2:'u',
3:'s',
4:'c',
5:'b',
6:'t', # think about this one later

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

-24:'W-'
}



particleT = Template(
"""
create ThePEG::ParticleData $name
# values set to 999999 are recalculated later from other model parameters
setup $name $pdg_code $name $mass $width $wcut $ctau $charge $color $spin 0
"""
)


class ParticleConverter:
    'Convert a FR particle to extract the information ThePEG needs.'
    def __init__(self,p,parmsubs,modelparameters):
        self.name = p.name
        self.pdg_code = p.pdg_code
        self.spin = p.spin
        self.color = p.color
        if self.color == 1:
            self.color = 0
        self.selfconjugate = 0

        self.mass = parmsubs[str(p.mass)]
        if type(self.mass) == str:
            value = modelparameters[self.mass]
            try:
                value = value.real
            except:
                pass
            newname = '%s_ABS' % self.mass
            self.mass = '${%s}' % newname
            modelparameters[newname] = abs(value)
        else:
            try:
                self.mass = self.mass.real
            except:
                pass
            self.mass = 999999. # abs(self.mass)

        hbarc = 197.3269631e-15 # GeV mm (I hope ;-) )
        self.width = parmsubs[str(p.width)]
        if type(self.width) == str:
            width = modelparameters[self.width]
            ctau = (hbarc / width) if width != 0 else 0
            newname = '%s_CTAU' % self.width
            self.ctau = '${%s}' % newname
            modelparameters[newname] = ctau

            wcut = 10 * width
            newname = '%s_WCUT' % self.width
            self.wcut = '${%s}' % newname
            modelparameters[newname] = wcut

            self.width = '${%s}' % self.width
        else:
            self.ctau = 999999. # (hbarc / self.width) if self.width != 0 else 0
            self.wcut = 999999. #10.0 * self.width
            self.width = 999999. # was blank line before

        self.charge = int(3 * p.charge)

    def subs(self):
        return self.__dict__

def check_effective_vertex(FR,p,ig) :
    for vertex in FR.all_vertices:
        if(len(vertex.particles) != 3) : continue
        if(p not in vertex.particles ) : continue
        ng=0
        for part in vertex.particles :
            if(part.pdg_code==ig) : ng+=1
        if(ng==2) :
            return False
    return True

# finds all dim-Dimention vectices in FR that involve pIn
def sort_vertices(FR,pIn,dim):
    possibleVertices = []
    for V in FR.all_vertices:
        # only keep 1 -> dim-1 vetrices
        if len(V.particles) != dim :
            continue
        # keep vertices with pIn
        if pIn == str(V.particles[0]) or pIn == str(V.particles[1]) or pIn == str(V.particles[2]):
            possibleVertices.append(V)
    return possibleVertices

def load_parameters(FR):
    if os.path.isfile("FR_Parameters.py") :
        return True
    else :
        file = open("FR_Parameters.py", "a")
        header = \
"""
import cmath
all_functions = []
class Function(object):
    def __init__(self, name, arguments, expression):
        global all_functions
        all_functions.append(self)

        self.name = name
        self.arguments = arguments
        self.expr = expression
    def __call__(self, *opt):
        for i, arg in enumerate(self.arguments):
            exec('%s = %s' % (arg, opt[i] ))
        return eval(self.expr)
complexconjugate = Function(name = 'complexconjugate',
                            arguments = ('z',),
                            expression = 'z.conjugate()')
re = Function(name = 're',
              arguments = ('z',),
              expression = 'z.real')
im = Function(name = 'im',
              arguments = ('z',),
              expression = 'z.imag')
sec = Function(name = 'sec',
             arguments = ('z',),
             expression = '1./cmath.cos(z.real)')
asec = Function(name = 'asec',
             arguments = ('z',),
             expression = 'cmath.acos(1./(z.real))')
csc = Function(name = 'csc',
             arguments = ('z',),
             expression = '1./cmath.sin(z.real)')
acsc = Function(name = 'acsc',
             arguments = ('z',),
             expression = 'cmath.asin(1./(z.real))')

cot = Function(name = 'cot',
               arguments = ('z',),
               expression = '1./cmath.tan(z.real)')
theta_function = Function(name = 'theta_function',
             arguments = ('x','y','z'),
             expression = 'y if x else z')
cond = Function(name = 'cond',
                arguments = ('condition','ExprTrue','ExprFalse'),
                expression = '(ExprTrue if condition==0.0 else ExprFalse)')
reglog = Function(name = 'reglog',
                arguments = ('z'),
                expression = '(0.0 if z==0.0 else cmath.log(z.real))')
"""
        file.write(header)
        for par in FR.all_parameters :
            name  = par.name
            value = par.value
            file.write(name+" = "+str(value)+"\n")
        file.close()
        return True

# extracts all possible splittings for incoming particle p
def sort_splittings(FR,Vertices,p):
    pSplittings = []
    for V in Vertices:
        # calculate the total coupling value for this splitting
        if load_parameters(FR) :
            from FR_Parameters import *
        coup = V.couplings
        keys = coup.keys()
        coupling_value = 0.
        for k in keys :
            coupling_value += eval(coup[k].value)
        # extract splitting format as p -> p1, p2
        p0set = False
        p1set = False
        p2set = False
        for particle in V.particles:
            if particle == p and not p0set:
                p0set = True
            elif not p1set :
                p1 = particle
                p1set = True
            elif not p2set :
                p2 = particle
                p2set = True
        if not p0set or not p1set or not p2set :
            continue
        id1 = abs(p1.pdg_code)
        id2 = abs(p2.pdg_code)
        # TODO need to improve this forbidden list assortment
        forbidden = [250, 251, 9000001, 9000002, 9000003, 9000004]
        if id1 in forbidden or id2 in forbidden:
            continue
        # put the bigger spin last
        if p1.spin > p2.spin :
            p1, p2 = p2, p1
        pp1p2 = [p,p1,p2,coupling_value]
        pp2p1 = [p,p2,p1,coupling_value]
        if pp1p2 not in pSplittings and pp2p1 not in pSplittings:
            pSplittings.append(pp1p2)
    return pSplittings

def extract_mass(FR,Vertex) :
    if load_parameters(FR) :
        from FR_Parameters import *
    m = [0.,0.,0.]
    m[0] = Vertex[0].mass.value
    m[1] = Vertex[1].mass.value
    m[2] = Vertex[2].mass.value
    for i in range(len(m)) :
        if isinstance(m[i],str) :
            m[i] = eval(m[i])
    return m

def split_name(Vertex,split=False) :
    p = [Vertex[0].name, Vertex[1].name, Vertex[2].name]
    for i in range(0,3) :
        if Vertex[i].pdg_code in SMPARTICLES :
            p[i] = SMPARTICLES[Vertex[i].pdg_code]
        else :
            p[i] = Vertex[i].name
    if not split :
        splitname = p[0] + p[1] + p[2]
        splitname = splitname.replace("+", "p")
        splitname = splitname.replace("-", "m")
        return splitname
    else :
        return p[0], p[1], p[2]

def isQuark(particle) :
    if abs(particle.pdg_code) >= 1 and abs(particle.pdg_code) <= 6 :
        return True
    else :
        return False

def isLepton(particle) :
    if abs(particle.pdg_code) >= 11 and abs(particle.pdg_code) <= 16 :
        return True
    else :
        return False

def isScalar(particle) :
    if particle.spin==1 :
        return True
    else :
        return False

def isGVB(particle) :
    if abs(particle.pdg_code)==24 :
        return True
    elif particle.pdg_code==23 :
        return True
    elif particle.pdg_code==22 :
        return True
    elif particle.pdg_code==21 :
        return True
    else :
        return False

def antiparticle(FR,particle) :
    for anti in FR.all_particles :
        if particle.pdg_code == -anti.pdg_code :
            return anti
    return particle

def thepeg_particles(FR,parameters,modelname,modelparameters,forbidden_names,hw_higgs,allow_fcnc) :
    plist = []
    antis = {}
    names = []
    splittings = []
    done_splitting_QCD = []
    done_splitting_QED = []
    done_splitting_EW  = []

    for p in FR.all_particles:
        if p.spin == -1:
            continue

        gsnames = ['goldstone',
                   'goldstoneboson',
                   'GoldstoneBoson']

        def gstest(name):
            try:
                return getattr(p,name)
            except AttributeError:
                return False

        if any(map(gstest, gsnames)):
            continue

        if p.pdg_code in SMPARTICLES:
            continue

        if p.pdg_code == 25 and not hw_higgs:
            plist.append(
"""
set /Herwig/Particles/h0:Mass_generator NULL
set /Herwig/Particles/h0:Width_generator NULL
rm /Herwig/Masses/HiggsMass
rm /Herwig/Widths/hWidth
"""
)
        if p.name in forbidden_names:
            print('RENAMING PARTICLE',p.name,'as ',p.name+'_UFO')
            p.name +="_UFO"
        subs = ParticleConverter(p,parameters,modelparameters).subs()
        if not (p.pdg_code == 25 and hw_higgs) :
            plist.append( particleT.substitute(subs) )

        pdg, name = subs['pdg_code'],  subs['name']
        names.append(name)
        if -pdg in antis:
            plist.append( 'makeanti %s %s\n' % (antis[-pdg], name) )

        elif not (p.pdg_code == 25 and hw_higgs) :
            plist.append( 'insert /Herwig/NewPhysics/NewModel:DecayParticles 0 %s\n' % name )
            plist.append( 'insert /Herwig/Shower/ShowerHandler:DecayInShower 0 %s #  %s' % (abs(pdg), name) )
            antis[pdg] = name
            selfconjugate = 1

        class SkipMe(Exception):
            pass

        def spin_name(s):
            spins = { 1 : 'Zero',
                      2 : 'Half',
                      3 : 'One' }
            if s not in spins:
                raise SkipMe()
            else:
                return spins[s]

        def col_name(c):
            cols = { 3 : 'Triplet',
                     6 : 'Sextet',
                     8 : 'Octet' }
            return cols[c]

        try:
            # QCD splitting functions
            if p.color in [3,6,8] and abs(pdg) not in done_splitting_QCD: # which colors?
                done_splitting_QCD.append(abs(pdg))
                splitname = '{name}SplitFnQCD'.format(name=p.name)
                sudname = '{name}SudakovQCD'.format(name=p.name)
                splittings.append(
"""
create Herwig::{s}{s}OneSplitFn {name}
set {name}:InteractionType QCD
set {name}:ColourStructure {c}{c}Octet
cp /Herwig/Shower/SudakovCommon {sudname}
set {sudname}:SplittingFunction {name}
do /Herwig/Shower/SplittingGenerator:AddFinalSplitting {pname}->{pname},g; {sudname}
""".format(s=spin_name(p.spin), name=splitname,
           c=col_name(p.color), pname=p.name, sudname=sudname)
                )
        except SkipMe:
            pass
        # QED splitting functions
        try:
            if p.charge != 0 and abs(pdg) not in done_splitting_QED:
                done_splitting_QED.append(abs(pdg))
                splitname = '{name}SplitFnQED'.format(name=p.name)
                sudname = '{name}SudakovQED'.format(name=p.name)
                splittings.append(
"""
create Herwig::{s}{s}OneSplitFn {name}
set {name}:InteractionType QED
set {name}:ColourStructure ChargedChargedNeutral
cp /Herwig/Shower/SudakovCommon {sudname}
set {sudname}:SplittingFunction {name}
set {sudname}:Alpha /Herwig/Shower/AlphaQED
do /Herwig/Shower/SplittingGenerator:AddFinalSplitting {pname}->{pname},gamma; {sudname}
""".format(s=spin_name(p.spin), name=splitname, pname=p.name, sudname=sudname)
                )
        except SkipMe:
            pass

        # EW and BSM splitting functions
        try:
            Vertices = sort_vertices(FR,str(p.name),3)
            pSplittings = sort_splittings(FR,Vertices,p)
            for Vertex in pSplittings :
                # do not do anything for CouplingValue < 1e-6
                if abs(Vertex[3].real + Vertex[3].imag) < 1e-6 :
                    continue
                # do not include QCD splittings
                if Vertex[1].pdg_code == 21 or Vertex[2].pdg_code == 21 :
                    continue
                # do not include QED splittings
                if Vertex[2].pdg_code == 22 and abs(Vertex[0].pdg_code) == abs(Vertex[1].pdg_code) or\
                   Vertex[1].pdg_code == 22 and abs(Vertex[0].pdg_code) == abs(Vertex[2].pdg_code):
                    continue
                # skip lepton vertices
                if isLepton(Vertex[1]) or isLepton(Vertex[2]) :
                    continue
                # loop over all possible configurations in the splitting
                for pos in range(0,3) :
                    # rearrange to all possible cases
                    V=[Vertex[0], Vertex[1], Vertex[2], Vertex[3]]
                    if pos==0:
                        V[0], V[1], V[2] = Vertex[0], Vertex[1], Vertex[2]
                    elif pos==1:
                        V[0], V[1], V[2] = Vertex[1], Vertex[2], Vertex[0]
                    else :
                        V[0], V[1], V[2] = Vertex[2], Vertex[0], Vertex[1]
                    # don't allow photon as progenitor
                    if V[0].pdg_code == 22 :
                        continue
                    # for a generic splitting m0 < m1+m2, otherwise it's a decay
                    m = extract_mass(FR,V)
                    # filter out m+i*zero
                    for ix in range(len(m)) :
                        if isinstance(m[ix], complex) :
                            m[ix] = m[ix].real
                    if m[0] > m[1] + m[2] :
                        continue

                    """
                    possible EW splittings:
                        - V > V  H   with SM GVBs and BSM (charged and nuetral) Scalars
                        - q > q' H   with SM quarks and BSM (charged and nuetral) Scalars, may include FCNC-inducing splittings
                        - H > H' V   with BSM (charged and nuetral) Scalars and SM GVBs (W,Z,photon)
                        - V > H  H'  with SM GVBs and BSM (charged and nuetral) Scalars
                        - H > H' H'' Higgs splittings
                    """

                    if isGVB(V[0]) :
                        # allow V > V H
                        if isGVB(V[1]) and isScalar(V[2]) :
                            pass
                        elif isGVB(V[2]) and isScalar(V[1]) :
                            V[0], V[1], V[2] = V[0], V[2], V[1]
                        # allow V > H H'
                        elif isScalar(V[1]) and isScalar(V[2]):
                            pass
                        # nothing else with a GVB progenitor
                        else :
                            continue
                    elif isQuark(V[0]) :
                        # allow q > q' H (including FCNC splittings)
                        if isQuark(V[1]) and isScalar(V[2]) :
                            pass
                        elif isQuark(V[2]) and isScalar(V[1]) :
                            V[0], V[1], V[2] = V[0], V[2], V[1]
                        # nothing else with a quark progenitor
                        else :
                            continue
                    elif isScalar(V[0]) :
                        # allow H > H' H''
                        if isScalar(V[1]) and isScalar(V[2]) :
                            pass
                        # allow H > H' V
                        elif isScalar(V[1]) and isGVB(V[2]) :
                            pass
                        elif isScalar(V[2]) and isGVB(V[1]) :
                            V[0], V[1], V[2] = V[0], V[2], V[1]
                        # nothing else with a scalar progenitor
                        else :
                            continue
                    # nothing else
                    else :
                        print(V[0]," > ",V[1], V[2], "splitting cannot be handled.")
                        continue

                    # getting the electric charge right
                    if V[0].charge != V[1].charge+V[2].charge :
                        if V[0].pdg_code*V[1].pdg_code < 0 :
                            if V[0].pdg_code > 0 :
                                V[1] = antiparticle(FR,V[1])
                            else :
                                V[0] = antiparticle(FR,V[0])
                        elif V[0].pdg_code*V[2].pdg_code < 0 :
                            if V[0].pdg_code > 0 :
                                V[2] = antiparticle(FR,V[2])
                            else :
                                V[0] = antiparticle(FR,V[0])
                    if V[0].charge != V[1].charge+V[2].charge :
                        if abs(V[0].charge-V[1].charge-V[2].charge) < 1e-10 :
                            pass
                        else :
                            continue

                    # deal with FCNC inducing splittings
                    if not allow_fcnc :
                        if isQuark(V[0]) and isQuark(V[1]) and isScalar(V[2]) and\
                           V[2].charge==0. and V[0].pdg_code!=V[1].pdg_code :
                            print("Omitting",V[0],"->",V[1],",",V[2], "FCNC-inducing splitting.",\
                                 "Use --allow-fcnc flag if you wish to keep this.")
                            continue

                    # set up the splitting
                    SPname = split_name(V,False)
                    # no scalar > quark, antiquark splitting
                    s = [spin_name(V[0].spin),spin_name(V[1].spin),spin_name(V[2].spin)]

                    if SPname not in done_splitting_EW and V[3].imag!=0.:
                        done_splitting_EW.append(SPname)
                        splitname = '{name}SplitFnEW'.format(name=SPname)
                        sudname = '{name}SudakovEW'.format(name=SPname)
                        p0name, p1name, p2name = split_name(V,True)
                        splittings.append(
"""
create Herwig::{s0}{s1}{s2}EWSplitFn {name}
set {name}:InteractionType EW
set {name}:ColourStructure EW
set {name}:CouplingValue {i}
cp /Herwig/Shower/SudakovCommon {sudname}
set {sudname}:SplittingFunction {name}
set {sudname}:Alpha /Herwig/Shower/AlphaEW
do /Herwig/Shower/SplittingGenerator:AddFinalSplitting {p0}->{p1},{p2}; {sudname}
""".format(s0=s[0],s1=s[1],s2=s[2],name=splitname,p0=p0name,\
          p1=p1name,p2=p2name,sudname=sudname,i=V[3].imag)
                    )
        except SkipMe:
            pass


        if p.charge == 0 and p.color == 1 and p.spin == 1 and not (p.pdg_code == 25 and hw_higgs) :
            if(check_effective_vertex(FR,p,21)) :
                plist.append(
"""
insert /Herwig/{ModelName}/V_GenericHGG:Bosons 0 {pname}
""".format(pname=p.name, ModelName=modelname)
                )
            if(check_effective_vertex(FR,p,22)) :
                plist.append(
"""
insert /Herwig/{ModelName}/V_GenericHPP:Bosons 0 {pname}
""".format(pname=p.name, ModelName=modelname)
                )

    return ''.join(plist)+''.join(splittings), names
