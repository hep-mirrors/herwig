from string import Template

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

-24:'W-',

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



def thepeg_particles(FR,parameters,modelname,modelparameters,forbidden_names):
    plist = []
    antis = {}
    names = []
    splittings = []
    done_splitting=[]
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

        if p.pdg_code == 25:
            plist.append(
"""
set /Herwig/Particles/h0:Mass_generator NULL
set /Herwig/Particles/h0:Width_generator NULL
rm /Herwig/Masses/HiggsMass
rm /Herwig/Widths/hWidth
"""
)
        if p.name in forbidden_names:
            print 'RENAMING PARTICLE',p.name,'as ',p.name+'_UFO'
            p.name +="_UFO"

        subs = ParticleConverter(p,parameters,modelparameters).subs()
        plist.append( particleT.substitute(subs) )

        pdg, name = subs['pdg_code'],  subs['name']
        names.append(name)
        if -pdg in antis:
            plist.append( 'makeanti %s %s\n' % (antis[-pdg], name) )

        else:
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
            if p.color in [3,6,8] and abs(pdg) not in done_splitting: # which colors?
                done_splitting.append(abs(pdg))
                splitname = '{name}SplitFn'.format(name=p.name)
                sudname = '{name}Sudakov'.format(name=p.name)
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

        if p.charge == 0 and p.color == 1 and p.spin == 1:
            plist.append(
"""
insert /Herwig/{ModelName}/V_GenericHPP:Bosons 0 {pname}
insert /Herwig/{ModelName}/V_GenericHGG:Bosons 0 {pname}
""".format(pname=p.name, ModelName=modelname)
            )
    return ''.join(plist)+''.join(splittings), names
