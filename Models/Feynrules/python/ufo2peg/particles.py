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
setup $name $pdg_code $name $mass $width $wcut $ctau $charge $color $spin 0
insert /Herwig/NewPhysics/NewModel:DecayParticles 0 $name
"""
)


class ParticleConverter:
    'Convert a FR particle to extract the information ThePEG needs.'
    def __init__(self,p,parmsubs):
        self.name = p.name
        self.pdg_code = p.pdg_code
        self.spin = p.spin
        self.color = p.color
        if self.color == 1:
            self.color = 0
        self.selfconjugate = 0
        self.mass = parmsubs[str(p.mass)]
        self.width = parmsubs[str(p.width)]
        try:
            self.mass = self.mass.real
        except:
            pass
        self.mass = abs(self.mass)
        hbarc = 197.3269631e-15 # GeV mm (I hope ;-) )
        self.ctau = (hbarc / self.width) if self.width != 0 else 0
        self.wcut = 10 * self.width
        self.charge = int(3 * p.charge)

    def subs(self):
        return self.__dict__



def thepeg_particles(FR,parameters):
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
            subs = ParticleConverter(p,parameters).subs()
            plist += particleT.substitute(subs)

            pdg, name = subs['pdg_code'],  subs['name']
            if -pdg in antis:
                plist += 'makeanti %s %s\n' % (antis[-pdg], name)
                
            else:
                antis[pdg] = name
                selfconjugate = 1
    return plist
