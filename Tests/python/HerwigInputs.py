# select the template to load
# collider
KNOWN_COLLIDERS = [
    "EE-Gamma",
    "BFactory",
    "EE",
    "DIS",
    "TVT",
    "LHC-GammaGamma",
    "LHC",
    "ISR",
    "Fermilab",
    "SppS",
    "Star",
    "SPS",
    "GammaGamma",
]

def identifyCollider(name) :
    """ Work out name of collider
    """
    collider = ""
    for cand_collider in KNOWN_COLLIDERS:
        if cand_collider in name:
            collider = cand_collider
            break
    del cand_collider
    if "EHS" in name : collider="SPS"
    assert collider
    have_hadronic_collider = collider in ["TVT","LHC","ISR","SppS","Star","SPS","Fermilab"]
    return (collider,have_hadronic_collider)

class StringBuilder(object):
    """
    Avoid expensive string additions until the end
    by building up a list first.

    This helper class avoids rewriting all the += lower down
    to list operations.
    """
    def __init__(self, init = None):
        self.lines = [] if init is None else [init]

    def __iadd__(self, line):
        self.lines.append(line)
        return self

    def __str__(self):
        return '\n'.join(self.lines)

def identifySimulation(name,collider,have_hadronic_collider) :
    """ Work out the type of simulation

    Identify the parton shower and source of the matrix elements
    """
    parameters = { 
        'shower'  : '',
        'bscheme' : '',
    }
    simulation=""
    # istart determines how many name parts need to be skipped
    istart = 1
    # Dipole shower with Matchbox Powheg
    if "Dipole-Matchbox-Powheg" in name :
        istart = 4
        simulation="Matchbox"
        parameters["shower"]  = "read Matchbox/Powheg-DipoleShower.in\n"
        
        # Dipole shower with internal Powheg - Todo: Finish modifying template files.
        '''
        elif "Dipole-Powheg" in name :
        istart = 3
        simulation="Powheg"
        parameters["shower"]  = "set /Herwig/EventHandlers/EventHandler:CascadeHandler /Herwig/DipoleShower/DipoleShowerHandler\nread snippets/Dipole_AutoTunes_gss.in\n"
        '''
        
    # Dipole shower with MCatNLO
    elif "Dipole-MCatNLO" in name :
        istart = 3
        simulation="Matchbox"
        parameters["shower"]  = "read Matchbox/MCatNLO-DipoleShower.in\n" 

    # Dipole shower with Matchbox LO
    elif "Dipole-Matchbox-LO" in name :
        istart = 4
        simulation="Matchbox"
        parameters["shower"]  = "read Matchbox/LO-DipoleShower.in\n" 

    # Dipole shower with internal LO
    elif "Dipole" in name :
        istart = 2
        simulation=""
        parameters["shower"]  = "set /Herwig/EventHandlers/EventHandler:CascadeHandler /Herwig/DipoleShower/DipoleShowerHandler\nread snippets/Dipole_AutoTunes_gss.in\n"
    
    # AO shower with Matchbox Powheg
    elif "Matchbox-Powheg" in name :
        istart = 3
        simulation="Matchbox"
        parameters["shower"] = "read Matchbox/Powheg-DefaultShower.in\n"

    # AO shower with MCatNLO
    elif "Matchbox" in name :
        istart = 2
        simulation="Matchbox"
        parameters["shower"] = "read Matchbox/MCatNLO-DefaultShower.in\n"

    # AO shower with internal Powheg    
    elif "Powheg" in name :
        istart = 2
        simulation="Powheg"

    # Dipole shower with merging    
    elif "Merging" in name :
        istart = 2
        simulation="Merging"
        thefactory="MergingFactory"
    
    # Flavour settings for Matchbox    
    if simulation=="Matchbox" :
        parameters["bscheme"] = "read Matchbox/FiveFlavourScheme.in\n"
    
        if "Dipole" in parameters["shower"] :
            parameters["bscheme"] += "read Matchbox/FiveFlavourNoBMassScheme.in\n"
        
        if collider not in ['DIS','EE'] :
            parameters["nlo"] = "read Matchbox/MadGraph-OpenLoops.in\n"

    # Flavour settings for dipole shower with internal ME
    if simulation=="" and "Dipole" in parameters["shower"] :
        parameters["bscheme"] = "read snippets/DipoleShowerFiveFlavours.in"
    
    # find the template
    if simulation=="" :
        if collider=="LHC-GammaGamma" :
            istart += 1
            templateName="Hadron-Gamma.in"
        elif have_hadronic_collider :
            templateName="Hadron.in"
        elif collider=="EE-Gamma" :
            istart+=1
            if("Direct" in name) :
                templateName="EE-Gamma-Direct.in"
            elif("Single-Resolved" in name) :
                templateName="EE-Gamma-Single-Resolved.in"
            elif("Double-Resolved" in name) :
                templateName="EE-Gamma-Double-Resolved.in"
            else :
                templateName="EE.in"
        elif collider=="GammaGamma" :
            templateName="GammaGamma.in"
        elif collider != "BFactory" :
            templateName= "%s.in" % collider
        else :
            templateName= "EE.in"
    else :
        if have_hadronic_collider :
            templateName= "Hadron-%s.in" % simulation 
        elif collider != "BFactory" :
            templateName= "%s-%s.in" % (collider,simulation) 
        else :
            templateName= "EE-%s.in" % simulation
    # work out the name of the parameter file
    parameterName="-".join(name.split("-")[istart:])
    del istart
    return (simulation,templateName,parameterName,parameters)

def addLeptonPairCut(minmass,maxmass):
    return "set /Herwig/Cuts/LeptonPairMassCut:MinMass %s*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass %s*GeV\n" %(minmass,maxmass)

def setHardProcessWidthToZero(list1):
  res=""
  for i in list1:
    res+="set /Herwig/Particles/"+i+":HardProcessWidth 0.\n"
  return res

def jet_kt_cut(ktmin,ktmax=-1.):
    output = "set /Herwig/Cuts/JetKtCut:MinKT {E}*GeV\n".format(E=ktmin)
    if ktmax>0 :
        output += "set /Herwig/Cuts/JetKtCut:MaxKT {E}*GeV\n".format(E=ktmax)
    return output

def mhat_cut(mmin,mmax=-1):
    output = "set /Herwig/Cuts/Cuts:MHatMin {E}*GeV\n".format(E=mmin)
    if mmax>0. :
        output += "set /Herwig/Cuts/Cuts:MHatMax {E}*GeV\n".format(E=mmax)
    return output

def mhat_minm_maxm(e1,e2,e3):
    return """\
set /Herwig/Cuts/Cuts:MHatMin {e1}*GeV
set /Herwig/Cuts/MassCut:MinM {e2}*GeV
set /Herwig/Cuts/MassCut:MaxM {e3}*GeV
""".format(**locals())

def collider_lumi(energy):
    return "set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy {E}*GeV\n".format(E=energy)

def insert_ME(me,process=None,ifname='Process',subprocess="SubProcess"):
    result = "insert /Herwig/MatrixElements/{subprocess}:MatrixElements 0 /Herwig/MatrixElements/{me}\n".format(**locals())
    if process is not None:
        result += "set /Herwig/MatrixElements/{me}:{ifname} {process}".format(**locals())
        
    return result

def particlegroup(factory,name,*particles):
    directory="MatrixElements/Matchbox"
    if(factory!="Factory") : directory="Merging"
    result = ["do /Herwig/{dir}/{fact}:StartParticleGroup {n}".format(n=name,fact=factory,dir=directory)]
    for p in particles:
        result.append(
            "insert /Herwig/{dir}/{fact}:ParticleGroup 0 /Herwig/Particles/{p}".format(p=p,fact=factory,dir=directory)
        )
    result.append("do /Herwig/{dir}/{fact}:EndParticleGroup".format(fact=factory,dir=directory))
    return '\n'.join(result)

didaddfirstjet=False
def addFirstJet(ptcut,ptmax=""):
    global didaddfirstjet
    if(didaddfirstjet):
      logging.error("Can only add jetcut once.")
      sys.exit(1)
  
    res="set  /Herwig/Cuts/Cuts:JetFinder  /Herwig/Cuts/JetFinder\n"
    res+="insert  /Herwig/Cuts/Cuts:MultiCuts 0  /Herwig/Cuts/JetCuts\n"
    res+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/FirstJet\n"
    if(ptcut!=""):
        res+="set /Herwig/Cuts/FirstJet:PtMin %s*GeV\n" % ptcut
    if(ptmax!=""):
        res+="set /Herwig/Cuts/FirstJet:PtMax %s*GeV\n" % ptmax
    didaddfirstjet=True
    return res

didaddsecondjet=False
def addSecondJet(ptcut):
    global didaddsecondjet
    if(didaddsecondjet):
      logging.error("Can only add second jetcut once.")
      sys.exit(1)
    res="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
    res+="set /Herwig/Cuts/SecondJet:PtMin "+ptcut+"*GeV\n"
    didaddsecondjet=True
    return res

didaddjetpair=False
def addJetPairCut(minmass,maxmass=""):
  global didaddjetpair
  if(didaddjetpair):
      logging.error("Can only add second jetcut once.")
      sys.exit(1)
  res="""\
create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so
set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet
set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet
insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass
set /Herwig/Cuts/JetPairMass:MassMin {mm}*GeV
""".format(mm=minmass)
  if maxmass != "" :
      res+= "set /Herwig/Cuts/JetPairMass:MassMin %s*GeV\n" %maxmass
  didaddjetpair=True
  return res
