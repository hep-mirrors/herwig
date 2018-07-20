#! /usr/bin/env python
import logging,sys,os
from string import strip, Template

import sys
if sys.version_info[:3] < (2,4,0):
    print "rivet scripts require Python version >= 2.4.0... exiting"
    sys.exit(1)

if __name__ == "__main__":
    import logging
    from optparse import OptionParser, OptionGroup
    parser = OptionParser(usage="%prog name [...]")


simulation=""

numberOfAddedProcesses=0
def addProcess(thefactory,theProcess,Oas,Oew,scale,mergedlegs,NLOprocesses):
    global numberOfAddedProcesses
    global simulation
    numberOfAddedProcesses+=1
    res ="set "+thefactory+":OrderInAlphaS "+Oas+"\n"
    res+="set "+thefactory+":OrderInAlphaEW "+Oew+"\n"
    res+="do "+thefactory+":Process "+theProcess+" "
    if ( mergedlegs != 0 ):
      if simulation!="Merging":
          print "simulation is not Merging, trying to add merged legs."
          sys.exit(1)
      res+="["
      for j in range(mergedlegs):
        res+=" j "
      res+="]"
    res+="\n"
    if (NLOprocesses!=0):
       if simulation!="Merging":
          print "simulation is not Merging, trying to add NLOProcesses."
          sys.exit(1)
       res+="set MergingFactory:NLOProcesses %s \n" % NLOprocesses
    if ( scale != "" ):
      res+="set "+thefactory+":ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/"+scale+"\n"
    return res


def addLeptonPairCut(minmass,maxmass):
    return "set /Herwig/Cuts/LeptonPairMassCut:MinMass "+minmass+"*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass "+maxmass+"*GeV\n"

didaddfirstjet=False
def addFirstJet(ptcut):
    global didaddfirstjet
    if(didaddfirstjet):
      logging.error("Can only add jetcut once.")
      sys.exit(1)
  
    res="set  /Herwig/Cuts/Cuts:JetFinder  /Herwig/Cuts/JetFinder\n"
    res+="insert  /Herwig/Cuts/Cuts:MultiCuts 0  /Herwig/Cuts/JetCuts\n"
    res+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/FirstJet\n"
    if(ptcut!=""):
        res+="set /Herwig/Cuts/FirstJet:PtMin "+ptcut+".*GeV\n"
    didaddfirstjet=True
    return res



didaddsecondjet=False
def addSecondJet(ptcut):
    global didaddsecondjet
    if(didaddsecondjet):
      logging.error("Can only add second jetcut once.")
      sys.exit(1)
    res="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
    res+="set /Herwig/Cuts/SecondJet:PtMin "+ptcut+".*GeV\n"
    didaddsecondjet=True
    return res


didaddjetpair=False
def addJetPairCut(minmass):
  global didaddjetpair
  if(didaddjetpair):
      logging.error("Can only add second jetcut once.")
      sys.exit(1)
  res="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
  res+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
  res+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
  res+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
  res+="set /Herwig/Cuts/JetPairMass:MassMin "+minmass+".*GeV\n"
  didaddjetpair=True
  return res

addedBRReweighter=False
def addBRReweighter():
  global addedBRReweighter
  if(addedBRReweighter):
    logging.error("Can only add BRReweighter once.")
    sys.exit(1)
  res="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
  res+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
  addedBRReweighter=True
  return res

def setHardProcessWidthToZero(list1):
  res=""
  for i in list1:
    res+="set /Herwig/Particles/"+i+":HardProcessWidth 0.\n"
  return res

selecteddecaymode=False
def selectDecayMode(particle,decaymodes):
  global selecteddecaymode
  res="do /Herwig/Particles/"+particle+":SelectDecayModes"
  for decay in decaymodes:
    res+=" /Herwig/Particles/"+particle+"/"+decay
  res+="\n"
  selecteddecaymode=True
  return res



(opts, args) = parser.parse_args()
## Check args
if len(args) != 1:
    logging.error("Must specify at least input file")
    sys.exit(1)

name=args[0]

# settings for four flavour scheme
fourFlavour="read Matchbox/FourFlavourScheme.in\ndo /Herwig/MatrixElements/Matchbox/Factory:StartParticleGroup bjet\ninsert /Herwig/MatrixElements/Matchbox/Factory:ParticleGroup 0 /Herwig/Particles/b\ninsert /Herwig/MatrixElements/Matchbox/Factory:ParticleGroup 0 /Herwig/Particles/bbar\ninsert /Herwig/MatrixElements/Matchbox/Factory:ParticleGroup 0 /Herwig/Particles/c\ninsert /Herwig/MatrixElements/Matchbox/Factory:ParticleGroup 0 /Herwig/Particles/cbar\ninsert /Herwig/MatrixElements/Matchbox/Factory:ParticleGroup 0 /Herwig/Particles/s\ninsert /Herwig/MatrixElements/Matchbox/Factory:ParticleGroup 0 /Herwig/Particles/sbar\ninsert /Herwig/MatrixElements/Matchbox/Factory:ParticleGroup 0 /Herwig/Particles/d\ninsert /Herwig/MatrixElements/Matchbox/Factory:ParticleGroup 0 /Herwig/Particles/dbar\ninsert /Herwig/MatrixElements/Matchbox/Factory:ParticleGroup 0 /Herwig/Particles/u\ninsert /Herwig/MatrixElements/Matchbox/Factory:ParticleGroup 0 /Herwig/Particles/ubar\ninsert /Herwig/MatrixElements/Matchbox/Factory:ParticleGroup 0 /Herwig/Particles/g\ndo /Herwig/MatrixElements/Matchbox/Factory:EndParticleGroup\nset /Herwig/Cuts/MatchboxJetMatcher:Group bjet\n"

collider=""
# select the template to load
# collider
parameters = {}
if(name.find("BFactory")==0) :
    collider="BFactory"
elif(name.find("LEP")==0) :
    collider="LEP"
elif(name.find("DIS")==0) :
    collider="DIS"
elif(name.find("TVT")==0) :
    collider="TVT"
elif(name.find("LHC-GammaGamma")==0) :
    collider="LHC-GammaGamma"
elif(name.find("LHC")==0) :
    collider="LHC"
elif(name.find("ISR")==0) :
    collider="ISR"
elif(name.find("SppS")==0) :
    collider="SppS"
elif(name.find("Star")==0) :
    collider="Star"

thefactory="Factory"

istart = 1
print name

parameters["shower"]=""
parameters["bscheme"]=""

# Dipole shower with Matchbox Powheg
if(name.find("Dipole-Matchbox-Powheg")>0) :
    istart = 4
    simulation="Matchbox"
    parameters["shower"]  = "read Matchbox/Powheg-DipoleShower.in\n"

    # Dipole shower with internal Powheg - Todo: Finish modifying template files.
    '''
    elif(name.find("Dipole-Powheg")>0) :
    istart = 3
    simulation="Powheg"
    parameters["shower"]  = "set /Herwig/EventHandlers/EventHandler:CascadeHandler /Herwig/DipoleShower/DipoleShowerHandler\nread snippets/Dipole_AutoTune_prel.in\n"
    '''

# Dipole shower with MCatNLO
elif(name.find("Dipole-MCatNLO")>0) :
    istart = 3
    simulation="Matchbox"
    parameters["shower"]  = "read Matchbox/MCatNLO-DipoleShower.in\n" 

# Dipole shower with Matchbox LO
elif(name.find("Dipole-Matchbox-LO")>0) :
    istart = 4
    simulation="Matchbox"
    parameters["shower"]  = "read Matchbox/LO-DipoleShower.in\n" 

# Dipole shower with internal LO
elif(name.find("Dipole")>0) :
    istart = 2
    simulation=""
    parameters["shower"]  = "set /Herwig/EventHandlers/EventHandler:CascadeHandler /Herwig/DipoleShower/DipoleShowerHandler\nread snippets/Dipole_AutoTune_prel.in\n"
    
# AO shower with Matchbox Powheg
elif(name.find("Matchbox-Powheg")>0) :
    istart = 3
    simulation="Matchbox"
    parameters["shower"] = "read Matchbox/Powheg-DefaultShower.in\n"

# AO shower with MCatNLO
elif(name.find("Matchbox")>0) :
    istart = 2
    simulation="Matchbox"
    parameters["shower"] = "read Matchbox/MCatNLO-DefaultShower.in\n"

# AO shower with inernal Powheg    
elif(name.find("Powheg")>0) :
    istart = 2
    simulation="Powheg"

# Dipole shower with merging    
elif(name.find("Merging")>0) :
    istart = 2
    simulation="Merging"
    thefactory="MergingFactory"
    
# Flavour settings for Matchbox    
if(simulation=="Matchbox") :
    parameters["bscheme"] = "read Matchbox/FiveFlavourScheme.in\n"
    
    if(parameters["shower"].find("Dipole")>=0) :
        parameters["bscheme"] += "read Matchbox/FiveFlavourNoBMassScheme.in\n"
        
    if(collider.find("DIS")<0 and collider.find("LEP")<0 ) :
        parameters["nlo"] = "read Matchbox/MadGraph-OpenLoops.in\n"

# Flavour settings for dipole shower with internal ME
if ( simulation == "" and parameters["shower"].find("Dipole")>=0 ) :
    parameters["bscheme"] = "read snippets/DipoleShowerFiveFlavours.in"
    

if(collider=="") :
    logging.error("Can\'t find collider")
    sys.exit(1)

# find the template
if(simulation=="") :
    if(collider.find("LHC-GammaGamma") >=0) :
        istart += 1
        templateName="Hadron-Gamma.in"
    elif(collider.find("TVT")>=0 or collider.find("LHC") >=0 or
       collider.find("ISR")>=0 or collider.find("SppS")>=0 or
       collider.find("Star")>=0) :
        templateName="Hadron.in"
    elif(collider.find("BFactory")<0) :
        templateName= "%s.in" % (collider)
    else :
        templateName= "LEP.in"
else :
    if(collider.find("TVT")>=0 or collider.find("LHC") >=0 or
       collider.find("ISR")>=0 or collider.find("SppS")>=0 or
       collider.find("Star")>=0) :
        templateName= "Hadron-%s.in" % (simulation) 
    elif(collider.find("BFactory")<0) :
        templateName= "%s-%s.in" % (collider,simulation) 
    else :
        templateName= "LEP-%s.in" % (simulation) 
with open(os.path.join("Rivet/Templates",templateName), 'r') as f:
    templateText = f.read()
template = Template( templateText )
# work out the name of the parameter file
nameSplit=name.split("-")
parameterName=nameSplit[istart]
for i in range(istart+1,len(nameSplit)) :
    parameterName += "-%s" % nameSplit[i] 

# work out the process and parameters
process=""
# Bfactory
if(collider=="BFactory") :
    if(simulation=="") :
        if(parameterName=="10.58-res") :
            process += "\ncreate Herwig::MEee2VectorMeson /Herwig/MatrixElements/MEUpsilon HwMELepton.so\nset /Herwig/MatrixElements/MEUpsilon:VectorMeson /Herwig/Particles/Upsilon(4S)\nset /Herwig/MatrixElements/MEUpsilon:Coupling 0.0004151809\ninsert /Herwig/MatrixElements/SubProcess:MatrixElements 0 /Herwig/MatrixElements/MEUpsilon"
        elif(parameterName=="10.58") :
            process += "\ncreate Herwig::MEee2VectorMeson /Herwig/MatrixElements/MEUpsilon HwMELepton.so\nset /Herwig/MatrixElements/MEUpsilon:VectorMeson /Herwig/Particles/Upsilon(4S)\nset /Herwig/MatrixElements/MEUpsilon:Coupling 0.0004151809\ninsert /Herwig/MatrixElements/SubProcess:MatrixElements 0 /Herwig/MatrixElements/MEUpsilon\n"
            process += "set /Herwig/MatrixElements/MEee2gZ2qq:MaximumFlavour 4\n"
        else :
            process+="insert  /Herwig/MatrixElements/SubProcess:MatrixElements 0 /Herwig/MatrixElements/MEee2gZ2qq\n"
            process+= "set /Herwig/MatrixElements/MEee2gZ2qq:MaximumFlavour 4\n"
    elif(simulation=="Powheg") :
        process = "set /Herwig/MatrixElements/PowhegMEee2gZ2qq:MaximumFlavour 4\n"
    elif(simulation=="Matchbox" ) :
        process =addProcess(thefactory,"e- e+ -> u ubar","0","2","",0,0)
        process+=addProcess(thefactory,"e- e+ -> d dbar","0","2","",0,0)
        process+=addProcess(thefactory,"e- e+ -> c cbar","0","2","",0,0)
        process+=addProcess(thefactory,"e- e+ -> s sbar","0","2","",0,0)
    elif(simulation=="Merging" ) :
        logging.warning("BFactory not explicitly tested for %s " % simulation)
        sys.exit(0)

# DIS
elif(collider=="DIS") :
    if(simulation=="") :
        if(parameterName.find("NoME")>=0) :
            process = "set /Herwig/Shower/ShowerHandler:HardEmission None"
            parameterName=parameterName.replace("NoME-","")
        else :
            process = ""
    elif(simulation=="Powheg") :
        process = ""
    elif(simulation=="Matchbox" ) :
      if(parameterName.find("e-")>=0) :
            process=addProcess(thefactory,"e- p -> e- j","0","2","",0,0)
      else :
            process=addProcess(thefactory,"e+ p -> e+ j","0","2","",0,0)
    elif(simulation=="Merging" ) :
        if(parameterName.find("e-")>=0) :
            process=addProcess(thefactory,"e- p -> e- j","0","2","",2,2)
        else :
            process=addProcess(thefactory,"e+ p -> e+ j","0","2","",2,2)

# LEP
elif(collider=="LEP") :
    if(simulation=="") :
        if(parameterName.find("gg")>=0) :
            process ="create Herwig::MEee2Higgs2SM /Herwig/MatrixElements/MEee2Higgs2SM\n"
            process+="insert /Herwig/MatrixElements/SubProcess:MatrixElements[0] /Herwig/MatrixElements/MEee2Higgs2SM\n"
            process+="set /Herwig/MatrixElements/MEee2Higgs2SM:Allowed Gluon\n"

        else :
            process="insert  /Herwig/MatrixElements/SubProcess:MatrixElements 0 /Herwig/MatrixElements/MEee2gZ2qq\n"
            if(parameterName=="10") :
                process+="set /Herwig/MatrixElements/MEee2gZ2qq:MaximumFlavour 4"
    elif(simulation=="Powheg") :
        process=""
        if(parameterName=="10") :
            process="set /Herwig/MatrixElements/PowhegMEee2gZ2qq:MaximumFlavour 4"
    elif(simulation=="Matchbox" ) :
        if(parameterName=="10") :
            process =addProcess(thefactory,"e- e+ -> u ubar","0","2","",0,0)
            process+=addProcess(thefactory,"e- e+ -> d dbar","0","2","",0,0)
            process+=addProcess(thefactory,"e- e+ -> c cbar","0","2","",0,0)
            process+=addProcess(thefactory,"e- e+ -> s sbar","0","2","",0,0)
        else :
            process=addProcess(thefactory,"e- e+ -> j j","0","2","",0,0)

    elif(simulation=="Merging" ) :
        if(parameterName=="10") :
          process=addProcess(thefactory,"e- e+ -> j j","0","2","",2,2)
          process+="read Matchbox/FourFlavourScheme.in"
        else :
          process=addProcess(thefactory,"e- e+ -> j j","0","2","",2,2)
# TVT
elif(collider=="TVT") :
    process="set /Herwig/Generators/EventGenerator:EventHandler:BeamB /Herwig/Particles/pbar-\n"
    if(parameterName.find("Run-II")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 1960.0\n"
    elif(parameterName.find("Run-I")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 1800.0\n"
    elif(parameterName.find("900")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 900.0\n"
    elif(parameterName.find("630")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 630.0\n"
    elif(parameterName.find("300")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 300.0\n"

    if(simulation=="") :
        if(parameterName.find("PromptPhoton")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEGammaJet\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 15.\n"
        elif(parameterName.find("DiPhoton-GammaGamma")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEGammaGamma\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaGamma","")
        elif(parameterName.find("DiPhoton-GammaJet")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEGammaJet\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaJet","")
        elif(parameterName.find("UE")>=0) :
            if (parameters["shower"].find("Dipole")>=0):
                process+="read snippets/MB-DipoleShower.in\n"
            else:
                process+="read snippets/MB.in\n"
            process+="read snippets/Diffraction.in\n"
                
            process += "set /Herwig/Decays/DecayHandler:LifeTimeOption 0\n"
            process += "set /Herwig/Decays/DecayHandler:MaxLifeTime 10*mm\n"
        elif(parameterName.find("Jets")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEQCD2to2\n"
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("Run-II-Jets-10")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/Cuts:MHatMin 500.*GeV\n"
            elif(parameterName.find("Run-II-Jets-11")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/Cuts:MHatMin 900.*GeV\n"
            elif(parameterName.find("Run-I-Jets-1")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
            elif(parameterName.find("Run-I-Jets-2")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 40.\n"
            elif(parameterName.find("Run-I-Jets-3")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 65.\n"
            elif(parameterName.find("Run-I-Jets-4")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 90.\n"
            elif(parameterName.find("Run-I-Jets-5")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 160.\n"
            elif(parameterName.find("Run-I-Jets-6")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/Cuts:MHatMin 100.*GeV\n"
            elif(parameterName.find("Run-I-Jets-7")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/Cuts:MHatMin 400.*GeV\n"
            elif(parameterName.find("Run-I-Jets-8")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/Cuts:MHatMin 700.*GeV\n"
            elif(parameterName.find("Run-II-Jets-0")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 15.\n"
            elif(parameterName.find("Run-II-Jets-1")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 25.\n"
            elif(parameterName.find("Run-II-Jets-2")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 40.\n"
            elif(parameterName.find("Run-II-Jets-3")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 60.\n"
            elif(parameterName.find("Run-II-Jets-4")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 85.\n"
            elif(parameterName.find("Run-II-Jets-5")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 110.\n"
            elif(parameterName.find("Run-II-Jets-6")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 160.\n"
            elif(parameterName.find("Run-II-Jets-7")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 250.\n"
            elif(parameterName.find("Run-II-Jets-8")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/Cuts:MHatMin 100.*GeV\n"
            elif(parameterName.find("Run-II-Jets-9")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/Cuts:MHatMin 300.*GeV\n"
            elif(parameterName.find("900-Jets-1")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 10.\n"
            elif(parameterName.find("300-Jets-1")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 6.\n"
            elif(parameterName.find("630-Jets-1")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
            elif(parameterName.find("630-Jets-2")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 40.\n"
            elif(parameterName.find("630-Jets-3")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 75.\n"
            elif(parameterName.find("900-Jets-1")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 10.\n"
        elif(parameterName.find("Run-I-WZ")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Electron\ninsert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Run-I-W")>=0 or parameterName.find("Run-II-W")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Electron\n"
        elif(parameterName.find("Run-I-Z")>=0 or parameterName.find("Run-II-Z-e")>=0) :
            process +="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Run-II-Z-LowMass-mu")>=0) :
            process +="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
            process+=addLeptonPairCut("25","70")
        elif(parameterName.find("Run-II-Z-HighMass-mu")>=0) :
            process +="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
            process+=addLeptonPairCut("150","600")
        elif(parameterName.find("Run-II-Z-mu")>=0) :
            process +="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
    elif(simulation=="Powheg") :
        if(parameterName.find("Run-I-WZ")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2W2ff\nset PowhegMEqq2W2ff:Process Electron\ninsert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Run-I-W")>=0 or parameterName.find("Run-II-W")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2W2ff\nset PowhegMEqq2W2ff:Process Electron\n"
        elif(parameterName.find("Run-I-Z")>=0 or parameterName.find("Run-II-Z-e")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Run-II-Z-LowMass-mu")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
            process+=addLeptonPairCut("25","70")
        elif(parameterName.find("Run-II-Z-HighMass-mu")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
            process+=addLeptonPairCut("150","600")
        elif(parameterName.find("Run-II-Z-mu")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("DiPhoton-GammaGamma")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEGammaGammaPowheg\n"
            process+="set MEGammaGammaPowheg:Process GammaGamma\n"
            process+="insert SubProcess:MatrixElements[0] MEGammaGamma\n"
            process+="set MEGammaGamma:Process gg\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaGamma","")
        elif(parameterName.find("DiPhoton-GammaJet")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEGammaGammaPowheg\n"
            process+="set MEGammaGammaPowheg:Process VJet\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaJet","")
    elif(simulation=="Matchbox" or simulation=="Merging" ) :
        if(parameterName.find("Jets")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p -> j j","2","0","MaxJetPtScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p p -> j j","2","0","MaxJetPtScale",1,0)
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("Run-II-Jets-10")>=0) :
                process+=addFirstJet("30")
                process+=addSecondJet("25")
                process+=addJetPairCut("500")
            elif(parameterName.find("Run-II-Jets-11")>=0) :
                process+=addFirstJet("30")
                process+=addSecondJet("25")
                process+=addJetPairCut("900")
            elif(parameterName.find("Run-II-Jets-12")>=0) :
                process+=addFirstJet("30")
                process+=addSecondJet("25")
                process+=addJetPairCut("300")
            elif(parameterName.find("Run-I-Jets-1")>=0) :
                process+=addFirstJet("20")
            elif(parameterName.find("Run-I-Jets-2")>=0) :
                process+=addFirstJet("40")
            elif(parameterName.find("Run-I-Jets-3")>=0) :
                process+=addFirstJet("65")
            elif(parameterName.find("Run-I-Jets-4")>=0) :
                process+=addFirstJet("90")
            elif(parameterName.find("Run-I-Jets-5")>=0) :
                process+=addFirstJet("160")
            elif(parameterName.find("Run-I-Jets-6")>=0) :
                process+=addFirstJet("30")
                process+=addSecondJet("25")
                process+=addJetPairCut("100")
            elif(parameterName.find("Run-I-Jets-7")>=0) :
                process+=addFirstJet("30")
                process+=addSecondJet("25")
                process+=addJetPairCut("400")
            elif(parameterName.find("Run-I-Jets-8")>=0) :
                process+=addFirstJet("30")
                process+=addSecondJet("25")
                process+=addJetPairCut("700")
            elif(parameterName.find("Run-II-Jets-0")>=0) :
                process+=addFirstJet("15")
            elif(parameterName.find("Run-II-Jets-1")>=0) :
                process+=addFirstJet("25")
            elif(parameterName.find("Run-II-Jets-2")>=0) :
                process+=addFirstJet("40")
            elif(parameterName.find("Run-II-Jets-3")>=0) :
                process+=addFirstJet("60")
            elif(parameterName.find("Run-II-Jets-4")>=0) :
                process+=addFirstJet("85")
            elif(parameterName.find("Run-II-Jets-5")>=0) :
                process+=addFirstJet("110")
            elif(parameterName.find("Run-II-Jets-6")>=0) :
                process+=addFirstJet("160")
            elif(parameterName.find("Run-II-Jets-7")>=0) :
                process+=addFirstJet("250")
            elif(parameterName.find("Run-II-Jets-8")>=0) :
                process+=addFirstJet("30")
                process+=addSecondJet("25")
                process+=addJetPairCut("100")
            elif(parameterName.find("Run-II-Jets-9")>=0) :
                process+=addFirstJet("30")
                process+=addSecondJet("25")
                process+=addJetPairCut("300")
            elif(parameterName.find("900-Jets-1")>=0) :
                process+=addFirstJet("10")
            elif(parameterName.find("300-Jets-1")>=0) :
                process+=addFirstJet("6")
            elif(parameterName.find("630-Jets-1")>=0) :
                process+=addFirstJet("20")
            elif(parameterName.find("630-Jets-2")>=0) :
                process+=addFirstJet("40")
            elif(parameterName.find("630-Jets-3")>=0) :
                process+=addFirstJet("75")
            elif(parameterName.find("900-Jets-1")>=0) :
                process+=addFirstJet("10")
            else :
                logging.error("Exit 00007")
                sys.exit(1)
        elif(parameterName.find("Run-I-WZ")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p pbar e+ e-","0","2","LeptonPairMassScale",0,0)
                process+=addProcess(thefactory,"p pbar e+ nu","0","2","LeptonPairMassScale",0,0)
                process+=addProcess(thefactory,"p pbar e- nu","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
                process+="do "+thefactory+":StartParticleGroup epm\n"
                process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/e+\n"
                process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/e-\n"
                process+="do "+thefactory+":EndParticleGroup\n"
                process+="do "+thefactory+":StartParticleGroup epmnu\n"
                process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/e+\n"
                process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/e-\n"
                process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/nu_e\n"
                process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/nu_ebar\n"
                process+="do "+thefactory+":EndParticleGroup\n"
                process+=addProcess(thefactory,"p pbar epm epmnu","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("Run-I-W")>=0 or parameterName.find("Run-II-W")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p pbar e+ nu","0","2","LeptonPairMassScale",0,0)
                process+=addProcess(thefactory,"p pbar e- nu","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
                process+="do "+thefactory+":StartParticleGroup epm\n"
                process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/e+\n"
                process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/e-\n"
                process+="do "+thefactory+":EndParticleGroup\n"
                process+=addProcess(thefactory,"p pbar epm nu","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("Run-I-Z")>=0 or parameterName.find("Run-II-Z-e")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p pbar e+ e-","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p pbar e+ e-","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("Run-II-Z-LowMass-mu")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p pbar mu+ mu-","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p pbar mu+ mu-","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("25","70")
        elif(parameterName.find("Run-II-Z-HighMass-mu")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p pbar mu+ mu-","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p pbar mu+ mu-","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("150","600")
        elif(parameterName.find("Run-II-Z-mu")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p pbar mu+ mu-","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p pbar mu+ mu-","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("60","120")
# Star
elif(collider=="Star" ) :
    process = "set /Herwig/Decays/DecayHandler:LifeTimeOption 0\n"
    process+= "set /Herwig/Decays/DecayHandler:MaxLifeTime 10*mm\n"
    process+= "set /Herwig/Generators/EventGenerator:EventHandler:BeamB /Herwig/Particles/p+\n"
    process+= "set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 200.0\n"
    process+= "set /Herwig/Cuts/Cuts:X2Min 0.01\n"
    if(simulation=="") :
        if(parameterName.find("UE")>=0) :
            if (parameters["shower"].find("Dipole")>=0):
                process+="read snippets/MB-DipoleShower.in\n"
            else:
                process+="read snippets/MB.in\n"    
            process+="read snippets/Diffraction.in\n"
            
            
        else :
            process+="insert SubProcess:MatrixElements[0] MEQCD2to2\n"
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("Jets-1")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 2.\n"
            elif(parameterName.find("Jets-2")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 5.\n"
            elif(parameterName.find("Jets-3")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
            elif(parameterName.find("Jets-4")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 25.\n"
    else :
        logging.error("Star not supported for %s " % simulation)
        sys.exit(1)
# ISR and SppS
elif(collider=="ISR" or collider =="SppS" ) :
    process="set /Herwig/Decays/DecayHandler:LifeTimeOption 0\n"
    process+="set /Herwig/Decays/DecayHandler:MaxLifeTime 10*mm\n"
    if(collider=="SppS") :
        process ="set /Herwig/Generators/EventGenerator:EventHandler:BeamB /Herwig/Particles/pbar-\n"
    if(parameterName.find("30")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 30.4\n"
    elif(parameterName.find("44")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 44.4\n"
    elif(parameterName.find("53")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 53.0\n"
    elif(parameterName.find("62")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 62.2\n"
    elif(parameterName.find("63")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 63.0\n"
    elif(parameterName.find("200")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 200.0\n"
    elif(parameterName.find("500")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 500.0\n"
    elif(parameterName.find("546")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 546.0\n"
    elif(parameterName.find("900")>=0) :
        process+="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 900.0\n"
    if(simulation=="") :
        if (parameters["shower"].find("Dipole")>=0):
            process+="read snippets/MB-DipoleShower.in\n"
        else:
            process+="read snippets/MB.in\n"
        process+="read snippets/Diffraction.in\n"
    else :
        logging.error(" SppS and ISR not supported for %s " % simulation)
        sys.exit(1)
# LHC
elif(collider=="LHC") :
    if(parameterName.find("7-")==0) :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 7000.0\n"
    elif(parameterName.find("8-")==0) :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 8000.0\n"
    elif(parameterName.find("13-")==0) :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 13000.0\n"
    elif(parameterName.find("900")==0) :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 900.0\n"
    elif(parameterName.find("2360")==0) :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 2360.0\n"
    elif(parameterName.find("2760")==0) :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 2760.0\n"
    else :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 7000.0\n"
    if(simulation=="") :
        if(parameterName.find("8-VBF")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEPP2HiggsVBF\n"
        elif(parameterName.find("VBF")>=0) :
            process+=selectDecayMode("h0",["h0->tau-,tau+;"])
            addedBRReweighter = True
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+="insert SubProcess:MatrixElements[0] MEPP2HiggsVBF\n"
        elif(parameterName.find("ggHJet")>=0) :
            process+=selectDecayMode("h0",["h0->tau-,tau+;"])
            addedBRReweighter = True
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+="insert SubProcess:MatrixElements[0] MEHiggsJet\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
        elif(parameterName.find("8-ggH")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEHiggs\n"
            process+="insert SubProcess:MatrixElements[0] MEHiggsJet\n"
            process+="set MEHiggsJet:Process qqbar\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("ggH")>=0) :
            process+=selectDecayMode("h0",["h0->tau-,tau+;"])
            addedBRReweighter = True
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+="insert SubProcess:MatrixElements[0] MEHiggs\n"
            process+="insert SubProcess:MatrixElements[0] MEHiggsJet\n"
            process+="set MEHiggsJet:Process qqbar\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("PromptPhoton")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEGammaJet\n"
            if(parameterName.find("PromptPhoton-1")>=0) :
                process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            elif(parameterName.find("PromptPhoton-2")>=0) :
                process+="set /Herwig/Cuts/PhotonKtCut:MinKT 25.\n"
            elif(parameterName.find("PromptPhoton-3")>=0) :
                process+="set /Herwig/Cuts/PhotonKtCut:MinKT 80.\n"
            elif(parameterName.find("PromptPhoton-4")>=0) :
                process+="set /Herwig/Cuts/PhotonKtCut:MinKT 150.\n"
        elif(parameterName.find("DiPhoton-GammaGamma")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEGammaGamma\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaGamma","")
        elif(parameterName.find("DiPhoton-GammaJet")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEGammaJet\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaJet","")
        elif(parameterName.find("8-WH")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEPP2WH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("8-ZH")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEPP2ZH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("WH")>=0) :
            process+=selectDecayMode("h0",["h0->b,bbar;"])
            process+=selectDecayMode("W+",["W+->nu_e,e+;",
                                           "W+->nu_mu,mu+;"])
            addedBRReweighter = True
            process+="insert SubProcess:MatrixElements[0] MEPP2WH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("ZH")>=0) :
            process+=selectDecayMode("h0",["h0->b,bbar;"])
            process+=selectDecayMode("Z0",["Z0->e-,e+;",
                                           "Z0->mu-,mu+;"])
            addedBRReweighter = True
            process+="insert SubProcess:MatrixElements[0] MEPP2ZH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("UE")>=0) :
            if (parameters["shower"].find("Dipole")>=0):
                process+="read snippets/MB-DipoleShower.in\n"
            else:
                process+="set /Herwig/Shower/ShowerHandler:IntrinsicPtGaussian 2.2*GeV\n"                
                process+="read snippets/MB.in\n"
            process+="read snippets/Diffraction.in\n"
            if(parameterName.find("Long")>=0) :
                process += "set /Herwig/Decays/DecayHandler:MaxLifeTime 100*mm\n"
        elif(parameterName.find("7-DiJets")>=0 or parameterName.find("8-DiJets")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEQCD2to2\n"
            process+="set MEQCD2to2:MaximumFlavour 5\n"
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("-A")>=0) :
               process+="set /Herwig/Cuts/JetKtCut:MinKT 45.\n"
               process+="set /Herwig/Cuts/JetKtCut:MinEta -3.\n"
               process+="set /Herwig/Cuts/JetKtCut:MaxEta  3.\n"
            elif(parameterName.find("-B")>=0) :
               process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
               process+="set /Herwig/Cuts/JetKtCut:MinEta -2.7\n"
               process+="set /Herwig/Cuts/JetKtCut:MaxEta  2.7\n"
            elif(parameterName.find("-C")>=0) :
               process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
               process+="set /Herwig/Cuts/JetKtCut:MinEta -4.8\n"
               process+="set /Herwig/Cuts/JetKtCut:MaxEta  4.8\n"
            if(parameterName.find("DiJets-1")>=0) :
                process+="set /Herwig/Cuts/Cuts:MHatMin 90.*GeV\n"
            elif(parameterName.find("DiJets-2")>=0) :
                process+="set /Herwig/Cuts/Cuts:MHatMin 200.*GeV\n"
            elif(parameterName.find("DiJets-3")>=0) :
                process+="set /Herwig/Cuts/Cuts:MHatMin 450.*GeV\n"
            elif(parameterName.find("DiJets-4")>=0) :
                process+="set /Herwig/Cuts/Cuts:MHatMin 750.*GeV\n"
            elif(parameterName.find("DiJets-5")>=0) :
                process+="set /Herwig/Cuts/Cuts:MHatMin 950.*GeV\n"
            elif(parameterName.find("DiJets-6")>=0) :
                process+="set /Herwig/Cuts/Cuts:MHatMin 1550.*GeV\n"
            elif(parameterName.find("DiJets-7")>=0) :
                process+="set /Herwig/Cuts/Cuts:MHatMin 2150.*GeV\n"
            elif(parameterName.find("DiJets-8")>=0) :
                process+="set /Herwig/Cuts/Cuts:MHatMin 2750.*GeV\n"
        elif(parameterName.find("7-Jets")>=0 or parameterName.find("8-Jets")>=0 or \
             parameterName.find("13-Jets")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEQCD2to2\n"
            process+="set MEQCD2to2:MaximumFlavour 5\n"
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("Jets-10")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 1800.\n"
            elif(parameterName.find("Jets-0")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 5.\n"
            elif(parameterName.find("Jets-1")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 10.\n"
            elif(parameterName.find("Jets-2")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
            elif(parameterName.find("Jets-3")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 40.\n"
            elif(parameterName.find("Jets-4")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 70.\n"
            elif(parameterName.find("Jets-5")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 150.\n"
            elif(parameterName.find("Jets-6")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 200.\n"
            elif(parameterName.find("Jets-7")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 300.\n"
            elif(parameterName.find("Jets-8")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 500.\n"
            elif(parameterName.find("Jets-9")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 800.\n"
        elif(parameterName.find("7-Charm")>=0 or \
             parameterName.find("7-Bottom")>=0) :
            if(parameterName.find("7-Bottom")>=0) :
                process+="cp MEHeavyQuark MEBottom\n" 
                process+="set MEBottom:QuarkType Bottom\n"
                process+="insert SubProcess:MatrixElements[0] MEBottom\n"
            else : 
                process+="cp MEHeavyQuark MECharm\n" 
                process+="set MECharm:QuarkType Charm\n"
                process+="insert SubProcess:MatrixElements[0] MECharm\n"
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("-0")>=0) :
                if(parameterName.find("7-Bottom")>=0) :
                    process+="set MEBottom:Process Pair\n" 
                process+="set /Herwig/Cuts/JetKtCut:MinKT 0.\n"
            elif(parameterName.find("-1")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 5.\n"
            elif(parameterName.find("-2")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
            elif(parameterName.find("-3")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 50.\n"
            elif(parameterName.find("-4")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 80.\n"
            elif(parameterName.find("-5")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 110.\n"
            elif(parameterName.find("-6")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/Cuts:MHatMin 90.*GeV\n"
            elif(parameterName.find("-7")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/Cuts:MHatMin 340.*GeV\n"
            elif(parameterName.find("-8")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/Cuts:MHatMin 500.*GeV\n"
        elif(parameterName.find("Top-L")>=0) :
            process+="set MEHeavyQuark:QuarkType Top\n"
            process+="insert SubProcess:MatrixElements[0] MEHeavyQuark\n"
            process+=selectDecayMode("t",["t->nu_e,e+,b;",
                                          "t->nu_mu,mu+,b;"])
            process+=addBRReweighter()
            
        elif(parameterName.find("Top-SL")>=0) :
            process+="set MEHeavyQuark:QuarkType Top\n"
            process+="insert SubProcess:MatrixElements[0] MEHeavyQuark\n"
            process+="set /Herwig/Particles/t:Synchronized Not_synchronized\n"
            process+="set /Herwig/Particles/tbar:Synchronized Not_synchronized\n"
            process+=selectDecayMode("t",["t->nu_e,e+,b;","t->nu_mu,mu+,b;"])
            process+=selectDecayMode("tbar",["tbar->b,bbar,cbar;",
                                             "tbar->bbar,cbar,d;",
                                             "tbar->bbar,cbar,s;",
                                             "tbar->bbar,s,ubar;",
                                             "tbar->bbar,ubar,d;"])
            process+=addBRReweighter()
            
        elif(parameterName.find("Top-All")>=0) :
            process+="set MEHeavyQuark:QuarkType Top\n"
            process+="insert SubProcess:MatrixElements[0] MEHeavyQuark\n"
        elif(parameterName.find("WZ")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process WZ\n"
            process+=selectDecayMode("W+",["W+->nu_e,e+;",
                                           "W+->nu_mu,mu+;"])
            process+=selectDecayMode("W-",["W-->nu_ebar,e-;",
                                           "W-->nu_mubar,mu-;"])
            process+=selectDecayMode("Z0",["Z0->e-,e+;",
                                           "Z0->mu-,mu+;"])
            addedBRReweighter = True
            
        elif(parameterName.find("WW-emu")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process WW\n"
            process+="set /Herwig/Particles/W+:Synchronized 0\n"
            process+="set /Herwig/Particles/W-:Synchronized 0\n"
            process+=selectDecayMode("W+",["W+->nu_e,e+;"])
            process+=selectDecayMode("W-",["W-->nu_mubar,mu-;"])
            addedBRReweighter = True
            
        elif(parameterName.find("WW-ll")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process WW\n"
            process+=selectDecayMode("W+",["W+->nu_e,e+;","W+->nu_mu,mu+;","W+->nu_tau,tau+;"])
            addedBRReweighter = True
            
        elif(parameterName.find("ZZ-ll")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process ZZ\n"
            process+=selectDecayMode("Z0",["Z0->e-,e+;",
                                           "Z0->mu-,mu+;",
                                           "Z0->tau-,tau+;"])
            addedBRReweighter = True

        elif(parameterName.find("ZZ-lv")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process ZZ\n"
            process+=selectDecayMode("Z0",["Z0->e-,e+;",
                                           "Z0->mu-,mu+;",
                                           "Z0->tau-,tau+;",
                                           "Z0->nu_e,nu_ebar;",
                                           "Z0->nu_mu,nu_mubar;",
                                           "Z0->nu_tau,nu_taubar;"])
            addedBRReweighter = True
            
        elif(parameterName.find("W-Z-e")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
            process+="insert SubProcess:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Electron\n"

        elif(parameterName.find("W-Z-mu")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
            process+="insert SubProcess:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Muon\n"
        elif(parameterName.find("W-e")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Electron\n"
        elif(parameterName.find("W-mu")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Muon\n"
        elif(parameterName.find("Z-e")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Z-mu")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("Z-LowMass-e")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
            process+="set /Herwig/Cuts/Cuts:MHatMin 20.*GeV\nset /Herwig/Cuts/MassCut:MinM 20.*GeV\nset /Herwig/Cuts/MassCut:MaxM 70.*GeV\n"
        elif(parameterName.find("Z-MedMass-e")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
            process+="set /Herwig/Cuts/Cuts:MHatMin 40.*GeV\nset /Herwig/Cuts/MassCut:MinM 40.*GeV\nset /Herwig/Cuts/MassCut:MaxM 130.*GeV\n"
        elif(parameterName.find("Z-LowMass-mu")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
            process+="set /Herwig/Cuts/Cuts:MHatMin 10.*GeV\nset /Herwig/Cuts/MassCut:MinM 10.*GeV\nset /Herwig/Cuts/MassCut:MaxM 70.*GeV\n"
        elif(parameterName.find("Z-Mass1")>=0) :
            process+="set /Herwig/Cuts/Cuts:MHatMin 10.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MinM 10.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MaxM 35.*GeV\n"
            if(parameterName.find("-e")>=0) :
                process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
            else :
                process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("Z-Mass2")>=0) :
            process+="set /Herwig/Cuts/Cuts:MHatMin 25.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MinM 25.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MaxM 70.*GeV\n"
            if(parameterName.find("-e")>=0) :
                process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
            else :
                process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("Z-Mass3")>=0) :
            process+="set /Herwig/Cuts/Cuts:MHatMin 60.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MinM 60.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MaxM 120.*GeV\n"
            if(parameterName.find("-e")>=0) :
                process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
            else :
                process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("Z-Mass4")>=0) :
            process+="set /Herwig/Cuts/Cuts:MHatMin 110.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MinM 110.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MaxM 8000.*GeV\n"
            if(parameterName.find("-e")>=0) :
                process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
            else :
                process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("W-Jet")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEWJet\nset MEWJet:WDecay Electron\n"
            if(parameterName.find("W-Jet-1-e")>=0) :
                process+="set /Herwig/Cuts/WBosonKtCut:MinKT 100.0*GeV\n"
                parameterName=parameterName.replace("W-Jet-1-e","W-Jet-e")
            elif(parameterName.find("W-Jet-2-e")>=0) :
                process+="set /Herwig/Cuts/WBosonKtCut:MinKT 190.0*GeV\n"
                parameterName=parameterName.replace("W-Jet-2-e","W-Jet-e")
            elif(parameterName.find("W-Jet-3-e")>=0) :
                process+="set /Herwig/Cuts/WBosonKtCut:MinKT 270.0*GeV\n"
                parameterName=parameterName.replace("W-Jet-3-e","W-Jet-e")
        elif(parameterName.find("Z-Jet")>=0) :
            if(parameterName.find("-e")>=0) :
                process+="insert SubProcess:MatrixElements[0] MEZJet\nset MEZJet:ZDecay Electron\n"
                if(parameterName.find("Z-Jet-0-e")>=0) :
                    process+="set /Herwig/Cuts/ZBosonKtCut:MinKT 35.0*GeV\n"
                    parameterName=parameterName.replace("Z-Jet-0-e","Z-Jet-e")
                elif(parameterName.find("Z-Jet-1-e")>=0) :
                    process+="set /Herwig/Cuts/ZBosonKtCut:MinKT 100.0*GeV\n"
                    parameterName=parameterName.replace("Z-Jet-1-e","Z-Jet-e")
                elif(parameterName.find("Z-Jet-2-e")>=0) :
                    process+="set /Herwig/Cuts/ZBosonKtCut:MinKT 190.0*GeV\n"
                    parameterName=parameterName.replace("Z-Jet-2-e","Z-Jet-e")
                elif(parameterName.find("Z-Jet-3-e")>=0) :
                    process+="set /Herwig/Cuts/ZBosonKtCut:MinKT 270.0*GeV\n"
                    parameterName=parameterName.replace("Z-Jet-3-e","Z-Jet-e")
            else :
                process+="insert SubProcess:MatrixElements[0] MEZJet\nset MEZJet:ZDecay Muon\n"
                process+="set /Herwig/Cuts/ZBosonKtCut:MinKT 35.0*GeV\n"
                parameterName=parameterName.replace("Z-Jet-0-mu","Z-Jet-mu")
        elif(parameterName.find("WGamma")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEPP2VGamma\nset MEPP2VGamma:Process 1\nset MEPP2VGamma:MassOption 1\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 10.\n"
            
            
            if(parameterName.find("-e")>=0) :
                process+=selectDecayMode("W+",["W+->nu_e,e+;"])
                process+=addBRReweighter()
            else :
                process+=selectDecayMode("W+",["W+->nu_mu,mu+;"])
                process+=addBRReweighter()
        elif(parameterName.find("ZGamma")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEPP2VGamma\nset MEPP2VGamma:Process 2\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 10.\n"
            if(parameterName.find("-e")>=0) :
                process+=selectDecayMode("Z0",["Z0->e-,e+;"])
                process+=addBRReweighter()
            else :
                process+=selectDecayMode("Z0",["Z0->mu-,mu+;"])
                process+=addBRReweighter()
        else :
            logging.error(" Process %s not supported for internal matrix elements" % name)
            sys.exit(1)
    elif(simulation=="Powheg") :
        if(parameterName.find("8-VBF")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEPP2HiggsVBF\n"
        elif(parameterName.find("VBF")>=0) :
            process+=selectDecayMode("h0",["h0->tau-,tau+;"])
            addedBRReweighter = True
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+="insert SubProcess:MatrixElements[0] PowhegMEPP2HiggsVBF\n"
        elif(parameterName.find("ggHJet")>=0) :
            logging.error(" Process %s not supported for POWHEG matrix elements" % name)
            sys.exit(1)
        elif(parameterName.find("8-ggH")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEHiggs\n"
        elif(parameterName.find("ggH")>=0) :
            process+=selectDecayMode("h0",["h0->tau-,tau+;"])
            addedBRReweighter = True
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+="insert SubProcess:MatrixElements[0] PowhegMEHiggs\n"
        elif(parameterName.find("8-WH")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEPP2WH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("8-ZH")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEPP2ZH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("WH")>=0) :
            process+=selectDecayMode("h0",["h0->b,bbar;"])
            process+=selectDecayMode("W+",["W+->nu_e,e+;",
                                           "W+->nu_mu,mu+;"])
            addedBRReweighter = True
            process+="insert SubProcess:MatrixElements[0] PowhegMEPP2WH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("ZH")>=0) :
            process+=selectDecayMode("h0",["h0->b,bbar;"])
            process+=selectDecayMode("Z0",["Z0->e-,e+;",
                                           "Z0->mu-,mu+;"])
            addedBRReweighter = True
            process+="insert SubProcess:MatrixElements[0] PowhegMEPP2ZH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("UE")>=0) :
            logging.error(" Process %s not supported for powheg matrix elements" % name)
            sys.exit(1)
        elif(parameterName.find("WZ")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="set /Herwig/Shower/ShowerHandler:SplitHardProcess No\n";
            process+="set /Herwig/Decays/ZDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/ZPowhegDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/WDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/WPowhegDecayer:PhotonGenerator NULL\n";
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SubProcess:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process WZ\n"
            process+=selectDecayMode("W+",["W+->nu_e,e+;",
                                           "W+->nu_mu,mu+;"])
            process+=selectDecayMode("W-",["W-->nu_ebar,e-;",
                                           "W-->nu_mubar,mu-;"])
            process+=selectDecayMode("Z0",["Z0->e-,e+;",
                                           "Z0->mu-,mu+;"])
            addedBRReweighter = True
            
        elif(parameterName.find("WW-emu")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="set /Herwig/Shower/ShowerHandler:SplitHardProcess No\n";
            process+="set /Herwig/Decays/ZDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/ZPowhegDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/WDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/WPowhegDecayer:PhotonGenerator NULL\n";
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SubProcess:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process WW\n"
            process+="set /Herwig/Particles/W+:Synchronized 0\n"
            process+="set /Herwig/Particles/W-:Synchronized 0\n"
            process+=selectDecayMode("W+",["W+->nu_e,e+;"])
            process+=selectDecayMode("W-",["W-->nu_mubar,mu-;"])
            addedBRReweighter = True
            
        elif(parameterName.find("WW-ll")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="set /Herwig/Shower/ShowerHandler:SplitHardProcess No\n";
            process+="set /Herwig/Decays/ZDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/ZPowhegDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/WDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/WPowhegDecayer:PhotonGenerator NULL\n";
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SubProcess:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process WW\n"
            process+=selectDecayMode("W+",["W+->nu_e,e+;",
                                           "W+->nu_mu,mu+;",
                                           "W+->nu_tau,tau+;"])
            addedBRReweighter = True
            
        elif(parameterName.find("ZZ-ll")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="set /Herwig/Shower/ShowerHandler:SplitHardProcess No\n";
            process+="set /Herwig/Decays/ZDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/ZPowhegDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/WDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/WPowhegDecayer:PhotonGenerator NULL\n";
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SubProcess:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process ZZ\n"
            process+=selectDecayMode("Z0",["Z0->e-,e+;",
                                           "Z0->mu-,mu+;",
                                           "Z0->tau-,tau+;"])
            addedBRReweighter = True
            
        elif(parameterName.find("ZZ-lv")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="set /Herwig/Shower/ShowerHandler:SplitHardProcess No\n";
            process+="set /Herwig/Decays/ZDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/ZPowhegDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/WDecayer:PhotonGenerator NULL\n";
            process+="set /Herwig/Decays/WPowhegDecayer:PhotonGenerator NULL\n";
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SubProcess:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process ZZ\n"
            process+=selectDecayMode("Z0",["Z0->e-,e+;",
                                           "Z0->mu-,mu+;",
                                           "Z0->tau-,tau+;",
                                           "Z0->nu_e,nu_ebar;",
                                           "Z0->nu_mu,nu_mubar;",
                                           "Z0->nu_tau,nu_taubar;"])
            addedBRReweighter = True
            
        elif(parameterName.find("W-Z-e")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2W2ff\nset PowhegMEqq2W2ff:Process Electron\n"
        elif(parameterName.find("W-Z-mu")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
            process+="insert SubProcess:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Muon\n"
        elif(parameterName.find("W-e")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2W2ff\nset PowhegMEqq2W2ff:Process Electron\n"
        elif(parameterName.find("W-mu")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2W2ff\nset PowhegMEqq2W2ff:Process Muon\n"
        elif(parameterName.find("Z-e")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Z-mu")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("Z-LowMass-e")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
            process+="set /Herwig/Cuts/Cuts:MHatMin 20.*GeV\nset /Herwig/Cuts/MassCut:MinM 20.*GeV\nset /Herwig/Cuts/MassCut:MaxM 70.*GeV\n"
        elif(parameterName.find("Z-MedMass-e")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
            process+="set /Herwig/Cuts/Cuts:MHatMin 40.*GeV\nset /Herwig/Cuts/MassCut:MinM 40.*GeV\nset /Herwig/Cuts/MassCut:MaxM 130.*GeV\n"
        elif(parameterName.find("Z-LowMass-mu")>=0) :
            process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
            process+="set /Herwig/Cuts/Cuts:MHatMin 10.*GeV\nset /Herwig/Cuts/MassCut:MinM 10.*GeV\nset /Herwig/Cuts/MassCut:MaxM 70.*GeV\n"
        elif(parameterName.find("Z-Mass1")>=0) :
            process+="set /Herwig/Cuts/Cuts:MHatMin 10.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MinM 10.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MaxM 35.*GeV\n"
            if(parameterName.find("-e")>=0) :
                process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
            else :
                process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("Z-Mass2")>=0) :
            process+="set /Herwig/Cuts/Cuts:MHatMin 25.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MinM 25.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MaxM 70.*GeV\n"
            if(parameterName.find("-e")>=0) :
                process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
            else :
                process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("Z-Mass3")>=0) :
            process+="set /Herwig/Cuts/Cuts:MHatMin 60.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MinM 60.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MaxM 120.*GeV\n"
            if(parameterName.find("-e")>=0) :
                process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
            else :
                process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("Z-Mass4")>=0) :
            process+="set /Herwig/Cuts/Cuts:MHatMin 110.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MinM 110.*GeV\n"
            process+="set /Herwig/Cuts/MassCut:MaxM 8000.*GeV\n"
            if(parameterName.find("-e")>=0) :
                process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
            else :
                process+="insert SubProcess:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("DiPhoton-GammaGamma")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEGammaGammaPowheg\n"
            process+="set MEGammaGammaPowheg:Process GammaGamma\n"
            process+="insert SubProcess:MatrixElements[0] MEGammaGamma\n"
            process+="set MEGammaGamma:Process gg\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaGamma","")
        elif(parameterName.find("DiPhoton-GammaJet")>=0) :
            process+="insert SubProcess:MatrixElements[0] MEGammaGammaPowheg\n"
            process+="set MEGammaGammaPowheg:Process VJet\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaJet","")
        else :
            logging.error(" Process %s not supported for internal POWHEG matrix elements" % name)
            sys.exit(1)
            
    elif( simulation=="Matchbox" or simulation=="Merging" ) :
        if(parameterName.find("8-VBF")>=0) :
            parameters["nlo"] = "read Matchbox/VBFNLO.in\n"
            if(simulation=="Merging"):
                process+="cd /Herwig/Merging/\n"
            process+="insert "+thefactory+":DiagramGenerator:RestrictLines 0 /Herwig/Particles/Z0\n"
            process+="insert "+thefactory+":DiagramGenerator:RestrictLines 0 /Herwig/Particles/W+\n"
            process+="insert "+thefactory+":DiagramGenerator:RestrictLines 0 /Herwig/Particles/W-\n"
            process+="insert "+thefactory+":DiagramGenerator:RestrictLines 0 /Herwig/Particles/gamma\n"
            process+="do "+thefactory+":DiagramGenerator:TimeLikeRange 0 0\n"
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p h0 j j","0","3","FixedScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p p h0 j j","0","3","FixedScale",1,1)
            process+=setHardProcessWidthToZero(["h0"])
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
            if(parameterName.find("GammaGamma")>=0) :
               process+=selectDecayMode("h0",["h0->gamma,gamma;"])
               process+=addBRReweighter()
               
        elif(parameterName.find("VBF")>=0) :
            process+=selectDecayMode("h0",["h0->tau-,tau+;"])
            process+=addBRReweighter()
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            parameters["nlo"] = "read Matchbox/VBFNLO.in\n"
            if(simulation=="Merging"):
                process+="cd /Herwig/Merging/\n"
            process+="insert "+thefactory+":DiagramGenerator:RestrictLines 0 /Herwig/Particles/Z0\n"
            process+="insert "+thefactory+":DiagramGenerator:RestrictLines 0 /Herwig/Particles/W+\n"
            process+="insert "+thefactory+":DiagramGenerator:RestrictLines 0 /Herwig/Particles/W-\n"
            process+="insert "+thefactory+":DiagramGenerator:RestrictLines 0 /Herwig/Particles/gamma\n"
            process+="do "+thefactory+":DiagramGenerator:TimeLikeRange 0 0\n"
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p h0 j j","0","3","FixedScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p p h0 j j","0","3","FixedScale",1,1)
            process+=setHardProcessWidthToZero(["h0"])
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
        elif(parameterName.find("ggHJet")>=0) :
            if(simulation=="Merging"):
               logging.warning("ggHJet not explicitly tested for %s " % simulation)
               sys.exit(0)
            parameters["nlo"] = "read Matchbox/MadGraph-GoSam.in\nread Matchbox/HiggsEffective.in\n"
            process+=selectDecayMode("h0",["h0->tau-,tau+;"])
            process+=addBRReweighter()
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+=setHardProcessWidthToZero(["h0"])
            process+=addProcess(thefactory,"p p h0 j","3","1","FixedScale",0,0)
            process+=addFirstJet("20")
            process+="set "+thefactory+":ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
        elif(parameterName.find("8-ggH")>=0) :
            parameters["nlo"] = "read Matchbox/MadGraph-GoSam.in\nread Matchbox/HiggsEffective.in\n"
            if(simulation=="Merging"):
                process+= "cd /Herwig/MatrixElements/Matchbox/Amplitudes\nset OpenLoops:HiggsEff On\nset MadGraph:Model heft\n"
                process+="cd /Herwig/Merging/\n"
            process+=setHardProcessWidthToZero(["h0"])
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p h0","2","1","FixedScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p p h0","2","1","FixedScale",2,2)
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
            if(parameterName.find("GammaGamma")>=0) :
               process+=selectDecayMode("h0",["h0->gamma,gamma;"])
               process+=addBRReweighter()
               
        elif(parameterName.find("ggH")>=0) :
            parameters["nlo"] = "read Matchbox/MadGraph-GoSam.in\nread Matchbox/HiggsEffective.in\n"
            if(simulation=="Merging"):
                process+= "cd /Herwig/MatrixElements/Matchbox/Amplitudes\nset OpenLoops:HiggsEff On\nset MadGraph:Model heft\n"
                process+="cd /Herwig/Merging/\n"
            process+=selectDecayMode("h0",["h0->tau-,tau+;"])
            process+=addBRReweighter()
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+=setHardProcessWidthToZero(["h0"])
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p h0","2","1","FixedScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p p h0","2","1","FixedScale",2,2)

            process+="set "+thefactory+":ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
        elif(parameterName.find("8-WH")>=0) :
            if(simulation=="Merging"):
              logging.warning("8-WH not explicitly tested for %s " % simulation)
              sys.exit(0)
            process+=setHardProcessWidthToZero(["h0","W+","W-"])
            process+=addProcess(thefactory,"p p W+ h0","0","2","FixedScale",0,0)
            process+=addProcess(thefactory,"p p W- h0","0","2","FixedScale",0,0)
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
            if(parameterName.find("GammaGamma")>=0) :
               process+=selectDecayMode("h0",["h0->gamma,gamma;"])
               process+=addBRReweighter()
               
        elif(parameterName.find("8-ZH")>=0) :
            if(simulation=="Merging"):
              logging.warning("8-ZH not explicitly tested for %s " % simulation)
              sys.exit(0)
            process+=setHardProcessWidthToZero(["h0","Z0"])
            process+=addProcess(thefactory,"p p Z0 h0","0","2","FixedScale",0,0)
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
            if(parameterName.find("GammaGamma")>=0) :
               process+=selectDecayMode("h0",["h0->gamma,gamma;"])
               process+=addBRReweighter()
               
        elif(parameterName.find("WH")>=0) :
            if(simulation=="Merging"):
              logging.warning("WH not explicitly tested for %s " % simulation)
              sys.exit(0)
            process+=selectDecayMode("h0",["h0->b,bbar;"])
            process+=addBRReweighter()
            process+=setHardProcessWidthToZero(["h0"])
            process+=addProcess(thefactory,"p p e+ nu h0","0","3","LeptonPairMassScale",0,0)
            process+=addProcess(thefactory,"p p e- nu h0","0","3","LeptonPairMassScale",0,0)
            process+=addProcess(thefactory,"p p mu+ nu h0","0","3","LeptonPairMassScale",0,0)
            process+=addProcess(thefactory,"p p mu- nu h0","0","3","LeptonPairMassScale",0,0)
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("ZH")>=0) :
            if(simulation=="Merging"):
              logging.warning("ZH not explicitly tested for %s " % simulation)
              sys.exit(0)
            process+=selectDecayMode("h0",["h0->b,bbar;"])
            process+=addBRReweighter()
            process+=setHardProcessWidthToZero(["h0"])
            process+=addProcess(thefactory,"p p e+ e- h0","0","3","LeptonPairMassScale",0,0)
            process+=addProcess(thefactory,"p p mu+ mu- h0","0","3","LeptonPairMassScale",0,0)
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("UE")>=0) :
            logging.error(" Process %s not supported for Matchbox matrix elements" % name)
            sys.exit(1)
        elif(parameterName.find("7-DiJets")>=0 or parameterName.find("8-DiJets")>=0) :
            if(simulation=="Matchbox"):
              process+=addProcess(thefactory,"p p j j","2","0","MaxJetPtScale",0,0)
            elif(simulation=="Merging"):
              process+=addProcess(thefactory,"p p j j","2","0","MaxJetPtScale",1,1)
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("-A")>=0) :
                 process+=addFirstJet("45")
                 process+=addSecondJet("25")
                 process+="set /Herwig/Cuts/FirstJet:YRange  -3. 3.\n"
                 process+="set /Herwig/Cuts/SecondJet:YRange -3. 3.\n"
            elif(parameterName.find("-B")>=0) :
                 process+=addFirstJet("20")
                 process+=addSecondJet("15")
                 process+="set /Herwig/Cuts/FirstJet:YRange  -2.7 2.7\n"
                 process+="set /Herwig/Cuts/SecondJet:YRange -2.7 2.7\n"
            elif(parameterName.find("-C")>=0) :
                 process+=addFirstJet("20")
                 process+=addSecondJet("15")
                 process+="set /Herwig/Cuts/FirstJet:YRange  -4.8 4.8\n"
                 process+="set /Herwig/Cuts/SecondJet:YRange -4.8 4.8\n"
            else :
                 logging.error("Exit 00001")
                 sys.exit(1)
            if(parameterName.find("DiJets-1")>=0) :
                process+=addJetPairCut("90")
            elif(parameterName.find("DiJets-2")>=0) :
                process+=addJetPairCut("200")
            elif(parameterName.find("DiJets-3")>=0) :
                process+=addJetPairCut("450")
            elif(parameterName.find("DiJets-4")>=0) :
                process+=addJetPairCut("750")
            elif(parameterName.find("DiJets-5")>=0) :
                process+=addJetPairCut("950")
            elif(parameterName.find("DiJets-6")>=0) :
                process+=addJetPairCut("1550")
            elif(parameterName.find("DiJets-7")>=0) :
                process+=addJetPairCut("2150")
            elif(parameterName.find("DiJets-8")>=0) :
                process+=addJetPairCut("2750")
            else :
                logging.error("Exit 00002")
                sys.exit(1)


        elif(parameterName.find("7-Jets")>=0 or parameterName.find("8-Jets")>=0 or \
             parameterName.find("13-Jets")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p j j","2","0","MaxJetPtScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p p j j","2","0","MaxJetPtScale",1,1)
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("Jets-10")>=0) :
                process+=addFirstJet("1800")
            elif(parameterName.find("Jets-0")>=0) :
                process+=addFirstJet("5")
            elif(parameterName.find("Jets-1")>=0) :
                process+=addFirstJet("10")
            elif(parameterName.find("Jets-2")>=0) :
                process+=addFirstJet("20")
            elif(parameterName.find("Jets-3")>=0) :
                process+=addFirstJet("40")
            elif(parameterName.find("Jets-4")>=0) :
                process+=addFirstJet("70")
            elif(parameterName.find("Jets-5")>=0) :
                process+=addFirstJet("150")
            elif(parameterName.find("Jets-6")>=0) :
                process+=addFirstJet("200")
            elif(parameterName.find("Jets-7")>=0) :
                process+=addFirstJet("300")
            elif(parameterName.find("Jets-8")>=0) :
                process+=addFirstJet("500")
            elif(parameterName.find("Jets-9")>=0) :
                process+=addFirstJet("800")
            else :
                logging.error("Exit 00003")
                sys.exit(1)
        elif(parameterName.find("7-Charm")>=0 or \
             parameterName.find("7-Bottom")>=0) :
            parameters["bscheme"]=fourFlavour
            process+="set /Herwig/Particles/b:HardProcessMass 4.2*GeV\n"
            process+="set /Herwig/Particles/bbar:HardProcessMass 4.2*GeV\n"
            
            
            if(parameterName.find("7-Bottom")>=0) :
                if(simulation=="Matchbox"):
                    process+=addProcess(thefactory,"p p b bbar","2","0","MaxJetPtScale",0,0)
                elif(simulation=="Merging"):
                    process+=addProcess(thefactory,"p p b bbar","2","0","MaxJetPtScale",1,0)
            else:
                if(simulation=="Matchbox"):
                    process+=addProcess(thefactory,"p p c cbar","2","0","MaxJetPtScale",0,0)
                elif(simulation=="Merging"):
                    process+=addProcess(thefactory,"p p c cbar","2","0","MaxJetPtScale",1,0)

            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("-0")>=0) :
                process+=addFirstJet("0")
            elif(parameterName.find("-1")>=0) :
                process+=addFirstJet("5")
            elif(parameterName.find("-2")>=0) :
                process+=addFirstJet("20")
            elif(parameterName.find("-3")>=0) :
                process+=addFirstJet("50")
            elif(parameterName.find("-4")>=0) :
                process+=addFirstJet("80")
            elif(parameterName.find("-5")>=0) :
                process+=addFirstJet("110")
            elif(parameterName.find("-6")>=0) :
                process+=addFirstJet("30")
                process+=addSecondJet("25")
                process+=addJetPairCut("90")
            elif(parameterName.find("-7")>=0) :
                process+=addFirstJet("30")
                process+=addSecondJet("25")
                process+=addJetPairCut("340")
            elif(parameterName.find("-8")>=0) :
                process+=addFirstJet("30")
                process+=addSecondJet("25")
                process+=addJetPairCut("500")
            else :
                logging.error("Exit 00004")
                sys.exit(1)
                  
        elif(parameterName.find("Top-L")>=0) :
            process+=setHardProcessWidthToZero(["t","tbar"])
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p t tbar","2","0","TopPairMTScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p p t tbar","2","0","TopPairMTScale",2,2)
            process+=selectDecayMode("t",["t->nu_e,e+,b;",
                                          "t->nu_mu,mu+,b;"])
            process+=addBRReweighter()
            
        elif(parameterName.find("Top-SL")>=0) :
            process+=setHardProcessWidthToZero(["t","tbar"])
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p t tbar","2","0","TopPairMTScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p p t tbar","2","0","TopPairMTScale",2,2)
            process+="set /Herwig/Particles/t:Synchronized Not_synchronized\n"
            process+="set /Herwig/Particles/tbar:Synchronized Not_synchronized\n"
            process+=selectDecayMode("t",["t->nu_e,e+,b;",
                                          "t->nu_mu,mu+,b;"])
            process+=selectDecayMode("tbar",["tbar->b,bbar,cbar;",
                                             "tbar->bbar,cbar,d;",
                                             "tbar->bbar,cbar,s;",
                                             "tbar->bbar,s,ubar;",
                                             "tbar->bbar,ubar,d;"])
            process+=addBRReweighter()
            
        elif(parameterName.find("Top-All")>=0) :
            process+=setHardProcessWidthToZero(["t","tbar"])
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p t tbar","2","0","TopPairMTScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p p t tbar","2","0","TopPairMTScale",2,2)
        elif(parameterName.find("WZ")>=0) :
            if(simulation=="Merging"):
              logging.warning("WZ not explicitly tested for %s " % simulation)
              sys.exit(0)
            process+=setHardProcessWidthToZero(["W+","W-","Z0"])
            process+=addProcess(thefactory,"p p W+ Z0","0","2","FixedScale",0,0)
            process+=addProcess(thefactory,"p p W- Z0","0","2","FixedScale",0,0)
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 171.6*GeV\n\n"
            process+=selectDecayMode("W+",["W+->nu_e,e+;",
                                           "W+->nu_mu,mu+;"])
            process+=selectDecayMode("W-",["W-->nu_ebar,e-;",
                                           "W-->nu_mubar,mu-;"])
            process+=selectDecayMode("Z0",["Z0->e-,e+;",
                                           "Z0->mu-,mu+;"])
            process+=addBRReweighter()
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("WW-emu")>=0) :
            if(simulation=="Merging"):
              logging.warning("WW-emu not explicitly tested for %s " % simulation)
              sys.exit(0)
            
            process+=setHardProcessWidthToZero(["W+","W-","Z0"])
            process+=addProcess(thefactory,"p p W+ W-","0","2","FixedScale",0,0)
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 160.8*GeV\n"
            process+="set /Herwig/Particles/W+:Synchronized 0\n"
            process+="set /Herwig/Particles/W-:Synchronized 0\n"
            process+=selectDecayMode("W+",["W+->nu_e,e+;"])
            process+=selectDecayMode("W-",["W-->nu_mubar,mu-;"])
            process+=addBRReweighter()
            parameters["bscheme"] = "read Matchbox/FourFlavourScheme.in\n"
            
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("WW-ll")>=0) :
            if(simulation=="Merging"):
              logging.warning("WW-ll not explicitly tested for %s " % simulation)
              sys.exit(0)
            process+=setHardProcessWidthToZero(["W+","W-","Z0"])
            process+=addProcess(thefactory,"p p W+ W-","0","2","FixedScale",0,0)
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 160.8*GeV\n"
            process+=selectDecayMode("W+",["W+->nu_e,e+;",
                                           "W+->nu_mu,mu+;",
                                           "W+->nu_tau,tau+;"])
            process+=addBRReweighter()
            process+=addLeptonPairCut("60","120")
            parameters["bscheme"] = "read Matchbox/FourFlavourScheme.in\n"

        elif(parameterName.find("ZZ-ll")>=0) :
            if(simulation=="Merging"):
              logging.warning("ZZ-ll not explicitly tested for %s " % simulation)
              sys.exit(0)
            process+=setHardProcessWidthToZero(["W+","W-","Z0"])
            process+=addProcess(thefactory,"p p Z0 Z0","0","2","FixedScale",0,0)
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 182.2*GeV\n"
            process+=selectDecayMode("Z0",["Z0->e-,e+;",
                                           "Z0->mu-,mu+;",
                                           "Z0->tau-,tau+;"])
            process+=addBRReweighter()
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("ZZ-lv")>=0) :
            if(simulation=="Merging"):
              logging.warning("ZZ-lv not explicitly tested for %s " % simulation)
              sys.exit(0)
            process+=setHardProcessWidthToZero(["W+","W-","Z0"])
            process+=addProcess(thefactory,"p p Z0 Z0","0","2","FixedScale",0,0)
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 182.2*GeV\n"
            process+=selectDecayMode("Z0",["Z0->e-,e+;",
                                           "Z0->mu-,mu+;",
                                           "Z0->tau-,tau+;",
                                           "Z0->nu_e,nu_ebar;",
                                           "Z0->nu_mu,nu_mubar;",
                                           "Z0->nu_tau,nu_taubar;"])
            process+=addBRReweighter()
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("W-Z-e")>=0) :
            if(simulation=="Matchbox"):
              process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",0,0)
              process+=addProcess(thefactory,"p p e+ nu","0","2","LeptonPairMassScale",0,0)
              process+=addProcess(thefactory,"p p e- nu","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
              process+="do "+thefactory+":StartParticleGroup epm\n"
              process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/e+\n"
              process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/e-\n"
              process+="do "+thefactory+":EndParticleGroup\n"
              process+="do "+thefactory+":StartParticleGroup epmnu\n"
              process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/e+\n"
              process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/e-\n"
              process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/nu_e\n"
              process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/nu_ebar\n"
              process+="do "+thefactory+":EndParticleGroup\n"
              process+=addProcess(thefactory,"p p epm epmnu","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("W-Z-mu")>=0) :
            if(simulation=="Matchbox"):
              process+=addProcess(thefactory,"p p mu+ mu-","0","2","LeptonPairMassScale",0,0)
              process+=addProcess(thefactory,"p p mu+ nu","0","2","LeptonPairMassScale",0,0)
              process+=addProcess(thefactory,"p p mu- nu","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
              process+="do "+thefactory+":StartParticleGroup mupm\n"
              process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/mu+\n"
              process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/mu-\n"
              process+="do "+thefactory+":EndParticleGroup\n"
              process+="do "+thefactory+":StartParticleGroup mupmnu\n"
              process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/mu+\n"
              process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/mu-\n"
              process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/nu_mu\n"
              process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/nu_mubar\n"
              process+="do "+thefactory+":EndParticleGroup\n"
              process+=addProcess(thefactory,"p p mupm mupmnu","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("W-e")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p e+ nu","0","2","LeptonPairMassScale",0,0)
                process+=addProcess(thefactory,"p p e- nu","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
                process+="do "+thefactory+":StartParticleGroup epm\n"
                process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/e+\n"
                process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/e-\n"
                process+="do "+thefactory+":EndParticleGroup\n"
                process+=addProcess(thefactory,"p p epm nu","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("60","120")

        elif(parameterName.find("W-mu")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p mu+ nu","0","2","LeptonPairMassScale",0,0)
                process+=addProcess(thefactory,"p p mu- nu","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
                process+="do "+thefactory+":StartParticleGroup mupm\n"
                process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/mu+\n"
                process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/mu-\n"
                process+="do "+thefactory+":EndParticleGroup\n"
                process+=addProcess(thefactory,"p p mupm nu","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("60","120")

        elif(parameterName.find("Z-e")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("Z-mu")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p mu+ mu-","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p p mu+ mu-","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("Z-jj")>=0) :
            if(simulation=="Merging"):
              logging.warning("Z-jj not explicitly tested for %s " % simulation)
              sys.exit(0)
            process+=addProcess(thefactory,"p p e+ e- j j","2","2","LeptonPairMassScale",0,0)
            process+=addFirstJet("40")
            process+=addSecondJet("30")
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("Z-LowMass-e")>=0) :
            if(simulation=="Matchbox"):
               process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
               process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("20","70")
        elif(parameterName.find("Z-MedMass-e")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("40","130")
        elif(parameterName.find("Z-LowMass-mu")>=0) :
            if(simulation=="Matchbox"):
                process+=addProcess(thefactory,"p p mu+ mu-","0","2","LeptonPairMassScale",0,0)
            elif(simulation=="Merging"):
                process+=addProcess(thefactory,"p p mu+ mu-","0","2","LeptonPairMassScale",2,2)
            process+=addLeptonPairCut("10","70")
        elif(parameterName.find("Z-Mass1")>=0) :
            process+=addLeptonPairCut("10","35")
            if(parameterName.find("-e")>=0) :
                if(simulation=="Matchbox"):
                    process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",0,0)
                elif(simulation=="Merging"):
                    process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",2,2)
            else :
                if(simulation=="Matchbox"):
                    process+=addProcess(thefactory,"p p mu+ mu-","0","2","LeptonPairMassScale",0,0)
                elif(simulation=="Merging"):
                    process+=addProcess(thefactory,"p p mu+ mu-","0","2","LeptonPairMassScale",2,2)
        elif(parameterName.find("Z-Mass2")>=0) :
            process+=addLeptonPairCut("25","70")
            if(parameterName.find("-e")>=0) :
                if(simulation=="Matchbox"):
                    process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",0,0)
                elif(simulation=="Merging"):
                    process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",2,2)
            else :
                if(simulation=="Matchbox"):
                    process+=addProcess(thefactory,"p p mu+ mu-","0","2","LeptonPairMassScale",0,0)
                elif(simulation=="Merging"):
                    process+=addProcess(thefactory,"p p mu+ mu-","0","2","LeptonPairMassScale",2,2)
        elif(parameterName.find("Z-Mass3")>=0) :
            process+=addLeptonPairCut("60","120")
            if(parameterName.find("-e")>=0) :
                if(simulation=="Matchbox"):
                    process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",0,0)
                elif(simulation=="Merging"):
                    process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",2,2)
            else :
                if(simulation=="Matchbox"):
                    process+=addProcess(thefactory,"p p mu+ mu-","0","2","LeptonPairMassScale",0,0)
                elif(simulation=="Merging"):
                    process+=addProcess(thefactory,"p p mu+ mu-","0","2","LeptonPairMassScale",2,2)
        elif(parameterName.find("Z-Mass4")>=0) :
            process+=addLeptonPairCut("115","8000")
            if(parameterName.find("-e")>=0) :
                if(simulation=="Matchbox"):
                  process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",0,0)
                elif(simulation=="Merging"):
                  process+=addProcess(thefactory,"p p e+ e-","0","2","LeptonPairMassScale",2,2)
            else :
                if(simulation=="Matchbox"):
                  process+=addProcess(thefactory,"p p mu+ mu-","0","2","LeptonPairMassScale",0,0)
                elif(simulation=="Merging"):
                  process+=addProcess(thefactory,"p p mu+ mu-","0","2","LeptonPairMassScale",2,2)
        elif(parameterName.find("W-Jet")>=0) :
            if(simulation=="Merging"):
              logging.warning("W-Jet not explicitly tested for %s " % simulation)
              sys.exit(0)
            
            process+=addProcess(thefactory,"p p e+ nu j","1","2","HTScale",0,0)
            process+=addProcess(thefactory,"p p e- nu j","1","2","HTScale",0,0)
            
            process+=addLeptonPairCut("60","120")
            if(parameterName.find("W-Jet-1-e")>=0) :
                process+=addFirstJet("100")
                parameterName=parameterName.replace("W-Jet-1-e","W-Jet-e")
            elif(parameterName.find("W-Jet-2-e")>=0) :
                process+=addFirstJet("190")
                parameterName=parameterName.replace("W-Jet-2-e","W-Jet-e")
            elif(parameterName.find("W-Jet-3-e")>=0) :
                process+=addFirstJet("270")
                parameterName=parameterName.replace("W-Jet-3-e","W-Jet-e")
            else :
                logging.error("Exit 00005")
                sys.exit(1)
        elif(parameterName.find("Z-Jet")>=0) :
            if(simulation=="Merging"):
              logging.warning("Z-Jet not explicitly tested for %s " % simulation)
              sys.exit(0)
            
            
            if(parameterName.find("-e")>=0) :
                process+=addProcess(thefactory,"p p e+ e- j","1","2","HTScale",0,0)
                if(parameterName.find("Z-Jet-0-e")>=0) :
                    process+=addFirstJet("35")
                    parameterName=parameterName.replace("Z-Jet-0-e","Z-Jet-e")
                elif(parameterName.find("Z-Jet-1-e")>=0) :
                    process+=addFirstJet("100")
                    parameterName=parameterName.replace("Z-Jet-1-e","Z-Jet-e")
                elif(parameterName.find("Z-Jet-2-e")>=0) :
                    process+=addFirstJet("190")
                    parameterName=parameterName.replace("Z-Jet-2-e","Z-Jet-e")
                elif(parameterName.find("Z-Jet-3-e")>=0) :
                    process+=addFirstJet("270")
                    parameterName=parameterName.replace("Z-Jet-3-e","Z-Jet-e")
                else :
                    logging.error("Exit 00006")
                    sys.exit(1)
            else :
                process+=addProcess(thefactory,"p p mu+ mu- j","1","2","HTScale",0,0)
                process+=addFirstJet("35")
                parameterName=parameterName.replace("Z-Jet-0-mu","Z-Jet-mu")
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("Z-bb")>=0) :
            if(simulation=="Merging"):
              logging.warning("Z-bb not explicitly tested for %s " % simulation)
              sys.exit(0)
            parameters["bscheme"]=fourFlavour
            process+="set /Herwig/Particles/b:HardProcessMass 4.2*GeV\nset /Herwig/Particles/bbar:HardProcessMass 4.2*GeV\n"
            process+=addProcess(thefactory,"p p e+ e- b bbar","2","2","FixedScale",0,0)
            process+=addLeptonPairCut("66","116")
            process+=addFirstJet("18")
            process+=addSecondJet("15")
            process+=addLeptonPairCut("60","120")
        elif(parameterName.find("Z-b")>=0) :
            if(simulation=="Merging"):
              logging.warning("Z-b not explicitly tested for %s " % simulation)
              sys.exit(0)
            process+="do "+thefactory+":StartParticleGroup bjet\n"
            process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/b\n"
            process+="insert "+thefactory+":ParticleGroup 0 /Herwig/Particles/bbar\n"
            process+="do "+thefactory+":EndParticleGroup\n"
            process+=addProcess(thefactory,"p p e+ e- bjet","1","2","FixedScale",0,0)
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 91.2*GeV\n"
            process+=addLeptonPairCut("60","120")
            process+=addFirstJet("15")
        elif(parameterName.find("W-b")>=0) :
            if(simulation=="Merging"):
              logging.warning("W-b not explicitly tested for %s " % simulation)
              sys.exit(0)
            parameters["bscheme"]=fourFlavour
            process += "set /Herwig/Particles/b:HardProcessMass 4.2*GeV\nset /Herwig/Particles/bbar:HardProcessMass 4.2*GeV\n"
            process+=addProcess(thefactory,"p p e-  nu b bbar","2","2","FixedScale",0,0)
            process+=addProcess(thefactory,"p p mu+ nu b bbar","2","2","FixedScale",0,0)
            process += "set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 80.4*GeV\n"
            process+=addFirstJet("30")
            process+=addLeptonPairCut("60","120")
        else :
            logging.error(" Process %s not supported for Matchbox matrix elements" % name)
            sys.exit(1)
# LHC-GammaGamma
elif(collider=="LHC-GammaGamma" ) :
    if(parameterName.find("-7-")>=0) :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 7000.0\n"
    elif(parameterName.find("-8-")>=0) :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 8000.0\n"
    else :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 7000.0\n"
    if(simulation=="") :
        if(parameterName.find("7")>=0) :
            process += "insert SubProcess:MatrixElements 0 /Herwig/MatrixElements/MEgg2ff\n"
            process += "set /Herwig/MatrixElements/MEgg2ff:Process Muon\n"
        else :
            logging.error(" Process %s not supported for default matrix elements" % name)
            sys.exit(1)
    else :
        logging.error("LHC-GammaGamma not supported for %s " % simulation)
        sys.exit(1)


parameters['parameterFile'] = os.path.join(collider,collider+"-"+parameterName+".in")
parameters['runname'] = name
parameters['process'] = process


#check if selecteddecaymode and addedBRReweighter is consistent

if selecteddecaymode and not addedBRReweighter:
    logging.error("Decaymode was selected but no BRReweighter was added.")
    sys.exit(1)

if addedBRReweighter and not selecteddecaymode:
    logging.error("BRReweighter was added but no Decaymode was selected.")
    sys.exit(1)

# check that we only add one process if in merging mode:

if numberOfAddedProcesses > 1 and simulation =="Merging":
    logging.error("In Merging only one process is allowed at the moment. See ticket #403.")
    sys.exit(1)

# Check if a process was added for Merging or Matchbox:

if numberOfAddedProcesses == 0 and (simulation =="Merging" or simulation =="Matchbox"):
    logging.error("No process was selected.")
    sys.exit(1)

# write the file

with open(os.path.join("Rivet",name+".in") ,'w') as f:
        f.write( template.substitute(parameters))






