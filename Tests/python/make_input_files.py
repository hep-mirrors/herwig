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
simulation=""
istart = 1
print name
if(name.find("Matchbox-Powheg")>0) :
    istart = 3
    simulation="Matchbox"
    parameters["shower"] = "read Matchbox/Powheg-DefaultShower.in\n"
elif(name.find("Matchbox")>0) :
    istart = 2
    simulation="Matchbox"
    parameters["shower"] = "read Matchbox/MCatNLO-DefaultShower.in\n"
elif(name.find("Dipole")>0) :
    istart = 2
    simulation="Matchbox"
    parameters["shower"]  = "read Matchbox/MCatNLO-DipoleShower.in\n"
elif(name.find("Powheg")>0) :
    istart = 2
    simulation="Powheg"

if(simulation=="Matchbox") :
    parameters["bscheme"] = "read Matchbox/FiveFlavourScheme.in\n"
    if(parameters["shower"].find("Dipole")>=0) :
        parameters["bscheme"] += "read Matchbox/FiveFlavourNoBMassScheme.in\n"
    if(collider.find("DIS")<0) :
        parameters["nlo"] = "read Matchbox/MadGraph-OpenLoops.in\n"

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
        process = "set /Herwig/MatrixElements/MEee2gZ2qq:MaximumFlavour 4"
        if(parameterName=="10.58") :
            process += "\ncreate Herwig::MEee2VectorMeson /Herwig/MatrixElements/MEUpsilon HwMELepton.so\nset /Herwig/MatrixElements/MEUpsilon:VectorMeson /Herwig/Particles/Upsilon(4S)\nset /Herwig/MatrixElements/MEUpsilon:Coupling 0.0004151809\ninsert /Herwig/MatrixElements/SimpleEE:MatrixElements 0 /Herwig/MatrixElements/MEUpsilon"
    elif(simulation=="Powheg") :
        process = "set /Herwig/MatrixElements/PowhegMEee2gZ2qq:MaximumFlavour 4"
        if(parameterName=="10.58") :
            process += "\ncreate Herwig::MEee2VectorMeson /Herwig/MatrixElements/MEUpsilon HwMELepton.so\nset /Herwig/MatrixElements/MEUpsilon:VectorMeson /Herwig/Particles/Upsilon(4S)\nset /Herwig/MatrixElements/MEUpsilon:Coupling 0.0004151809\ninsert /Herwig/MatrixElements/SimpleEE:MatrixElements 0 /Herwig/MatrixElements/MEUpsilon"
    elif(simulation=="Matchbox" ) :
        process = "do Factory:Process e- e+ -> u ubar\ndo Factory:Process e- e+ -> d dbar\ndo Factory:Process e- e+ -> c cbar\ndo Factory:Process e- e+ -> s sbar"
        if(parameterName=="10.58") :
            process += "\ninsert /Herwig/Generators/EventGenerator:EventHandler:SubProcessHandlers 0 /Herwig/MatrixElements/SimpleEE\ncreate Herwig::MEee2VectorMeson /Herwig/MatrixElements/MEUpsilon HwMELepton.so\nset /Herwig/MatrixElements/MEUpsilon:VectorMeson /Herwig/Particles/Upsilon(4S)\nset /Herwig/MatrixElements/MEUpsilon:Coupling 0.0004151809\ninsert /Herwig/MatrixElements/SimpleEE:MatrixElements 0 /Herwig/MatrixElements/MEUpsilon\n"
# DIS
elif(collider=="DIS") :
    if(simulation=="") :
        if(parameterName.find("NoME")>=0) :
            process = "set /Herwig/Shower/Evolver:MECorrMode 0"
            parameterName=parameterName.replace("NoME-","")
        else :
            process = ""
    elif(simulation=="Powheg") :
        process = ""
    elif(simulation=="Matchbox" ) :
        if(parameterName.find("e-")>=0) :
            process="do Factory:Process e- p -> e- j"
        else :
            process="do Factory:Process e+ p -> e+ j"
# LEP
elif(collider=="LEP") :
    if(simulation=="") :
        process=""
        if(parameterName=="10") :
            process="set /Herwig/MatrixElements/MEee2gZ2qq:MaximumFlavour 4"
    elif(simulation=="Powheg") :
        process=""
        if(parameterName=="10") :
            process="set /Herwig/MatrixElements/PowhegMEee2gZ2qq:MaximumFlavour 4"
    elif(simulation=="Matchbox" ) :
        if(parameterName=="10") :
            process="do Factory:Process e- e+ -> u ubar\ndo Factory:Process e- e+ -> d dbar\ndo Factory:Process e- e+ -> c cbar\ndo Factory:Process e- e+ -> s sbar"
        else :
            process="do Factory:Process e- e+ -> j j"
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
            process+="insert SimpleQCD:MatrixElements[0] MEGammaJet\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 15.\n"
        elif(parameterName.find("DiPhoton-GammaGamma")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEGammaGamma\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaGamma","")
        elif(parameterName.find("DiPhoton-GammaJet")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEGammaJet\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaJet","")
        elif(parameterName.find("UE")>=0) :
            process += "insert SimpleQCD:MatrixElements[0] MEMinBias\n"
            process += "set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            process += "set /Herwig/Generators/EventGenerator:EventHandler:Cuts /Herwig/Cuts/MinBiasCuts\n"
            process += "create Herwig::MPIXSecReweighter /Herwig/Generators/MPIXSecReweighter\n"
            process += "insert /Herwig/Generators/EventGenerator:EventHandler:PostSubProcessHandlers 0 /Herwig/Generators/MPIXSecReweighter\n"
            process += "set /Herwig/Decays/DecayHandler:LifeTimeOption 0\n"
            process += "set /Herwig/Decays/DecayHandler:MaxLifeTime 10*mm\n"
        elif(parameterName.find("Jets")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEQCD2to2\n"
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("Run-II-Jets-10")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 500.*GeV\n"
            elif(parameterName.find("Run-II-Jets-11")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 900.*GeV\n"
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
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 100.*GeV\n"
            elif(parameterName.find("Run-I-Jets-7")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 400.*GeV\n"
            elif(parameterName.find("Run-I-Jets-8")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 700.*GeV\n"
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
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 100.*GeV\n"
            elif(parameterName.find("Run-II-Jets-9")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 300.*GeV\n"
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
            process+="insert SimpleQCD:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Electron\ninsert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Run-I-W")>=0 or parameterName.find("Run-II-W")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Electron\n"
        elif(parameterName.find("Run-I-Z")>=0 or parameterName.find("Run-II-Z-e")>=0) :
            process +="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Run-II-Z-LowMass-mu")>=0) :
            process +="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 25*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 70*GeV\n"
        elif(parameterName.find("Run-II-Z-HighMass-mu")>=0) :
            process +="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 150*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 600*GeV\n"
        elif(parameterName.find("Run-II-Z-mu")>=0) :
            process +="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
    elif(simulation=="Powheg") :
        if(parameterName.find("Run-I-WZ")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2W2ff\nset PowhegMEqq2W2ff:Process Electron\ninsert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Run-I-W")>=0 or parameterName.find("Run-II-W")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2W2ff\nset PowhegMEqq2W2ff:Process Electron\n"
        elif(parameterName.find("Run-I-Z")>=0 or parameterName.find("Run-II-Z-e")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Run-II-Z-LowMass-mu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 25*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 70*GeV\n"
        elif(parameterName.find("Run-II-Z-HighMass-mu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 150*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 600*GeV\n"
        elif(parameterName.find("Run-II-Z-mu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("DiPhoton-GammaGamma")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEGammaGammaPowheg\n"
            process+="set MEGammaGammaPowheg:Process GammaGamma\n"
            process+="insert SimpleQCD:MatrixElements[0] MEGammaGamma\n"
            process+="set MEGammaGamma:Process gg\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaGamma","")
        elif(parameterName.find("DiPhoton-GammaJet")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEGammaGammaPowheg\n"
            process+="set MEGammaGammaPowheg:Process VJet\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaJet","")
    elif(simulation=="Matchbox" ) :
        if(parameterName.find("Jets")>=0) :
            process+="set Factory:OrderInAlphaS 2\nset Factory:OrderInAlphaEW 0\n"
            process+="do Factory:Process p p j j\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/MaxJetPtScale\n"
            process+="set  /Herwig/Cuts/Cuts:JetFinder  /Herwig/Cuts/JetFinder\n"
            process+="insert  /Herwig/Cuts/Cuts:MultiCuts 0  /Herwig/Cuts/JetCuts\n"
            process+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/FirstJet\n"
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("Run-II-Jets-10")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 30.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 25.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 500.*GeV\n"
            elif(parameterName.find("Run-II-Jets-11")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 30.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 25.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 900.*GeV\n"
            elif(parameterName.find("Run-II-Jets-12")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 30.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 25.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 300.*GeV\n"
            elif(parameterName.find("Run-I-Jets-1")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 20.\n"
            elif(parameterName.find("Run-I-Jets-2")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 40.\n"
            elif(parameterName.find("Run-I-Jets-3")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 65.\n"
            elif(parameterName.find("Run-I-Jets-4")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 90.\n"
            elif(parameterName.find("Run-I-Jets-5")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 160.\n"
            elif(parameterName.find("Run-I-Jets-6")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 30.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 25.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 100.*GeV\n"
            elif(parameterName.find("Run-I-Jets-7")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 30.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 25.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 400.*GeV\n"
            elif(parameterName.find("Run-I-Jets-8")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 30.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 25.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 700.*GeV\n"
            elif(parameterName.find("Run-II-Jets-0")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 15.\n"
            elif(parameterName.find("Run-II-Jets-1")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 25.\n"
            elif(parameterName.find("Run-II-Jets-2")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 40.\n"
            elif(parameterName.find("Run-II-Jets-3")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 60.\n"
            elif(parameterName.find("Run-II-Jets-4")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 85.\n"
            elif(parameterName.find("Run-II-Jets-5")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 110.\n"
            elif(parameterName.find("Run-II-Jets-6")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 160.\n"
            elif(parameterName.find("Run-II-Jets-7")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 250.\n"
            elif(parameterName.find("Run-II-Jets-8")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 30.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 25.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 100.*GeV\n"
            elif(parameterName.find("Run-II-Jets-9")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 30.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 25.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 300.*GeV\n"
            elif(parameterName.find("900-Jets-1")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 10.\n"
            elif(parameterName.find("300-Jets-1")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 6.\n"
            elif(parameterName.find("630-Jets-1")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 20.\n"
            elif(parameterName.find("630-Jets-2")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 40.\n"
            elif(parameterName.find("630-Jets-3")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 75.\n"
            elif(parameterName.find("900-Jets-1")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 10.\n"
        elif(parameterName.find("Run-I-WZ")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p pbar e+ e-\ndo Factory:Process p pbar e+ nu\ndo Factory:Process p pbar e- nu\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("Run-I-W")>=0 or parameterName.find("Run-II-W")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p pbar e+ nu\ndo Factory:Process p pbar e- nu\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("Run-I-Z")>=0 or parameterName.find("Run-II-Z-e")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p pbar e+ e-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("Run-II-Z-LowMass-mu")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p pbar mu+ mu-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 25*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 70*GeV\n"
        elif(parameterName.find("Run-II-Z-HighMass-mu")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p pbar mu+ mu-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 150.*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 600*GeV\n"
        elif(parameterName.find("Run-II-Z-mu")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p pbar mu+ mu-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
# Star
elif(collider=="Star" ) :
    process = "set /Herwig/Decays/DecayHandler:LifeTimeOption 0\n"
    process+= "set /Herwig/Decays/DecayHandler:MaxLifeTime 10*mm\n"
    process+= "set /Herwig/Generators/EventGenerator:EventHandler:BeamB /Herwig/Particles/p+\n"
    process+= "set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 200.0\n"
    process+= "set /Herwig/Cuts/QCDCuts:X2Min 0.01\n"
    if(simulation=="") :
        if(parameterName.find("UE")>=0) :
            process += "insert SimpleQCD:MatrixElements[0] MEMinBias\n"
            process += "set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            process += "set /Herwig/Generators/EventGenerator:EventHandler:Cuts /Herwig/Cuts/MinBiasCuts\n"
            process += "create Herwig::MPIXSecReweighter /Herwig/Generators/MPIXSecReweighter\n"
            process += "insert /Herwig/Generators/EventGenerator:EventHandler:PostSubProcessHandlers 0 /Herwig/Generators/MPIXSecReweighter\n"
        else :
            process+="insert SimpleQCD:MatrixElements[0] MEQCD2to2\n"
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
        process += "insert SimpleQCD:MatrixElements[0] MEMinBias\n"
        process += "set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
        process += "set /Herwig/Generators/EventGenerator:EventHandler:Cuts /Herwig/Cuts/MinBiasCuts\n"
        process += "create Herwig::MPIXSecReweighter /Herwig/Generators/MPIXSecReweighter\n"
        process += "insert /Herwig/Generators/EventGenerator:EventHandler:PostSubProcessHandlers 0 /Herwig/Generators/MPIXSecReweighter\n"
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
            process+="insert SimpleQCD:MatrixElements[0] MEPP2HiggsVBF\n"
        elif(parameterName.find("VBF")>=0) :
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->tau-,tau+;\n"
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+="insert SimpleQCD:MatrixElements[0] MEPP2HiggsVBF\n"
        elif(parameterName.find("ggHJet")>=0) :
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->tau-,tau+;\n"
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+="insert SimpleQCD:MatrixElements[0] MEHiggsJet\n"
        elif(parameterName.find("8-ggH")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEHiggs\n"
            process+="insert SimpleQCD:MatrixElements[0] MEHiggsJet\n"
            process+="set MEHiggsJet:Process qqbar\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("ggH")>=0) :
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->tau-,tau+;\n"
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+="insert SimpleQCD:MatrixElements[0] MEHiggs\n"
            process+="insert SimpleQCD:MatrixElements[0] MEHiggsJet\n"
            process+="set MEHiggsJet:Process qqbar\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("PromptPhoton")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEGammaJet\n"
            if(parameterName.find("PromptPhoton-1")>=0) :
                process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            elif(parameterName.find("PromptPhoton-2")>=0) :
                process+="set /Herwig/Cuts/PhotonKtCut:MinKT 25.\n"
            elif(parameterName.find("PromptPhoton-3")>=0) :
                process+="set /Herwig/Cuts/PhotonKtCut:MinKT 80.\n"
            elif(parameterName.find("PromptPhoton-4")>=0) :
                process+="set /Herwig/Cuts/PhotonKtCut:MinKT 150.\n"
        elif(parameterName.find("DiPhoton-GammaGamma")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEGammaGamma\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaGamma","")
        elif(parameterName.find("DiPhoton-GammaJet")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEGammaJet\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaJet","")
        elif(parameterName.find("8-WH")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEPP2WH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("8-ZH")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEPP2ZH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("WH")>=0) :
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->b,bbar;\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes W+->nu_e,e+; W+->nu_mu,mu+;\n"
            process+="insert SimpleQCD:MatrixElements[0] MEPP2WH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("ZH")>=0) :
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->b,bbar;\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes Z0->e-,e+; Z0->mu-,mu+;\n"
            process+="insert SimpleQCD:MatrixElements[0] MEPP2ZH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("UE")>=0) :
            process += "insert SimpleQCD:MatrixElements[0] MEMinBias\n"
            process += "set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            process += "set /Herwig/Generators/EventGenerator:EventHandler:Cuts /Herwig/Cuts/MinBiasCuts\n"
            process += "create Herwig::MPIXSecReweighter /Herwig/Generators/MPIXSecReweighter\n"
            process += "insert /Herwig/Generators/EventGenerator:EventHandler:PostSubProcessHandlers 0 /Herwig/Generators/MPIXSecReweighter\n"
            if(parameterName.find("Long")>=0) :
                process += "set /Herwig/Decays/DecayHandler:MaxLifeTime 100*mm\n"
        elif(parameterName.find("7-Jets")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEQCD2to2\n"
            process+="set MEQCD2to2:MaximumFlavour 5\n"
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("7-Jets-0")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 5.\n"
            elif(parameterName.find("7-Jets-10")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 200.*GeV\n"
            elif(parameterName.find("7-Jets-11")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 600.*GeV\n"
            elif(parameterName.find("7-Jets-12")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 1000.*GeV\n"
            elif(parameterName.find("7-Jets-13")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 1600.*GeV\n"
            elif(parameterName.find("7-Jets-14")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 2200.*GeV\n"
            elif(parameterName.find("7-Jets-15")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 2800.*GeV\n"
            elif(parameterName.find("7-Jets-1")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 10.\n"
            elif(parameterName.find("7-Jets-2")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
            elif(parameterName.find("7-Jets-3")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 40.\n"
            elif(parameterName.find("7-Jets-4")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 70.\n"
            elif(parameterName.find("7-Jets-5")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 150.\n"
            elif(parameterName.find("7-Jets-6")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 200.\n"
            elif(parameterName.find("7-Jets-7")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 300.\n"
            elif(parameterName.find("7-Jets-8")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 500.\n"
            elif(parameterName.find("7-Jets-9")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 20.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 90.*GeV\n"
        elif(parameterName.find("7-Charm")>=0 or \
             parameterName.find("7-Bottom")>=0) :
            if(parameterName.find("7-Bottom")>=0) :
                process+="cp MEHeavyQuark MEBottom\n" 
                process+="set MEBottom:QuarkType Bottom\n"
                process+="insert SimpleQCD:MatrixElements[0] MEBottom\n"
            else : 
                process+="cp MEHeavyQuark MECharm\n" 
                process+="set MECharm:QuarkType Charm\n"
                process+="insert SimpleQCD:MatrixElements[0] MECharm\n"
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("7-Heavy-0")>=0) :
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
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 90.*GeV\n"
            elif(parameterName.find("-7")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 340.*GeV\n"
            elif(parameterName.find("-8")>=0) :
                process+="set /Herwig/Cuts/JetKtCut:MinKT 30.\n"
                process+="set /Herwig/Cuts/QCDCuts:MHatMin 500.*GeV\n"
        elif(parameterName.find("7-Top-L")>=0) :
            process+="set MEHeavyQuark:QuarkType Top\n"
            process+="insert SimpleQCD:MatrixElements[0] MEHeavyQuark\n"
            process+="do /Herwig/Particles/t:SelectDecayModes t->nu_e,e+,b; t->nu_mu,mu+,b;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("7-Top-SL")>=0) :
            process+="set MEHeavyQuark:QuarkType Top\n"
            process+="insert SimpleQCD:MatrixElements[0] MEHeavyQuark\n"
            process+="set /Herwig/Particles/t:Synchronized Not_synchronized\n"
            process+="set /Herwig/Particles/tbar:Synchronized Not_synchronized\n"
            process+="do /Herwig/Particles/t:SelectDecayModes t->nu_e,e+,b; t->nu_mu,mu+,b;\n"
            process+="do /Herwig/Particles/tbar:SelectDecayModes tbar->b,bbar,cbar; tbar->bbar,cbar,d; tbar->bbar,cbar,s; tbar->bbar,s,ubar; tbar->bbar,ubar,d;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("7-Top-All")>=0) :
            process+="set MEHeavyQuark:QuarkType Top\n"
            process+="insert SimpleQCD:MatrixElements[0] MEHeavyQuark\n"
        elif(parameterName.find("WZ")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process WZ\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+; /Herwig/Particles/W+/W+->nu_mu,mu+;\n"
            process+="do /Herwig/Particles/W-:SelectDecayModes /Herwig/Particles/W-/W-->nu_ebar,e-; /Herwig/Particles/W-/W-->nu_mubar,mu-;\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("WW-emu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process WW\n"
            process+="set /Herwig/Particles/W+:Synchronized 0\n"
            process+="set /Herwig/Particles/W-:Synchronized 0\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+;\n"
            process+="do /Herwig/Particles/W-:SelectDecayModes /Herwig/Particles/W-/W-->nu_mubar,mu-;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("WW-ll")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process WW\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+; /Herwig/Particles/W+/W+->nu_mu,mu+; /Herwig/Particles/W+/W+->nu_tau,tau+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("ZZ-ll")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process ZZ\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+; /Herwig/Particles/Z0/Z0->tau-,tau+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("ZZ-lv")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process ZZ\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+; /Herwig/Particles/Z0/Z0->tau-,tau+; /Herwig/Particles/Z0/Z0->nu_e,nu_ebar; /Herwig/Particles/Z0/Z0->nu_mu,nu_mubar; /Herwig/Particles/Z0/Z0->nu_tau,nu_taubar;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("W-Z-e")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
            process+="insert SimpleQCD:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Electron\n"
        elif(parameterName.find("W-Z-mu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
            process+="insert SimpleQCD:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Muon\n"
        elif(parameterName.find("W-e")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Electron\n"
        elif(parameterName.find("W-mu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Muon\n"
        elif(parameterName.find("Z-e")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Z-mu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("Z-LowMass-e")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
            process+="set /Herwig/Cuts/QCDCuts:MHatMin 20.*GeV\nset /Herwig/Cuts/MassCut:MinM 20.*GeV\nset /Herwig/Cuts/MassCut:MaxM 70.*GeV\n"
        elif(parameterName.find("Z-MedMass-e")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
            process+="set /Herwig/Cuts/QCDCuts:MHatMin 40.*GeV\nset /Herwig/Cuts/MassCut:MinM 40.*GeV\nset /Herwig/Cuts/MassCut:MaxM 130.*GeV\n"
        elif(parameterName.find("Z-LowMass-mu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
            process+="set /Herwig/Cuts/QCDCuts:MHatMin 10.*GeV\nset /Herwig/Cuts/MassCut:MinM 10.*GeV\nset /Herwig/Cuts/MassCut:MaxM 70.*GeV\n"
        elif(parameterName.find("W-Jet")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEWJet\nset MEWJet:WDecay Electron\n"
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
                process+="insert SimpleQCD:MatrixElements[0] MEZJet\nset MEZJet:ZDecay Electron\n"
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
                process+="insert SimpleQCD:MatrixElements[0] MEZJet\nset MEZJet:ZDecay Muon\n"
                process+="set /Herwig/Cuts/ZBosonKtCut:MinKT 35.0*GeV\n"
                parameterName=parameterName.replace("Z-Jet-0-mu","Z-Jet-mu")
        else :
            logging.error(" Process %s not supported for internal matrix elements" % name)
            sys.exit(1)
    elif(simulation=="Powheg") :
        if(parameterName.find("8-VBF")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2HiggsVBF\n"
        elif(parameterName.find("VBF")>=0) :
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->tau-,tau+;\n"
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2HiggsVBF\n"
        elif(parameterName.find("ggHJet")>=0) :
            logging.error(" Process %s not supported for POWHEG matrix elements" % name)
            sys.exit(1)
        elif(parameterName.find("8-ggH")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEHiggs\n"
        elif(parameterName.find("ggH")>=0) :
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->tau-,tau+;\n"
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEHiggs\n"
        elif(parameterName.find("8-WH")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2WH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("8-ZH")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2ZH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("WH")>=0) :
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->b,bbar;\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes W+->nu_e,e+; W+->nu_mu,mu+;\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2WH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("ZH")>=0) :
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->b,bbar;\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes Z0->e-,e+; Z0->mu-,mu+;\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2ZH\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 0.0*GeV\n"
        elif(parameterName.find("UE")>=0) :
            logging.error(" Process %s not supported for powheg matrix elements" % name)
            sys.exit(1)
        elif(parameterName.find("WZ")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process WZ\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+; /Herwig/Particles/W+/W+->nu_mu,mu+;\n"
            process+="do /Herwig/Particles/W-:SelectDecayModes /Herwig/Particles/W-/W-->nu_ebar,e-; /Herwig/Particles/W-/W-->nu_mubar,mu-;\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("WW-emu")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process WW\n"
            process+="set /Herwig/Particles/W+:Synchronized 0\n"
            process+="set /Herwig/Particles/W-:Synchronized 0\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+;\n"
            process+="do /Herwig/Particles/W-:SelectDecayModes /Herwig/Particles/W-/W-->nu_mubar,mu-;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("WW-ll")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process WW\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+; /Herwig/Particles/W+/W+->nu_mu,mu+; /Herwig/Particles/W+/W+->nu_tau,tau+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("ZZ-ll")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process ZZ\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+; /Herwig/Particles/Z0/Z0->tau-,tau+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("ZZ-lv")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process ZZ\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+; /Herwig/Particles/Z0/Z0->tau-,tau+; /Herwig/Particles/Z0/Z0->nu_e,nu_ebar; /Herwig/Particles/Z0/Z0->nu_mu,nu_mubar; /Herwig/Particles/Z0/Z0->nu_tau,nu_taubar;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("W-Z-e")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2W2ff\nset PowhegMEqq2W2ff:Process Electron\n"
        elif(parameterName.find("W-Z-mu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
            process+="insert SimpleQCD:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Muon\n"
        elif(parameterName.find("W-e")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2W2ff\nset PowhegMEqq2W2ff:Process Electron\n"
        elif(parameterName.find("W-mu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2W2ff\nset PowhegMEqq2W2ff:Process Muon\n"
        elif(parameterName.find("Z-e")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Z-mu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
        elif(parameterName.find("Z-LowMass-e")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
            process+="set /Herwig/Cuts/QCDCuts:MHatMin 20.*GeV\nset /Herwig/Cuts/MassCut:MinM 20.*GeV\nset /Herwig/Cuts/MassCut:MaxM 70.*GeV\n"
        elif(parameterName.find("Z-MedMass-e")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
            process+="set /Herwig/Cuts/QCDCuts:MHatMin 40.*GeV\nset /Herwig/Cuts/MassCut:MinM 40.*GeV\nset /Herwig/Cuts/MassCut:MaxM 130.*GeV\n"
        elif(parameterName.find("Z-LowMass-mu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
            process+="set /Herwig/Cuts/QCDCuts:MHatMin 10.*GeV\nset /Herwig/Cuts/MassCut:MinM 10.*GeV\nset /Herwig/Cuts/MassCut:MaxM 70.*GeV\n"
        elif(parameterName.find("DiPhoton-GammaGamma")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEGammaGammaPowheg\n"
            process+="set MEGammaGammaPowheg:Process GammaGamma\n"
            process+="insert SimpleQCD:MatrixElements[0] MEGammaGamma\n"
            process+="set MEGammaGamma:Process gg\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaGamma","")
        elif(parameterName.find("DiPhoton-GammaJet")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEGammaGammaPowheg\n"
            process+="set MEGammaGammaPowheg:Process VJet\n"
            process+="set /Herwig/Cuts/PhotonKtCut:MinKT 5.\n"
            process+="set /Herwig/Cuts/JetKtCut:MinKT 5.\n"
            parameterName=parameterName.replace("-GammaJet","")
        else :
            logging.error(" Process %s not supported for internal POWHEG matrix elements" % name)
            sys.exit(1)
            
    elif(simulation=="Matchbox" ) :
        if(parameterName.find("8-VBF")>=0) :
            parameters["nlo"] = "read Matchbox/VBFNLO.in\n"
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 3\n"
            process+="insert Factory:DiagramGenerator:RestrictLines 0 /Herwig/Particles/Z0\n"
            process+="insert Factory:DiagramGenerator:RestrictLines 0 /Herwig/Particles/W+\n"
            process+="insert Factory:DiagramGenerator:RestrictLines 0 /Herwig/Particles/W-\n"
            process+="insert Factory:DiagramGenerator:RestrictLines 0 /Herwig/Particles/gamma\n"
            process+="do Factory:DiagramGenerator:TimeLikeRange 0 0\n"
            process+="do Factory:Process p p h0 j j\n"
            process+="set /Herwig/Particles/h0:HardProcessWidth 0.\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
        elif(parameterName.find("VBF")>=0) :
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->tau-,tau+;\n"
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            parameters["nlo"] = "read Matchbox/VBFNLO.in\n"
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 3\n"
            process+="insert Factory:DiagramGenerator:RestrictLines 0 /Herwig/Particles/Z0\n"
            process+="insert Factory:DiagramGenerator:RestrictLines 0 /Herwig/Particles/W+\n"
            process+="insert Factory:DiagramGenerator:RestrictLines 0 /Herwig/Particles/W-\n"
            process+="insert Factory:DiagramGenerator:RestrictLines 0 /Herwig/Particles/gamma\n"
            process+="do Factory:DiagramGenerator:TimeLikeRange 0 0\n"
            process+="do Factory:Process p p h0 j j\n"
            process+="set /Herwig/Particles/h0:HardProcessWidth 0.\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
        elif(parameterName.find("ggHJet")>=0) :
            parameters["nlo"] = "read Matchbox/MadGraph-GoSam.in\nread Matchbox/HiggsEffective.in\n"
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->tau-,tau+;\n"
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+="set Factory:OrderInAlphaS 3\nset Factory:OrderInAlphaEW 1\n"
            process+="set /Herwig/Particles/h0:HardProcessWidth 0.\n"
            process+="do Factory:Process p p h0 j\n"
            process+="set /Herwig/Cuts/Cuts:JetFinder /Herwig/Cuts/JetFinder\n"
            process+="insert /Herwig/Cuts/Cuts:MultiCuts 0 /Herwig/Cuts/JetCuts\n"
            process+="insert /Herwig/Cuts/JetCuts:JetRegions 0 /Herwig/Cuts/FirstJet\n"
            process+="set /Herwig/Cuts/FirstJet:PtMin 20.\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
        elif(parameterName.find("8-ggH")>=0) :
            parameters["nlo"] = "read Matchbox/MadGraph-GoSam.in\nread Matchbox/HiggsEffective.in\n"
            process+="set Factory:OrderInAlphaS 2\nset Factory:OrderInAlphaEW 1\n"
            process+="set /Herwig/Particles/h0:HardProcessWidth 0.\n"
            process+="do Factory:Process p p h0\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
        elif(parameterName.find("ggH")>=0) :
            parameters["nlo"] = "read Matchbox/MadGraph-GoSam.in\nread Matchbox/HiggsEffective.in\n"
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->tau-,tau+;\n"
            process+="set /Herwig/Particles/tau-:Stable Stable\n"
            process+="set Factory:OrderInAlphaS 2\nset Factory:OrderInAlphaEW 1\n"
            process+="set /Herwig/Particles/h0:HardProcessWidth 0.\n"
            process+="do Factory:Process p p h0\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
        elif(parameterName.find("8-WH")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\n"
            process+="set /Herwig/Particles/h0:HardProcessWidth 0.\n"
            process+="do Factory:Process p p W+ h0\n"
            process+="do Factory:Process p p W- h0\n"
            process+="set /Herwig/Particles/W+:HardProcessWidth 0.\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
        elif(parameterName.find("8-ZH")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\n"
            process+="set /Herwig/Particles/h0:HardProcessWidth 0.\n"
            process+="set /Herwig/Particles/Z0:HardProcessWidth 0.\n"
            process+="do Factory:Process p p Z0 h0\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 125.7\n"
        elif(parameterName.find("WH")>=0) :
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->b,bbar;\n"
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 3\n"
            process+="set /Herwig/Particles/h0:HardProcessWidth 0.\n"
            process+="do Factory:Process p p e+ nu h0\n"
            process+="do Factory:Process p p e- nu h0\n"
            process+="do Factory:Process p p mu+ nu h0\n"
            process+="do Factory:Process p p mu- nu h0\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("ZH")>=0) :
            process+="do /Herwig/Particles/h0:SelectDecayModes h0->b,bbar;\n"
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 3\n"
            process+="set /Herwig/Particles/h0:HardProcessWidth 0.\n"
            process+="do Factory:Process p p e+ e- h0\n"
            process+="do Factory:Process p p mu+ mu- h0\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("UE")>=0) :
            logging.error(" Process %s not supported for Matchbox matrix elements" % name)
            sys.exit(1)
        elif(parameterName.find("7-Jets")>=0) :
            process+="set Factory:OrderInAlphaS 2\nset Factory:OrderInAlphaEW 0\n"
            process+="do Factory:Process p p j j\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/MaxJetPtScale\n"
            process+="set  /Herwig/Cuts/Cuts:JetFinder  /Herwig/Cuts/JetFinder\n"
            process+="insert  /Herwig/Cuts/Cuts:MultiCuts 0  /Herwig/Cuts/JetCuts\n"
            process+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/FirstJet\n"
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("7-Jets-0")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 5.\n"
            elif(parameterName.find("7-Jets-10")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 20.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 15.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 200.*GeV\n"
            elif(parameterName.find("7-Jets-11")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 20.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 15.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 600.*GeV\n"
            elif(parameterName.find("7-Jets-12")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 20.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 15.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 1000.*GeV\n"
            elif(parameterName.find("7-Jets-13")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 20.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 15.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 1600.*GeV\n"
            elif(parameterName.find("7-Jets-14")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 20.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 15.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 2200.*GeV\n"
            elif(parameterName.find("7-Jets-15")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 20.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 15.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 2800.*GeV\n"
            elif(parameterName.find("7-Jets-1")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 10.\n"
            elif(parameterName.find("7-Jets-2")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 20.\n"
            elif(parameterName.find("7-Jets-3")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 40.\n"
            elif(parameterName.find("7-Jets-4")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 70.\n"
            elif(parameterName.find("7-Jets-5")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 150.\n"
            elif(parameterName.find("7-Jets-6")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 200.\n"
            elif(parameterName.find("7-Jets-7")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 300.\n"
            elif(parameterName.find("7-Jets-8")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 500.\n"
            elif(parameterName.find("7-Jets-9")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 20.*GeV\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 15.*GeV\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 90.*GeV\n"
        elif(parameterName.find("7-Charm")>=0 or \
             parameterName.find("7-Bottom")>=0) :
            parameters["bscheme"]=fourFlavour
            process+="set /Herwig/Particles/b:HardProcessMass 4.2*GeV\n"
            process+="set /Herwig/Particles/bbar:HardProcessMass 4.2*GeV\n"
            process+="set Factory:OrderInAlphaS 2\nset Factory:OrderInAlphaEW 0\n"
            if(parameterName.find("7-Bottom")>=0) :
                process+="do Factory:Process p p b bbar\n"
            else:
                process+="do Factory:Process p p c cbar\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/MaxJetPtScale\n"
            process+="set  /Herwig/Cuts/Cuts:JetFinder  /Herwig/Cuts/JetFinder\n"
            process+="insert  /Herwig/Cuts/Cuts:MultiCuts 0  /Herwig/Cuts/JetCuts\n"
            process+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/FirstJet\n"
            process+="set /Herwig/UnderlyingEvent/MPIHandler:IdenticalToUE 0\n"
            if(parameterName.find("-0")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 0.\n"
            elif(parameterName.find("-1")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 5.\n"
            elif(parameterName.find("-2")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 20.\n"
            elif(parameterName.find("-3")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 50.\n"
            elif(parameterName.find("-4")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 80.\n"
            elif(parameterName.find("-5")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 110.\n"
            elif(parameterName.find("-6")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 30.\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 25.\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 90.*GeV\n"
            elif(parameterName.find("-7")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 30.\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 25.\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 340.*GeV\n"
            elif(parameterName.find("-8")>=0) :
                process+="insert /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 30.\n"
                process+="set /Herwig/Cuts/SecondJet:PtMin 25.\n"
                process+="create ThePEG::JetPairRegion /Herwig/Cuts/JetPairMass JetCuts.so\n"
                process+="set /Herwig/Cuts/JetPairMass:FirstRegion /Herwig/Cuts/FirstJet\n"
                process+="set /Herwig/Cuts/JetPairMass:SecondRegion /Herwig/Cuts/SecondJet\n"
                process+="insert /Herwig/Cuts/JetCuts:JetPairRegions 0  /Herwig/Cuts/JetPairMass\n"
                process+="set /Herwig/Cuts/JetPairMass:MassMin 500.*GeV\n"
        elif(parameterName.find("7-Top-L")>=0) :
            process+="set /Herwig/Particles/t:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/tbar:HardProcessWidth 0.*GeV\n"
            process+="set Factory:OrderInAlphaS 2\nset Factory:OrderInAlphaEW 0\n"
            process+="do Factory:Process p p t tbar\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/TopPairMTScale\n"
            process+="do /Herwig/Particles/t:SelectDecayModes t->nu_e,e+,b; t->nu_mu,mu+,b;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("7-Top-SL")>=0) :
            process+="set /Herwig/Particles/t:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/tbar:HardProcessWidth 0.*GeV\n"
            process+="set Factory:OrderInAlphaS 2\nset Factory:OrderInAlphaEW 0\n"
            process+="do Factory:Process p p t tbar\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/TopPairMTScale\n"
            process+="set /Herwig/Particles/t:Synchronized Not_synchronized\n"
            process+="set /Herwig/Particles/tbar:Synchronized Not_synchronized\n"
            process+="do /Herwig/Particles/t:SelectDecayModes t->nu_e,e+,b; t->nu_mu,mu+,b;\n"
            process+="do /Herwig/Particles/tbar:SelectDecayModes tbar->b,bbar,cbar; tbar->bbar,cbar,d; tbar->bbar,cbar,s; tbar->bbar,s,ubar; tbar->bbar,ubar,d;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("7-Top-All")>=0) :
            process+="set /Herwig/Particles/t:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/tbar:HardProcessWidth 0.*GeV\n"
            process+="set Factory:OrderInAlphaS 2\nset Factory:OrderInAlphaEW 0\n"
            process+="do Factory:Process p p t tbar\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/TopPairMTScale\n"
        elif(parameterName.find("WZ")>=0) :
            process+="set /Herwig/Particles/W+:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/W-:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/Z0:HardProcessWidth 0.*GeV\n"
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p W+ Z0\ndo Factory:Process p p W- Z0\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 171.6*GeV\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+; /Herwig/Particles/W+/W+->nu_mu,mu+;\n"
            process+="do /Herwig/Particles/W-:SelectDecayModes /Herwig/Particles/W-/W-->nu_ebar,e-; /Herwig/Particles/W-/W-->nu_mubar,mu-;\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("WW-emu")>=0) :
            process+="set /Herwig/Particles/W+:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/W-:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/Z0:HardProcessWidth 0.*GeV\n"
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p W+ W-\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 160.8*GeV\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/Particles/W+:Synchronized 0\n"
            process+="set /Herwig/Particles/W-:Synchronized 0\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+;\n"
            process+="do /Herwig/Particles/W-:SelectDecayModes /Herwig/Particles/W-/W-->nu_mubar,mu-;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("WW-ll")>=0) :
            process+="set /Herwig/Particles/W+:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/W-:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/Z0:HardProcessWidth 0.*GeV\n"
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p W+ W-\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 160.8*GeV\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+; /Herwig/Particles/W+/W+->nu_mu,mu+; /Herwig/Particles/W+/W+->nu_tau,tau+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("ZZ-ll")>=0) :
            process+="set /Herwig/Particles/W+:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/W-:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/Z0:HardProcessWidth 0.*GeV\n"
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p Z0 Z0\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 182.2*GeV\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+; /Herwig/Particles/Z0/Z0->tau-,tau+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("ZZ-lv")>=0) :
            process+="set /Herwig/Particles/W+:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/W-:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/Z0:HardProcessWidth 0.*GeV\n"
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p Z0 Z0\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 182.2*GeV\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+; /Herwig/Particles/Z0/Z0->tau-,tau+; /Herwig/Particles/Z0/Z0->nu_e,nu_ebar; /Herwig/Particles/Z0/Z0->nu_mu,nu_mubar; /Herwig/Particles/Z0/Z0->nu_tau,nu_taubar;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generators/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("W-Z-e")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\n"
            process+="do Factory:Process p p e+ e-\ndo Factory:Process p p e+ nu\ndo Factory:Process p p e- nu\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("W-Z-mu")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\n"
            process+="do Factory:Process p p mu+ mu-\ndo Factory:Process p p mu+ nu\ndo Factory:Process p p mu- nu\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("W-e")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ nu\ndo Factory:Process p p e- nu\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("W-mu")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p mu+ nu\ndo Factory:Process p p mu- nu\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("Z-e")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ e-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("Z-mu")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p mu+ mu-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("Z-jj")>=0) :
            process+="set Factory:OrderInAlphaS 2\nset Factory:OrderInAlphaEW 2\n"
            process+="do Factory:Process p p e+ e- j j\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set  /Herwig/Cuts/Cuts:JetFinder  /Herwig/Cuts/JetFinder\n"
            process+="insert  /Herwig/Cuts/Cuts:MultiCuts 0  /Herwig/Cuts/JetCuts\n"
            process+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/FirstJet\n"
            process+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
            process+="set /Herwig/Cuts/FirstJet:PtMin 40.*GeV\n"
            process+="set /Herwig/Cuts/SecondJet:PtMin 30.*GeV\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("Z-LowMass-e")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ e-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 20*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 70*GeV\n"
        elif(parameterName.find("Z-MedMass-e")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ e-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 40*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 130*GeV\n"
        elif(parameterName.find("Z-LowMass-mu")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p mu+ mu-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 10*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 70*GeV\n"
        elif(parameterName.find("W-Jet")>=0) :
            process+="set Factory:OrderInAlphaS 1\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ nu j\ndo Factory:Process p p e- nu j\n\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/HTScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
            process+="set /Herwig/Cuts/Cuts:JetFinder /Herwig/Cuts/JetFinder\n"
            process+="insert /Herwig/Cuts/Cuts:MultiCuts 0 /Herwig/Cuts/JetCuts\n"
            process+="insert /Herwig/Cuts/JetCuts:JetRegions 0 /Herwig/Cuts/FirstJet\n"
            if(parameterName.find("W-Jet-1-e")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 100.*GeV\n"
                parameterName=parameterName.replace("W-Jet-1-e","W-Jet-e")
            elif(parameterName.find("W-Jet-2-e")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 190.0*GeV\n"
                parameterName=parameterName.replace("W-Jet-2-e","W-Jet-e")
            elif(parameterName.find("W-Jet-3-e")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 270.0*GeV\n"
                parameterName=parameterName.replace("W-Jet-3-e","W-Jet-e")
        elif(parameterName.find("Z-Jet")>=0) :
            process+="set Factory:OrderInAlphaS 1\nset Factory:OrderInAlphaEW 2\n"
            if(parameterName.find("-e")>=0) :
                process+="do Factory:Process p p e+ e- j\n"
                if(parameterName.find("Z-Jet-0-e")>=0) :
                    process+="set /Herwig/Cuts/FirstJet:PtMin 35.*GeV\n"
                    parameterName=parameterName.replace("Z-Jet-0-e","Z-Jet-e")
                elif(parameterName.find("Z-Jet-1-e")>=0) :
                    process+="set /Herwig/Cuts/FirstJet:PtMin 100.*GeV\n"
                    parameterName=parameterName.replace("Z-Jet-1-e","Z-Jet-e")
                elif(parameterName.find("Z-Jet-2-e")>=0) :
                    process+="set /Herwig/Cuts/FirstJet:PtMin 190.0*GeV\n"
                    parameterName=parameterName.replace("Z-Jet-2-e","Z-Jet-e")
                elif(parameterName.find("Z-Jet-3-e")>=0) :
                    process+="set /Herwig/Cuts/FirstJet:PtMin 270.0*GeV\n"
                    parameterName=parameterName.replace("Z-Jet-3-e","Z-Jet-e")
            else :
                process+="do Factory:Process p p mu+ mu- j\n"
                process+="set /Herwig/Cuts/FirstJet:PtMin 35.*GeV\n"
                parameterName=parameterName.replace("Z-Jet-0-mu","Z-Jet-mu")
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/HTScale\n"
            process+="set /Herwig/Cuts/Cuts:JetFinder /Herwig/Cuts/JetFinder\n"
            process+="insert /Herwig/Cuts/Cuts:MultiCuts 0 /Herwig/Cuts/JetCuts\n"
            process+="insert /Herwig/Cuts/JetCuts:JetRegions 0 /Herwig/Cuts/FirstJet\n"
        elif(parameterName.find("Z-bb")>=0) :
            parameters["bscheme"]=fourFlavour
            process+="set /Herwig/Particles/b:HardProcessMass 4.2*GeV\nset /Herwig/Particles/bbar:HardProcessMass 4.2*GeV\n"
            process+="set Factory:OrderInAlphaS 2\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ e- b bbar\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 91.2*GeV\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 66*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 116*GeV\n"
            process+="set /Herwig/Cuts/Cuts:JetFinder /Herwig/Cuts/JetFinder\n"
            process+="insert  /Herwig/Cuts/Cuts:MultiCuts 0  /Herwig/Cuts/JetCuts\n"
            process+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/FirstJet\n"
            process+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
            process+="set  /Herwig/Cuts/FirstJet:PtMin 18.*GeV\n"
            process+="set  /Herwig/Cuts/SecondJet:PtMin 15.*GeV\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("Z-b")>=0) :
            process+="do Factory:StartParticleGroup bjet\n"
            process+="insert Factory:ParticleGroup 0 /Herwig/Particles/b\n"
            process+="insert Factory:ParticleGroup 0 /Herwig/Particles/bbar\n"
            process+="do Factory:EndParticleGroup\n"
            process+="set Factory:OrderInAlphaS 1\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ e- bjet\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 91.2*GeV\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
            process+="set /Herwig/Cuts/Cuts:JetFinder /Herwig/Cuts/JetFinder\n"
            process+="insert  /Herwig/Cuts/Cuts:MultiCuts 0  /Herwig/Cuts/JetCuts\n"
            process+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/FirstJet\n"
            process+="set  /Herwig/Cuts/FirstJet:PtMin 15.*GeV\n"
        elif(parameterName.find("W-b")>=0) :
            parameters["bscheme"]=fourFlavour
            process += "set /Herwig/Particles/b:HardProcessMass 4.2*GeV\nset /Herwig/Particles/bbar:HardProcessMass 4.2*GeV\n"
            process += "set Factory:OrderInAlphaS 2\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ nu b bbar\ndo Factory:Process p p e- nu b bbar\n"
            process += "set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 80.4*GeV\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/Cuts/Cuts:JetFinder /Herwig/Cuts/JetFinder\n"
            process+="insert  /Herwig/Cuts/Cuts:MultiCuts 0  /Herwig/Cuts/JetCuts\n"
            process+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/FirstJet\n"
            process+="set  /Herwig/Cuts/FirstJet:PtMin 30.*GeV\n"
            process+="set /Herwig/Cuts/LeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/LeptonPairMassCut:MaxMass 120*GeV\n"
        else :
            logging.error(" Process %s not supported for Matchbox matrix elements" % name)
            sys.exit(1)
# Star
elif(collider=="LHC-GammaGamma" ) :
    if(parameterName.find("-7-")>=0) :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 7000.0\n"
    elif(parameterName.find("-8-")>=0) :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 8000.0\n"
    else :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 7000.0\n"
    if(simulation=="") :
        if(parameterName.find("7")>=0) :
            process += "insert SimpleQCD:MatrixElements 0 /Herwig/MatrixElements/MEgg2ff\n"
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

# write the file
if(simulation=="Matchbox" ) :
    with open(os.path.join("Rivet",name+".in") ,'w') as f:
        f.write( template.substitute(parameters))
else :
    with open(os.path.join("Rivet",name+".in") ,'w') as f:
        f.write( template.substitute(parameters))






