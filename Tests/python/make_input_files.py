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

collider=""
# select the template to load
# collider
if(name.find("BFactory")==0) :
    collider="BFactory"
elif(name.find("LEP")==0) :
    collider="LEP"
elif(name.find("DIS")==0) :
    collider="DIS"
elif(name.find("TVT")==0) :
    collider="TVT"
elif(name.find("LHC")==0) :
    collider="LHC"
simulation=""
if(name.find("Matchbox")>0) :
    simulation="Matchbox"
elif(name.find("Dipole")>0) :
    simulation="Dipole"
elif(name.find("Powheg")>0) :
    simulation="Powheg"

if(collider=="") :
    logging.error("Can\'t find collider")
    sys.exit(1)

# find the template
if(simulation=="") :
    if(collider.find("TVT")>=0 or collider.find("LHC")>=0) :
        templateName="Hadron.in"
    elif(collider.find("BFactory")<0) :
        templateName= "%s.in" % (collider)
    else :
        templateName= "LEP.in"
else :
    if(collider.find("TVT")>=0 or collider.find("LHC")>=0) :
        templateName= "Hadron-%s.in" % (simulation) 
    elif(collider.find("BFactory")<0) :
        templateName= "%s-%s.in" % (collider,simulation) 
    else :
        templateName= "LEP-%s.in" % (simulation) 
with open(os.path.join("Rivet/Templates",templateName), 'r') as f:
    templateText = f.read()
template = Template( templateText )

# work out the name of the parameter file
if(simulation=="") : 
    istart = 1
else :
    istart = 2
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
    elif(simulation=="Matchbox" or simulation == "Dipole" ) :
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
    elif(simulation=="Matchbox" or simulation == "Dipole" ) :
        if(parameterName.find("e-")>=0) :
            process="do Factory:Process e- p -> e- j"
        else :
            process="do Factory:Process e+ p -> e+ j"
# LEP
elif(collider=="LEP") :
    if(simulation=="") :
        process=""
        if(parameterName=="-10") :
            process="set /Herwig/MatrixElements/MEee2gZ2qq:MaximumFlavour 4"
    elif(simulation=="Powheg") :
        process=""
        if(parameterName=="-10") :
            process="set /Herwig/MatrixElements/PowhegMEee2gZ2qq:MaximumFlavour 4"
    elif(simulation=="Matchbox" or simulation == "Dipole" ) :
        if(parameterName=="-10") :
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

    if(simulation=="") :
        if(parameterName.find("Run-I-WZ")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Electron\ninsert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Run-I-W")>=0 or parameterName.find("Run-II-W")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Electron\n"
        elif(parameterName.find("Run-I-Z")>=0 or parameterName.find("Run-II-Z-e")>=0) :
            process +="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Run-II-Z-mu")>=0) :
            process +="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Muon\n"
    elif(simulation=="Powheg") :
        if(parameterName.find("Run-I-WZ")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2W2ff\nset PowhegMEqq2W2ff:Process Electron\ninsert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Run-I-W")>=0 or parameterName.find("Run-II-W")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2W2ff\nset PowhegMEqq2W2ff:Process Electron\n"
        elif(parameterName.find("Run-I-Z")>=0 or parameterName.find("Run-II-Z-e")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
        elif(parameterName.find("Run-II-Z-mu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Muon\n"
    elif(simulation=="Matchbox" or simulation == "Dipole" ) :
        if(parameterName.find("Run-I-WZ")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p pbar e+ e-\ndo Factory:Process p pbar e+ nu\ndo Factory:Process p pbar e- nu\n"
            process+="set /Herwig/Cuts/ChargedLeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/ChargedLeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("Run-I-W")>=0 or parameterName.find("Run-II-W")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p pbar e+ nu\ndo Factory:Process p pbar e- nu\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
        elif(parameterName.find("Run-I-Z")>=0 or parameterName.find("Run-II-Z-e")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p pbar e+ e-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/ChargedLeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/ChargedLeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("Run-II-Z-mu")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p pbar mu+ mu-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/ChargedLeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/ChargedLeptonPairMassCut:MaxMass 120*GeV\n"


# LHC
elif(collider=="LHC") :
    if(parameterName.find("-7-")>=0) :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 7000.0\n"
    elif(parameterName.find("-8-")>=0) :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 8000.0\n"
    else :
        process="set /Herwig/Generators/EventGenerator:EventHandler:LuminosityFunction:Energy 7000.0\n"

    if(simulation=="") :
        if(parameterName.find("WZ")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process WZ\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+; /Herwig/Particles/W+/W+->nu_mu,mu+;\n"
            process+="do /Herwig/Particles/W-:SelectDecayModes /Herwig/Particles/W-/W-->nu_ebar,e-; /Herwig/Particles/W-/W-->nu_mubar,mu-;\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/EventHandlers/LHCHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("WW-emu")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process WW\n"
            process+="set /Herwig/Particles/W+:Synchronized 0\n"
            process+="set /Herwig/Particles/W-:Synchronized 0\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+;\n"
            process+="do /Herwig/Particles/W-:SelectDecayModes /Herwig/Particles/W-/W-->nu_mubar,mu-;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/EventHandlers/LHCHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("WW-ll")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process WW\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+; /Herwig/Particles/W+/W+->nu_mu,mu+; /Herwig/Particles/W+/W+->nu_tau,tau+\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/EventHandlers/LHCHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("ZZ-ll")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process ZZ\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+; /Herwig/Particles/Z0/Z0->tau-,tau+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/EventHandlers/LHCHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("ZZ-lv")>=0) :
            process+="insert SimpleQCD:MatrixElements[0] MEPP2VV\nset MEPP2VV:Process ZZ\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+; /Herwig/Particles/Z0/Z0->tau-,tau+; /Herwig/Particles/Z0/Z0->nu_e,nu_ebar; /Herwig/Particles/Z0/Z0->nu_mu,nu_mubar; /Herwig/Particles/Z0/Z0->nu_tau,nu_taubar;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/EventHandlers/LHCHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("LHC-W-Z-e")>0) :
            process+="insert SimpleQCD:MatrixElements[0] MEqq2gZ2ff\nset MEqq2gZ2ff:Process Electron\n"
            process+="insert SimpleQCD:MatrixElements[0] MEqq2W2ff\nset MEqq2W2ff:Process Electron\n"
        elif(parameterName.find("LHC-W-Z-mu")>0) :
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
            process+="set /Herwig/Cuts/QCDCuts:MHatMin 20.*GeV\nset /Herwig/Cuts/MassCut:MinM 10.*GeV\nset /Herwig/Cuts/MassCut:MaxM 70.*GeV\n"
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
            process+="insert SimpleQCD:MatrixElements[0] MEZJet\nset MEZJet:ZDecay Electron\n"
            if(parameterName.find("Z-Jet-1-e")>=0) :
                process+="set /Herwig/Cuts/ZBosonKtCut:MinKT 100.0*GeV\n"
                parameterName=parameterName.replace("Z-Jet-1-e","Z-Jet-e")
            elif(parameterName.find("Z-Jet-2-e")>=0) :
                process+="set /Herwig/Cuts/ZBosonKtCut:MinKT 190.0*GeV\n"
                parameterName=parameterName.replace("Z-Jet-2-e","Z-Jet-e")
            elif(parameterName.find("Z-Jet-3-e")>=0) :
                process+="set /Herwig/Cuts/ZBosonKtCut:MinKT 270.0*GeV\n"
                parameterName=parameterName.replace("Z-Jet-3-e","Z-Jet-e")
        else :
            logging.error(" Process %s not supported for internal matrix elements" % name)
            sys.exit(1)
    elif(simulation=="Powheg") :
        if(parameterName.find("WZ")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generator/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process WZ\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+; /Herwig/Particles/W+/W+->nu_mu,mu+;\n"
            process+="do /Herwig/Particles/W-:SelectDecayModes /Herwig/Particles/W-/W-->nu_ebar,e-; /Herwig/Particles/W-/W-->nu_mubar,mu-;\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generator/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("WW-emu")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generator/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process WW\n"
            process+="set /Herwig/Particles/W+:Synchronized 0\n"
            process+="set /Herwig/Particles/W-:Synchronized 0\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+;\n"
            process+="do /Herwig/Particles/W-:SelectDecayModes /Herwig/Particles/W-/W-->nu_mubar,mu-;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/EventHandlers/LHCHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("WW-ll")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generator/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process WW\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+; /Herwig/Particles/W+/W+->nu_mu,mu+; /Herwig/Particles/W+/W+->nu_tau,tau+\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/EventHandlers/LHCHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("ZZ-ll")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generator/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process ZZ\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+; /Herwig/Particles/Z0/Z0->tau-,tau+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/EventHandlers/LHCHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("ZZ-lv")>=0) :
            process+="create Herwig::HwDecayHandler /Herwig/NewPhysics/DecayHandler\n"
            process+="set /Herwig/NewPhysics/DecayHandler:NewStep No\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 0 /Herwig/Particles/tau-\n"
            process+="insert /Herwig/NewPhysics/DecayHandler:Excluded 1 /Herwig/Particles/tau+\n"
            process+="insert /Herwig/Generator/EventGenerator:EventHandler:PreCascadeHandlers 0 /Herwig/NewPhysics/DecayHandler\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEPP2VV\nset PowhegMEPP2VV:Process ZZ\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+; /Herwig/Particles/Z0/Z0->tau-,tau+; /Herwig/Particles/Z0/Z0->nu_e,nu_ebar; /Herwig/Particles/Z0/Z0->nu_mu,nu_mubar; /Herwig/Particles/Z0/Z0->nu_tau,nu_taubar;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/EventHandlers/LHCHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("LHC-W-Z-e")>0) :
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2gZ2ff\nset PowhegMEqq2gZ2ff:Process Electron\n"
            process+="insert SimpleQCD:MatrixElements[0] PowhegMEqq2W2ff\nset PowhegMEqq2W2ff:Process Electron\n"
        elif(parameterName.find("LHC-W-Z-mu")>0) :
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
            process+="set /Herwig/Cuts/QCDCuts:MHatMin 20.*GeV\nset /Herwig/Cuts/MassCut:MinM 10.*GeV\nset /Herwig/Cuts/MassCut:MaxM 70.*GeV\n"
        else :
            logging.error(" Process %s not supported for internal POWHEG matrix elements" % name)
            sys.exit(1)
            
    elif(simulation=="Matchbox" or simulation == "Dipole" ) :
        if(parameterName.find("WZ")>=0) :
            process+="set /Herwig/Particles/W+:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/W-:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/Z0:HardProcessWidth 0.*GeV\n"
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p W+ Z0\ndo Factory:Process p p W- Z0\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 171.6*GeV\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n\n"
            process+="set /Herwig/Cuts/MassCut:MinM 66*GeV\nset /Herwig/Cuts/MassCut:MaxM 116*GeV\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+; /Herwig/Particles/W+/W+->nu_mu,mu+;\n"
            process+="do /Herwig/Particles/W-:SelectDecayModes /Herwig/Particles/W-/W-->nu_ebar,e-; /Herwig/Particles/W-/W-->nu_mubar,mu-;\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/Generator/EventGenerator:EventHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("WW-emu")>=0) :
            process+="set /Herwig/Particles/W+:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/W-:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/Z0:HardProcessWidth 0.*GeV\n"
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process pnob pnob W+ W-\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 160.8*GeV\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/Particles/W+:Synchronized 0\n"
            process+="set /Herwig/Particles/W-:Synchronized 0\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+;\n"
            process+="do /Herwig/Particles/W-:SelectDecayModes /Herwig/Particles/W-/W-->nu_mubar,mu-;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/EventHandlers/LHCHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("WW-ll")>=0) :
            process+="set /Herwig/Particles/W+:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/W-:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/Z0:HardProcessWidth 0.*GeV\n"
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process pnob pnob W+ W-\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 160.8*GeV\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="do /Herwig/Particles/W+:SelectDecayModes /Herwig/Particles/W+/W+->nu_e,e+; /Herwig/Particles/W+/W+->nu_mu,mu+; /Herwig/Particles/W+/W+->nu_tau,tau+\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/EventHandlers/LHCHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("ZZ-ll")>=0) :
            process+="set /Herwig/Particles/W+:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/W-:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/Z0:HardProcessWidth 0.*GeV\n"
            process+="set PPFactory:OrderInAlphaS 0\nset PPFactory:OrderInAlphaEW 2\ndo PPFactory:Process p p Z0 Z0\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 182.2*GeV\nset PPFactory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+; /Herwig/Particles/Z0/Z0->tau-,tau+;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/EventHandlers/LHCHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("ZZ-lv")>=0) :
            process+="set /Herwig/Particles/W+:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/W-:HardProcessWidth 0.*GeV\n"
            process+="set /Herwig/Particles/Z0:HardProcessWidth 0.*GeV\n"
            process+="set PPFactory:OrderInAlphaS 0\nset PPFactory:OrderInAlphaEW 2\ndo PPFactory:Process p p Z0 Z0\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 182.2*GeV\nset PPFactory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="do /Herwig/Particles/Z0:SelectDecayModes /Herwig/Particles/Z0/Z0->e-,e+; /Herwig/Particles/Z0/Z0->mu-,mu+; /Herwig/Particles/Z0/Z0->tau-,tau+; /Herwig/Particles/Z0/Z0->nu_e,nu_ebar; /Herwig/Particles/Z0/Z0->nu_mu,nu_mubar; /Herwig/Particles/Z0/Z0->nu_tau,nu_taubar;\n"
            process+="create Herwig::BranchingRatioReweighter /Herwig/Generators/BRReweighter\n"
            process+="insert /Herwig/EventHandlers/LHCHandler:PostHadronizationHandlers 0 /Herwig/Generators/BRReweighter\n"
        elif(parameterName.find("LHC-W-Z-e")>0) :
            process+="set PPFactory:OrderInAlphaS 0\nset PPFactory:OrderInAlphaEW 2\n"
            process+="do PPFactory:Process p p e+ e-\ndo PPFactory:Process p p e+ nu\ndo PPFactory:Process p p e- nu\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/ChargedLeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/ChargedLeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("LHC-W-Z-mu")>0) :
            process+="set PPFactory:OrderInAlphaS 0\nset PPFactory:OrderInAlphaEW 2\n"
            process+="do PPFactory:Process p p mu+ mu-\ndo PPFactory:Process p p mu+ nu\ndo PPFactory:Process p p mu- nu\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/ChargedLeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/ChargedLeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("W-e")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ nu\nFactory:Process p p e- nu\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
        elif(parameterName.find("W-mu")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p mu+ nu\nFactory:Process p p mu- nu\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
        elif(parameterName.find("Z-e")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ e-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
        elif(parameterName.find("Z-mu")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p mu+ mu-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/ChargedLeptonPairMassCut:MinMass 60*GeV\nset /Herwig/Cuts/ChargedLeptonPairMassCut:MaxMass 120*GeV\n"
        elif(parameterName.find("Z-LowMass-e")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ e-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/ChargedLeptonPairMassCut:MinMass 20*GeV\nset /Herwig/Cuts/ChargedLeptonPairMassCut:MaxMass 70*GeV\n"
        elif(parameterName.find("Z-MedMass-e")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ e-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/ChargedLeptonPairMassCut:MinMass 40*GeV\nset /Herwig/Cuts/ChargedLeptonPairMassCut:MaxMass 130*GeV\n"
        elif(parameterName.find("Z-LowMass-mu")>=0) :
            process+="set Factory:OrderInAlphaS 0\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p mu+ mu-\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/LeptonPairMassScale\n"
            process+="set /Herwig/Cuts/ChargedLeptonPairMassCut:MinMass 20*GeV\nset /Herwig/Cuts/ChargedLeptonPairMassCut:MaxMass 70*GeV\n"
        elif(parameterName.find("W-Jet")>=0) :
            process+="set Factory:OrderInAlphaS 1\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ nu j\ndo Factory:Process p p e- nu j\n\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/HTScale\n"
            process+="set /Herwig/Cuts/QCDCuts:JetFinder /Herwig/Cuts/JetFinder\n"
            process+="insert /Herwig/Cuts/QCDCuts:MultiCuts 0 /Herwig/Cuts/JetCuts\n"
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
            process+="set PPFactory:OrderInAlphaS 1\nset PPFactory:OrderInAlphaEW 2\ndo PPFactory:Process p p e+ e- j\n"
            process+="set Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/HTScale\n"
            process+="set /Herwig/Cuts/QCDCuts:JetFinder /Herwig/Cuts/JetFinder\n"
            process+="insert /Herwig/Cuts/QCDCuts:MultiCuts 0 /Herwig/Cuts/JetCuts\n"
            process+="insert /Herwig/Cuts/JetCuts:JetRegions 0 /Herwig/Cuts/FirstJet\n"
            if(parameterName.find("Z-Jet-1-e")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 100.*GeV\n"
                parameterName=parameterName.replace("Z-Jet-1-e","Z-Jet-e")
            elif(parameterName.find("Z-Jet-2-e")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 190.0*GeV\n"
                parameterName=parameterName.replace("Z-Jet-2-e","Z-Jet-e")
            elif(parameterName.find("Z-Jet-3-e")>=0) :
                process+="set /Herwig/Cuts/FirstJet:PtMin 270.0*GeV\n"
                parameterName=parameterName.replace("Z-Jet-3-e","Z-Jet-e")
        elif(parameterName.find("Z-bb")>=0) :
            process+="set /Herwig/Particles/b:HardProcessMass 4.2*GeV\nset /Herwig/Particles/bbar:HardProcessMass 4.2*GeV\n"
            process+="set Factory:OrderInAlphaS 2\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ e- b bbar\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 91.2*GeV\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/Cuts/MassCut:MinM 66*GeV\nset /Herwig/Cuts/MassCut:MaxM 116*GeV\n"
            process+="set /Herwig/CutsQCDCuts:JetFinder /Herwig/Cuts/JetFinder\n"
            process+="insert  /Herwig/Cuts/QCDCuts:MultiCuts 0  /Herwig/Cuts/JetCuts\n"
            process+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/FirstJet\n"
            process+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/SecondJet\n"
            process+="set  /Herwig/Cuts/FirstJet:PtMin 18.*GeV\n"
            process+="set  /Herwig/Cuts/SecondJet:PtMin 15.*GeV\n"
        elif(parameterName.find("Z-b")>=0) :
            process+="set Factory:OrderInAlphaS 1\nset Factory:OrderInAlphaEW 2\ndo Factory:Process p p e+ e- b\ndo Factory:Process p p e+ e- bbar\n"
            process+="set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 91.2*GeV\nset Factory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/Cuts/MassCut:MinM 66*GeV\nset /Herwig/Cuts/MassCut:MaxM 116*GeV\n"
            process+="set /Herwig/CutsQCDCuts:JetFinder /Herwig/Cuts/JetFinder\n"
            process+="insert  /Herwig/Cuts/QCDCuts:MultiCuts 0  /Herwig/Cuts/JetCuts\n"
            process+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/FirstJet\n"
            process+="set  /Herwig/Cuts/FirstJet:PtMin 15.*GeV\n"
        elif(parameterName.find("W-b")>=0) :
            process += "set /Herwig/Particles/b:HardProcessMass 4.2*GeV\nset /Herwig/Particles/bbar:HardProcessMass 4.2*GeV\n"
            process += "set PPFactory:OrderInAlphaS 2\nset PPFactory:OrderInAlphaEW 2\ndo PPFactory:Process p p e+ nu b bbar\ndo PPFactory:Process p p e- nu b bbar\n"
            process += "set /Herwig/MatrixElements/Matchbox/Scales/FixedScale:FixedScale 80.4*GeV\nset PPFactory:ScaleChoice /Herwig/MatrixElements/Matchbox/Scales/FixedScale\n"
            process+="set /Herwig/CutsQCDCuts:JetFinder /Herwig/Cuts/JetFinder\n"
            process+="insert  /Herwig/Cuts/QCDCuts:MultiCuts 0  /Herwig/Cuts/JetCuts\n"
            process+="insert  /Herwig/Cuts/JetCuts:JetRegions 0  /Herwig/Cuts/FirstJet\n"
            process+="set  /Herwig/Cuts/FirstJet:PtMin 30.*GeV\n"
        else :
            logging.error(" Process %s not supported for Matchbox matrix elements" % name)
            sys.exit(1)
# write the file
with open(os.path.join("Rivet",name+".in") ,'w') as f:
    f.write( template.substitute({ 'process' : process,
                                   'runname' : name,
                                   'parameterFile' : os.path.join(collider,collider+"-"+parameterName+".in") }))






