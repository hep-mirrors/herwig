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
    logging.error("Must specify at least one AIDA histogram file")
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
    if(collider.find("BFactory")<0) :
        templateName= "%s.in" % (collider)
    else :
        templateName= "%s.in" % ("LEP") 
else :
    if(collider.find("BFactory")<0) :
        templateName= "%s-%s.in" % (collider,simulation) 
    else :
        templateName= "%s-%s.in" % ("LEP",simulation) 
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
    process="set EventGenerator:EventHandler:BeamB /Herwig/Particles/pbar-\n"
    if(parameterName.find("Run-II")>=0) :
        process+="set EventGenerator:EventHandler:LuminosityFunction:Energy 1960.0\n"
    elif(parameterName.find("Run-II")>=0) :
        process+="set EventGenerator:EventHandler:LuminosityFunction:Energy 1800.0\n"

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

    process+="set EventGenerator:EventHandler:LuminosityFunction:Energy 7000.0\n"

# write the file
with open(os.path.join("Rivet",name+".in") ,'w') as f:
    f.write( template.substitute({ 'process' : process,
                                   'runname' : name,
                                   'parameterFile' : os.path.join(collider,collider+"-"+parameterName+".in") }))






