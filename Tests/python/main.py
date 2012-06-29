import compare
# directories for the comparison
directory1="trunk/"
directory2="recon/"
# location for the plots
plotLocation="http://projects.hepforge.org/herwig/private/images/peter/comparison1"
# file for the wiki
wiki = open("wiki.info",'w') 
wiki.write("= Comparision of Heriwg++ results =\n")
# lepton-lepton processes
wiki.write("=== Lepton-Lepton ===\n")
wiki.write("|| Process || Cross Section/nb || Error/nb || Cross Section/nb || Error/nb || Fractional Difference sigma || Average Fractional Difference Distribution  || Topdraw || Postscript ||\n")
# quarks
compare.compareLEPQuarks(directory1,directory2,wiki,plotLocation)
# leptons
compare.compareLEPLeptons(directory1,directory2,wiki,plotLocation)
# VH
compare.compareLEPVH(directory1,directory2,wiki,plotLocation)
# VV
compare.compareLEPVV(directory1,directory2,wiki,plotLocation)
# VBF
compare.compareLEPVBF(directory1,directory2,wiki,plotLocation)
# top decay
compare.compareTopDecay(directory1,directory2,wiki,plotLocation)
# charm event shapes
compare.compareCharmShapes(directory1,directory2,wiki,plotLocation) 
# LEP event shapes
compare.compareLEPShapes(directory1,directory2,wiki,plotLocation) 
# hadron-hadron processes at the LHC
wiki.write("=== LHC ===\n")
wiki.write("|| Process || Cross Section/nb || Error/nb || Cross Section/nb || Error/nb || Fractional Difference sigma || Average Fractional Difference Distribution  || Topdraw || Postscript ||\n")
# compare W production
compare.compareW(directory1,directory2,wiki,plotLocation)
# compare Z production
compare.compareZ(directory1,directory2,wiki,plotLocation)
# compare W+jet production
compare.compareWJet(directory1,directory2,wiki,plotLocation)
# compare Z+jet production
compare.compareZJet(directory1,directory2,wiki,plotLocation)
# compare Higgs production
compare.compareHiggs(directory1,directory2,wiki,plotLocation)
# compare Higgs +jet production
compare.compareHiggsJet(directory1,directory2,wiki,plotLocation)
# compare VBF WW
compare.compareWWVBF(directory1,directory2,wiki,plotLocation)
# compare VBF ZZ
compare.compareZZVBF(directory1,directory2,wiki,plotLocation)
# compare VBF all
compare.compareVBF(directory1,directory2,wiki,plotLocation)
# compare WW production
compare.compareWW(directory1,directory2,wiki,plotLocation)
# compare WZ production
compare.compareWZ(directory1,directory2,wiki,plotLocation)
# compare ZZ production
compare.compareZZ(directory1,directory2,wiki,plotLocation)
# compare WGamma production
compare.compareWGamma(directory1,directory2,wiki,plotLocation)
# compare ZGamma production
compare.compareZGamma(directory1,directory2,wiki,plotLocation)
# compare photon pair production
compare.compareGammaGamma(directory1,directory2,wiki,plotLocation)
# compare photon + jet production
compare.compareGammaJet(directory1,directory2,wiki,plotLocation)
# compare WH production
compare.compareWH(directory1,directory2,wiki,plotLocation)
# compare ZH production
compare.compareZH(directory1,directory2,wiki,plotLocation)
# compare QCD production
compare.compareQCD(directory1,directory2,wiki,plotLocation)
# compare QCDFast production
compare.compareQCDFast(directory1,directory2,wiki,plotLocation)
# compare top pair production
compare.compareTop(directory1,directory2,wiki,plotLocation)
# compare bottom pair +Higgs production
compare.comparebbH(directory1,directory2,wiki,plotLocation)
# compare top pair +Higgs production
compare.comparettH(directory1,directory2,wiki,plotLocation)
# compare shower in W
compare.compareWShower(directory1,directory2,wiki,plotLocation)
# compare shower in Z
compare.compareZShower(directory1,directory2,wiki,plotLocation)
# compare shower in W Powheg
compare.compareWPowheg(directory1,directory2,wiki,plotLocation)
# compare shower in Z Powheg
compare.compareZPowheg(directory1,directory2,wiki,plotLocation)
# compare shower in H
compare.compareHJet(directory1,directory2,wiki,plotLocation)
# compare shower in H Powheg
compare.compareHPowheg(directory1,directory2,wiki,plotLocation)
# compare shower in WH
compare.compareWHJet(directory1,directory2,wiki,plotLocation)
# compare shower in ZH
compare.compareZHJet(directory1,directory2,wiki,plotLocation)
# compare shower in WH
compare.compareWHJetPowheg(directory1,directory2,wiki,plotLocation)
# compare shower in ZH
compare.compareZHJetPowheg(directory1,directory2,wiki,plotLocation)
# DIS processes at HERA
wiki.write("=== DIS ===\n")
wiki.write("|| Process || Cross Section/nb || Error/nb || Cross Section/nb || Error/nb || Fractional Difference sigma || Average Fractional Difference Distribution  || Topdraw || Postscript ||\n")
# neutral current
compare.compareNeutralCurrent(directory1,directory2,wiki,plotLocation)
# charged current
compare.compareChargedCurrent(directory1,directory2,wiki,plotLocation)
# photon initiated
wiki.write("=== Photon Initiated ===\n")
wiki.write("|| Process || Cross Section/nb || Error/nb || Cross Section/nb || Error/nb || Fractional Difference sigma || Average Fractional Difference Distribution  || Topdraw || Postscript ||\n")
# fermion production
compare.compareGammaFF(directory1,directory2,wiki,plotLocation)
# W production
compare.compareGammaWW(directory1,directory2,wiki,plotLocation)
# gamma P
compare.compareGammaP(directory1,directory2,wiki,plotLocation)
# close the wiki file
wiki.close()	   
