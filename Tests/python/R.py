#! /usr/bin/env python
from __future__ import print_function
import yoda,os,math,subprocess,optparse
from string import Template
# get the path for the rivet data
p = subprocess.Popen(["rivet-config", "--datadir"],stdout=subprocess.PIPE)
path=p.communicate()[0].strip()
thresholds = [2.*1.87,2.*5.28]
#Define the arguments
op = optparse.OptionParser(usage=__doc__)
op.add_option("--process"         , dest="processes"       , default=[], action="append")
op.add_option("--path"            , dest="path"            , default=path)
op.add_option("--non-perturbative", dest="nonPerturbative" , default=False, action="store_true")
op.add_option("--perturbative"    , dest="perturbative"    , default=False, action="store_true")
op.add_option("--dest"            , dest="dest"            , default="Rivet")
op.add_option("--list"            , dest="list"            , default=False, action="store_true")
op.add_option("--max-energy"      , dest="maxEnergy"       , default=11.5, type="float", action="store")
op.add_option("--plots"           , dest="plot"           , default=False, action="store_true")
opts, args = op.parse_args()
path=opts.path
# list of analyses
analyses={}
analyses["CLEO_2007_I753556"] = ["d01-x01-y01"]
analyses["DASP_1978_I129715"] = ["d01-x01-y01"]
analyses["CLEO_1984_I193577"] = ["d01-x01-y01"]
analyses["AMY_1990_I294525" ] = ["d01-x01-y01"]
analyses["BABAR_2009_I797507"] = ["d01-x01-y01"]
analyses["BELLE_2015_I1336624"] = ["d02-x01-y01","d01-x01-y01","d01-x01-y02","d01-x01-y03"]
analyses["TASSO_1982_I176887"] = ["d01-x01-y01","d02-x01-y01","d03-x01-y01"]
analyses["CRYSTAL_BALL_1986_I238081"] = ["d01-x01-y01"]
analyses["CRYSTAL_BALL_1990_I294419"] = ["d01-x01-y01","d02-x01-y01"]
analyses["CRYSTAL_BALL_1988_I261078"] = ["d01-x01-y01"]
analyses["TOPAZ_1990_I283003"] = ["d01-x01-y01"]
analyses["TOPAZ_1993_I353845"] = ["d02-x01-y01"]
analyses["TOPAZ_1995_I381777"] = ["d01-x01-y01"]
analyses["VENUS_1987_I251274"] = ["d01-x01-y01"]
analyses["VENUS_1990_I283774"] = ["d01-x01-y01"]
analyses["VENUS_1990_I296392"] = ["d03-x01-y01"]
analyses["VENUS_1999_I500179"] = ["d01-x01-y01"]
analyses["MD1_1986_I364141"] = ["d01-x01-y01"]
analyses["MUPI_1972_I84978"] = ["d01-x01-y01"]
analyses["MUPI_1973_I95215"] = ["d01-x01-y01"]
analyses["NMD_1974_I745"   ] = ["d01-x01-y01","d01-x01-y02"]
analyses["GAMMAGAMMA_1979_I141722"] = ["d01-x01-y01","d02-x01-y01","d02-x01-y02"]
analyses["MARKI_1975_I100733"] = ["d01-x01-y01","d02-x01-y01"]
analyses["MARKI_1977_I119979"] = ["d01-x01-y01","d02-x01-y01"]
analyses["MARKI_1976_I108144"] = ["d01-x01-y01"]
analyses["DASP_1982_I178613"] = ["d02-x01-y01"]
analyses["DESY147_1980_I153896"] = ["d01-x01-y01"]
analyses["DESY147_1978_I131524"] = ["d01-x01-y01","d02-x01-y01"]
analyses["CLEO_1983_I188805"] = ["d01-x01-y01"]
analyses["CLEO_1998_I445351"] = ["d01-x01-y01"]
analyses["CLEO_1983_I188803"] = ["d02-x01-y01"]
analyses["CLEOC_2008_I777917"] = ["d06-x01-y01","d06-x01-y02"]
analyses["CLEO_2006_I700665"] = ["d01-x01-y01"]
analyses["CUSB_1982_I180613"] = ["d03-x01-y01"]
analyses["MAC_1985_I206052"] = ["d01-x01-y01"]
analyses["MARKII_1991_I295286"] = ["d01-x01-y01"]
analyses["MARKJ_1979_I141976"] = ["d01-x01-y01"]
analyses["MARKJ_1980_I158857"] = ["d01-x01-y01"]
analyses["MARKJ_1982_I166369"] = ["d01-x01-y01"]
analyses["MARKJ_1983_I182337"] = ["d04-x01-y01"]
analyses["MARKJ_1984_I196567"] = ["d01-x01-y01"]
analyses["MARKJ_1986_I230297"] = ["d01-x01-y01"]
analyses["LENA_1982_I179431"] = ["d01-x01-y01","d02-x01-y01","d03-x01-y01"]
analyses["MARKII_1979_I143939"] = ["d02-x01-y01","d03-x01-y01"]
analyses["ARGUS_1992_I319102"] = ["d01-x01-y01"]
analyses["CELLO_1981_I166365"] = ["d01-x01-y01"]
analyses["CELLO_1984_I202783"] = ["d02-x01-y01"]
analyses["CELLO_1987_I236981"] = ["d01-x01-y01"]
analyses["TASSO_1979_I140303"] = ["d01-x01-y01"]
analyses["TASSO_1980_I143690"] = ["d01-x01-y01"]
analyses["TASSO_1984_I199468"] = ["d01-x01-y01","d02-x01-y01"]
analyses["JADE_1987_I234905"] = ["d01-x01-y01"]
analyses["PLUTO_1982_I166799"] = ["d01-x01-y01","d02-x01-y01"]
analyses["PLUTO_1979_I142517"] = ["d01-x01-y01"]
analyses["PLUTO_1979_I140818"] = ["d02-x01-y01"]
analyses["PLUTO_1979_I140294"] = ["d01-x01-y01","d02-x01-y01"]
analyses["PLUTO_1980_I152291"] = ["d01-x01-y01","d02-x01-y01"]
analyses["BBAR_1980_I152630"] = ["d01-x01-y01"]
analyses["BESII_2000_I505323"] = ["d01-x01-y01"]
analyses["BESII_2002_I552757"] = ["d01-x01-y01"]
analyses["BES_1995_I39870"] = ["d01-x01-y01"]
analyses["BESII_2009_I814778"] = ["d01-x01-y01"]
analyses["BESII_2006_I717665"] = ["d02-x01-y01"]
analyses["GAMMAGAMMA_1979_I133588"] = ["d01-x01-y01","d01-x01-y02"]
analyses["GAMMAGAMMA_1975_I100016"] = ["d01-x01-y01","d01-x01-y02"]
analyses["GAMMAGAMMA_1973_I84794"] = ["d03-x01-y01","d04-x01-y01"]
analyses["GAMMAGAMMA_1981_I158474"] = ["d02-x01-y01","d02-x01-y02","d01-x01-y05","d01-x01-y06"]
analyses["TASSO_1984_I195333"] = ["d01-x01-y01","d04-x01-y01"]
analyses["PLUTO_1977_I110272"]=["d02-x01-y01"]
analyses["FENICE_1996_I426675"]=["d01-x01-y01"]
analyses["PLUTO_1981_I165122"]=["d01-x01-y01","d02-x01-y01","d03-x01-y01"]
analyses["KEDR_2019_I1673357"]=["d01-x01-y01","d01-x01-y02"]
analyses["BELLE_2011_I878228"] = ["d02-x01-y01"]
# list analyses if needed
if(opts.list) :
    print (" ".join(analyses.keys()))
    quit()
if(opts.plot) :
    output=""
    for analysis in analyses.keys():
        for plot in analyses[analysis]:
            output+= " -m/%s/%s" % (analysis,plot)
    for i in range(1,7) :
        output += " -m/BESII_2004_I622224/d0%s-x01-y01" % i
    print (output)
    quit()

energies={}
def nearestEnergy(en) :
    Emin=0
    delta=1e30
    anals=[]
    for val in energies :
        if(abs(val-en)<delta) :
            delta = abs(val-en)
            Emin = val
            anals=energies[val]
    return (Emin,delta,anals)

for analysis in analyses :
    aos=yoda.read(os.path.join(os.path.join(os.getcwd(),path.decode('UTF-8')),analysis+".yoda"))
    if(len(aos)==0) : continue
    for plot in analyses[analysis] :
        histo = aos["/REF/%s/%s" %(analysis,plot)]
        for point in histo.points() :
            energy = point.x()
            if(analysis=="FENICE_1996_I426675") :
                energy = math.sqrt(energy)
            if(energy>200) :
                energy *= 0.001
            emin,delta,anals = nearestEnergy(energy)
            if(energy in energies) :
                if(analysis not in energies[energy]) :
                    energies[energy].append(analysis)
            elif(delta<1e-7) :
                if(analysis not in anals) :
                    anals.append(analysis)
            else :
                energies[energy]=[analysis]

# add the bes spectra
for val in [2.2,2.6,3.0,3.2,4.6,4.8] :
    if val in energies :
        energies[val].append("BESII_2004_I622224")
    else:
        energies[val] = ["BESII_2004_I622224"]
        
# set up the templates
with open("python/LowEnergy-EE-Perturbative.in", 'r') as f:
    templateText = f.read()
perturbative=Template( templateText )
with open("python/LowEnergy-EE-NonPerturbative.in", 'r') as f:
    templateText = f.read()
nonPerturbative=Template( templateText )


# lepton matrix element
lepton_me="insert SubProcess:MatrixElements 0 MEee2gZ2ll\nset MEee2gZ2ll:Allowed Muon\n"
# low energy matrix element
mes = ["MEee2Pions", "MEee2Kaons", "MEee3Pions", "MEee4Pions", "MEee2EtaPiPi",
       "MEee2EtaPrimePiPi", "MEee2EtaPhi", "MEee2EtaOmega", "MEee2OmegaPi",
       "MEee2OmegaPiPi", "MEee2PhiPi", "MEee2PiGamma", "MEee2EtaGamma",
       "MEee2PPbar", "MEee2LL"   , "MEee2KKPi" ]
proc=lepton_me
for matrix in mes :
    proc+="insert SubProcess:MatrixElements 0 %s\n" % matrix

targets=""
for energy in sorted(energies) :
    anal=""
    for analysis in energies[energy]: 
        anal+= "insert /Herwig/Analysis/Rivet:Analyses 0 %s\n" %analysis
    # input file for perturbative QCD
    maxflavour =5
    if(energy<thresholds[0]) :
        maxflavour=3
    elif(energy<thresholds[1]) :
        maxflavour=4
    if(opts.perturbative) :
        inputPerturbative = perturbative.substitute({"ECMS" : "%8.6f" % energy, "ANALYSES" : anal,
                                                     "lepton" : lepton_me, 'maxflavour' : maxflavour})
        with open(opts.dest+"/Rivet-LowEnergy-EE-Perturbative-%8.6f.in" % energy ,'w') as f:
            f.write(inputPerturbative)
        targets += "Rivet-LowEnergy-EE-Perturbative-%8.6f.yoda " % energy
    # input file for currents
    if(opts.nonPerturbative and energy <= opts.maxEnergy) :
        inputNonPerturbative = nonPerturbative.substitute({"ECMS" : "%8.6f" % energy, "ANALYSES" : anal,
                                                           "processes" : proc})
        with open(opts.dest+"/Rivet-LowEnergy-EE-NonPerturbative-%8.6f.in" % energy ,'w') as f:
            f.write(inputNonPerturbative)
        targets += "Rivet-LowEnergy-EE-NonPerturbative-%8.6f.yoda " % energy
print (targets)
