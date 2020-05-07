#! /usr/bin/env python
from __future__ import print_function
import yoda,os,math,subprocess,optparse
from string import Template
# get the path for the rivet data
p = subprocess.Popen(["rivet-config", "--datadir"],stdout=subprocess.PIPE)
path=p.communicate()[0].strip().decode("UTF-8")
#Define the arguments
op = optparse.OptionParser(usage=__doc__)
op.add_option("--process"         , dest="processes"       , default=[], action="append")
op.add_option("--path"            , dest="path"            , default=path)
op.add_option("--non-perturbative", dest="nonPerturbative" , default=False, action="store_true")
op.add_option("--perturbative"    , dest="perturbative"    , default=False, action="store_true")
op.add_option("--dest"            , dest="dest"            , default="Rivet")
op.add_option("--list"            , dest="list"            , default=False, action="store_true")
op.add_option("--flavour"         , dest="flavour"         , default="All"  )
op.add_option("--plots"           , dest="plot"           , default=False, action="store_true")
opts, args = op.parse_args()
path=opts.path
thresholds = [0.7,2.*.5,2.*1.87,2.*5.28]
# the list of analyses and processes
analyses = { 'KK'           : {}, 'PiPi'      : {}, 'PPbar'   : {}, "3Pi"      : {},
             "EtaprimePiPi" : {}, "4Pi"       : {}, "EtaPhi"  : {}, "EtaOmega" : {},
             "2K1Pi"        : {}, "2K2Pi"     : {}, "4K"      : {}, "6m"       : {},
             "EtaPiPi"      : {}, "OmegaPi"   : {}, "PiGamma" : {}, "EtaGamma" : {},
             "PhiPi"        : {}, "OmegaPiPi" : {}, "DD"      : {}, "BB"       : {},
             "5Pi"          : {}, "LL"        : {}, "Baryon"  : {} }
# pi+pi-
analyses["PiPi"]["KLOE_2009_I797438"   ] = ["d02-x01-y01"]
analyses["PiPi"]["KLOE_2005_I655225"   ] = ["d02-x01-y01"]
analyses["PiPi"]["KLOE2_2017_I1634981" ] = ["d01-x01-y01"]
analyses["PiPi"]["BABAR_2009_I829441"  ] = ["d01-x01-y01"]
analyses["PiPi"]["DM1_1978_I134061"    ] = ["d01-x01-y01"]
analyses["PiPi"]["DM2_1989_I267118"    ] = ["d01-x01-y01"]
analyses["PiPi"]["CMD2_2007_I728302"   ] = ["d02-x01-y01"]
analyses["PiPi"]["CMD2_2006_I728191"   ] = ["d03-x01-y01"]
analyses["PiPi"]["BESIII_2016_I1385603"] = ["d01-x01-y01"]
analyses["PiPi"]["SND_2005_I686349"    ] = ["d01-x01-y01"]
analyses["PiPi"]["CMD_1985_I221309"    ] = ["d01-x01-y01","d02-x01-y01"]
analyses["PiPi"]["CMD2_2002_I568807"   ] = ["d01-x01-y02"]
analyses["PiPi"]["CMD2_1999_I498859"   ] = ["d01-x01-y01"]
analyses['PiPi']["CLEOC_2005_I693873"  ] = ["d01-x01-y01"]
analyses['PiPi']["ND_1991_I321108"     ] = ["d11-x01-y01"]
analyses['PiPi']["OLYA_1984_I208231"   ] = ["d01-x01-y01"]
# K+K- and K_S^0 K_L^0
analyses['KK']["BESIII_2018_I1704558"] = ["d01-x01-y01"]
analyses['KK']["BABAR_2013_I1238807" ] = ["d01-x01-y01"]
analyses['KK']["DM1_1981_I156053"    ] = ["d01-x01-y01"]
analyses['KK']["DM1_1981_I156054"    ] = ["d01-x01-y01"]
analyses['KK']["CLEOC_2005_I693873"  ] = ["d01-x01-y02"]
analyses['KK']["BABAR_2015_I1383130" ] = ["d01-x01-y04"]
analyses['KK']["DM2_1988_I262690"    ] = ["d01-x01-y01"]
analyses['KK']["SND_2007_I755881"    ] = ["d01-x01-y01"]
analyses['KK']["SND_2001_I533574"    ] = ["d01-x01-y01","d01-x01-y02","d01-x01-y03",
                                          "d02-x01-y01","d02-x01-y02","d02-x01-y03"]
analyses['KK']["SND_2006_I720035"    ] = ["d01-x01-y01"]
analyses['KK']["BABAR_2014_I1287920" ] = ["d09-x01-y01"]
analyses['KK']["CMD2_2003_I601222"   ] = ["d01-x01-y01"]
analyses['KK']["CMD3_2016_I1444990"  ] = ["d01-x01-y06"]
analyses['KK']["CMD2_1995_I406880"   ] = ["d01-x01-y01","d01-x01-y02"]
analyses['KK']["CMD2_1999_I502164"   ] = ["d01-x01-y01","d02-x01-y01",
                                          "d03-x01-y01","d04-x01-y01"]
analyses['KK']["CMD2_2008_I782516"   ] = ["d01-x01-y01","d02-x01-y01"]
analyses['KK']["ND_1991_I321108"     ] = ["d12-x01-y01","d13-x01-y01"]
analyses['KK']["OLYA_1981_I173076"   ] = ["d01-x01-y01"]
analyses['KK']["SND_2016_I1484677"   ] = ["d01-x01-y01","d02-x01-y01"]
# proton-antiproton
analyses['PPbar']["BESIII_2019_I1736235"] = ["d01-x01-y01"]
analyses['PPbar']["BESIII_2019_I1718337"] = ["d01-x01-y01"]
analyses['PPbar']["BESIII_2015_I1358937"] = ["d01-x01-y05"]
analyses['PPbar']["BABAR_2013_I1217421" ] = ["d01-x01-y01"]
analyses['PPbar']["BABAR_2013_I1247058" ] = ["d01-x01-y01"]
analyses['PPbar']["SND_2014_I1321689"   ] = ["d01-x01-y01","d02-x01-y01"]
analyses['PPbar']["CMD3_2016_I1385598"  ] = ["d01-x01-y06"]
analyses['PPbar']["CLEOC_2005_I693873"  ] = ["d01-x01-y03"]
analyses['PPbar']["BABAR_2006_I700020"  ] = ["d01-x01-y01","d02-x01-y01"]
analyses['PPbar']["DM2_1983_I190558"    ] = ["d01-x01-y01"]
analyses["PPbar"]["DM2_1990_I297706"    ] = ["d01-x01-y01"]
analyses["PPbar"]["DM1_1979_I141565"    ] = ["d01-x01-y01"]
analyses["PPbar"]["FENICE_1998_I471263" ] = ["d01-x01-y01"]
analyses["PPbar"]["FENICE_1994_I377833" ] = ["d01-x01-y01"]
analyses['PPbar']["BESII_2005_I685906"  ] = ["d01-x01-y01"]
analyses['PPbar']["BESIII_2014_I1286898"] = ["d01-x01-y06"]
# pi0 gamma
analyses["PiGamma"]["SND_2018_I1694988"] = ["d01-x01-y01"]
analyses["PiGamma"]["SND_2016_I1418483"] = ["d01-x01-y05"]
analyses["PiGamma"]["SND_2003_I612867" ] = ["d01-x01-y01"]
analyses["PiGamma"]["CMD2_2005_I658856"] = ["d02-x01-y01"]
analyses["PiGamma"]["SND_2000_I524221" ] = ["d01-x01-y02"]
# eta gamma
analyses["EtaGamma"]["CMD2_2005_I658856" ] = ["d01-x01-y01"]
analyses["EtaGamma"]["SND_2006_I717778"  ] = ["d01-x01-y01","d02-x01-y01"]
analyses["EtaGamma"]["SND_2014_I1275333" ] = ["d01-x01-y01"]
analyses["EtaGamma"]["SND_2000_I524221"  ] = ["d01-x01-y01"]
analyses["EtaGamma"]["CMD2_1999_I503154" ] = ["d01-x01-y01"]
analyses["EtaGamma"]["CMD2_2001_I554522" ] = ["d01-x01-y01"]
analyses['EtaGamma']["CMD2_1995_I406880" ] = ["d01-x01-y04"]
analyses['EtaGamma']["BABAR_2006_I716277"] = ["d01-x01-y01"]
# 3 pion
analyses["3Pi"]["BABAR_2004_I656680"     ] = ["d01-x01-y01"]
analyses["3Pi"]["BESIII_2019_I1773081"   ] = ["d01-x01-y01"]
analyses["3Pi"]["SND_2002_I582183"       ] = ["d01-x01-y01"]
analyses["3Pi"]["SND_2003_I619011"       ] = ["d01-x01-y01"]
analyses["3Pi"]["SND_1999_I508003"       ] = ["d01-x01-y01"]
analyses["3Pi"]["SND_2001_I533574"       ] = ["d01-x01-y04","d02-x01-y04"]
analyses["3Pi"]["CMD2_2000_I523691"      ] = ["d01-x01-y01"]
analyses["3Pi"]["CMD2_1998_I480170"      ] = ["d01-x01-y01"]
analyses['3Pi']["CMD2_1995_I406880"      ] = ["d01-x01-y03"]
analyses['3Pi']["DM2_1992_I339265"       ] = ["d01-x01-y01"]
analyses['3Pi']["DM1_1980_I140174"       ] = ["d01-x01-y01"]
analyses['3Pi']["ND_1991_I321108"        ] = ["d05-x01-y01","d10-x01-y04"]
analyses['3Pi']["GAMMAGAMMA_1981_I158474"] = ["d01-x01-y01"]
analyses["3Pi"]["CLEO_2006_I691720"      ] = ["d01-x01-y01"]
analyses["3Pi"]["SND_2015_I1389908"      ] = ["d01-x01-y01"]
# eta pipi
analyses["EtaPiPi"]["BABAR_2018_I1700745"] = ["d01-x01-y01","d02-x01-y01"]
analyses["EtaPiPi"]["BABAR_2018_I1647139"] = ["d01-x01-y01"]
analyses["EtaPiPi"]["SND_2015_I1332929"  ] = ["d01-x01-y01"]
analyses["EtaPiPi"]["SND_2018_I1638368"  ] = ["d01-x01-y01"]
analyses["EtaPiPi"]["BABAR_2007_I758568" ] = ["d01-x01-y01","d02-x01-y01"]
analyses["EtaPiPi"]["CMD2_2000_I532970"  ] = ["d02-x01-y01"]
analyses["EtaPiPi"]["DM2_1988_I264144"   ] = ["d01-x01-y01"]
analyses['EtaPiPi']["ND_1991_I321108"    ] = ["d06-x01-y01","d14-x01-y01"]
analyses['EtaPiPi']["CMD3_2019_I1744510" ] = ["d02-x01-y01"]
# eta' pipi
analyses["EtaprimePiPi"]["BABAR_2007_I758568"] = ["d05-x01-y01","d06-x01-y01"]
# Eta Phi
analyses["EtaPhi"]["BABAR_2008_I765258"  ] = ["d04-x01-y01","d05-x01-y01"]
analyses["EtaPhi"]["BABAR_2007_I758568"  ] = ["d08-x01-y01","d09-x01-y01"]
analyses["EtaPhi"]["SND_2018_I1693737"   ] = ["d01-x01-y01"]
analyses["EtaPhi"]["BABAR_2017_I1511276" ] = ["d03-x01-y01"]
analyses["EtaPhi"]["SND_2019_I1726419"   ] = ["d01-x01-y01","d01-x01-y03"]
analyses["EtaPhi"]["CMD3_2019_I1740541"  ] = ["d01-x01-y06","d02-x01-y06","d03-x01-y06"]
analyses["EtaPhi"]["CMD3_2017_I1606078"  ] = ["d01-x01-y01"]
analyses["EtaPhi"]["BABAR_2006_I709730"  ] = ["d02-x01-y01"]
analyses["EtaPhi"]["BESII_2008_I801210"  ] = ["d01-x01-y03"]
analyses["EtaPhi"]["BABAR_2006_I731865"  ] = ["d01-x01-y02"]
analyses["EtaPhi"]["BELLE_2009_I823878"  ] = ["d01-x01-y01"]
# Eta Omega
analyses["EtaOmega"]["SND_2016_I1473343" ] = ["d01-x01-y01"]
analyses["EtaOmega"]["BABAR_2006_I709730"] = ["d02-x01-y01"]
analyses["EtaOmega"]["SND_2019_I1726419" ] = ["d01-x01-y01","d01-x01-y02"]
analyses["EtaOmega"]["CMD3_2017_I1606078"] = ["d01-x01-y01","d01-x01-y02"]
analyses["EtaOmega"]["BESII_2008_I801210"  ] = ["d01-x01-y03"]
# 4 pions
analyses["4Pi"]["BABAR_2017_I1621593"    ] = ["d01-x01-y01","d02-x01-y01"]
analyses["4Pi"]["BABAR_2012_I1086164"    ] = ["d01-x01-y01"]
analyses["4Pi"]["CMD2_2000_I531667"      ] = ["d01-x01-y01"]
analyses["4Pi"]["CMD2_2004_I648023"      ] = ["d01-x01-y01"]
analyses["4Pi"]["BABAR_2005_I676691"     ] = ["d01-x01-y01"]
analyses["4Pi"]["CMD2_2000_I511375"      ] = ["d01-x01-y01"]
analyses["4Pi"]["CMD2_1999_I483994"      ] = ["d01-x01-y01","d02-x01-y01","d03-x01-y01"]
analyses["4Pi"]["BESII_2008_I801210"     ] = ["d01-x01-y01"]
analyses["4Pi"]["KLOE_2008_I791841"      ] = ["d01-x01-y01"]
analyses['4Pi']["ND_1991_I321108"        ] = ["d07-x01-y01","d08-x01-y01","d10-x01-y01","d10-x01-y02",
                                              "d01-x01-y01","d02-x01-y01","d03-x01-y01","d04-x01-y01","d10-x01-y03"]
analyses['4Pi']["BESII_2007_I750713"     ] = ["d01-x01-y03"]
analyses['4Pi']["SND_2001_I579319"       ] = ["d01-x01-y01","d02-x01-y01"]
analyses['4Pi']["DM1_1982_I168552"       ] = ["d01-x01-y01"]
analyses['4Pi']["DM1_1979_I132828"       ] = ["d01-x01-y01"]
analyses['4Pi']["GAMMAGAMMA_1980_I153382"] = ["d01-x01-y01"]
analyses['4Pi']["GAMMAGAMMA_1981_I158474"] = ["d01-x01-y02"]
# (these are Omega(-> pi0 gamma) pi0)
analyses["OmegaPi"]["SND_2016_I1489182"  ] = ["d01-x01-y01"]
analyses["OmegaPi"]["SND_2000_I527752"   ] = ["d01-x01-y01"]
analyses["OmegaPi"]["SND_2000_I503946"   ] = ["d01-x01-y01"]
analyses["OmegaPi"]["CMD2_2003_I616446"  ] = ["d01-x01-y01"]
# non Omega
analyses["OmegaPi"]["SND_2002_I587084"  ] = ["d01-x01-y01"]
analyses["OmegaPi"]["CMD2_2004_I630009" ] = ["d01-x01-y01"]
analyses["OmegaPi"]["KLOE_2008_I791841" ] = ["d02-x01-y01"]
# from 4 Pion
analyses["OmegaPi"]["CMD2_1999_I483994" ] = ["d03-x01-y01"]
analyses['OmegaPi']["ND_1991_I321108"   ] = ["d01-x01-y01","d02-x01-y01","d03-x01-y01",
                                             "d04-x01-y01","d10-x01-y03"]
# 5 pion and related
analyses["OmegaPiPi"]["DM1_1981_I166964"   ] = ["d01-x01-y01"]
analyses["OmegaPiPi"]["DM2_1992_I339265"   ] = ["d02-x01-y01"]
analyses["OmegaPiPi"]["CMD2_2000_I532970"  ] = ["d01-x01-y01"]
analyses["OmegaPiPi"]["BABAR_2018_I1700745"] = ["d01-x01-y01","d03-x01-y01"]
analyses["OmegaPiPi"]["BABAR_2007_I758568" ] = ["d01-x01-y01","d03-x01-y01","d04-x01-y01"]
analyses['OmegaPiPi']["ND_1991_I321108"    ] = ["d14-x01-y01"]
analyses["5Pi"]["CMD2_2000_I532970"        ] = ["d03-x01-y01"]
analyses["5Pi"]["BABAR_2007_I758568"       ] = ["d01-x01-y01"]
analyses['5Pi']["ND_1991_I321108"          ] = ["d14-x01-y01"]
analyses['5Pi']["GAMMAGAMMA_1981_I158474"  ] = ["d01-x01-y03"]
analyses["5Pi"]["BABAR_2018_I1700745"      ] = ["d01-x01-y01"]
# 2K 1 pi
analyses["2K1Pi"]["BABAR_2008_I765258"  ] = ["d01-x01-y01","d02-x01-y01"]
analyses["2K1Pi"]["DM1_1982_I176801"    ] = ["d01-x01-y01"]
analyses["2K1Pi"]["DM2_1991_I318558"    ] = ["d01-x01-y01","d02-x01-y01"]
analyses["2K1Pi"]["BESII_2008_I801208"  ] = ["d01-x01-y01"]
analyses["2K1Pi"]["SND_2018_I1637194"   ] = ["d01-x01-y01"]
analyses["2K1Pi"]["BESIII_2018_I1691798"] = ["d01-x01-y01"]
analyses["2K1Pi"]["BABAR_2017_I1511276" ] = ["d01-x01-y01"]
analyses["PhiPi"]["BABAR_2017_I1511276" ] = ["d01-x01-y01","d02-x01-y01"]
analyses["PhiPi"]["BABAR_2008_I765258"  ] = ["d02-x01-y01","d03-x01-y01"]
# 2K 2 pi
analyses["2K2Pi"]["DM1_1982_I169382"    ] = ["d01-x01-y01"]
analyses["2K2Pi"]["BABAR_2005_I676691"  ] = ["d02-x01-y01"]
analyses["2K2Pi"]["BABAR_2014_I1287920" ] = ["d09-x01-y01","d10-x01-y01","d11-x01-y01"]
analyses["2K2Pi"]["BABAR_2012_I892684"  ] = ["d01-x01-y01","d02-x01-y01","d03-x01-y01",
                                             "d04-x01-y01","d05-x01-y01",
                                             "d06-x01-y01","d07-x01-y01"]
analyses["2K2Pi"]["BABAR_2007_I747875"  ] = ["d01-x01-y01","d02-x01-y01","d03-x01-y01",
                                             "d04-x01-y01","d05-x01-y01","d07-x01-y01"]
analyses["2K2Pi"]["BESII_2008_I801210"  ] = ["d01-x01-y02"]
analyses["2K2Pi"]["BESII_2008_I801208"  ] = ["d01-x01-y02"]
analyses["2K2Pi"]["BELLE_2009_I809630"  ] = ["d01-x01-y01"]
analyses["2K2Pi"]["CMD3_2016_I1395968"  ] = ["d01-x01-y01"]
analyses['2K2Pi']["BESII_2007_I750713"  ] = ["d01-x01-y04"]
analyses["2K2Pi"]["BABAR_2017_I1511276" ] = ["d03-x01-y01","d04-x01-y01"]
analyses["2K2Pi"]["BABAR_2017_I1591716" ] = ["d01-x01-y01","d02-x01-y01"]
analyses['2K2Pi']["BESIII_2018_I1699641"] = ["d01-x01-y01","d02-x01-y01"]
analyses['2K2Pi']["CMD3_2019_I1770428"  ] = ["d01-x01-y06"]
# 4K
analyses["4K"]["BESIII_2019_I1743841"] = ["d01-x01-y01","d02-x01-y01"]
analyses["4K"]["BABAR_2005_I676691"  ] = ["d03-x01-y01"]
analyses["4K"]["BABAR_2014_I1287920" ] = ["d12-x01-y01"]
analyses["4K"]["BABAR_2012_I892684"  ] = ["d08-x01-y01"]
analyses["4K"]["BABAR_2007_I747875"  ] = ["d07-x01-y01"]
analyses['4K']["BESII_2007_I750713"  ] = ["d01-x01-y06","d01-x01-y07"]
# 6 mesons
analyses["6m"]["CMD3_2013_I1217420" ] = ["d01-x01-y01"]
analyses["6m"]["SND_2019_I1726419"  ] = ["d01-x01-y01","d01-x01-y04"]
analyses["6m"]["CMD3_2017_I1606078" ] = ["d01-x01-y03","d01-x01-y04"]
analyses["6m"]["CMD3_2019_I1720610" ] = ["d01-x01-y01","d01-x01-y02","d01-x01-y03"]
analyses["6m"]["BABAR_2018_I1700745"] = ["d04-x01-y01","d05-x01-y01"]
analyses["6m"]["SND_2016_I1471515"  ] = ["d01-x01-y06"]
analyses["6m"]["DM1_1981_I166353"   ] = ["d01-x01-y01"]
analyses["6m"]["BABAR_2006_I709730" ] = ["d01-x01-y01","d02-x01-y01","d03-x01-y01"]
analyses["6m"]["BABAR_2007_I758568" ] = ["d05-x01-y01","d07-x01-y01",
                                        "d08-x01-y01","d09-x01-y01","d10-x01-y01","d11-x01-y01"]
analyses["6m"]["BESII_2007_I763880" ] = ["d01-x01-y01","d01-x01-y02","d01-x01-y03","d01-x01-y04",
                                         "d01-x01-y05","d01-x01-y06","d01-x01-y07"]
analyses["6m"]["BESII_2007_I762901" ] = ["d01-x01-y01","d01-x01-y02","d01-x01-y03","d01-x01-y04",
                                         "d01-x01-y06","d01-x01-y07","d01-x01-y08","d01-x01-y09","d01-x01-y10"]
analyses["6m"]["CLEO_2006_I691720"  ] = ["d01-x01-y02","d01-x01-y03","d01-x01-y04","d01-x01-y05",
                                         "d01-x01-y07","d01-x01-y08","d01-x01-y09","d01-x01-y10","d01-x01-y11",
                                         "d01-x01-y12","d01-x01-y13","d01-x01-y14","d01-x01-y15","d01-x01-y17"]
analyses["6m"]["BESII_2008_I801210" ] = ["d01-x01-y03","d01-x01-y04","d01-x01-y05"]
analyses["6m"]["BESII_2008_I801208" ] = ["d01-x01-y03","d01-x01-y04","d01-x01-y05","d01-x01-y06"]
analyses["6m"]["MARKI_1982_I169326" ] = ["d06-x01-y01"]
analyses["6m"]["MARKI_1975_I100592" ] = ["d01-x01-y01","d02-x01-y01"]
analyses['6m']["BESII_2007_I750713" ] = ["d01-x01-y08","d01-x01-y09","d01-x01-y11",
                                        "d01-x01-y12","d01-x01-y13","d01-x01-y14",
                                         "d01-x01-y15","d01-x01-y16","d01-x01-y17","d01-x01-y18"]
analyses['6m']["SND_2016_I1473343"  ] = ["d01-x01-y01"]
# other baryon processes
analyses['Baryon']["BESIII_2017_I1509241"  ] = ["d01-x01-y01"]
# DD
analyses["DD"]["BELLE_2007_I723333"       ] = ["d01-x01-y01","d02-x01-y01"]
analyses["DD"]["BELLE_2007_I756012"       ] = ["d01-x01-y01"]
analyses["DD"]["BELLE_2007_I756643"       ] = ["d01-x01-y01"]
analyses["DD"]["BELLE_2008_I757220"       ] = ["d01-x01-y01","d02-x01-y01"]
analyses["DD"]["BELLE_2008_I759073"       ] = ["d01-x01-y01"]
analyses["DD"]["BABAR_2008_I776519"       ] = ["d01-x01-y01","d01-x01-y02"]
analyses["DD"]["BELLE_2008_I791660"       ] = ["d01-x01-y01"]
analyses["DD"]["BELLE_2013_I1225975"      ] = ["d01-x01-y01"]
analyses["DD"]["BELLE_2014_I1282602"      ] = ["d01-x01-y01"]
analyses["DD"]["BELLE_2015_I1324785"      ] = ["d01-x01-y01"]
analyses["DD"]["BESIII_2016_I1457597"     ] = ["d01-x01-y07"]
analyses["DD"]["BESIII_2015_I1355215"     ] = ["d01-x01-y10"]
analyses["DD"]["BESIII_2015_I1377204"     ] = ["d01-x01-y10"]
analyses["DD"]["BESIII_2016_I1495838"     ] = ["d01-x01-y01","d02-x01-y01"]
analyses["DD"]["CRYSTAL_BALL_1986_I238081"] = ["d02-x01-y01"]
analyses["DD"]["CLEOC_2008_I777917"       ] = ["d01-x01-y01","d01-x01-y02","d01-x01-y03",
                                               "d02-x01-y01","d02-x01-y02","d02-x01-y03",
                                               "d03-x01-y01","d03-x01-y02","d03-x01-y03",
                                               "d04-x01-y01","d04-x01-y02",
                                               "d05-x01-y01","d05-x01-y02"]
analyses["DD"]["BELLE_2017_I1613517"      ] = ["d01-x01-y01","d01-x01-y02"]
analyses["DD"]["BESIII_2014_I1323621"     ] = ["d01-x01-y01"]
analyses["DD"]["BESIII_2015_I1406939"     ] = ["d02-x01-y06","d03-x01-y06"]
analyses["DD"]["BESIII_2017_I1628093"     ] = ["d01-x01-y01"]
analyses["DD"]["BESIII_2019_I1723934"     ] = ["d01-x01-y01"]
analyses["DD"]["BESIII_2019_I1756876"     ] = ["d01-x01-y09","d01-x01-y10"]
analyses["DD"]["BABAR_2007_I729388"       ] = ["d01-x01-y01"]
analyses["DD"]["BESIII_2015_I1329785"     ] = ["d01-x01-y08","d02-x01-y08","d03-x01-y08"]
analyses["DD"]["BESIII_2017_I1494065"     ] = ["d01-x01-y01","d02-x01-y01"]
analyses["DD"]["BESIII_2017_I1596897"     ] = ["d01-x01-y01"]
analyses["DD"]["BESIII_2018_I1653121"     ] = ["d01-x01-y01","d01-x01-y02"]
analyses["DD"]["BESIII_2020_I1762922"     ] = ["d01-x01-y01"]
analyses["DD"]["BESIII_2018_I1633425"     ] = ["d01-x01-y01"]
analyses["DD"]["BESIII_2018_I1685535"     ] = ["d01-x01-y01","d02-x01-y01"]
analyses["DD"]["BELLE_2011_I878228"       ] = ["d01-x01-y01","d01-x01-y02","d01-x01-y03"]
analyses["DD"]["BABAR_2010_I864027"       ] = ["d01-x01-y01","d01-x01-y02","d01-x01-y03"]
analyses["DD"]["BABAR_2009_I815035"       ] = ["d01-x01-y01","d01-x01-y02","d01-x01-y03","d02-x01-y01"]
analyses["DD"]["BES_1999_I508349"         ] = ["d01-x01-y01","d01-x01-y02","d01-x01-y03","d01-x01-y04"]
# BB
analyses["BB"]["BELLE_2016_I1389855" ] = ["d01-x01-y02","d01-x01-y03"]
analyses["BB"]["BELLE_2008_I764099"  ] = ["d01-x01-y01","d02-x01-y01",
                                          "d03-x01-y01","d04-x01-y01"]
analyses["BB"]["CLEO_1999_I474676"   ] = ["d01-x01-y01","d01-x01-y02"]
analyses["BB"]["CUSB_1991_I325661"   ] = ["d01-x01-y01"]
analyses["BB"]["CLEO_1991_I29927"    ] = ["d01-x01-y01"]
# hyperons
analyses["LL"]["BESIII_2018_I1627871"] = ["d01-x01-y01"]
analyses["LL"]["DM2_1990_I297706"    ] = ["d02-x01-y01"]
analyses["LL"]["BESIII_2019_I1758883"] = ["d01-x01-y05"]
analyses["LL"]["BESIII_2019_I1726357"] = ["d01-x01-y01"]
analyses["LL"]["BABAR_2007_I760730"  ] = ["d01-x01-y01","d02-x01-y01","d03-x01-y01"]
# list the analysis if required and quit()
allProcesses=False
if "All" in opts.processes :
    allProcesses=True
    processes = sorted(list(analyses.keys()))
else :
    processes = sorted(list(set(opts.processes)))
if(opts.list) :
    for process in processes :
        print (" ".join(analyses[process]))
    quit()
if(opts.plot) :
    output=""
    for process in processes:
        for analysis in analyses[process] :
            if(analysis=="CMD3_2019_I1770428") :
                for iy in range(1,3) :
                    output+= " -m/%s/%s" % (analysis,"d02-x01-y0%s"%iy)
            elif(analysis=="BES_1999_I508349") :
                for ix in range(2,4) :
                    for iy in range(1,3) :
                        output+= " -m/%s/%s" % (analysis,"d0%s-x01-y0%s"%(ix,iy))
            elif(analysis=="BESIII_2019_I1726357") :
                for ix in range(2,4) :
                    output+= " -m/%s/%s" % (analysis,"d0%s-x01-y01"% ix)
            for plot in analyses[process][analysis]:
                output+= " -m/%s/%s" % (analysis,plot)
    print (output)
    quit()
# mapping of process to me to use
me = { "PiPi"         : "MEee2Pions",
       "KK"           : "MEee2Kaons",
       "3Pi"          : "MEee3Pions",
       "4Pi"          : "MEee4Pions",
       "EtaPiPi"      : "MEee2EtaPiPi",
       "EtaprimePiPi" : "MEee2EtaPrimePiPi",
       "EtaPhi"       : "MEee2EtaPhi",
       "EtaOmega"     : "MEee2EtaOmega",
       "OmegaPi"      : "MEee2OmegaPi",
       "OmegaPiPi"    : "MEee2OmegaPiPi",
       "PhiPi"        : "MEee2PhiPi",
       "PiGamma"      : "MEee2PiGamma",
       "EtaGamma"     : "MEee2EtaGamma",
       "PPbar"        : "MEee2PPbar",
       "LL"           : "MEee2LL"   ,
       "2K1Pi"        : "MEee2KKPi" }

# get the particle masses from Herwig
particles = { "pi+" : 0., "pi0" : 0. ,"eta" : 0. ,"eta'" : 0. ,"phi" : 0. ,"omega" : 0. ,"p+" : 0. ,"K+" : 0.}
for val in particles :
    tempTxt = "get /Herwig/Particles/%s:NominalMass\nget /Herwig/Particles/%s:WidthLoCut\n" % (val,val)
    with open("temp.in",'w') as f:
        f.write(tempTxt)
    p = subprocess.Popen(["../src/Herwig", "read","temp.in"],stdout=subprocess.PIPE)
    vals = p.communicate()[0].split()
    mass = float(vals[0])-float(vals[1])
    particles[val]=mass
    os.remove("temp.in")
# minimum CMS energies for specific processes
minMass = { "PiPi"         : 2.*particles["pi+"],
            "KK"           : 2.*particles["K+"],
            "3Pi"          : 2.*particles["pi+"]+particles["pi0"],
            "4Pi"          : 2.*particles["pi+"]+2.*particles["pi0"],
            "EtaPiPi"      : particles["eta"]+2.*particles["pi+"],
            "EtaprimePiPi" : particles["eta'"]+2.*particles["pi+"],
            "EtaPhi"       : particles["phi"]+particles["eta"],
            "EtaOmega"     : particles["omega"]+particles["eta"],
            "OmegaPi"      : particles["omega"]+particles["pi0"],
            "OmegaPiPi"    : particles["omega"]+2.*particles["pi0"],
            "PhiPi"        : particles["phi"]+particles["pi0"],
            "PiGamma"      : particles["pi0"],
            "EtaGamma"     : particles["eta"],
            "PPbar"        : 2.*particles["p+"],
            "LL"           : 0.,
            "2K1Pi"        : 2.*particles["K+"]+particles["pi0"] }
# energies we need
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

for process in processes:
    if(process not in analyses) : continue
    matrix=""
    if( process in me ) :
        matrix = me[process]
    for analysis in analyses[process] :
        aos=yoda.read(os.path.join(os.path.join(os.getcwd(),path),analysis+".yoda"))
        if(len(aos)==0) : continue
        for plot in analyses[process][analysis] :
            histo = aos["/REF/%s/%s" %(analysis,plot)]
            for point in histo.points() :
                energy = point.x()
                if(analysis=="KLOE_2009_I797438" or
                   analysis=="KLOE_2005_I655225" or
                   analysis=="KLOE2_2017_I1634981" or
                   analysis=="FENICE_1994_I377833") :
                    energy = math.sqrt(energy)
                if(energy>200) :
                    energy *= 0.001
                emin,delta,anals = nearestEnergy(energy)
                if(energy in energies) :
                    if(analysis not in energies[energy][1]) :
                        energies[energy][1].append(analysis)
                    if(matrix!="" and matrix not in energies[energy][0] and
                       minMass[process]<=energy) :
                        energies[energy][0].append(matrix)
                elif(delta<1e-7) :
                    if(analysis not in anals[1]) :
                        anals[1].append(analysis)
                    if(matrix!="" and matrix not in anals[0] and
                       minMass[process]<=energy) :
                        anals[0].append(matrix)
                else :
                    if(matrix=="") :
                        energies[energy]=[[],[analysis]]
                    elif(minMass[process]<=energy) :
                        energies[energy]=[[matrix],[analysis]]

with open("python/LowEnergy-EE-Perturbative.in", 'r') as f:
    templateText = f.read()
perturbative=Template( templateText )
with open("python/LowEnergy-EE-NonPerturbative.in", 'r') as f:
    templateText = f.read()
nonPerturbative=Template( templateText )

targets=""
for energy in sorted(energies) :
    anal=""
    for analysis in energies[energy][1]: 
        anal+= "insert /Herwig/Analysis/Rivet:Analyses 0 %s\n" %analysis
    proc=""
    matrices = energies[energy][0]
    if(allProcesses) : matrices = me.values()
    for matrix in  matrices:
        proc+="insert SubProcess:MatrixElements 0 %s\n" % matrix
        proc+="set %s:Flavour %s\n" % (matrix,opts.flavour)
    maxflavour =5
    if(energy<thresholds[1]) :
        maxflavour=2
    elif(energy<thresholds[2]) :
        maxflavour=3
    elif(energy<thresholds[3]) :
        maxflavour=4
    # input file for perturbative QCD
    if(opts.perturbative and energy> thresholds[0]) :
        inputPerturbative = perturbative.substitute({"ECMS" : "%8.6f" % energy, "ANALYSES" : anal,
                                                     "lepton" : "", "maxflavour" : maxflavour})
        with open(opts.dest+"/Rivet-LowEnergy-EE-Perturbative-%8.6f.in" % energy ,'w') as f:
            f.write(inputPerturbative)
        targets += "Rivet-LowEnergy-EE-Perturbative-%8.6f.yoda " % energy
    # input file for currents
    if(opts.nonPerturbative and proc!="") :
        inputNonPerturbative = nonPerturbative.substitute({"ECMS" : "%8.6f" % energy, "ANALYSES" : anal,
                                                           "processes" : proc})
        with open(opts.dest+"/Rivet-LowEnergy-EE-NonPerturbative-%8.6f.in" % energy ,'w') as f:
            f.write(inputNonPerturbative)
        targets += "Rivet-LowEnergy-EE-NonPerturbative-%8.6f.yoda " % energy
print (targets)
