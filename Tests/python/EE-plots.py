#! /usr/bin/env python
import glob,os
directory="Rivet-EE"

header="""<html>
<head>
<title>{title}</title>
<style>
      html {{ font-family: sans-serif; }}
      img {{ border: 0; }}
      a {{ text-decoration: none; font-weight: bold; }}
    </style>
            <script type="text/x-mathjax-config">
        MathJax.Hub.Config({{
          tex2jax: {{inlineMath: [["$","$"]]}}
        }});
        </script>
        <script type="text/javascript"
          src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
        </script>
        </head>
<body><center><h1>{title}</h1></center>"""


analyses={ "HadronDecays"     : { },
           "TauDecays" : {},
           "Charged" : {"TotalChargedMult" : { 0 : {}, 1 : {}, 4 : {}, 5 : {}, 51 : {}, 41 : {} , "C" : {} },
                        "ChargedSpectrum" : { 0 : {}, 1 : {}, 2 : {}, 4 : {}, 5 : {} },
                        "ChargedRapidityThrust"    : { },
                        "ChargedpTInThrust"        : { },
                        "ChargedpTOutThrust"       : { },
                        "ChargedRapiditySphericity": { },
                        "ChargedpTInSphericity"    : { },
                        "ChargedpTOutSphericity"   : { },
                        "ChargedpLSphericity"      : { },
                        "ChargedpTSphericity"      : { },
                        "ChargedpT2Sphericity"     : { },
                        "DistChargedMult"  : { 0 : {}, 1 : {}, 2 : {}, 4 : {}, 5: {}, 21 : {}, "C" : {} }},
           "IdentifiedParticle"  : { 22      : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     111     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     211     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     221     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     331     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     223     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     333     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     321     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     311     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     313     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     323     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     2212    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     413     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     423     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     3122    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     3212    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     2224    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     3312    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     3222    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     "3224B" : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     431     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     433     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     3112    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     3224    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     3114    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     3324    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     3124    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     443     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     9010221 : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     9000211 : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     225     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     335     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     113     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     213     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     421     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     411     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     425     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     511     : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     4122    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     3334    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     4332    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     4132    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     4112    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     4114    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     4124    : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}},
                                     14122   : { "x" : {}, "p" : {}, "xi" : {}, "Other" : {}, "Ratio" : {}}},
           "IdentifiedParticleFlavour" : {111  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} }, 211  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} },
                                          321  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} }, 311  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} },
                                          313  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} }, 333  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} },
                                          2212 : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} }, 3122 : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} },
                                          413  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} }}, 
           "MultiplicityFlavour" : {111  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} },  211   : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} },
                                    321  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} },  2212  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} },
                                    413  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} },  3122  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} },
                                    313  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} },  333   : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} },
                                    311  : { 1 : {}, 4 : {}, 5 : {}, 41 : {}, 51 : {} },
           }, 
           "Multiplicity"  : { 22: {}, 111 : {}, 211 : {}, 221 : {}, 331 : {}, 113 : {},9010221 : {},9000211 : {}, 213 : {},
                               223 : {}, 333 : {}, 321  : {}, 311 : {}, 411 : {}, 421 : {}, 431 : {}, 521 : {}, 531 : {}, 511 : {},
                               313 : {}, 323 : {}, 2212 : {}, 413 : {}, 423 : {}, 433 : {}, 513 : {}, 515 : {},
                               3122  : {}, 3312 : {}, 3212 : {}, 3112 : {}, 3114 : {}, 3324: {}, 3334 : {},
                               443 : {}, 100443 : {}, 553 : {}, 20223 : {}, 20333 : {}, 20443 : {}, 225 :{}, 335 : {}, 20431 : {}, 435 : {}, 315 : {}, 325 : {},
                               3222 : {}, 2224 : {},  3224 : {}, 3114 : {}, 4122 : {}, 5122 :{}, 3124 :{}, 4222 : {}, "3222B" : {}, "3224B" : {},
           },
           "EventShapes" : { "T" : {}, "S" : {}, "D" : {}, "O" : {}, "Minor" : {}, "Major" : {},
                             "y12_dur"  : {}, "y23_dur"  : {}, "y34_dur"  : {}, "y45_dur"  : {}, "y56_dur"  : {},
                             "y12_jade" : {}, "y23_jade" : {}, "y34_jade" : {}, "y45_jade" : {}, "y56_jade" : {},
                             "HeavyJetMass" : {} , "JetMassDifference" : {} ,
                             "LightJetMass" : {}, "TotalJetMass" : {} , "EEC" : {}, "AEEC" : {} ,
                             "P" : {}, "A" : {} , "BW" : {}, "BT" : {}, "BN" : {},
                             "Bdiff" : {}, "C" : {},
                             "1jet_dur"  : {}, "2jet_dur"  : {}, "3jet_dur"  : {}, "4jet_dur"  : {}, "5jet_dur"  : {}, "6jet_dur"  : {},
                             "1jet_jade" : {}, "2jet_jade" : {}, "3jet_jade" : {}, "4jet_jade" : {}, "5jet_jade" : {}, "6jet_jade" : {},
                             "Moment_T":{}, "Moment_H":{},"Moment_C":{},"Moment_S":{},"Moment_L":{},"Moment_y":{},"Moment_BW":{},
                             "Moment_BN":{},"Moment_BT":{},"Moment_O":{},"Moment_M":{},"Moment_m":{},},
           "EventShapesFlavour" : { "T" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},
                                    "S" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},
                                    "D" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},
                                    "O" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},
                                    "Minor" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},
                                    "Major" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},
                             "y12" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},"y23" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},"y34" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},"y45" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},"y56" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}}, "HeavyJetMass" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}} , "JetMassDifference" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}} ,
                             "LightJetMass" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}}, "TotalJetMass" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}} , "EEC" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}}, "AEEC" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}} ,
                             "P" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}}, "A" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}} , "BW" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}}, "BT" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}}, "BN" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},
                             "Bdiff" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}}, "C" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},
                                    "1jet" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},"2jet" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},"3jet" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},"4jet" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},"5jet" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},"6jet" : { 1 : {}, 2 : {}, 4 : {}, 5 : {}},},

           "QED" : {},
}
# hadron decays
analyses["HadronDecays"][223] = ["/A2_2017_I1486671/d01-x01-y01","/MC_OmegaPhia1_3Pion_Decay/dalitz_1",
                                 "/MC_OmegaPhia1_3Pion_Decay/m0_1","/MC_OmegaPhia1_3Pion_Decay/mminus_1","/MC_OmegaPhia1_3Pion_Decay/mplus_1",
                                 "/MC_OmegaPhia1_3Pion_Decay/xhist_1","/MC_OmegaPhia1_3Pion_Decay/yhist_1",]
analyses["HadronDecays"][221] = ["/A2_2017_I1486671/d02-x01-y01","/KLOE2_2016_I1416990/d01-x01-y01","/KLOE2_2016_I1416990/d02-x01-y01",
                                 "/KLOE2_2016_I1416990/d02-x01-y02","/KLOE2_2016_I1416990/d02-x01-y03","/KLOE2_2016_I1416990/d02-x01-y04",
                                 "/KLOE2_2016_I1416990/d02-x01-y05","/KLOE2_2016_I1416990/d02-x01-y06","/KLOE2_2016_I1416990/d02-x01-y07",
                                 "/KLOE2_2016_I1416990/d02-x01-y08","/KLOE2_2016_I1416990/d02-x01-y09","/KLOE2_2016_I1416990/d02-x01-y10",
                                 "/KLOE2_2016_I1416990/d02-x01-y11","/KLOE2_2016_I1416990/d02-x01-y12","/KLOE2_2016_I1416990/d02-x01-y13",
                                 "/KLOE2_2016_I1416990/d02-x01-y14","/KLOE2_2016_I1416990/d02-x01-y15","/KLOE2_2016_I1416990/d02-x01-y16",
                                 "/KLOE2_2016_I1416990/d02-x01-y17","/MC_Eta_Decay/dpi0pi0_0","/MC_Eta_Decay/dpi0pip_0",
                                 "/MC_Eta_Decay/dpippim_0","/MC_Eta_Decay/mgammagamma_0","/MC_Eta_Decay/mpi0gamma_0",
                                 "/MC_Eta_Decay/mpimgamma_0","/MC_Eta_Decay/mpipgamma_0","/MC_Eta_Decay/mpippim_0",
                                 "/MC_Eta_Decay/photonenergy_0",]

analyses["HadronDecays"][331] = ["/BESIII_2015_I1364494/d01-x01-y03","/BESIII_2018_I1641075/d01-x01-y05",
                                 "/MC_Eta_Decay/dpi0eta","/MC_Eta_Decay/dpi0pi0_1","/MC_Eta_Decay/dpi0pi0_2",
                                 "/MC_Eta_Decay/dpi0pim_0","/MC_Eta_Decay/dpi0pim_1","/MC_Eta_Decay/dpi0pip_1",
                                 "/MC_Eta_Decay/dpimeta","/MC_Eta_Decay/dpipeta","/MC_Eta_Decay/dpippim_1",
                                "/MC_Eta_Decay/dpippim_2","/MC_Eta_Decay/mgammagamma_1","/MC_Eta_Decay/mpi0gamma_1",
                                 "/MC_Eta_Decay/mpimgamma_1","/MC_Eta_Decay/mpipgamma_1","/MC_Eta_Decay/mpippim_1",
                                 "/MC_Eta_Decay/photonenergy_1",]
analyses["HadronDecays"][333] = ["/KLOE2_2014_I1317236/d01-x01-y01","/KLOE2_2016_I1416825/d01-x01-y01",
                                 "/KLOE_2002_I585183/d01-x01-y01","/KLOE_2009_I818106/d01-x01-y01",
                                 "/SND_2000_I525398/d01-x01-y01","/SND_2000_I527094/d01-x01-y01",
                                 "/SND_2001_I558279/d01-x01-y01","/SND_2001_I558279/d02-x01-y01",
                                 "/MC_OmegaPhia1_3Pion_Decay/dalitz_2",
                                 "/MC_OmegaPhia1_3Pion_Decay/m0_2","/MC_OmegaPhia1_3Pion_Decay/mminus_2","/MC_OmegaPhia1_3Pion_Decay/mplus_2",
                                 "/MC_OmegaPhia1_3Pion_Decay/xhist_2","/MC_OmegaPhia1_3Pion_Decay/yhist_2",]
analyses["HadronDecays"][411] = ["/BESIII_2017_I1519425/d01-x01-y01","/BESIII_2017_I1519425/d02-x01-y01",
                                 "/MC_D_Dalitz/dalitz3","/MC_D_Dalitz/dalitz4","/MC_D_Dalitz/dalitz5",
                                 "/MC_D_Dalitz/h_Kpi04","/MC_D_Dalitz/h_Kpiall3","/MC_D_Dalitz/h_Kpihigh3",
                                 "/MC_D_Dalitz/h_Kpilow3","/MC_D_Dalitz/h_Kpip4","/MC_D_Dalitz/h_kppim5",
                                 "/MC_D_Dalitz/h_kppip5","/MC_D_Dalitz/h_pipi3","/MC_D_Dalitz/h_pipi4","/MC_D_Dalitz/h_pippim5",
                                 "/MC_Semi_Leptonic_Decay/h_411m_10313p_13p_energy","/MC_Semi_Leptonic_Decay/h_411m_10313p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411m_111p_13p_energy","/MC_Semi_Leptonic_Decay/h_411m_111p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411m_113p_13p_energy","/MC_Semi_Leptonic_Decay/h_411m_113p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411m_221p_13p_energy","/MC_Semi_Leptonic_Decay/h_411m_221p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411m_223p_13p_energy","/MC_Semi_Leptonic_Decay/h_411m_223p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411m_311p_13p_energy","/MC_Semi_Leptonic_Decay/h_411m_311p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411m_313p_13p_energy","/MC_Semi_Leptonic_Decay/h_411m_313p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411m_315p_13p_energy","/MC_Semi_Leptonic_Decay/h_411m_315p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411m_331p_13p_energy","/MC_Semi_Leptonic_Decay/h_411m_331p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411p_10313m_11m_energy","/MC_Semi_Leptonic_Decay/h_411p_10313m_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411p_111p_11m_energy","/MC_Semi_Leptonic_Decay/h_411p_111p_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411p_113p_11m_energy","/MC_Semi_Leptonic_Decay/h_411p_113p_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411p_221p_11m_energy","/MC_Semi_Leptonic_Decay/h_411p_221p_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411p_223p_11m_energy","/MC_Semi_Leptonic_Decay/h_411p_223p_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411p_311m_11m_energy","/MC_Semi_Leptonic_Decay/h_411p_311m_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411p_313m_11m_energy","/MC_Semi_Leptonic_Decay/h_411p_313m_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411p_315m_11m_energy","/MC_Semi_Leptonic_Decay/h_411p_315m_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_411p_331p_11m_energy","/MC_Semi_Leptonic_Decay/h_411p_331p_11m_scale"]
analyses["HadronDecays"][421] = ["/BESIII_2015_I1391138/d01-x01-y03","/BESIII_2015_I1391138/d02-x01-y03",
                                 "/MC_D_Dalitz/dalitz1","/MC_D_Dalitz/dalitz2","/MC_D_Dalitz/h_minus1","/MC_D_Dalitz/h_minus2",
                                 "/MC_D_Dalitz/h_neutral2","/MC_D_Dalitz/h_pipi1","/MC_D_Dalitz/h_pipi2","/MC_D_Dalitz/h_plus1",
                                 "/MC_Semi_Leptonic_Decay/h_421m_10323p_13p_energy","/MC_Semi_Leptonic_Decay/h_421m_10323p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_421m_211p_13p_energy","/MC_Semi_Leptonic_Decay/h_421m_211p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_421m_213p_13p_energy","/MC_Semi_Leptonic_Decay/h_421m_213p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_421m_321p_13p_energy","/MC_Semi_Leptonic_Decay/h_421m_321p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_421m_323p_13p_energy","/MC_Semi_Leptonic_Decay/h_421m_323p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_421m_325p_13p_energy","/MC_Semi_Leptonic_Decay/h_421m_325p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_421p_10323m_11m_energy","/MC_Semi_Leptonic_Decay/h_421p_10323m_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_421p_211m_11m_energy","/MC_Semi_Leptonic_Decay/h_421p_211m_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_421p_213m_11m_energy","/MC_Semi_Leptonic_Decay/h_421p_213m_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_421p_321m_11m_energy","/MC_Semi_Leptonic_Decay/h_421p_321m_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_421p_323m_11m_energy","/MC_Semi_Leptonic_Decay/h_421p_323m_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_421p_325m_11m_energy","/MC_Semi_Leptonic_Decay/h_421p_325m_11m_scale",
                                 "/BABAR_2015_I1334693/d01-x01-y01"]
analyses["HadronDecays"][431] = ["/MC_D_Dalitz/dalitz6","/MC_D_Dalitz/h_kppim6","/MC_D_Dalitz/h_kppip6","/MC_D_Dalitz/h_pippim6",
                                 "/MC_Semi_Leptonic_Decay/h_431m_221p_13p_energy","/MC_Semi_Leptonic_Decay/h_431m_221p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_431m_311p_13p_energy","/MC_Semi_Leptonic_Decay/h_431m_311p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_431m_331p_13p_energy","/MC_Semi_Leptonic_Decay/h_431m_331p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_431m_333p_13p_energy","/MC_Semi_Leptonic_Decay/h_431m_333p_13p_scale",
                                 "/MC_Semi_Leptonic_Decay/h_431p_221p_11m_energy","/MC_Semi_Leptonic_Decay/h_431p_221p_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_431p_311m_11m_energy","/MC_Semi_Leptonic_Decay/h_431p_311m_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_431p_311p_11m_energy","/MC_Semi_Leptonic_Decay/h_431p_311p_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_431p_313p_11m_energy","/MC_Semi_Leptonic_Decay/h_431p_313p_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_431p_331p_11m_energy","/MC_Semi_Leptonic_Decay/h_431p_331p_11m_scale",
                                 "/MC_Semi_Leptonic_Decay/h_431p_333p_11m_energy","/MC_Semi_Leptonic_Decay/h_431p_333p_11m_scale"]

analyses["HadronDecays"][511]=["/BELLE_2011_I878990/d01-x01-y01","/BELLE_2013_I1238273/d02-x01-y01", 
                               "/BELLE_2013_I1238273/d04-x01-y01","/BELLE_2013_I1238273/d05-x01-y01",
                               "/BELLE_2015_I1330289/d01-x01-y02",
                               "/BELLE_2015_I1397632/d01-x01-y01","/BELLE_2015_I1397632/d01-x01-y02",
                               "/BELLE_2015_I1397632/d02-x01-y01","/BELLE_2015_I1397632/d02-x01-y02","/ARGUS_1994_I354224/d01-x01-y01",
                               "/MC_Semi_Leptonic_Decay/h_511m_10411p_13p_energy","/MC_Semi_Leptonic_Decay/h_511m_10411p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_511m_10413p_13p_energy","/MC_Semi_Leptonic_Decay/h_511m_10413p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_511m_20413p_13p_energy","/MC_Semi_Leptonic_Decay/h_511m_20413p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_511m_213p_13p_energy","/MC_Semi_Leptonic_Decay/h_511m_213p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_511m_411p_13p_energy","/MC_Semi_Leptonic_Decay/h_511m_411p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_511m_413p_13p_energy","/MC_Semi_Leptonic_Decay/h_511m_413p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_511m_415p_13p_energy","/MC_Semi_Leptonic_Decay/h_511m_415p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_511p_10411m_11m_energy","/MC_Semi_Leptonic_Decay/h_511p_10411m_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_511p_10413m_11m_energy","/MC_Semi_Leptonic_Decay/h_511p_10413m_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_511p_20413m_11m_energy","/MC_Semi_Leptonic_Decay/h_511p_20413m_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_511p_211m_11m_energy","/MC_Semi_Leptonic_Decay/h_511p_211m_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_511p_213m_11m_energy","/MC_Semi_Leptonic_Decay/h_511p_213m_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_511p_411m_11m_energy","/MC_Semi_Leptonic_Decay/h_511p_411m_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_511p_413m_11m_energy","/MC_Semi_Leptonic_Decay/h_511p_413m_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_511p_415m_11m_energy","/MC_Semi_Leptonic_Decay/h_511p_415m_11m_scale"]


analyses["HadronDecays"][521]=["/BELLE_2013_I1238273/d01-x01-y01","/BELLE_2013_I1238273/d03-x01-y01",
                               "/BELLE_2017_I1512299/d01-x01-y01","/BELLE_2017_I1512299/d02-x01-y01",
                               "/BELLE_2017_I1512299/d03-x01-y01","/BELLE_2017_I1512299/d04-x01-y01",
                               "/MC_Semi_Leptonic_Decay/h_521m_10421p_13p_energy","/MC_Semi_Leptonic_Decay/h_521m_10421p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_521m_10423p_13p_energy","/MC_Semi_Leptonic_Decay/h_521m_10423p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_521m_20423p_13p_energy","/MC_Semi_Leptonic_Decay/h_521m_20423p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_521m_331p_13p_energy","/MC_Semi_Leptonic_Decay/h_521m_331p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_521m_421p_13p_energy","/MC_Semi_Leptonic_Decay/h_521m_421p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_521m_423p_13p_energy","/MC_Semi_Leptonic_Decay/h_521m_423p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_521m_425p_13p_energy","/MC_Semi_Leptonic_Decay/h_521m_425p_13p_scale",
                               "/MC_Semi_Leptonic_Decay/h_521p_10421m_11m_energy","/MC_Semi_Leptonic_Decay/h_521p_10421m_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_521p_10423m_11m_energy","/MC_Semi_Leptonic_Decay/h_521p_10423m_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_521p_113p_11m_energy","/MC_Semi_Leptonic_Decay/h_521p_113p_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_521p_20423m_11m_energy","/MC_Semi_Leptonic_Decay/h_521p_20423m_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_521p_223p_11m_energy","/MC_Semi_Leptonic_Decay/h_521p_223p_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_521p_331p_11m_energy","/MC_Semi_Leptonic_Decay/h_521p_331p_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_521p_421m_11m_energy","/MC_Semi_Leptonic_Decay/h_521p_421m_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_521p_423m_11m_energy","/MC_Semi_Leptonic_Decay/h_521p_423m_11m_scale",
                               "/MC_Semi_Leptonic_Decay/h_521p_425m_11m_energy","/MC_Semi_Leptonic_Decay/h_521p_425m_11m_scale",
                               "/BABAR_2013_I1116411/d01-x01-y01"]
analyses["HadronDecays"][443] = ["/BESIII_2018_I1697377/d01-x01-y01",
                                 "/MC_Meson_Meson_Leptons_Decay/h2_443p_22p_11_mVf",
                                 "/MC_Meson_Meson_Leptons_Decay/h2_443p_22p_11_mVfbar",
                                 "/MC_Meson_Meson_Leptons_Decay/h2_443p_22p_11_mff",]
analyses["HadronDecays"][441] = ["/BESIII_2019_I1724880/d01-x01-y01"]
analyses["HadronDecays"][553] = ["/ARGUS_1988_I251097/d01-x01-y01","/ARGUS_1988_I251097/d01-x01-y02",
                                 "/ARGUS_1988_I251097/d01-x01-y03","/ARGUS_1988_I251097/d01-x01-y04",
                                 "/ARGUS_1988_I251097/d01-x01-y05","/ARGUS_1988_I251097/d01-x01-y06",
                                 "/ARGUS_1988_I251097/d01-x01-y07","/ARGUS_1988_I251097/d03-x01-y01",
                                 "/ARGUS_1988_I251097/d07-x01-y01","/ARGUS_1989_I262415/d03-x01-y01",
                                 "/ARGUS_1989_I262551/d02-x01-y01","/ARGUS_1989_I262551/d03-x01-y01",
                                 "/ARGUS_1989_I276860/d05-x01-y01","/ARGUS_1989_I276860/d06-x01-y01",
                                 "/ARGUS_1989_I276860/d07-x01-y01","/ARGUS_1989_I276860/d08-x01-y01",
                                 "/ARGUS_1989_I276860/d09-x01-y01","/ARGUS_1989_I276860/d10-x01-y01",
                                 "/ARGUS_1989_I276860/d11-x01-y01","/ARGUS_1989_I276860/d12-x01-y01",
                                 "/ARGUS_1990_I278933/d01-x01-y02","/ARGUS_1990_I278933/d02-x01-y02",
                                 "/ARGUS_1990_I278933/d04-x01-y01","/ARGUS_1990_I278933/d06-x01-y01",
                                 "/ARGUS_1993_S2669951/d03-x01-y01","/ARGUS_1993_S2789213/d02-x01-y01",
                                 "/ARGUS_1993_S2789213/d02-x01-y02","/ARGUS_1993_S2789213/d02-x01-y03",
                                 "/ARGUS_1993_S2789213/d02-x01-y04","/ARGUS_1993_S2789213/d02-x01-y05",
                                 "/ARGUS_1993_S2789213/d05-x01-y01","/ARGUS_1993_S2789213/d08-x01-y01",
                                 "/ARGUS_1993_S2789213/d11-x01-y01","/ARGUS_1993_S2789213/d14-x01-y01",
                                 "/PLUTO_1981_I165122/d06-x01-y01",
                                 "/ARGUS_1989_I276860/d01-x01-y01","/ARGUS_1989_I276860/d01-x01-y02",
                                 "/ARGUS_1989_I276860/d02-x01-y01","/ARGUS_1989_I276860/d03-x01-y01",
                                 "/ARGUS_1989_I276860/d04-x01-y01","/ARGUS_1989_I276860/d04-x01-y02"]
analyses["HadronDecays"][100553] = ["/ARGUS_1988_I251097/d04-x01-y01","/ARGUS_1988_I251097/d08-x01-y01",
                                    "/ARGUS_1989_I262551/d02-x01-y02","/ARGUS_1989_I262551/d04-x01-y01",
                                    "/ARGUS_1990_I278933/d01-x01-y03","/ARGUS_1990_I278933/d02-x01-y03",
                                    "/ARGUS_1990_I278933/d04-x01-y02","/ARGUS_1990_I278933/d06-x01-y02",
                                    "/ARGUS_1993_S2669951/d04-x01-y01",
                                    "/MC_Onium_PiPi_Decay/h_100443_443_hpi0pi0","/MC_Onium_PiPi_Decay/h_100443_443_hpippim",
                                    "/MC_Onium_PiPi_Decay/h_100443_443_mpi0pi0","/MC_Onium_PiPi_Decay/h_100443_443_mpippim",]
analyses["HadronDecays"][300553] = ["/ARGUS_1993_S2789213/d03-x01-y01","/ARGUS_1993_S2789213/d03-x01-y02",
                                    "/ARGUS_1993_S2789213/d03-x01-y03","/ARGUS_1993_S2789213/d03-x01-y04",
                                    "/ARGUS_1993_S2789213/d03-x01-y05","/ARGUS_1993_S2789213/d06-x01-y01",
                                    "/ARGUS_1993_S2789213/d09-x01-y01","/ARGUS_1992_I319102/d03-x01-y01",
                                    "/ARGUS_1993_S2789213/d12-x01-y01",
                                    "/ARGUS_1993_S2653028/d01-x01-y01","/ARGUS_1993_S2653028/d02-x01-y01",
                                    "/ARGUS_1993_S2653028/d03-x01-y01","/ARGUS_1993_S2653028/d04-x01-y01",
                                    "/ARGUS_1993_S2653028/d05-x01-y01","/ARGUS_1993_S2653028/d06-x01-y01",
                                    "/ARGUS_1993_S2653028/d07-x01-y01","/ARGUS_1993_S2653028/d08-x01-y01",
                                    "/ARGUS_1993_S2653028/d09-x01-y01","/ARGUS_1993_S2653028/d10-x01-y01",
                                    "/ARGUS_1993_S2653028/d11-x01-y01","/BELLE_2001_S4598261/d01-x01-y01",
                                    "/BELLE_2001_S4598261/d02-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d104-x01-y01","/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d106-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d110-x01-y01","/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d113-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d116-x01-y01","/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d29-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d30-x01-y01","/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d31-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d32-x01-y01","/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d33-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d34-x01-y01","/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d48-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d50-x01-y01","/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d51-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d53-x01-y01","/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d60-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d61-x01-y01","/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d62-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d63-x01-y01","/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d64-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d65-x01-y01","/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d87-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d88-x01-y01","/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d89-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d90-x01-y01","/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d92-x01-y01",
                                    "/PDG_Upsilon_4S_HADRON_MULTIPLICITIES/d96-x01-y01",
                                    "/BABAR_2007_S6895344/d03-x01-y01","/BABAR_2007_S6895344/d04-x01-y01",
                                    "/BABAR_2003_I593379/d01-x01-y01","/BABAR_2003_I593379/d01-x01-y02",
                                    "/BABAR_2003_I593379/d01-x01-y03","/BABAR_2003_I593379/d01-x01-y04",
                                    "/BABAR_2003_I593379/d01-x01-y05","/BABAR_2003_I593379/d01-x01-y06",
                                    "/BABAR_2003_I593379/d01-x01-y07","/BABAR_2003_I593379/d06-x01-y01",
                                    "/BABAR_2003_I593379/d07-x01-y01","/BABAR_2003_I593379/d07-x01-y02",
                                    "/BABAR_2003_I593379/d08-x01-y01","/BABAR_2003_I593379/d10-x01-y01"]
analyses["HadronDecays"][20113] = ["/MC_OmegaPhia1_3Pion_Decay/dalitz0","/MC_OmegaPhia1_3Pion_Decay/dalitz2",
                                   "/MC_OmegaPhia1_3Pion_Decay/hist0"  ,"/MC_OmegaPhia1_3Pion_Decay/hist2A",
                                   "/MC_OmegaPhia1_3Pion_Decay/hist2B" ,"/MC_OmegaPhia1_3Pion_Decay/hist2C",]
analyses["HadronDecays"][20213] = ["/MC_OmegaPhia1_3Pion_Decay/dalitz1","/MC_OmegaPhia1_3Pion_Decay/dalitz3",
                                   "/MC_OmegaPhia1_3Pion_Decay/hist1A" ,"/MC_OmegaPhia1_3Pion_Decay/hist1B",
                                   "/MC_OmegaPhia1_3Pion_Decay/hist3A" ,"/MC_OmegaPhia1_3Pion_Decay/hist3B",]
# charged multiplicity (total)
analyses["Charged"]["TotalChargedMult"][0][12.0 ] = ["/JADE_1983_I190818/d01-x01-y01"]
analyses["Charged"]["TotalChargedMult"][0][14.0 ] = ["/TASSO_1989_I277658/d02-x01-y01"]
analyses["Charged"]["TotalChargedMult"][0][21.65] = ["/PLUTO_1980_I154270/d01-x01-y01"]
analyses["Charged"]["TotalChargedMult"][0][22.0 ] = ["/TASSO_1980_I143691/d01-x01-y01"]
analyses["Charged"]["TotalChargedMult"][0][29.0 ] = ["/TPC_1987_I235694/d05-x01-y04","/HRS_1986_I18502/d03-x01-y01"]
analyses["Charged"]["TotalChargedMult"][0][50.0 ] = ["/AMY_1990_I295160/d02-x01-y01"]
analyses["Charged"]["TotalChargedMult"][0][55.7 ] = ["/AMY_1990_I295160/d02-x02-y01"]
analyses["Charged"]["TotalChargedMult"][0][91.2 ] = ["/ALEPH_1991_S2435284/d02-x01-y01","/OPAL_1992_I321190/d05-x01-y01",
                                                     "/DELPHI_1991_I301657/d04-x01-y01","/ALEPH_2004_S5765862/d01-x01-y01",
                                                     "/DELPHI_1996_S3430090/d35-x01-y01","/DELPHI_1998_I473409/d01-x01-y01",
                                                     "/OPAL_1998_S3780481/d09-x01-y04","/ALEPH_1996_S3486095/d19-x01-y01"]
analyses["Charged"]["TotalChargedMult"][5][91.2] = ["/OPAL_2002_S5361494/d01-x01-y01","/OPAL_1998_S3780481/d09-x01-y03",
                                                    "/DELPHI_1998_I473409/d02-x01-y01","/SLD_2004_S5693039/d08-x02-y03",
                                                    "/SLD_1996_S3398250/d01-x01-y01"]
analyses["Charged"]["TotalChargedMult"][4][91.2] = ["/OPAL_2002_S5361494/d01-x01-y02","/OPAL_1998_S3780481/d09-x01-y02",
                                                    "/SLD_2004_S5693039/d08-x02-y02","/SLD_1996_S3398250/d02-x01-y01"]
analyses["Charged"]["TotalChargedMult"][1][91.2] = ["/OPAL_2002_S5361494/d01-x01-y03","/OPAL_1998_S3780481/d09-x01-y01",
                                                    "/DELPHI_1998_I473409/d03-x01-y01","/SLD_2004_S5693039/d08-x02-y01",
                                                    "/SLD_1996_S3398250/d03-x01-y01"]
analyses["Charged"]["TotalChargedMult"][51][91.2] = ["/OPAL_2002_S5361494/d01-x01-y04","/SLD_2004_S5693039/d08-x03-y03",
                                                     "/SLD_1996_S3398250/d05-x01-y01"]
analyses["Charged"]["TotalChargedMult"][41][91.2] = ["/SLD_2004_S5693039/d08-x03-y02","/SLD_1996_S3398250/d04-x01-y01"]

analyses["Charged"]["TotalChargedMult"][1][29.0] = ["/TPC_1987_I235694/d04-x01-y04"]
analyses["Charged"]["TotalChargedMult"][4][29.0] = ["/TPC_1987_I235694/d03-x01-y04"]
analyses["Charged"]["TotalChargedMult"][5][29.0] = ["/TPC_1987_I235694/d02-x01-y04"]
analyses["Charged"]["TotalChargedMult"][1 ][195.0] = ["/DELPHI_2000_S4328825/d01-x01-y03"]
analyses["Charged"]["TotalChargedMult"][4 ][195.0] = ["/DELPHI_2000_S4328825/d01-x01-y02"]
analyses["Charged"]["TotalChargedMult"][5 ][195.0] = ["/DELPHI_2000_S4328825/d01-x01-y01"]
analyses["Charged"]["TotalChargedMult"][51][195.0] = ["/DELPHI_2000_S4328825/d01-x01-y04"]
# charged multiplicity (dist)
analyses["Charged"]["DistChargedMult"][0][10.47] = ["/ARGUS_1992_I319102/d02-x01-y01"]
analyses["Charged"]["DistChargedMult"][0][14.0] = ["/TASSO_1989_I277658/d05-x01-y01"]
analyses["Charged"]["DistChargedMult"][0][22.0] = ["/TASSO_1989_I277658/d05-x01-y02"]
analyses["Charged"]["DistChargedMult"][0][29.0] = ["/HRS_1986_I18502/d01-x01-y01"]
analyses["Charged"]["DistChargedMult"][0][34.8] = ["/TASSO_1989_I277658/d05-x01-y03"]
analyses["Charged"]["DistChargedMult"][0][35.0] = ["/TASSO_1988_I263859/d06-x01-y01"]
analyses["Charged"]["DistChargedMult"][0][43.6] = ["/TASSO_1989_I277658/d05-x01-y04"]
analyses["Charged"]["DistChargedMult"][0][50.0] = ["/AMY_1990_I295160/d01-x01-y01"]
analyses["Charged"]["DistChargedMult"][0][52.0] = ["/AMY_1990_I295160/d01-x01-y02"]
analyses["Charged"]["DistChargedMult"][0][55.0] = ["/AMY_1990_I295160/d01-x01-y03"]
analyses["Charged"]["DistChargedMult"][0][56.0] = ["/AMY_1990_I295160/d01-x01-y04"]
analyses["Charged"]["DistChargedMult"][0][57.0] = ["/AMY_1990_I295160/d01-x01-y05"]
analyses["Charged"]["DistChargedMult"][0][60.0] = ["/AMY_1990_I295160/d01-x01-y06"]
analyses["Charged"]["DistChargedMult"][0][60.8] = ["/AMY_1990_I295160/d01-x01-y07"]
analyses["Charged"]["DistChargedMult"][0][61.4] = ["/AMY_1990_I295160/d01-x01-y08"]
analyses["Charged"]["DistChargedMult"][0][55.7] = ["/AMY_1990_I295160/d01-x01-y09"]
analyses["Charged"]["DistChargedMult"][0][91.2] = ["/ALEPH_1991_S2435284/d01-x01-y01","/OPAL_1992_I321190/d01-x01-y01",
                                                   "/DELPHI_1991_I301657/d02-x01-y01","/L3_2004_I652683/d59-x01-y01",
                                                   "/ALEPH_1996_S3486095/d18-x01-y01"]
analyses["Charged"]["DistChargedMult"][2][91.2] = ["/L3_2004_I652683/d59-x01-y02"]
analyses["Charged"]["DistChargedMult"][5][91.2] = ["/L3_2004_I652683/d59-x01-y03"]

analyses["Charged"]["DistChargedMult"][0][130.1]=["/L3_2004_I652683/d60-x01-y01"]
analyses["Charged"]["DistChargedMult"][0][136.3]=["/L3_2004_I652683/d60-x01-y02"]
analyses["Charged"]["DistChargedMult"][0][161.3]=["/L3_2004_I652683/d60-x01-y03"]
analyses["Charged"]["DistChargedMult"][0][172.3]=["/L3_2004_I652683/d61-x01-y01"]
analyses["Charged"]["DistChargedMult"][0][182.8]=["/L3_2004_I652683/d61-x01-y02"]
analyses["Charged"]["DistChargedMult"][0][188.6]=["/L3_2004_I652683/d61-x01-y03"]
analyses["Charged"]["DistChargedMult"][0][194.4]=["/L3_2004_I652683/d62-x01-y01"]
analyses["Charged"]["DistChargedMult"][0][200.2]=["/L3_2004_I652683/d62-x01-y02"]
analyses["Charged"]["DistChargedMult"][0][206.2]=["/L3_2004_I652683/d62-x01-y03"]


analyses["Charged"]["DistChargedMult"]["C"][91.2]=["/DELPHI_1991_I324035/d01-x01-y01","/DELPHI_1991_I324035/d02-x01-y01",
                                                   "/DELPHI_1991_I324035/d03-x01-y01","/DELPHI_1991_I324035/d04-x01-y01",
                                                   "/DELPHI_1991_I324035/d05-x01-y01","/DELPHI_1991_I324035/d06-x01-y01",
                                                   "/DELPHI_1991_I324035/d07-x01-y01","/DELPHI_1991_I324035/d08-x01-y01",
                                                   "/DELPHI_1991_I324035/d09-x01-y01","/DELPHI_1991_I324035/d10-x01-y01",
                                                   "/DELPHI_1991_I324035/d11-x01-y01","/DELPHI_1991_I324035/d12-x01-y01",
                                                   "/DELPHI_1991_I324035/d13-x01-y01",
                                                   "/DELPHI_1992_I334948/d01-x01-y01","/DELPHI_1992_I334948/d01-x01-y02",
                                                   "/DELPHI_1992_I334948/d01-x01-y03","/DELPHI_1992_I334948/d02-x01-y01",
                                                   "/DELPHI_1992_I334948/d02-x01-y02","/DELPHI_1992_I334948/d02-x01-y03",
                                                   "/DELPHI_1992_I334948/d03-x01-y01","/DELPHI_1992_I334948/d03-x01-y02",
                                                   "/DELPHI_1992_I334948/d03-x01-y03",]

analyses["Charged"]["DistChargedMult"][21][ 5.25] = ["/OPAL_2004_I631361/d01-x01-y01"]
analyses["Charged"]["DistChargedMult"][21][ 5.98] = ["/OPAL_2004_I631361/d01-x01-y02"]
analyses["Charged"]["DistChargedMult"][21][ 6.98] = ["/OPAL_2004_I631361/d01-x01-y03"]
analyses["Charged"]["DistChargedMult"][21][ 8.43] = ["/OPAL_2004_I631361/d02-x01-y01"]
analyses["Charged"]["DistChargedMult"][21][10.92] = ["/OPAL_2004_I631361/d02-x01-y02"]
analyses["Charged"]["DistChargedMult"][21][14.24] = ["/OPAL_2004_I631361/d03-x01-y01"]
analyses["Charged"]["DistChargedMult"][21][17.72] = ["/OPAL_2004_I631361/d03-x01-y02"]

analyses["Charged"]["ChargedSpectrum"][0][2.2  ] = ["/BESII_2004_I622224/d01-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][2.6  ] = ["/BESII_2004_I622224/d02-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][3.0  ] = ["/BESII_2004_I622224/d03-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][3.2  ] = ["/BESII_2004_I622224/d04-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][4.6  ] = ["/BESII_2004_I622224/d05-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][4.8  ] = ["/BESII_2004_I622224/d06-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][12.0 ] = ["/TASSO_1980_I153511/d05-x01-y01","/TASSO_1982_I177174/d02-x01-y01",
                                                    "/TASSO_1982_I177174/d03-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][13.0 ] = ["/TASSO_1980_I143691/d05-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][14.0 ] = ["/TASSO_1982_I177174/d01-x01-y01","/TASSO_1982_I177174/d02-x01-y02","/TASSO_1982_I177174/d03-x01-y02"]
analyses["Charged"]["ChargedSpectrum"][0][19.5 ] = ["/TASSO_1980_I143691/d06-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][22.0 ] = ["/TASSO_1982_I177174/d01-x01-y02","/TASSO_1982_I177174/d02-x01-y03","/TASSO_1982_I177174/d03-x01-y03"]
analyses["Charged"]["ChargedSpectrum"][0][25.0 ] = ["/TASSO_1982_I177174/d02-x01-y04","/TASSO_1982_I177174/d03-x01-y04"]
analyses["Charged"]["ChargedSpectrum"][0][29.0 ] = ["/TPC_1988_I262143/d01-x01-y04","/HRS_1985_I201482/d10-x01-y01",
                                                    "/HRS_1985_I201482/d11-x01-y01","/HRS_1985_I201482/d12-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][30.0 ] = ["/TASSO_1982_I177174/d02-x01-y05","/TASSO_1982_I177174/d03-x01-y05","/TASSO_1980_I143691/d07-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][30.8 ] = ["/TASSO_1980_I153511/d06-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][33.0 ] = ["/TASSO_1982_I177174/d01-x01-y03"]
analyses["Charged"]["ChargedSpectrum"][0][34.0 ] = ["/TASSO_1982_I177174/d02-x01-y06","/TASSO_1982_I177174/d03-x01-y06"]
analyses["Charged"]["ChargedSpectrum"][0][35.0 ] = ["/TASSO_1982_I177174/d02-x01-y07","/TASSO_1982_I177174/d03-x01-y07",
                                                    "/TASSO_1988_I263859/d10-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][55.2 ] = ["/AMY_1990_I283337/d02-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][58.0 ] = ["/TOPAZ_1995_I381900/d01-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][91.2 ] = ["/DELPHI_1998_I473409/d16-x01-y01","/DELPHI_1998_I473409/d17-x01-y01",
                                                    "/ALEPH_1996_S3486095/d09-x01-y01","/OPAL_1998_S3780481/d04-x01-y01",
                                                    "/OPAL_1998_S3780481/d08-x01-y01","/SLD_2004_S5693039/d01-x01-y01",
                                                    "/L3_2004_I652683/d65-x01-y01","/SLD_1999_S3743934/d04-x01-y01",
                                                    "/ALEPH_1996_S3486095/d17-x01-y01","/DELPHI_1996_S3430090/d07-x01-y01",
                                                    "/DELPHI_1996_S3430090/d08-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][133.0] = ["/ALEPH_2004_S5765862/d02-x01-y01","/ALEPH_2004_S5765862/d11-x01-y01",
                                                    "/ALEPH_2004_S5765862/d19-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][161.0] = ["/ALEPH_2004_S5765862/d03-x01-y01","/ALEPH_2004_S5765862/d12-x01-y01",
                                                    "/ALEPH_2004_S5765862/d20-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][172.0] = ["/ALEPH_2004_S5765862/d04-x01-y01","/ALEPH_2004_S5765862/d13-x01-y01",
                                                    "/ALEPH_2004_S5765862/d21-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][183.0] = ["/DELPHI_2003_I620250/d32-x01-y01","/ALEPH_2004_S5765862/d05-x01-y01",
                                                    "/ALEPH_2004_S5765862/d14-x01-y01","/ALEPH_2004_S5765862/d22-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][189.0] = ["/DELPHI_2003_I620250/d32-x01-y02","/ALEPH_2004_S5765862/d06-x01-y01",
                                                    "/ALEPH_2004_S5765862/d15-x01-y01","/ALEPH_2004_S5765862/d23-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][192.0] = ["/DELPHI_2003_I620250/d32-x01-y03"]
analyses["Charged"]["ChargedSpectrum"][0][196.0] = ["/DELPHI_2003_I620250/d32-x01-y04","/ALEPH_2004_S5765862/d07-x01-y01",
                                                    "/ALEPH_2004_S5765862/d16-x01-y01","/ALEPH_2004_S5765862/d24-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][200.0] = ["/DELPHI_2003_I620250/d33-x01-y01","/ALEPH_2004_S5765862/d08-x01-y01",
                                                    "/ALEPH_2004_S5765862/d17-x01-y01","/ALEPH_2004_S5765862/d25-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][200.5] = ["/OPAL_2003_I595335/d04-x01-y01","/OPAL_2003_I595335/d05-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][202.0] = ["/DELPHI_2003_I620250/d33-x01-y02"]
analyses["Charged"]["ChargedSpectrum"][0][205.0] = ["/DELPHI_2003_I620250/d33-x01-y03"]
analyses["Charged"]["ChargedSpectrum"][0][206.0] = ["/ALEPH_2004_S5765862/d09-x01-y01","/ALEPH_2004_S5765862/d18-x01-y01",
                                                    "/ALEPH_2004_S5765862/d26-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][207.0] = ["/DELPHI_2003_I620250/d33-x01-y04"]



analyses["Charged"]["ChargedSpectrum"][0][130.1]=["/L3_2004_I652683/d66-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][136.3]=["/L3_2004_I652683/d66-x01-y02"]
analyses["Charged"]["ChargedSpectrum"][0][161.3]=["/L3_2004_I652683/d66-x01-y03"]
analyses["Charged"]["ChargedSpectrum"][0][172.3]=["/L3_2004_I652683/d67-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][182.8]=["/L3_2004_I652683/d67-x01-y02"]
analyses["Charged"]["ChargedSpectrum"][0][188.6]=["/L3_2004_I652683/d67-x01-y03"]
analyses["Charged"]["ChargedSpectrum"][0][194.4]=["/L3_2004_I652683/d68-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][0][200.2]=["/L3_2004_I652683/d68-x01-y02"]
analyses["Charged"]["ChargedSpectrum"][0][206.2]=["/L3_2004_I652683/d68-x01-y03"]


analyses["Charged"]["ChargedRapidityThrust"][13.0 ] = ["/TASSO_1980_I143691/d02-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][19.5 ] = ["/TASSO_1980_I143691/d03-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][30.0 ] = ["/TASSO_1980_I143691/d04-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][55.2 ] = ["/AMY_1990_I283337/d01-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][91.2 ] = ["/DELPHI_1996_S3430090/d05-x01-y01","/ALEPH_1996_S3486095/d10-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][133.0] = ["/ALEPH_2004_S5765862/d36-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][161.0] = ["/ALEPH_2004_S5765862/d37-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][172.0] = ["/ALEPH_2004_S5765862/d38-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][183.0] = ["/DELPHI_2003_I620250/d30-x01-y01","/ALEPH_2004_S5765862/d39-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][189.0] = ["/DELPHI_2003_I620250/d30-x01-y02","/ALEPH_2004_S5765862/d40-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][192.0] = ["/DELPHI_2003_I620250/d30-x01-y03"]
analyses["Charged"]["ChargedRapidityThrust"][196.0] = ["/DELPHI_2003_I620250/d30-x01-y04","/ALEPH_2004_S5765862/d41-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][200.0] = ["/DELPHI_2003_I620250/d31-x01-y01","/ALEPH_2004_S5765862/d42-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][200.5] = ["/OPAL_2003_I595335/d03-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][202.0] = ["/DELPHI_2003_I620250/d31-x01-y02"]
analyses["Charged"]["ChargedRapidityThrust"][205.0] = ["/DELPHI_2003_I620250/d31-x01-y03"]
analyses["Charged"]["ChargedRapidityThrust"][206.0] = ["/ALEPH_2004_S5765862/d43-x01-y01"]
analyses["Charged"]["ChargedRapidityThrust"][207.0] = ["/DELPHI_2003_I620250/d31-x01-y04"]

analyses["Charged"]["ChargedpTInThrust" ][91.2 ] = ["/DELPHI_1996_S3430090/d01-x01-y01","/ALEPH_1996_S3486095/d11-x01-y01"]
analyses["Charged"]["ChargedpTInThrust" ][133.0] = ["/ALEPH_2004_S5765862/d27-x01-y01"]
analyses["Charged"]["ChargedpTInThrust" ][161.0] = ["/ALEPH_2004_S5765862/d28-x01-y01"]
analyses["Charged"]["ChargedpTInThrust" ][172.0] = ["/ALEPH_2004_S5765862/d29-x01-y01"]
analyses["Charged"]["ChargedpTInThrust" ][183.0] = ["/DELPHI_2003_I620250/d34-x01-y01","/ALEPH_2004_S5765862/d30-x01-y01"]
analyses["Charged"]["ChargedpTInThrust" ][189.0] = ["/DELPHI_2003_I620250/d34-x01-y02","/ALEPH_2004_S5765862/d31-x01-y01"]
analyses["Charged"]["ChargedpTInThrust" ][192.0] = ["/DELPHI_2003_I620250/d34-x01-y03"]
analyses["Charged"]["ChargedpTInThrust" ][196.0] = ["/DELPHI_2003_I620250/d34-x01-y04","/ALEPH_2004_S5765862/d32-x01-y01"]
analyses["Charged"]["ChargedpTInThrust" ][200.0] = ["/DELPHI_2003_I620250/d35-x01-y01","/ALEPH_2004_S5765862/d33-x01-y01"]
analyses["Charged"]["ChargedpTInThrust" ][200.5] = ["/OPAL_2003_I595335/d01-x01-y01"]
analyses["Charged"]["ChargedpTInThrust" ][202.0] = ["/DELPHI_2003_I620250/d35-x01-y02"]
analyses["Charged"]["ChargedpTInThrust" ][205.0] = ["/DELPHI_2003_I620250/d35-x01-y03"]
analyses["Charged"]["ChargedpTInThrust" ][206.0] = ["/ALEPH_2004_S5765862/d34-x01-y01"]
analyses["Charged"]["ChargedpTInThrust" ][207.0] = ["/DELPHI_2003_I620250/d35-x01-y04"]
analyses["Charged"]["ChargedpTInSphericity" ][35.0]=["/TASSO_1988_I263859/d07-x01-y01"]
analyses["Charged"]["ChargedpTInSphericity" ][55.2 ]    = ["/AMY_1990_I283337/d06-x01-y01"]
analyses["Charged"]["ChargedpTInSphericity" ][91.2 ]    = ["/DELPHI_1996_S3430090/d03-x01-y01"]

analyses["Charged"]["ChargedpTOutSphericity"][35.0 ]    = ["/TASSO_1988_I263859/d08-x01-y01"]
analyses["Charged"]["ChargedpTOutSphericity"][55.2 ]    = ["/AMY_1990_I283337/d07-x01-y01"]
analyses["Charged"]["ChargedpTOutSphericity"][91.2 ]    = ["/DELPHI_1996_S3430090/d04-x01-y01"]
analyses["Charged"]["ChargedRapiditySphericity"][35.0 ] = ["/TASSO_1988_I263859/d11-x01-y01"]
analyses["Charged"]["ChargedRapiditySphericity"][91.2 ] = ["/DELPHI_1996_S3430090/d06-x01-y01"]
analyses["Charged"]["ChargedRapiditySphericity"][133.0] = ["/ALEPH_2004_S5765862/d44-x01-y01"]
analyses["Charged"]["ChargedRapiditySphericity"][161.0] = ["/ALEPH_2004_S5765862/d45-x01-y01"]
analyses["Charged"]["ChargedRapiditySphericity"][172.0] = ["/ALEPH_2004_S5765862/d46-x01-y01"]
analyses["Charged"]["ChargedRapiditySphericity"][183.0] = ["/ALEPH_2004_S5765862/d47-x01-y01"]
analyses["Charged"]["ChargedRapiditySphericity"][189.0] = ["/ALEPH_2004_S5765862/d48-x01-y01"]
analyses["Charged"]["ChargedRapiditySphericity"][196.0] = ["/ALEPH_2004_S5765862/d49-x01-y01"]
analyses["Charged"]["ChargedRapiditySphericity"][200.0] = ["/ALEPH_2004_S5765862/d50-x01-y01"]
analyses["Charged"]["ChargedRapiditySphericity"][206.0] = ["/ALEPH_2004_S5765862/d51-x01-y01"]

analyses["Charged"]["ChargedpTOutThrust"][91.2 ] = ["/DELPHI_1996_S3430090/d02-x01-y01","/ALEPH_1996_S3486095/d12-x01-y01"]
analyses["Charged"]["ChargedpTOutThrust"][183.0] = ["/DELPHI_2003_I620250/d36-x01-y01"]
analyses["Charged"]["ChargedpTOutThrust"][189.0] = ["/DELPHI_2003_I620250/d36-x01-y02"]
analyses["Charged"]["ChargedpTOutThrust"][192.0] = ["/DELPHI_2003_I620250/d36-x01-y03"]
analyses["Charged"]["ChargedpTOutThrust"][196.0] = ["/DELPHI_2003_I620250/d36-x01-y04"]
analyses["Charged"]["ChargedpTOutThrust"][200.0] = ["/DELPHI_2003_I620250/d37-x01-y01"]
analyses["Charged"]["ChargedpTOutThrust"][200.5] = ["/OPAL_2003_I595335/d02-x01-y01"]
analyses["Charged"]["ChargedpTOutThrust"][202.0] = ["/DELPHI_2003_I620250/d37-x01-y02"]
analyses["Charged"]["ChargedpTOutThrust"][205.0] = ["/DELPHI_2003_I620250/d37-x01-y03"]
analyses["Charged"]["ChargedpTOutThrust"][206.0] = ["/ALEPH_2004_S5765862/d35-x01-y01"]
analyses["Charged"]["ChargedpTOutThrust"][207.0] = ["/DELPHI_2003_I620250/d37-x01-y04"]
# identified particle (flavour sep)
analyses["IdentifiedParticleFlavour"][111 ][5][91.2]=["/DELPHI_1996_I401100/d03-x01-y01","/SLD_2004_S5693039/d05-x01-y03"]
analyses["IdentifiedParticleFlavour"][211 ][5][91.2]=["/DELPHI_1998_I473409/d26-x01-y01","/DELPHI_1998_I473409/d27-x01-y01",
                                                      "/SLD_1999_S3743934/d10-x01-y03"]
analyses["IdentifiedParticleFlavour"][321 ][5][91.2]=["/DELPHI_1998_I473409/d28-x01-y01","/DELPHI_1998_I473409/d29-x01-y01",
                                                      "/SLD_2004_S5693039/d06-x01-y03","/SLD_1999_S3743934/d12-x01-y03"]
analyses["IdentifiedParticleFlavour"][2212][5][91.2]=["/DELPHI_1998_I473409/d30-x01-y01","/DELPHI_1998_I473409/d31-x01-y01",
                                                      "/SLD_2004_S5693039/d07-x01-y03","/SLD_1999_S3743934/d16-x01-y03"]

analyses["IdentifiedParticleFlavour"][211 ][4][91.2]=["/SLD_2004_S5693039/d05-x01-y02","/SLD_1999_S3743934/d10-x01-y02"]
                                                      
analyses["IdentifiedParticleFlavour"][321 ][4][91.2]=["/SLD_2004_S5693039/d06-x01-y02","/SLD_1999_S3743934/d12-x01-y02"]
analyses["IdentifiedParticleFlavour"][2212][4][91.2]=["/SLD_2004_S5693039/d07-x01-y02","/SLD_1999_S3743934/d16-x01-y02"]

analyses["IdentifiedParticleFlavour"][211 ][1][91.2]=["/DELPHI_1998_I473409/d34-x01-y01","/DELPHI_1998_I473409/d35-x01-y01",
                                                      "/SLD_2004_S5693039/d05-x01-y01","/SLD_1999_S3743934/d10-x01-y01"]
analyses["IdentifiedParticleFlavour"][321 ][1][91.2]=["/DELPHI_1998_I473409/d36-x01-y01","/DELPHI_1998_I473409/d37-x01-y01",
                                                      "/SLD_2004_S5693039/d05-x01-y03","/SLD_2004_S5693039/d06-x01-y01",
                                                      "/SLD_1999_S3743934/d12-x01-y01"]
analyses["IdentifiedParticleFlavour"][2212][1][91.2]=["/DELPHI_1998_I473409/d38-x01-y01","/DELPHI_1998_I473409/d39-x01-y01",
                                                      "/SLD_2004_S5693039/d07-x01-y01","/SLD_1999_S3743934/d16-x01-y01"]

analyses["IdentifiedParticleFlavour"][413][5][91.2]=["/OPAL_1995_I382219/d04-x01-y01"]
analyses["IdentifiedParticleFlavour"][413][4][91.2]=["/OPAL_1995_I382219/d05-x01-y01"]

analyses["IdentifiedParticleFlavour"][313 ][1][91.2]=["/SLD_1999_S3743934/d14-x01-y01"]
analyses["IdentifiedParticleFlavour"][313 ][4][91.2]=["/SLD_1999_S3743934/d14-x01-y02"]
analyses["IdentifiedParticleFlavour"][313 ][5][91.2]=["/SLD_1999_S3743934/d14-x01-y03"]
analyses["IdentifiedParticleFlavour"][3122][1][91.2]=["/SLD_1999_S3743934/d18-x01-y01"]
analyses["IdentifiedParticleFlavour"][3122][4][91.2]=["/SLD_1999_S3743934/d18-x01-y02"]
analyses["IdentifiedParticleFlavour"][3122][5][91.2]=["/SLD_1999_S3743934/d18-x01-y03"]
analyses["IdentifiedParticleFlavour"][311 ][1][91.2]=["/SLD_1999_S3743934/d20-x01-y01"]
analyses["IdentifiedParticleFlavour"][311 ][4][91.2]=["/SLD_1999_S3743934/d20-x01-y02"]
analyses["IdentifiedParticleFlavour"][311 ][5][91.2]=["/SLD_1999_S3743934/d20-x01-y03"]
analyses["IdentifiedParticleFlavour"][333 ][1][91.2]=["/SLD_1999_S3743934/d22-x01-y01"]
analyses["IdentifiedParticleFlavour"][333 ][4][91.2]=["/SLD_1999_S3743934/d22-x01-y02"]
analyses["IdentifiedParticleFlavour"][333 ][5][91.2]=["/SLD_1999_S3743934/d22-x01-y03"]

analyses["IdentifiedParticleFlavour"][211 ][41][91.2]=["/SLD_1999_S3743934/d11-x01-y01"]
analyses["IdentifiedParticleFlavour"][211 ][51][91.2]=["/SLD_1999_S3743934/d11-x01-y02"]
analyses["IdentifiedParticleFlavour"][321 ][41][91.2]=["/SLD_1999_S3743934/d13-x01-y01"]
analyses["IdentifiedParticleFlavour"][321 ][51][91.2]=["/SLD_1999_S3743934/d13-x01-y02"]
analyses["IdentifiedParticleFlavour"][313 ][41][91.2]=["/SLD_1999_S3743934/d15-x01-y01"]
analyses["IdentifiedParticleFlavour"][313 ][51][91.2]=["/SLD_1999_S3743934/d15-x01-y02"]
analyses["IdentifiedParticleFlavour"][2212][41][91.2]=["/SLD_1999_S3743934/d17-x01-y01"]
analyses["IdentifiedParticleFlavour"][2212][51][91.2]=["/SLD_1999_S3743934/d17-x01-y02"]
analyses["IdentifiedParticleFlavour"][3122][41][91.2]=["/SLD_1999_S3743934/d19-x01-y01"]
analyses["IdentifiedParticleFlavour"][3122][51][91.2]=["/SLD_1999_S3743934/d19-x01-y02"]
analyses["IdentifiedParticleFlavour"][311 ][41][91.2]=["/SLD_1999_S3743934/d21-x01-y01"]
analyses["IdentifiedParticleFlavour"][311 ][51][91.2]=["/SLD_1999_S3743934/d21-x01-y02"]
analyses["IdentifiedParticleFlavour"][333 ][41][91.2]=["/SLD_1999_S3743934/d23-x01-y01"]
analyses["IdentifiedParticleFlavour"][333 ][51][91.2]=["/SLD_1999_S3743934/d23-x01-y02"]

analyses["MultiplicityFlavour"][211 ][41][91.2]=["/SLD_1999_S3743934/d25-x01-y01"]
analyses["MultiplicityFlavour"][211 ][51][91.2]=["/SLD_1999_S3743934/d25-x01-y02"]
analyses["MultiplicityFlavour"][321 ][41][91.2]=["/SLD_1999_S3743934/d25-x02-y01"]
analyses["MultiplicityFlavour"][321 ][51][91.2]=["/SLD_1999_S3743934/d25-x02-y02"]
analyses["MultiplicityFlavour"][311 ][41][91.2]=["/SLD_1999_S3743934/d25-x03-y01"]
analyses["MultiplicityFlavour"][311 ][51][91.2]=["/SLD_1999_S3743934/d25-x03-y02"]
analyses["MultiplicityFlavour"][313 ][41][91.2]=["/SLD_1999_S3743934/d25-x04-y01"]
analyses["MultiplicityFlavour"][313 ][51][91.2]=["/SLD_1999_S3743934/d25-x04-y02"]
analyses["MultiplicityFlavour"][333 ][41][91.2]=["/SLD_1999_S3743934/d25-x05-y01"]
analyses["MultiplicityFlavour"][333 ][51][91.2]=["/SLD_1999_S3743934/d25-x05-y02"]
analyses["MultiplicityFlavour"][2212][41][91.2]=["/SLD_1999_S3743934/d25-x06-y01"]
analyses["MultiplicityFlavour"][2212][51][91.2]=["/SLD_1999_S3743934/d25-x06-y02"]
analyses["MultiplicityFlavour"][3122][41][91.2]=["/SLD_1999_S3743934/d25-x07-y01"]
analyses["MultiplicityFlavour"][3122][51][91.2]=["/SLD_1999_S3743934/d25-x07-y02"]

analyses["MultiplicityFlavour"][211 ][1][91.2]=["/SLD_2004_S5693039/d05-x02-y01","/SLD_1999_S3743934/d24-x01-y02",
                                                "/DELPHI_1998_I473409/d03-x01-y02"]
analyses["MultiplicityFlavour"][211 ][4][91.2]=["/SLD_2004_S5693039/d05-x02-y02","/SLD_1999_S3743934/d24-x01-y03"]
analyses["MultiplicityFlavour"][211 ][5][91.2]=["/SLD_2004_S5693039/d05-x02-y03","/SLD_1999_S3743934/d24-x01-y04",
                                                "/DELPHI_1998_I473409/d02-x01-y02"]
analyses["MultiplicityFlavour"][321 ][1][91.2]=["/SLD_2004_S5693039/d06-x02-y01","/SLD_1999_S3743934/d24-x02-y02",
                                                "/DELPHI_1998_I473409/d03-x01-y03"]
analyses["MultiplicityFlavour"][321 ][4][91.2]=["/SLD_2004_S5693039/d06-x02-y02","/SLD_1999_S3743934/d24-x02-y03"]
analyses["MultiplicityFlavour"][321 ][5][91.2]=["/SLD_2004_S5693039/d06-x02-y03","/SLD_1999_S3743934/d24-x02-y04",
                                                "/DELPHI_1998_I473409/d02-x01-y03"]
analyses["MultiplicityFlavour"][311 ][1][91.2]=["/SLD_1999_S3743934/d24-x03-y02"]
analyses["MultiplicityFlavour"][311 ][4][91.2]=["/SLD_1999_S3743934/d24-x03-y03"]
analyses["MultiplicityFlavour"][311 ][5][91.2]=["/SLD_1999_S3743934/d24-x03-y04"]
analyses["MultiplicityFlavour"][313 ][1][91.2]=["/SLD_1999_S3743934/d24-x04-y02"]
analyses["MultiplicityFlavour"][313 ][4][91.2]=["/SLD_1999_S3743934/d24-x04-y03"]
analyses["MultiplicityFlavour"][313 ][5][91.2]=["/SLD_1999_S3743934/d24-x04-y04"]
analyses["MultiplicityFlavour"][333 ][1][91.2]=["/SLD_1999_S3743934/d24-x05-y02"]
analyses["MultiplicityFlavour"][333 ][4][91.2]=["/SLD_1999_S3743934/d24-x05-y03"]
analyses["MultiplicityFlavour"][333 ][5][91.2]=["/SLD_1999_S3743934/d24-x05-y04"]
analyses["MultiplicityFlavour"][2212][1][91.2]=["/SLD_2004_S5693039/d07-x02-y01","/SLD_1999_S3743934/d24-x06-y02",
                                                "/DELPHI_1998_I473409/d03-x01-y04"]
analyses["MultiplicityFlavour"][2212][4][91.2]=["/SLD_2004_S5693039/d07-x02-y02","/SLD_1999_S3743934/d24-x06-y03"]
analyses["MultiplicityFlavour"][2212][5][91.2]=["/SLD_2004_S5693039/d07-x02-y03","/SLD_1999_S3743934/d24-x06-y04",
                                                "/DELPHI_1998_I473409/d02-x01-y04"]
analyses["MultiplicityFlavour"][3122][1][91.2]=["/SLD_1999_S3743934/d24-x07-y02"]
analyses["MultiplicityFlavour"][3122][4][91.2]=["/SLD_1999_S3743934/d24-x07-y03"]
analyses["MultiplicityFlavour"][3122][5][91.2]=["/SLD_1999_S3743934/d24-x07-y04"]

analyses["Charged"]["ChargedSpectrum"][1][91.2]=["/DELPHI_1998_I473409/d32-x01-y01","/DELPHI_1998_I473409/d33-x01-y01",
                                                 "/DELPHI_1997_I428178/d01-x01-y03","/OPAL_1998_S3780481/d01-x01-y01",
                                                 "/OPAL_1998_S3780481/d05-x01-y01","/SLD_2004_S5693039/d08-x01-y01"]
analyses["Charged"]["ChargedSpectrum"][2][91.2]=["/L3_2004_I652683/d65-x01-y02"]
analyses["Charged"]["ChargedSpectrum"][4][91.2]=["/DELPHI_1997_I428178/d01-x01-y02","/OPAL_1998_S3780481/d02-x01-y01",
                                                 "/OPAL_1998_S3780481/d06-x01-y01","/SLD_2004_S5693039/d08-x01-y02"]
analyses["Charged"]["ChargedSpectrum"][5][91.2]=["/DELPHI_1998_I473409/d24-x01-y01","/DELPHI_1998_I473409/d25-x01-y01",
                                                 "/DELPHI_1997_I428178/d01-x01-y01","/OPAL_1998_S3780481/d03-x01-y01",
                                                 "/OPAL_1998_S3780481/d07-x01-y01","/SLD_2004_S5693039/d08-x01-y03","/L3_2004_I652683/d65-x01-y03"]


# identified particle distributions
# photons
# x_E
analyses["IdentifiedParticle"][22  ]["x" ][14.0]=["/CELLO_1983_I191415/d01-x01-y01"]
analyses["IdentifiedParticle"][22  ]["x" ][22.0]=["/CELLO_1983_I191415/d02-x01-y01"]
analyses["IdentifiedParticle"][22  ]["x" ][29.0]=["/TPC_1985_I205868/d01-x01-y01" ]
analyses["IdentifiedParticle"][22  ]["x" ][34.0]=["/CELLO_1983_I191415/d03-x01-y01"]
analyses["IdentifiedParticle"][22  ]["x" ][35.0]=["/CELLO_1989_I276764/d02-x01-y01","/JADE_1990_I282847/d01-x01-y01"]
analyses["IdentifiedParticle"][22  ]["x" ][44.0]=["/JADE_1990_I282847/d02-x01-y01"]
analyses["IdentifiedParticle"][22  ]["x" ][91.2]=["/OPAL_1998_S3749908/d02-x01-y01"]
# xi
analyses["IdentifiedParticle"][22  ]["xi"][91.2]=["/ALEPH_1996_S3486095/d28-x01-y01","/OPAL_1998_S3749908/d03-x01-y01"]
# charged pions
# x
analyses["IdentifiedParticle"][211 ]["x" ][10.0 ] = ["/ARGUS_1989_I276860/d09-x01-y02"]
analyses["IdentifiedParticle"][211 ]["x" ][10.52] = ["/BELLE_2013_I1216515/d01-x01-y01"]
analyses["IdentifiedParticle"][211 ]["x" ][12.0 ] = ["/TASSO_1980_I153656/d02-x01-y02"]
analyses["IdentifiedParticle"][211 ]["x" ][14.0 ] = ["/TASSO_1983_I181470/d20-x01-y01"]
analyses["IdentifiedParticle"][211 ]["x" ][22.0 ] = ["/TASSO_1983_I181470/d22-x01-y01"]
analyses["IdentifiedParticle"][211 ]["x" ][29.0 ] = ["/TPC_1988_I262143/d01-x01-y01","/TPC_1988_I262143/d05-x01-y01"]
analyses["IdentifiedParticle"][211 ]["x" ][30.0 ] = ["/TASSO_1980_I153656/d05-x01-y02"]
analyses["IdentifiedParticle"][211 ]["x" ][34.0 ] = ["/TASSO_1989_I267755/d07-x01-y01","/TASSO_1983_I181470/d24-x01-y01"]
analyses["IdentifiedParticle"][211 ]["x" ][44.0 ] = ["/TASSO_1989_I267755/d10-x01-y01"]
analyses["IdentifiedParticle"][211 ]["x" ][91.2 ] = ["/ALEPH_1995_I382179/d01-x01-y01","/DELPHI_1998_I473409/d19-x01-y01",
                                                     "/ALEPH_1996_S3486095/d25-x01-y01","/SLD_2004_S5693039/d02-x01-y02",
                                                     "/SLD_1999_S3743934/d01-x01-y02"]
analyses["IdentifiedParticle"][211 ]["Other" ][91.2 ] = ["/SLD_1999_S3743934/d26-x01-y01","/SLD_1999_S3743934/d26-x01-y02",
                                                         "/SLD_1999_S3743934/d27-x01-y01","/SLD_2004_S5693039/d09-x01-y01",
                                                         "/SLD_2004_S5693039/d09-x01-y02","/SLD_2004_S5693039/d09-x01-y03",]
# p
analyses["IdentifiedParticle"][211 ]["p" ][10.0 ] = ["/ARGUS_1989_I276860/d05-x01-y02"]
analyses["IdentifiedParticle"][211 ]["p" ][10.54] = ["/BABAR_2013_I1238276/d01-x01-y01","/BABAR_2013_I1238276/d02-x01-y01"]
analyses["IdentifiedParticle"][211 ]["p" ][12.0 ] = ["/TASSO_1980_I153656/d02-x01-y01"]
analyses["IdentifiedParticle"][211 ]["p" ][14.0 ] = ["/TASSO_1983_I181470/d19-x01-y01"]
analyses["IdentifiedParticle"][211 ]["p" ][22.0 ] = ["/TASSO_1983_I181470/d25-x01-y01"]
analyses["IdentifiedParticle"][211 ]["p" ][30.0 ] = ["/TASSO_1980_I153656/d05-x01-y01"]
analyses["IdentifiedParticle"][211 ]["p" ][34.0 ] = ["/TASSO_1983_I181470/d13-x01-y01"]
analyses["IdentifiedParticle"][211 ]["p" ][91.2 ] = ["/DELPHI_1998_I473409/d18-x01-y01","/OPAL_1994_S2927284/d01-x01-y01"]
# xi
analyses["IdentifiedParticle"][211 ]["xi"][58.0 ] = ["/TOPAZ_1995_I381900/d02-x01-y01"]
# ratios
analyses["IdentifiedParticle"][211 ]["Ratio"][12.0 ] = ["/TASSO_1980_I153656/d08-x01-y01"]
analyses["IdentifiedParticle"][211 ]["Ratio"][29.0 ] = ["/TPC_1988_I262143/d06-x01-y01"]
analyses["IdentifiedParticle"][211 ]["Ratio"][30.0 ] = ["/TASSO_1980_I153656/d11-x01-y01"]
analyses["IdentifiedParticle"][211 ]["Ratio"][34.0 ] = ["/TASSO_1989_I267755/d01-x01-y01"]
analyses["IdentifiedParticle"][211 ]["Ratio"][44.0 ] = ["/TASSO_1989_I267755/d04-x01-y01"]
analyses["IdentifiedParticle"][211 ]["Ratio"][91.2 ] = ["/DELPHI_1998_I473409/d04-x01-y01","/SLD_1999_S3743934/d01-x01-y01"]
# neutral pions
# x
analyses["IdentifiedParticle"][111 ]["x" ][10.0]=["/ARGUS_1990_I278933/d03-x01-y01","/ARGUS_1990_I278933/d03-x01-y02"]
analyses["IdentifiedParticle"][111 ]["x" ][14.0]=["/TASSO_1982_I168232/d02-x03-y03","/CELLO_1983_I191415/d04-x01-y01"]
analyses["IdentifiedParticle"][111 ]["x" ][22.0]=["/CELLO_1983_I191415/d05-x01-y01"]
analyses["IdentifiedParticle"][111 ]["x" ][29.0]=["/TPC_1985_I205868/d02-x01-y01" ]
analyses["IdentifiedParticle"][111 ]["x" ][34.0]=["/TASSO_1982_I168232/d03-x03-y03","/TASSO_1986_I230950/d02-x01-y01",
                                                  "/CELLO_1983_I191415/d06-x01-y01"]
analyses["IdentifiedParticle"][111 ]["x" ][35.0]=["/CELLO_1989_I276764/d03-x01-y01","/CELLO_1989_I276764/d04-x01-y01",
                                                  "/JADE_1990_I282847/d03-x01-y01"]
analyses["IdentifiedParticle"][111 ]["x" ][44.0]=["/TASSO_1989_I267755/d13-x01-y01","/JADE_1990_I282847/d04-x01-y01"]
analyses["IdentifiedParticle"][111 ]["x" ][91.2]=["/DELPHI_1996_I401100/d01-x01-y01","/ALEPH_1996_S3486095/d29-x01-y01",
                                                  "/OPAL_1998_S3749908/d04-x01-y01",]
# p/E
analyses["IdentifiedParticle"][111 ]["p" ][14.0]=["/TASSO_1982_I168232/d02-x01-y01","/TASSO_1982_I168232/d02-x02-y02"]
analyses["IdentifiedParticle"][111 ]["p" ][34.0]=["/TASSO_1982_I168232/d03-x01-y01","/TASSO_1982_I168232/d03-x02-y02",
                                                  "/TASSO_1986_I230950/d01-x01-y01"]
# xi
analyses["IdentifiedParticle"][111 ]["xi"][91.2]=["/OPAL_1998_S3749908/d05-x01-y01"]
# eta
#x
analyses["IdentifiedParticle"][221 ]["x" ][10.0]=["/ARGUS_1990_I278933/d05-x01-y01","/ARGUS_1990_I278933/d05-x01-y02"]
analyses["IdentifiedParticle"][221 ]["x" ][29.0]=["/HRS_1988_I250824/d01-x01-y01"]
analyses["IdentifiedParticle"][221 ]["x" ][35.0]=["/CELLO_1989_I276764/d05-x01-y01" ,"/JADE_1990_I282847/d05-x01-y01"]
analyses["IdentifiedParticle"][221 ]["x" ][91.2]=["/ALEPH_2002_S4823664/d02-x01-y02","/L3_1992_I336180/d01-x01-y01",
                                                  "/ALEPH_1996_S3486095/d30-x01-y01","/OPAL_1998_S3749908/d06-x01-y01",]
# xi
analyses["IdentifiedParticle"][221 ]["xi"][91.2]=["/L3_1992_I336180/d02-x01-y01","/OPAL_1998_S3749908/d07-x01-y01"]
# eta'
# x
analyses["IdentifiedParticle"][331 ]["x" ][91.2]=["/L3_1997_I427107/d07-x01-y01"    ,"/L3_1997_I427107/d09-x01-y01",
                                                  "/ALEPH_1996_S3486095/d31-x01-y01","/OPAL_1998_S3749908/d12-x01-y01"]
# xi
analyses["IdentifiedParticle"][331 ]["xi"][91.2]=["/L3_1997_I427107/d08-x01-y01","/L3_1997_I427107/d10-x01-y01",
                                                  "/OPAL_1998_S3749908/d13-x01-y01"]
# rho +/-
analyses["IdentifiedParticle"][213 ]["x" ][91.2] = ["/OPAL_1998_S3749908/d08-x01-y01"]
analyses["IdentifiedParticle"][213 ]["xi"][91.2] = ["/OPAL_1998_S3749908/d09-x01-y01"]
# rho0
analyses["IdentifiedParticle"][113 ]["x" ][10.0] = ["/ARGUS_1993_S2789213/d10-x01-y01"]
analyses["IdentifiedParticle"][113 ]["x" ][35.0] = ["/JADE_1984_I203145/d02-x01-y01"]
analyses["IdentifiedParticle"][113 ]["x" ][91.2] = ["/DELPHI_1999_S3960137/d01-x01-y01","/ALEPH_1996_S3486095/d37-x01-y01"]
# omega
analyses["IdentifiedParticle"][223 ]["x" ][10.0] = ["/ARGUS_1993_S2789213/d13-x01-y01"]
analyses["IdentifiedParticle"][223 ]["x" ][91.2] = ["/ALEPH_2002_S4823664/d03-x01-y02","/L3_1997_I427107/d05-x01-y01",
                                                    "/ALEPH_1996_S3486095/d38-x01-y01","/OPAL_1998_S3749908/d10-x01-y01",]
analyses["IdentifiedParticle"][223 ]["xi"][91.2] = ["/L3_1997_I427107/d06-x01-y01","/OPAL_1998_S3749908/d11-x01-y01"]
# phi
analyses["IdentifiedParticle"][333 ]["x"    ][10.0] = ["/ARGUS_1989_I262551/d01-x01-y01"]
analyses["IdentifiedParticle"][333 ]["x"    ][29.0] = ["/TPC_1984_I200105/d01-x01-y01"]
analyses["IdentifiedParticle"][333 ]["x"    ][91.2] = ["/DELPHI_1996_I420528/d03-x01-y01","/ALEPH_1996_S3486095/d40-x01-y01",
                                                       "/OPAL_1998_S3702294/d02-x01-y03","/SLD_1999_S3743934/d09-x01-y01"]
analyses["IdentifiedParticle"][333 ]["Other"][29.0] = ["/TPC_1984_I200105/d03-x01-y01"]
# f_2
analyses["IdentifiedParticle"][225]["x"    ][91.2]=["/DELPHI_1999_S3960137/d01-x01-y03","/OPAL_1998_S3702294/d02-x01-y02"]
# f_2'
analyses["IdentifiedParticle"][335]["x"    ][91.2]=["/DELPHI_1996_I416741/d01-x01-y01"]
# f_0(980)
analyses["IdentifiedParticle"][9010221]["x"    ][10.0]=["/ARGUS_1993_S2669951/d02-x01-y01"]
analyses["IdentifiedParticle"][9010221]["x"    ][91.2]=["/DELPHI_1999_S3960137/d01-x01-y02","/OPAL_1998_S3702294/d02-x01-y01"]
# a_0 =/-
analyses["IdentifiedParticle"][9000211]["x"    ][91.2]=["/OPAL_1998_S3749908/d14-x01-y01"]
analyses["IdentifiedParticle"][9000211]["xi"   ][91.2]=["/OPAL_1998_S3749908/d15-x01-y01"]
# strange mesons
# K0
# x
analyses["IdentifiedParticle"][311 ]["x"    ][3.63] = ["/PLUTO_1977_I118873/d02-x01-y01"]
analyses["IdentifiedParticle"][311 ]["x"    ][4.03] = ["/PLUTO_1977_I118873/d03-x01-y01"]
analyses["IdentifiedParticle"][311 ]["x"    ][4.5]  = ["/PLUTO_1977_I118873/d04-x01-y01"]
analyses["IdentifiedParticle"][311 ]["x"    ][9.4]  = ["/PLUTO_1981_I165122/d05-x01-y01"]
analyses["IdentifiedParticle"][311 ]["x"    ][10.0] = ["/ARGUS_1989_I276860/d11-x01-y02"]
analyses["IdentifiedParticle"][311 ]["x"    ][14.0] = ["/TASSO_1980_I153341/d04-x01-y01","/TASSO_1985_I205119/d01-x01-y01"]
analyses["IdentifiedParticle"][311 ]["x"    ][22.0] = ["/TASSO_1985_I205119/d02-x01-y01"]
analyses["IdentifiedParticle"][311 ]["x"    ][29.0] = ["/TPC_1984_I205869/d04-x01-y01","/HRS_1990_I280958/d03-x01-y01"]
analyses["IdentifiedParticle"][311 ]["x"    ][30.0] = ["/PLUTO_1981_I165122/d04-x01-y01"]
analyses["IdentifiedParticle"][311 ]["x"    ][34.0] = ["/TASSO_1985_I205119/d03-x01-y01"]
analyses["IdentifiedParticle"][311 ]["x"    ][34.5] = ["/TASSO_1990_I284251/d01-x01-y03"]
analyses["IdentifiedParticle"][311 ]["x"    ][35.0] = ["/TASSO_1990_I284251/d01-x01-y02","/CELLO_1990_I283026/d01-x01-y01"]
analyses["IdentifiedParticle"][311 ]["x"    ][42.6] = ["/TASSO_1990_I284251/d01-x01-y01"]
analyses["IdentifiedParticle"][311 ]["x"    ][91.2] = ["/OPAL_2000_S4418603/d03-x01-y01","/DELPHI_1995_I377487/d08-x01-y01",
                                                       "/ALEPH_1996_S3486095/d32-x01-y01","/SLD_1999_S3743934/d05-x01-y01"]
# p
analyses["IdentifiedParticle"][311 ]["p"    ][10.0] = ["/ARGUS_1989_I276860/d07-x01-y02"]
analyses["IdentifiedParticle"][311 ]["p"    ][14.0] = ["/TASSO_1980_I153341/d02-x01-y01","/TASSO_1985_I205119/d07-x01-y01"]
analyses["IdentifiedParticle"][311 ]["p"    ][22.0] = ["/TASSO_1985_I205119/d08-x01-y01"]
analyses["IdentifiedParticle"][311 ]["p"    ][34.0] = ["/TASSO_1985_I205119/d09-x01-y01"]
# xi
analyses["IdentifiedParticle"][311 ]["xi"   ][58.0] = ["/TOPAZ_1995_I381900/d03-x01-y01"]
analyses["IdentifiedParticle"][311 ]["xi"   ][91.2] = ["/DELPHI_1995_I377487/d09-x01-y01"]
# other
analyses["IdentifiedParticle"][311 ]["Other"][14.8 ] = ["/TASSO_1990_I284251/d06-x01-y01","/TASSO_1990_I284251/d06-x01-y02"]
analyses["IdentifiedParticle"][311 ]["Other"][21.5 ] = ["/TASSO_1990_I284251/d07-x01-y01","/TASSO_1990_I284251/d07-x01-y02"]
analyses["IdentifiedParticle"][311 ]["Other"][29.0 ] = ["/HRS_1990_I280958/d04-x01-y01"]
analyses["IdentifiedParticle"][311 ]["Other"][35.0 ] = ["/TASSO_1990_I284251/d05-x01-y03","/TASSO_1990_I284251/d05-x01-y04",""]
analyses["IdentifiedParticle"][311 ]["Other"][42.6 ] = ["/TASSO_1990_I284251/d05-x01-y01","/TASSO_1990_I284251/d05-x01-y02"]
# K+/-
# x
analyses["IdentifiedParticle"][321 ]["x"    ][10.0 ] = ["/ARGUS_1989_I276860/d10-x01-y02"]
analyses["IdentifiedParticle"][321 ]["x"    ][10.52] = ["/BELLE_2013_I1216515/d01-x01-y02"]
analyses["IdentifiedParticle"][321 ]["x"    ][12.0 ] = ["/TASSO_1980_I153656/d03-x01-y02"]
analyses["IdentifiedParticle"][321 ]["x"    ][14.0 ] = ["/TASSO_1983_I181470/d26-x01-y01"]
analyses["IdentifiedParticle"][321 ]["x"    ][22.0 ] = ["/TASSO_1983_I181470/d10-x01-y01"]
analyses["IdentifiedParticle"][321 ]["x"    ][29.0 ] = ["/TPC_1988_I262143/d01-x01-y02","/TPC_1988_I262143/d05-x01-y02"]
analyses["IdentifiedParticle"][321 ]["x"    ][30.0 ] = ["/TASSO_1980_I153656/d06-x01-y02"]
analyses["IdentifiedParticle"][321 ]["x"    ][34.0 ] = ["/TASSO_1989_I267755/d08-x01-y01","/TASSO_1983_I181470/d12-x01-y01"]
analyses["IdentifiedParticle"][321 ]["x"    ][44.0 ] = ["/TASSO_1989_I267755/d11-x01-y01"]
analyses["IdentifiedParticle"][321 ]["x"    ][91.2 ] = ["/ALEPH_1995_I382179/d02-x01-y01","/DELPHI_1995_I394052/d05-x01-y01",
                                                        "/DELPHI_1998_I473409/d21-x01-y01","/SLD_1999_S3743934/d02-x01-y02",
                                                        "/ALEPH_1996_S3486095/d26-x01-y01","/SLD_2004_S5693039/d03-x01-y02"]
# p
analyses["IdentifiedParticle"][321 ]["p"    ][10.0 ] = ["/ARGUS_1989_I276860/d06-x01-y02"]
analyses["IdentifiedParticle"][321 ]["p"    ][10.54] = ["/BABAR_2013_I1238276/d01-x01-y02","/BABAR_2013_I1238276/d02-x01-y02"]
analyses["IdentifiedParticle"][321 ]["p"    ][12.0 ] = ["/TASSO_1980_I153656/d03-x01-y01"]
analyses["IdentifiedParticle"][321 ]["p"    ][14.0 ] = ["/TASSO_1983_I181470/d21-x01-y01"]
analyses["IdentifiedParticle"][321 ]["p"    ][22.0 ] = ["/TASSO_1983_I181470/d27-x01-y01"]
analyses["IdentifiedParticle"][321 ]["p"    ][30.0 ] = ["/TASSO_1980_I153656/d06-x01-y01"]
analyses["IdentifiedParticle"][321 ]["p"    ][34.0 ] = ["/TASSO_1983_I181470/d15-x01-y01"]
analyses["IdentifiedParticle"][321 ]["p"    ][91.2 ] = ["/DELPHI_1995_I394052/d03-x01-y01","/DELPHI_1998_I473409/d20-x01-y01",
                                                        "/OPAL_1994_S2927284/d02-x01-y01"]
# xi
analyses["IdentifiedParticle"][321 ]["xi"   ][58.0 ] = ["/TOPAZ_1995_I381900/d02-x01-y02"]
# ratio
analyses["IdentifiedParticle"][321 ]["Ratio"][12.0 ] = ["/TASSO_1980_I153656/d09-x01-y01"]
analyses["IdentifiedParticle"][321 ]["Ratio"][29.0 ] = ["/TPC_1988_I262143/d07-x01-y01","/TPC_1988_I262143/d06-x01-y02"]
analyses["IdentifiedParticle"][321 ]["Ratio"][30.0 ] = ["/TASSO_1980_I153656/d12-x01-y01"]
analyses["IdentifiedParticle"][321 ]["Ratio"][34.0 ] = ["/TASSO_1989_I267755/d02-x01-y01"]
analyses["IdentifiedParticle"][321 ]["Ratio"][44.0 ] = ["/TASSO_1989_I267755/d05-x01-y01"]
analyses["IdentifiedParticle"][321 ]["Ratio"][91.2 ] = ["/DELPHI_1998_I473409/d05-x01-y01","/SLD_1999_S3743934/d02-x01-y01"]
# other
analyses["IdentifiedParticle"][321 ]["Other"][91.2 ] = ["/SLD_1999_S3743934/d30-x01-y01","/SLD_1999_S3743934/d30-x01-y02",
                                                        "/SLD_1999_S3743934/d31-x01-y01","/SLD_2004_S5693039/d10-x01-y01",
                                                        "/SLD_2004_S5693039/d10-x01-y02","/SLD_2004_S5693039/d10-x01-y03"]
# K*0
analyses["IdentifiedParticle"][313 ]["x"    ][10.0 ] = ["/ARGUS_1993_S2789213/d07-x01-y01"]
analyses["IdentifiedParticle"][313 ]["x"    ][29.0 ] = ["/TPC_1984_I205869/d03-x01-y01"]
analyses["IdentifiedParticle"][313 ]["x"    ][91.2 ] = ["/DELPHI_1996_I420528/d01-x01-y01","/ALEPH_1996_S3486095/d39-x01-y01",
                                                        "/OPAL_1997_S3608263/d01-x01-y01","/SLD_1999_S3743934/d08-x01-y01"]
analyses["IdentifiedParticle"][313 ]["Other"][91.2 ] = ["/SLD_1999_S3743934/d28-x01-y01","/SLD_1999_S3743934/d28-x01-y02",
                                                        "/SLD_1999_S3743934/d29-x01-y01"]
# K* +/-
analyses["IdentifiedParticle"][323 ]["x"    ][10.0 ] = ["/ARGUS_1993_S2789213/d04-x01-y01"]
analyses["IdentifiedParticle"][323 ]["x"    ][14.8 ] = ["/TASSO_1990_I284251/d02-x01-y01"]
analyses["IdentifiedParticle"][323 ]["x"    ][21.5 ] = ["/TASSO_1990_I284251/d03-x01-y01"]
analyses["IdentifiedParticle"][323 ]["x"    ][34.5 ] = ["/TASSO_1990_I284251/d08-x01-y03"]
analyses["IdentifiedParticle"][323 ]["x"    ][35.0 ] = ["/TASSO_1990_I284251/d08-x01-y02","/CELLO_1990_I283026/d02-x01-y01",
                                                        "/JADE_1984_I203145/d03-x01-y01"]
analyses["IdentifiedParticle"][323 ]["x"    ][42.6 ] = ["/TASSO_1990_I284251/d08-x01-y01"]
analyses["IdentifiedParticle"][323 ]["x"    ][91.2 ] = ["/OPAL_1993_I342766/d01-x01-y01","/DELPHI_1995_I377487/d10-x01-y01",
                                                        "/ALEPH_1996_S3486095/d43-x01-y01"]
analyses["IdentifiedParticle"][323 ]["Other"][35.0 ] = ["/TASSO_1990_I284251/d10-x01-y01","/TASSO_1990_I284251/d10-x01-y02"]
# charm
# D+/-
analyses["IdentifiedParticle"][421 ]["x"    ][10.5 ] = ["/CLEO_2004_S5809304/d03-x01-y01","/CLEO_2004_S5809304/d04-x01-y01"]
analyses["IdentifiedParticle"][421 ]["x"    ][29.0 ] = ["/HRS_1988_I23360/d02-x01-y01"]
# D0
analyses["IdentifiedParticle"][411 ]["x"    ][10.5 ] = ["/CLEO_2004_S5809304/d02-x01-y01","/CLEO_2004_S5809304/d09-x01-y01"]
analyses["IdentifiedParticle"][411 ]["x"    ][29.0 ] = ["/HRS_1988_I23360/d02-x01-y02"]
# D* 0
analyses["IdentifiedParticle"][423 ]["x"    ][10.5]=["/CLEO_2004_S5809304/d07-x01-y01","/CLEO_2004_S5809304/d08-x01-y01"]
# D* +/-
analyses["IdentifiedParticle"][413 ]["x"    ][29.0 ] = ["/TPC_1986_I217416/d01-x01-y01","/TPC_1986_I217416/d01-x01-y02",
                                                        "/HRS_1988_I23360/d01-x01-y01","/HRS_1988_I23360/d01-x01-y02"]
analyses["IdentifiedParticle"][413 ]["x"    ][34.4 ] = ["/JADE_1984_I202785/d01-x01-y01"]
analyses["IdentifiedParticle"][413 ]["x"    ][36.2 ] = ["/TASSO_1989_I278856/d01-x01-y01","/TASSO_1989_I278856/d01-x01-y02",
                                                        "/TASSO_1989_I278856/d02-x01-y01","/TASSO_1989_I278856/d02-x01-y02"]
analyses["IdentifiedParticle"][413 ]["x"    ][91.2 ] = ["/ALEPH_1999_S4193598/d01-x01-y01"]
analyses["IdentifiedParticle"][413 ]["x"    ][10.5 ] = ["/CLEO_2004_S5809304/d05-x01-y01","/CLEO_2004_S5809304/d06-x01-y01"]
analyses["IdentifiedParticle"][413 ]["Other"][34.4 ] = ["/JADE_1984_I202785/d03-x01-y01"]
# D_2
analyses["IdentifiedParticle"][425 ]["x"    ][10.0 ] = ["/ARGUS_1989_I268577/d02-x01-y01"]
# D_s+
analyses["IdentifiedParticle"][431 ]["x"    ][10.5 ] = ["/CLEO_2000_I526554/d02-x01-y01","/CLEO_2000_I526554/d04-x01-y01"]
# D_s*+
analyses["IdentifiedParticle"][433 ]["x"    ][10.5 ] = ["/CLEO_2000_I526554/d01-x01-y01","/CLEO_2000_I526554/d03-x01-y01"]
# charmonium
analyses["IdentifiedParticle"][443 ]["x"    ][91.2 ] = ["/OPAL_1996_S3257789/d01-x01-y01"]
#
#  Baryons
#
# light unflavoured
# proton
# x
analyses["IdentifiedParticle"][2212]["x"    ][10.0 ] = ["/ARGUS_1989_I276860/d12-x01-y02"]
analyses["IdentifiedParticle"][2212]["x"    ][12.0 ] = ["/TASSO_1980_I153656/d04-x01-y02"]
analyses["IdentifiedParticle"][2212]["x"    ][14.0 ] = ["/TASSO_1983_I181470/d14-x01-y01"]
analyses["IdentifiedParticle"][2212]["x"    ][22.0 ] = ["/TASSO_1983_I181470/d16-x01-y01"]
analyses["IdentifiedParticle"][2212]["x"    ][29.0 ] = ["/TPC_1988_I262143/d01-x01-y03","/TPC_1988_I262143/d05-x01-y03"]
analyses["IdentifiedParticle"][2212]["x"    ][30.0 ] = ["/TASSO_1980_I153656/d07-x01-y02"]
analyses["IdentifiedParticle"][2212]["x"    ][34.0 ] = ["/TASSO_1989_I267755/d09-x01-y01","/TASSO_1983_I181470/d18-x01-y01"]
analyses["IdentifiedParticle"][2212]["x"    ][44.0 ] = ["/TASSO_1989_I267755/d12-x01-y01"]
analyses["IdentifiedParticle"][2212]["x"    ][91.2 ] = ["/ALEPH_1995_I382179/d03-x01-y01","/DELPHI_1995_I394052/d06-x01-y01",
                                                        "/DELPHI_1998_I473409/d23-x01-y01","/ALEPH_1996_S3486095/d27-x01-y01",
                                                        "/SLD_2004_S5693039/d04-x01-y02","/SLD_1999_S3743934/d03-x01-y02"]
# p
analyses["IdentifiedParticle"][2212]["p"    ][10.0 ] = ["/ARGUS_1989_I276860/d08-x01-y02"]
analyses["IdentifiedParticle"][2212]["p"    ][10.54] = ["/BABAR_2013_I1238276/d01-x01-y03","/BABAR_2013_I1238276/d02-x01-y03"] 
analyses["IdentifiedParticle"][2212]["p"    ][12.0 ] = ["/TASSO_1980_I153656/d04-x01-y01"]
analyses["IdentifiedParticle"][2212]["p"    ][14.0 ] = ["/TASSO_1983_I181470/d23-x01-y01"]
analyses["IdentifiedParticle"][2212]["p"    ][22.0 ] = ["/TASSO_1983_I181470/d11-x01-y01"]
analyses["IdentifiedParticle"][2212]["p"    ][30.0 ] = ["/TASSO_1980_I153656/d07-x01-y01"]
analyses["IdentifiedParticle"][2212]["p"    ][34.0 ] = ["/TASSO_1989_I267755/d03-x01-y01","/JADE_1981_I166363/d01-x01-y01",
                                                        "/TASSO_1983_I181470/d17-x01-y01"]
analyses["IdentifiedParticle"][2212]["p"    ][44.0 ] = ["/TASSO_1989_I267755/d06-x01-y01"]
analyses["IdentifiedParticle"][2212]["p"    ][91.2 ] = ["/DELPHI_1995_I394052/d04-x01-y01","/DELPHI_1998_I473409/d22-x01-y01",
                                                        "/OPAL_1994_S2927284/d03-x01-y01",]
# xi
analyses["IdentifiedParticle"][2212]["xi"   ][58.0 ] = ["/TOPAZ_1995_I381900/d02-x01-y03"]
# ratio
analyses["IdentifiedParticle"][2212]["Ratio"][12.0 ] = ["/TASSO_1980_I153656/d10-x01-y01"]
analyses["IdentifiedParticle"][2212]["Ratio"][29.0 ] = ["/TPC_1988_I262143/d06-x01-y03","/TPC_1988_I262143/d07-x01-y02",
                                                        "/TPC_1988_I262143/d07-x01-y03"]
analyses["IdentifiedParticle"][2212]["Ratio"][30.0 ] = ["/TASSO_1980_I153656/d13-x01-y01"]
analyses["IdentifiedParticle"][2212]["Ratio"][91.2 ] = ["/SLD_1999_S3743934/d03-x01-y01","/DELPHI_1998_I473409/d06-x01-y01"]
analyses["IdentifiedParticle"][2212]["Other"][91.2 ] = ["/SLD_1999_S3743934/d32-x01-y01","/SLD_1999_S3743934/d32-x01-y02",
                                                        "/SLD_1999_S3743934/d33-x01-y01","/SLD_2004_S5693039/d11-x01-y01",
                                                        "/SLD_2004_S5693039/d11-x01-y02","/SLD_2004_S5693039/d11-x01-y03"]
# Delta++
analyses["IdentifiedParticle"][2224]["x"    ][91.2 ] = ["/OPAL_1995_S3198391/d01-x01-y01","/DELPHI_1995_I399737/d01-x01-y01"]
# hyperons
# lambda0
# x
analyses["IdentifiedParticle"][3122]["x"    ][10.0 ] = ["/ARGUS_1988_I251097/d05-x01-y01","/ARGUS_1988_I251097/d06-x01-y01"]
analyses["IdentifiedParticle"][3122]["x"    ][10.52] = ["/BELLE_2017_I1606201/d01-x01-y01"]
analyses["IdentifiedParticle"][3122]["x"    ][14.0 ] = ["/TASSO_1985_I205119/d04-x01-y01"]
analyses["IdentifiedParticle"][3122]["x"    ][22.0 ] = ["/TASSO_1985_I205119/d05-x01-y01"]
analyses["IdentifiedParticle"][3122]["x"    ][29.0 ] = ["/HRS_1992_I339573/d01-x01-y01"]
analyses["IdentifiedParticle"][3122]["x"    ][34.0 ] = ["/TASSO_1985_I205119/d06-x01-y01"] 
analyses["IdentifiedParticle"][3122]["x"    ][34.8 ] = ["/TASSO_1989_I266893/d08-x01-y01"]
analyses["IdentifiedParticle"][3122]["x"    ][35.0 ] = ["/CELLO_1990_I283026/d03-x01-y01"]
analyses["IdentifiedParticle"][3122]["x"    ][91.2 ] = ["/OPAL_1997_S3396100/d01-x01-y01","/ALEPH_1996_S3486095/d33-x01-y01",
                                                        "/DELPHI_1993_I360638/d01-x01-y01","/SLD_1999_S3743934/d07-x01-y01"]
# p
analyses["IdentifiedParticle"][3122]["p"    ][14.0 ] = ["/TASSO_1985_I205119/d10-x01-y01"]
analyses["IdentifiedParticle"][3122]["p"    ][22.0 ] = ["/TASSO_1985_I205119/d11-x01-y01"]
analyses["IdentifiedParticle"][3122]["p"    ][34.0 ] = ["/JADE_1981_I166363/d02-x01-y01","/TASSO_1985_I205119/d12-x01-y01"]
analyses["IdentifiedParticle"][3122]["p"    ][34.8 ] = ["/TASSO_1989_I266893/d03-x01-y01"]
# xi
analyses["IdentifiedParticle"][3122]["xi"   ][91.2 ] = ["/OPAL_1997_S3396100/d02-x01-y01"]
# other
analyses["IdentifiedParticle"][3122]["Other"][34.8 ] = ["/TASSO_1989_I266893/d04-x01-y01","/TASSO_1989_I266893/d05-x01-y01",
                                                        "/TASSO_1989_I266893/d06-x01-y01","/TASSO_1989_I266893/d07-x01-y01",
                                                        "/TASSO_1989_I266893/d15-x01-y01","/TASSO_1989_I266893/d15-x01-y02",
                                                        "/TASSO_1989_I266893/d15-x01-y03"]
analyses["IdentifiedParticle"][3122]["Other"][91.2 ] = ["/DELPHI_1993_I360638/d03-x01-y01","/DELPHI_1993_I360638/d04-x01-y01",
                                                        "/DELPHI_1993_I360638/d05-x01-y01","/DELPHI_1993_I360638/d06-x01-y01",
                                                        "/SLD_1999_S3743934/d34-x01-y01","/SLD_1999_S3743934/d34-x01-y02",
                                                        "/SLD_1999_S3743934/d35-x01-y01"]
# Sigma+
analyses["IdentifiedParticle"][3222]["x"    ][91.2 ] = ["/OPAL_1997_I421977/d01-x01-y01"]
# sigma0
analyses["IdentifiedParticle"][3212]["x"    ][10.52] = ["/BELLE_2017_I1606201/d02-x01-y01"]
# sigma-
analyses["IdentifiedParticle"][3112]["x"    ][91.2 ] = ["/OPAL_1997_I421977/d02-x01-y01","/DELPHI_2000_I524694/d01-x01-y01"]
# Sigma*+
analyses["IdentifiedParticle"][3224   ]["x"    ][10.52] = ["/BELLE_2017_I1606201/d03-x01-y01"]
analyses["IdentifiedParticle"][3224   ]["x"    ][91.2 ] = ["/DELPHI_1995_S3137023/d03-x01-y01","/OPAL_1997_S3396100/d05-x01-y01"]
analyses["IdentifiedParticle"][3224   ]["xi"   ][91.2 ] = ["/OPAL_1997_S3396100/d06-x01-y01"]
analyses["IdentifiedParticle"]["3224B"]["x"    ][91.2 ] = ["/ALEPH_1996_S3486095/d35-x01-y01"]
# sigma*-
analyses["IdentifiedParticle"][3114   ]["x"    ][91.2 ] = ["/OPAL_1997_S3396100/d07-x01-y01"]
analyses["IdentifiedParticle"][3114   ]["xi"   ][91.2 ] = ["/OPAL_1997_S3396100/d08-x01-y01"]
# xi-
analyses["IdentifiedParticle"][3312]["x"    ][10.0 ] = ["/ARGUS_1988_I251097/d09-x01-y01"]
analyses["IdentifiedParticle"][3312]["x"    ][10.52] = ["/BELLE_2017_I1606201/d05-x01-y01"]
analyses["IdentifiedParticle"][3312]["x"    ][34.4 ] = ["/TASSO_1983_I192072/d02-x01-y01"]
analyses["IdentifiedParticle"][3312]["x"    ][34.8 ] = ["/TASSO_1989_I266893/d23-x01-y01"]
analyses["IdentifiedParticle"][3312]["p"    ][34.8 ] = ["/TASSO_1989_I266893/d18-x01-y01"]
analyses["IdentifiedParticle"][3312]["x"    ][91.2 ] = ["/OPAL_1997_S3396100/d03-x01-y01","/DELPHI_1995_S3137023/d02-x01-y01",
                                                        "/ALEPH_1996_S3486095/d34-x01-y01"]
analyses["IdentifiedParticle"][3312]["xi"   ][91.2 ] = ["/OPAL_1997_S3396100/d04-x01-y01","/DELPHI_2006_I719387/d01-x03-y01",]
analyses["IdentifiedParticle"][3312]["Other"][34.8 ] = ["/TASSO_1989_I266893/d19-x01-y01","/TASSO_1989_I266893/d20-x01-y01",
                                                        "/TASSO_1989_I266893/d21-x01-y01","/TASSO_1989_I266893/d22-x01-y01"]
# xi*0
analyses["IdentifiedParticle"][3324]["x"    ][10.52] = ["/BELLE_2017_I1606201/d07-x01-y01"]
analyses["IdentifiedParticle"][3324]["x"    ][91.2 ] = ["/OPAL_1997_S3396100/d09-x01-y01","/ALEPH_1996_S3486095/d36-x01-y01"]
analyses["IdentifiedParticle"][3324]["xi"   ][91.2 ] = ["/OPAL_1997_S3396100/d10-x01-y01"]
# omega
analyses["IdentifiedParticle"][3334]["x"    ][10.52] = ["/BELLE_2017_I1606201/d06-x01-y01"]
# lambda 1520
analyses["IdentifiedParticle"][3124]["x"    ][10.0 ] = ["/ARGUS_1989_I262415/d04-x01-y01"]
analyses["IdentifiedParticle"][3124]["x"    ][10.52] = ["/BELLE_2017_I1606201/d04-x01-y01"]
analyses["IdentifiedParticle"][3124]["x"    ][91.2 ] = ["/OPAL_1997_S3396100/d11-x01-y01","/DELPHI_2000_I524694/d03-x01-y01"]
analyses["IdentifiedParticle"][3124]["xi"   ][91.2 ] = ["/OPAL_1997_S3396100/d12-x01-y01"]
# charm baryons
# lambda_c
analyses["IdentifiedParticle"][4122 ]["x"    ][10.52] = ["/BELLE_2017_I1606201/d08-x01-y01"]
analyses["IdentifiedParticle"][4122 ]["x"    ][10.54] = ["/BABAR_2007_S6895344/d01-x01-y01"]
analyses["IdentifiedParticle"][4122 ]["Other"][10.5 ] = ["/CLEO_2001_I552541/d03-x01-y01","/CLEO_2001_I552541/d03-x01-y02",
                                                         "/CLEO_2001_I552541/d03-x01-y03","/CLEO_2001_I552541/d03-x01-y04",
                                                         "/CLEO_2001_I552541/d04-x01-y01","/CLEO_2001_I552541/d04-x01-y02",
                                                         "/CLEO_2001_I552541/d04-x01-y03","/CLEO_2001_I552541/d04-x01-y04",]
# sigma_c0
analyses["IdentifiedParticle"][4112 ]["x"    ][10.52] = ["/BELLE_2017_I1606201/d11-x01-y01"]
# sigma_c*0
analyses["IdentifiedParticle"][4114 ]["x"    ][10.52] = ["/BELLE_2017_I1606201/d12-x01-y01"]
# xi_c
analyses["IdentifiedParticle"][4132 ]["x"    ][10.52] = ["/BELLE_2017_I1606201/d14-x01-y01","/BELLE_2017_I1606201/d15-x01-y01"]
# omega_c
analyses["IdentifiedParticle"][4332 ]["x"    ][10.52] = ["/BELLE_2017_I1606201/d13-x01-y01"]
# lambda_c(2595)
analyses["IdentifiedParticle"][14122]["x"    ][10.52] = ["/BELLE_2017_I1606201/d09-x01-y01"]
# lambda_c(2625)
analyses["IdentifiedParticle"][4124 ]["x"    ][10.52] = ["/BELLE_2017_I1606201/d10-x01-y01"]
# b fragmentation
analyses["IdentifiedParticle"][511]["weak"     ] = ["/DELPHI_2011_I890503/d01-x01-y01","/SLD_2002_S4869273/d01-x01-y01",
                                                    "/ALEPH_2001_S4656318/d01-x01-y01","/OPAL_2003_I599181/d01-x01-y01"]
analyses["IdentifiedParticle"][511]["weak_mean"] = ["/DELPHI_2011_I890503/d02-x01-y01","/ALEPH_2001_S4656318/d07-x01-y01",
                                                    "/OPAL_2003_I599181/d02-x01-y01"]
analyses["IdentifiedParticle"][511]["lead"     ] = ["/ALEPH_2001_S4656318/d01-x01-y02"]
analyses["IdentifiedParticle"][511]["lead_mean"] = ["/ALEPH_2001_S4656318/d07-x01-y02"]
# multiplcities
# mesons
analyses["Multiplicity"][211 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d01-x01-y01"]
analyses["Multiplicity"][211 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d01-x01-y02"]
analyses["Multiplicity"][211 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d01-x01-y03","/DELPHI_1996_S3430090/d36-x01-y01",
                                         "/DELPHI_1998_I473409/d01-x01-y02","/SLD_2004_S5693039/d02-x02-y02",
                                         "/SLD_1999_S3743934/d24-x01-y01"]
analyses["Multiplicity"][211 ][165.0] = ["/PDG_HADRON_MULTIPLICITIES/d01-x01-y04"] 
analyses["Multiplicity"][111 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d02-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d02-x01-y01","/ARGUS_1990_I278933/d01-x01-y01"]
analyses["Multiplicity"][111 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d02-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d02-x01-y02"]
analyses["Multiplicity"][111 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d02-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d02-x01-y03",
                                         "/DELPHI_1996_S3430090/d36-x01-y02","/ALEPH_1996_S3486095/d44-x01-y02"]

analyses["Multiplicity"][321 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d03-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d03-x01-y01"]
analyses["Multiplicity"][321 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d03-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d03-x01-y02"]
analyses["Multiplicity"][321 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d03-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d03-x01-y03",
                                         "/DELPHI_1996_S3430090/d36-x01-y03","/DELPHI_1998_I473409/d01-x01-y03",
                                         "/SLD_2004_S5693039/d03-x02-y02","/SLD_1999_S3743934/d24-x02-y01"]
analyses["Multiplicity"][321 ][165.0] = ["/PDG_HADRON_MULTIPLICITIES/d03-x01-y04","/PDG_HADRON_MULTIPLICITIES_RATIOS/d03-x01-y04"]
analyses["Multiplicity"][311 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d04-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d04-x01-y01","/PLUTO_1981_I165122/d02-x01-y01"]
analyses["Multiplicity"][311 ][30.  ] = ["/TASSO_1990_I284251/d04-x01-y01"]
analyses["Multiplicity"][311 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d04-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d04-x01-y02"]
analyses["Multiplicity"][311 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d04-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d04-x01-y03",
                                         "/DELPHI_1996_S3430090/d36-x01-y04","/ALEPH_1996_S3486095/d44-x01-y05","/SLD_1999_S3743934/d24-x03-y01"]
analyses["Multiplicity"][311 ][165.0] = ["/PDG_HADRON_MULTIPLICITIES/d04-x01-y04","/PDG_HADRON_MULTIPLICITIES_RATIOS/d04-x01-y04"]
analyses["Multiplicity"][221 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d05-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d05-x01-y01","/ARGUS_1990_I278933/d02-x01-y01"]
analyses["Multiplicity"][221 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d05-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d05-x01-y02"]
analyses["Multiplicity"][221 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d05-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d05-x01-y03",
                                         "/DELPHI_1996_S3430090/d36-x01-y05","/ALEPH_1996_S3486095/d44-x01-y03"]
analyses["Multiplicity"][331 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d06-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d06-x01-y01"]
analyses["Multiplicity"][331 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d06-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d06-x01-y02"]
analyses["Multiplicity"][331 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d06-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d06-x01-y03",
                                         "/DELPHI_1996_S3430090/d36-x01-y06","/ALEPH_1996_S3486095/d44-x01-y04"]
analyses["Multiplicity"][411 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d07-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d07-x01-y01"]
analyses["Multiplicity"][411 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d07-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d07-x01-y02"]
analyses["Multiplicity"][411 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d07-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d07-x01-y03","/DELPHI_1996_S3430090/d36-x01-y07"]
analyses["Multiplicity"][421 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d08-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d08-x01-y01"]
analyses["Multiplicity"][421 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d08-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d08-x01-y02"]
analyses["Multiplicity"][421 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d08-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d08-x01-y03","/DELPHI_1996_S3430090/d36-x01-y08"]
analyses["Multiplicity"][431 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d09-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d09-x01-y01"]
analyses["Multiplicity"][431 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d09-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d09-x01-y02"]
analyses["Multiplicity"][431 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d09-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d09-x01-y03"]
analyses["Multiplicity"][511 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d10-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d10-x01-y01","/DELPHI_1996_S3430090/d36-x01-y09"]
analyses["Multiplicity"][521 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d11-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d11-x01-y01"]
analyses["Multiplicity"][531 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d12-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d12-x01-y01"]
analyses["Multiplicity"][9010221][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d13-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d13-x01-y01"]
analyses["Multiplicity"][9010221][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d13-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d13-x01-y02"]
analyses["Multiplicity"][9010221][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d13-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d13-x01-y03","/DELPHI_1996_S3430090/d37-x01-y01"]
analyses["Multiplicity"][9000211][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d14-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d14-x01-y01"]
analyses["Multiplicity"][113 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d15-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d15-x01-y01","/ARGUS_1993_S2789213/d01-x01-y02"]
analyses["Multiplicity"][113 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d15-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d15-x01-y02"]
analyses["Multiplicity"][113 ][34.0 ] = ["/TASSO_1982_I179022/d01-x01-y01"]
analyses["Multiplicity"][113 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d15-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d15-x01-y03",
                                         "/ALEPH_1996_S3486095/d44-x01-y06","/DELPHI_1996_S3430090/d38-x01-y01"]
analyses["Multiplicity"][213 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d16-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d16-x01-y01"]
analyses["Multiplicity"][223 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d17-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d17-x01-y01","/ARGUS_1993_S2789213/d01-x01-y01"]
analyses["Multiplicity"][223 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d17-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d17-x01-y02",
                                         "/ALEPH_1996_S3486095/d44-x01-y07"]
analyses["Multiplicity"][323 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d18-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d18-x01-y01","/ARGUS_1993_S2789213/d01-x01-y04"]
analyses["Multiplicity"][323 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d18-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d18-x01-y02"]
analyses["Multiplicity"][323 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d18-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d18-x01-y03",
                                         "/OPAL_1993_I342766/d02-x01-y01","/DELPHI_1996_S3430090/d38-x01-y02","/ALEPH_1996_S3486095/d44-x01-y09"]
analyses["Multiplicity"][313 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d19-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d19-x01-y01","/ARGUS_1993_S2789213/d01-x01-y03"]
analyses["Multiplicity"][313 ][30.  ] = ["/TASSO_1990_I284251/d09-x01-y01"]
analyses["Multiplicity"][313 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d19-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d19-x01-y02"]
analyses["Multiplicity"][313 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d19-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d19-x01-y03",
                                         "/DELPHI_1996_S3430090/d38-x01-y03","/ALEPH_1996_S3486095/d44-x01-y10","/SLD_1999_S3743934/d24-x04-y01"]
analyses["Multiplicity"][333 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d20-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d20-x01-y01","/ARGUS_1993_S2789213/d01-x01-y05"]
analyses["Multiplicity"][333 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d20-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d20-x01-y02"]
analyses["Multiplicity"][333 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d20-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d20-x01-y03",
                                         "/DELPHI_1996_S3430090/d38-x01-y04","/ALEPH_1996_S3486095/d44-x01-y08","/SLD_1999_S3743934/d24-x05-y01"]
analyses["Multiplicity"][413 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d21-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d21-x01-y01"]
analyses["Multiplicity"][413 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d21-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d21-x01-y02"]
analyses["Multiplicity"][413 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d21-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d21-x01-y03",
                                         "/DELPHI_1996_S3430090/d38-x01-y05","/OPAL_1995_I382219/d03-x01-y01"]
analyses["Multiplicity"][423 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d22-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d22-x01-y01"]
analyses["Multiplicity"][423 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d22-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d22-x01-y02"]
analyses["Multiplicity"][433 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d23-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d23-x01-y01"]
analyses["Multiplicity"][433 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d23-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d23-x01-y02"]
analyses["Multiplicity"][513 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d24-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d24-x01-y01"]
analyses["Multiplicity"][443    ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d25-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d25-x01-y01"]
analyses["Multiplicity"][443    ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d25-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d25-x01-y02",
                                            "/OPAL_1996_S3257789/d02-x01-y01"]
analyses["Multiplicity"][100443 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d26-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d26-x01-y01",
                                            "/OPAL_1996_S3257789/d02-x01-y02"]
analyses["Multiplicity"][553    ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d27-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d27-x01-y01"]
analyses["Multiplicity"][20223  ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d28-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d28-x01-y01"]
analyses["Multiplicity"][20333  ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d29-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d29-x01-y01"]
analyses["Multiplicity"][20443  ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d30-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d30-x01-y01"]
analyses["Multiplicity"][225 ][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d31-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d31-x01-y01"]
analyses["Multiplicity"][225 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d31-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d31-x01-y02"]
analyses["Multiplicity"][225 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d31-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d31-x01-y03","/DELPHI_1996_S3430090/d39-x01-y01"]
analyses["Multiplicity"][335 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d32-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d32-x01-y01"]
analyses["Multiplicity"][325 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d33-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d33-x01-y01"]
analyses["Multiplicity"][315 ][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d34-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d34-x01-y01"]
analyses["Multiplicity"][315 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d34-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d34-x01-y02","/DELPHI_1996_S3430090/d39-x01-y02"]
analyses["Multiplicity"][515 ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d35-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d35-x01-y01"]
analyses["Multiplicity"][20431][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d36-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d36-x01-y01"]
analyses["Multiplicity"][435  ][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d37-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d37-x01-y01"]
#baryons
analyses["Multiplicity"][2212][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d38-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d38-x01-y01"]
analyses["Multiplicity"][2212][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d38-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d38-x01-y02"]
analyses["Multiplicity"][2212][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d38-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d38-x01-y03",
                                         "/DELPHI_1996_S3430090/d40-x01-y01","/DELPHI_1998_I473409/d01-x01-y04",
                                         "/SLD_2004_S5693039/d04-x02-y02","/SLD_1999_S3743934/d24-x06-y01"]
analyses["Multiplicity"][2212][165.0] = ["/PDG_HADRON_MULTIPLICITIES/d38-x01-y04","/PDG_HADRON_MULTIPLICITIES_RATIOS/d38-x01-y04"]
analyses["Multiplicity"][3122][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d39-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d39-x01-y01","/ARGUS_1988_I251097/d02-x01-y01"]
analyses["Multiplicity"][3122][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d39-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d39-x01-y02"]
analyses["Multiplicity"][3122][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d39-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d39-x01-y03",
                                         "/DELPHI_1996_S3430090/d40-x01-y02","/ALEPH_1996_S3486095/d44-x01-y11","/DELPHI_1993_I360638/d02-x01-y01",
                                         "/SLD_1999_S3743934/d24-x07-y01"]
analyses["Multiplicity"][3122][165.0] = ["/PDG_HADRON_MULTIPLICITIES/d39-x01-y04","/PDG_HADRON_MULTIPLICITIES_RATIOS/d39-x01-y04"]
analyses["Multiplicity"][3212][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d40-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d40-x01-y01","/ARGUS_1988_I251097/d02-x01-y03"]
analyses["Multiplicity"][3212][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d40-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d40-x01-y02"]
analyses["Multiplicity"][3112][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d41-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d41-x01-y01"]
analyses["Multiplicity"][3222][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d42-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d42-x01-y01","/ALEPH_1996_S3486095/d44-x01-y12"]
analyses["Multiplicity"][3312][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d44-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d44-x01-y01","/ARGUS_1988_I251097/d02-x01-y02"]
analyses["Multiplicity"][3312][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d44-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d44-x01-y02"]
analyses["Multiplicity"][3312][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d44-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d44-x01-y03",
                                         "/DELPHI_1996_S3430090/d40-x01-y03","/ALEPH_1996_S3486095/d44-x01-y13"]
analyses["Multiplicity"][2224][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d45-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d45-x01-y01"]
analyses["Multiplicity"][2224][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d45-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d45-x01-y02","/DELPHI_1996_S3430090/d40-x01-y05"]
analyses["Multiplicity"][3114][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d46-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d46-x01-y01","/ARGUS_1988_I251097/d02-x01-y04"]
analyses["Multiplicity"][3114][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d46-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d46-x01-y02"]
analyses["Multiplicity"][3114][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d46-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d46-x01-y03"]
analyses["Multiplicity"][3224][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d47-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d47-x01-y01","/ARGUS_1988_I251097/d02-x01-y05"]
analyses["Multiplicity"][3224][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d47-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d47-x01-y02"]
analyses["Multiplicity"][3224][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d47-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d47-x01-y03"]
analyses["Multiplicity"][3324][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d49-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d49-x01-y01","/ARGUS_1988_I251097/d02-x01-y06"]
analyses["Multiplicity"][3324][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d49-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d49-x01-y02",
                                         "/DELPHI_1996_S3430090/d40-x01-y07","/ALEPH_1996_S3486095/d44-x01-y15"]
analyses["Multiplicity"][3334][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d50-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d50-x01-y01","/ARGUS_1988_I251097/d02-x01-y07"]
analyses["Multiplicity"][3334][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d50-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d50-x01-y02"]
analyses["Multiplicity"][3334][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d50-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d50-x01-y03",
                                         "/DELPHI_1996_S3430090/d40-x01-y04","/ALEPH_1996_S3486095/d44-x01-y16"]
analyses["Multiplicity"][4122][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d51-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d51-x01-y01",
                                         "/BABAR_2007_S6895344/d02-x01-y01"]
analyses["Multiplicity"][4122][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d51-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d51-x01-y02"]
analyses["Multiplicity"][4122][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d51-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d51-x01-y03"]
analyses["Multiplicity"][5122][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d52-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d52-x01-y01","/DELPHI_1996_S3430090/d40-x01-y08"]
analyses["Multiplicity"][4222][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d53-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d53-x01-y01"]
analyses["Multiplicity"][3124][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d54-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d54-x01-y01"]
analyses["Multiplicity"][3124][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d54-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d54-x01-y02"]
#
analyses["Multiplicity"]["3222B"][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d43-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d43-x01-y01"]
analyses["Multiplicity"]["3224B"][10.  ] = ["/PDG_HADRON_MULTIPLICITIES/d48-x01-y01","/PDG_HADRON_MULTIPLICITIES_RATIOS/d48-x01-y01"]
analyses["Multiplicity"]["3224B"][32.0 ] = ["/PDG_HADRON_MULTIPLICITIES/d48-x01-y02","/PDG_HADRON_MULTIPLICITIES_RATIOS/d48-x01-y02"]
analyses["Multiplicity"]["3224B"][91.2 ] = ["/PDG_HADRON_MULTIPLICITIES/d48-x01-y03","/PDG_HADRON_MULTIPLICITIES_RATIOS/d48-x01-y03",
                                            "/DELPHI_1996_S3430090/d40-x01-y06","/ALEPH_1996_S3486095/d44-x01-y14"]
# event shapes
# thrust based
analyses["EventShapes"]["T"][14.0 ] = ["/TASSO_1990_S2148048/d08-x01-y01"]
analyses["EventShapes"]["T"][22.0 ] = ["/TASSO_1990_S2148048/d08-x01-y02"]
analyses["EventShapes"]["T"][29.0 ] = ["/HRS_1985_I201482/d03-x01-y01","/HRS_1985_I201482/d04-x01-y01"]
analyses["EventShapes"]["T"][35.0 ] = ["/TASSO_1990_S2148048/d08-x01-y03","/TASSO_1988_I263859/d03-x01-y01",
                                       "/JADE_1998_S3612880/d06-x01-y01"]
analyses["EventShapes"]["T"][44.0 ] = ["/TASSO_1990_S2148048/d08-x01-y04","/JADE_1998_S3612880/d02-x01-y01"]
analyses["EventShapes"]["T"][45.0 ] = ["/DELPHI_2003_I620250/d01-x01-y01"]
analyses["EventShapes"]["T"][55.2 ] = ["/AMY_1990_I283337/d12-x01-y01"]
analyses["EventShapes"]["T"][58.0 ] = ["/TOPAZ_1993_I361661/d01-x01-y01"]
analyses["EventShapes"]["T"][66.0 ] = ["/DELPHI_2003_I620250/d01-x01-y02"]
analyses["EventShapes"]["T"][76.0 ] = ["/DELPHI_2003_I620250/d01-x01-y03"]
analyses["EventShapes"]["T"][91.2 ] = ["/DELPHI_1996_S3430090/d11-x01-y01","/ALEPH_1996_S3486095/d03-x01-y01",
                                       "/OPAL_2004_S6132243/d01-x01-y01","/ALEPH_2004_S5765862/d54-x01-y01"]
analyses["EventShapes"]["T"][133.0] = ["/ALEPH_2004_S5765862/d55-x01-y01","/OPAL_2004_S6132243/d01-x01-y02"]
analyses["EventShapes"]["T"][161.0] = ["/ALEPH_2004_S5765862/d56-x01-y01"]
analyses["EventShapes"]["T"][172.0] = ["/ALEPH_2004_S5765862/d57-x01-y01"]
analyses["EventShapes"]["T"][177.0] = ["/OPAL_2004_S6132243/d01-x01-y03"]
analyses["EventShapes"]["T"][183.0] = ["/DELPHI_2003_I620250/d38-x01-y01","/ALEPH_2004_S5765862/d58-x01-y01"]
analyses["EventShapes"]["T"][189.0] = ["/DELPHI_2003_I620250/d38-x01-y02","/ALEPH_2004_S5765862/d59-x01-y01"]
analyses["EventShapes"]["T"][192.0] = ["/DELPHI_2003_I620250/d38-x01-y03"]
analyses["EventShapes"]["T"][196.0] = ["/DELPHI_2003_I620250/d38-x01-y04"]
analyses["EventShapes"]["T"][200.0] = ["/DELPHI_2003_I620250/d39-x01-y01","/ALEPH_2004_S5765862/d60-x01-y01"]
analyses["EventShapes"]["T"][197.0] = ["/OPAL_2004_S6132243/d01-x01-y04"]
analyses["EventShapes"]["T"][202.0] = ["/DELPHI_2003_I620250/d39-x01-y02"]
analyses["EventShapes"]["T"][205.0] = ["/DELPHI_2003_I620250/d39-x01-y03"]
analyses["EventShapes"]["T"][206.0] = ["/ALEPH_2004_S5765862/d61-x01-y01"]
analyses["EventShapes"]["T"][207.0] = ["/DELPHI_2003_I620250/d39-x01-y04"]
analyses["EventShapes"]["T"][41.4 ] = ["/L3_2004_I652683/d21-x01-y01"]
analyses["EventShapes"]["T"][55.3 ] = ["/L3_2004_I652683/d21-x01-y02"]
analyses["EventShapes"]["T"][65.4 ] = ["/L3_2004_I652683/d21-x01-y03"]
analyses["EventShapes"]["T"][75.7 ] = ["/L3_2004_I652683/d22-x01-y01"]
analyses["EventShapes"]["T"][82.3 ] = ["/L3_2004_I652683/d22-x01-y02"]
analyses["EventShapes"]["T"][85.1 ] = ["/L3_2004_I652683/d22-x01-y03"]
analyses["EventShapes"]["T"][130.1] = ["/L3_2004_I652683/d23-x01-y01"]
analyses["EventShapes"]["T"][136.3] = ["/L3_2004_I652683/d23-x01-y02"]
analyses["EventShapes"]["T"][161.3] = ["/L3_2004_I652683/d23-x01-y03"]
analyses["EventShapes"]["T"][172.3] = ["/L3_2004_I652683/d24-x01-y01"]
analyses["EventShapes"]["T"][182.8] = ["/L3_2004_I652683/d24-x01-y02"]
analyses["EventShapes"]["T"][188.6] = ["/L3_2004_I652683/d24-x01-y03"]
analyses["EventShapes"]["T"][194.4] = ["/L3_2004_I652683/d25-x01-y01"]
analyses["EventShapes"]["T"][200.2] = ["/L3_2004_I652683/d25-x01-y02"]
analyses["EventShapes"]["T"][206.2] = ["/L3_2004_I652683/d25-x01-y03"]
analyses["EventShapes"]["Moment_T"][91.2 ] = ["/OPAL_2004_S6132243/d15-x01-y01"]
analyses["EventShapes"]["Moment_T"][133.0] = ["/OPAL_2004_S6132243/d15-x01-y02"]
analyses["EventShapes"]["Moment_T"][177.0] = ["/OPAL_2004_S6132243/d15-x01-y03"]
analyses["EventShapes"]["Moment_T"][197.0] = ["/OPAL_2004_S6132243/d15-x01-y04"]

analyses["EventShapesFlavour"]["T"][2][91.2] = ["/L3_2004_I652683/d47-x01-y01"]
analyses["EventShapesFlavour"]["T"][5][91.2] = ["/L3_2004_I652683/d47-x01-y02"]
analyses["EventShapesFlavour"]["HeavyJetMass"][2][91.2] = ["/L3_2004_I652683/d48-x01-y01"]
analyses["EventShapesFlavour"]["HeavyJetMass"][5][91.2] = ["/L3_2004_I652683/d48-x01-y02"]
analyses["EventShapesFlavour"]["BT"][2][91.2] = ["/L3_2004_I652683/d49-x01-y01"]
analyses["EventShapesFlavour"]["BT"][5][91.2] = ["/L3_2004_I652683/d49-x01-y02"]
analyses["EventShapesFlavour"]["BW"][2][91.2] = ["/L3_2004_I652683/d50-x01-y01"]
analyses["EventShapesFlavour"]["BW"][5][91.2] = ["/L3_2004_I652683/d50-x01-y02"]
analyses["EventShapesFlavour"]["C"][2][91.2] = ["/L3_2004_I652683/d51-x01-y01"]
analyses["EventShapesFlavour"]["C"][5][91.2] = ["/L3_2004_I652683/d51-x01-y02"]
analyses["EventShapesFlavour"]["D"][2][91.2] = ["/L3_2004_I652683/d52-x01-y01"]
analyses["EventShapesFlavour"]["D"][5][91.2] = ["/L3_2004_I652683/d52-x01-y02"]

analyses["EventShapes"]["Moment_H" ][91.2 ] = ["/OPAL_2004_S6132243/d16-x01-y01"]
analyses["EventShapes"]["Moment_H" ][133.0] = ["/OPAL_2004_S6132243/d16-x01-y02"]
analyses["EventShapes"]["Moment_H" ][177.0] = ["/OPAL_2004_S6132243/d16-x01-y03"]
analyses["EventShapes"]["Moment_H" ][197.0] = ["/OPAL_2004_S6132243/d16-x01-y04"]
analyses["EventShapes"]["Moment_C" ][91.2 ] = ["/OPAL_2004_S6132243/d17-x01-y01"]
analyses["EventShapes"]["Moment_C" ][133.0] = ["/OPAL_2004_S6132243/d17-x01-y02"]
analyses["EventShapes"]["Moment_C" ][177.0] = ["/OPAL_2004_S6132243/d17-x01-y03"]
analyses["EventShapes"]["Moment_C" ][197.0] = ["/OPAL_2004_S6132243/d17-x01-y04"]
analyses["EventShapes"]["Moment_BT"][91.2 ] = ["/OPAL_2004_S6132243/d18-x01-y01"]
analyses["EventShapes"]["Moment_BT"][133.0] = ["/OPAL_2004_S6132243/d18-x01-y02"]
analyses["EventShapes"]["Moment_BT"][177.0] = ["/OPAL_2004_S6132243/d18-x01-y03"]
analyses["EventShapes"]["Moment_BT"][197.0] = ["/OPAL_2004_S6132243/d18-x01-y04"]
analyses["EventShapes"]["Moment_BW"][91.2 ] = ["/OPAL_2004_S6132243/d19-x01-y01"]
analyses["EventShapes"]["Moment_BW"][133.0] = ["/OPAL_2004_S6132243/d19-x01-y02"]
analyses["EventShapes"]["Moment_BW"][177.0] = ["/OPAL_2004_S6132243/d19-x01-y03"]
analyses["EventShapes"]["Moment_BW"][197.0] = ["/OPAL_2004_S6132243/d19-x01-y04"]
analyses["EventShapes"]["Moment_y" ][91.2 ] = ["/OPAL_2004_S6132243/d20-x01-y01"]
analyses["EventShapes"]["Moment_y" ][133.0] = ["/OPAL_2004_S6132243/d20-x01-y02"]
analyses["EventShapes"]["Moment_y" ][177.0] = ["/OPAL_2004_S6132243/d20-x01-y03"]
analyses["EventShapes"]["Moment_y" ][197.0] = ["/OPAL_2004_S6132243/d20-x01-y04"]
analyses["EventShapes"]["Moment_M" ][91.2 ] = ["/OPAL_2004_S6132243/d21-x01-y01"]
analyses["EventShapes"]["Moment_M" ][133.0] = ["/OPAL_2004_S6132243/d21-x01-y02"]
analyses["EventShapes"]["Moment_M" ][177.0] = ["/OPAL_2004_S6132243/d21-x01-y03"]
analyses["EventShapes"]["Moment_M" ][197.0] = ["/OPAL_2004_S6132243/d21-x01-y04"]
analyses["EventShapes"]["Moment_m" ][91.2 ] = ["/OPAL_2004_S6132243/d22-x01-y01"]
analyses["EventShapes"]["Moment_m" ][133.0] = ["/OPAL_2004_S6132243/d22-x01-y02"]
analyses["EventShapes"]["Moment_m" ][177.0] = ["/OPAL_2004_S6132243/d22-x01-y03"]
analyses["EventShapes"]["Moment_m" ][197.0] = ["/OPAL_2004_S6132243/d22-x01-y04"]
analyses["EventShapes"]["Moment_S" ][91.2 ] = ["/OPAL_2004_S6132243/d23-x01-y01"]
analyses["EventShapes"]["Moment_S" ][133.0] = ["/OPAL_2004_S6132243/d23-x01-y02"]
analyses["EventShapes"]["Moment_S" ][177.0] = ["/OPAL_2004_S6132243/d23-x01-y03"]
analyses["EventShapes"]["Moment_S" ][197.0] = ["/OPAL_2004_S6132243/d23-x01-y04"]
analyses["EventShapes"]["Moment_O" ][91.2 ] = ["/OPAL_2004_S6132243/d24-x01-y01"]
analyses["EventShapes"]["Moment_O" ][133.0] = ["/OPAL_2004_S6132243/d24-x01-y02"]
analyses["EventShapes"]["Moment_O" ][177.0] = ["/OPAL_2004_S6132243/d24-x01-y03"]
analyses["EventShapes"]["Moment_O" ][197.0] = ["/OPAL_2004_S6132243/d24-x01-y04"]
analyses["EventShapes"]["Moment_L" ][91.2 ] = ["/OPAL_2004_S6132243/d25-x01-y01"]
analyses["EventShapes"]["Moment_L" ][133.0] = ["/OPAL_2004_S6132243/d25-x01-y02"]
analyses["EventShapes"]["Moment_L" ][177.0] = ["/OPAL_2004_S6132243/d25-x01-y03"]
analyses["EventShapes"]["Moment_L" ][197.0] = ["/OPAL_2004_S6132243/d25-x01-y04"]
analyses["EventShapes"]["Moment_BN"][91.2 ] = ["/OPAL_2004_S6132243/d26-x01-y01"]
analyses["EventShapes"]["Moment_BN"][133.0] = ["/OPAL_2004_S6132243/d26-x01-y02"]
analyses["EventShapes"]["Moment_BN"][177.0] = ["/OPAL_2004_S6132243/d26-x01-y03"]
analyses["EventShapes"]["Moment_BN"][197.0] = ["/OPAL_2004_S6132243/d26-x01-y04"]

analyses["EventShapes"]["Major"][45.0 ] = ["/DELPHI_2003_I620250/d02-x01-y01"]
analyses["EventShapes"]["Major"][55.2 ] = ["/AMY_1990_I283337/d13-x01-y01"]
analyses["EventShapes"]["Major"][66.0 ] = ["/DELPHI_2003_I620250/d02-x01-y02"]
analyses["EventShapes"]["Major"][76.0 ] = ["/DELPHI_2003_I620250/d02-x01-y03"]
analyses["EventShapes"]["Major"][91.2 ] = ["/DELPHI_1996_S3430090/d12-x01-y01","/OPAL_2004_S6132243/d07-x01-y01",
                                           "/ALEPH_2004_S5765862/d94-x01-y01"]
analyses["EventShapes"]["Major"][133.0] = ["/ALEPH_2004_S5765862/d95-x01-y01","/OPAL_2004_S6132243/d07-x01-y02"]
analyses["EventShapes"]["Major"][161.0] = ["/ALEPH_2004_S5765862/d96-x01-y01"]
analyses["EventShapes"]["Major"][172.0] = ["/ALEPH_2004_S5765862/d97-x01-y01"]
analyses["EventShapes"]["Major"][177.0] = ["/OPAL_2004_S6132243/d07-x01-y03"]
analyses["EventShapes"]["Major"][183.0] = ["/DELPHI_2003_I620250/d40-x01-y01","/ALEPH_2004_S5765862/d98-x01-y01"]
analyses["EventShapes"]["Major"][189.0] = ["/DELPHI_2003_I620250/d40-x01-y02","/ALEPH_2004_S5765862/d99-x01-y01"]
analyses["EventShapes"]["Major"][192.0] = ["/DELPHI_2003_I620250/d40-x01-y03"]
analyses["EventShapes"]["Major"][196.0] = ["/DELPHI_2003_I620250/d40-x01-y04"]
analyses["EventShapes"]["Major"][197.0] = ["/OPAL_2004_S6132243/d07-x01-y04"]
analyses["EventShapes"]["Major"][200.0] = ["/DELPHI_2003_I620250/d41-x01-y01","/ALEPH_2004_S5765862/d100-x01-y01"]
analyses["EventShapes"]["Major"][202.0] = ["/DELPHI_2003_I620250/d41-x01-y02"]
analyses["EventShapes"]["Major"][205.0] = ["/DELPHI_2003_I620250/d41-x01-y03"]
analyses["EventShapes"]["Major"][206.0] = ["/ALEPH_2004_S5765862/d101-x01-y01"]
analyses["EventShapes"]["Major"][207.0] = ["/DELPHI_2003_I620250/d41-x01-y04"]
analyses["EventShapes"]["Minor"][45.0 ] = ["/DELPHI_2003_I620250/d03-x01-y01"]
analyses["EventShapes"]["Minor"][55.2 ] = ["/AMY_1990_I283337/d14-x01-y01"]
analyses["EventShapes"]["Minor"][66.0 ] = ["/DELPHI_2003_I620250/d03-x01-y02"]
analyses["EventShapes"]["Minor"][76.0 ] = ["/DELPHI_2003_I620250/d03-x01-y03"]
analyses["EventShapes"]["Minor"][91.2 ] = ["/DELPHI_1996_S3430090/d13-x01-y01","/ALEPH_2004_S5765862/d102-x01-y01",
                                           "/ALEPH_1996_S3486095/d04-x01-y01","/OPAL_2004_S6132243/d08-x01-y01"]
analyses["EventShapes"]["Minor"][133.0] = ["/ALEPH_2004_S5765862/d103-x01-y01","/OPAL_2004_S6132243/d08-x01-y02"]
analyses["EventShapes"]["Minor"][161.0] = ["/ALEPH_2004_S5765862/d104-x01-y01"]
analyses["EventShapes"]["Minor"][172.0] = ["/ALEPH_2004_S5765862/d105-x01-y01"]
analyses["EventShapes"]["Minor"][177.0] = ["/OPAL_2004_S6132243/d08-x01-y03"]
analyses["EventShapes"]["Minor"][183.0] = ["/DELPHI_2003_I620250/d42-x01-y01","/ALEPH_2004_S5765862/d106-x01-y01"]
analyses["EventShapes"]["Minor"][189.0] = ["/DELPHI_2003_I620250/d42-x01-y02","/ALEPH_2004_S5765862/d107-x01-y01"]
analyses["EventShapes"]["Minor"][192.0] = ["/DELPHI_2003_I620250/d42-x01-y03"]
analyses["EventShapes"]["Minor"][196.0] = ["/DELPHI_2003_I620250/d42-x01-y04"]
analyses["EventShapes"]["Minor"][197.0] = ["/OPAL_2004_S6132243/d08-x01-y04"]
analyses["EventShapes"]["Minor"][200.0] = ["/DELPHI_2003_I620250/d43-x01-y01","/ALEPH_2004_S5765862/d108-x01-y01"]
analyses["EventShapes"]["Minor"][202.0] = ["/DELPHI_2003_I620250/d43-x01-y02"]
analyses["EventShapes"]["Minor"][205.0] = ["/DELPHI_2003_I620250/d43-x01-y03"]
analyses["EventShapes"]["Minor"][206.0] = ["/ALEPH_2004_S5765862/d109-x01-y01"]
analyses["EventShapes"]["Minor"][207.0] = ["/DELPHI_2003_I620250/d43-x01-y04"]
# jet broadenings
analyses["EventShapes"]["BW"][45.0 ] = ["/DELPHI_2003_I620250/d13-x01-y01"]
analyses["EventShapes"]["BW"][66.0 ] = ["/DELPHI_2003_I620250/d13-x01-y02"]
analyses["EventShapes"]["BW"][76.0 ] = ["/DELPHI_2003_I620250/d13-x01-y03"]
analyses["EventShapes"]["BW"][91.2 ] = ["/DELPHI_1996_S3430090/d23-x01-y01","/ALEPH_2004_S5765862/d78-x01-y01","/OPAL_2004_S6132243/d05-x01-y01"]
analyses["EventShapes"]["BW"][133.0] = ["/OPAL_2004_S6132243/d05-x01-y02","/ALEPH_2004_S5765862/d79-x01-y01"]
analyses["EventShapes"]["BW"][161.0] = ["/ALEPH_2004_S5765862/d80-x01-y01"]
analyses["EventShapes"]["BW"][172.0] = ["/ALEPH_2004_S5765862/d81-x01-y01"]
analyses["EventShapes"]["BW"][177.0] = ["/OPAL_2004_S6132243/d05-x01-y03"]
analyses["EventShapes"]["BW"][183.0] = ["/DELPHI_2003_I620250/d46-x01-y01","/ALEPH_2004_S5765862/d82-x01-y01"]
analyses["EventShapes"]["BW"][189.0] = ["/DELPHI_2003_I620250/d46-x01-y02","/ALEPH_2004_S5765862/d83-x01-y01"]
analyses["EventShapes"]["BW"][192.0] = ["/DELPHI_2003_I620250/d46-x01-y03"]
analyses["EventShapes"]["BW"][196.0] = ["/DELPHI_2003_I620250/d46-x01-y04"]
analyses["EventShapes"]["BW"][197.0] = ["/OPAL_2004_S6132243/d05-x01-y04"]
analyses["EventShapes"]["BW"][200.0] = ["/DELPHI_2003_I620250/d47-x01-y01","/ALEPH_2004_S5765862/d84-x01-y01"]
analyses["EventShapes"]["BW"][202.0] = ["/DELPHI_2003_I620250/d47-x01-y02"]
analyses["EventShapes"]["BW"][205.0] = ["/DELPHI_2003_I620250/d47-x01-y03"]
analyses["EventShapes"]["BW"][206.0] = ["/ALEPH_2004_S5765862/d85-x01-y01"]
analyses["EventShapes"]["BW"][207.0] = ["/DELPHI_2003_I620250/d47-x01-y04"]
analyses["EventShapes"]["BW"][41.4 ] = ["/L3_2004_I652683/d36-x01-y01"]
analyses["EventShapes"]["BW"][55.3 ] = ["/L3_2004_I652683/d36-x01-y02"]
analyses["EventShapes"]["BW"][65.4 ] = ["/L3_2004_I652683/d36-x01-y03"]
analyses["EventShapes"]["BW"][75.7 ] = ["/L3_2004_I652683/d37-x01-y01"]
analyses["EventShapes"]["BW"][82.3 ] = ["/L3_2004_I652683/d37-x01-y02"]
analyses["EventShapes"]["BW"][85.1 ] = ["/L3_2004_I652683/d37-x01-y03"]
analyses["EventShapes"]["BW"][130.1] = ["/L3_2004_I652683/d38-x01-y01"]
analyses["EventShapes"]["BW"][136.3] = ["/L3_2004_I652683/d38-x01-y02"]
analyses["EventShapes"]["BW"][161.3] = ["/L3_2004_I652683/d38-x01-y03"]
analyses["EventShapes"]["BW"][172.3] = ["/L3_2004_I652683/d39-x01-y01"]
analyses["EventShapes"]["BW"][182.8] = ["/L3_2004_I652683/d39-x01-y02"]
analyses["EventShapes"]["BW"][188.6] = ["/L3_2004_I652683/d39-x01-y03"]
analyses["EventShapes"]["BW"][194.4] = ["/L3_2004_I652683/d40-x01-y01"]
analyses["EventShapes"]["BW"][200.2] = ["/L3_2004_I652683/d40-x01-y02"]
analyses["EventShapes"]["BW"][206.2] = ["/L3_2004_I652683/d40-x01-y03"]
analyses["EventShapes"]["BW"][35.0] = ["/JADE_1998_S3612880/d09-x01-y01"]
analyses["EventShapes"]["BW"][44.0] = ["/JADE_1998_S3612880/d05-x01-y01"]
analyses["EventShapes"]["BN"][45.0 ] = ["/DELPHI_2003_I620250/d14-x01-y01"]
analyses["EventShapes"]["BN"][66.0 ] = ["/DELPHI_2003_I620250/d14-x01-y02"]
analyses["EventShapes"]["BN"][76.0 ] = ["/DELPHI_2003_I620250/d14-x01-y03"]
analyses["EventShapes"]["BN"][91.2 ] = ["/DELPHI_1996_S3430090/d24-x01-y01","/OPAL_2004_S6132243/d13-x01-y01"]
analyses["EventShapes"]["BN"][133.0] = ["/OPAL_2004_S6132243/d13-x01-y02"]
analyses["EventShapes"]["BN"][177.0] = ["/OPAL_2004_S6132243/d13-x01-y03"]
analyses["EventShapes"]["BN"][197.0] = ["/OPAL_2004_S6132243/d13-x01-y04"]
analyses["EventShapes"]["BT"][41.4 ] = ["/L3_2004_I652683/d31-x01-y01"]
analyses["EventShapes"]["BT"][55.3 ] = ["/L3_2004_I652683/d31-x01-y02"]
analyses["EventShapes"]["BT"][65.4 ] = ["/L3_2004_I652683/d31-x01-y03"]
analyses["EventShapes"]["BT"][75.7 ] = ["/L3_2004_I652683/d32-x01-y01"]
analyses["EventShapes"]["BT"][82.3 ] = ["/L3_2004_I652683/d32-x01-y02"]
analyses["EventShapes"]["BT"][85.1 ] = ["/L3_2004_I652683/d32-x01-y03"]
analyses["EventShapes"]["BT"][91.2 ] = ["/DELPHI_1996_S3430090/d25-x01-y01","/OPAL_2004_S6132243/d04-x01-y01",
                                        "/ALEPH_2004_S5765862/d70-x01-y01"]
analyses["EventShapes"]["BT"][130.1] = ["/L3_2004_I652683/d33-x01-y01"]
analyses["EventShapes"]["BT"][133.0] = ["/OPAL_2004_S6132243/d04-x01-y02","/ALEPH_2004_S5765862/d71-x01-y01"]
analyses["EventShapes"]["BT"][136.3] = ["/L3_2004_I652683/d33-x01-y02"]
analyses["EventShapes"]["BT"][161.3] = ["/L3_2004_I652683/d33-x01-y03"]
analyses["EventShapes"]["BT"][172.3] = ["/L3_2004_I652683/d34-x01-y01"]
analyses["EventShapes"]["BT"][177.0] = ["/OPAL_2004_S6132243/d04-x01-y03"]
analyses["EventShapes"]["BT"][182.8] = ["/L3_2004_I652683/d34-x01-y02"]
analyses["EventShapes"]["BT"][188.6] = ["/L3_2004_I652683/d34-x01-y03"]
analyses["EventShapes"]["BT"][194.4] = ["/L3_2004_I652683/d35-x01-y01"]
analyses["EventShapes"]["BT"][197.0] = ["/OPAL_2004_S6132243/d04-x01-y04"]
analyses["EventShapes"]["BT"][200.2] = ["/L3_2004_I652683/d35-x01-y02"]
analyses["EventShapes"]["BT"][206.2] = ["/L3_2004_I652683/d35-x01-y03"]
analyses["EventShapes"]["BT"][45.0 ] = ["/DELPHI_2003_I620250/d15-x01-y01"]
analyses["EventShapes"]["BT"][66.0 ] = ["/DELPHI_2003_I620250/d15-x01-y02"]
analyses["EventShapes"]["BT"][76.0 ] = ["/DELPHI_2003_I620250/d15-x01-y03"]
analyses["EventShapes"]["BT"][161.0] = ["/ALEPH_2004_S5765862/d72-x01-y01"]
analyses["EventShapes"]["BT"][172.0] = ["/ALEPH_2004_S5765862/d73-x01-y01"]
analyses["EventShapes"]["BT"][183.0] = ["/DELPHI_2003_I620250/d48-x01-y01","/ALEPH_2004_S5765862/d74-x01-y01"]
analyses["EventShapes"]["BT"][189.0] = ["/DELPHI_2003_I620250/d48-x01-y02","/ALEPH_2004_S5765862/d75-x01-y01"]
analyses["EventShapes"]["BT"][192.0] = ["/DELPHI_2003_I620250/d48-x01-y03"]
analyses["EventShapes"]["BT"][196.0] = ["/DELPHI_2003_I620250/d48-x01-y04"]
analyses["EventShapes"]["BT"][200.0] = ["/DELPHI_2003_I620250/d49-x01-y01","/ALEPH_2004_S5765862/d76-x01-y01"]
analyses["EventShapes"]["BT"][202.0] = ["/DELPHI_2003_I620250/d49-x01-y02"]
analyses["EventShapes"]["BT"][205.0] = ["/DELPHI_2003_I620250/d49-x01-y03"]
analyses["EventShapes"]["BT"][206.0] = ["/ALEPH_2004_S5765862/d77-x01-y01"]
analyses["EventShapes"]["BT"][207.0] = ["/DELPHI_2003_I620250/d49-x01-y04"]
analyses["EventShapes"]["BT"][35.0] = ["/JADE_1998_S3612880/d08-x01-y01"]
analyses["EventShapes"]["BT"][44.0] = ["/JADE_1998_S3612880/d04-x01-y01"]
analyses["EventShapes"]["Bdiff"][45.0 ] = ["/DELPHI_2003_I620250/d16-x01-y01"]
analyses["EventShapes"]["Bdiff"][66.0 ] = ["/DELPHI_2003_I620250/d16-x01-y02"]
analyses["EventShapes"]["Bdiff"][76.0 ] = ["/DELPHI_2003_I620250/d16-x01-y03"]
analyses["EventShapes"]["Bdiff"][91.2 ] = ["/DELPHI_1996_S3430090/d26-x01-y01"]
analyses["EventShapes"]["Bdiff"][183.0] = ["/DELPHI_2003_I620250/d50-x01-y01"]
analyses["EventShapes"]["Bdiff"][189.0] = ["/DELPHI_2003_I620250/d50-x01-y02"]
analyses["EventShapes"]["Bdiff"][192.0] = ["/DELPHI_2003_I620250/d50-x01-y03"]
analyses["EventShapes"]["Bdiff"][196.0] = ["/DELPHI_2003_I620250/d50-x01-y04"]
analyses["EventShapes"]["Bdiff"][200.0] = ["/DELPHI_2003_I620250/d51-x01-y01"]
analyses["EventShapes"]["Bdiff"][202.0] = ["/DELPHI_2003_I620250/d51-x01-y02"]
analyses["EventShapes"]["Bdiff"][205.0] = ["/DELPHI_2003_I620250/d51-x01-y03"]
analyses["EventShapes"]["Bdiff"][207.0] = ["/DELPHI_2003_I620250/d51-x01-y04"]
# C and D
analyses["EventShapes"]["C"][45.0 ] = ["/DELPHI_2003_I620250/d17-x01-y01"]
analyses["EventShapes"]["C"][66.0 ] = ["/DELPHI_2003_I620250/d17-x01-y02"]
analyses["EventShapes"]["C"][76.0 ] = ["/DELPHI_2003_I620250/d17-x01-y03"]
analyses["EventShapes"]["C"][91.2 ] = ["/DELPHI_1996_S3430090/d18-x01-y01","/ALEPH_1996_S3486095/d07-x01-y01",
                                       "/ALEPH_2004_S5765862/d86-x01-y01","/OPAL_2004_S6132243/d03-x01-y01"]
analyses["EventShapes"]["C"][133.0] = ["/OPAL_2004_S6132243/d03-x01-y02","/ALEPH_2004_S5765862/d87-x01-y01"]
analyses["EventShapes"]["C"][161.0] = ["/ALEPH_2004_S5765862/d88-x01-y01"]
analyses["EventShapes"]["C"][172.0] = ["/ALEPH_2004_S5765862/d89-x01-y01"]
analyses["EventShapes"]["C"][177.0] = ["/OPAL_2004_S6132243/d03-x01-y03"]
analyses["EventShapes"]["C"][183.0] = ["/DELPHI_2003_I620250/d52-x01-y01","/ALEPH_2004_S5765862/d90-x01-y01"]
analyses["EventShapes"]["C"][189.0] = ["/DELPHI_2003_I620250/d52-x01-y02","/ALEPH_2004_S5765862/d91-x01-y01"]
analyses["EventShapes"]["C"][192.0] = ["/DELPHI_2003_I620250/d52-x01-y03"]
analyses["EventShapes"]["C"][196.0] = ["/DELPHI_2003_I620250/d52-x01-y04"]
analyses["EventShapes"]["C"][197.0] = ["/OPAL_2004_S6132243/d03-x01-y04"]
analyses["EventShapes"]["C"][200.0] = ["/DELPHI_2003_I620250/d53-x01-y01","/ALEPH_2004_S5765862/d92-x01-y01"]
analyses["EventShapes"]["C"][202.0] = ["/DELPHI_2003_I620250/d53-x01-y02"]
analyses["EventShapes"]["C"][205.0] = ["/DELPHI_2003_I620250/d53-x01-y03"]
analyses["EventShapes"]["C"][206.0] = ["/ALEPH_2004_S5765862/d93-x01-y01"]
analyses["EventShapes"]["C"][207.0] = ["/DELPHI_2003_I620250/d53-x01-y04"]
analyses["EventShapes"]["C"][130.1] = ["/L3_2004_I652683/d41-x01-y01"]
analyses["EventShapes"]["C"][136.3] = ["/L3_2004_I652683/d41-x01-y02"]
analyses["EventShapes"]["C"][161.3] = ["/L3_2004_I652683/d41-x01-y03"]
analyses["EventShapes"]["C"][172.3] = ["/L3_2004_I652683/d42-x01-y01"]
analyses["EventShapes"]["C"][182.8] = ["/L3_2004_I652683/d42-x01-y02"]
analyses["EventShapes"]["C"][188.6] = ["/L3_2004_I652683/d42-x01-y03"]
analyses["EventShapes"]["C"][194.4] = ["/L3_2004_I652683/d43-x01-y01"]
analyses["EventShapes"]["C"][200.2] = ["/L3_2004_I652683/d43-x01-y02"]
analyses["EventShapes"]["C"][206.2] = ["/L3_2004_I652683/d43-x01-y03"]

analyses["EventShapes"]["D"][91.2 ] = ["/DELPHI_1996_S3430090/d19-x01-y01","/OPAL_2004_S6132243/d14-x01-y01"]
analyses["EventShapes"]["D"][130.1] = ["/L3_2004_I652683/d44-x01-y01"]
analyses["EventShapes"]["D"][133.0] = ["/OPAL_2004_S6132243/d14-x01-y02"]
analyses["EventShapes"]["D"][136.3] = ["/L3_2004_I652683/d44-x01-y02"]
analyses["EventShapes"]["D"][161.3] = ["/L3_2004_I652683/d44-x01-y03"]
analyses["EventShapes"]["D"][172.3] = ["/L3_2004_I652683/d45-x01-y01"]
analyses["EventShapes"]["D"][177.0] = ["/OPAL_2004_S6132243/d14-x01-y03"]
analyses["EventShapes"]["D"][182.8] = ["/L3_2004_I652683/d45-x01-y02"]
analyses["EventShapes"]["D"][188.6] = ["/L3_2004_I652683/d45-x01-y03"]
analyses["EventShapes"]["D"][194.4] = ["/L3_2004_I652683/d46-x01-y01"]
analyses["EventShapes"]["D"][197.0] = ["/OPAL_2004_S6132243/d14-x01-y04"]
analyses["EventShapes"]["D"][200.2] = ["/L3_2004_I652683/d46-x01-y02"]
analyses["EventShapes"]["D"][206.2] = ["/L3_2004_I652683/d46-x01-y03"]
analyses["EventShapes"]["D"][183.0] = ["/DELPHI_2003_I620250/d54-x01-y01"]
analyses["EventShapes"]["D"][189.0] = ["/DELPHI_2003_I620250/d54-x01-y02"]
analyses["EventShapes"]["D"][192.0] = ["/DELPHI_2003_I620250/d54-x01-y03"]
analyses["EventShapes"]["D"][196.0] = ["/DELPHI_2003_I620250/d54-x01-y04"]
analyses["EventShapes"]["D"][200.0] = ["/DELPHI_2003_I620250/d55-x01-y01"]
analyses["EventShapes"]["D"][202.0] = ["/DELPHI_2003_I620250/d55-x01-y02"]
analyses["EventShapes"]["D"][205.0] = ["/DELPHI_2003_I620250/d55-x01-y03"]
analyses["EventShapes"]["D"][207.0] = ["/DELPHI_2003_I620250/d55-x01-y04"]
# hemispheres
analyses["EventShapes"]["HeavyJetMass"][14.0] = ["/TASSO_1989_I279165/d02-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][22.0] = ["/TASSO_1989_I279165/d02-x01-y02"]
analyses["EventShapes"]["HeavyJetMass"][34.8] = ["/TASSO_1989_I279165/d02-x01-y03"]
analyses["EventShapes"]["HeavyJetMass"][35.0] = ["/JADE_1998_S3612880/d07-x01-y01"]


analyses["EventShapes"]["HeavyJetMass"][43.5] = ["/TASSO_1989_I279165/d02-x01-y04"]
analyses["EventShapes"]["HeavyJetMass"][44.0] = ["/JADE_1998_S3612880/d03-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][45.0] = ["/DELPHI_2003_I620250/d07-x01-y01","/DELPHI_2003_I620250/d08-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][55.2] = ["/AMY_1990_I283337/d21-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][58.0] = ["/TOPAZ_1993_I361661/d02-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][66.0] = ["/DELPHI_2003_I620250/d07-x01-y02","/DELPHI_2003_I620250/d08-x01-y02"]
analyses["EventShapes"]["HeavyJetMass"][76.0] = ["/DELPHI_2003_I620250/d07-x01-y03","/DELPHI_2003_I620250/d08-x01-y03"]
analyses["EventShapes"]["HeavyJetMass"][91.2] = ["/DELPHI_1996_S3430090/d20-x01-y01","/ALEPH_1996_S3486095/d06-x01-y01","/OPAL_2004_S6132243/d02-x01-y01","/ALEPH_2004_S5765862/d62-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][133.0] = ["/OPAL_2004_S6132243/d02-x01-y02","/ALEPH_2004_S5765862/d63-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][161.0] = ["/ALEPH_2004_S5765862/d64-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][172.0] = ["/ALEPH_2004_S5765862/d65-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][177.0] = ["/OPAL_2004_S6132243/d02-x01-y03"]
analyses["EventShapes"]["HeavyJetMass"][183.0] = ["/DELPHI_2003_I620250/d56-x01-y01","/DELPHI_2003_I620250/d58-x01-y01",
                                                  "/DELPHI_2003_I620250/d60-x01-y01","/ALEPH_2004_S5765862/d66-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][189.0] = ["/DELPHI_2003_I620250/d56-x01-y02","/DELPHI_2003_I620250/d58-x01-y02",
                                                  "/DELPHI_2003_I620250/d60-x01-y02","/ALEPH_2004_S5765862/d67-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][192.0] = ["/DELPHI_2003_I620250/d56-x01-y03","/DELPHI_2003_I620250/d58-x01-y03",
                                                  "/DELPHI_2003_I620250/d60-x01-y03"]
analyses["EventShapes"]["HeavyJetMass"][196.0] = ["/DELPHI_2003_I620250/d56-x01-y04","/DELPHI_2003_I620250/d58-x01-y04",
                                                  "/DELPHI_2003_I620250/d60-x01-y04"]
analyses["EventShapes"]["HeavyJetMass"][197.0] = ["/OPAL_2004_S6132243/d02-x01-y04"]
analyses["EventShapes"]["HeavyJetMass"][200.0] = ["/DELPHI_2003_I620250/d57-x01-y01","/DELPHI_2003_I620250/d59-x01-y01",
                                                  "/DELPHI_2003_I620250/d61-x01-y01","/ALEPH_2004_S5765862/d68-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][202.0] = ["/DELPHI_2003_I620250/d57-x01-y02","/DELPHI_2003_I620250/d59-x01-y02",
                                                  "/DELPHI_2003_I620250/d61-x01-y02"]
analyses["EventShapes"]["HeavyJetMass"][205.0] = ["/DELPHI_2003_I620250/d57-x01-y03","/DELPHI_2003_I620250/d59-x01-y03",
                                                  "/DELPHI_2003_I620250/d61-x01-y03"]
analyses["EventShapes"]["HeavyJetMass"][206.0] = ["/ALEPH_2004_S5765862/d69-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][207.0] = ["/DELPHI_2003_I620250/d57-x01-y04","/DELPHI_2003_I620250/d59-x01-y04",
                                                  "/DELPHI_2003_I620250/d61-x01-y04"]
analyses["EventShapes"]["HeavyJetMass"][41.4 ] = ["/L3_2004_I652683/d26-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][55.3 ] = ["/L3_2004_I652683/d26-x01-y02"]
analyses["EventShapes"]["HeavyJetMass"][65.4 ] = ["/L3_2004_I652683/d26-x01-y03"]
analyses["EventShapes"]["HeavyJetMass"][75.7 ] = ["/L3_2004_I652683/d27-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][82.3 ] = ["/L3_2004_I652683/d27-x01-y02"]
analyses["EventShapes"]["HeavyJetMass"][85.1 ] = ["/L3_2004_I652683/d27-x01-y03"]
analyses["EventShapes"]["HeavyJetMass"][130.1] = ["/L3_2004_I652683/d28-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][136.3] = ["/L3_2004_I652683/d28-x01-y02"]
analyses["EventShapes"]["HeavyJetMass"][161.3] = ["/L3_2004_I652683/d28-x01-y03"]
analyses["EventShapes"]["HeavyJetMass"][172.3] = ["/L3_2004_I652683/d29-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][182.8] = ["/L3_2004_I652683/d29-x01-y02"]
analyses["EventShapes"]["HeavyJetMass"][188.6] = ["/L3_2004_I652683/d29-x01-y03"]
analyses["EventShapes"]["HeavyJetMass"][194.4] = ["/L3_2004_I652683/d30-x01-y01"]
analyses["EventShapes"]["HeavyJetMass"][200.2] = ["/L3_2004_I652683/d30-x01-y02"]
analyses["EventShapes"]["HeavyJetMass"][206.2] = ["/L3_2004_I652683/d30-x01-y03"]

analyses["EventShapes"]["JetMassDifference"][14.0] = ["/TASSO_1989_I279165/d01-x01-y01"]
analyses["EventShapes"]["JetMassDifference"][22.0] = ["/TASSO_1989_I279165/d01-x01-y02"]
analyses["EventShapes"]["JetMassDifference"][34.8] = ["/TASSO_1989_I279165/d01-x01-y03"]
analyses["EventShapes"]["JetMassDifference"][43.5] = ["/TASSO_1989_I279165/d01-x01-y04"]
analyses["EventShapes"]["JetMassDifference"][45.0] = ["/DELPHI_2003_I620250/d10-x01-y01"]
analyses["EventShapes"]["JetMassDifference"][55.2] = ["/AMY_1990_I283337/d22-x01-y01"]
analyses["EventShapes"]["JetMassDifference"][66.0] = ["/DELPHI_2003_I620250/d10-x01-y02"]
analyses["EventShapes"]["JetMassDifference"][76.0] = ["/DELPHI_2003_I620250/d10-x01-y03"]
analyses["EventShapes"]["JetMassDifference"][91.2] = ["/DELPHI_1996_S3430090/d22-x01-y01","/ALEPH_2004_S5765862/d110-x01-y01"]
analyses["EventShapes"]["JetMassDifference"][133.0] = ["/ALEPH_2004_S5765862/d111-x01-y01"]
analyses["EventShapes"]["JetMassDifference"][161.0] = ["/ALEPH_2004_S5765862/d112-x01-y01"]
analyses["EventShapes"]["JetMassDifference"][172.0] = ["/ALEPH_2004_S5765862/d113-x01-y01"]
analyses["EventShapes"]["JetMassDifference"][183.0] = ["/DELPHI_2003_I620250/d64-x01-y01","/ALEPH_2004_S5765862/d114-x01-y01"]
analyses["EventShapes"]["JetMassDifference"][189.0] = ["/DELPHI_2003_I620250/d64-x01-y02","/ALEPH_2004_S5765862/d115-x01-y01"]
analyses["EventShapes"]["JetMassDifference"][192.0] = ["/DELPHI_2003_I620250/d64-x01-y03"]
analyses["EventShapes"]["JetMassDifference"][196.0] = ["/DELPHI_2003_I620250/d64-x01-y04"]
analyses["EventShapes"]["JetMassDifference"][200.0] = ["/DELPHI_2003_I620250/d65-x01-y01","/ALEPH_2004_S5765862/d116-x01-y01"]
analyses["EventShapes"]["JetMassDifference"][202.0] = ["/DELPHI_2003_I620250/d65-x01-y02"]
analyses["EventShapes"]["JetMassDifference"][205.0] = ["/DELPHI_2003_I620250/d65-x01-y03"]
analyses["EventShapes"]["JetMassDifference"][206.0] = ["/ALEPH_2004_S5765862/d117-x01-y01"]
analyses["EventShapes"]["JetMassDifference"][207.0] = ["/DELPHI_2003_I620250/d65-x01-y04"]
analyses["EventShapes"]["LightJetMass"][14.0] = ["/TASSO_1989_I279165/d03-x01-y01"]
analyses["EventShapes"]["LightJetMass"][22.0] = ["/TASSO_1989_I279165/d03-x01-y02"]
analyses["EventShapes"]["LightJetMass"][34.8] = ["/TASSO_1989_I279165/d03-x01-y03"]
analyses["EventShapes"]["LightJetMass"][43.5] = ["/TASSO_1989_I279165/d03-x01-y04"]
analyses["EventShapes"]["LightJetMass"][45.0] = ["/DELPHI_2003_I620250/d09-x01-y01"]
analyses["EventShapes"]["LightJetMass"][55.2] = ["/AMY_1990_I283337/d20-x01-y01"]
analyses["EventShapes"]["LightJetMass"][66.0] = ["/DELPHI_2003_I620250/d09-x01-y02"]
analyses["EventShapes"]["LightJetMass"][76.0] = ["/DELPHI_2003_I620250/d09-x01-y03"]
analyses["EventShapes"]["LightJetMass"][91.2] = ["/DELPHI_1996_S3430090/d21-x01-y01","/OPAL_2004_S6132243/d12-x01-y01"]
analyses["EventShapes"]["LightJetMass"][133.0] = ["/OPAL_2004_S6132243/d12-x01-y02"]
analyses["EventShapes"]["LightJetMass"][177.0] = ["/OPAL_2004_S6132243/d12-x01-y03"]
analyses["EventShapes"]["LightJetMass"][183.0] = ["/DELPHI_2003_I620250/d62-x01-y01"]
analyses["EventShapes"]["LightJetMass"][189.0] = ["/DELPHI_2003_I620250/d62-x01-y02"]
analyses["EventShapes"]["LightJetMass"][192.0] = ["/DELPHI_2003_I620250/d62-x01-y03"]
analyses["EventShapes"]["LightJetMass"][196.0] = ["/DELPHI_2003_I620250/d62-x01-y04"]
analyses["EventShapes"]["LightJetMass"][197.0] = ["/OPAL_2004_S6132243/d12-x01-y04"]
analyses["EventShapes"]["LightJetMass"][200.0] = ["/DELPHI_2003_I620250/d63-x01-y01"]
analyses["EventShapes"]["LightJetMass"][202.0] = ["/DELPHI_2003_I620250/d63-x01-y02"]
analyses["EventShapes"]["LightJetMass"][205.0] = ["/DELPHI_2003_I620250/d63-x01-y03"]
analyses["EventShapes"]["LightJetMass"][207.0] = ["/DELPHI_2003_I620250/d63-x01-y04"]
analyses["EventShapes"]["TotalJetMass"][45.0] = ["/DELPHI_2003_I620250/d11-x01-y01","/DELPHI_2003_I620250/d12-x01-y01"]
analyses["EventShapes"]["TotalJetMass"][66.0] = ["/DELPHI_2003_I620250/d11-x01-y02","/DELPHI_2003_I620250/d12-x01-y02"]
analyses["EventShapes"]["TotalJetMass"][76.0] = ["/DELPHI_2003_I620250/d11-x01-y03","/DELPHI_2003_I620250/d12-x01-y03"]
# jets
# y12
analyses["EventShapes"]["y12_dur"][91.2 ] = ["/ALEPH_2004_S5765862/d149-x01-y01"]
analyses["EventShapes"]["y12_dur"][133.0] = ["/ALEPH_2004_S5765862/d150-x01-y01"]
analyses["EventShapes"]["y12_dur"][161.0] = ["/ALEPH_2004_S5765862/d151-x01-y01"]
analyses["EventShapes"]["y12_dur"][172.0] = ["/ALEPH_2004_S5765862/d152-x01-y01"]
analyses["EventShapes"]["y12_dur"][183.0] = ["/ALEPH_2004_S5765862/d153-x01-y01"]
analyses["EventShapes"]["y12_dur"][189.0] = ["/ALEPH_2004_S5765862/d154-x01-y01"]
analyses["EventShapes"]["y12_dur"][200.0] = ["/ALEPH_2004_S5765862/d155-x01-y01"]
analyses["EventShapes"]["y12_dur"][206.0] = ["/ALEPH_2004_S5765862/d156-x01-y01"]
# y23
analyses["EventShapes"]["y23_dur"][22.0] = ["/JADE_1998_S3612880/d12-x01-y01"]
analyses["EventShapes"]["y23_dur"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d24-x01-y01","/JADE_1998_S3612880/d11-x01-y01"]
analyses["EventShapes"]["y23_dur"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d25-x01-y01","/JADE_1998_S3612880/d10-x01-y01"]
analyses["EventShapes"]["y23_dur"][58.0 ] = ["/TOPAZ_1993_I361661/d03-x01-y01"]
analyses["EventShapes"]["y23_dur"][91.2 ] = ["/ALEPH_2004_S5765862/d157-x01-y01","/ALEPH_1996_S3486095/d05-x01-y01"
                                             "/OPAL_2004_S6132243/d06-x01-y01","/JADE_OPAL_2000_S4300807/d26-x01-y01",
                                             "/DELPHI_1996_S3430090/d27-x01-y01"]
analyses["EventShapes"]["y23_dur"][133.0] = ["/ALEPH_2004_S5765862/d158-x01-y01","/OPAL_2004_S6132243/d06-x01-y02",
                                             "/JADE_OPAL_2000_S4300807/d27-x01-y01"]
analyses["EventShapes"]["y23_dur"][161.0] = ["/ALEPH_2004_S5765862/d159-x01-y01","/JADE_OPAL_2000_S4300807/d28-x01-y01"]
analyses["EventShapes"]["y23_dur"][172.0] = ["/ALEPH_2004_S5765862/d160-x01-y01","/JADE_OPAL_2000_S4300807/d29-x01-y01"]
analyses["EventShapes"]["y23_dur"][177.0] = ["/OPAL_2004_S6132243/d06-x01-y03"]
analyses["EventShapes"]["y23_dur"][183.0] = ["/ALEPH_2004_S5765862/d161-x01-y01","/JADE_OPAL_2000_S4300807/d30-x01-y01"]
analyses["EventShapes"]["y23_dur"][189.0] = ["/ALEPH_2004_S5765862/d162-x01-y01","/JADE_OPAL_2000_S4300807/d31-x01-y01"]
analyses["EventShapes"]["y23_dur"][197.0] = ["/OPAL_2004_S6132243/d06-x01-y04"]
analyses["EventShapes"]["y23_dur"][200.0] = ["/ALEPH_2004_S5765862/d163-x01-y01"]
analyses["EventShapes"]["y23_dur"][206.0] = ["/ALEPH_2004_S5765862/d164-x01-y01"]
# y34                           
analyses["EventShapes"]["y34_dur"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d24-x01-y02"]
analyses["EventShapes"]["y34_dur"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d25-x01-y02"]
analyses["EventShapes"]["y34_dur"][91.2 ] = ["/ALEPH_2004_S5765862/d165-x01-y01","/JADE_OPAL_2000_S4300807/d26-x01-y02",
                                             "/DELPHI_1996_S3430090/d29-x01-y01"]
analyses["EventShapes"]["y34_dur"][133.0] = ["/ALEPH_2004_S5765862/d166-x01-y01","/JADE_OPAL_2000_S4300807/d27-x01-y02"]
analyses["EventShapes"]["y34_dur"][161.0] = ["/ALEPH_2004_S5765862/d167-x01-y01","/JADE_OPAL_2000_S4300807/d28-x01-y02"]
analyses["EventShapes"]["y34_dur"][172.0] = ["/ALEPH_2004_S5765862/d168-x01-y01","/JADE_OPAL_2000_S4300807/d29-x01-y02"]
analyses["EventShapes"]["y34_dur"][183.0] = ["/ALEPH_2004_S5765862/d169-x01-y01","/JADE_OPAL_2000_S4300807/d30-x01-y02"]
analyses["EventShapes"]["y34_dur"][189.0] = ["/ALEPH_2004_S5765862/d170-x01-y01","/JADE_OPAL_2000_S4300807/d31-x01-y02"]
analyses["EventShapes"]["y34_dur"][206.0] = ["/ALEPH_2004_S5765862/d172-x01-y01"]
# y45                           
analyses["EventShapes"]["y45_dur"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d24-x01-y03"]
analyses["EventShapes"]["y45_dur"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d25-x01-y03"]
analyses["EventShapes"]["y45_dur"][91.2 ] = ["/ALEPH_2004_S5765862/d173-x01-y01","/JADE_OPAL_2000_S4300807/d26-x01-y03",
                                             "/DELPHI_1996_S3430090/d31-x01-y01"]
analyses["EventShapes"]["y45_dur"][133.0] = ["/ALEPH_2004_S5765862/d174-x01-y01","/JADE_OPAL_2000_S4300807/d27-x01-y03"]
analyses["EventShapes"]["y45_dur"][161.0] = ["/ALEPH_2004_S5765862/d175-x01-y01","/JADE_OPAL_2000_S4300807/d28-x01-y03"]
analyses["EventShapes"]["y45_dur"][172.0] = ["/ALEPH_2004_S5765862/d176-x01-y01","/JADE_OPAL_2000_S4300807/d29-x01-y03"]
analyses["EventShapes"]["y45_dur"][183.0] = ["/ALEPH_2004_S5765862/d177-x01-y01","/JADE_OPAL_2000_S4300807/d30-x01-y03"]
analyses["EventShapes"]["y45_dur"][189.0] = ["/ALEPH_2004_S5765862/d178-x01-y01","/JADE_OPAL_2000_S4300807/d31-x01-y03"]
analyses["EventShapes"]["y45_dur"][200.0] = ["/ALEPH_2004_S5765862/d179-x01-y01"]
# y56                           
analyses["EventShapes"]["y56_dur"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d24-x01-y04"]
analyses["EventShapes"]["y56_dur"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d25-x01-y04"]
analyses["EventShapes"]["y56_dur"][91.2 ] = ["/ALEPH_2004_S5765862/d180-x01-y01","/JADE_OPAL_2000_S4300807/d26-x01-y04"]
analyses["EventShapes"]["y56_dur"][133.0] = ["/ALEPH_2004_S5765862/d181-x01-y01","/JADE_OPAL_2000_S4300807/d27-x01-y04"]
analyses["EventShapes"]["y56_dur"][161.0] = ["/ALEPH_2004_S5765862/d182-x01-y01","/JADE_OPAL_2000_S4300807/d28-x01-y04"]
analyses["EventShapes"]["y56_dur"][172.0] = ["/ALEPH_2004_S5765862/d183-x01-y01","/JADE_OPAL_2000_S4300807/d29-x01-y04"]
analyses["EventShapes"]["y56_dur"][183.0] = ["/ALEPH_2004_S5765862/d184-x01-y01","/JADE_OPAL_2000_S4300807/d30-x01-y04"]
analyses["EventShapes"]["y56_dur"][189.0] = ["/ALEPH_2004_S5765862/d185-x01-y01","/JADE_OPAL_2000_S4300807/d31-x01-y04"]
analyses["EventShapes"]["y56_dur"][200.0] = ["/ALEPH_2004_S5765862/d186-x01-y01"]
# jade scheme
analyses["EventShapes"]["y23_jade"][57.7 ] = ["/AMY_1995_I406129/d02-x01-y01","/AMY_1995_I406129/d03-x01-y01",
                                              "/AMY_1995_I406129/d04-x01-y01","/AMY_1995_I406129/d06-x01-y01"]
analyses["EventShapes"]["y23_jade"][91.2 ] = ["/DELPHI_1996_S3430090/d28-x01-y01"]
analyses["EventShapes"]["y34_jade"][91.2 ] = ["/DELPHI_1996_S3430090/d30-x01-y01"]
analyses["EventShapes"]["y45_jade"][91.2 ] = ["/DELPHI_1996_S3430090/d32-x01-y01"]
# jet fractions
# 1 jet
analyses["EventShapes"]["1jet_dur"][91.2 ] = ["/ALEPH_2004_S5765862/d187-x01-y01"]
analyses["EventShapes"]["1jet_dur"][133.0] = ["/ALEPH_2004_S5765862/d188-x01-y01"]
analyses["EventShapes"]["1jet_dur"][161.0] = ["/ALEPH_2004_S5765862/d189-x01-y01"]
analyses["EventShapes"]["1jet_dur"][172.0] = ["/ALEPH_2004_S5765862/d190-x01-y01"]
analyses["EventShapes"]["1jet_dur"][183.0] = ["/ALEPH_2004_S5765862/d191-x01-y01"]
analyses["EventShapes"]["1jet_dur"][189.0] = ["/ALEPH_2004_S5765862/d192-x01-y01"]
analyses["EventShapes"]["1jet_dur"][200.0] = ["/ALEPH_2004_S5765862/d193-x01-y01"]
analyses["EventShapes"]["1jet_dur"][206.0] = ["/ALEPH_2004_S5765862/d194-x01-y01"]
# 2 jet
analyses["EventShapes"]["2jet_dur"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d16-x01-y01"]
analyses["EventShapes"]["2jet_dur"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d17-x01-y01"]
analyses["EventShapes"]["2jet_dur"][91.2 ] = ["/ALEPH_2004_S5765862/d195-x01-y01","/JADE_OPAL_2000_S4300807/d18-x01-y01"]
analyses["EventShapes"]["2jet_dur"][133.0] = ["/ALEPH_2004_S5765862/d196-x01-y01","/JADE_OPAL_2000_S4300807/d19-x01-y01"]
analyses["EventShapes"]["2jet_dur"][161.0] = ["/ALEPH_2004_S5765862/d197-x01-y01","/JADE_OPAL_2000_S4300807/d20-x01-y01"]
analyses["EventShapes"]["2jet_dur"][172.0] = ["/ALEPH_2004_S5765862/d198-x01-y01","/JADE_OPAL_2000_S4300807/d21-x01-y01"]
analyses["EventShapes"]["2jet_dur"][183.0] = ["/ALEPH_2004_S5765862/d199-x01-y01","/JADE_OPAL_2000_S4300807/d22-x01-y01"]
analyses["EventShapes"]["2jet_dur"][189.0] = ["/ALEPH_2004_S5765862/d200-x01-y01","/JADE_OPAL_2000_S4300807/d23-x01-y01"]
analyses["EventShapes"]["2jet_dur"][200.0] = ["/ALEPH_2004_S5765862/d201-x01-y01"]
analyses["EventShapes"]["2jet_dur"][206.0] = ["/ALEPH_2004_S5765862/d202-x01-y01"]
# 3 jet
analyses["EventShapes"]["3jet_dur"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d16-x01-y02"]
analyses["EventShapes"]["3jet_dur"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d17-x01-y02"]
analyses["EventShapes"]["3jet_dur"][91.2 ] = ["/ALEPH_2004_S5765862/d203-x01-y01","/JADE_OPAL_2000_S4300807/d18-x01-y02"]
analyses["EventShapes"]["3jet_dur"][133.0] = ["/ALEPH_2004_S5765862/d204-x01-y01","/JADE_OPAL_2000_S4300807/d19-x01-y02"]
analyses["EventShapes"]["3jet_dur"][161.0] = ["/ALEPH_2004_S5765862/d205-x01-y01","/JADE_OPAL_2000_S4300807/d20-x01-y02"]
analyses["EventShapes"]["3jet_dur"][172.0] = ["/ALEPH_2004_S5765862/d206-x01-y01","/JADE_OPAL_2000_S4300807/d21-x01-y02"]
analyses["EventShapes"]["3jet_dur"][183.0] = ["/ALEPH_2004_S5765862/d207-x01-y01","/JADE_OPAL_2000_S4300807/d22-x01-y02"]
analyses["EventShapes"]["3jet_dur"][189.0] = ["/ALEPH_2004_S5765862/d208-x01-y01","/JADE_OPAL_2000_S4300807/d23-x01-y02"]
analyses["EventShapes"]["3jet_dur"][200.0] = ["/ALEPH_2004_S5765862/d209-x01-y01"]
analyses["EventShapes"]["3jet_dur"][206.0] = ["/ALEPH_2004_S5765862/d210-x01-y01"]
# 4 jet
analyses["EventShapes"]["4jet_dur"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d16-x01-y03"]
analyses["EventShapes"]["4jet_dur"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d17-x01-y03"]
analyses["EventShapes"]["4jet_dur"][91.2 ] = ["/ALEPH_2004_S5765862/d211-x01-y01","/JADE_OPAL_2000_S4300807/d18-x01-y03"]
analyses["EventShapes"]["4jet_dur"][133.0] = ["/ALEPH_2004_S5765862/d212-x01-y01","/JADE_OPAL_2000_S4300807/d19-x01-y03"]
analyses["EventShapes"]["4jet_dur"][161.0] = ["/ALEPH_2004_S5765862/d213-x01-y01","/JADE_OPAL_2000_S4300807/d20-x01-y03"]
analyses["EventShapes"]["4jet_dur"][172.0] = ["/ALEPH_2004_S5765862/d214-x01-y01","/JADE_OPAL_2000_S4300807/d21-x01-y03"]
analyses["EventShapes"]["4jet_dur"][183.0] = ["/ALEPH_2004_S5765862/d215-x01-y01","/JADE_OPAL_2000_S4300807/d22-x01-y03"]
analyses["EventShapes"]["4jet_dur"][189.0] = ["/ALEPH_2004_S5765862/d216-x01-y01","/JADE_OPAL_2000_S4300807/d23-x01-y03"]
analyses["EventShapes"]["4jet_dur"][200.0] = ["/ALEPH_2004_S5765862/d217-x01-y01"]
analyses["EventShapes"]["4jet_dur"][206.0] = ["/ALEPH_2004_S5765862/d218-x01-y01"]
# 5 jet
analyses["EventShapes"]["5jet_dur"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d16-x01-y04"]
analyses["EventShapes"]["5jet_dur"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d17-x01-y04"]
analyses["EventShapes"]["5jet_dur"][91.2 ] = ["/ALEPH_2004_S5765862/d219-x01-y01","/JADE_OPAL_2000_S4300807/d18-x01-y04"]
analyses["EventShapes"]["5jet_dur"][133.0] = ["/ALEPH_2004_S5765862/d220-x01-y01","/JADE_OPAL_2000_S4300807/d19-x01-y04"]
analyses["EventShapes"]["5jet_dur"][161.0] = ["/ALEPH_2004_S5765862/d221-x01-y01","/JADE_OPAL_2000_S4300807/d20-x01-y04"]
analyses["EventShapes"]["5jet_dur"][172.0] = ["/ALEPH_2004_S5765862/d222-x01-y01","/JADE_OPAL_2000_S4300807/d21-x01-y04"]
analyses["EventShapes"]["5jet_dur"][183.0] = ["/ALEPH_2004_S5765862/d223-x01-y01","/JADE_OPAL_2000_S4300807/d22-x01-y04"]
analyses["EventShapes"]["5jet_dur"][189.0] = ["/ALEPH_2004_S5765862/d224-x01-y01","/JADE_OPAL_2000_S4300807/d23-x01-y04"]
analyses["EventShapes"]["5jet_dur"][200.0] = ["/ALEPH_2004_S5765862/d225-x01-y01"]
analyses["EventShapes"]["5jet_dur"][206.0] = ["/ALEPH_2004_S5765862/d226-x01-y01"]
# 6 jet                                       
analyses["EventShapes"]["6jet_dur"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d16-x01-y05"]
analyses["EventShapes"]["6jet_dur"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d17-x01-y05"]
analyses["EventShapes"]["6jet_dur"][91.2 ] = ["/ALEPH_2004_S5765862/d227-x01-y01","/JADE_OPAL_2000_S4300807/d18-x01-y05"]
analyses["EventShapes"]["6jet_dur"][133.0] = ["/ALEPH_2004_S5765862/d228-x01-y01","/JADE_OPAL_2000_S4300807/d19-x01-y05"]
analyses["EventShapes"]["6jet_dur"][161.0] = ["/ALEPH_2004_S5765862/d229-x01-y01","/JADE_OPAL_2000_S4300807/d20-x01-y05"]
analyses["EventShapes"]["6jet_dur"][172.0] = ["/ALEPH_2004_S5765862/d230-x01-y01","/JADE_OPAL_2000_S4300807/d21-x01-y05"]
analyses["EventShapes"]["6jet_dur"][183.0] = ["/ALEPH_2004_S5765862/d231-x01-y01","/JADE_OPAL_2000_S4300807/d22-x01-y05"]
analyses["EventShapes"]["6jet_dur"][189.0] = ["/ALEPH_2004_S5765862/d232-x01-y01","/JADE_OPAL_2000_S4300807/d23-x01-y05"]
analyses["EventShapes"]["6jet_dur"][200.0] = ["/ALEPH_2004_S5765862/d233-x01-y01"]
analyses["EventShapes"]["6jet_dur"][206.0] = ["/ALEPH_2004_S5765862/d234-x01-y01"]
# EEC
analyses["EventShapes"]["EEC" ][7.7 ] = ["/PLUTO_1981_I156315/d01-x01-y01"]
analyses["EventShapes"]["EEC" ][9.4 ] = ["/PLUTO_1981_I156315/d01-x01-y02"]
analyses["EventShapes"]["EEC" ][12.0] = ["/PLUTO_1981_I156315/d01-x01-y03"]
analyses["EventShapes"]["EEC" ][13.0] = ["/PLUTO_1981_I156315/d01-x01-y04"]
analyses["EventShapes"]["EEC" ][14.0] = ["/TASSO_1987_I248660/d01-x01-y01","/JADE_1984_I202784/d01-x01-y01"]
analyses["EventShapes"]["EEC" ][17.0] = ["/PLUTO_1981_I156315/d01-x01-y05"]
analyses["EventShapes"]["EEC" ][22.0] = ["/TASSO_1987_I248660/d02-x01-y01","/JADE_1984_I202784/d01-x01-y02",
                                         "/CELLO_1982_I12010/d01-x01-y01","/PLUTO_1981_I156315/d01-x01-y06"]
analyses["EventShapes"]["EEC" ][27.6] = ["/PLUTO_1981_I156315/d01-x01-y07"]
analyses["EventShapes"]["EEC" ][29.0] = ["/MAC_1985_I202924/d01-x01-y01","/MAC_1985_I202924/d01-x01-y02"]
analyses["EventShapes"]["EEC" ][30.8] = ["/PLUTO_1981_I156315/d01-x01-y08"]
analyses["EventShapes"]["EEC" ][34.0] = ["/JADE_1984_I202784/d01-x01-y03","/CELLO_1982_I12010/d01-x01-y02"]
analyses["EventShapes"]["EEC" ][34.6] = ["/PLUTO_1985_I215869/d01-x01-y01"]
analyses["EventShapes"]["EEC" ][34.8] = ["/TASSO_1987_I248660/d03-x01-y01"]
analyses["EventShapes"]["EEC" ][43.5] = ["/TASSO_1987_I248660/d04-x01-y01"]
analyses["EventShapes"]["EEC" ][53.3] = ["/TOPAZ_1989_I279575/d01-x01-y01","/TOPAZ_1989_I279575/d01-x01-y02"]
analyses["EventShapes"]["EEC" ][91.2] = ["/DELPHI_1996_S3430090/d33-x01-y01"]
analyses["EventShapes"]["EEC" ][59.5] = ["/TOPAZ_1989_I279575/d02-x01-y01","/TOPAZ_1989_I279575/d02-x01-y02"]
# AEEC 
analyses["EventShapes"]["AEEC"][8.65] = ["/PLUTO_1981_I156315/d04-x01-y01"]
analyses["EventShapes"]["AEEC"][14.0] = ["/JADE_1984_I202784/d02-x01-y01"]
analyses["EventShapes"]["AEEC"][22.0] = ["/JADE_1984_I202784/d02-x01-y02","/CELLO_1982_I12010/d03-x01-y01"]
analyses["EventShapes"]["AEEC"][29.0] = ["/MAC_1985_I202924/d01-x01-y03"]
analyses["EventShapes"]["AEEC"][30.8] = ["/PLUTO_1981_I156315/d05-x01-y01"]
analyses["EventShapes"]["AEEC"][34.0] = ["/JADE_1984_I202784/d02-x01-y03","/CELLO_1982_I12010/d03-x01-y02"]
analyses["EventShapes"]["AEEC"][53.3] = ["/TOPAZ_1989_I279575/d01-x01-y03"]
analyses["EventShapes"]["AEEC"][59.5] = ["/TOPAZ_1989_I279575/d02-x01-y03"]
analyses["EventShapes"]["AEEC"][91.2] = ["/DELPHI_1996_S3430090/d34-x01-y01"]
# sphericity based 
analyses["EventShapes"]["S"][12.0 ] = ["/TASSO_1980_I153511/d01-x01-y01"]
analyses["EventShapes"]["S"][14.0 ] = ["/TASSO_1990_S2148048/d06-x01-y01"]
analyses["EventShapes"]["S"][22.0 ] = ["/TASSO_1990_S2148048/d06-x01-y02"]
analyses["EventShapes"]["S"][29.0 ] = ["/HRS_1985_I201482/d01-x01-y01"]
analyses["EventShapes"]["S"][30.8 ] = ["/TASSO_1980_I153511/d02-x01-y01"]
analyses["EventShapes"]["S"][35.0 ] = ["/TASSO_1990_S2148048/d06-x01-y03","/TASSO_1988_I263859/d01-x01-y01"]
analyses["EventShapes"]["S"][44.0 ] = ["/TASSO_1990_S2148048/d06-x01-y04"]
analyses["EventShapes"]["S"][45.0 ] = ["/DELPHI_2003_I620250/d04-x01-y01"]
analyses["EventShapes"]["S"][55.2 ] = ["/AMY_1990_I283337/d16-x01-y01"]
analyses["EventShapes"]["S"][66.0 ] = ["/DELPHI_2003_I620250/d04-x01-y02"]
analyses["EventShapes"]["S"][76.0 ] = ["/DELPHI_2003_I620250/d04-x01-y03"]
analyses["EventShapes"]["S"][91.2 ] = ["/ALEPH_2004_S5765862/d141-x01-y01","/DELPHI_1996_S3430090/d15-x01-y01",
                                       "/ALEPH_1996_S3486095/d01-x01-y01","/OPAL_2004_S6132243/d10-x01-y01"]
analyses["EventShapes"]["S"][133.0] = ["/ALEPH_2004_S5765862/d142-x01-y01","/OPAL_2004_S6132243/d10-x01-y02"]
analyses["EventShapes"]["S"][161.0] = ["/ALEPH_2004_S5765862/d143-x01-y01"]
analyses["EventShapes"]["S"][172.0] = ["/ALEPH_2004_S5765862/d144-x01-y01"]
analyses["EventShapes"]["S"][177.0] = ["/OPAL_2004_S6132243/d10-x01-y03"]
analyses["EventShapes"]["S"][183.0] = ["/DELPHI_2003_I620250/d66-x01-y01","/ALEPH_2004_S5765862/d145-x01-y01"]
analyses["EventShapes"]["S"][189.0] = ["/DELPHI_2003_I620250/d66-x01-y02","/ALEPH_2004_S5765862/d146-x01-y01"]
analyses["EventShapes"]["S"][192.0] = ["/DELPHI_2003_I620250/d66-x01-y03"]
analyses["EventShapes"]["S"][196.0] = ["/DELPHI_2003_I620250/d66-x01-y04"]
analyses["EventShapes"]["S"][197.0] = ["/OPAL_2004_S6132243/d10-x01-y04"]
analyses["EventShapes"]["S"][200.0] = ["/DELPHI_2003_I620250/d67-x01-y01","/ALEPH_2004_S5765862/d147-x01-y01"]
analyses["EventShapes"]["S"][202.0] = ["/DELPHI_2003_I620250/d67-x01-y02"]
analyses["EventShapes"]["S"][205.0] = ["/DELPHI_2003_I620250/d67-x01-y03"]
analyses["EventShapes"]["S"][206.0] = ["/ALEPH_2004_S5765862/d148-x01-y01"]
analyses["EventShapes"]["S"][207.0] = ["/DELPHI_2003_I620250/d67-x01-y04"]

analyses["QED"] = ["/ALEPH_1996_S3196992/d03-x01-y01","/ALEPH_1996_S3196992/d04-x01-y01",
                   "/ALEPH_1996_S3196992/d01-x01-y01","/ALEPH_1996_S3196992/d02-x01-y01",
                   "/ALEPH_1996_S3196992/d05-x01-y01","/ALEPH_1996_S3196992/d06-x01-y01",
                   "/ALEPH_1996_S3196992/d07-x01-y01","/ALEPH_1996_S3196992/d08-x01-y01",
                   "/OPAL_1993_S2692198/d01-x01-y01","/OPAL_1993_S2692198/d02-x01-y01",
                   "/OPAL_1993_S2692198/d03-x01-y01","/OPAL_1993_S2692198/d03-x01-y02",
                   "/OPAL_1993_S2692198/d03-x01-y03","/OPAL_1993_S2692198/d03-x01-y04",
                   "/OPAL_1993_S2692198/d04-x01-y01","/OPAL_1993_S2692198/d04-x01-y02",
                   "/OPAL_1993_S2692198/d04-x01-y03","/OPAL_1993_S2692198/d04-x01-y04",]

analyses["EventShapes"]["P"][45.0 ] = ["/DELPHI_2003_I620250/d05-x01-y01"]
analyses["EventShapes"]["P"][66.0 ] = ["/DELPHI_2003_I620250/d05-x01-y02"]
analyses["EventShapes"]["P"][76.0 ] = ["/DELPHI_2003_I620250/d05-x01-y03"]
analyses["EventShapes"]["P"][91.2 ] = ["/DELPHI_1996_S3430090/d17-x01-y01"]
analyses["EventShapes"]["P"][133.0] = ["/ALEPH_2004_S5765862/d126-x01-y01"]
analyses["EventShapes"]["P"][161.0] = ["/ALEPH_2004_S5765862/d127-x01-y01"]
analyses["EventShapes"]["P"][172.0] = ["/ALEPH_2004_S5765862/d128-x01-y01"]
analyses["EventShapes"]["P"][183.0] = ["/DELPHI_2003_I620250/d68-x01-y01","/ALEPH_2004_S5765862/d129-x01-y01"]
analyses["EventShapes"]["P"][189.0] = ["/DELPHI_2003_I620250/d68-x01-y02","/ALEPH_2004_S5765862/d130-x01-y01"]
analyses["EventShapes"]["P"][192.0] = ["/DELPHI_2003_I620250/d68-x01-y03"]
analyses["EventShapes"]["P"][196.0] = ["/DELPHI_2003_I620250/d68-x01-y04"]
analyses["EventShapes"]["P"][200.0] = ["/DELPHI_2003_I620250/d69-x01-y01","/ALEPH_2004_S5765862/d131-x01-y01"]
analyses["EventShapes"]["P"][202.0] = ["/DELPHI_2003_I620250/d69-x01-y02"]
analyses["EventShapes"]["P"][205.0] = ["/DELPHI_2003_I620250/d69-x01-y03"]
analyses["EventShapes"]["P"][206.0] = ["/ALEPH_2004_S5765862/d132-x01-y01"]
analyses["EventShapes"]["P"][207.0] = ["/DELPHI_2003_I620250/d69-x01-y04"]

analyses["EventShapes"]["A"][12.0 ] = ["/TASSO_1980_I153511/d03-x01-y01"]
analyses["EventShapes"]["A"][14.0 ] = ["/TASSO_1990_S2148048/d07-x01-y01"]
analyses["EventShapes"]["A"][22.0 ] = ["/TASSO_1990_S2148048/d07-x01-y02","/HRS_1985_I201482/d06-x01-y01"]
analyses["EventShapes"]["A"][30.8 ] = ["/TASSO_1980_I153511/d04-x01-y01"]
analyses["EventShapes"]["A"][35.0 ] = ["/TASSO_1990_S2148048/d07-x01-y03","/TASSO_1988_I263859/d02-x01-y01"]
analyses["EventShapes"]["A"][44.0 ] = ["/TASSO_1990_S2148048/d07-x01-y04"]
analyses["EventShapes"]["A"][55.2 ] = ["/AMY_1990_I283337/d17-x01-y01"]
analyses["EventShapes"]["A"][91.2 ] = ["/ALEPH_2004_S5765862/d118-x01-y01","/DELPHI_1996_S3430090/d16-x01-y01",
                                       "/ALEPH_1996_S3486095/d02-x01-y01","/OPAL_2004_S6132243/d09-x01-y01"]
analyses["EventShapes"]["A"][133.0] = ["/ALEPH_2004_S5765862/d119-x01-y01","/OPAL_2004_S6132243/d09-x01-y02"]
analyses["EventShapes"]["A"][161.0] = ["/ALEPH_2004_S5765862/d120-x01-y01"]
analyses["EventShapes"]["A"][172.0] = ["/ALEPH_2004_S5765862/d121-x01-y01"]
analyses["EventShapes"]["A"][177.0] = ["/OPAL_2004_S6132243/d09-x01-y03"]
analyses["EventShapes"]["A"][183.0] = ["/DELPHI_2003_I620250/d70-x01-y01","/ALEPH_2004_S5765862/d122-x01-y01"]
analyses["EventShapes"]["A"][189.0] = ["/DELPHI_2003_I620250/d70-x01-y02","/ALEPH_2004_S5765862/d123-x01-y01"]
analyses["EventShapes"]["A"][192.0] = ["/DELPHI_2003_I620250/d70-x01-y03"]
analyses["EventShapes"]["A"][196.0] = ["/DELPHI_2003_I620250/d70-x01-y04"]
analyses["EventShapes"]["A"][197.0] = ["/OPAL_2004_S6132243/d09-x01-y04"]
analyses["EventShapes"]["A"][200.0] = ["/DELPHI_2003_I620250/d71-x01-y01","/ALEPH_2004_S5765862/d124-x01-y01"]
analyses["EventShapes"]["A"][202.0] = ["/DELPHI_2003_I620250/d71-x01-y02"]
analyses["EventShapes"]["A"][205.0] = ["/DELPHI_2003_I620250/d71-x01-y03"]
analyses["EventShapes"]["A"][206.0] = ["/ALEPH_2004_S5765862/d125-x01-y01"]
analyses["EventShapes"]["A"][207.0] = ["/DELPHI_2003_I620250/d71-x01-y04"]

analyses["EventShapes"]["O"][45.0 ] = ["/DELPHI_2003_I620250/d06-x01-y01"]
analyses["EventShapes"]["O"][55.2 ] = ["/AMY_1990_I283337/d15-x01-y01"]
analyses["EventShapes"]["O"][66.0 ] = ["/DELPHI_2003_I620250/d06-x01-y02"]
analyses["EventShapes"]["O"][76.0 ] = ["/DELPHI_2003_I620250/d06-x01-y03"]
analyses["EventShapes"]["O"][91.2 ] = ["/ALEPH_2004_S5765862/d133-x01-y01","/DELPHI_1996_S3430090/d14-x01-y01",
                                       "/ALEPH_1996_S3486095/d08-x01-y01","/OPAL_2004_S6132243/d11-x01-y01"]
analyses["EventShapes"]["O"][133.0] = ["/ALEPH_2004_S5765862/d134-x01-y01","/OPAL_2004_S6132243/d11-x01-y02"]
analyses["EventShapes"]["O"][161.0] = ["/ALEPH_2004_S5765862/d135-x01-y01"]
analyses["EventShapes"]["O"][172.0] = ["/ALEPH_2004_S5765862/d136-x01-y01"]
analyses["EventShapes"]["O"][177.0] = ["/OPAL_2004_S6132243/d11-x01-y03"]
analyses["EventShapes"]["O"][183.0] = ["/DELPHI_2003_I620250/d44-x01-y01","/ALEPH_2004_S5765862/d137-x01-y01"]
analyses["EventShapes"]["O"][189.0] = ["/DELPHI_2003_I620250/d44-x01-y02","/ALEPH_2004_S5765862/d138-x01-y01"]
analyses["EventShapes"]["O"][192.0] = ["/DELPHI_2003_I620250/d44-x01-y03"]
analyses["EventShapes"]["O"][196.0] = ["/DELPHI_2003_I620250/d44-x01-y04"]
analyses["EventShapes"]["O"][197.0] = ["/OPAL_2004_S6132243/d11-x01-y04"]
analyses["EventShapes"]["O"][200.0] = ["/DELPHI_2003_I620250/d45-x01-y01","/ALEPH_2004_S5765862/d139-x01-y01"]
analyses["EventShapes"]["O"][202.0] = ["/DELPHI_2003_I620250/d45-x01-y02"]
analyses["EventShapes"]["O"][205.0] = ["/DELPHI_2003_I620250/d45-x01-y03"]
analyses["EventShapes"]["O"][206.0] = ["/ALEPH_2004_S5765862/d140-x01-y01"]
analyses["EventShapes"]["O"][207.0] = ["/DELPHI_2003_I620250/d45-x01-y04"]

analyses["TauDecays"]["KK"   ] = ["/BABAR_2018_I1679886/d01-x01-y01","/MC_TAU_Decay/h_2B_m2KK","/MC_TAU_Decay/h_2B_mKK"]
analyses["TauDecays"]["Kpi"  ] = ["/BELLE_2007_I753243/d01-x01-y01","/MC_TAU_Decay/h_2B_m2KpiA","/MC_TAU_Decay/h_2B_m2KpiB",
                                  "/MC_TAU_Decay/h_2B_mKpiA","/MC_TAU_Decay/h_2B_mKpiB"]
analyses["TauDecays"]["2pi"  ] = ["/BELLE_2008_I786560/d01-x01-y01","/ALEPH_2014_I1267648/d01-x01-y01","/CLEO_1999_I508944/d01-x01-y01",
                                  "/MC_TAU_Decay/h_2B_m2pipi","/MC_TAU_Decay/h_2B_mpipi"]
analyses["TauDecays"]["3pi"  ] = ["/BELLE_2010_I841618/d01-x01-y01","/BABAR_2007_S7266081/d01-x01-y01",
                                  "/BABAR_2007_S7266081/d02-x01-y01","/BABAR_2007_S7266081/d11-x01-y01",
                                  "/BABAR_2005_S6181155/d01-x01-y01","/ALEPH_2014_I1267648/d02-x01-y01","/ALEPH_2014_I1267648/d04-x01-y01",
                                  "/MC_TAU_Decay/h_3B_pi0pi0pim_1","/MC_TAU_Decay/h_3B_pi0pi0pim_2","/MC_TAU_Decay/h_3B_pi0pi0pim_3",
                                  "/MC_TAU_Decay/h_3B_pippimpim_1","/MC_TAU_Decay/h_3B_pippimpim_2","/MC_TAU_Decay/h_3B_pippimpim_3"]
analyses["TauDecays"]["Kpipi"] = ["/BELLE_2010_I841618/d02-x01-y01" ,"/BABAR_2007_S7266081/d03-x01-y01",
                                  "/BABAR_2007_S7266081/d04-x01-y01","/BABAR_2007_S7266081/d05-x01-y01",
                                  "/BABAR_2007_S7266081/d12-x01-y01",
                                  "/MC_TAU_Decay/h_3B_pi0pi0km_1","/MC_TAU_Decay/h_3B_pi0pi0km_2","/MC_TAU_Decay/h_3B_pi0pi0km_3",
                                  "/MC_TAU_Decay/h_3B_pimk0pi0_1","/MC_TAU_Decay/h_3B_pimk0pi0_2","/MC_TAU_Decay/h_3B_pimk0pi0_3",
                                  "/MC_TAU_Decay/h_3B_pimk0pi0_4"]
analyses["TauDecays"]["KKpi" ] = ["/BELLE_2010_I841618/d03-x01-y01","/BABAR_2007_S7266081/d06-x01-y01",
                                  "/BABAR_2007_S7266081/d07-x01-y01","/BABAR_2007_S7266081/d08-x01-y01",
                                  "/BABAR_2007_S7266081/d13-x01-y01",
                                  "/MC_TAU_Decay/h_3B_klpimkl_1","/MC_TAU_Decay/h_3B_klpimkl_2","/MC_TAU_Decay/h_3B_klpimkl_3",
                                  "/MC_TAU_Decay/h_3B_kmpi0k0_1","/MC_TAU_Decay/h_3B_kmpi0k0_2","/MC_TAU_Decay/h_3B_kmpi0k0_3",
                                  "/MC_TAU_Decay/h_3B_kmpi0k0_4","/MC_TAU_Decay/h_3B_kmpimkp_1","/MC_TAU_Decay/h_3B_kmpimkp_2",
                                  "/MC_TAU_Decay/h_3B_kmpimkp_3","/MC_TAU_Decay/h_3B_kmpimkp_4","/MC_TAU_Decay/h_3B_kmpimpip_1",
                                  "/MC_TAU_Decay/h_3B_kmpimpip_2","/MC_TAU_Decay/h_3B_kmpimpip_3","/MC_TAU_Decay/h_3B_kmpimpip_4",
                                  "/MC_TAU_Decay/h_3B_kspimkl_1","/MC_TAU_Decay/h_3B_kspimkl_2","/MC_TAU_Decay/h_3B_kspimkl_3",
                                  "/MC_TAU_Decay/h_3B_kspimkl_4","/MC_TAU_Decay/h_3B_kspimks_1","/MC_TAU_Decay/h_3B_kspimks_2",
                                  "/MC_TAU_Decay/h_3B_kspimks_3"]
analyses["TauDecays"]["3K"   ] = ["/BELLE_2010_I841618/d04-x01-y01","/BABAR_2007_S7266081/d09-x01-y01",
                                  "/BABAR_2007_S7266081/d10-x01-y01","/BABAR_2007_S7266081/d14-x01-y01"]
analyses["TauDecays"]["Keta"] = ["/MC_TAU_Decay/h_2B_m2Keta","/MC_TAU_Decay/h_2B_mKeta"]
analyses["TauDecays"]["lnu" ] = ["/MC_TAU_Decay/h_2B_m2enu","/MC_TAU_Decay/h_2B_m2munu",
                                 "/MC_TAU_Decay/h_2B_menu","/MC_TAU_Decay/h_2B_mmunu"]
analyses["TauDecays"]["2pieta"] = ["/MC_TAU_Decay/h_3B_pimpi0eta_1","/MC_TAU_Decay/h_3B_pimpi0eta_2",
                                   "/MC_TAU_Decay/h_3B_pimpi0eta_3","/MC_TAU_Decay/h_3B_pimpi0eta_4"]
analyses["TauDecays"]["2pigamma"] = ["/MC_TAU_Decay/h_3B_pimpi0gamma_1","/MC_TAU_Decay/h_3B_pimpi0gamma_2",
                                     "/MC_TAU_Decay/h_3B_pimpi0gamma_3","/MC_TAU_Decay/h_3B_pimpi0gamma_4"]
analyses["TauDecays"]["4pi"] = ["/ALEPH_2014_I1267648/d03-x01-y01","/ALEPH_2014_I1267648/d05-x01-y01",
                                "/MC_TAU_Decay/h_4B_pipi_1","/MC_TAU_Decay/h_4B_pipi_2",
                                "/MC_TAU_Decay/h_4B_pipi_3","/MC_TAU_Decay/h_4B_pipi_4",
                                "/MC_TAU_Decay/h_4B_pipi_5","/MC_TAU_Decay/h_4B_pipi_6",
                                "/MC_TAU_Decay/h_4B_pipipi_1","/MC_TAU_Decay/h_4B_pipipi_2",
                                "/MC_TAU_Decay/h_4B_pipipi_3","/MC_TAU_Decay/h_4B_pipipi_4",
                                "/MC_TAU_Decay/h_4B_pipipi_5","/MC_TAU_Decay/h_4B_pipipipi_1",
                                "/MC_TAU_Decay/h_4B_pipipipi_2"]
analyses["TauDecays"]["5pi"] = ["/MC_TAU_Decay/h_5B_pipi1_1","/MC_TAU_Decay/h_5B_pipi1_2",
                                "/MC_TAU_Decay/h_5B_pipi1_3","/MC_TAU_Decay/h_5B_pipi1_4",
                                "/MC_TAU_Decay/h_5B_pipi1_5","/MC_TAU_Decay/h_5B_pipi2_1",
                                "/MC_TAU_Decay/h_5B_pipi2_2","/MC_TAU_Decay/h_5B_pipi3_1",
                                "/MC_TAU_Decay/h_5B_pipi3_2","/MC_TAU_Decay/h_5B_pipi3_3",
                                "/MC_TAU_Decay/h_5B_pipipi1_1","/MC_TAU_Decay/h_5B_pipipi1_2",
                                "/MC_TAU_Decay/h_5B_pipipi1_3","/MC_TAU_Decay/h_5B_pipipi1_4",
                                "/MC_TAU_Decay/h_5B_pipipi1_5","/MC_TAU_Decay/h_5B_pipipi2_1",
                                "/MC_TAU_Decay/h_5B_pipipi2_2","/MC_TAU_Decay/h_5B_pipipi3_1",
                                "/MC_TAU_Decay/h_5B_pipipi3_2","/MC_TAU_Decay/h_5B_pipipi3_3",
                                "/MC_TAU_Decay/h_5B_pipipipi1_1","/MC_TAU_Decay/h_5B_pipipipi1_2",
                                "/MC_TAU_Decay/h_5B_pipipipi1_3","/MC_TAU_Decay/h_5B_pipipipi2_1",
                                "/MC_TAU_Decay/h_5B_pipipipi2_2","/MC_TAU_Decay/h_5B_pipipipi3_1",
                                "/MC_TAU_Decay/h_5B_pipipipi3_2","/MC_TAU_Decay/h_5B_q1",
                                "/MC_TAU_Decay/h_5B_q2","/MC_TAU_Decay/h_5B_q3"]
analyses["EventShapes"]["2jet_jade"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d07-x01-y01"]
analyses["EventShapes"]["3jet_jade"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d07-x01-y02"]
analyses["EventShapes"]["4jet_jade"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d07-x01-y03"]
analyses["EventShapes"]["5jet_jade"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d07-x01-y04"]
analyses["EventShapes"]["6jet_jade"][35.0 ] = ["/JADE_OPAL_2000_S4300807/d07-x01-y05"]
analyses["EventShapes"]["2jet_jade"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d08-x01-y01"]
analyses["EventShapes"]["3jet_jade"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d08-x01-y02"]
analyses["EventShapes"]["4jet_jade"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d08-x01-y03"]
analyses["EventShapes"]["5jet_jade"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d08-x01-y04"]
analyses["EventShapes"]["6jet_jade"][44.0 ] = ["/JADE_OPAL_2000_S4300807/d08-x01-y05"]
analyses["EventShapes"]["2jet_jade"][91.2 ] = ["/JADE_OPAL_2000_S4300807/d09-x01-y01"]
analyses["EventShapes"]["3jet_jade"][91.2 ] = ["/JADE_OPAL_2000_S4300807/d09-x01-y02"]
analyses["EventShapes"]["4jet_jade"][91.2 ] = ["/JADE_OPAL_2000_S4300807/d09-x01-y03"]
analyses["EventShapes"]["5jet_jade"][91.2 ] = ["/JADE_OPAL_2000_S4300807/d09-x01-y04"]
analyses["EventShapes"]["6jet_jade"][91.2 ] = ["/JADE_OPAL_2000_S4300807/d09-x01-y05"]
analyses["EventShapes"]["2jet_jade"][133.0] = ["/JADE_OPAL_2000_S4300807/d10-x01-y01"]
analyses["EventShapes"]["3jet_jade"][133.0] = ["/JADE_OPAL_2000_S4300807/d10-x01-y02"]
analyses["EventShapes"]["4jet_jade"][133.0] = ["/JADE_OPAL_2000_S4300807/d10-x01-y03"]
analyses["EventShapes"]["5jet_jade"][133.0] = ["/JADE_OPAL_2000_S4300807/d10-x01-y04"]
analyses["EventShapes"]["6jet_jade"][133.0] = ["/JADE_OPAL_2000_S4300807/d10-x01-y05"]
analyses["EventShapes"]["2jet_jade"][161.0] = ["/JADE_OPAL_2000_S4300807/d11-x01-y01"]
analyses["EventShapes"]["3jet_jade"][161.0] = ["/JADE_OPAL_2000_S4300807/d11-x01-y02"]
analyses["EventShapes"]["4jet_jade"][161.0] = ["/JADE_OPAL_2000_S4300807/d11-x01-y03"]
analyses["EventShapes"]["5jet_jade"][161.0] = ["/JADE_OPAL_2000_S4300807/d11-x01-y04"]
analyses["EventShapes"]["6jet_jade"][161.0] = ["/JADE_OPAL_2000_S4300807/d11-x01-y05"]
analyses["EventShapes"]["2jet_jade"][172.0] = ["/JADE_OPAL_2000_S4300807/d12-x01-y01"]
analyses["EventShapes"]["3jet_jade"][172.0] = ["/JADE_OPAL_2000_S4300807/d12-x01-y02"]
analyses["EventShapes"]["4jet_jade"][172.0] = ["/JADE_OPAL_2000_S4300807/d12-x01-y03"]
analyses["EventShapes"]["5jet_jade"][172.0] = ["/JADE_OPAL_2000_S4300807/d12-x01-y04"]
analyses["EventShapes"]["6jet_jade"][172.0] = ["/JADE_OPAL_2000_S4300807/d12-x01-y05"]
analyses["EventShapes"]["2jet_jade"][183.0] = ["/JADE_OPAL_2000_S4300807/d13-x01-y01"]
analyses["EventShapes"]["3jet_jade"][183.0] = ["/JADE_OPAL_2000_S4300807/d13-x01-y02"]
analyses["EventShapes"]["4jet_jade"][183.0] = ["/JADE_OPAL_2000_S4300807/d13-x01-y03"]
analyses["EventShapes"]["5jet_jade"][183.0] = ["/JADE_OPAL_2000_S4300807/d13-x01-y04"]
analyses["EventShapes"]["6jet_jade"][183.0] = ["/JADE_OPAL_2000_S4300807/d13-x01-y05"]
analyses["EventShapes"]["2jet_jade"][189.0] = ["/JADE_OPAL_2000_S4300807/d14-x01-y01"]
analyses["EventShapes"]["3jet_jade"][189.0] = ["/JADE_OPAL_2000_S4300807/d14-x01-y02"]
analyses["EventShapes"]["4jet_jade"][189.0] = ["/JADE_OPAL_2000_S4300807/d14-x01-y03"]
analyses["EventShapes"]["5jet_jade"][189.0] = ["/JADE_OPAL_2000_S4300807/d14-x01-y04"]
analyses["EventShapes"]["6jet_jade"][189.0] = ["/JADE_OPAL_2000_S4300807/d14-x01-y05"]


# /OPAL_2001_S4553896/d03-x01-y01
# /OPAL_2001_S4553896/d04-x01-y01
# /OPAL_2001_S4553896/d05-x01-y01
# /OPAL_2001_S4553896/d06-x01-y01
# /OPAL_2004_I631361/d05-x01-y01
# /OPAL_2004_I631361/d05-x01-y02
# /OPAL_2004_I648738/d06-x01-y01
# /OPAL_2004_I648738/d06-x01-y02
# /OPAL_2004_I648738/d06-x01-y03
# /OPAL_2004_I648738/d07-x01-y01
# /OPAL_2004_I648738/d07-x01-y02
# /OPAL_2004_I648738/d07-x01-y03
# /OPAL_2004_I648738/d08-x01-y01
# /OPAL_2004_I648738/d08-x01-y02
# /OPAL_2004_I648738/d08-x01-y03
# /OPAL_2004_I648738/d09-x01-y01
# /OPAL_2004_I648738/d09-x01-y02
# /OPAL_2004_I648738/d09-x01-y03
# /OPAL_2004_I648738/d10-x01-y01
# /OPAL_2004_I648738/d10-x01-y02
# /OPAL_2004_I648738/d11-x01-y01
# /OPAL_2004_I648738/d11-x01-y02



figures=glob.glob("%s/*/*.dat" % directory)

plotOutput="""<div style="float:left; font-size:smaller; font-weight:bold;">
    <a href="#{name}">&#9875;</a><a href="{dat}">&#8984</a> {name}: {energy} GeV<br/>
    <a name="{name}"><a href="{pdf}">
      <img src="{png}">
    </a></a>
  </div>"""
    
def writePlots(plots,output) :
    global figures
    output.write("<div style=\"float:none; overflow:auto; width:100%\">\n")
    for energy in sorted(plots.keys()) :
        try :
            float(energy)
        except:
            continue
        for name in sorted(plots[energy]) :
            dat=name[1:] +".dat"
            figName = ("%s/%s" %(directory,dat))
            if(figName not in figures) : continue
            del figures[figures.index(figName)]
            pdf=name[1:] +".pdf"
            png=name[1:] +".png"
            output.write(plotOutput.format(name=name[1:],pdf=pdf,png=png,dat=dat,energy=energy))
            output.write("\n")
    output.write("</div>")
    
def writePlots2(plots,output) :
    global figures
    output.write("<div style=\"float:none; overflow:auto; width:100%\">\n")
    for name in sorted(plots) :
        dat=name[1:] +".dat"
        figName = ("%s/%s" %(directory,dat))
        if(figName not in figures) : continue
        del figures[figures.index(figName)]
        pdf=name[1:] +".pdf"
        png=name[1:] +".png"
        output.write(plotOutput.format(name=name[1:],pdf=pdf,png=png,dat=dat,energy=""))
        output.write("\n")
    output.write("</div>")

particleNames = {
    22 : "$\\gamma$",
    111 : "$\\pi^0$", 211 : "$\\pi^\\pm$", 221 : "$\\eta$", 331 : "$\\eta^\\prime$",
    113 : "$\\rho^0$", 213 : "$\\rho^\\pm$", 223 : "$\\omega$", 333 : "$\\phi$",
    115 : "$a_2^0$", 215 : "$a_2^\pm$", 225 : "$f_2$", 335 : "$f^\\prime_2$",
    20223 : "$f_1$", 20333 : "$f^\\prime_1$", 20113 : "$a^0_1$", 20213 : "$a^\\pm_1$",
    9010221 : "$f_0(980)$", 9000211 : "$a^\\pm_0(980)$",
    
    311 : "$K^0,\\bar{K}^0$", 321 : "$K^\\pm$",  313 : "$K^{*0},\\bar{K}^{*0}$", 323 : "$K^{*\\pm}$",
    315 : "$K^0_2,\\bar{K}^0_2$", 325 : "$K^\\pm_2$",
    411 : "$D^\\pm$", 421 : "$D^0,\\bar{D}^0$", 413: "$D^{*\\pm}$", 415 : "$D^\\pm_2$", 425 : "$D^0_2, \\bar{D}^0_2$", 423: "$D^{*0},\\bar{D}^{*0}$",
    431 : "$D_s^\\pm$", 435 : "$D^\\pm_{s2}$", 20431 : "$D^\\pm_{s1}$", 433 : "$D_s^{*\\pm}$",

    511 : "$B^0,\\bar{B}^0, B^\\pm$", 521 : "$B^\\pm$", 531 : "$B^0_s,\\bar{B}^0_s$", 513 : "$B^*$",515 : "$B^{**}$",
    
    443 : "$J/\\psi$" , 100443 : "$\\psi(2S)$", 553 : "$\Upsilon(1S)$", 20443 : "$\\chi_{c1}(1P)$", 441  : "$\\eta_c$", 100553 : "$\Upsilon(2S)$", 300553 : "$\Upsilon(4S)$",
    
    2212 : "$p,\\bar{p}$", 2224 : "$\\Delta^{++},\\bar{\\Delta}^{--}$",

    3122 : "$\\Lambda^0,\\bar{\\Lambda}^0$",
    3222 : "$\\Sigma^+,\\bar{\Sigma}^-$", "3222B" : "$\\Sigma^\\pm,\\bar{\Sigma}^\\pm$",
    3212 : "$\\Sigma^0,\\bar{\Sigma}^0$", 3112 : "$\\Sigma^-,\\bar{\Sigma}^+$",
    3224 : "$\\Sigma^{*+},\\bar{\Sigma}^{*-}$", "3224B" : "$\\Sigma^{*\\pm},\\bar{\Sigma}^{*\\pm}$",
    3214 : "$\\Sigma^{*0},\\bar{\Sigma}^{*0}$", 3114 : "$\\Sigma^{*-},\\bar{\Sigma}^{*+}$",
    3322 : "$\\Xi^0,\\bar{\\Xi}^0$", 3312 : "$\\Xi^-,\\bar{\\Xi}^+$", 3324 : "$\\Xi^{*0},\\bar{\\Xi}^{*0}$", 3314 : "$\\Xi^{*-},\\bar{\\Xi}^{*+}$",
    3334 : "$\\Omega^-,\\bar{\\Omega}^+$", 3124 : "$\\Lambda^0(1520),\\bar{\\Lambda}^0(1520)$",

    4122 : "$\\Lambda_c^+,\\bar{\\Lambda}^+_c$", 4222 : "$\\Sigma_c^{++}, \\Sigma_c^{0}, \\bar{\\Sigma}_c^{++}, \\bar{\\Sigma}_c^{0}$",
    
    5122 : "$\\Lambda_b^0,\\bar{\\Lambda}^0_b$", 14122 : "$\\Lambda_c(2595)^+,\\bar{\\Lambda}(2595)_c^+$",
    4124 : "$\\Lambda_c(2625)^+,\\bar{\\Lambda}(2595)_c^+$", 4112 : "$\\Sigma_c^0,\\bar{\\Sigma}_c^0$", 4114 : "$\\Sigma_c^{*0},\\bar{\\Sigma}_c^{*0}$",
    4332 : "$\\Omega_c^0,\\bar{\\Omega}_c^0$", 4132 : "$\\Xi_c^0,\\bar{\\Xi}_c^0$"
    }

def writeMisc() :
    global figures
    misc=open(os.path.join(directory,"misc.html"),'w')
    misc.write(header.format(title="Comparisions of Herwig7 and Miscellaneous $e^+e^-$ Data"))
    misc.write("<div style=\"float:none; overflow:auto; width:100%\">\n")
    for val in sorted(figures) :
        name=val.replace("Rivet-EE","").replace(".dat","")
        dat=name[1:] +".dat"
        figName = ("%s/%s" %(directory,dat))
        if(figName not in figures) : continue
        del figures[figures.index(figName)]
        pdf=name[1:] +".pdf"
        png=name[1:] +".png"
        misc.write(plotOutput.format(name=name[1:],pdf=pdf,png=png,dat=dat,energy=""))
        misc.write("\n")
    misc.write("</div>")
    # footer
    misc.write("</body>\n</html>")
    misc.close()
    
def writeQED() :
    global figures
    qed=open(os.path.join(directory,"qed.html"),'w')
    qed.write(header.format(title="Comparisions of Herwig7 and Photon Radiation Data"))
    writePlots2(analyses["QED"],qed)
    # footer
    qed.write("</body>\n</html>")
    qed.close()

def writeTauDecays() :
    global figures
    decays=open(os.path.join(directory,"taus.html"),'w')
    decays.write(header.format(title="Comparisions of Herwig7 and Tau Decay Data"))
    decays.write("<h2 id=\"2pi\">$\\tau\\to\\nu_\\tau\\pi^-\\pi^0$</h2>\n")
    writePlots2(analyses["TauDecays"]["2pi"],decays)
    decays.write("<h2 id=\"Kpi\">$\\tau\\to\\nu_\\tau K\\pi$</h2>\n")
    writePlots2(analyses["TauDecays"]["Kpi"],decays)
    decays.write("<h2 id=\"KK\">$\\tau\\to\\nu_\\tau KK$</h2>\n")
    writePlots2(analyses["TauDecays"]["KK"],decays)

    decays.write("<h2 id=\"lnu\">$\\tau\\to\\nu_\\tau \\ell\\nu$</h2>\n")
    writePlots2(analyses["TauDecays"]["lnu"],decays)
    decays.write("<h2 id=\"Keta\">$\\tau\\to\\nu_\\tau K\\eta$</h2>\n")
    writePlots2(analyses["TauDecays"]["Keta"],decays)
    
    decays.write("<h2 id=\"3pi\">$\\tau\\to\\nu_\\tau\\pi\\pi\\pi$</h2>\n")
    writePlots2(analyses["TauDecays"]["3pi"],decays)
    decays.write("<h2 id=\"Kpipi\">$\\tau\\to\\nu_\\tau K\\pi\\pi$</h2>\n")
    writePlots2(analyses["TauDecays"]["Kpipi"],decays)
    decays.write("<h2 id=\"KKpi\">$\\tau\\to\\nu_\\tau KK\\pi$</h2>\n")
    writePlots2(analyses["TauDecays"]["KKpi"],decays)
    decays.write("<h2 id=\"3K\">$\\tau\\to\\nu_\\tau KKK$</h2>\n")
    writePlots2(analyses["TauDecays"]["3K"],decays)
    decays.write("<h2 id=\"2pieta\">$\\tau\\to\\nu_\\tau\\eta\\pi\\pi$</h2>\n")
    writePlots2(analyses["TauDecays"]["2pieta"],decays)
    decays.write("<h2 id=\"2pigamma\">$\\tau\\to\\nu_\\tau\\gamma\\pi\\pi$</h2>\n")
    writePlots2(analyses["TauDecays"]["2pigamma"],decays)

    
    decays.write("<h2 id=\"4pi\">$\\tau\\to\\nu_\\tau4\\pi$</h2>\n")
    writePlots2(analyses["TauDecays"]["4pi"],decays)
    decays.write("<h2 id=\"5pi\">$\\tau\\to\\nu_\\tau5\\pi$</h2>\n")
    writePlots2(analyses["TauDecays"]["5pi"],decays)
    # footer
    decays.write("</body>\n</html>")
    decays.close()


def writeDecays() :
    global figures
    decays=open(os.path.join(directory,"decays.html"),'w')
    decays.write(header.format(title="Comparisions of Herwig7 and Hadronic Decay Data"))
    # mesons
    decays.write("<h2 id=\"MESONS\">Meson</h2>\n")
    decays.write("<h3 id=\"m_light\">Light, Unflavoured</h3>\n")
    for val in [221,331,223,333,20113,20213] : 
        decays.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots2(analyses["HadronDecays"][val],decays)
        decays.write("</div>\n")
        
    decays.write("<h3 id=\"m_charm\">Charm Mesons</h3>\n")
    for val in [411,421,431] : 
        decays.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots2(analyses["HadronDecays"][val],decays)
        decays.write("</div>\n")
        
    decays.write("<h3 id=\"m_charm\">Bottom Mesons</h3>\n")
    for val in [511,521] : 
        decays.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots2(analyses["HadronDecays"][val],decays)
        decays.write("</div>\n")
        
    decays.write("<h3 id=\"m_ccbar\">Charmonium</h3>\n")
    for val in [441,443] : 
        decays.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots2(analyses["HadronDecays"][val],decays)
        decays.write("</div>\n")
        
    decays.write("<h3 id=\"m_bbbar\">Bottomonium</h3>\n")
    for val in [553,100553,300553] : 
        decays.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots2(analyses["HadronDecays"][val],decays)
        decays.write("</div>\n")
    
    # footer
    decays.write("</body>\n</html>")
    decays.close()
    
def writeMult() :
    global figures
    mult=open(os.path.join(directory,"mult.html"),'w')
    mult.write(header.format(title="Comparisions of Herwig7 and $e^+e^-$ Mult Particle Spectra"))
    # mesons
    mult.write("<h2 id=\"MESONS\">Meson</h2>\n")
    mult.write("<h3 id=\"m_light\">Light, Unflavoured</h3>\n")
    for val in [211,111,221,331,213,113,223,333,225,335,20223,20333,9010221, 9000211] : 
        mult.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots(analyses["Multiplicity"][val],mult)
        mult.write("</div>\n")
    mult.write("<h3 id=\"m_strange\">Strange</h3>\n")
    for val in [311,321,313,323,315,325] :
        mult.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots(analyses["Multiplicity"][val],mult)
        mult.write("</div>\n")
    mult.write("<h3 id=\"m_charm\">Charm</h3>\n")
    for val in [411,421,413,423,431,433,435,20431] : 
        mult.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots(analyses["Multiplicity"][val],mult)
        mult.write("</div>\n")
    mult.write("<h3 id=\"m_bottom\">Bottom</h3>\n")
    for val in [511,521,531,513,515] : 
        mult.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots(analyses["Multiplicity"][val],mult)
        mult.write("</div>\n")
    mult.write("<h3 id=\"m_ccbar\">$c\\bar{c}$</h3>\n")
    for val in [443,100443,20443] : 
        mult.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots(analyses["Multiplicity"][val],mult)
        mult.write("</div>\n")
    mult.write("<h3 id=\"m_bbbar\">$b\\bar{b}$</h3>\n")
    for val in [553] : 
        mult.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots(analyses["Multiplicity"][val],mult)
        mult.write("</div>\n")
    # baryons
    mult.write("<h2 id=\"BARYONS\">Baryons</h2>\n")
    mult.write("<h3 id=\"b_light\">Light, Unflavoured</h3>\n")
    for val in [2212,2224] :
        mult.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots(analyses["Multiplicity"][val],mult)
        mult.write("</div>\n")
    mult.write("<h3 id=\"b_strange\">Hyperons</h3>\n")
    for val in [3122,3222,"3222B",3212,3112,3114,3224,"3224B",3312,3324,3334,3124] :
        mult.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots(analyses["Multiplicity"][val],mult)
        mult.write("</div>\n")
    mult.write("<h3 id=\"b_charm\">Charm Baryons</h3>\n")
    for val in [4122,4222] :
        mult.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots(analyses["Multiplicity"][val],mult)
        mult.write("</div>\n")
    mult.write("<h3 id=\"b_bottom\">Bottom Baryons</h3>\n")
    for val in [5122] :
        mult.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (val,particleNames[val]))
        writePlots(analyses["Multiplicity"][val],mult)
        mult.write("</div>\n")
    # footer
    mult.write("</body>\n</html>")
    mult.close()

# output the identified particle plots
def writeIdentified(index) :
    ident=open(os.path.join(directory,"identified.html"),'w')
    latexNames = { "x"     : "Scaled Momentum/Energy",
                   "xi"    : "Log of Scaled Momentum/Energy",
                   "p"     : "Momentum/Energy",
                   "Ratio" : "Ratios of particle multiplicities",
                   "Other" : "Other distributions" }
    # header for page
    ident.write(header.format(title="Comparisions of Herwig7 and $e^+e^-$ Identified Particle Spectra"))
    # line for index
    index.write("<li> <a href=\"identified.html\">Identified Particle Spectra</a>\n")
    index.write("<ul>\n")
    # photons
    # page
    ident.write("<h2 id=\"gamma\">Photons</h2>\n")
    ident.write("<div style=\"float:none; overflow:auto; \">\n<h3 id=\"%s\">%s</h3>\n" % (22,particleNames[22]))
    # index
    index.write("<li> <a href=\"identified.html#gamma\">Photons:<a/>\n")
    index.write(" <a href=\"identified.html#22\">%s,<a/>\n" % particleNames[22])
    # plots
    for val in ["x","p","xi"] :
        if(len(analyses["IdentifiedParticle"][22][val])==0) : continue
        ident.write("<div style=\"float:none; overflow:auto; \">\n<h3 id=\"%s_%s\">%s</h3>\n" % (22,val,latexNames[val]))
        writePlots(analyses["IdentifiedParticle"][22][val],ident)
    ident.write("</div>\n")
    ident.write("<h2 id=\"MESONS\">Mesons</h2>\n")
    # Light unflavoured mesons
    # page
    ident.write("<h3 id=\"m_light\">Light, Unflavoured</h3>\n")
    # index
    index.write("<li><a href=\"identified.html#m_light\">Light unflavoured mesons:<a/>\n")
    # loop over particles
    for pdgId in [211,111,221,331,213,113,223,333,225,335, 9010221, 9000211] :
        if(pdgId not in analyses["IdentifiedParticle"]) : continue
        sumL = len(analyses["IdentifiedParticle"][pdgId]["x"])+len(analyses["IdentifiedParticle"][pdgId]["p"])+\
               len(analyses["IdentifiedParticle"][pdgId]["xi"])+len(analyses["IdentifiedParticle"][pdgId]["Ratio"])\
               +len(analyses["IdentifiedParticle"][pdgId]["Other"])
        if(sumL==0) : continue
        # lines in html
        ident.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (pdgId,particleNames[pdgId]))
        index.write(" <a href=\"identified.html#%s\">%s,<a/>\n" % (pdgId,particleNames[pdgId]))
        # plots
        for val in ["x","p","xi","Ratio","Other"] :
            if(len(analyses["IdentifiedParticle"][pdgId][val])==0) : continue
            ident.write("<div style=\"float:none; overflow:auto; \">\n<h3 id=\"%s_%s\">%s</h3>\n" % (pdgId,val,latexNames[val]))
            writePlots(analyses["IdentifiedParticle"][pdgId][val],ident)
    # Strange mesons
    # page
    ident.write("<h3 id=\"m_strange\">Strange</h3>\n")
    # index
    index.write("<li><a href=\"identified.html#m_strange\">Strange mesons:<a/>\n")
    # loop over particles
    for pdgId in [311,321,313,323] :
        if(pdgId not in analyses["IdentifiedParticle"]) : continue
        sumL = len(analyses["IdentifiedParticle"][pdgId]["x"])+len(analyses["IdentifiedParticle"][pdgId]["p"])+\
               len(analyses["IdentifiedParticle"][pdgId]["xi"])+len(analyses["IdentifiedParticle"][pdgId]["Ratio"])\
               +len(analyses["IdentifiedParticle"][pdgId]["Other"])
        if(sumL==0) : continue
        # lines in html
        ident.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (pdgId,particleNames[pdgId]))
        index.write(" <a href=\"identified.html#%s\">%s,<a/>\n" % (pdgId,particleNames[pdgId]))
        # plots
        for val in ["x","p","xi","Ratio","Other"] :
            if(len(analyses["IdentifiedParticle"][pdgId][val])==0) : continue
            ident.write("<div style=\"float:none; overflow:auto; \">\n<h3 id=\"%s_%s\">%s</h3>\n" % (pdgId,val,latexNames[val]))
            writePlots(analyses["IdentifiedParticle"][pdgId][val],ident)
    # charm meson
    # page
    ident.write("<h3 id=\"m_charm\">Charm</h3>\n")
    # index
    index.write("<li><a href=\"identified.html#m_charm\">Charm mesons:<a/>\n")
    # loop over particles
    for pdgId in [411,421,413,423,425,431,433] :
        if(pdgId not in analyses["IdentifiedParticle"]) : continue
        sumL = len(analyses["IdentifiedParticle"][pdgId]["x"])+len(analyses["IdentifiedParticle"][pdgId]["p"])+\
               len(analyses["IdentifiedParticle"][pdgId]["xi"])+len(analyses["IdentifiedParticle"][pdgId]["Ratio"])\
               +len(analyses["IdentifiedParticle"][pdgId]["Other"])
        if(sumL==0) : continue
        # lines in html
        ident.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (pdgId,particleNames[pdgId]))
        index.write(" <a href=\"identified.html#%s\">%s,<a/>\n" % (pdgId,particleNames[pdgId]))
        # plots
        for val in ["x","p","xi","Ratio","Other"] :
            if(len(analyses["IdentifiedParticle"][pdgId][val])==0) : continue
            ident.write("<div style=\"float:none; overflow:auto; \">\n<h3 id=\"%s_%s\">%s</h3>\n" % (pdgId,val,latexNames[val]))
            writePlots(analyses["IdentifiedParticle"][pdgId][val],ident)
    # charmonium
    # page
    ident.write("<h3 id=\"m_charm\">Charmonium</h3>\n")
    # index
    index.write("<li><a href=\"identified.html#m_charm\">Charmonium:<a/>\n")
    # loop over particles
    for pdgId in [443] :
        if(pdgId not in analyses["IdentifiedParticle"]) : continue
        sumL = len(analyses["IdentifiedParticle"][pdgId]["x"])+len(analyses["IdentifiedParticle"][pdgId]["p"])+\
               len(analyses["IdentifiedParticle"][pdgId]["xi"])+len(analyses["IdentifiedParticle"][pdgId]["Ratio"])\
               +len(analyses["IdentifiedParticle"][pdgId]["Other"])
        if(sumL==0) : continue
        # lines in html
        ident.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (pdgId,particleNames[pdgId]))
        index.write(" <a href=\"identified.html#%s\">%s,<a/>\n" % (pdgId,particleNames[pdgId]))
        # plots
        for val in ["x","p","xi","Ratio","Other"] :
            if(len(analyses["IdentifiedParticle"][pdgId][val])==0) : continue
            ident.write("<div style=\"float:none; overflow:auto; \">\n<h3 id=\"%s_%s\">%s</h3>\n" % (pdgId,val,latexNames[val]))
            writePlots(analyses["IdentifiedParticle"][pdgId][val],ident)
    # Baryons
    ident.write("<h2 id=\"BARYONS\">Baryons</h2>\n")
    # light baryons
    ident.write("<h3 id=\"b_light\">Light, Unflavoured</h3>\n")
    index.write("<li><a href=\"identified.html#b_light\">Light unflavoured baryons:<a/>\n")
    # loop over particles
    for pdgId in [2212,2224] :
        if(pdgId not in analyses["IdentifiedParticle"]) : continue
        sumL = len(analyses["IdentifiedParticle"][pdgId]["x"])+len(analyses["IdentifiedParticle"][pdgId]["p"])+\
               len(analyses["IdentifiedParticle"][pdgId]["xi"])+len(analyses["IdentifiedParticle"][pdgId]["Ratio"])\
               +len(analyses["IdentifiedParticle"][pdgId]["Other"])
        if(sumL==0) : continue
        # lines in html
        ident.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (pdgId,particleNames[pdgId]))
        index.write(" <a href=\"identified.html#%s\">%s,<a/>\n" % (pdgId,particleNames[pdgId]))
        # plots
        for val in ["x","p","xi","Ratio","Other"] :
            if(len(analyses["IdentifiedParticle"][pdgId][val])==0) : continue
            ident.write("<div style=\"float:none; overflow:auto; \">\n<h3 id=\"%s_%s\">%s</h3>\n" % (pdgId,val,latexNames[val]))
            writePlots(analyses["IdentifiedParticle"][pdgId][val],ident)
    # hyperons
    ident.write("<h3 id=\"b_strange\">Hyperons</h3>\n")
    index.write("<li><a href=\"identified.html#b_strange\">Hyperons:<a/>\n")
    # loop over particles
    for pdgId in [3122,3222,3212,3112,3224,"3224B",3114,3312,3322,3314,3324,3334,3124] :
        if(pdgId not in analyses["IdentifiedParticle"]) : continue
        sumL = len(analyses["IdentifiedParticle"][pdgId]["x"])+len(analyses["IdentifiedParticle"][pdgId]["p"])+\
               len(analyses["IdentifiedParticle"][pdgId]["xi"])+len(analyses["IdentifiedParticle"][pdgId]["Ratio"])\
               +len(analyses["IdentifiedParticle"][pdgId]["Other"])
        if(sumL==0) : continue
        # lines in html
        ident.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (pdgId,particleNames[pdgId]))
        index.write(" <a href=\"identified.html#%s\">%s,<a/>\n" % (pdgId,particleNames[pdgId]))
        # plots
        for val in ["x","p","xi","Ratio","Other"] :
            if(len(analyses["IdentifiedParticle"][pdgId][val])==0) : continue
            ident.write("<div style=\"float:none; overflow:auto; \">\n<h3 id=\"%s_%s\">%s</h3>\n" % (pdgId,val,latexNames[val]))
            writePlots(analyses["IdentifiedParticle"][pdgId][val],ident)
    # charm baryons
    ident.write("<h3 id=\"b_charm\">Charm Baryons</h3>\n")
    index.write("<li><a href=\"identified.html#b_charm\">Charm baryons:<a/>\n")
    # loop over particles
    for pdgId in [4122,4112,4114,4332,4132,14122,4124] :
        if(pdgId not in analyses["IdentifiedParticle"]) : continue
        sumL = len(analyses["IdentifiedParticle"][pdgId]["x"])+len(analyses["IdentifiedParticle"][pdgId]["p"])+\
               len(analyses["IdentifiedParticle"][pdgId]["xi"])+len(analyses["IdentifiedParticle"][pdgId]["Ratio"])\
               +len(analyses["IdentifiedParticle"][pdgId]["Other"])
        if(sumL==0) : continue
        # lines in html
        ident.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"%s\">%s</h4>\n" % (pdgId,particleNames[pdgId]))
        index.write(" <a href=\"identified.html#%s\">%s,<a/>\n" % (pdgId,particleNames[pdgId]))
        # plots
        for val in ["x","p","xi","Ratio","Other"] :
            if(len(analyses["IdentifiedParticle"][pdgId][val])==0) : continue
            ident.write("<div style=\"float:none; overflow:auto; \">\n<h3 id=\"%s_%s\">%s</h3>\n" % (pdgId,val,latexNames[val]))
            writePlots(analyses["IdentifiedParticle"][pdgId][val],ident)
    # bottom fragmentation
    ident.write("<h2 id=\"b_frag\">Bottom Fragmentation Function</h2>\n")
    index.write("<li><a href=\"identified.html#b_frag\">Bottom Fragmentation Function:<a/>\n")
    for val in ["weak","lead","weak_mean","lead_mean"] :
        if val=="weak" :
            name="Weakly Decaying B hadron"
            name2="Weakly Decaying"
        elif val=="lead" :
            name="Leading Decaying B hadron"
            name2="Leading"
        elif val=="weak_mean":
            name="Weakly Decaying B hadron (average)"
            name2="Weakly Decaying (average)"
        elif val=="lead_mean" :
            name="Leading Decaying B hadron (average)"
            name2="Leading (average)"
        ident.write("<div style=\"float:none; overflow:auto; \">\n<h3 id=\"%s\">%s</h3>\n" % (val,name))
        index.write("<a href=\"identified.html#%s\">%s,<a/>\n" % (val,name2))
        writePlots2(analyses["IdentifiedParticle"][511][val],ident)
        ident.write("</div>\n")
    # footer
    index.write(" </ul>\n")    
    ident.write("</body>\n</html>")
    ident.close()
    
def writeFlavour() :
    global figures
    flavour=open(os.path.join(directory,"flavour.html"),'w')
    flavour.write(header.format(title="Comparisions of Herwig7 and Flavour Separated $e^+e^-$ Data"))
    # total multiplicity
    flavour.write("<h2 id=\"mult\">Charged Particle Multiplicity</h2>\n")
    for flav in [1,4,5] :
        if(flav==1) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_1\">Light</h3>\n")
        elif(flav==4) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_4\">Charm</h3>\n")
        elif(flav==5) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_5\">Bottom</h3>\n")
        writePlots(analyses["Charged"]["TotalChargedMult"][flav],flavour)
        flavour.write("</div>\n")
        
    flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_51\">Bottom-Light Difference</h3>\n")
    writePlots(analyses["Charged"]["TotalChargedMult"][51],flavour)
    flavour.write("</div>\n")
    flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_41\">Charm-Light Difference</h3>\n")
    writePlots(analyses["Charged"]["TotalChargedMult"][41],flavour)
    flavour.write("</div>\n")
    # multiplicity dist
    flavour.write("<h2 id=\"multdist\">Charged Particle Multiplicity Distribution</h2>\n")
    for flav in [2,4,5] :
        if(flav==1) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_1\">Light</h3>\n")
        elif(flav==2) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"spectrum_2\">u, d, s, c</h3>\n")
        elif(flav==4) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_4\">Charm</h3>\n")
        elif(flav==5) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_5\">Bottom</h3>\n")
        writePlots(analyses["Charged"]["DistChargedMult"][flav],flavour)
        flavour.write("</div>\n")
    flavour.write("</div>\n")
    # event shapes
    flavour.write("<h2 id=\"EVENT\">Event Shapes</h2>\n")
    # thrust
    flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"thrust\">Thrust</h3>\n")
    for flav in [2,5] :
        if(flav==1) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_1\">Light</h3>\n")
        elif(flav==2) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"spectrum_2\">u, d, s, c</h3>\n")
        elif(flav==4) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_4\">Charm</h3>\n")
        elif(flav==5) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_5\">Bottom</h3>\n")
        writePlots(analyses["EventShapesFlavour"]["T"][flav],flavour)
        flavour.write("</div>\n")
    flavour.write("</div>\n")
    # heavy jet mass
    flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"heavyjetmass\">Heavy Jet Mass</h3>\n")
    for flav in [2,5] :
        if(flav==1) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_1\">Light</h3>\n")
        elif(flav==2) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"spectrum_2\">u, d, s, c</h3>\n")
        elif(flav==4) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_4\">Charm</h3>\n")
        elif(flav==5) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_5\">Bottom</h3>\n")
        writePlots(analyses["EventShapesFlavour"]["HeavyJetMass"][flav],flavour)
        flavour.write("</div>\n")
    flavour.write("</div>\n")
    # BT
    flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"bt\">Total Jet Broadening</h3>\n")
    for flav in [2,5] :
        if(flav==1) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_1\">Light</h3>\n")
        elif(flav==2) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"spectrum_2\">u, d, s, c</h3>\n")
        elif(flav==4) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_4\">Charm</h3>\n")
        elif(flav==5) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_5\">Bottom</h3>\n")
        writePlots(analyses["EventShapesFlavour"]["BT"][flav],flavour)
        flavour.write("</div>\n")
    flavour.write("</div>\n")
    # BW
    flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"bw\">Wide Jet Broadening</h3>\n")
    for flav in [2,5] :
        if(flav==1) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_1\">Light</h3>\n")
        elif(flav==2) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"spectrum_2\">u, d, s, c</h3>\n")
        elif(flav==4) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_4\">Charm</h3>\n")
        elif(flav==5) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_5\">Bottom</h3>\n")
        writePlots(analyses["EventShapesFlavour"]["BW"][flav],flavour)
        flavour.write("</div>\n")
    flavour.write("</div>\n")
    # C
    flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"C\">C</h3>\n")
    for flav in [2,5] :
        if(flav==1) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_1\">Light</h3>\n")
        elif(flav==2) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"spectrum_2\">u, d, s, c</h3>\n")
        elif(flav==4) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_4\">Charm</h3>\n")
        elif(flav==5) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_5\">Bottom</h3>\n")
        writePlots(analyses["EventShapesFlavour"]["C"][flav],flavour)
        flavour.write("</div>\n")
    flavour.write("</div>\n")
    # D
    flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"D\">D</h3>\n")
    for flav in [2,5] :
        if(flav==1) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_1\">Light</h3>\n")
        elif(flav==2) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"spectrum_2\">u, d, s, c</h3>\n")
        elif(flav==4) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_4\">Charm</h3>\n")
        elif(flav==5) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"mult_5\">Bottom</h3>\n")
        writePlots(analyses["EventShapesFlavour"]["D"][flav],flavour)
        flavour.write("</div>\n")
    flavour.write("</div>\n")
    
    
    # spectrum
    flavour.write("<h2 id=\"spectrum\">Charged Particle Spectrum</h2>\n")
    for flav in [1,2,4,5] :
        if(flav==1) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"spectrum_1\">Light</h3>\n")
        elif(flav==2) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"spectrum_2\">u, d, s, c</h3>\n")
        elif(flav==4) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"spectrum_4\">Charm</h3>\n")
        elif(flav==5) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"spectrum_5\">Bottom</h3>\n")
        writePlots(analyses["Charged"]["ChargedSpectrum"][flav],flavour)
        flavour.write("</div>\n")
    # identified particle specrtra
    flavour.write("<h2 id=\"spectrum\">Identified Particle Spectra</h2>\n")
    for flav in [1,4,5,41,51] :
        if(flav==1) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"identf_1\">Light</h3>\n")
        elif(flav==4) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"identf_4\">Charm</h3>\n")
        elif(flav==5) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"identf_5\">Bottom</h3>\n")
        elif(flav==41) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"identf_41\">Ratio Charm to Light</h3>\n")
        elif(flav==51) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"identf_51\">Ratio Bottom to Light</h3>\n")
        for val in [111,211,311,321,313,333,2212,3122,413] :
            if(len(analyses["IdentifiedParticleFlavour"][val][flav])==0) : continue
            flavour.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"f_%s_%s\">%s</h4>\n" % (flav,val,particleNames[val]))
            writePlots(analyses["IdentifiedParticleFlavour"][val][flav],flavour)
            flavour.write("</div>\n")
    # multiplicities
    flavour.write("<h2 id=\"mults\">Identified Particle Multiplicities</h2>\n")
    for flav in [1,4,5,41,51] :
        if(flav==1) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"identm_1\">Light</h3>\n")
        elif(flav==4) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"identm_4\">Charm</h3>\n")
        elif(flav==5) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"identm_5\">Bottom</h3>\n")
        elif(flav==41) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"identm_41\">Charm-Light</h3>\n")
        elif(flav==51) :
            flavour.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"identm_51\">Bottom-Light</h3>\n")
        for val in [211,311,321,313,333,2212,3122] :
            if(len(analyses["MultiplicityFlavour"][val][flav])==0) : continue
            flavour.write("<div style=\"float:none; overflow:auto; \">\n<h4 id=\"f_%s_%s\">%s</h4>\n" % (flav,val,particleNames[val]))
            writePlots(analyses["MultiplicityFlavour"][val][flav],flavour)
            flavour.write("</div>\n")

    
    flavour.write("</body>\n</html>")
    flavour.close()
    
def writeCharged() :
    global figures
    charged=open(os.path.join(directory,"charged.html"),'w')
    charged.write(header.format(title="Comparisions of Herwig7 and $e^+e^-$ Data on Charged Particles"))
    # total multiplicity
    charged.write("<h2 id=\"mult\">Charged Particle Multiplicity</h2>\n")
    charged.write("<div style=\"float:none; overflow:auto; width:100%\">\n")
    writePlots(analyses["Charged"]["TotalChargedMult"][0],charged)
    charged.write("</div>\n")
    # multiplicity dist
    charged.write("<h2 id=\"multdist\">Charged Particle Multiplicity Distribution</h2>\n")
    charged.write("<div style=\"float:none; overflow:auto; width:100%\">\n")
    writePlots(analyses["Charged"]["DistChargedMult"][0],charged)
    charged.write("</div>\n")
    charged.write("<h2 id=\"multdist_c\">Charged Particle Multiplicity Distributions, with cuts</h2>\n")
    charged.write("<div style=\"float:none; overflow:auto; width:100%\">\n")
    writePlots(analyses["Charged"]["DistChargedMult"]["C"],charged)
    charged.write("</div>\n")
    # spectra
    charged.write("<h2 id=\"spectra\">Charged Particle Spectra</h2>\n")
    charged.write("<div style=\"float:none; overflow:auto; width:100%\">\n")
    writePlots(analyses["Charged"]["ChargedSpectrum"][0],charged)
    charged.write("</div>\n")
    # dist w.r.t thrust
    charged.write("<h2 id=\"T\">Charged Particle Distribution w.r.t the Thrust axis</h2>\n")
    # rap
    charged.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"Trap\">Rapidity</h3>\n")
    writePlots(analyses["Charged"]["ChargedRapidityThrust"],charged)
    charged.write("</div>\n")
    # ptin
    charged.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"Tptin\">$p_\\perp^{\\text{in}}$</h3>\n")
    writePlots(analyses["Charged"]["ChargedpTInThrust"],charged)
    charged.write("</div>\n")
    # ptout
    charged.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"Tptout\">$p_\\perp^{\\text{out}}$</h3>\n")
    writePlots(analyses["Charged"]["ChargedpTOutThrust"],charged)
    charged.write("</div>\n")
    charged.write("</body>\n</html>")
    # dist w.r.t sphericity
    charged.write("<h2 id=\"T\">Charged Particle Distribution w.r.t the Sphericity axis</h2>\n")
    # rap
    charged.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"Srap\">Rapidity</h3>\n")
    writePlots(analyses["Charged"]["ChargedRapiditySphericity"],charged)
    charged.write("</div>\n")
    # ptin
    charged.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"Sptin\">$p_\\perp^{\\text{in}}$</h3>\n")
    writePlots(analyses["Charged"]["ChargedpTInSphericity"],charged)
    charged.write("</div>\n")
    # ptout
    charged.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"Sptout\">$p_\\perp^{\\text{out}}$</h3>\n")
    writePlots(analyses["Charged"]["ChargedpTOutSphericity"],charged)
    charged.write("</div>\n")
    charged.write("</body>\n</html>")
    charged.close()
    
def writeJets(index) :
    jets=open(os.path.join(directory,"jets.html"),'w')
    jets.write(header.format(title="Comparisions of Herwig7 and $e^+e^-$ Jet Data"))
    index.write("<li> <a href=\"jets.html\">Jets</a>\n")
    index.write("<ul>\n")
    # fractions
    jets.write("<h2 id=\"FRACTION\">Jet Fractions</h2>\n")
    index.write("<li> <a href=\"jets.html#FRACTION\">Jet Fractions:</a> \n")
    lFormat="""<div style=\"float:none; overflow:auto; width:100%\">\n<{hlevel} id=\"{tag}\">{name}</{hlevel}>\n"""
    for i in range(1,7) :
        if(i<6) : index.write(" <a href=\"jets.html#%sjet\">%s jet,<a/>\n" % (i,i))
        else    :  index.write(" <a href=\"jets.html#%sjet\">%s jet.<a/>\n" % (i,i))
        jets .write(lFormat.format(hlevel="h3",tag   ="%sjet" % i,name  ="%s Jet Fraction" % i))
        jets .write(lFormat.format(hlevel="h4",tag   ="%sjet_dur" % i,name  ="%s Jet Fraction (Durham)" % i))
        writePlots(analyses["EventShapes"]["%sjet_dur" % i],jets)
        jets.write("</div>\n")
        jets .write(lFormat.format(hlevel="h4",tag   ="%sjet_jade" % i,name  ="%s Jet Fraction (JADE)" % i))
        writePlots(analyses["EventShapes"]["%sjet_jade" % i],jets)
        jets.write("</div>\n")
        jets.write("</div>\n")
    # differential jet rates
    index.write("<li> <a href=\"jets.html#DRATE\">Differential jet rates:</a> \n")
    jets.write("<h2 id=\"DRATE\">Differential Jet Rates</h2>\n")
    for i in range(1,6) :
        yval="%s%s" % (i,i+1)
        if(i!=5) : index.write(" <a href=\"jets.html#y%s\">$y_{%s}$,<a/>\n" %(yval,yval))
        else     : index.write(" <a href=\"jets.html#y%s\">$y_{%s}$.<a/>\n" %(yval,yval))
        jets .write(lFormat.format(hlevel="h3",tag   ="y%s" % yval,name  ="$y_{%s}$" % yval))
        if(len(analyses["EventShapes"]["y%s_dur" % yval])!=0) :
            jets .write(lFormat.format(hlevel="h4",tag   ="y%s_dur" % yval,name  ="$y_{%s}$ (Durham)" % yval))
            writePlots(analyses["EventShapes"]["y%s_dur" % yval],jets)
            jets.write("</div>\n")
        if(len(analyses["EventShapes"]["y%s_jade" % yval])!=0) :
            jets .write(lFormat.format(hlevel="h4",tag   ="y%s_jade" % yval,name  ="$y_{%s}$ (JADE)" % yval))
            writePlots(analyses["EventShapes"]["y%s_jade" % yval],jets)
            jets.write("</div>\n")


        jets.write("</div>\n")
    # footer
    index.write(" </ul>\n")
    jets.write("</div>\n")
    jets.write("</body>\n</html>")
    jets.close()

def writeEventShapes() :
    global figures
    event=open(os.path.join(directory,"event.html"),'w')
    event.write(header.format(title="Comparisions of Herwig7 and $e^+e^-$ Event Shape Data"))
    # thrust and related
    event.write("<h2 id=\"THRUST\">Thrust and Related Variables</h2>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"thrust\">Thrust</h3>\n")
    writePlots(analyses["EventShapes"]["T"],event)
    event.write("</div>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"major\">Thrust Major</h3>\n")
    writePlots(analyses["EventShapes"]["Major"],event)
    event.write("</div>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"minor\">Thrust Minor</h3>\n")
    writePlots(analyses["EventShapes"]["Minor"],event)
    event.write("</div>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"oblateness\">Oblateness</h3>\n")
    writePlots(analyses["EventShapes"]["O"],event)
    event.write("</div>\n")
    # sphericity and related
    event.write("<h2 id=\"SPHERICITY\">Sphericity and Related Variables</h2>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"sphericity\">Sphericity</h3>\n")
    writePlots(analyses["EventShapes"]["S"],event)
    event.write("</div>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"planarity\">Planarity</h3>\n")
    writePlots(analyses["EventShapes"]["P"],event)
    event.write("</div>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"aplanarity\">Aplanarity</h3>\n")
    writePlots(analyses["EventShapes"]["A"],event)
    event.write("</div>\n")
    # jet masses
    event.write("<h2 id=\"MASSES\">Jet Masses</h2>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"heavyjetmass\">Heavy Jet Mass</h3>\n")
    writePlots(analyses["EventShapes"]["HeavyJetMass"],event)
    event.write("</div>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"lightjetmass\">Light Jet Mass</h3>\n")
    writePlots(analyses["EventShapes"]["LightJetMass"],event)
    event.write("</div>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"totaljetmass\">Total Jet Mass</h3>\n")
    writePlots(analyses["EventShapes"]["TotalJetMass"],event)
    event.write("</div>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"jetmassdifference\">Jet Mass Difference</h3>\n")
    writePlots(analyses["EventShapes"]["JetMassDifference"],event)
    event.write("</div>\n")
    event.write("<h2 id=\"EEC\">Energy-Energy Correlations</h2>\n")
    # EEC and AEEC
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"eec\">Energy-Energy Correlation</h3>\n")
    writePlots(analyses["EventShapes"]["EEC"],event)
    event.write("</div>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"aeec\">Asymmetry of the Energy-Energy Correlation</h3>\n")
    writePlots(analyses["EventShapes"]["AEEC"],event)
    event.write("</div>\n")
    # jet broadening
    event.write("<h2 id=\"JB\">Jet Broadenings</h2>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"bt\">Total Jet Broadening</h3>\n")
    writePlots(analyses["EventShapes"]["BT"],event)
    event.write("</div>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"bw\">Wide Jet Broadening</h3>\n")
    writePlots(analyses["EventShapes"]["BW"],event)
    event.write("</div>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"bn\">Narrow Jet Broadening</h3>\n")
    writePlots(analyses["EventShapes"]["BN"],event)
    event.write("</div>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"bdiff\">Difference of Jet Broadenings</h3>\n")
    writePlots(analyses["EventShapes"]["Bdiff"],event)
    event.write("</div>\n")
    # C and D
    event.write("<h2 id=\"CD\">C and D parameters</h2>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"C\">C</h3>\n")
    writePlots(analyses["EventShapes"]["C"],event)
    event.write("</div>\n")
    event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"D\">D</h3>\n")
    writePlots(analyses["EventShapes"]["D"],event)
    event.write("</div>\n")
    # moments of event shapes
    event.write("<h2 id=\"MOMENT\">Moments of Event Shapes</h2>\n")
    for val in ["Moment_T","Moment_M","Moment_m","Moment_O","Moment_H","Moment_L","Moment_BT","Moment_BW","Moment_BN","Moment_C","Moment_S","Moment_y"] :
        if(val=="Moment_T") :
            event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"m_thrust\">Thrust</h3>\n")
        elif(val=="Moment_M") :
            event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"m_major\">Thrust Major</h3>\n")
        elif(val=="Moment_m") :
            event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"m_minor\">Thrust Minor</h3>\n")
        elif(val=="Moment_O") :
            event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"m_O\">Oblateness</h3>\n")
        elif(val=="Moment_H") :
            event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"m_heavy\">Heavy Jet Mass</h3>\n")
        elif(val=="Moment_L") :
            event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"m_light\">Light Jet Mass</h3>\n")
        elif(val=="Moment_BT") :
            event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"m_bt\">Total Jet Broadening</h3>\n")
        elif(val=="Moment_BW") :
            event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"m_bw\">Wide Jet Broadening</h3>\n")
        elif(val=="Moment_BN") :
            event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"m_bn\">Narrow Jet Broadening</h3>\n")
        elif(val=="Moment_C") :
            event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"m_C\">C-parameter</h3>\n")
        elif(val=="Moment_S") :
            event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"m_S\">Sphericity</h3>\n")
        elif(val=="Moment_y") :
            event.write("<div style=\"float:none; overflow:auto; width:100%\">\n<h3 id=\"m_y\">$y_{23}$</h3>\n")
        writePlots(analyses["EventShapes"][val],event)
        event.write("</div>\n")
    event.write("</body>\n</html>")
    event.close()

#  ,"Moment_C":{},"Moment_S":{} "Moment_y":{},
    
def writeGluon() :
    global figures
    gluon=open(os.path.join(directory,"gluon.html"),'w')
    gluon.write(header.format(title="Comparisions of Herwig7 and $e^+e^-$ Data on Gluon Jets"))
    # charged particles dists
    gluon.write("<h2 id=\"multdist\">Charged Particle Multiplicity in Gluon Jets</h2>\n")
    gluon.write("<div style=\"float:none; overflow:auto; width:100%\">\n")
    writePlots(analyses["Charged"]["DistChargedMult"][21],gluon)
    gluon.write("</div>\n")
    
    # footer
    gluon.write("</body>\n</html>")
    gluon.close()
    
print 'Total no of figures',len(figures)
# output the event shapes    
writeEventShapes()
writeCharged()
writeMult()
writeFlavour()
writeDecays()
writeTauDecays()
writeQED()
writeGluon()


index=open(os.path.join(directory,"herwig.html"),'w')
index.write(header.format(title="Comparisions of Herwig7 and $e^+e^-$ Data"))
# event shapes
index.write("<ul>\n")
index.write("<li> <a href=\"event.html\">Event Shapes</a>\n")
index.write("<ul>\n")
index.write("<li> <a href=\"event.html#THRUST\">Thrust related:</a> \n")
index.write(" <a href=\"event.html#thrust\">thrust,<a/>\n")
index.write(" <a href=\"event.html#major\">major,<a/>\n")
index.write(" <a href=\"event.html#minor\">minor,<a/>\n")
index.write(" <a href=\"event.html#oblateness\">oblateness.<a/>\n")
index.write("<li> <a href=\"event.html#SPHERICITY\">Sphericity related:</a> \n")
index.write(" <a href=\"event.html#sphericity\">sphericity,<a/>\n")
index.write(" <a href=\"event.html#planarity\">planarity,<a/>\n")
index.write(" <a href=\"event.html#aplanarity\">aplanarity.<a/>\n")
index.write("<li> <a href=\"event.html#MASSES\">Jet Masses:</a> \n")
index.write(" <a href=\"event.html#heavyjetmass\">heavy jet mass,<a/>\n")
index.write(" <a href=\"event.html#lightjetmass\">light jet mass,<a/>\n")
index.write(" <a href=\"event.html#totaljetmass\">total jet mass,<a/>\n")
index.write(" <a href=\"event.html#jetmassdifference\">jet mass difference,<a/>\n")
index.write("<li> <a href=\"event.html#EEC\">Energy-Energy Correlations:</a> \n")
index.write(" <a href=\"event.html#eec\">energy-energy correlation,<a/>\n")
index.write(" <a href=\"event.html#aeec\">asymmetry.<a/>\n")
index.write("<li> <a href=\"event.html#JB\">Jet Broadening:</a> \n")
index.write(" <a href=\"event.html#bt\">total jet broadening,<a/>\n")
index.write(" <a href=\"event.html#bw\">wide jet broadening,<a/>\n")
index.write(" <a href=\"event.html#bn\">narrow jet broadening,<a/>\n")
index.write(" <a href=\"event.html#bdiff\">difference of jet broadenings.<a/>\n")
index.write("<li> <a href=\"event.html#CD\">C and D parameters:</a> \n")
index.write(" <a href=\"event.html#C\">C,<a/>\n")
index.write(" <a href=\"event.html#D\">D.<a/>\n")
index.write(" </ul>\n")
# charged particles
index.write("<li> <a href=\"charged.html\">Charged Particles</a>\n")
index.write("<ul>\n")
index.write("<li> <a href=\"charged.html#mult\">Total Charged Multiplicity</a> \n")
index.write("<li> <a href=\"charged.html#multdist\">Charged Multiplicity Distribution</a> \n")
index.write("<li> <a href=\"charged.html#spectra\">Charged Spectra</a> \n")
index.write("<li> <a href=\"charged.html#T\">Distributions w.r.t thrust axis:</a> \n")
index.write(" <a href=\"charged.html#Trap\">rapidity,<a/>\n")
index.write(" <a href=\"charged.html#Tptin\">$p_\\perp^{\\text{in}}$,<a/>\n")
index.write(" <a href=\"charged.html#Tptout\">$p_\\perp^{\\text{out}}$.<a/>\n")
index.write("<li> <a href=\"charged.html#T\">Distributions w.r.t sphericity axis:</a> \n")
index.write(" <a href=\"charged.html#Srap\">rapidity,<a/>\n")
index.write(" <a href=\"charged.html#Sptin\">$p_\\perp^{\\text{in}}$,<a/>\n")
index.write(" <a href=\"charged.html#Sptout\">$p_\\perp^{\\text{out}}$.<a/>\n")
index.write(" </ul>\n")
# jets
writeJets(index)
# identified particle spectra
writeIdentified(index)
# identified particle multiplicity
# mesons
index.write("<li> <a href=\"mult.html\">Identified Particle Multiplicities</a>\n")
index.write("<ul>\n")
index.write("<li><a href=\"mult.html#m_light\">Light unflavoured mesons:<a/>\n")
for val in [211,111,221,331,213,113,223,333,225,335,20223,20333,9010221, 9000211] : 
    index.write(" <a href=\"mult.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
index.write("<li><a href=\"mult.html#m_strange\">Strange mesons:<a/>\n")
for val in [311,321,313,323,315,325] :
    index.write(" <a href=\"mult.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
index.write("<li><a href=\"mult.html#m_charm\">Charm mesons:<a/>\n")
for val in [411,421,413,423,431,433,435,20431] : 
    index.write(" <a href=\"mult.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
index.write("<li><a href=\"mult.html#m_bottom\">Bottom mesons:<a/>\n")
for val in [511,521,531,513,515] : 
    index.write(" <a href=\"mult.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
index.write("<li><a href=\"mult.html#m_ccbar\">$c\\bar{c}$ mesons:<a/>\n")
for val in  [443,100443,20443] : 
    index.write(" <a href=\"mult.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
index.write("<li><a href=\"mult.html#m_bbbar\">$b\\bar{b}$ mesons:<a/>\n")
for val in  [553] : 
    index.write(" <a href=\"mult.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
# baryons
index.write("<li><a href=\"mult.html#b_light\">Light unflavoured baryons:<a/>\n")
for val in [2212,2224] :
    index.write(" <a href=\"mult.html#%s\">%s,<a/>\n" % (val,particleNames[val]))    
index.write("<li><a href=\"mult.html#b_strange\">Hyperons:<a/>\n")
for val in [3122,3222,"3222B",3212,3112,3114,3224,"3224B",3312,3324,3334,3124] :
    index.write(" <a href=\"mult.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
index.write("<li><a href=\"mult.html#b_charm\">Charm:<a/>\n")
for val in [4122,4222] :
    index.write(" <a href=\"mult.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
index.write("<li><a href=\"mult.html#b_bottom\">Bottom:<a/>\n")
for val in  [5122] :
    index.write(" <a href=\"mult.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
index.write(" </ul>\n")
# flavour
index.write("<li> <a href=\"flavour.html\">Flavour Separated</a>\n")
index.write(" <ul>\n")

index.write("<li> <a href=\"flavour.html#mult\">Total Charged Multiplicity: </a> \n")
index.write(" <a href=\"flavour.html#mult_1\">light,<a/>\n")
index.write(" <a href=\"flavour.html#mult_4\">charm,<a/>\n")
index.write(" <a href=\"flavour.html#mult_5\">bottom.<a/>\n")

index.write("<li> <a href=\"flavour.html#spectrum\">Charged Particle Spectrum: </a> \n")
index.write(" <a href=\"flavour.html#spectrum_1\">light,<a/>\n")
index.write(" <a href=\"flavour.html#spectrum_4\">charm,<a/>\n")
index.write(" <a href=\"flavour.html#spectrum_5\">bottom.<a/>\n")

index.write("</ul>\n")
# hadron decays
index.write("<li> <a href=\"decays.html\">Hadron Decays</a>\n")
index.write(" <ul>\n")
index.write("<li><a href=\"decays.html#m_light\">Light unflavoured mesons:<a/>\n")
for val in [221,331,223,333] : 
    index.write(" <a href=\"decays.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
index.write("<li><a href=\"decays.html#m_charm\">Charm mesons:<a/>\n")
for val in [411,421,431] : 
    index.write(" <a href=\"decays.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
index.write("<li><a href=\"decays.html#m_charm\">Bottom mesons:<a/>\n")
for val in [511,521] : 
    index.write(" <a href=\"decays.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
index.write("<li><a href=\"decays.html#m_ccbar\">Charmonium:<a/>\n")
for val in [441,443] : 
    index.write(" <a href=\"decays.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
index.write("<li><a href=\"decays.html#m_bbbar\">Bottomonium:<a/>\n")
for val in [553,100553,300553] : 
    index.write(" <a href=\"decays.html#%s\">%s,<a/>\n" % (val,particleNames[val]))
index.write(" </ul>\n")
# tau decays
index.write("<li> <a href=\"taus.html\">Tau Decays</a>\n")

print 'Unused figures',len(figures)
writeMisc()

    # decays.write("<h2 id=\"2pi\">$\\tau\\to\\nu_\\tau\\pi^-\\pi^0$</h2>\n")
    # writePlots2(analyses["TauDecays"]["2pi"],decays)
    # decays.write("<h2 id=\"Kpi\">$\\tau\\to\\nu_\\tau K\\pi$</h2>\n")
    # writePlots2(analyses["TauDecays"]["Kpi"],decays)
    # decays.write("<h2 id=\"KK\">$\\tau\\to\\nu_\\tau KK$</h2>\n")
    # writePlots2(analyses["TauDecays"]["KK"],decays)
    # decays.write("<h2 id=\"3pi\">$\\tau\\to\\nu_\\tau\\pi\\pi\\pi$</h2>\n")
    # writePlots2(analyses["TauDecays"]["3pi"],decays)
    # decays.write("<h2 id=\"Kpipi\">$\\tau\\to\\nu_\\tau K\\pi\\pi$</h2>\n")
    # writePlots2(analyses["TauDecays"]["Kpipi"],decays)
    # decays.write("<h2 id=\"KKpi\">$\\tau\\to\\nu_\\tau KK\\pi$</h2>\n")
    # writePlots2(analyses["TauDecays"]["KKpi"],decays)
    # decays.write("<h2 id=\"3K\">$\\tau\\to\\nu_\\tau KKK$</h2>\n")
    # writePlots2(analyses["TauDecays"]["3K"],decays)
# qed
index.write("<li> <a href=\"qed.html\">QED Radiation</a>\n")
index.write("<li> <a href=\"gluon.html\">Gluons</a>\n")
index.write("<li> <a href=\"misc.html\">Other Plots</a>\n")
# footer
index.write("</ul>\n")
index.write("</body>\n</html>")
index.close()
