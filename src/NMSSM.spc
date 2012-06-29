# NMHDECAY OUTPUT IN SLHA FORMAT
# Info about spectrum calculator
BLOCK SPINFO        # Program information
     1   NMHDECAY   # spectrum calculator
     2   3          # version number
# Input parameters
BLOCK MODSEL
    3     1         # NMSSM PARTICLE CONTENT
BLOCK SMINPUTS
     1     1.27920000E+02   # ALPHA_EM^-1(MZ)
     2     1.16639000E-05   # GF
     3     1.17200000E-01   # ALPHA_S(MZ)
     4     9.11870000E+01   # MZ
     5     4.21400000E+00   # MB(MB)
     6     1.71400000E+02   # MTOP (POLE MASS)
     7     1.77700000E+00   # MTAU
# SMINPUTS Beyond SLHA:
# MW:     0.80420000E+02
# MS:     0.19000000E+00
# MC:     0.14000000E+01
# VUS:     0.22000000E+00
# VCB:     0.40000000E-01
# VUB:     0.40000000E-02
BLOCK MINPAR
     3     2.00000000E+00   # TANBETA
BLOCK EXTPAR
     1     1.50000000E+02   # M1
     2     3.00000000E+02   # M2
     3     1.00000000E+03   # M3
    11     2.50000000E+03   # ATOP
    12     2.50000000E+03   # ABOTTOM
    13     2.50000000E+03   # ATAU
   124     1.32099188E+03   # MA AT Q_STSB
    31     1.00000000E+06   # LEFT SELECTRON
    32     1.00000000E+06   # LEFT SMUON
    33     1.00000000E+06   # LEFT STAU
    34     1.00000000E+06   # RIGHT SELECTRON
    35     1.00000000E+06   # RIGHT SMUON
    36     1.00000000E+06   # RIGHT STAU
    41     1.00000000E+06   # LEFT 1ST GEN. SQUARKS
    42     1.00000000E+06   # LEFT 2ND GEN. SQUARKS
    43     1.00000000E+06   # LEFT 3RD GEN. SQUARKS
    44     1.00000000E+06   # RIGHT U-SQUARKS
    45     1.00000000E+06   # RIGHT C-SQUARKS
    46     1.00000000E+06   # RIGHT T-SQUARKS
    47     1.00000000E+06   # RIGHT D-SQUARKS
    48     1.00000000E+06   # RIGHT S-SQUARKS
    49     1.00000000E+06   # RIGHT B-SQUARKS
    61     7.02000000E-01   # LAMBDA
    62     4.90000000E-02   # KAPPA
    63     1.28000000E+03   # A_LAMBDA
    64     1.00000000E+01   # A_KAPPA
    65     5.30000000E+02   # MU_EFF
# 
BLOCK MASS   # Mass spectrum 
#  PDG Code     mass             particle 
        25     1.39836374E+02   # lightest neutral scalar
        35     1.51443222E+02   # second neutral scalar
        45     1.30169034E+03   # third neutral scalar
        36     4.89195961E+01   # lightest pseudoscalar
        46     1.30475321E+03   # second pseudoscalar
        37     1.30252610E+03   # charged Higgs
   1000001     1.03794721E+03   #  ~d_L
   2000001     1.03715965E+03   #  ~d_R
   1000002     1.03619464E+03   #  ~u_L
   2000002     1.03663055E+03   #  ~u_R
   1000003     1.03794721E+03   #  ~s_L
   2000003     1.03715965E+03   #  ~s_R
   1000004     1.03619464E+03   #  ~c_L
   2000004     1.03663055E+03   #  ~c_R
   1000005     1.03589101E+03   #  ~b_1
   2000005     1.03921740E+03   #  ~b_2
   1000006     8.77891045E+02   #  ~t_1
   2000006     1.18786303E+03   #  ~t_2
   1000011     1.00063432E+03   #  ~e_L
   2000011     1.00054844E+03   #  ~e_R
   1000012     9.98816192E+02   #  ~nue_L
   1000013     1.00063432E+03   #  ~mu_L
   2000013     1.00054844E+03   #  ~mu_R
   1000014     9.98816192E+02   #  ~numu_L
   1000015     9.99312737E+02   #  ~tau_1
   2000015     1.00187154E+03   #  ~tau_2
   1000016     9.98816192E+02   #  ~nutau_L
   1000021     1.07855652E+03   #  ~g
   1000022     9.12928865E+01   # neutralino(1)
   1000023     1.46429979E+02   # neutralino(2)
   1000025     2.94547359E+02   # neutralino(3)
   1000035    -5.56163202E+02   # neutralino(4)
   1000045     5.64325430E+02   # neutralino(5)
   1000024     2.92959564E+02   # chargino(1)
   1000037     5.58677258E+02   # chargino(2)
# Low energy observables
BLOCK LOWEN
# Exp. 2 Sigma: 3.07E-04 < BR(b -> s gamma) < 4.07E-04:
         1     3.41523755E-04   # BR(b -> s gamma)
        11     3.72101862E-04   # (BR(b -> s gamma)+Theor.Err.)
        12     2.91032383E-04   # (BR(b -> s gamma)-Theor.Err.)
# Exp. 2 Sigma: 4.99E-01 < Delta M_d < 5.15E-01:
         2     6.07354601E-01   # Delta M_d in ps^-1
        21     1.05776810E+00   # Delta M_d +Theor.Err.
        22     1.63385106E-01   # Delta M_d -Theor.Err.
# Exp. 2 Sigma: 1.753E+01 < Delta Ms < 1.801E+01:
         3     2.10452493E+01   # Delta M_s in ps^-1
        31     2.76961825E+01   # Delta M_s +Theor.Err.
        32     1.44881882E+01   # Delta M_s -Theor.Err.
# Exp. 95% C.L.: BR(Bs->mu+mu-) < 5.8E-08:
         4     3.54098427E-09   # BR(Bs -> mu+mu-)
        41     6.01429512E-09   # BR(Bs -> mu+mu-)+Theor.Err.
        42     1.71891137E-09   # BR(Bs -> mu+mu-)-Theor.Err.
# Exp. 2 Sigma: 3.40E-05 < BR(B+ > tau+ + nu_tau) < 2.30E-04:
         5     1.31775781E-04   # BR(B+ -> tau+ + nu_tau)
        51     2.63603107E-04   # BR(B+ -> tau+ + nu_tau) + Theor.Err.
        52     5.68712649E-05   # BR(B+ -> tau+ + nu_tau) - Theor.Err.
# 
BLOCK HMIX AT Q=  1.00000000E+03 # (STOP/SBOTTOM MASSES)
     1     5.30000000E+02   # MU_EFF
     2     2.00000402E+00   # TAN(BETA)
     3     1.71023505E+02   # V(Q)
     4     1.74501955E+06   # MA**2
# 3*3 Higgs mixing
BLOCK NMHMIX
  1  1     4.49872460E-01   # S_(1,1)
  1  2     8.91681851E-01   # S_(1,2)
  1  3     5.01821419E-02   # S_(1,3)
  2  1     2.88136045E-02   # S_(2,1)
  2  2    -7.06510966E-02   # S_(2,2)
  2  3     9.97084850E-01   # S_(2,3)
  3  1     8.92627888E-01   # S_(3,1)
  3  2    -4.47115086E-01   # S_(3,2)
  3  3    -5.74765508E-02   # S_(3,3)
# 3*3 Pseudoscalar Higgs mixing
BLOCK NMAMIX
  1  1    -7.53734413E-02   # P_(1,1)
  1  2    -3.76867207E-02   # P_(1,2)
  1  3     9.96442951E-01   # P_(1,3)
  2  1     8.91245670E-01   # P_(2,1)
  2  2     4.45622835E-01   # P_(2,2)
  2  3     8.42700692E-02   # P_(2,3)
# 3rd generation sfermion mixing
BLOCK STOPMIX  # Stop mixing matrix
  1  1    -7.07590438E-01   # Rst_(1,1)
  1  2     7.06622794E-01   # Rst_(1,2)
  2  1    -7.06622794E-01   # Rst_(2,1)
  2  2    -7.07590438E-01   # Rst_(2,2)
BLOCK SBOTMIX  # Sbottom mixing matrix
  1  1    -6.22498693E-01   # Rsb_(1,1)
  1  2     7.82620839E-01   # Rsb_(1,2)
  2  1    -7.82620839E-01   # Rsb_(2,1)
  2  2    -6.22498693E-01   # Rsb_(2,2)
BLOCK STAUMIX  # Stau mixing matrix
  1  1    -6.95140010E-01   # Rsl_(1,1)
  1  2     7.18874374E-01   # Rsl_(1,2)
  2  1    -7.18874374E-01   # Rsl_(2,1)
  2  2    -6.95140010E-01   # Rsl_(2,2)
# Gaugino-Higgsino mixing
BLOCK NMNMIX  # 5*5 Neutralino Mixing Matrix
  1  1     1.13458283E-01   # N_(1,1)
  1  2    -5.54357071E-02   # N_(1,2)
  1  3    -4.98782346E-02   # N_(1,3)
  1  4    -1.98077653E-01   # N_(1,4)
  1  5     9.70737609E-01   # N_(1,5)
  2  1     9.86083701E-01   # N_(2,1)
  2  2    -3.84507802E-02   # N_(2,2)
  2  3     9.97300301E-02   # N_(2,3)
  2  4    -4.07516099E-02   # N_(2,4)
  2  5    -1.20638712E-01   # N_(2,5)
  3  1     7.54509718E-02   # N_(3,1)
  3  2     9.53764002E-01   # N_(3,2)
  3  3    -2.24716204E-01   # N_(3,3)
  3  4     1.71385532E-01   # N_(3,4)
  3  5     6.90724660E-02   # N_(3,5)
  4  1    -1.87054069E-02   # N_(4,1)
  4  2     2.83352907E-02   # N_(4,2)
  4  3     6.97656324E-01   # N_(4,3)
  4  4     6.92377943E-01   # N_(4,4)
  4  5     1.80929977E-01   # N_(4,5)
  5  1     9.33992005E-02   # N_(5,1)
  5  2    -2.91512880E-01   # N_(5,2)
  5  3    -6.71077018E-01   # N_(5,3)
  5  4     6.71076995E-01   # N_(5,4)
  5  5     7.48874698E-02   # N_(5,5)
BLOCK UMIX  # Chargino U Mixing Matrix
  1  1     9.44173996E-01   # U_(1,1)
  1  2    -3.29447211E-01   # U_(1,2)
  2  1     3.29447211E-01   # U_(2,1)
  2  2     9.44173996E-01   # U_(2,2)
BLOCK VMIX  # Chargino V Mixing Matrix
  1  1     9.62164136E-01   # V_(1,1)
  1  2    -2.72470503E-01   # V_(1,2)
  2  1     2.72470503E-01   # V_(2,1)
  2  2     9.62164136E-01   # V_(2,2)
# Higgs reduced couplings
# (as compared to a SM Higgs with same mass)
BLOCK REDCOUP
# H1
  1  1     9.96930616E-01   # U-type fermions
  1  2     1.00594540E+00   # D-type fermions
  1  3     9.98733573E-01   # W,Z bosons
  1  4     9.82183861E-01   # Gluons
  1  5     1.01232522E+00   # Photons
# H2
  2  1    -7.89903273E-02   # U-type fermions
  2  2     6.44291783E-02   # D-type fermions
  2  3    -5.03064262E-02   # W,Z bosons
  2  4     8.42729081E-02   # Gluons
  2  5     8.53137747E-02   # Photons
# H3
  3  1    -4.99889863E-01   # U-type fermions
  3  2     1.99597664E+00   # D-type fermions
  3  3    -7.16563065E-04   # W,Z bosons
  3  4     4.97960461E-01   # Gluons
  3  5     3.79551399E-01   # Photons
# A1
  4  1    -4.21350346E-02   # U-type fermions
  4  2    -1.68540138E-01   # D-type fermions
  4  3     0.00000000E+00   # W,Z bosons
  4  4     6.09163550E-02   # Gluons
  4  5     1.54258690E-01   # Photons
# A2
  5  1     4.98221476E-01   # U-type fermions
  5  2     1.99288590E+00   # D-type fermions
  5  3     0.00000000E+00   # W,Z bosons
  5  4     5.02901556E-01   # Gluons
  5  5     6.47453065E-01   # Photons
# 
# GAUGE AND YUKAWA COUPLINGS AT THE SUSY SCALE
BLOCK GAUGE Q=  1.00000000E+03 # (SUSY SCALE)
         1     3.62224894E-01   # g1(Q,DR_bar)
         2     6.42793227E-01   # g2(Q,DR_bar)
         3     1.05981813E+00   # g3(Q,DR_bar)
BLOCK YU Q=  1.00000000E+03 # (SUSY SCALE)
  3  3     9.51283729E-01   # HTOP(Q,DR_bar)
BLOCK YD Q=  1.00000000E+03 # (SUSY SCALE)
  3  3     3.09237206E-02   # HBOT(Q,DR_bar)
BLOCK YE Q=  1.00000000E+03 # (SUSY SCALE)
  3  3     2.23241912E-02   # HTAU(Q,DR_bar)
BLOCK L/K Q=  1.00000000E+03 # (SUSY SCALE)
         1     7.02000000E-01   # LAMBDA(Q,DR_bar)
         2     4.90000000E-02   # KAPPA(Q,DR_bar)
# 
# SOFT TRILINEAR COUPLINGS AT THE SUSY SCALE
BLOCK AU Q=  1.00000000E+03 # (SUSY SCALE)
  3  3     2.50000000E+03   # ATOP
BLOCK AD Q=  1.00000000E+03 # (SUSY SCALE)
  3  3     2.50000000E+03   # ABOT
BLOCK AE Q=  1.00000000E+03 # (SUSY SCALE)
  3  3     2.50000000E+03   # ATAU
BLOCK AL/AK Q=  1.00000000E+03 # (SUSY SCALE)
         1     1.28000000E+03   # ALAMBDA
         2     1.00000000E+01   # AKAPPA
# 
# SOFT MASSES AT THE SUSY SCALE
BLOCK MSOFT Q=  1.00000000E+03 # (SUSY SCALE)
         1     1.50000000E+02   # M1
         2     3.00000000E+02   # M2
         3     1.00000000E+03   # M3
        21     1.10587466E+06   # M_HD^2
        22     9.23131746E+04   # M_HU^2
        23    -2.79072679E+03   # M_S^2
        31     1.00000000E+06   # M_eL^2
        32     1.00000000E+06   # M_muL^2
        33     1.00000000E+06   # M_tauL^2
        34     1.00000000E+06   # M_eR^2
        35     1.00000000E+06   # M_muR^2
        36     1.00000000E+06   # M_tauR^2
        41     1.00000000E+06   # M_q1L^2
        42     1.00000000E+06   # M_q2L^2
        43     1.00000000E+06   # M_q3L^2
        44     1.00000000E+06   # M_uR^2
        45     1.00000000E+06   # M_cR^2
        46     1.00000000E+06   # M_tR^2
        47     1.00000000E+06   # M_dR^2
        48     1.00000000E+06   # M_sR^2
        49     1.00000000E+06   # M_bR^2
# 
# GAUGE AND YUKAWA COUPLINGS AT THE GUT SCALE
BLOCK GAUGEGUT MGUT=  2.10086949E+16 # (GUT SCALE)
         1     7.10170399E-01   # g1(MGUT,DR_bar), GUT normalization
         2     7.10170401E-01   # g2(MGUT,DR_bar)
         3     7.02117834E-01   # g3(MGUT,DR_bar)
BLOCK YUGUT MGUT=  2.10086949E+16 # (GUT SCALE)
  3  3     1.38184598E+00   # HTOP(MGUT,DR_bar)
BLOCK YDGUT MGUT=  2.10086949E+16 # (GUT SCALE)
  3  3     1.80681615E-02   # HBOT(MGUT,DR_bar)
BLOCK YEGUT MGUT=  2.10086949E+16 # (GUT SCALE)
  3  3     2.20124415E-02   # HTAU(MGUT,DR_bar)
BLOCK LGUT/KGUT MGUT=  2.10086949E+16 # (GUT SCALE)
         1     3.48198696E+00   # LAMBDA(MGUT,DR_bar)
         2     4.18288468E-01   # KAPPA(MGUT,DR_bar)
# 
# SOFT TRILINEAR COUPLINGS AT THE GUT SCALE
BLOCK AUGUT MGUT=  2.10086949E+16 # (GUT SCALE)
  3  3     5.97876433E+04   # ATOP
BLOCK ADGUT MGUT=  2.10086949E+16 # (GUT SCALE)
  3  3     3.26140851E+04   # ABOT
BLOCK AEGUT MGUT=  2.10086949E+16 # (GUT SCALE)
  3  3     2.63837779E+04   # ATAU
BLOCK ALGUT/AKGUT MGUT=  2.10086949E+16 # (GUT SCALE)
         1     1.13506016E+05   # ALAMBDA
         2     1.43426913E+05   # AKAPPA
# 
# SOFT MASSES SQUARED AT THE GUT SCALE
BLOCK MSOFTGUT MGUT=  2.10086949E+16 # (GUT SCALE)
         1     5.43544490E+02   # M1
         2     6.46206574E+02   # M2
         3     5.20293372E+02   # M3
        21     2.82594533E+09   # M_HD^2
        22     3.64457715E+09   # M_HU^2
        23     5.74138558E+09   # M_S^2
        31     1.64121131E+06   # M_eL^2
        32     1.64121131E+06   # M_muL^2
        33     1.66388725E+06   # M_tauL^2
        34     4.22386848E+05   # M_eR^2
        35     4.22386848E+05   # M_muR^2
        36     4.68676867E+05   # M_tauR^2
        41     5.81471927E+05   # M_q1L^2
        42     5.81471927E+05   # M_q2L^2
        43     2.53008569E+08   # M_q3L^2
        44     9.63907122E+05   # M_uR^2
        45     9.63907122E+05   # M_cR^2
        46     5.25968483E+08   # M_tR^2
        47     2.02153850E+05   # M_dR^2
        48     2.02153850E+05   # M_sR^2
        49     2.48848188E+05   # M_bR^2
