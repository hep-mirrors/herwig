// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPV class.
//

#include "RPV.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void RPV::persistentOutput(PersistentOStream & os ) const {
  os << lambdaLLE_ << lambdaLQD_ << lambdaUDD_
     << LLEVertex_ << LQDVertex_ << UDDVertex_;
}

void RPV::persistentInput(PersistentIStream & is, int) {
  is >> lambdaLLE_ >> lambdaLQD_ >> lambdaUDD_
     >> LLEVertex_ >> LQDVertex_ >> UDDVertex_;
}

ClassDescription<RPV> RPV::initRPV;
// Definition of the static class description member.

void RPV::Init() {

  static ClassDocumentation<RPV> documentation
    ("The RPV class is the base class for the implementation of the"
     " R-parity violating MSSM.");

  static Reference<RPV,AbstractFFSVertex> interfaceLLEVertex
    ("Vertex/LLE",
     "The vertex for the trillinear LLE interaction",
     &RPV::LLEVertex_, false, false, true, false, false);

  static Reference<RPV,AbstractFFSVertex> interfaceLQDVertex
    ("Vertex/LQD",
     "The vertex for the trillinear LQD interaction",
     &RPV::LQDVertex_, false, false, true, false, false);

  static Reference<RPV,AbstractFFSVertex> interfaceUDDVertex
    ("Vertex/UDD",
     "The vertex for the trillinear UDD interaction",
     &RPV::UDDVertex_, false, false, true, false, false);

}

void RPV::extractParameters(bool checkmodel) {
  MSSM::extractParameters(false);
  if(checkmodel) {
    map<string,ParamMap>::const_iterator pit;
    pit = parameters().find("modsel");
    if(pit == parameters().end()) return;
    ParamMap::const_iterator it;
    // nmssm or mssm
    it = pit->second.find(3);
    int inmssm = (it != pit->second.end()) ? int(it->second) : 0;
    if(inmssm != 0) 
      throw Exception() << "R-parity violating MSSM model"
			<< " used but NMSSM read in." << Exception::runerror; 
    // RPV
    it = pit->second.find(4);
    int irpv = (it != pit->second.end()) ? int(it->second) : 0;
    if(irpv != 1) throw Exception() << "RPV model used but no RPV in input file"
				    << Exception::runerror; 
    // CPV
    it = pit->second.find(5);
    int icpv = (it != pit->second.end()) ? int(it->second) : 0;
    if(icpv != 0) throw Exception() << "RPV model does not support CPV" 
				  << Exception::runerror; 
    // flavour violation
    it = pit->second.find(6);
    int ifv = (it != pit->second.end()) ? int(it->second) : 0;
    if(ifv != 0) throw Exception() << "RPV model does not support "
				 << "flavour violation"
				 << Exception::runerror;
  }
  // get the RPV parameters
  map<string,ParamMap>::const_iterator pit;
  pit=parameters().find("rvlamlle");
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first==-1) continue;
      int i = it->first/100-1;
      int k = it->first%10-1;
      int j = (it->first%100)/10-1;
      lambdaLLE_[i][j][k] = it->second;
    }
  }
  pit=parameters().find("rvlamlqd");
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first==-1) continue;
      int i = it->first/100-1;
      int k = it->first%10-1;
      int j = (it->first%100)/10-1;
      lambdaLQD_[i][j][k] = it->second;
    }
  }
  pit=parameters().find("rvlamudd");
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first==-1) continue;
      int i = it->first/100-1;
      int k = it->first%10-1;
      int j = (it->first%100)/10-1;
      lambdaUDD_[i][j][k] = it->second;
    }
  }
}

void RPV::createMixingMatrices() {
//   map<string,pair<MatrixSize, MixingVector> >::const_iterator it;
//   for(it=mixings().begin();it!=mixings().end();++it) {
//     string name=it->first;
//     // pseudo-scalar higgs mixing
//     if (name == "nmamix") {
//       createMixingMatrix(theHiggsAMix,name,it->second.second,it->second.first);
//     }
//   }
  // base class for neutralinos and charginos
  MSSM::createMixingMatrices();
}
// Block MINPAR  # SUSY breaking input parameters
//      3    1.000000000000000e+01   # tanb
//      4    1.000000000000000e+00   # sign(mu)
//      1    1.000000000000000e+02   # m0
//      2    2.500000000000000e+02   # m12
//      5   -1.000000000000000e+02   # A0


// # Higgs mixing
// Block alpha   # Effective Higgs mixing parameter
//           -1.149203839391468e-01   # alpha

// Block stopmix  # stop mixing matrix
//   1  1     5.567136252574828e-01   # O_{11}
//   1  2     8.307044838284376e-01   # O_{12}
//   2  1     8.307044838284376e-01   # O_{21}
//   2  2    -5.567136252574828e-01   # O_{22}

// Block sbotmix  # sbottom mixing matrix
//   1  1     9.494595966791360e-01   # O_{11}
//   1  2     3.138892707212089e-01   # O_{12}
//   2  1    -3.138892707212089e-01   # O_{21}
//   2  2     9.494595966791360e-01   # O_{22}

// Block staumix  # stau mixing matrix
//   1  1     2.666587767723425e-01   # O_{11}
//   1  2     9.637910026402394e-01   # O_{12}
//   2  1     9.637910026402394e-01   # O_{21}
//   2  2    -2.666587767723425e-01   # O_{22}

// Block nmix  # neutralino mixing matrix
//   1  1     9.851357669950006e-01   # N_{1,1}
//   1  2    -5.771661350693690e-02   # N_{1,2}
//   1  3     1.517355210244914e-01   # N_{1,3}
//   1  4    -5.614841735871662e-02   # N_{1,4}
//   2  1     1.073787916649285e-01   # N_{2,1}
//   2  2     9.402548967081881e-01   # N_{2,2}
//   2  3    -2.795671827746073e-01   # N_{2,3}
//   2  4     1.619651648729555e-01   # N_{2,4}
//   3  1    -6.112790189200206e-02   # N_{3,1}
//   3  2     9.103171953480492e-02   # N_{3,2}
//   3  3     6.945942965173417e-01   # N_{3,3}
//   3  4     7.109960399990973e-01   # N_{3,4}
//   4  1    -1.193343843912277e-01   # N_{4,1}
//   4  2     3.229593593319475e-01   # N_{4,2}
//   4  3     6.452575340284490e-01   # N_{4,3}
//   4  4    -6.819818705078531e-01   # N_{4,4}

// Block Umix  # chargino U mixing matrix 
//   1  1     9.140759450766424e-01   # U_{1,1}
//   1  2    -4.055430515151790e-01   # U_{1,2}
//   2  1     4.055430515151790e-01   # U_{2,1}
//   2  2     9.140759450766424e-01   # U_{2,2}

// Block Vmix  # chargino V mixing matrix 
//   1  1     9.711123714587470e-01   # V_{1,1}
//   1  2    -2.386226351370898e-01   # V_{1,2}
//   2  1     2.386226351370898e-01   # V_{2,1}
//   2  2     9.711123714587470e-01   # V_{2,2}

// Block RVT Q= 9.118760000000000e+01 # R-Parity violating LLE soft terms 
//   1 1 1    0.000000000000000e+00   # T_{111}
//   1 1 2    0.000000000000000e+00   # T_{112}
//   1 1 3    0.000000000000000e+00   # T_{113}
//   1 2 1    0.000000000000000e+00   # T_{121}
//   1 2 2    0.000000000000000e+00   # T_{122}
//   1 2 3    4.243543810455972e+01   # T_{123}
//   1 3 1    0.000000000000000e+00   # T_{131}
//   1 3 2    0.000000000000000e+00   # T_{132}
//   1 3 3    0.000000000000000e+00   # T_{133}
//   2 1 1    0.000000000000000e+00   # T_{211}
//   2 1 2    0.000000000000000e+00   # T_{212}
//   2 1 3   -4.243543810455972e+01   # T_{213}
//   2 2 1    0.000000000000000e+00   # T_{221}
//   2 2 2    0.000000000000000e+00   # T_{222}
//   2 2 3    0.000000000000000e+00   # T_{223}
//   2 3 1    0.000000000000000e+00   # T_{231}
//   2 3 2    0.000000000000000e+00   # T_{232}
//   2 3 3    0.000000000000000e+00   # T_{233}
//   3 1 1    0.000000000000000e+00   # T_{311}
//   3 1 2    0.000000000000000e+00   # T_{312}
//   3 1 3    0.000000000000000e+00   # T_{313}
//   3 2 1    0.000000000000000e+00   # T_{321}
//   3 2 2    0.000000000000000e+00   # T_{322}
//   3 2 3    0.000000000000000e+00   # T_{323}
//   3 3 1    0.000000000000000e+00   # T_{331}
//   3 3 2    0.000000000000000e+00   # T_{332}
//   3 3 3    0.000000000000000e+00   # T_{333}
// Block RVTP Q= 9.118760000000000e+01 # R-Parity violating LQD soft terms 
//   1 1 1    0.000000000000000e+00   # T'_{111}
//   1 1 2    0.000000000000000e+00   # T'_{112}
//   1 1 3    0.000000000000000e+00   # T'_{113}
//   1 2 1    0.000000000000000e+00   # T'_{121}
//   1 2 2    0.000000000000000e+00   # T'_{122}
//   1 2 3    0.000000000000000e+00   # T'_{123}
//   1 3 1    0.000000000000000e+00   # T'_{131}
//   1 3 2    0.000000000000000e+00   # T'_{132}
//   1 3 3    0.000000000000000e+00   # T'_{133}
//   2 1 1    0.000000000000000e+00   # T'_{211}
//   2 1 2    0.000000000000000e+00   # T'_{212}
//   2 1 3    0.000000000000000e+00   # T'_{213}
//   2 2 1    0.000000000000000e+00   # T'_{221}
//   2 2 2    0.000000000000000e+00   # T'_{222}
//   2 2 3    0.000000000000000e+00   # T'_{223}
//   2 3 1    0.000000000000000e+00   # T'_{231}
//   2 3 2    0.000000000000000e+00   # T'_{232}
//   2 3 3    0.000000000000000e+00   # T'_{233}
//   3 1 1    0.000000000000000e+00   # T'_{311}
//   3 1 2    0.000000000000000e+00   # T'_{312}
//   3 1 3    0.000000000000000e+00   # T'_{313}
//   3 2 1    0.000000000000000e+00   # T'_{321}
//   3 2 2    0.000000000000000e+00   # T'_{322}
//   3 2 3    0.000000000000000e+00   # T'_{323}
//   3 3 1    0.000000000000000e+00   # T'_{331}
//   3 3 2    0.000000000000000e+00   # T'_{332}
//   3 3 3    0.000000000000000e+00   # T'_{333}
// Block RVTPP Q= 9.118760000000000e+01 # R-Parity violating UDD soft terms 
//   1 1 1    0.000000000000000e+00   # T''_{111}
//   1 1 2    0.000000000000000e+00   # T''_{112}
//   1 1 3    0.000000000000000e+00   # T''_{113}
//   1 2 1    0.000000000000000e+00   # T''_{121}
//   1 2 2    0.000000000000000e+00   # T''_{122}
//   1 2 3    0.000000000000000e+00   # T''_{123}
//   1 3 1    0.000000000000000e+00   # T''_{131}
//   1 3 2    0.000000000000000e+00   # T''_{132}
//   1 3 3    0.000000000000000e+00   # T''_{133}
//   2 1 1    0.000000000000000e+00   # T''_{211}
//   2 1 2    0.000000000000000e+00   # T''_{212}
//   2 1 3    0.000000000000000e+00   # T''_{213}
//   2 2 1    0.000000000000000e+00   # T''_{221}
//   2 2 2    0.000000000000000e+00   # T''_{222}
//   2 2 3    0.000000000000000e+00   # T''_{223}
//   2 3 1    0.000000000000000e+00   # T''_{231}
//   2 3 2    0.000000000000000e+00   # T''_{232}
//   2 3 3    0.000000000000000e+00   # T''_{233}
//   3 1 1    0.000000000000000e+00   # T''_{311}
//   3 1 2    0.000000000000000e+00   # T''_{312}
//   3 1 3    0.000000000000000e+00   # T''_{313}
//   3 2 1    0.000000000000000e+00   # T''_{321}
//   3 2 2    0.000000000000000e+00   # T''_{322}
//   3 2 3    0.000000000000000e+00   # T''_{323}
//   3 3 1    0.000000000000000e+00   # T''_{331}
//   3 3 2    0.000000000000000e+00   # T''_{332}
//   3 3 3    0.000000000000000e+00   # T''_{333}
// Block RVKAPPA Q= 9.118760000000000e+01 # R-Parity violating kappa 
//      1    0.000000000000000e+00   # kappa_{1}
//      2    0.000000000000000e+00   # kappa_{2}
//      3    0.000000000000000e+00   # kappa_{3}
// Block RVD Q= 9.118760000000000e+01 # R-Parity violating D 
//      1    0.000000000000000e+00   # D_{1}
//      2    0.000000000000000e+00   # D_{2}
//      3    0.000000000000000e+00   # D_{3}
// Block RVSNVEV Q= 9.118760000000000e+01 # sneutrino VEVs D 
//      1    0.000000000000000e+00   # SneutrinoVev_{1}
//      2    0.000000000000000e+00   # SneutrinoVev_{2}
//      3    0.000000000000000e+00   # SneutrinoVev_{3}
// Block RVM2LH1 Q= 9.118760000000000e+01 # M2LH1 
//      1    0.000000000000000e+00   # M2LH1_{1}
//      2    0.000000000000000e+00   # M2LH1_{2}
//      3    0.000000000000000e+00   # M2LH1_{3}
// Block RVNMIX Q= 9.118760000000000e+01 # neutrino-neutralino mixing matrix 
//   1 1    1.000000000000000e+00   # N_{11}
//   1 2    0.000000000000000e+00   # N_{12}
//   1 3    0.000000000000000e+00   # N_{13}
//   1 4    0.000000000000000e+00   # N_{14}
//   1 5    0.000000000000000e+00   # N_{15}
//   1 6    0.000000000000000e+00   # N_{16}
//   1 7    0.000000000000000e+00   # N_{17}
//   2 1    0.000000000000000e+00   # N_{21}
//   2 2    1.000000000000000e+00   # N_{22}
//   2 3    0.000000000000000e+00   # N_{23}
//   2 4    0.000000000000000e+00   # N_{24}
//   2 5    0.000000000000000e+00   # N_{25}
//   2 6    0.000000000000000e+00   # N_{26}
//   2 7    0.000000000000000e+00   # N_{27}
//   3 1    0.000000000000000e+00   # N_{31}
//   3 2    0.000000000000000e+00   # N_{32}
//   3 3    1.000000000000000e+00   # N_{33}
//   3 4    0.000000000000000e+00   # N_{34}
//   3 5    0.000000000000000e+00   # N_{35}
//   3 6    0.000000000000000e+00   # N_{36}
//   3 7    0.000000000000000e+00   # N_{37}
//   4 1    0.000000000000000e+00   # N_{41}
//   4 2    0.000000000000000e+00   # N_{42}
//   4 3    0.000000000000000e+00   # N_{43}
//   4 4    9.847565328754548e-01   # N_{44}
//   4 5    1.103472257924572e-01   # N_{45}
//   4 6   -6.109421671307692e-02   # N_{46}
//   4 7   -1.197729410310909e-01   # N_{47}
//   5 1    0.000000000000000e+00   # N_{51}
//   5 2    0.000000000000000e+00   # N_{52}
//   5 3    0.000000000000000e+00   # N_{53}
//   5 4   -6.117909738186333e-02   # N_{54}
//   5 5    9.417013337872188e-01   # N_{55}
//   5 6    9.139986887133542e-02   # N_{56}
//   5 7    3.179650609064091e-01   # N_{57}
//   6 1    0.000000000000000e+00   # N_{61}
//   6 2    0.000000000000000e+00   # N_{62}
//   6 3    0.000000000000000e+00   # N_{63}
//   6 4    1.526096816059464e-01   # N_{64}
//   6 5   -2.757226329309796e-01   # N_{65}
//   6 6    6.951143481105823e-01   # N_{66}
//   6 7    6.461449975203246e-01   # N_{67}
//   7 1    0.000000000000000e+00   # N_{71}
//   7 2    0.000000000000000e+00   # N_{72}
//   7 3    0.000000000000000e+00   # N_{73}
//   7 4   -5.676243549025390e-02   # N_{74}
//   7 5    1.581110919350379e-01   # N_{75}
//   7 6    7.104432445349301e-01   # N_{76}
//   7 7   -6.834100561295584e-01   # N_{77}
// Block hmix Q= 4.658779932434547e+02  # Higgs mixing parameters
//      1     3.523836232180990e+02   # mu(Q)MSSM DRbar
//      2     9.750797435881633e+00   # tan beta(Q)MSSM DRbar
//      3     2.450549911158079e+02   # higgs vev(Q)MSSM DRbar
//      4     1.540022571859735e+05   # mA^2(Q)MSSM DRbar
// Block msoft Q= 4.658779932434547e+02 # MSSM DRbar SUSY breaking parameters
//      1     1.019307676894905e+02   # M_1(Q)
//      2     1.918845353043126e+02   # M_2(Q)
//      3     5.864985652616149e+02   # M_3(Q)
//     21     3.224468700882639e+04   # mH1^2(Q)
//     22    -1.265998776222425e+05   # mH2^2(Q)
//     31     1.929790251783819e+02   # meL(Q)
//     32     1.929762221011005e+02   # mmuL(Q)
//     33     1.943010234240223e+02   # mtauL(Q)
//     34     1.359061860029394e+02   # meR(Q)
//     35     1.358981170295141e+02   # mmuR(Q)
//     36     1.270729528972383e+02   # mtauR(Q)
//     41     5.456782453568811e+02   # mqL1(Q)
//     42     5.456765542407577e+02   # mqL2(Q)
//     43     4.973670810960051e+02   # mqL3(Q)
//     44     5.277735920289771e+02   # muR(Q)
//     45     5.277718568028397e+02   # mcR(Q)
//     46     4.228101648527092e+02   # mtR(Q)
//     47     5.256849074418041e+02   # mdR(Q)
//     48     5.256830975517062e+02   # msR(Q)
//     49     5.224091548560939e+02   # mbR(Q)
// Block au Q= 4.658779932434547e+02  
//   1  1    -6.800325328454909e+02   # Au(Q)MSSM DRbar
//   2  2    -6.800290625507030e+02   # Ac(Q)MSSM DRbar
//   3  3    -5.004599927360533e+02   # At(Q)MSSM DRbar
// Block ad Q= 4.658779932434547e+02  
//   1  1    -8.573073511584474e+02   # Ad(Q)MSSM DRbar
//   2  2    -8.573041034412643e+02   # As(Q)MSSM DRbar
//   3  3    -7.942004406414295e+02   # Ab(Q)MSSM DRbar
// Block ae Q= 4.658779932434547e+02  
//   1  1    -2.510763383610091e+02   # Ae(Q)MSSM DRbar
//   2  2    -2.510705792217225e+02   # Amu(Q)MSSM DRbar
//   3  3    -2.478683315008348e+02   # Atau(Q)MSSM DRbar

void RPV::doinit() {
  MSSM::doinit();
  addVertex(LLEVertex_);
  addVertex(LQDVertex_);
  addVertex(UDDVertex_);
}
