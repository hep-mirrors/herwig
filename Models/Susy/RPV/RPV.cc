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
  os << lambdaLLE_ << lambdaLQD_ << lambdaUDD_ << ounit(vnu_,GeV)
     << LLEVertex_ << LQDVertex_ << UDDVertex_
     << HiggsAMix_ << HiggsPMix_;
}

void RPV::persistentInput(PersistentIStream & is, int) {
  is >> lambdaLLE_ >> lambdaLQD_ >> lambdaUDD_ >> iunit(vnu_,GeV)
     >> LLEVertex_ >> LQDVertex_ >> UDDVertex_
     >> HiggsAMix_ >> HiggsPMix_;
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
  // lambda
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
  // lambda'
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
  // lambda''
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
  // sneutrino vevs
  pit=parameters().find("rvsnvev");
  vnu_.resize(3);
  if( pit != parameters().end() ) {
    for(ParamMap::const_iterator it = pit->second.begin();
	it!=pit->second.end();++it) {
      if(it->first>0) {
	assert(it->first>=1&&it->first<=3);
	vnu_[it->first-1] = it->second*GeV;
      }
    }
  }
}

void RPV::createMixingMatrices() {
  map<string,pair<MatrixSize, MixingVector> >::const_iterator it;
  for(it=mixings().begin();it!=mixings().end();++it) {
    string name=it->first;
    cerr << "testing in mixings loop " << name << "\n";
    // pseudo-scalar higgs mixing
    if (name == "rvamix") {
      cerr << "testing in mixing create A\n";
      createMixingMatrix(HiggsAMix_,name,it->second.second,it->second.first);
    }
    else if (name == "rvlmix") {
      cerr << "testing in mixing create C\n";
      createMixingMatrix(HiggsPMix_,name,it->second.second,it->second.first);
    }
  }
  // base class for neutralinos and charginos
  MSSM::createMixingMatrices();
}













// # Higgs mixing
// Block alpha   # Effective Higgs mixing parameter
//           -1.149203839391468e-01   # alpha

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
