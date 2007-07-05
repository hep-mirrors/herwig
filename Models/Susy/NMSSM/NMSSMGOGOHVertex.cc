// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMGOGOHVertex class.
//

#include "NMSSMGOGOHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMGOGOHVertex::NMSSMGOGOHVertex() {
  _q2last=0.*GeV2;
  _couplast=0.;
  int ieven[3]={25,35,45};
  int iodd [2]={36,46};
  int ichar[2]={1000024,1000037};
  int ineut[5]={1000022,1000023,1000025,1000035,1000045};
  vector<int> first,second,third;
  // CP-even charginos
  for(unsigned int ix=0;ix<2;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      for(unsigned int iz=0;iz<3;++iz) {
	first .push_back(-ichar[ix]);
	second.push_back( ichar[iy]);
	third .push_back( ieven[iz]);
      }
    }
  }
  // CP-odd charginos
  for(unsigned int ix=0;ix<2;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      for(unsigned int iz=0;iz<2;++iz) {
	first .push_back(-ichar[ix]);
	second.push_back( ichar[iy]);
	third .push_back( iodd [iz]);
      }
    }
  }
  // CP-even neutralinos
  for(unsigned int ix=0;ix<5;++ix) {
    for(unsigned int iy=0;iy<5;++iy) {
      for(unsigned int iz=0;iz<3;++iz) {
	first .push_back( ineut[ix]);
	second.push_back( ineut[iy]);
	third .push_back( ieven[iz]);
      }
    }
  }
  // CP-odd  neutralinos
  for(unsigned int ix=0;ix<5;++ix) {
    for(unsigned int iy=0;iy<5;++iy) {
      for(unsigned int iz=0;iz<2;++iz) {
	first .push_back( ineut[ix]);
	second.push_back( ineut[iy]);
	third .push_back( iodd[iz]);
      }
    }
  }
  // charged higgs
  for(unsigned int ix=0;ix<5;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      first .push_back(-ichar[iy]);
      second.push_back( ineut[ix]);
      third .push_back(    37    );
    }
  }
}

void NMSSMGOGOHVertex::persistentOutput(PersistentOStream & os) const {
  os << _mixV << _mixU << _mixN << _mixS << _mixP << _lambda << _kappa << _sinb
     << _cosb << _sw << _cw << _theSM;
}

void NMSSMGOGOHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _mixV >> _mixU >> _mixN >> _mixS >> _mixP >> _lambda >> _kappa >> _sinb
     >> _cosb >> _sw >> _cw >> _theSM;
}

void NMSSMGOGOHVertex::doinit() throw(InitException) {
  // SM parameters
  _theSM = generator()->standardModel();
  _sw = _theSM->sin2ThetaW();
  _cw=sqrt(1.-_sw);
  _sw=sqrt(_sw);
  // NMSSM parameters
  tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(_theSM);
  if(!model) 
    throw InitException() << "Must have the NMSSM Model in NMSSMGOGOHVertex::doinit()"
			  << Exception::runerror;
  // get the mixing matrices
  // higgs
  _mixS=model->CPevenHiggsMix();
  if(!_mixS) throw InitException() << "Mixing matrix for CP-even neutral Higgs"
				   << " bosons is not set in NMSSMGOGOHVertex::doinit()" 
				   << Exception::runerror;
  _mixP=model->CPoddHiggsMix();
  if(!_mixP) throw InitException() << "Mixing matrix for CP-odd neutral Higgs"
				   << " bosons is not set in NMSSMGOGOHVertex::doinit()" 
				   << Exception::runerror;
  // charginos
  _mixU = model->charginoUMix();
  _mixV = model->charginoVMix();
  if(!_mixU || !_mixV)
    throw InitException() << "NMSSMGOGOHVertex::doinit - "
			  << "A mixing matrix pointer is null.  U: " 
			  << _mixU << "  V: " << _mixV
			  << Exception::abortnow;
  // neutralinos
  _mixN  = model->neutralinoMix();
  if(!_mixN)
    throw InitException() << "NMSSMGOGOHVertex::doinit - The neutralino "
			  << "mixing matrix pointer is null." 
			  << Exception::abortnow;
  // kappa and lambda couplings
  _lambda = model->lambda();
  _kappa  = model->kappa();
  // sin and cos beta
  double beta = atan(model->tanBeta());
  _sinb=sin(beta);
  _cosb=cos(beta);
  // order in the couplings
  orderInGem(1);
  orderInGs(0);
  FFSVertex::doinit();
}

ClassDescription<NMSSMGOGOHVertex> NMSSMGOGOHVertex::initNMSSMGOGOHVertex;
// Definition of the static class description member.

void NMSSMGOGOHVertex::Init() {

  static ClassDocumentation<NMSSMGOGOHVertex> documentation
    ("The NMSSMGOGOHVertex class implements the couplings of the Higgs bosons"
     " of the NMSSM and the electroweak gauginos");

}

void NMSSMGOGOHVertex::setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
				   tcPDPtr part3,int ioff) {
  int ihigg = part3->id();
  int ig1 = part1->id();
  int ig2 = part2->id();
  // set the prefactor to 1
  setNorm(1.);
  // electromagentic coupling
  if(q2!=_q2last) {
    double alpha = _theSM->alphaEM(q2);
    _couplast = sqrt(4.0*Constants::pi*alpha);
  }
  double rt = sqrt(0.5);
  // CP-even neutral higgs
  if(ihigg==25||ihigg==35||ihigg==45) {
    int iloc = (ihigg-25)/10;
    // chargino
    if(abs(ig1)==1000024||abs(ig1)==1000037) {
      int ic1 = abs(ig1)==1000024 ? 0 : 1;
      int ic2 = abs(ig2)==1000024 ? 0 : 1;
      Complex Q = _lambda*rt*(*_mixS)(iloc,2)*(*_mixU)(ic1,1)*(*_mixV)(ic2,1)
	+_couplast/_sw*rt*((*_mixS)(iloc,0)*(*_mixU)(ic1,1)*(*_mixV)(ic2,0)+
			   (*_mixS)(iloc,1)*(*_mixU)(ic1,0)*(*_mixV)(ic2,1));
    }
    // neutralino
    else {
      int in1 = ig1<1000024 ? (ig1-1000022) : (ig1-1000005)/10; 
      int in2 = ig2<1000024 ? (ig2-1000022) : (ig2-1000005)/10; 
      Complex Qdp;
    }
  }
  // CP-odd  neutral higgs
  else if(ihigg==36||ihigg==46) {
    int iloc = (ihigg-36)/10;
    // chargino
    if(abs(ig1)==1000024||abs(ig1)==1000037) {
      int ic1 = abs(ig1)==1000024 ? 0 : 1;
      int ic2 = abs(ig2)==1000024 ? 0 : 1;
      Complex R = _lambda*rt*(*_mixP)(iloc,2)*(*_mixU)(ic1,1)*(*_mixV)(ic2,1)
	-_couplast/_sw*rt*((*_mixP)(iloc,0)*(*_mixU)(ic1,1)*(*_mixV)(ic2,0)+
			   (*_mixP)(iloc,1)*(*_mixU)(ic1,0)*(*_mixV)(ic2,1));
    }
    // neutralino
    else {
      int in1 = ig1<1000024 ? (ig1-1000022) : (ig1-1000005)/10; 
      int in2 = ig2<1000024 ? (ig2-1000022) : (ig2-1000005)/10; 
      Complex Rdp;
    }
  }
  // charged higgs
  else {
    int in,ic;
    if(abs(ig1)==1000024||abs(ig1)==1000037) {
      in = ig1<1000024 ? (ig1-1000022) : (ig1-1000005)/10; 
      ic = abs(ig2)==1000024 ? 0 : 1;
    }
    else {
      in = ig2<1000024 ? (ig2-1000022) : (ig2-1000005)/10; 
      ic = abs(ig1)==1000024 ? 0 : 1;
    }
    Complex QpL = _couplast/_sw*_cosb*((*_mixV)(ic,0)*(*_mixN)(in,3)
				       +rt*(*_mixV)(ic,1)*((*_mixN)(in,1)
							   +_sw/_cw*(*_mixN)(in,0)))
      +_lambda*_sinb*(*_mixV)(ic,1)*(*_mixN)(in,4);
    Complex QpR = _couplast/_sw*_sinb*((*_mixU)(ic,0)*(*_mixN)(in,2)
				       -rt*(*_mixU)(ic,1)*((*_mixN)(in,1)
							   +_sw/_cw*(*_mixN)(in,0)))
      +_lambda*_cosb*(*_mixU)(ic,1)*(*_mixN)(in,4);
    if(ihigg>0) {
      setLeft (-conj(QpL));
      setRight(-QpR);
    }
    else {
      setLeft (-QpL);
      setRight(-conj(QpR));
    }
  }
}
