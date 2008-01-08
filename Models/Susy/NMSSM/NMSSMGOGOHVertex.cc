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

NMSSMGOGOHVertex::NMSSMGOGOHVertex() : _lambda(0.), _kappa(0.), _sinb(0.),
				       _cosb(0.), _sw(0.), _cw(0.), 
				       _q2last(0.*MeV2), _couplast(0.) {
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
      second.push_back(ineut[ix]);
      third .push_back(37);

      first .push_back(ichar[iy]);
      second.push_back(ineut[ix]);
      third .push_back(-37);

    }
  }
  setList(first, second, third);
}

void NMSSMGOGOHVertex::persistentOutput(PersistentOStream & os) const {
   os << _mixV << _mixU << _mixN << _mixS << _mixP << _lambda << _kappa << _sinb
      << _cosb << _sw << _cw;
}

void NMSSMGOGOHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _mixV >> _mixU >> _mixN >> _mixS >> _mixP >> _lambda >> _kappa >> _sinb
     >> _cosb >> _sw >> _cw;
}

void NMSSMGOGOHVertex::doinit() throw(InitException) {
  // SM parameters
  generator()->standardModel() = generator()->standardModel();
  _sw = generator()->standardModel()->sin2ThetaW();
  _cw=sqrt(1.-_sw);
  _sw=sqrt(_sw);
  // NMSSM parameters
  tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must have the NMSSM Model in "
			  << "NMSSMGOGOHVertex::doinit()"
			  << Exception::runerror;
  // get the mixing matrices
  // higgs
  _mixS=model->CPevenHiggsMix();
  if(!_mixS) 
    throw InitException() << "Mixing matrix for CP-even neutral Higgs"
			  << " bosons is not set in NMSSMGOGOHVertex::doinit()" 
			  << Exception::runerror;
  _mixP=model->CPoddHiggsMix();
  if(!_mixP) 
    throw InitException() << "Mixing matrix for CP-odd neutral Higgs"
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
				   tcPDPtr part3,int ) {
  long ihigg = part3->id();
  long ig1 = part1->id();
  long ig2 = part2->id();
  // set the prefactor to 1
  setNorm(1.);
  // electromagentic coupling
  if(q2!=_q2last) {
    _couplast = electroMagneticCoupling(q2);
    _q2last = q2;
  }
  double rt = sqrt(0.5);
  // CP-even neutral higgs
  if(ihigg == 25 || ihigg == 35 || ihigg == 45) {
    int iloc = (ihigg - 25)/10;
    // chargino
    if(abs(ig1) == 1000024 || abs(ig1) == 1000037) {
      int ic1 = (abs(ig1) == 1000024) ? 0 : 1;
      int ic2 = (abs(ig2) == 1000024) ? 0 : 1;
      Complex Q = _lambda*rt*(*_mixS)(iloc,2)*(*_mixU)(ic1,1)*(*_mixV)(ic2,1)
	-_couplast/_sw*rt*((*_mixS)(iloc,0)*(*_mixU)(ic1,1)*(*_mixV)(ic2,0)+
			   (*_mixS)(iloc,1)*(*_mixU)(ic1,0)*(*_mixV)(ic2,1));
      setLeft(conj(Q));
      setRight(Q);
    }
    // neutralino
    else {
      int in1 = (ig1 < 1000024) ? (ig1 - 1000022) : (ig1 - 1000005)/10; 
      int in2 = (ig2 < 1000024) ? (ig2 - 1000022) : (ig2 - 1000005)/10;
      Complex us1 = (*_mixS)(iloc, 0); Complex us2 = (*_mixS)(iloc, 1);
      Complex us3 = (*_mixS)(iloc, 2); Complex ni2 = (*_mixN)(in1,1);
      Complex nj2 = (*_mixN)(in2,1); Complex ni3 = (*_mixN)(in1,2);
      Complex nj3 = (*_mixN)(in2,2); Complex ni4 = (*_mixN)(in1,3);
      Complex nj4 = (*_mixN)(in2,3); Complex ni5 = (*_mixN)(in1,4);
      Complex nj5 = (*_mixN)(in2,4);
      Complex Qdp = 
	0.5*( (us1*_cosb + us2*_sinb) * 
	      ( _couplast/_sw/_cw*(ni2*conj(nj3) + nj2*conj(ni3)) 
		+ sqrt(2)*_lambda*(ni5*conj(nj4) + nj5*conj(ni4)) ) 
	      + (us1*_sinb - us2*_cosb)*
	      ( _couplast/_sw/_cw*(ni2*conj(nj4) + nj2*conj(ni4)) 
		- sqrt(2)*_lambda*(ni5*conj(nj3) + nj5*conj(ni3)) ) )
	- sqrt(2)*_kappa*us3*(ni5*conj(nj5) + nj5*conj(ni5));
      setLeft(-Qdp);
      setRight(-conj(Qdp));
    }
  }
  // CP-odd  neutral higgs
  else if(ihigg==36||ihigg==46) {
    int iloc = (ihigg-36)/10;
    
    // chargino
    if(abs(ig1)==1000024||abs(ig1)==1000037) {
      int ic1 = abs(ig1)==1000024 ? 0 : 1;
      int ic2 = abs(ig2)==1000024 ? 0 : 1;
      Complex R = -_lambda*rt*(*_mixP)(iloc,2)*(*_mixU)(ic1,1)*(*_mixV)(ic2,1)
	-_couplast/_sw*rt*((*_mixP)(iloc,0)*(*_mixU)(ic1,1)*(*_mixV)(ic2,0)+
			   (*_mixP)(iloc,1)*(*_mixU)(ic1,0)*(*_mixV)(ic2,1));
      R *= Complex(0., 1.);
      setLeft(conj(R));
      setRight(R);
    }
    // neutralino
    else {
      int in1 = (ig1 < 1000024) ? (ig1 - 1000022) : (ig1 - 1000005)/10; 
      int in2 = (ig2 < 1000024) ? (ig2 - 1000022) : (ig2 - 1000005)/10; 
      Complex up1 = (*_mixP)(iloc, 0); Complex up2 = (*_mixP)(iloc, 1);
      Complex up3 = (*_mixP)(iloc, 2); Complex ni2 = (*_mixN)(in1,1);
      Complex nj2 = (*_mixN)(in2,1); Complex ni3 = (*_mixN)(in1,2);
      Complex nj3 = (*_mixN)(in2,2); Complex ni4 = (*_mixN)(in1,3);
      Complex nj4 = (*_mixN)(in2,3); Complex ni5 = (*_mixN)(in1,4);
      Complex nj5 = (*_mixN)(in2,4);

      Complex Rdp = 
	-0.5*( (up1*_cosb + up2*_sinb) * 
	      ( _couplast/_sw/_cw*(ni2*conj(nj3) + nj2*conj(ni3)) 
		- sqrt(2)*_lambda*(ni5*conj(nj4) + nj5*conj(ni4)) ) 
	      + (up1*_sinb - up2*_cosb)*
	      ( _couplast/_sw/_cw*(ni2*conj(nj4) + nj2*conj(ni4)) 
		+ sqrt(2)*_lambda*(ni5*conj(nj3) + nj5*conj(ni3)) ) )
	- sqrt(2)*_kappa*up3*(ni5*conj(nj5) + nj5*conj(ni5));
      Rdp *= Complex(0., 1.);
      setLeft(-Rdp);
      setRight(-conj(Rdp));
    }
  }
  // charged higgs
  else {
    int in,ic;
    if(abs(ig1)==1000024||abs(ig1)==1000037) {
      in = (ig1 < 1000024) ? (ig1-1000022) : (ig1-1000005)/10; 
      ic = (abs(ig2) == 1000024) ? 0 : 1;
    }
    else {
      in = (ig2 < 1000024) ? (ig2-1000022) : (ig2-1000005)/10; 
      ic = (abs(ig1) == 1000024) ? 0 : 1;
    }
    Complex QpL = _couplast/_sw*_cosb * 
      ( (*_mixV)(ic,0)*((*_mixN)(in,3)*_cosb - (*_mixN)(in,2)*_sinb)
	+ rt*(*_mixV)(ic,1)*( 2.*_sw*(*_mixN)(in,0) 
			      + (sqr(_cw) - sqr(_sw))*(*_mixN)(in,1)/_cw) )
      - _lambda*_sinb*(*_mixV)(ic,1)*(*_mixN)(in,4);
    
    Complex QpR = _couplast/_sw*_sinb * 
      ( (*_mixU)(ic,0)*((*_mixN)(in,3)*_sinb + (*_mixN)(in,2)*_cosb)
	- rt*(*_mixU)(ic,1)*( 2.*_sw*(*_mixN)(in,0) 
			      + (sqr(_cw) - sqr(_sw))*(*_mixN)(in,1)/_cw) )
      - _lambda*_cosb*(*_mixU)(ic,1)*(*_mixN)(in,4);
    
    if(ihigg > 0) {
      setLeft (-conj(QpL));
      setRight(-QpR);
    }
    else {
      setLeft (-QpL);
      setRight(-conj(QpR));
    }
  }
}
