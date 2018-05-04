// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMGOGOHVertex class.
//

#include "NMSSMGOGOHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMGOGOHVertex::NMSSMGOGOHVertex() : _lambda(0.), _kappa(0.), _sinb(0.),
				       _cosb(0.), _sw(0.), _cw(0.),
				       _q2last(0.*MeV2), _couplast(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void NMSSMGOGOHVertex::persistentOutput(PersistentOStream & os) const {
   os << _mixV << _mixU << _mixN << _mixS << _mixP << _lambda << _kappa << _sinb
      << _cosb << _sw << _cw;
}

void NMSSMGOGOHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _mixV >> _mixU >> _mixN >> _mixS >> _mixP >> _lambda >> _kappa >> _sinb
     >> _cosb >> _sw >> _cw;
}

void NMSSMGOGOHVertex::doinit() {
  int ieven[3]={25,35,45};
  int iodd [2]={36,46};
  long ichar[2]={1000024,1000037};
  long ineut[5]={1000022,1000023,1000025,1000035,1000045};
  // CP-even charginos
  for(unsigned int ix=0;ix<2;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      for(unsigned int iz=0;iz<3;++iz) {
	addToList(-ichar[ix], ichar[iy], ieven[iz]);
      }
    }
  }
  // CP-odd charginos
  for(unsigned int ix=0;ix<2;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      for(unsigned int iz=0;iz<2;++iz) {
	addToList(-ichar[ix], ichar[iy], iodd [iz]);

      }
    }
  }
  // CP-even neutralinos
  for(unsigned int ix=0;ix<5;++ix) {
    for(unsigned int iy=0;iy<5;++iy) {
      for(unsigned int iz=0;iz<3;++iz) {
	addToList( ineut[ix], ineut[iy], ieven[iz]);
        }
      }
    }
  // CP-odd  neutralinos
  for(unsigned int ix=0;ix<5;++ix) {
    for(unsigned int iy=0;iy<5;++iy) {
      for(unsigned int iz=0;iz<2;++iz) {
	addToList( ineut[ix], ineut[iy], iodd[iz]);
        }
      }
    }

  // charged higgs
  for(unsigned int ix=0;ix<5;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      addToList(ineut[ix], -ichar[iy], 37);

      addToList(ineut[ix], ichar[iy], -37);

    }
  }

   tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(generator()->standardModel());
     // SM parameters
  // sin theta_W

  double sw2=sin2ThetaW();
  _cw=sqrt(1.0 - sw2);
  _sw=sqrt(sw2);
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
  FFSVertex::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSMGOGOHVertex,FFSVertex>
describeHerwigNMSSMGOGOHVertex("Herwig::NMSSMGOGOHVertex", "HwSusy.so HwNMSSM.so");

void NMSSMGOGOHVertex::Init() {

  static ClassDocumentation<NMSSMGOGOHVertex> documentation
    ("The NMSSMGOGOHVertex class implements the couplings of the Higgs bosons"
     " of the NMSSM and the electroweak gauginos");

}

void NMSSMGOGOHVertex::setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
				   tcPDPtr part3) {
  long id1(part1->id()), id2(part2->id()), 
    id3(part3->id()), ihigg(0), ig1(0), ig2(0);
  if( abs(id1) == 25 || abs(id1) == 35 || abs(id1) == 45 || 
      abs(id1) == 36 || abs(id1) == 46 || abs(id1) == 37 ) {
    ihigg = id1;
    ig1 = id2;
    ig2 = id3;
  }
  else if( abs(id2) == 25 || abs(id2) == 35 || 
	   abs(id2) == 45 ||abs(id2) == 36 ||abs(id2) == 46 || abs(id2) == 37  ) {
    ihigg = id2;
    ig1 = id1;
    ig2 = id3;
  }
  else if( abs(id3) ==25 || abs(id3) == 35 || 
	   abs(id3) == 45 ||abs(id3) == 36 ||abs(id3) == 46 || abs(id3) == 37 ) {
    ihigg = id3;
    ig1 = id1;
    ig2 = id2;
  }
  else {
    throw HelicityConsistencyError() 
      << "NMSSMGOGOHVertex::setCoupling - There is no higgs particle in "
      << "this vertex. Particles: " << id1 << " " << id2 << " " << id3
      << Exception::runerror;
    return;
  }
  
  // weak coupling
  if(q2!=_q2last) {
    _couplast = weakCoupling(q2);
    _q2last = q2;
  }
  double rt = sqrt(0.5);
  // CP-even neutral higgs
  if(ihigg == 25 || ihigg == 35 || ihigg == 45) {
    int iloc = (ihigg - 25)/10;
    // chargino
    if(abs(ig1) == 1000024 || abs(ig1) == 1000037) {
      if( ig1 < 0 ) swap(ig1, ig2);
      int ic1 = abs(ig1)==1000024 ? 0 : 1;
      int ic2 = abs(ig2)==1000024 ? 0 : 1;
      Complex coupL = -_lambda*rt*conj((*_mixS)(iloc,2)*(*_mixU)(ic1,1)*(*_mixV)(ic2,1))
	-_couplast*rt*(conj((*_mixS)(iloc,0)*(*_mixU)(ic1,1)*(*_mixV)(ic2,0) +
			    (*_mixS)(iloc,1)*(*_mixU)(ic1,0)*(*_mixV)(ic2,1)));
      Complex coupR = -_lambda*rt*(*_mixS)(iloc,2)*(*_mixU)(ic2,1)*(*_mixV)(ic1,1)
	-_couplast*rt*((*_mixS)(iloc,0)*(*_mixU)(ic2,1)*(*_mixV)(ic1,0)+
		       (*_mixS)(iloc,1)*(*_mixU)(ic2,0)*(*_mixV)(ic1,1));
      left(coupL);
      right(coupR);
      norm(1.0);
    }
    // neutralino
    else  {
      int in1 = (ig1 < 1000024) ? (ig1 - 1000022) : (ig1 - 1000005)/10; 
      int in2 = (ig2 < 1000024) ? (ig2 - 1000022) : (ig2 - 1000005)/10;
      Complex us1 = (*_mixS)(iloc, 0), us2 = (*_mixS)(iloc, 1);
      Complex us3 = (*_mixS)(iloc, 2);
      Complex ni1 = (*_mixN)(in1,0), nj1 = (*_mixN)(in2,0);	  
      Complex ni2 = (*_mixN)(in1,1), nj2 = (*_mixN)(in2,1);
      Complex ni3 = (*_mixN)(in1,3), nj3 = (*_mixN)(in2,3); 
      Complex ni4 = (*_mixN)(in1,2), nj4 = (*_mixN)(in2,2); 
      Complex ni5 = (*_mixN)(in1,4), nj5 = (*_mixN)(in2,4);
      Complex YL =  
	- _lambda*rt*(us2*(ni4*nj5 + ni5*nj4) + 
		      us1*(ni3*nj5 + ni5*nj3) +
		      us3*(ni3*nj4 + ni4*nj3))
	+ sqrt(2.)*_kappa*us3*ni5*nj5 
	- _couplast*0.5*(us2*(ni2*nj3 + ni3*nj2) -
			 us1*(ni2*nj4 + ni4*nj2))
	+ _couplast*0.5*_sw*(us2*(ni1*nj3 + ni3*nj1) - 
			     us1*(ni1*nj4 + ni4*nj1) )/_cw;
      left(-conj(YL));
      right(-YL);
      norm(1.0);
    }
  }
  // CP-odd  neutral higgs
  else if(ihigg==36||ihigg==46) {
    int iloc = (ihigg-36)/10;
    // chargino
    if(abs(ig1)==1000024||abs(ig1)==1000037) {
	if( ig1 < 0 ) swap(ig1, ig2);
	int ic1 = abs(ig1)==1000024 ? 0 : 1;
	int ic2 = abs(ig2)==1000024 ? 0 : 1;
	Complex QL = Complex(0,-1.0)*
	  (_lambda*rt*conj((*_mixP)(iloc,2)*(*_mixU)(ic1,1)*(*_mixV)(ic2,1))
	   -_couplast*rt*conj(((*_mixP)(iloc,0)*(*_mixU)(ic1,1)*(*_mixV)(ic2,0) +
			       (*_mixP)(iloc,1)*(*_mixU)(ic1,0)*(*_mixV)(ic2,1))));
	Complex QR = Complex(0,-1.0)*
	  (_lambda*rt*(*_mixP)(iloc,2)*(*_mixU)(ic2,1)*(*_mixV)(ic1,1)
	   -_couplast*rt*((*_mixP)(iloc,0)*(*_mixU)(ic2,1)*(*_mixV)(ic1,0) +
			  (*_mixP)(iloc,1)*(*_mixU)(ic2,0)*(*_mixV)(ic1,1)));
	left(QL);
	right(-QR);
	norm(1.);
    }
    // neutralino
    else {
      int in1 = (ig1 < 1000024) ? (ig1 - 1000022) : (ig1 - 1000005)/10; 
      int in2 = (ig2 < 1000024) ? (ig2 - 1000022) : (ig2 - 1000005)/10;
      Complex up1 = (*_mixP)(iloc, 0), up2 = (*_mixP)(iloc, 1);
      Complex up3 = (*_mixP)(iloc, 2);
      Complex ni1 = (*_mixN)(in1,0), nj1 = (*_mixN)(in2,0);	  
      Complex ni2 = (*_mixN)(in1,1), nj2 = (*_mixN)(in2,1); 
      Complex ni3 = (*_mixN)(in1,2), nj3 = (*_mixN)(in2,2); 
      Complex ni4 = (*_mixN)(in1,3), nj4 = (*_mixN)(in2,3); 
      Complex ni5 = (*_mixN)(in1,4), nj5 = (*_mixN)(in2,4);
      Complex AL = 
	_lambda*rt*(up2*(ni3*nj5 + ni5*nj3) +
		    up1*(ni4*nj5 + ni5*nj4) + 
		    up3*(ni3*nj4 + ni4*nj3))
	- sqrt(2.)*_kappa*up3*ni5*nj5
	- _couplast*0.5*(up2*(ni2*nj4 + ni4*nj2) - 
			 up1*(ni2*nj3 + ni3*nj2))
	+ _couplast*0.5*_sw*(up2*(ni1*nj4 + ni4*nj1) - 
			     up1*(ni1*nj3 + ni3*nj1))/_cw;
      AL *= Complex(0.0, -1.0);
      left(conj(AL));
      right(AL);
      norm(1.);
    }
  }
  // charged higgs
  else {
    if (abs(ig1) == 1000024 || abs(ig1) == 1000037) swap (ig1,ig2);
    int in = (abs(ig1) < 1000024) ? (ig1-1000022) : (ig1-1000005)/10; 
    int ic = (abs(ig2) == 1000024) ? 0 : 1;
    Complex QpR = _lambda*_cosb*(*_mixU)(ic,1)*(*_mixN)(in,4)
      -_sinb*_couplast*(rt*(*_mixU)(ic,1)*(_sw*(*_mixN)(in,0)/_cw + (*_mixN)(in,1))
		 	 - (*_mixU)(ic,0)*(*_mixN)(in,2));
    Complex QpL = _lambda*_sinb*(*_mixV)(ic,1)*(*_mixN)(in,4)
      + _couplast*_cosb*(rt*(*_mixV)(ic,1)
			 *(_sw*(*_mixN)(in,0)/_cw + (*_mixN)(in,1))
			 + (*_mixV)(ic,0)*(*_mixN)(in,3));
    QpL = conj(QpL);
    if(ihigg > 0) {
      left (QpL);
      right(QpR);
      norm(-1.);
    }
    else {
      left (conj(QpR));
      right(conj(QpL));
      norm(-1.);
    }
  }
}
