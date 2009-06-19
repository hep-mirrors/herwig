// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMGGHVertex class.
//

#include "NMSSMGGHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/Susy/NMSSM/NMSSM.h"

using namespace Herwig;

NMSSMGGHVertex::NMSSMGGHVertex() : 
  _sw(0.), _cw(0.), _mw(0.*MeV), _mz(0.*MeV),_lambdaVEV(0.*MeV), _lambda(0.),
   _v1(0.*MeV), _v2(0.*MeV), _triTp(0.*MeV), _triBt(0.*MeV),
  _triTa(0.*MeV),_sb(0.), _cb(0.),
  _masslast(make_pair(0.*MeV,0.*MeV)), _q2last(0.*MeV2), _couplast(0.),_coup(0.),
  _hlast(0), _recalc(true) {
  
  //PDG codes for particles at vertices
  vector<long> first(5,21),second(5,21), third(5);
  third[0] = 25;
  third[1] = 35;
  third[2] = 36;
  third[3] = 45;
  third[4] = 46;
  setList(first, second, third);
}

void NMSSMGGHVertex::doinit() {
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if( !_theSM ) {
    throw InitException() << "NMSSMGGHVertex::doinit - The SM pointer is null!"
			  << Exception::abortnow;
  }
  // SM parameters
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt(1. - _theSM->sin2ThetaW());
  _mw = getParticleData(24)->mass();
  _mz = getParticleData(23)->mass();
  _top = getParticleData(6);
  _bt = getParticleData(5);
  _charm = getParticleData(4);
  _tau = getParticleData(15);
  
  //NMSSM parameters
  tcNMSSMPtr nmssm = dynamic_ptr_cast<tcNMSSMPtr>(_theSM);
  _mixS = nmssm->CPevenHiggsMix();
  _mixP = nmssm->CPoddHiggsMix();
  _mixQt = nmssm->stopMix();
  _mixQb = nmssm->sbottomMix();
  _mixQta = nmssm->stauMix();
  
  double beta = atan(nmssm->tanBeta());
  _sb = sin(beta);
  _cb = cos(beta);
  
 _v1 = sqrt(2.)*_mw*_cb;
 _v2 = sqrt(2.)*_mw*_sb;
 
   _lambda = nmssm->lambda();
  _lambdaVEV = nmssm->lambdaVEV();
 
 
	_triTp = nmssm->topTrilinear();
	_triBt = nmssm->bottomTrilinear();
    _triTa = nmssm->tauTrilinear();



  // resize vectors here and use setNParticles method
  // to the set the actual number in the loop.
  // Also only the top mass hass to be calculated at runtime
  masses.resize(9, Energy());
  masses[3] = getParticleData(1000005)->mass();
  masses[4] = getParticleData(2000005)->mass();
  masses[5] = getParticleData(1000006)->mass();
  masses[6] = getParticleData(2000006)->mass();
  masses[7] = getParticleData(1000015)->mass();  
  masses[8] = getParticleData(2000015)->mass();  
  
  type.resize(9, PDT::Spin0);
  type[0] = PDT::Spin1Half;
  type[1] = PDT::Spin1Half;
  type[2] = PDT::Spin1Half;
  couplings.resize(9);

  orderInGem(1);
  orderInGs(2);

  VVSLoopVertex::doinit();
}

void NMSSMGGHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << _sw << _cw << ounit(_mw, GeV) << ounit(_mz, GeV)
   << ounit(_lambdaVEV,GeV) << _lambda << ounit(_v1,GeV) << ounit(_v2,GeV)
    << ounit(_triTp,GeV) << ounit(_triBt,GeV) << ounit(_triTa,GeV) 
     << _top << _bt << _charm << _tau << _mixS << _mixP << _mixQt << _mixQb
	  << _mixQta << _sb << _cb; 
}


void NMSSMGGHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> _sw >> _cw >> iunit(_mw, GeV) >> iunit(_mz, GeV)
     >> iunit(_lambdaVEV,GeV) >> _lambda >> iunit(_v1,GeV) >> iunit(_v2,GeV)
    >> iunit(_triTp,GeV) >> iunit(_triBt,GeV) >> iunit(_triTa,GeV) 
     >> _top >> _bt >> _charm >> _tau >> _mixS >> _mixP >> _mixQt >>
	  _mixQb >> _mixQta >> _sb >> _cb; 
}

ClassDescription<NMSSMGGHVertex> NMSSMGGHVertex::initNMSSMGGHVertex;
// Definition of the static class description member.

void NMSSMGGHVertex::Init() {

  static ClassDocumentation<NMSSMGGHVertex> documentation
    ("The effective coupling of a higgs to a pair of gluons in the "
     "NMSSM.");

}

void NMSSMGGHVertex::setCoupling(Energy2 q2, tcPDPtr p1, tcPDPtr p2,
				 tcPDPtr p3) {			 
  long hid(p3->id());
  if( q2 != _q2last ) {
    _couplast = sqr(strongCoupling(q2));
	 _coup = weakCoupling(q2);
    _q2last = q2;
    _recalc = true;
  }
  setNorm(_couplast*_coup);

  if( hid != _hlast ) {
    _hlast = hid;
    _recalc = true;
    if( hid % 5 == 0 ) {
	
		Energy _mt= 0.*MeV;
Energy _mb= 0.*MeV;
Energy _mc= 0.*MeV;
Energy _mta= 0.*MeV;

      setNParticles(9);
      int iloc = (hid - 25)/10;
      //top quark
	   _mt = _theSM->mass(q2, _top);
      masses[0] = _theSM->mass(q2, _top);
      Complex cpl = -_mt*
	                (*_mixS)(iloc, 1)*0.5/_sb/_mw;
       couplings[0].first = cpl;
      couplings[0].second = cpl;
      //couplings[0].second = cpl;
	  //bottom quark
	  _mb = _theSM->mass(q2, _bt);
      masses[1] = _theSM->mass(q2, _bt);
      cpl = -_mb*(*_mixS)(iloc, 0)*0.5/_cb/_mw;
      couplings[1].first = cpl;
      couplings[1].second = cpl;	 
	  //charm quark
	   _mc = _theSM->mass(q2, _charm);
	  masses[2] = _theSM->mass(q2, _charm);
	  cpl = - _mc*(*_mixS)(iloc, 1)*0.5/_sb/_mw;
      couplings[2].first = cpl;
      couplings[2].second = cpl;
      //~b_1
     double f1 = _mb/_mw/_cb;
      Complex f2 = (0.5*_mz*(_sb*(*_mixS)(iloc,1) - _cb*(*_mixS)(iloc,0))/_cw)*UnitRemoval::InvE;
      
      cpl = -f2*( (1. - 2.*sqr(_sw)/3.)*(*_mixQb)(0, 0)*(*_mixQb)(0, 0)
	        + 2.*sqr(_sw)*(*_mixQb)(0, 1)*(*_mixQb)(0, 1)/3.)
			- f1*_mb*UnitRemoval::InvE*(*_mixS)(iloc,0)
			*((*_mixQb)(0, 0)*(*_mixQb)(0, 0) + (*_mixQb)(0, 1)*(*_mixQb)(0, 1)) 
	        - 0.5*f1*(_lambdaVEV*(*_mixS)(iloc,1) + _lambda*_v2*(*_mixS)(iloc,2)/_coup 
			+  _triBt*(*_mixS)(iloc,0))*((*_mixQb)(0, 1)*(*_mixQb)(0, 0)
			+ (*_mixQb)(0, 0)*(*_mixQb)(0, 1))* UnitRemoval::InvE;
			
		
		
      couplings[3].first = cpl; 
	  couplings[3].second = cpl; 
      //~b_2
      cpl =  -f2*((1. - 2.*sqr(_sw)/3.)*(*_mixQb)(1, 0)*(*_mixQb)(1, 0)
	         + 2.*sqr(_sw)*(*_mixQb)(1, 1)*(*_mixQb)(1, 1)/3.)
			- f1*_mb*UnitRemoval::InvE*(*_mixS)(iloc,0)*
			((*_mixQb)(1, 0)*(*_mixQb)(1, 0) + (*_mixQb)(1, 1)*(*_mixQb)(1, 1)) 
	        - 0.5*f1*(_lambdaVEV*(*_mixS)(iloc,1) + _lambda*_v2*(*_mixS)(iloc,2)/_coup
			+  _triBt*(*_mixS)(iloc,0))*((*_mixQb)(0, 1)*(*_mixQb)(0, 0)
			+ (*_mixQb)(0, 0)*(*_mixQb)(0, 1))*UnitRemoval::InvE;

			
      couplings[4].first = cpl; 
	  couplings[4].second = cpl; 
     
	   //~t_1
      f1 = _mt/_mw/_sb;

	  cpl   =-f2*( (1. - 4.*sqr(_sw)/3.)*(*_mixQt)(0, 0)*(*_mixQt)(0, 0)
	        + 4.*sqr(_sw)*(*_mixQt)(0, 1)*(*_mixQt)(0, 1)/3.)
	        - f1*_mt*UnitRemoval::InvE*(*_mixS)(iloc,1)
			*((*_mixQt)(0, 0)*(*_mixQt)(0, 0)
			+ (*_mixQt)(0, 1)*(*_mixQt)(0, 1))  
	        -  0.5*f1*(_lambdaVEV*(*_mixS)(iloc,0) + _lambda*_v1*(*_mixS)(iloc,2)/_coup
			+ _triTp*(*_mixS)(iloc,1))*((*_mixQt)(0, 1)*(*_mixQt)(0, 0)
			+ (*_mixQt)(0, 0)*(*_mixQt)(0, 1))*UnitRemoval::InvE;
	  
	   
		
      couplings[5].first = cpl;
	  couplings[5].second = cpl; 
      //~t_2
	 cpl =  -f2*( (1. - 4.*sqr(_sw)/3.)*(*_mixQt)(1, 0)*(*_mixQt)(1, 0)
	        + 4.*sqr(_sw)*(*_mixQt)(1, 1)*(*_mixQt)(1, 1)/3.)
	        - f1*_mt*UnitRemoval::InvE*(*_mixS)(iloc,1)
			*((*_mixQt)(1, 0)*(*_mixQt)(1, 0)
			+ (*_mixQt)(1, 1)*(*_mixQt)(1, 1))  
	        -  0.5*f1*(_lambdaVEV*(*_mixS)(iloc,0) + _lambda*_v1*(*_mixS)(iloc,2)/_coup
			+ _triTp*(*_mixS)(iloc,1))*((*_mixQt)(1, 1)*(*_mixQt)(1, 0)
			+ (*_mixQt)(1, 0)*(*_mixQt)(1, 1))*UnitRemoval::InvE;
		
      couplings[6].first = cpl;
	  couplings[6].second = cpl; 

	  //~tau_1
	    _mta = _theSM->mass(q2, _tau); 
	  f1 = _mta/_mw/_cb;
	  cpl =-f2*( (1. - 2.*sqr(_sw))*(*_mixQta)(0, 0)*(*_mixQta)(0, 0) 
		    + 2.*sqr(_sw)*(*_mixQta)(0, 1)*(*_mixQta)(0, 1))
		    - f1*_mta*UnitRemoval::InvE*(*_mixS)(iloc,0)
			*((*_mixQta)(0, 0)*(*_mixQta)(0, 0) 
			+ (*_mixQta)(0, 1)*(*_mixQta)(0, 1)) 
			- 0.5*f1*(_lambdaVEV*(*_mixS)(iloc,1) + _lambda*_v2*(*_mixS)(iloc,2)/_coup
			+  _triTa*(*_mixS)(iloc,0))*((*_mixQta)(0, 1)*(*_mixQta)(0, 0) 
		    + (*_mixQta)(0, 0)*(*_mixQta)(0, 1)) * UnitRemoval::InvE;
		
      couplings[7].first = cpl; 
	  couplings[7].second = cpl; 
      //~tau_2
	 cpl =  -f2*( (1. - 2.*sqr(_sw))*(*_mixQta)(1, 0)*(*_mixQta)(1, 0)
	        + 2.*sqr(_sw)*(*_mixQta)(1, 1)*(*_mixQta)(1, 1))
	        - f1*_mta*UnitRemoval::InvE*(*_mixS)(iloc,1)
			*((*_mixQta)(1, 0)*(*_mixQta)(1, 0)
			+ (*_mixQta)(1, 1)*(*_mixQta)(1, 1))  
	        -  0.5*f1*(_lambdaVEV*(*_mixS)(iloc,0) 
			+ _lambda*_v1*(*_mixS)(iloc,2)/_coup
			+ _triTa*(*_mixS)(iloc,1) )*((*_mixQta)(1, 1)*(*_mixQta)(1, 0)
			+ (*_mixQta)(1, 0)*(*_mixQta)(1, 1))*UnitRemoval::InvE;
		
		
      couplings[8].first = cpl; 
	  couplings[8].second = cpl; 	  
	  

    }
	
	
	
    else {
      setNParticles(9);
      int iloc = (hid - 36)/10;
	  	Energy _mt= 0.*MeV;
Energy _mb= 0.*MeV;
Energy _mc= 0.*MeV;
Energy _mta= 0.*MeV;
      //top quark
	   _mt =  _theSM->mass(q2, _top);
      masses[0] = _theSM->mass(q2, _top);
	  Complex cpl = Complex(0., -1.)*0.5*_mt* 
					(*_mixP)(iloc, 1)/_sb/_mw;
      couplings[0].first = cpl;
      couplings[0].second = -cpl;
	  //bottom quark
      masses[1] = _theSM->mass(q2, _bt);
	  _mb = _theSM->mass(q2, _bt);
       cpl = Complex(0., -1.)*0.5*_mb*
	        (*_mixP)(iloc, 0)/_cb/_mw;
      couplings[1].first = cpl;
      couplings[1].second = -cpl;	 
	  //charm quark
	  masses[2] = _theSM->mass(q2, _charm);
	   _mc = _theSM->mass(q2, _charm);
       cpl = Complex(0., -1.)*0.5*_mc*
	       (*_mixP)(iloc, 1)/_sb/_mw;
      couplings[2].first = cpl;
      couplings[2].second = -cpl;
        //sbottom_1
      double f1 = _mb/_mw/_cb;
      Complex f2 = 0.5*_mz*(-_cb*(*_mixP)(iloc,0) 
	               + _sb*(*_mixP)(iloc,1))/_cw*UnitRemoval::InvE;
      cpl = -f2*((1. - 2.*sqr(_sw)/3.)
	          *(*_mixQb)(0, 0)*(*_mixQb)(0, 0)
	        + 2.*sqr(_sw)*(*_mixQb)(0, 1)*(*_mixQb)(0, 1)/3.)
			- f1*_mb*UnitRemoval::InvE*(*_mixP)(iloc,0)
			*((*_mixQb)(0, 0)*(*_mixQb)(0, 0)+ (*_mixQb)(0, 1)*(*_mixQb)(0, 1)) 
	        + 0.5*f1*Complex((_lambdaVEV*(*_mixP)(iloc,1)
			 + _lambda*_v2*(*_mixP)(iloc,2)/_coup
			-  _triBt*(*_mixP)(iloc,0))*((*_mixQb)(0, 1)*(*_mixQb)(0, 0)
			+ (*_mixQb)(0, 0)*(*_mixQb)(0, 1)) * UnitRemoval::InvE);
			
		cpl *= Complex(0.,1.0);	 
		 couplings[3].first = cpl;
         couplings[3].second = -cpl;
		 
	
      //~b_2
      cpl =  -f2*( (1. - 2.*sqr(_sw)/3.)
	         *(*_mixQb)(1, 0)*(*_mixQb)(1, 0)
	         + 2.*sqr(_sw)*(*_mixQb)(1, 1)*(*_mixQb)(1, 1)/3.)
			- f1*_mb*UnitRemoval::InvE*(*_mixP)(iloc,0)
			 *((*_mixQb)(1, 0)*(*_mixQb)(1, 0) 
			+ (*_mixQb)(1, 1)*(*_mixQb)(1, 1)) 
			+ 0.5*f1*Complex((_lambdaVEV*(*_mixP)(iloc,1)
			+ _lambda*_v2*(*_mixP)(iloc,2)/_coup
			-  _triBt*(*_mixP)(iloc,0))*((*_mixQb)(1, 1)*(*_mixQb)(1, 0)
			+ (*_mixQb)(1, 0)*(*_mixQb)(1, 1)) * UnitRemoval::InvE);
		 
	cpl *=Complex(0.,1.0);			
			
      couplings[4].first = cpl; 
	  couplings[4].second = -cpl; 
     
	   //~t_1
      f1 = _mt/_mw/_sb;

	  cpl   = -f2*( (1. - 4.*sqr(_sw)/3.)
	          *(*_mixQt)(0, 0)*(*_mixQt)(0, 0)
		    + 4.*sqr(_sw)*(*_mixQt)(0, 1)*(*_mixQt)(0, 1)/3.)
			- f1*_mt*UnitRemoval::InvE*(*_mixP)(iloc,1)
			*((*_mixQt)(0, 0)*(*_mixQt)(0, 0)
	        + (*_mixQt)(0, 1)*(*_mixQt)(0, 1))  
	        +  0.5*f1*Complex((_lambdaVEV*(*_mixP)(iloc,0) 
			+ _lambda*_v1*(*_mixP)(iloc,2)/_coup
			- _triTp*(*_mixP)(iloc,1) )*((*_mixQt)(0, 1)*(*_mixQt)(0, 0) 
		    + (*_mixQt)(0, 0)*(*_mixQt)(0, 1))*UnitRemoval::InvE);
			
	cpl *= Complex(0.,1.0);		
		
      couplings[5].first = cpl;
	  couplings[5].second = -cpl; 
      //~t_2
	 cpl =  -f2*( (1. - 4.*sqr(_sw)/3.)
	        *(*_mixQt)(1, 0)*(*_mixQt)(1, 0)
	        + 4.*sqr(_sw)*(*_mixQt)(1, 1)*(*_mixQt)(1, 1)/3.)
	        - f1*_mt*UnitRemoval::InvE*(*_mixP)(iloc,1)
			*((*_mixQt)(1, 0)*(*_mixQt)(1, 0)
			+ (*_mixQt)(1, 1)*(*_mixQt)(1, 1))  
	        +  0.5*f1*Complex((_lambdaVEV*(*_mixP)(iloc,0) 
			+ _lambda*_v1*(*_mixS)(iloc,2)/_coup
			- _triTp*(*_mixP)(iloc,1) )*((*_mixQt)(1, 1)*(*_mixQt)(1, 0)
			+ (*_mixQt)(1, 0)*(*_mixQt)(1, 1))*UnitRemoval::InvE);
			
	cpl *= Complex(0.,1.0);
		
      couplings[6].first = cpl;
	  couplings[6].second = -cpl; 

	  //stau_1
	   _mta=_theSM->mass(q2, _tau); 
	  f1 = _mta/_mw/_cb;
	  cpl = -f2*( (1. - 2.*sqr(_sw))
	        *(*_mixQta)(0, 0)*(*_mixQta)(0, 0) 
		    + 2.*sqr(_sw)*(*_mixQta)(0, 1)*(*_mixQta)(0, 1))
		    - f1*_mta*UnitRemoval::InvE*(*_mixP)(iloc,0)
			*((*_mixQta)(0, 0)*(*_mixQta)(0, 0) 
			+ (*_mixQta)(0, 1)*(*_mixQta)(0, 1)) 
			+ 0.5*f1*Complex(( _lambdaVEV*(*_mixP)(iloc,1) 
			+ _lambda*_v2*(*_mixP)(iloc,2)/_coup
			-  _triTa*(*_mixP)(iloc,0))*((*_mixQta)(0, 1)*(*_mixQta)(0, 0) 
		    + (*_mixQta)(0, 0)*(*_mixQta)(0, 1))* UnitRemoval::InvE);
			
		cpl *=Complex(0.,1.0);
		
      couplings[7].first = cpl; 
	  couplings[7].second = -cpl; 
      //~tau_2
	 cpl =  -f2*( (1. - 2.*sqr(_sw))
	        *(*_mixQta)(1, 0)*(*_mixQta)(1, 0)
	        + 2.*sqr(_sw)*(*_mixQta)(1, 1)*(*_mixQta)(1, 1))
	        - f1*_mta*UnitRemoval::InvE*(*_mixP)(iloc,1)
			*((*_mixQta)(1, 0)*(*_mixQta)(1, 0)
			+ (*_mixQta)(1, 1)*(*_mixQta)(1, 1))  
	        +  0.5*f1*Complex(( _lambdaVEV*(*_mixP)(iloc,0) 
			+ _lambda*_v1*(*_mixP)(iloc,2)/_coup
			- _triTa*(*_mixP)(iloc,1) )*((*_mixQta)(1, 1)*(*_mixQta)(1, 0)
			+ (*_mixQta)(1, 0)*(*_mixQta)(1, 1))*UnitRemoval::InvE);
			
	cpl *= Complex(0.,1.0);
		
		
      couplings[8].first = cpl; 
	  couplings[8].second = -cpl; 		
		
    }
	
  }

  if( _recalc ) {
    VVSLoopVertex::setCoupling(q2, p1, p2, p3);
    _recalc = false;
  } 
}
