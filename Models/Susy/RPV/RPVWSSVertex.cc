// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVWSSVertex class.
//

#include "RPVWSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "RPV.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

RPVWSSVertex::RPVWSSVertex() :_sw(0.), _cw(0.), _s2w(0.), _c2w(0.),
			      _interactions(0), _q2last(), 
			      _ulast(0), _dlast(0), _gblast(0),
			      _factlast(0.), _couplast(0.) {
  orderInGs(0);
  orderInGem(1);
}				 

IBPtr RPVWSSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVWSSVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVWSSVertex::doinit() {
  // extract the model 
  tRPVPtr model = dynamic_ptr_cast<tRPVPtr>(generator()->standardModel());
  if( !model ) throw InitException() << "RPVWSSVertex::doinit() - The"
				     << " pointer to the RPV object is null!"
				     << Exception::abortnow;
  // get the mixing matrices
  MixingMatrixPtr mixH = model->CPevenHiggsMix() ;
  MixingMatrixPtr mixP = model->CPoddHiggsMix()  ;
  MixingMatrixPtr mixC = model->ChargedHiggsMix();
  // find the codes for the Higgs bosons
  vector<long> pseudo(1,36);
  vector<long> scalar(2);
  vector<long> charged(1,37);
  scalar[0] = 25;
  scalar[1] = 35;
  if(mixH&&mixH->size().first>2) {
    scalar.push_back(1000012);
    scalar.push_back(1000014);
    scalar.push_back(1000016);
  }
  if(mixP&&mixP->size().first>1) {
    pseudo.push_back(1000017);
    pseudo.push_back(1000018);
    pseudo.push_back(1000019);
  }
  if(mixC&&mixC->size().first>2) {
    charged.push_back(-1000011);
    charged.push_back(-1000013);
    charged.push_back(-1000015);
    charged.push_back(-2000011);
    charged.push_back(-2000013);
    charged.push_back(-2000015);
  }
  // sfermion interactions
  if(_interactions==0||_interactions==1) {
    // squarks
    //W-
    //LL-squarks
    for(long ix=1000001;ix<1000006;ix+=2) {
      addToList(-24,ix+1,-ix);
    }
    //1-2 stop sbottom
    addToList(-24,1000006,-2000005);
    //2-1 stop sbottom
    addToList(-24,2000006,-1000005);
    //2-2 stop sbottom
    addToList(-24,2000006,-2000005);
    //W+
    for(long ix=1000001;ix<1000006;ix+=2) {
      addToList(24,-(ix+1),ix);
    }
    //1-2 stop sbottom
    addToList(24,-1000006,2000005);
    //2-1 stop sbottom
    addToList(24,-2000006,1000005);
    //2-2 stop sbottom
    addToList(24,-2000006,2000005);
    //LL squarks
    for(long ix=1000001;ix<1000007;++ix) {
      addToList(23,ix,-ix);
    }
    //RR squarks
    for(long ix=2000001;ix<2000007;++ix) {
      addToList(23,ix,-ix);
    }
    //L-Rbar stop
    addToList(23,1000006,-2000006);
    //Lbar-R stop
    addToList(23,-1000006,2000006);
    //L-Rbar sbottom
    addToList(23,1000005,-2000005);
    //Lbar-R sbottom
    addToList(23,-1000005,2000005);
    // gamma
    //squarks
    for(long ix=1000001;ix<1000007;++ix) {
      addToList(22,ix,-ix);
    }
    for(long ix=2000001;ix<2000007;++ix) {
      addToList(22,ix,-ix);
    }
    //sleptons
    // gamma
    for(long ix=1000011;ix<1000016;ix+=2) {
      addToList(22,ix,-ix);
    }
    for(long ix=2000011;ix<2000016;ix+=2) {
      addToList(22,ix,-ix);
    }
    // Z
    //LL-sleptons
    for(long ix=1000011;ix<1000017;ix+=2) {
      addToList(23,ix,-ix);
    }
    // RR-sleptons
    for(long ix=2000011;ix<2000016;ix+=2) {
      addToList(23,ix,-ix);
    }
    //L-Rbar stau
    addToList(23,1000015,-2000015);
    //Lbar-R stau
    addToList(23,-1000015,2000015);
    if(!(mixH&&mixH->size().first>2)) {
      for(long ix=1000012;ix<1000017;ix+=2) {
	// sneutrinos
	addToList(23,ix,-ix);
      }
      //LL-sleptons
      for(long ix=1000011;ix<1000016;ix+=2) {
	addToList(-24,-ix,ix+1);
      }
      //2-L stau
      addToList(-24,-2000015,1000016);
      //LL-sleptons
      for(long ix=1000011;ix<1000016;ix+=2) {
	addToList(24,ix,-ix-1);
      }
      //2-L stau
      addToList(24,2000015,-1000016);
    }
  }
  if(_interactions==0||_interactions==2) {
    // charged Higgs and photon
    addToList(22,37,-37);
    // charged Higgs and Z
    for(unsigned int ix=0;ix<charged.size();++ix) {
      for(unsigned int iy=0;iy<charged.size();++iy) {
	if(abs(charged[ix])>1000000&&abs(charged[iy])>100000 && 
	   ( charged[ix]==charged[iy] ||
	     (abs(charged[ix])%1000000==15&&abs(charged[iy])%1000000==15)))
	  continue;
	addToList(23,charged[ix],-charged[iy]);
      }
    }
    // neutral Higgs and Z
    for(unsigned int ix=0;ix<scalar.size();++ix) {
      for(unsigned int iy=0;iy<pseudo.size();++iy) {
	addToList(23,scalar[ix],pseudo[iy]);
      }
    }
    // charged higss, scalar higgs and W
    for(unsigned int ix=0;ix<charged.size();++ix) {
      for(unsigned int iy=0;iy<scalar.size();++iy) {
	addToList( 24,-charged[ix],scalar[iy]);
	addToList(-24, charged[ix],scalar[iy]);
      }
    }
    // charged higss, pseudoscalar higgs and W
    for(unsigned int ix=0;ix<charged.size();++ix) {
      for(unsigned int iy=0;iy<pseudo.size();++iy) {
	addToList( 24,-charged[ix],pseudo[iy]);
	addToList(-24, charged[ix],pseudo[iy]);
      }
      for(unsigned int iy=0;iy<scalar.size();++iy) {
	addToList( 24,-charged[ix],scalar[iy]);
	addToList(-24, charged[ix],scalar[iy]);
      }
    }
  }
  VSSVertex::doinit();
  // sfermion mixing
  _stop = model->stopMix();
  _sbottom = model->sbottomMix();
  if(!_stop || !_sbottom)
    throw InitException() << "RPVWSSVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " stop: " << _stop << " sbottom: " << _sbottom
			  << Exception::abortnow;
  _stau = model->stauMix();
  if(!_stau && (!mixC || mixC->size().first<2))
    throw InitException() << "RPVWSSVertex::doinit() either the stau"
			  << " mixing matrix must be set or the stau"
			  << " included in mixing with the"
			  << " charged Higgs bosons" << Exception::abortnow;
  // weak mixing
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt( 1. - sin2ThetaW() );
  _s2w = 2.*_cw*_sw;
  _c2w = 1. - 2.*sin2ThetaW();

  // Coupling of Z to scalar and pseudoscalar
  Cijeo_.resize(pseudo.size(),vector<Complex>(scalar.size(),0.));
  for(unsigned int i=0;i<pseudo.size();++i) {
    for(unsigned int j=0;j<scalar.size();++j) {
      for(unsigned int ix=0;ix<mixH->size().second;++ix) {
	double sign = ix!=1 ? 1. : -1.;
	Cijeo_[i][j] += sign*(*mixP)(i,ix)*(*mixH)(j,ix);
      }
    }
  }
  // Coupling of W to scalar charged
  Cijec_.resize(scalar.size(),vector<Complex>(charged.size(),0.));
  for(unsigned int i=0;i<scalar.size();++i) {
    for(unsigned int j=0;j<charged.size();++j) {
      for(unsigned int ix=0;ix<mixH->size().second;++ix) {
	double sign = ix!=1 ? 1. : -1.;
	Cijec_[i][j] += sign*(*mixH)(i,ix)*(*mixC)(j,ix);
      }
    }
  }
  // Coupling of W to pseudopseudo charged
  Cijco_.resize(pseudo.size(),vector<Complex>(charged.size(),0.));
  for(unsigned int i=0;i<pseudo.size();++i) {
    for(unsigned int j=0;j<charged.size();++j) {
      for(unsigned int ix=0;ix<mixP->size().second;++ix) {
	// not sure about this need it to get agreement with SPheno
	//double sign = ix!=1 ? 1. : -1.;
	double sign = 1.;
	Cijco_[i][j] += sign*(*mixP)(i,ix)*(*mixC)(j,ix);
      }
    }
  }
  // Coupling of Z to charged Higgs
  Cijc_.resize(charged.size(),vector<Complex>(charged.size(),0.));
  for(unsigned int i=0;i<charged.size();++i) {
    for(unsigned int j=0;j<charged.size();++j) {
      for(unsigned int ix=5;ix<mixC->size().second;++ix) {
	Cijc_[i][j] += (*mixC)(i,ix)*(*mixC)(j,ix);
      }
    }
  }
}

void RPVWSSVertex::persistentOutput(PersistentOStream & os) const {
  os << _interactions << _sw  << _cw
     << _stau << _stop << _sbottom << _s2w << _c2w
     << Cijeo_ << Cijec_ << Cijco_ << Cijc_;
}

void RPVWSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> _interactions >> _sw >> _cw
     >> _stau >> _stop >> _sbottom >> _s2w >> _c2w
     >> Cijeo_ >> Cijec_ >> Cijco_ >> Cijc_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RPVWSSVertex,Helicity::VSSVertex>
describeHerwigRPVWSSVertex("Herwig::RPVWSSVertex", "HwSusy.so HwRPV.so");

void RPVWSSVertex::Init() {

  static ClassDocumentation<RPVWSSVertex> documentation
    ("There is no documentation for the RPVWSSVertex class");

  static Switch<RPVWSSVertex,unsigned int> interfaceInteractions
    ("Interactions",
     "Which interactions to include",
     &RPVWSSVertex::_interactions, 0, false, false);
  static SwitchOption interfaceInteractionsAll
    (interfaceInteractions,
     "All",
     "Include both the interactions which would have been sfermion"
     " and Higgs bosons with the gauge bosons in the MSSM",
     0);
  static SwitchOption interfaceInteractionsSfermions
    (interfaceInteractions,
     "Sfermions",
     "Include the sfermion interactions",
     1);
  static SwitchOption interfaceInteractionsHiggs
    (interfaceInteractions,
     "Higgs",
     "Include the Higgs boson interactions",
     2);

}


void RPVWSSVertex::setCoupling(Energy2 q2,tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3) {
  long gboson = part1->id();
  assert(     gboson  == ParticleID::Z0    ||
	      gboson  == ParticleID::gamma || 
	  abs(gboson) == ParticleID::Wplus );
  long h1ID = part2->id();
  long h2ID = part3->id();
  // squarks and sleptons
  if(((abs(h1ID)>=1000001&&abs(h1ID)<=1000006)||
      (abs(h1ID)>=2000001&&abs(h1ID)<=2000006)) ||
     (((abs(h1ID)>=1000011&&abs(h1ID)<=1000016)||
       (abs(h1ID)>=2000011&&abs(h1ID)<=2000016))&&
      Cijeo_.size()==1)) {
    long sf1(abs(part2->id())),sf2(abs(part3->id()));
    if( sf1 % 2 != 0 ) swap(sf1, sf2);
    if( sf1 != _ulast || sf2 != _dlast || gboson != _gblast) {
      _gblast = gboson;
      _ulast = sf1;
      _dlast = sf2;
      //photon is simplest
      if( gboson == ParticleID::gamma )
	_factlast = getParticleData(sf1)->charge()/eplus;
      else {
	//determine which helicity state
	unsigned int alpha(sf1/1000000 - 1), beta(sf2/1000000 - 1);
	//mixing factors
	Complex m1a(0.), m1b(0.);
	if( sf1 == ParticleID::SUSY_t_1 || sf1 == ParticleID::SUSY_t_2 )
	  m1a = (*_stop)(alpha, 0);
	else if( sf1 == ParticleID::SUSY_b_1 || sf1 == ParticleID::SUSY_b_2 )
	  m1a = (*_sbottom)(alpha, 0);
	else if( sf1 == ParticleID::SUSY_tau_1minus || 
		 sf1 == ParticleID::SUSY_tau_2minus )
	  m1a = (*_stau)(alpha, 0);
	else
	  m1a = (alpha == 0) ? Complex(1.) : Complex(0.);
	if( sf2 == ParticleID::SUSY_t_1 || sf2 == ParticleID::SUSY_t_2 )
	  m1b = (*_stop)(beta, 0);
	else if( sf2 == ParticleID::SUSY_b_1 || sf2 == ParticleID::SUSY_b_2 )
	  m1b = (*_sbottom)(beta, 0);
	else if( sf2 == ParticleID::SUSY_tau_1minus || 
		 sf2 == ParticleID::SUSY_tau_2minus )
	  m1b = (*_stau)(beta, 0);
	else
	  m1b = (beta == 0) ? Complex(1.) : Complex(0.);
	//W boson
	if( abs(gboson) == ParticleID::Wplus ) {
	  _factlast = m1a*m1b/sqrt(2)/_sw;
	}
	//Z boson
	else {
	  if( sf1 == ParticleID::SUSY_nu_eL || sf1 == ParticleID::SUSY_nu_muL ||
	      sf1 == ParticleID::SUSY_nu_tauL ) {
	    _factlast = 1./_cw/2./_sw;
	  }
	  else {
	    double lmda(1.);
	    if( sf2 % 2 == 0 ) lmda = -1.;
	    _factlast = lmda*m1a*m1b;
	    if( alpha == beta) {
	      double ef = getParticleData(sf1)->charge()/eplus;
	      _factlast += 2.*ef*sqr(_sw);
	    }
	    _factlast *= -1./2./_cw/_sw; 
	  }
	}
      }
    }
  }
  // Higgs bosons
  else {
    _gblast = gboson;
    _ulast = h1ID;
    _dlast = h2ID;
    _factlast = 0.;
    if( gboson == ParticleID::Z0 ) {
      if( part2->charged() ) {
	unsigned int c1 = abs(h1ID) < 1000000 ? 0 : 
	  (abs(h1ID) < 2000000 ? (abs(h1ID)-1000009)/2 : (abs(h1ID)-2000003)/2);
	unsigned int c2 = abs(h2ID) < 1000000 ? 0 : 
	  (abs(h2ID) < 2000000 ? (abs(h2ID)-1000009)/2 : (abs(h2ID)-2000003)/2);
	if(c1==c2) _factlast = (_c2w-Cijc_[c1][c2])/_s2w;
	else       _factlast = -Cijc_[c1][c2]/_s2w;
	if(part2->iCharge()<0) _factlast *= -1.;
      }
      else {
	if(h1ID == ParticleID::h0         || h1ID  == ParticleID::H0         ||
	   h1ID == ParticleID::SUSY_nu_eL || h1ID == ParticleID::SUSY_nu_muL ||
	   h1ID == ParticleID::SUSY_nu_tauL ) {
	  unsigned int is = h1ID < 1000000 ? (h1ID-25)/10 : (h1ID-1000008)/2;
	  unsigned int ip = h2ID < 1000000 ? 0 : (h2ID-1000016);
	  _factlast =  Complex(0.,1.)*Cijeo_[ip][is]/_s2w;
	}
	else {
	  unsigned int is = h2ID < 1000000 ? (h2ID-25)/10 : (h2ID-1000008)/2;
	  unsigned int ip = h1ID < 1000000 ? 0 : (h1ID-1000016);
	  _factlast = -Complex(0.,1.)*Cijeo_[ip][is]/_s2w;
	}

      }
    }
    else if( gboson == ParticleID::gamma ) {
      _factlast = part2->iCharge()/3;
    }
    else {
      long scalar  = part2->charged() ? h2ID : h1ID;
      long charged = part2->charged() ? h1ID : h2ID;
      unsigned int ic = abs(charged) < 1000000 ? 0 : 
	(abs(charged) < 2000000 ? (abs(charged)-1000009)/2 : (abs(charged)-2000003)/2);
      if(scalar == ParticleID::h0         || scalar  == ParticleID::H0         ||
	 scalar == ParticleID::SUSY_nu_eL || scalar == ParticleID::SUSY_nu_muL ||
	 scalar == ParticleID::SUSY_nu_tauL ) {
	unsigned int ih = scalar < 1000000 ? (scalar-25)/10 : (scalar-1000008)/2;
	_factlast = -0.5*Cijec_[ih][ic]/_sw;
	if(gboson<0) _factlast *= -1.;
      }
      else {
	unsigned int ih = scalar < 1000000 ? 0 : (scalar-1000016);
	_factlast =  Complex(0., 0.5)*Cijco_[ih][ic]/_sw;
      } 
      if(part3->charged()) _factlast *= -1.;
    }
  }
  if( q2 != _q2last || _couplast==0. ) {
    _q2last = q2;
    _couplast = electroMagneticCoupling(q2);
  }
  norm(_couplast*_factlast);
}
