// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVSSSVertex class.
//

#include "RPVSSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "RPV.h"
#include <cassert>

using namespace Herwig;
using namespace ThePEG::Helicity;

RPVSSSVertex::RPVSSSVertex() : interactions_(0), q2Last_(ZERO), gLast_(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr RPVSSSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVSSSVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVSSSVertex::persistentOutput(PersistentOStream & os) const {
  os << interactions_ << stau_ << sbottom_ << stop_ << ounit(scalarScalarScalar_,GeV)
     << ounit(scalarPseudoPseudo_,GeV) << ounit(scalarSup_,GeV) << ounit(scalarSdown_,GeV)
     << ounit(scalarSneutrino_,GeV) << ounit(scalarSlepton_,GeV)
     << ounit(chargedSquark_,GeV) << ounit(chargedSlepton_,GeV)
     << ounit(scalarChargedCharged_,GeV) << ounit(pseudoChargedCharged_,GeV)
     << ounit(pseudoSup_,GeV) << ounit(pseudoSdown_,GeV) << ounit(pseudoSlepton_,GeV);
}

void RPVSSSVertex::persistentInput(PersistentIStream & is, int) {
  is >> interactions_ >> stau_ >> sbottom_ >> stop_ >> iunit(scalarScalarScalar_,GeV)
     >> iunit(scalarPseudoPseudo_,GeV) >> iunit(scalarSup_,GeV) >> iunit(scalarSdown_,GeV)
     >> iunit(scalarSneutrino_,GeV) >> iunit(scalarSlepton_,GeV)
     >> iunit(chargedSquark_,GeV) >> iunit(chargedSlepton_,GeV)
     >> iunit(scalarChargedCharged_,GeV) >> iunit(pseudoChargedCharged_,GeV)
     >> iunit(pseudoSup_,GeV) >> iunit(pseudoSdown_,GeV) >> iunit(pseudoSlepton_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RPVSSSVertex,SSSVertex>
describeHerwigRPVSSSVertex("Herwig::RPVSSSVertex", "HwRPV.so");

void RPVSSSVertex::Init() {

  static ClassDocumentation<RPVSSSVertex> documentation
    ("The RPVSSSVertex class implements the coupling of three"
     " scalar particles in the RPV model");

  static Switch<RPVSSSVertex,unsigned int> interfaceInteractions
    ("Interactions",
     "Which interactions to include",
     &RPVSSSVertex::interactions_, 0, false, false);
  static SwitchOption interfaceInteractionsAll
    (interfaceInteractions,
     "All",
     "Include both triple Higgs and Higgs sfermion interactions",
     0);
  static SwitchOption interfaceInteractionsHiggsHiggsHiggs
    (interfaceInteractions,
     "HiggsHiggsHiggs",
     "Only include triple Higgs boson interactions",
     1);
  static SwitchOption interfaceInteractionsHiggsSfermions
    (interfaceInteractions,
     "HiggsSfermions",
     "Only include Higgs sfermion interactions",
     2);

}

void  RPVSSSVertex::doinit() {
  // extract the model 
  tRPVPtr model = dynamic_ptr_cast<tRPVPtr>(generator()->standardModel());
  if( !model ) throw InitException() << "RPVSSSVertex::doinit() - The"
				     << " pointer to the RPV object is null!"
				     << Exception::abortnow;
  // get the Higgs mixing matrices
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
  // triple Higgs interactions
  if(interactions_==0||interactions_==1) {
    // three neutral scalar bosons
    for(unsigned int i1=0;i1<scalar.size();++i1) {
      for(unsigned int i2=0;i2<=i1;++i2) {
	for(unsigned int i3=0;i3<=i2;++i3) {
	  addToList(scalar[i1],scalar[i2],scalar[i3]);
	}
      }
    }
    // one neutral scalar two pseudoscalars
    for(unsigned int i1=0;i1<scalar.size();++i1) {
      for(unsigned int i2=0;i2<pseudo.size();++i2) {
	for(unsigned int i3=0;i3<=i2;++i3) {
	  addToList(scalar[i1],pseudo[i2],pseudo[i3]);
	}
      }
    }
    // one neutral scalar two charged scalars
    for(unsigned int i1=0;i1<scalar.size();++i1) {
      for(unsigned int i2=0;i2<charged.size();++i2) {
	for(unsigned int i3=0;i3<charged.size();++i3) {
	  if(!(abs(charged[i2])>1000000&&abs(charged[i3])>1000000&&
	     abs(charged[i2])%1000000==abs(charged[i3])%1000000))
	    addToList(scalar[i1],charged[i2],-charged[i3]);
	}
      }
    }
    // one pseudo scalar two charged scalars
    if(charged.size()>1) {
      for(unsigned int i1=0;i1<pseudo.size();++i1) {
	for(unsigned int i2=0;i2<charged.size();++i2) {
	  for(unsigned int i3=0;i3<charged.size();++i3) {
	    if(i2==i3) continue;
	    if(!(abs(charged[i2])>1000000&&abs(charged[i3])>1000000&&
		 abs(charged[i2])%1000000==abs(charged[i3])%1000000))
	      addToList(pseudo[i1],charged[i2],-charged[i3]);
	  }
	}
      }
    }
  }
  // sfermion interactions
  if(interactions_==0||interactions_==2) {
    // scalar neutral higgs sfermion sfermion 
    for(unsigned int i = 0; i < scalar.size(); ++i) {
      // squarks
      for(unsigned int j = 1; j < 7; ++j) {
	long lj = 1000000 + j;
	long rj = 2000000 + j;
	//LLbar
	addToList(scalar[i],lj,-lj);
	//RRbar
	addToList(scalar[i],rj,-rj);
	//LRbar
	addToList(scalar[i],lj,-rj);
	//RLbar
	addToList(scalar[i],rj,-lj);
      }
      // sleptons
      for(unsigned int j = 11; j < 17; ++j) {
	long lj = 1000000 + j;
	long rj = 2000000 + j;
	if( j % 2 != 0) {
	  // LL
	  addToList(scalar[i],lj,-lj);
	  // RR
	  addToList(scalar[i],rj,-rj);
	  //LRbar
	  addToList(scalar[i],lj,-rj);
	  //RLbar
	  addToList(scalar[i],rj,-lj);
	}
	// sneutrino only if no mixing
	else if(scalar.size()==2) {
	  addToList(scalar[i],lj,-lj);
	}
      }
    }
    // pseudoscalar sfermion sfermion
    for(unsigned int i = 0; i < pseudo.size(); ++i) {
      // squarks
      for(unsigned int j = 1; j < 7; ++j) {
	long lj = 1000000 + j;
	long rj = 2000000 + j;
	//LRbar
	addToList(pseudo[i],lj,-rj);
	//RLbar
	addToList(pseudo[i],rj,-lj);
      }
      // sleptons
      if(scalar.size()==2) {
	for(unsigned int j = 11; j < 17; j += 2) {
	  long lj = 1000000 + j;
	  long rj = 2000000 + j;
	  addToList(36,lj,-rj);
	  addToList(36,rj,-lj);
	}
      }
    }
    //outgoing H+
    for(unsigned int i=0;i<charged.size();++i) {
      // squarks
      for(long ii = 2; ii < 7; ii += 2) {
	//LL
	addToList(charged[i], 999999 + ii, -1000000 - ii);
	//RR
	addToList(charged[i], 1999999 + ii, -2000000 - ii);
	//RL
	addToList(charged[i], 1999999 + ii, -1000000 - ii);
	//LR
	addToList(charged[i], 999999 + ii, -2000000 - ii);
      }
      if(scalar.size()==2) {
	for(long ii = 11; ii < 17; ii += 2) {
	  addToList(37, 1000000 + ii, -1000001 - ii);
	  addToList(37, 2000000 + ii, -1000001 - ii);
	}
      }
      //outgoing H-
      for(long ii = 2; ii < 7; ii += 2) {
	//LL
	addToList(-charged[i], 1000000 + ii, -999999 - ii);
	//RR
	addToList(-charged[i], 2000000 + ii, -1999999 - ii);
	//RL
	addToList(-charged[i], 1000000 + ii, -1999999 - ii);
	//LR
	addToList(-charged[i], 2000000 + ii, -999999 - ii);
      }
      if(scalar.size()==2) {
	for(long ii = 11; ii < 17; ii += 2) {
	  addToList(-37, 1000001 + ii, -1000000 - ii);
	  addToList(-37, 1000001 + ii, -2000000 - ii);
	}
      }
    }
  }
  SSSVertex::doinit();
  // extract the sfermion mixing matrices
  // sfermion mixing
  stop_ = model->stopMix();
  sbottom_ = model->sbottomMix();
  if(!stop_ || !sbottom_)
    throw InitException() << "RPVSSSVertex::doinit() - "
			  << "A mixing matrix pointer is null."
			  << " stop: " << stop_ << " sbottom: " << sbottom_
			  << Exception::abortnow;
  stau_ = model->stauMix();
  if(!stau_ && (!mixC || mixC->size().first<2))
    throw InitException() << "RPVSSSVertex::doinit() either the stau"
			  << " mixing matrix must be set or the stau"
			  << " included in mixing with the"
			  << " charged Higgs bosons" << Exception::abortnow;
  // couplings of the neutral Higgs bosons 
  Energy mw = getParticleData(ParticleID::Wplus)->mass();
  double sw = sqrt(sin2ThetaW());
  double cw = sqrt(1.-sin2ThetaW());
  // extract the vevs
  vector<Energy> vnu = model->sneutrinoVEVs();
  double g = electroMagneticCoupling(sqr(mw))/sw;
  Energy v = 2.*mw/g;
  double tanb = model->tanBeta();
  vector<Energy> vevs(5);
  vevs[0] = sqrt((sqr(v)-sqr(vnu[0])-sqr(vnu[1])-sqr(vnu[2]))/
		   (1.+sqr(tanb)));
  vevs[1] = vevs[0]*tanb;
  for(unsigned int ix=0;ix<vnu.size();++ix) vevs[2+ix] = vnu[ix];
  // bilinear RPV terms
  // coupling of three scalar Higgs bosons
  scalarScalarScalar_.resize(scalar.size(),vector<vector<complex<Energy > > >
			     (scalar.size(),vector<complex<Energy> >(scalar.size(),complex<Energy>(ZERO))));
  double delta[5]={1.,-1.,1.,1.,1.};
  if(vevs.size()>scalar.size()) vevs.resize(scalar.size());
  for(unsigned int i1=0;i1<scalar.size();++i1) {
    for(unsigned int i2=0;i2<scalar.size();++i2) {
      for(unsigned int i3=0;i3<scalar.size();++i3) {
	scalarScalarScalar_[i1][i2][i3] = ZERO;
	for(unsigned int ix=0;ix<vevs.size();++ix) {
	  for(unsigned int k=0;k<scalar.size();++k) {
	    scalarScalarScalar_[i1][i2][i3] += 
	      delta[ix]*vevs[ix]*(*mixH)(i1,ix)*delta[k]*(*mixH)(i2,k)*(*mixH)(i3,k)+
	      delta[ix]*vevs[ix]*(*mixH)(i2,ix)*delta[k]*(*mixH)(i1,k)*(*mixH)(i3,k)+
	      delta[ix]*vevs[ix]*(*mixH)(i3,ix)*delta[k]*(*mixH)(i1,k)*(*mixH)(i2,k);
	  }
	}
	scalarScalarScalar_[i1][i2][i3] *= -0.25*g/sqr(cw);
      }
    }
  }
  // coupling of scalar Higgs to 2 pseudoscalars 
  scalarPseudoPseudo_.resize(scalar.size(),vector<vector<complex<Energy > > >
			     (pseudo.size(),vector<complex<Energy> >(pseudo.size(),complex<Energy>(ZERO))));
  for(unsigned int i1=0;i1<scalar.size();++i1) {
    for(unsigned int i2=0;i2<pseudo.size();++i2) {
      for(unsigned int i3=0;i3<pseudo.size();++i3) {
	scalarPseudoPseudo_[i1][i2][i3] = ZERO;
	for(unsigned int ix=0;ix<vevs.size();++ix) {
	  for(unsigned int k=0;k<pseudo.size();++k) {
	    scalarPseudoPseudo_[i1][i2][i3] += 
	      delta[k]*delta[ix]*vevs[ix]*(*mixH)(i1,ix)*(*mixP)(i2,k)*(*mixP)(i3,k);
	  }
	}
	scalarPseudoPseudo_[i1][i2][i3] *= -0.25*g/sqr(cw);
      }
    }
  }
  // coupling of scalar Higgs to 2 charged Higgs
  double gp = g*sw/cw;
  Energy mu = model->muParameter();
  vector<Energy> eps = model->epsilon();
  scalarChargedCharged_.resize(scalar.size(),vector<vector<complex<Energy > > >
			       (charged.size(),vector<complex<Energy> >(charged.size(),complex<Energy>(ZERO))));
  double yl[3] = {sqrt(2.)*getParticleData(11)->mass()/vevs[1],
		  sqrt(2.)*getParticleData(13)->mass()/vevs[1],
		  sqrt(2.)*getParticleData(15)->mass()/vevs[1]};
  complex<Energy> Al[3] =  {complex<Energy>(ZERO),complex<Energy>(ZERO),yl[2]*model->tauTrilinear()};
  for(unsigned int i1=0;i1<scalar.size();++i1) {
    for(unsigned int i2=0;i2<charged.size();++i2) {
      for(unsigned int i3=0;i3<charged.size();++i3) {
	complex<Energy> umdi = vevs[0]*(*mixH)(i1,0)-vevs[1]*(*mixH)(i1,1)+vevs[2]*(*mixH)(i1,2)+
	  vevs[3]*(*mixH)(i1,3)+vevs[4]*(*mixH)(i1,4);
	scalarChargedCharged_[i1][i2][i3] =
	  // HH term
	  0.25*sqr(g)*(-vevs[1]*((*mixH)(i1,1)*((*mixC)(i2,1)*(*mixC)(i3,1)+(*mixC)(i2,2)*(*mixC)(i3,2))+
				 (*mixH)(i1,0)*((*mixC)(i2,2)*(*mixC)(i3,1)+(*mixC)(i2,1)*(*mixC)(i3,2)))
		       -vevs[0]*((*mixH)(i1,0)*((*mixC)(i2,1)*(*mixC)(i3,1)+(*mixC)(i2,2)*(*mixC)(i3,2))+
				 (*mixH)(i1,1)*((*mixC)(i2,2)*(*mixC)(i3,1)+(*mixC)(i2,1)*(*mixC)(i3,2)))
		       +(vevs[2]*(*mixH)(i1,2)+vevs[3]*(*mixH)(i1,3)+vevs[4]*(*mixH)(i1,4))*
		       ((*mixC)(i2,1)*(*mixC)(i3,1)-(*mixC)(i2,2)*(*mixC)(i3,2)))-
	  0.25*sqr(gp)*((*mixC)(i2,1)*(*mixC)(i3,1)-(*mixC)(i2,2)*(*mixC)(i3,2))*umdi
	  -0.5*(*mixC)(i2,1)*(*mixC)(i3,1)*(yl[0]*vevs[2]*(*mixH)(i1,2)+yl[1]*vevs[3]*(*mixH)(i1,3)+yl[2]*vevs[4]*(*mixH)(i1,4))
	  // LL term
	  +0.25*(sqr(g)-sqr(gp))*umdi*((*mixC)(i2,2)*(*mixC)(i3,2)+(*mixC)(i2,3)*(*mixC)(i3,3)+(*mixC)(i2,4)*(*mixC)(i3,4))
	  -vevs[0]*(*mixH)(i1,1)*(yl[0]*(*mixC)(i2,2)*(*mixC)(i3,2)+yl[1]*(*mixC)(i2,3)*(*mixC)(i3,3)+yl[2]*(*mixC)(i2,4)*(*mixC)(i3,4))
	  -0.25*sqr(g)*((vevs[2]*(*mixC)(i3,2)+vevs[3]*(*mixC)(i3,3)+vevs[4]*(*mixC)(i3,4))*
			((*mixH)(i1,2)*(*mixC)(i2,2)+(*mixH)(i1,3)*(*mixC)(i2,3)+(*mixH)(i1,4)*(*mixC)(i2,4))+
			(vevs[2]*(*mixC)(i2,2)+vevs[3]*(*mixC)(i2,3)+vevs[4]*(*mixC)(i2,4))*
			((*mixH)(i1,2)*(*mixC)(i3,2)+(*mixH)(i1,3)*(*mixC)(i3,3)+(*mixH)(i1,4)*(*mixC)(i3,4)))
	  // RR term
	  +0.5*sqr(gp)*umdi*((*mixC)(i2,5)*(*mixC)(i3,5)+(*mixC)(i2,6)*(*mixC)(i3,6)+(*mixC)(i2,7)*(*mixC)(i3,7))
	  -vevs[0]*(*mixH)(i1,1)*(yl[0]*(*mixC)(i2,5)*(*mixC)(i3,5)+yl[1]*(*mixC)(i2,6)*(*mixC)(i3,6)+yl[2]*(*mixC)(i2,7)*(*mixC)(i3,7))
	  -0.5*((vevs[2]*yl[0]*(*mixC)(i2,5)+vevs[3]*yl[1]*(*mixC)(i2,6)+vevs[4]*yl[2]*(*mixC)(i2,7))*
		(yl[0]*(*mixH)(i1,2)*(*mixC)(i3,5)+yl[1]*(*mixH)(i1,3)*(*mixC)(i3,6)+yl[2]*(*mixH)(i1,4)*(*mixC)(i3,7))+
		(vevs[2]*yl[0]*(*mixC)(i3,5)+vevs[3]*yl[1]*(*mixC)(i3,6)+vevs[4]*yl[2]*(*mixC)(i3,7))*
		(yl[0]*(*mixH)(i1,2)*(*mixC)(i2,5)+yl[1]*(*mixH)(i1,3)*(*mixC)(i2,6)+yl[2]*(*mixH)(i1,4)*(*mixC)(i2,7)))
	  // LR
	  -(*mixH)(i1,1)/sqrt(2.)*(Al[0]*(*mixC)(i2,2)*(*mixC)(i3,2)+Al[1]*(*mixC)(i2,3)*(*mixC)(i3,3)+Al[2]*(*mixC)(i2,4)*(*mixC)(i3,4))
	  +(*mixH)(i1,2)/sqrt(2.)*mu*(yl[0]*(*mixC)(i2,2)*(*mixC)(i3,2)+yl[1]*(*mixC)(i2,3)*(*mixC)(i3,3)+yl[2]*(*mixC)(i2,4)*(*mixC)(i3,4))
	  // RL
	  -(*mixH)(i1,1)/sqrt(2.)*(Al[0]*(*mixC)(i3,2)*(*mixC)(i2,2)+Al[1]*(*mixC)(i3,3)*(*mixC)(i2,3)+Al[2]*(*mixC)(i3,4)*(*mixC)(i2,4))
	  +(*mixH)(i1,2)/sqrt(2.)*mu*(yl[0]*(*mixC)(i3,2)*(*mixC)(i2,2)+yl[1]*(*mixC)(i3,3)*(*mixC)(i2,3)+yl[2]*(*mixC)(i3,4)*(*mixC)(i2,4))
	  // HL
	  -0.25*sqr(g)*(((*mixH)(i1,2)*(*mixC)(i3,2)+(*mixH)(i1,3)*(*mixC)(i3,3)+(*mixH)(i1,4)*(*mixC)(i3,4))*
			(vevs[0]*(*mixC)(i2,0)+vevs[1]*(*mixC)(i2,1))
			+(vevs[2]*(*mixC)(i3,2)+vevs[3]*(*mixC)(i3,3)+vevs[4]*(*mixC)(i3,4))*
			((*mixH)(i1,0)*(*mixC)(i2,0)+(*mixH)(i1,1)*(*mixC)(i2,1)))
	  +0.5*(*mixH)(i1,0)*(*mixC)(i2,0)*
	  (sqr(yl[0])*vevs[2]*(*mixC)(i3,2)+sqr(yl[1])*vevs[3]*(*mixC)(i3,3)+sqr(yl[2])*vevs[4]*(*mixC)(i3,4))
	  +0.5*vevs[0]*(*mixC)(i2,0)*(sqr(yl[0])*(*mixH)(i1,2)*(*mixC)(i3,2)+
				      sqr(yl[1])*(*mixH)(i1,3)*(*mixC)(i3,3)+
				      sqr(yl[2])*(*mixH)(i1,4)*(*mixC)(i3,4))
	  // LH
	  -0.25*sqr(g)*(((*mixH)(i1,2)*(*mixC)(i2,2)+(*mixH)(i1,3)*(*mixC)(i2,3)+(*mixH)(i1,4)*(*mixC)(i2,4))*
			(vevs[0]*(*mixC)(i3,0)+vevs[1]*(*mixC)(i3,1))
			+(vevs[2]*(*mixC)(i2,2)+vevs[3]*(*mixC)(i2,3)+vevs[4]*(*mixC)(i2,4))*
			((*mixH)(i1,0)*(*mixC)(i3,0)+(*mixH)(i1,1)*(*mixC)(i3,1)))
	  +0.5*(*mixH)(i1,0)*(*mixC)(i3,0)*
	  (sqr(yl[0])*vevs[2]*(*mixC)(i2,2)+sqr(yl[1])*vevs[3]*(*mixC)(i2,3)+sqr(yl[2])*vevs[4]*(*mixC)(i2,4))
	  +0.5*vevs[0]*(*mixC)(i3,0)*(sqr(yl[0])*(*mixH)(i1,2)*(*mixC)(i2,2)+
				      sqr(yl[1])*(*mixH)(i1,3)*(*mixC)(i2,3)+
				      sqr(yl[2])*(*mixH)(i1,4)*(*mixC)(i2,4))
	  // HR
	  +sqrt(0.5)*(eps[0]*yl[0]*(*mixC)(i3,5)+eps[1]*yl[1]*(*mixC)(i3,6)+eps[2]*yl[2]*(*mixC)(i3,7))*
	  ((*mixH)(i1,0)*(*mixC)(i2,1)+(*mixH)(i1,1)*(*mixC)(i2,0))
	  +sqrt(0.5)*(*mixC)(i2,0)*(Al[0]*(*mixH)(i1,2)*(*mixC)(i3,5)+
				    Al[1]*(*mixH)(i1,2)*(*mixC)(i3,6)+Al[2]*(*mixH)(i1,2)*(*mixC)(i3,7))
	  +sqrt(0.5)*mu*(*mixC)(i2,1)*(yl[0]*(*mixH)(i1,2)*(*mixC)(i3,5)+
				       yl[1]*(*mixH)(i1,3)*(*mixC)(i3,6)+
				       yl[2]*(*mixH)(i1,4)*(*mixC)(i3,7))
	  // RH
	  +sqrt(0.5)*(eps[0]*yl[0]*(*mixC)(i2,5)+eps[1]*yl[1]*(*mixC)(i2,6)+eps[2]*yl[2]*(*mixC)(i2,7))*
	  ((*mixH)(i1,0)*(*mixC)(i3,1)+(*mixH)(i1,1)*(*mixC)(i3,0))
	  +sqrt(0.5)*(*mixC)(i3,0)*(Al[0]*(*mixH)(i1,2)*(*mixC)(i2,5)+
				    Al[1]*(*mixH)(i1,2)*(*mixC)(i2,6)+Al[2]*(*mixH)(i1,2)*(*mixC)(i2,7))
	  +sqrt(0.5)*mu*(*mixC)(i3,1)*(yl[0]*(*mixH)(i1,2)*(*mixC)(i2,5)+
				       yl[1]*(*mixH)(i1,3)*(*mixC)(i2,6)+
				       yl[2]*(*mixH)(i1,4)*(*mixC)(i2,7));
	// normalization
	scalarChargedCharged_[i1][i2][i3] /= g;
      }
    }
  }
  // coupling of the pseudoscalar Higgs to two charged Higgs bosons
  pseudoChargedCharged_.resize(pseudo.size(),vector<vector<complex<Energy > > >
  			       (charged.size(),vector<complex<Energy> >(charged.size(),complex<Energy>(ZERO))));
  if(charged.size()<5) {
    for(unsigned int i1=0;i1<pseudo.size();++i1) {
      for(unsigned int i2=0;i2<charged.size();++i2) {
	for(unsigned int i3=0;i3<charged.size();++i3) {
	  pseudoChargedCharged_[i1][i2][i3] = ZERO;
	  for(unsigned int il=0;il<3;++il) {
	    pseudoChargedCharged_[i1][i2][i3] += (Al[il]*(*mixP)(i1,0)-mu*yl[il]*(*mixP)(i1,1))*
	      ((*mixC)(i2,il+2)*(*mixC)(i3,il+5)-(*mixC)(i2,il+5)*(*mixC)(i3,il+2));
	  }
	  // normalization // do not revert to *= , breaks with XCode 5.1
	  pseudoChargedCharged_[i1][i2][i3] = pseudoChargedCharged_[i1][i2][i3] * Complex(0.,-1.)*sqrt(0.5)/g;
	}
      }
    }
  }
  // couplings of the Higgs bosons to the up type squarks
  scalarSup_.resize(scalar.size(),vector<vector<vector<complex<Energy > > > >
		    (3,vector<vector<complex<Energy> > >(2,vector<complex<Energy> >(2,complex<Energy>(ZERO)))));
  for(unsigned int iq=0;iq<3;++iq) {
    double y = sqrt(2.)*getParticleData(2*long(iq+1))->mass()/vevs[1];
    complex<Energy> A = iq!=2 ? complex<Energy>(ZERO) : y*model->topTrilinear();
    Complex mixing[2][2];
    if(iq!=2) {
      mixing[0][0] = 1.;
      mixing[0][1] = 0.;
      mixing[1][0] = 0.;
      mixing[1][1] = 1.;
    }
    else {
      mixing[0][0] = (*stop_)(0,0);
      mixing[0][1] = (*stop_)(0,1);
      mixing[1][0] = (*stop_)(1,0);
      mixing[1][1] = (*stop_)(1,1);
    }
    // loop over Higgs
    for(unsigned int ix=0;ix<2;++ix) {
      for(unsigned int iy=0;iy<2;++iy) {
	for(unsigned int ih=0;ih<scalar.size();++ih) {
	  // first LL term
	  Complex LL1 =-(sqr(g)/4.- sqr(gp)/12.)*mixing[ix][0]*conj(mixing[iy][0]);
	  // first RR term 
	  Complex RR1 = -1./3.*sqr(gp)*mixing[ix][1]*conj(mixing[iy][1]);
	  complex<Energy> factors[5];
	  for(unsigned int iz=0;iz<5;++iz) factors[iz] = delta[iz]*vevs[iz]*(LL1+RR1);
	  factors[0] += y/sqrt(2.)*mu*( mixing[ix][0]*conj(mixing[iy][1]) + 
					mixing[ix][1]*conj(mixing[iy][0]) );
	  factors[1] += -vevs[1]*sqr(y)*( mixing[ix][0]*conj(mixing[iy][0]) + 
					  mixing[ix][1]*conj(mixing[iy][1]) )
	    -A/sqrt(2.)*(mixing[ix][0]*conj(mixing[iy][1]) + 
			 mixing[ix][1]*conj(mixing[iy][0]));
	  for(unsigned int iz=0;iz<3;++iz) {
	    factors[iz+2] += -y*eps[iz]/sqrt(2.)*( mixing[ix][0]*conj(mixing[iy][1]) + 
						    mixing[ix][1]*conj(mixing[iy][0]) );
	  }
	  for(unsigned int iz=0;iz<scalar.size();++iz) {
	    scalarSup_[ih][iq][ix][iy] += (*mixH)(ih,iz)*factors[iz];
	  }
	  scalarSup_[ih][iq][ix][iy] /= g;
	}
      }
    }
  }
  // couplings of the pseudoscalar Higgs bosons to the up type squarks
  pseudoSup_.resize(pseudo.size(),vector<complex<Energy > >(3,complex<Energy>(ZERO)));
  for(unsigned int iq=0;iq<3;++iq) {
    double y = sqrt(2.)*getParticleData(2*long(iq+1))->mass()/vevs[1];
    complex<Energy> A = iq!=2 ? complex<Energy>(ZERO) : y*model->topTrilinear();
    // loop over Higgs
    for(unsigned int ih=0;ih<pseudo.size();++ih) {
      pseudoSup_[ih][iq] = Complex(0.,-1.)*sqrt(0.5)*(A*(*mixP)(ih,1)-y*mu*(*mixP)(ih,0))/g;
    }
  }
  // couplings of the Higgs bosons to the down type squarks
  scalarSdown_.resize(scalar.size(),vector<vector<vector<complex<Energy > > > >
		      (3,vector<vector<complex<Energy> > >(2,vector<complex<Energy> >(2,complex<Energy>(ZERO)))));
  for(unsigned int iq=0;iq<3;++iq) {
    double y = sqrt(2.)*getParticleData(2*long(iq+1)-1)->mass()/vevs[0];
    Complex mixing[2][2];
    complex<Energy> A = iq!=2 ? complex<Energy>(ZERO) : y*model->bottomTrilinear();
    if(iq!=2) {
      mixing[0][0] = 1.;
      mixing[0][1] = 0.;
      mixing[1][0] = 0.;
      mixing[1][1] = 1.;
    }
    else {
      mixing[0][0] = (*sbottom_)(0,0);
      mixing[0][1] = (*sbottom_)(0,1);
      mixing[1][0] = (*sbottom_)(1,0);
      mixing[1][1] = (*sbottom_)(1,1);
    }
    for(unsigned int ix=0;ix<2;++ix) {
      for(unsigned int iy=0;iy<2;++iy) {
	for(unsigned int ih=0;ih<scalar.size();++ih) {
	  // first LL term
	  Complex LL1 = (sqr(g)/4.+sqr(gp)/12.)*mixing[ix][0]*conj(mixing[iy][0]);
	  // first RR term 
	  Complex RR1 = 1./6.*sqr(gp)*mixing[ix][1]*conj(mixing[iy][1]);
	  complex<Energy> factors[5];
	  for(unsigned int iz=0;iz<5;++iz) factors[iz] = delta[iz]*vevs[iz]*(LL1+RR1);
	  factors[0] += -vevs[0]*sqr(y)*(mixing[ix][0]*conj(mixing[iy][0])+
					 mixing[ix][1]*conj(mixing[iy][1]))
	    - A/sqrt(2.)*(mixing[ix][0]*conj(mixing[iy][1])+
			  mixing[ix][1]*conj(mixing[iy][0]));
	  factors[1] +=  y/sqrt(2.)*mu*(mixing[ix][0]*conj(mixing[iy][1])+
					mixing[ix][1]*conj(mixing[iy][0]));
	  for(unsigned int iz=0;iz<scalar.size();++iz) {
	    scalarSdown_[ih][iq][ix][iy] += (*mixH)(ih,iz)*factors[iz];
	  }
	  scalarSdown_[ih][iq][ix][iy] /= g;
	}
      }
    }
  }
  // couplings of the pseudoscalar Higgs bosons to the down type squarks
  pseudoSdown_.resize(pseudo.size(),vector<complex<Energy > >(3,complex<Energy>(ZERO)));
  for(unsigned int iq=0;iq<3;++iq) {
    double y = sqrt(2.)*getParticleData(2*long(iq+1)-1)->mass()/vevs[0];
    complex<Energy> A = iq!=2 ? complex<Energy>(ZERO) : y*model->bottomTrilinear();
    for(unsigned int ih=0;ih<pseudo.size();++ih) {
      pseudoSdown_[ih][iq] = Complex(0.,-1.)*sqrt(0.5)*(A*(*mixP)(ih,0)-mu*y*(*mixP)(ih,1))/g;
    }
  }
  // couplings of the scalar Higgs bosons to the sneutrinos
  if(scalar.size()==2) {
    scalarSneutrino_.resize(scalar.size(),vector<complex<Energy > >(3,complex<Energy>(ZERO)));
    for(unsigned int il=0;il<3;++il) {
      // loop over Higgs
      for(unsigned int ih=0;ih<scalar.size();++ih) {
	// first LL term
	Complex LL1 =-(sqr(g)/4.+sqr(gp)/4.);
	complex<Energy> factors[5];
	for(unsigned int iz=0;iz<5;++iz) factors[iz] = delta[iz]*vevs[iz]*LL1;
	for(unsigned int iz=0;iz<scalar.size();++iz) {
	  scalarSneutrino_[ih][il] += (*mixH)(ih,iz)*factors[iz];
	}
	scalarSneutrino_[ih][il] /= g;
      }
    }
  }
  // couplings of the Higgs bosons to the charged sleptons
  if(charged.size()==1) {
    scalarSlepton_.resize(scalar.size(),vector<vector<vector<complex<Energy > > > >
			  (3,vector<vector<complex<Energy> > >(2,vector<complex<Energy> >(2,complex<Energy>(ZERO)))));
    for(unsigned int il=0;il<3;++il) {
      double y = sqrt(2.)*getParticleData(2*long(il+6)-1)->mass()/vevs[0];
      Complex mixing[2][2];
      complex<Energy> A = il!=2 ? complex<Energy>(ZERO) : y*model->tauTrilinear();
      if(il!=2) {
	mixing[0][0] = 1.;
	mixing[0][1] = 0.;
	mixing[1][0] = 0.;
	mixing[1][1] = 1.;
      }
      else {
	mixing[0][0] = (*stau_)(0,0);
	mixing[0][1] = (*stau_)(0,1);
	mixing[1][0] = (*stau_)(1,0);
	mixing[1][1] = (*stau_)(1,1);
      }
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  for(unsigned int ih=0;ih<scalar.size();++ih) {
	    // first LL term
	    Complex LL1 = (sqr(g)/4.-sqr(gp)/4.)*mixing[ix][0]*conj(mixing[iy][0]);
	    // first RR term 
	    Complex RR1 = 0.5*sqr(gp)*mixing[ix][1]*conj(mixing[iy][1]);
	    complex<Energy> factors[5];
	    for(unsigned int iz=0;iz<5;++iz) factors[iz] = delta[iz]*vevs[iz]*(LL1+RR1);
	    factors[0] += -vevs[0]*sqr(y)*(mixing[ix][0]*conj(mixing[iy][0])+
					   mixing[ix][1]*conj(mixing[iy][1]))
	      - A/sqrt(2.)*(mixing[ix][0]*conj(mixing[iy][1])+
			    mixing[ix][1]*conj(mixing[iy][0]));
	    factors[1] +=  y/sqrt(2.)*mu*(mixing[ix][0]*conj(mixing[iy][1])+
					  mixing[ix][1]*conj(mixing[iy][0]));
	    for(unsigned int iz=0;iz<scalar.size();++iz) {
	      scalarSlepton_[ih][il][ix][iy] += (*mixH)(ih,iz)*factors[iz];
	    }
	    scalarSlepton_[ih][il][ix][iy] /= g;
	  }
	}
      }
    }
  }
  // couplings of the pseudoscalar Higgs bosons to the charged sleptons
  if(charged.size()==1) {
    pseudoSlepton_.resize(pseudo.size(),vector<complex<Energy> >(3,complex<Energy>(ZERO)));
    for(unsigned int il=0;il<3;++il) {
      double y = sqrt(2.)*getParticleData(2*long(il+6)-1)->mass()/vevs[0];
      complex<Energy> A = il!=2 ? complex<Energy>(ZERO) : y*model->tauTrilinear();
      for(unsigned int ih=0;ih<pseudo.size();++ih) {
	pseudoSlepton_[ih][il] = Complex(0.,-1.)*sqrt(0.5)*(A*(*mixP)(ih,0)-mu*y*(*mixP)(ih,1))/g;
      }
    }
  }
  // charged Higgs squarks
  chargedSquark_.resize(charged.size(),vector<vector<vector<complex<Energy > > > >
			(3,vector<vector<complex<Energy> > >(2,vector<complex<Energy> >(2,complex<Energy>(ZERO)))));
  for(unsigned int iq=0;iq<3;++iq) {
    double yd = sqrt(2.)*getParticleData(2*long(iq+1)-1)->mass()/vevs[0];
    double yu = sqrt(2.)*getParticleData(2*long(iq+1))->mass()/vevs[1];
    Complex mixd[2][2],mixu[2][2];
    complex<Energy> Ad = iq!=2 ? complex<Energy>(ZERO) : yd*model->bottomTrilinear();
    complex<Energy> Au = iq!=2 ? complex<Energy>(ZERO) : yu*model->topTrilinear();
    if(iq!=2) {
      mixd[0][0] = mixu[0][0] = 1.;
      mixd[0][1] = mixu[0][1] = 0.;
      mixd[1][0] = mixu[1][0] = 0.;
      mixd[1][1] = mixu[1][1] = 1.;
    }
    else {
      mixd[0][0] = (*sbottom_)(0,0);
      mixd[0][1] = (*sbottom_)(0,1);
      mixd[1][0] = (*sbottom_)(1,0);
      mixd[1][1] = (*sbottom_)(1,1);
      mixu[0][0] = (*stop_)(0,0);
      mixu[0][1] = (*stop_)(0,1);
      mixu[1][0] = (*stop_)(1,0);
      mixu[1][1] = (*stop_)(1,1);
    }
    for(unsigned int ih=0;ih<charged.size();++ih) {
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  // various terms (up first) (down second)
	  complex<Energy> LL(ZERO);
	  for(unsigned int iz=0;iz<vevs.size();++iz)
	    LL += 0.5*sqr(g)*vevs[iz]*(*mixC)(ih,iz);
	  LL += -vevs[0]*(*mixC)(ih,0)*sqr(yd)-vevs[1]*(*mixC)(ih,1)*sqr(yu);
	  LL *= -sqrt(0.5);
	  complex<Energy> RR = 
	    yu*yd*(vevs[1]*(*mixC)(ih,0) + vevs[0]*(*mixC)(ih,1))/sqrt(2.);
	  complex<Energy> LR = (*mixC)(ih,0)*Ad+(*mixC)(ih,1)*yd*mu;
	  complex<Energy> RL = yu*mu*(*mixC)(ih,0) + Au*(*mixC)(ih,1);
	  chargedSquark_[ih][iq][ix][iy] = 
	    mixu[ix][0]*mixd[iy][0]*LL + mixu[ix][1]*mixd[iy][1]*RR +
	    mixu[ix][0]*mixd[iy][1]*LR + mixu[ix][1]*mixd[iy][0]*RL;
	  chargedSquark_[ih][iq][ix][iy] /=g;
	}
      }
    }
  }
  // charged Higgs slepton
  if(charged.size()==1) {
    chargedSlepton_.resize(charged.size(),vector<vector<complex<Energy > > >
			   (3,vector<complex<Energy> >(2,complex<Energy>(ZERO))));
    for(unsigned int il=0;il<3;++il) {
      double y = sqrt(2.)*getParticleData(2*long(il+6)-1)->mass()/vevs[0];
      Complex mixd[2][2];
      complex<Energy> Al = il!=2 ? complex<Energy>(ZERO) : y*model->tauTrilinear();
      if(il!=2) {
	mixd[0][0] = 1.;
	mixd[0][1] = 0.;
	mixd[1][0] = 0.;
	mixd[1][1] = 1.;
      }
      else {
	mixd[0][0] = (*stau_)(0,0);
	mixd[0][1] = (*stau_)(0,1);
	mixd[1][0] = (*stau_)(1,0);
	mixd[1][1] = (*stau_)(1,1);
      }
      for(unsigned int ih=0;ih<charged.size();++ih) {
	for(unsigned int ix=0;ix<2;++ix) {
	  // various terms (charged lepton second)
	  complex<Energy> LL(ZERO);
	  for(unsigned int iz=0;iz<vevs.size();++iz)
	    LL += 0.5*sqr(g)*vevs[iz]*(*mixC)(ih,iz);
	  LL += -vevs[0]*(*mixC)(ih,0)*sqr(y);
	  LL *= -sqrt(0.5);
	  complex<Energy> LR = (*mixC)(ih,0)*Al+(*mixC)(ih,1)*y*mu;
	  chargedSlepton_[ih][il][ix] = mixd[ix][0]*LL + mixd[ix][1]*LR;
	  chargedSlepton_[ih][il][ix] /=g;
	}
      }
    }
  }
}

namespace {

int scalarEigenState(long id) {
  if(id<1000000)
    return (id-25)/10;
  else
    return (id-1000008)/2;
}

int pseudoEigenState(long id) {
  if(id<1000000)
    return 0;
  else
    return (id-1000016);
}

int chargedEigenState(long id) {
  if(abs(id)<1000000)
    return 0;
  else if(abs(id)<2000000)
    return (abs(id)-1000009)/2;
  else
    return (abs(id)-2000003)/2;
}

}

void RPVSSSVertex::setCoupling(Energy2 q2, tcPDPtr part1,
			      tcPDPtr part2,tcPDPtr part3) {
  // prefactor
  if( q2 != q2Last_ || gLast_==0. ) {
    q2Last_ = q2;
    gLast_ = weakCoupling(q2);
  }
  // some ordering of the particles
  // sort in order of PDG codes, smallest first
  if(abs(part1->id())>abs(part2->id())) swap(part1,part2);
  if(abs(part1->id())>abs(part3->id())) swap(part1,part3);
  if(abs(part2->id())>abs(part3->id())) swap(part2,part3);
  // make sure squarks 2nd and 3rd
  if(abs(part1->id())%1000000<=6) swap(part1,part2);
  if(abs(part1->id())%1000000<=6) swap(part1,part3);
  // extract particle ids
  long sca1(part1->id()), sca2(part2->id()), sca3(part3->id());
  // Higgs squark couplings
  if(abs(sca2)%1000000<=6) {
    // charged Higgs
    if( part1->charged() ) {
      int ih = chargedEigenState(sca1);
      if(abs(sca2)%2!=0) swap(sca2,sca3);
      unsigned int alpha = abs(sca2)/1000000-1;
      unsigned int beta  = abs(sca3)/1000000-1;
      unsigned int iq = (abs(sca2)%1000000-2)/2;
      norm(gLast_*chargedSquark_[ih][iq][alpha][beta]*UnitRemoval::InvE);
      return;
    }
    // neutral Higgs
    else {
      long sm(0);
      // get the left/right light/heavy state
      unsigned int alpha(abs(sca2)/1000000 - 1), beta(abs(sca3)/1000000 - 1);
      sm = abs(sca2)%1000000;
      bool pseudo = sca1==36 || (sca1>=1000017&&sca1<=1000019);
      int higgs = pseudo ? pseudoEigenState(sca1) : scalarEigenState(sca1);
      complex<Energy> coup;
      if(!pseudo) {
	if( sm % 2 == 0 ) {
	  coup = scalarSup_  [higgs][(sm-2)/2][alpha][beta];
	}
	else {
	  coup = scalarSdown_[higgs][(sm-1)/2][alpha][beta];
	}
      }
      else {
	if( sm % 2 == 0 ) {
	  coup = pseudoSup_  [higgs][(sm-2)/2];
	}
	else {
	  coup = pseudoSdown_[higgs][(sm-1)/2];
	}
	if((alpha==1&&sca2<0)||(beta==1&&sca3<0)) coup *=-1.;
      }
      // set coupling and return
      norm(gLast_*coup*UnitRemoval::InvE);
      return;
    }
  }
  // all neutral
  else if(!part1->charged()&&!part2->charged()&&!part3->charged()) {
    // neutral Higgs sneutrino
    if(scalarScalarScalar_.size()==2&&abs(sca2)>1000000) {
      assert(!(sca1==36 || (sca1>=1000017&&sca1<=1000019)));
      int il = (abs(sca2)-1000012)/2;
      norm(gLast_*scalarSneutrino_[scalarEigenState(sca1)][il]*UnitRemoval::InvE);
      return;
    }
    else {
      if(sca1==36 || (sca1>=1000017&&sca2<=1000019)) swap(sca1,sca2);
      if(sca1==36 || (sca1>=1000017&&sca2<=1000019)) swap(sca1,sca3);
      // 2 pseudoscalar 1 scale
      if(sca2==36 || (sca2>=1000017&&sca2<=1000019)) {
	norm(gLast_*scalarPseudoPseudo_[scalarEigenState(sca1)]
	     [pseudoEigenState(sca2)][pseudoEigenState(sca3)]*UnitRemoval::InvE);
	return;
      }
      // 3 scalars
      else {
	norm(gLast_*scalarScalarScalar_[scalarEigenState(sca1)]
	     [scalarEigenState(sca2)][scalarEigenState(sca3)]*UnitRemoval::InvE);
	return;
      }
    }
  }
  // two charged
  else {
    // put the charged last
    if(!part2->charged()) {
      swap(part1,part2);
      swap(sca1 ,sca2 );
    }
    if(!part3->charged()) {
      swap(part1,part3);
      swap(sca1 ,sca3 );
    }
    // sleptons
    if(scalarChargedCharged_[0].size()<5&&abs(sca2)>1000000) {
      // neutral Higgs charged sleptons
      if((abs(sca2)>=1000011&&abs(sca2)<=1000015&&abs(sca2)%2!=0)&&
	 (abs(sca3)>=1000011&&abs(sca3)<=1000015&&abs(sca3)%2!=0)) {
	long sm(0);
	// get the left/right light/heavy state
	unsigned int alpha(abs(sca2)/1000000 - 1), beta(abs(sca3)/1000000 - 1);
	sm = (abs(sca2)%1000000-11)/2;
	bool pseudo = sca1==36 || (sca1>=1000017&&sca1<=1000019);
	int higgs = pseudo ? pseudoEigenState(sca1) : scalarEigenState(sca1);
	complex<Energy> coup;
	if(!pseudo) {
	  coup = scalarSlepton_[higgs][(sm-1)/2][alpha][beta];
	}
	else {
    	  coup = pseudoSlepton_[higgs][(sm-1)/2];
	  if((alpha==1&&sca2<0)||(beta==1&&sca3<0)) coup *=-1.;
    	}
	// set coupling and return
	norm(gLast_*coup*UnitRemoval::InvE);
	return;
      }
      // charged Higgs
      else {
	if(abs(sca2)<1000000) {
	  swap(part1,part2);
	  swap(sca1 ,sca2 );
	}
	if(abs(sca3)<1000000) {
	  swap(part1,part3);
	  swap(sca1 ,sca3 );
	}
	if(part3->charged()) {
	  swap(part2,part3);
	  swap(sca2 ,sca3 );
	}
	unsigned int il    = (abs(sca2)%1000000-11)/2;
	unsigned int alpha = abs(sca2)/1000000-1;
	norm(gLast_*chargedSlepton_[chargedEigenState(sca1)][il][alpha]*UnitRemoval::InvE);
      }
    }
    else {
      // pseudoscalar charged charged
      if(sca1==36 || (sca1>=1000017&&sca1<=1000019)) {
	if(sca2>0) swap(sca2,sca3);
	norm(gLast_*pseudoChargedCharged_[pseudoEigenState(sca1)]
	     [chargedEigenState(sca2)][chargedEigenState(sca3)]*UnitRemoval::InvE);
      }
      // scalar charged charged
      else {
	norm(gLast_*scalarChargedCharged_[scalarEigenState(sca1)]
	     [chargedEigenState(sca2)][chargedEigenState(sca3)]*UnitRemoval::InvE);
      }
    }
    return;
  }	
}
