
#include "VBFNLOCommonBlocks.h"
#include "VBFNLOMEqq2hqqg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/Throw.h"
 
#include "Diagrams/DC_Hjjj_CC_qBqB.h"
#include "Diagrams/DC_Hjjj_CC_qBq.h"
#include "Diagrams/DC_Hjjj_CC_qqB.h"
#include "Diagrams/DC_Hjjj_CC_qq.h"
#include "Diagrams/DC_Hjjj_NC_qBqB.h"
#include "Diagrams/DC_Hjjj_NC_qBq.h"
#include "Diagrams/DC_Hjjj_NC_qqB.h"
#include "Diagrams/DC_Hjjj_NC_qq.h"

#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"

using namespace Herwig;


VBFNLOMEqq2hqqg::VBFNLOMEqq2hqqg():VBFNLOMEPP2hJetJetJet(){
}

VBFNLOMEqq2hqqg::VBFNLOMEqq2hqqg(int current, bool incoming1, bool incoming2, int decay, bool nw)
  : VBFNLOMEPP2hJetJetJet(current, incoming1, incoming2, decay, nw) {
}

VBFNLOMEqq2hqqg::~VBFNLOMEqq2hqqg(){}

IBPtr VBFNLOMEqq2hqqg::clone() const {
  return new_ptr(*this);
}

IBPtr VBFNLOMEqq2hqqg::fullclone() const {
  return new_ptr(*this);
}

void VBFNLOMEqq2hqqg::prepareMomenta(double (& pbar)[14][4],
					double (&qbar) [5],int) const{
  L5MomToDouble( (mePartonData()[0]->id() > 0 ? 
  		  lastMEMomenta()[0]  : 
  		  lastMEMomenta()[2]), &pbar[0][0]);

  L5MomToDouble( (mePartonData()[2]->id() > 0 ? 
  		  lastMEMomenta()[2] : 
  		  lastMEMomenta()[0]), &pbar[1][0]);

  L5MomToDouble( (mePartonData()[1]->id() > 0 ? 
  		  lastMEMomenta()[1] : 
  		  lastMEMomenta()[3]), &pbar[2][0]);

  L5MomToDouble( (mePartonData()[3]->id() > 0 ? 
  		  lastMEMomenta()[3] : 
  		  lastMEMomenta()[1]), &pbar[3][0]);

  L5MomToDouble( (lastMEMomenta()[4]), &qbar[0]);
  qbar[4]=0;

  pbar[4][0]=0;
  pbar[4][1]=0;
  pbar[4][2]=0;
  pbar[4][3]=0;

  if (NDecayProducts() == 1) {
    L5MomToDouble( (lastMEMomenta()[5]), &pbar[5][0]);
  }
  else if (NDecayProducts() == 2) {
    L5MomToDouble( (lastMEMomenta()[5]+lastMEMomenta()[6]), &pbar[5][0]);
  }
  else if (NDecayProducts() == 4) {
    L5MomToDouble( (lastMEMomenta()[5]+lastMEMomenta()[6]+lastMEMomenta()[7]+lastMEMomenta()[8]), &pbar[5][0]);
  }


  return;
} 

void VBFNLOMEqq2hqqg::initDiagramContainers() const {
  if ( theCurrent == neutral ){
    if ( theIncoming1 && theIncoming2 ) theDiagramContainer = RCPtr<DC_Hjjj_NC_qq>::Create(DC_Hjjj_NC_qq(const_cast<VBFNLOMEqq2hqqg*>(this)));
    if ( !theIncoming1 && theIncoming2 ) theDiagramContainer = RCPtr<DC_Hjjj_NC_qBq>::Create(DC_Hjjj_NC_qBq(const_cast<VBFNLOMEqq2hqqg*>(this)));
    if ( theIncoming1 && !theIncoming2 ) theDiagramContainer = RCPtr<DC_Hjjj_NC_qqB>::Create(DC_Hjjj_NC_qqB(const_cast<VBFNLOMEqq2hqqg*>(this)));
    if ( !theIncoming1 && !theIncoming2 ) theDiagramContainer = RCPtr<DC_Hjjj_NC_qBqB>::Create(DC_Hjjj_NC_qBqB(const_cast<VBFNLOMEqq2hqqg*>(this)));
  }
  else if ( theCurrent == charged ){
    if ( theIncoming1 && theIncoming2 ) theDiagramContainer = RCPtr<DC_Hjjj_CC_qq>::Create(DC_Hjjj_CC_qq(const_cast<VBFNLOMEqq2hqqg*>(this)));
    if ( !theIncoming1 && theIncoming2 ) theDiagramContainer = RCPtr<DC_Hjjj_CC_qBq>::Create(DC_Hjjj_CC_qBq(const_cast<VBFNLOMEqq2hqqg*>(this)));
    if ( theIncoming1 && !theIncoming2 ) theDiagramContainer = RCPtr<DC_Hjjj_CC_qqB>::Create(DC_Hjjj_CC_qqB(const_cast<VBFNLOMEqq2hqqg*>(this)));
    if ( !theIncoming1 && !theIncoming2 ) theDiagramContainer = RCPtr<DC_Hjjj_CC_qBqB>::Create(DC_Hjjj_CC_qBqB(const_cast<VBFNLOMEqq2hqqg*>(this)));
  }
  else
    assert(false);
}

void VBFNLOMEqq2hqqg::getDiagrams() const {
  if (theInteger == 0){
    if (!theDiagramContainer){
      initDiagramContainers();
    }
    allPossibleDiagrams = theDiagramContainer->getDiagrams();
    for (vector<DiagPtr>::const_iterator di = allPossibleDiagrams.begin();
	 di != allPossibleDiagrams.end() ; di++) {
      if ( allowedDiagram(*di) ) {
	add(*di);
      }
    }
  return;
  }

  bool useOldStyle = true;

  if (useOldStyle) {
    tcPDPtr Wplus = getParticleData(ParticleID::Wplus);
    tcPDPtr Wminus = getParticleData(ParticleID::Wminus);
    tcPDPtr Z0 = getParticleData(ParticleID::Z0);
    tcPDPtr h0 = getParticleData(ParticleID::h0);
    tcPDPtr g = getParticleData(ParticleID::g);

    PDVector QuarksAndAntiQuarks1,QuarksAndAntiQuarks2;
    for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	  q != theQuarkFlavours.end(); ++q ){
      if (requestedAsIncoming1(*q)) QuarksAndAntiQuarks1.push_back(*q);
      if (requestedAsIncoming1((**q).CC())) QuarksAndAntiQuarks1.push_back((**q).CC());
      if (requestedAsIncoming2(*q)) QuarksAndAntiQuarks2.push_back(*q);
      if (requestedAsIncoming2((**q).CC())) QuarksAndAntiQuarks2.push_back((**q).CC());
    }
    if (NDecayProducts() == 1 ) {
      if ( theCurrent == neutral ) {
	for ( PDVector::const_iterator q1 = QuarksAndAntiQuarks1.begin();
	      q1 != QuarksAndAntiQuarks1.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = QuarksAndAntiQuarks2.begin();
		q2 != QuarksAndAntiQuarks2.end(); ++q2 ) {
	    add(new_ptr((Tree2toNDiagram(5), *q1, *q1, Z0, Z0, *q2, 2, *q1,  4, *q2,  1, g,  3, h0,    -1)));
	    add(new_ptr((Tree2toNDiagram(5), *q1, Z0, Z0, *q2, *q2, 1, *q1,  3, *q2,  4, g,  2, h0,    -2)));
	    add(new_ptr((Tree2toNDiagram(4), *q1, Z0, Z0, *q2,  1, *q1,  5, *q1,  3, *q2,  5, g,  2, h0,  -3)));
	    add(new_ptr((Tree2toNDiagram(4), *q1, Z0, Z0, *q2,  3, *q2,  1, *q1,  5, *q2,  5, g,  2, h0,  -4)));
	  }
	}
      }

      else if ( theCurrent == charged ) {
	for ( PDVector::const_iterator q1 = QuarksAndAntiQuarks1.begin();
	      q1 != QuarksAndAntiQuarks1.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = QuarksAndAntiQuarks2.begin();
		q2 != QuarksAndAntiQuarks2.end(); ++q2 ) {
	    Charge chargediff = (**q1).charge() - (**q2).charge() ;
	    if ( chargediff != ZERO ) {
	      tcPDPtr q1prime = SU2Helper::SU2CC(*q1);
	      tcPDPtr q2prime = SU2Helper::SU2CC(*q2);
	      if ((**q1).charge()-(*q1prime).charge() == (*q2prime).charge()-(**q2).charge()){
		if ( chargediff > ZERO) {
		  add(new_ptr((Tree2toNDiagram(5), *q1, *q1, Wplus, Wplus, *q2, 2, q1prime,  4, q2prime,  1, g,  3, h0,    -1)));
		  add(new_ptr((Tree2toNDiagram(5), *q1, Wplus, Wplus, *q2, *q2, 1, q1prime,  3, q2prime,  4, g,  2, h0,    -2)));
		  add(new_ptr((Tree2toNDiagram(4), *q1, Wplus, Wplus, *q2,  1, q1prime,  5, q1prime,  3, q2prime,  5, g,  2, h0,  -3)));
		  add(new_ptr((Tree2toNDiagram(4), *q1, Wplus, Wplus, *q2,  3, q2prime,  1, q1prime,  5, q2prime,  5, g,  2, h0,  -4)));
		}
		else if ( chargediff < ZERO ) {
		  add(new_ptr((Tree2toNDiagram(5), *q1, *q1, Wminus, Wminus, *q2, 2, q1prime,  4, q2prime,  1, g,  3, h0,    -1)));
		  add(new_ptr((Tree2toNDiagram(5), *q1, Wminus, Wminus, *q2, *q2, 1, q1prime,  3, q2prime,  4, g,  2, h0,    -2)));
		  add(new_ptr((Tree2toNDiagram(4), *q1, Wminus, Wminus, *q2,  1, q1prime,  5, q1prime,  3, q2prime,  5, g,  2, h0,  -3)));
		  add(new_ptr((Tree2toNDiagram(4), *q1, Wminus, Wminus, *q2,  3, q2prime,  1, q1prime,  5, q2prime,  5, g,  2, h0,  -4)));
		}
	      }
	    }
	  }
	}
      }
      else {
	throw ThePEG::Exception() << "Please insert a pointer to W or Z "
				  << "boson as ExchangedBoson "
				  << "in the ME interface."
				  << ThePEG::Exception::abortnow;
      }
    }
    else if (NDecayProducts() == 2) {
      // get the decay products
      tcPDPtr dec1, dec2;

      if (DecayChannel() == 1) {
	dec1 = getParticleData(ParticleID::gamma);
	dec2 = getParticleData(ParticleID::gamma);
      }
      else if (DecayChannel() == 2) {
	dec1 = getParticleData(ParticleID::muplus);
	dec2 = getParticleData(ParticleID::muminus);
      }
      else if (DecayChannel() == 3) {
	dec1 = getParticleData(ParticleID::tauplus);
	dec2 = getParticleData(ParticleID::tauminus);
      }
      else if (DecayChannel() == 4) {
	dec1 = getParticleData(ParticleID::b);
	dec2 = getParticleData(ParticleID::bbar);
      }
      else   
	throw ThePEG::Exception() << "DecayChannel() doesn't fit to"
				  << "the result of NDecayProducts()"
				  << ThePEG::Exception::abortnow;
      if ( theCurrent == neutral ) {
	for ( PDVector::const_iterator q1 = QuarksAndAntiQuarks1.begin();
	      q1 != QuarksAndAntiQuarks1.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = QuarksAndAntiQuarks2.begin();
		q2 != QuarksAndAntiQuarks2.end(); ++q2 ) {
	    add(new_ptr((Tree2toNDiagram(5), *q1, *q1, Z0, Z0, *q2, 2, *q1,  4, *q2,  1, g,  3, h0,  9,dec1,  9,dec2,    -1)));
	    add(new_ptr((Tree2toNDiagram(5), *q1, Z0, Z0, *q2, *q2, 1, *q1,  3, *q2,  4, g,  2, h0,  9,dec1,  9,dec2,    -2)));
	    add(new_ptr((Tree2toNDiagram(4), *q1, Z0, Z0, *q2,  1, *q1,  5, *q1,  3, *q2,  5, g,  2, h0,  9,dec1,  9,dec2,  -3)));
	    add(new_ptr((Tree2toNDiagram(4), *q1, Z0, Z0, *q2,  3, *q2,  1, *q1,  5, *q2,  5, g,  2, h0,  9,dec1,  9,dec2,  -4)));
	  }
	}
      }

      else if ( theCurrent == charged ) {
	for ( PDVector::const_iterator q1 = QuarksAndAntiQuarks1.begin();
	      q1 != QuarksAndAntiQuarks1.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = QuarksAndAntiQuarks2.begin();
		q2 != QuarksAndAntiQuarks2.end(); ++q2 ) {
	    Charge chargediff = (**q1).charge() - (**q2).charge() ;
	    if ( chargediff != ZERO ) {
	      tcPDPtr q1prime = SU2Helper::SU2CC(*q1);
	      tcPDPtr q2prime = SU2Helper::SU2CC(*q2);
	      if ((**q1).charge()-(*q1prime).charge() == (*q2prime).charge()-(**q2).charge()){
		if ( chargediff > ZERO) {
		  add(new_ptr((Tree2toNDiagram(5), *q1, *q1, Wplus, Wplus, *q2, 2, q1prime,  4, q2prime,  1, g,  3, h0,  9,dec1,  9,dec2,    -1)));
		  add(new_ptr((Tree2toNDiagram(5), *q1, Wplus, Wplus, *q2, *q2, 1, q1prime,  3, q2prime,  4, g,  2, h0,  9,dec1,  9,dec2,    -2)));
		  add(new_ptr((Tree2toNDiagram(4), *q1, Wplus, Wplus, *q2,  1, q1prime,  5, q1prime,  3, q2prime,  5, g,  2, h0,  9,dec1,  9,dec2,  -3)));
		  add(new_ptr((Tree2toNDiagram(4), *q1, Wplus, Wplus, *q2,  3, q2prime,  1, q1prime,  5, q2prime,  5, g,  2, h0,  9,dec1,  9,dec2,  -4)));
		}
		else if ( chargediff < ZERO ) {
		  add(new_ptr((Tree2toNDiagram(5), *q1, *q1, Wminus, Wminus, *q2, 2, q1prime,  4, q2prime,  1, g,  3, h0,  9,dec1,  9,dec2,    -1)));
		  add(new_ptr((Tree2toNDiagram(5), *q1, Wminus, Wminus, *q2, *q2, 1, q1prime,  3, q2prime,  4, g,  2, h0,  9,dec1,  9,dec2,    -2)));
		  add(new_ptr((Tree2toNDiagram(4), *q1, Wminus, Wminus, *q2,  1, q1prime,  5, q1prime,  3, q2prime,  5, g,  2, h0,  9,dec1,  9,dec2,  -3)));
		  add(new_ptr((Tree2toNDiagram(4), *q1, Wminus, Wminus, *q2,  3, q2prime,  1, q1prime,  5, q2prime,  5, g,  2, h0,  9,dec1,  9,dec2,  -4)));
		}
	      }
	    }
	  }
	}
      }
      else {
	throw ThePEG::Exception() << "Please insert a pointer to W or Z "
				  << "boson as ExchangedBoson "
				  << "in the ME interface."
				  << ThePEG::Exception::abortnow;
      }
    }
    else if (NDecayProducts() == 4) {
      // get the decay products
      tcPDPtr VecBos1, VecBos2, dec1, dec2, dec3, dec4;

      if (DecayChannel() == 5) {
	VecBos1 = getParticleData(ParticleID::Wplus);
	VecBos2 = getParticleData(ParticleID::Wminus);
	dec1 = getParticleData(ParticleID::nu_e);
	dec2 = getParticleData(ParticleID::eplus);
	dec3 = getParticleData(ParticleID::nu_ebar);
	dec4 = getParticleData(ParticleID::eminus);
      }
      else if (DecayChannel() == 6) {
	VecBos1 = getParticleData(ParticleID::Z0);
	VecBos2 = getParticleData(ParticleID::Z0);
	dec1 = getParticleData(ParticleID::eminus);
	dec2 = getParticleData(ParticleID::eplus);
	dec3 = getParticleData(ParticleID::eminus);
	dec4 = getParticleData(ParticleID::eplus);
      }
      else if (DecayChannel() == 7) {
	VecBos1 = getParticleData(ParticleID::Z0);
	VecBos2 = getParticleData(ParticleID::Z0);
	dec1 = getParticleData(ParticleID::eminus);
	dec2 = getParticleData(ParticleID::eplus);
	dec3 = getParticleData(ParticleID::nu_e);
	dec4 = getParticleData(ParticleID::nu_ebar);
      }
      if ( theCurrent == neutral ) {
	for ( PDVector::const_iterator q1 = QuarksAndAntiQuarks1.begin();
	      q1 != QuarksAndAntiQuarks1.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = QuarksAndAntiQuarks2.begin();
		q2 != QuarksAndAntiQuarks2.end(); ++q2 ) {
	    add(new_ptr((Tree2toNDiagram(5), *q1, *q1, Z0, Z0, *q2, 2, *q1,  4, *q2,  1, g,  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4,    -1)));
	    add(new_ptr((Tree2toNDiagram(5), *q1, Z0, Z0, *q2, *q2, 1, *q1,  3, *q2,  4, g,  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4,    -2)));
	    add(new_ptr((Tree2toNDiagram(4), *q1, Z0, Z0, *q2,  1, *q1,  5, *q1,  3, *q2,  5, g,  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4,  -3)));
	    add(new_ptr((Tree2toNDiagram(4), *q1, Z0, Z0, *q2,  3, *q2,  1, *q1,  5, *q2,  5, g,  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4,  -4)));
	  }
	}
      }

      else if ( theCurrent == charged ) {
	for ( PDVector::const_iterator q1 = QuarksAndAntiQuarks1.begin();
	      q1 != QuarksAndAntiQuarks1.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = QuarksAndAntiQuarks2.begin();
		q2 != QuarksAndAntiQuarks2.end(); ++q2 ) {
	    Charge chargediff = (**q1).charge() - (**q2).charge() ;
	    if ( chargediff != ZERO ) {
	      tcPDPtr q1prime = SU2Helper::SU2CC(*q1);
	      tcPDPtr q2prime = SU2Helper::SU2CC(*q2);
	      if ((**q1).charge()-(*q1prime).charge() == (*q2prime).charge()-(**q2).charge()){
		if ( chargediff > ZERO) {
		  add(new_ptr((Tree2toNDiagram(5), *q1, *q1, Wplus, Wplus, *q2, 2, q1prime,  4, q2prime,  1, g,  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4,    -1)));
		  add(new_ptr((Tree2toNDiagram(5), *q1, Wplus, Wplus, *q2, *q2, 1, q1prime,  3, q2prime,  4, g,  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4,    -2)));
		  add(new_ptr((Tree2toNDiagram(4), *q1, Wplus, Wplus, *q2,  1, q1prime,  5, q1prime,  3, q2prime,  5, g,  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4,  -3)));
		  add(new_ptr((Tree2toNDiagram(4), *q1, Wplus, Wplus, *q2,  3, q2prime,  1, q1prime,  5, q2prime,  5, g,  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4,  -4)));
		}
		else if ( chargediff < ZERO ) {
		  add(new_ptr((Tree2toNDiagram(5), *q1, *q1, Wminus, Wminus, *q2, 2, q1prime,  4, q2prime,  1, g,  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4,    -1)));
		  add(new_ptr((Tree2toNDiagram(5), *q1, Wminus, Wminus, *q2, *q2, 1, q1prime,  3, q2prime,  4, g,  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4,    -2)));
		  add(new_ptr((Tree2toNDiagram(4), *q1, Wminus, Wminus, *q2,  1, q1prime,  5, q1prime,  3, q2prime,  5, g,  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4,  -3)));
		  add(new_ptr((Tree2toNDiagram(4), *q1, Wminus, Wminus, *q2,  3, q2prime,  1, q1prime,  5, q2prime,  5, g,  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4,  -4)));
		}
	      }
	    }
	  }
	}
      }
      else {
	throw ThePEG::Exception() << "Please insert a pointer to W or Z "
				  << "boson as ExchangedBoson "
				  << "in the ME interface."
				  << ThePEG::Exception::abortnow;
      }
    }
  }
  else { // ! useOldStyle

    PDPtr Wplus = getParticleData(ParticleID::Wplus); 
    PDPtr Wminus = getParticleData(ParticleID::Wminus); 
    PDPtr Z0 = getParticleData(ParticleID::Z0); 
    PDPtr h0 = getParticleData(ParticleID::h0); 
    PDPtr g = getParticleData(ParticleID::g); 
    
    for (int genU = 0; genU < 3; genU++) {
      for (int genL = 0; genL < 3; genL++) {
	PDPtr qu, qd, ku, kd;
	if (genU == 0){
	  qu = getParticleData(ParticleID::u);
	  qd = getParticleData(ParticleID::d);
	}
	else if (genU == 1){
	  qu = getParticleData(ParticleID::c);
	  qd = getParticleData(ParticleID::s);
	}
	else if (genU == 2){
	  qu = getParticleData(ParticleID::t);
	  qd = getParticleData(ParticleID::b);
	}
	if (genL == 0){
	  ku = getParticleData(ParticleID::u);
	  kd = getParticleData(ParticleID::d);
	}
	else if (genL == 1){
	  ku = getParticleData(ParticleID::c);
	  kd = getParticleData(ParticleID::s);
	}
	else if (genL == 2){
	  ku = getParticleData(ParticleID::t);
	  kd = getParticleData(ParticleID::b);
	}

	tcPDPtr quB = (*qu).CC();
	tcPDPtr qdB = (*qd).CC();

	tcPDPtr quC = SU2Helper::SU2CC(qu);
	tcPDPtr qdC = SU2Helper::SU2CC(qd);
	tcPDPtr quBC = SU2Helper::SU2CC(quB);
	tcPDPtr qdBC = SU2Helper::SU2CC(qdB);

	tcPDPtr kuB = (*ku).CC();
	tcPDPtr kdB = (*kd).CC();

	tcPDPtr kuC = SU2Helper::SU2CC(ku);
	tcPDPtr kdC = SU2Helper::SU2CC(kd);
	tcPDPtr kuBC = SU2Helper::SU2CC(kuB);
	tcPDPtr kdBC = SU2Helper::SU2CC(kdB);
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, Wminus, Wminus, kdB, kd, 1, qdBC, 3, kdC, 4, g, 2, h0, -1)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Wminus, Wminus, kd, 3, kdC, 1, qdBC, 5, kdC, 5, g, 2, h0, -2)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, qdB, Wminus, Wminus, kd, 2, qdBC, 4, kdC, 1, g, 3, h0, -3)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Wminus, Wminus, kd, 1, qdBC, 5, qdBC, 3, kdC, 5, g, 2, h0, -4)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, Wminus, Wminus, ku, kuB, 1, qdBC, 3, kuBC, 4, g, 2, h0, -5)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Wminus, Wminus, kuB, 3, kuBC, 1, qdBC, 5, kuBC, 5, g, 2, h0, -6)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, qdB, Wminus, Wminus, kuB, 2, qdBC, 4, kuBC, 1, g, 3, h0, -7)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Wminus, Wminus, kuB, 1, qdBC, 5, qdBC, 3, kuBC, 5, g, 2, h0, -8)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Wplus, Wplus, kd, kdB, 1, qdC, 3, kdBC, 4, g, 2, h0, -9)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Wplus, Wplus, kdB, 3, kdBC, 1, qdC, 5, kdBC, 5, g, 2, h0, -10)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, qd, Wplus, Wplus, kdB, 2, qdC, 4, kdBC, 1, g, 3, h0, -11)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Wplus, Wplus, kdB, 1, qdC, 5, qdC, 3, kdBC, 5, g, 2, h0, -12)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Wplus, Wplus, kuB, ku, 1, qdC, 3, kuC, 4, g, 2, h0, -13)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Wplus, Wplus, ku, 3, kuC, 1, qdC, 5, kuC, 5, g, 2, h0, -14)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, qd, Wplus, Wplus, ku, 2, qdC, 4, kuC, 1, g, 3, h0, -15)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Wplus, Wplus, ku, 1, qdC, 5, qdC, 3, kuC, 5, g, 2, h0, -16)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, Wplus, Wplus, kd, kdB, 1, quBC, 3, kdBC, 4, g, 2, h0, -17)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Wplus, Wplus, kdB, 3, kdBC, 1, quBC, 5, kdBC, 5, g, 2, h0, -18)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, quB, Wplus, Wplus, kdB, 2, quBC, 4, kdBC, 1, g, 3, h0, -19)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Wplus, Wplus, kdB, 1, quBC, 5, quBC, 3, kdBC, 5, g, 2, h0, -20)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, Wplus, Wplus, kuB, ku, 1, quBC, 3, kuC, 4, g, 2, h0, -21)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Wplus, Wplus, ku, 3, kuC, 1, quBC, 5, kuC, 5, g, 2, h0, -22)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, quB, Wplus, Wplus, ku, 2, quBC, 4, kuC, 1, g, 3, h0, -23)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Wplus, Wplus, ku, 1, quBC, 5, quBC, 3, kuC, 5, g, 2, h0, -24)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Wminus, Wminus, kdB, kd, 1, quC, 3, kdC, 4, g, 2, h0, -25)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Wminus, Wminus, kd, 3, kdC, 1, quC, 5, kdC, 5, g, 2, h0, -26)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, qu, Wminus, Wminus, kd, 2, quC, 4, kdC, 1, g, 3, h0, -27)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Wminus, Wminus, kd, 1, quC, 5, quC, 3, kdC, 5, g, 2, h0, -28)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Wminus, Wminus, ku, kuB, 1, quC, 3, kuBC, 4, g, 2, h0, -29)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Wminus, Wminus, kuB, 3, kuBC, 1, quC, 5, kuBC, 5, g, 2, h0, -30)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, qu, Wminus, Wminus, kuB, 2, quC, 4, kuBC, 1, g, 3, h0, -31)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Wminus, Wminus, kuB, 1, quC, 5, quC, 3, kuBC, 5, g, 2, h0, -32)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, Z0, Z0, kdB, kd, 1, qdB, 3, kd, 4, g, 2, h0, -33)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Z0, Z0, kd, 3, kd, 1, qdB, 5, kd, 5, g, 2, h0, -34)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, qdB, Z0, Z0, kd, 2, qdB, 4, kd, 1, g, 3, h0, -35)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Z0, Z0, kd, 1, qdB, 5, qdB, 3, kd, 5, g, 2, h0, -36)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, Z0, Z0, kd, kdB, 1, qdB, 3, kdB, 4, g, 2, h0, -37)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Z0, Z0, kdB, 3, kdB, 1, qdB, 5, kdB, 5, g, 2, h0, -38)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, qdB, Z0, Z0, kdB, 2, qdB, 4, kdB, 1, g, 3, h0, -39)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Z0, Z0, kdB, 1, qdB, 5, qdB, 3, kdB, 5, g, 2, h0, -40)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, Z0, Z0, kuB, ku, 1, qdB, 3, ku, 4, g, 2, h0, -41)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Z0, Z0, ku, 3, ku, 1, qdB, 5, ku, 5, g, 2, h0, -42)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, qdB, Z0, Z0, ku, 2, qdB, 4, ku, 1, g, 3, h0, -43)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Z0, Z0, ku, 1, qdB, 5, qdB, 3, ku, 5, g, 2, h0, -44)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, Z0, Z0, ku, kuB, 1, qdB, 3, kuB, 4, g, 2, h0, -45)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Z0, Z0, kuB, 3, kuB, 1, qdB, 5, kuB, 5, g, 2, h0, -46)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, qdB, Z0, Z0, kuB, 2, qdB, 4, kuB, 1, g, 3, h0, -47)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Z0, Z0, kuB, 1, qdB, 5, qdB, 3, kuB, 5, g, 2, h0, -48)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Z0, Z0, kdB, kd, 1, qd, 3, kd, 4, g, 2, h0, -49)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, kd, 3, kd, 1, qd, 5, kd, 5, g, 2, h0, -50)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, qd, Z0, Z0, kd, 2, qd, 4, kd, 1, g, 3, h0, -51)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, kd, 1, qd, 5, qd, 3, kd, 5, g, 2, h0, -52)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Z0, Z0, kd, kdB, 1, qd, 3, kdB, 4, g, 2, h0, -53)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, kdB, 3, kdB, 1, qd, 5, kdB, 5, g, 2, h0, -54)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, qd, Z0, Z0, kdB, 2, qd, 4, kdB, 1, g, 3, h0, -55)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, kdB, 1, qd, 5, qd, 3, kdB, 5, g, 2, h0, -56)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Z0, Z0, kuB, ku, 1, qd, 3, ku, 4, g, 2, h0, -57)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, ku, 3, ku, 1, qd, 5, ku, 5, g, 2, h0, -58)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, qd, Z0, Z0, ku, 2, qd, 4, ku, 1, g, 3, h0, -59)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, ku, 1, qd, 5, qd, 3, ku, 5, g, 2, h0, -60)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Z0, Z0, ku, kuB, 1, qd, 3, kuB, 4, g, 2, h0, -61)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, kuB, 3, kuB, 1, qd, 5, kuB, 5, g, 2, h0, -62)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, qd, Z0, Z0, kuB, 2, qd, 4, kuB, 1, g, 3, h0, -63)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, kuB, 1, qd, 5, qd, 3, kuB, 5, g, 2, h0, -64)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, Z0, Z0, kdB, kd, 1, quB, 3, kd, 4, g, 2, h0, -65)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Z0, Z0, kd, 3, kd, 1, quB, 5, kd, 5, g, 2, h0, -66)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, quB, Z0, Z0, kd, 2, quB, 4, kd, 1, g, 3, h0, -67)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Z0, Z0, kd, 1, quB, 5, quB, 3, kd, 5, g, 2, h0, -68)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, Z0, Z0, kd, kdB, 1, quB, 3, kdB, 4, g, 2, h0, -69)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Z0, Z0, kdB, 3, kdB, 1, quB, 5, kdB, 5, g, 2, h0, -70)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, quB, Z0, Z0, kdB, 2, quB, 4, kdB, 1, g, 3, h0, -71)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Z0, Z0, kdB, 1, quB, 5, quB, 3, kdB, 5, g, 2, h0, -72)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, Z0, Z0, kuB, ku, 1, quB, 3, ku, 4, g, 2, h0, -73)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Z0, Z0, ku, 3, ku, 1, quB, 5, ku, 5, g, 2, h0, -74)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, quB, Z0, Z0, ku, 2, quB, 4, ku, 1, g, 3, h0, -75)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Z0, Z0, ku, 1, quB, 5, quB, 3, ku, 5, g, 2, h0, -76)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, Z0, Z0, ku, kuB, 1, quB, 3, kuB, 4, g, 2, h0, -77)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Z0, Z0, kuB, 3, kuB, 1, quB, 5, kuB, 5, g, 2, h0, -78)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, quB, Z0, Z0, kuB, 2, quB, 4, kuB, 1, g, 3, h0, -79)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Z0, Z0, kuB, 1, quB, 5, quB, 3, kuB, 5, g, 2, h0, -80)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Z0, Z0, kdB, kd, 1, qu, 3, kd, 4, g, 2, h0, -81)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, kd, 3, kd, 1, qu, 5, kd, 5, g, 2, h0, -82)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, qu, Z0, Z0, kd, 2, qu, 4, kd, 1, g, 3, h0, -83)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, kd, 1, qu, 5, qu, 3, kd, 5, g, 2, h0, -84)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Z0, Z0, kd, kdB, 1, qu, 3, kdB, 4, g, 2, h0, -85)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, kdB, 3, kdB, 1, qu, 5, kdB, 5, g, 2, h0, -86)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, qu, Z0, Z0, kdB, 2, qu, 4, kdB, 1, g, 3, h0, -87)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, kdB, 1, qu, 5, qu, 3, kdB, 5, g, 2, h0, -88)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Z0, Z0, kuB, ku, 1, qu, 3, ku, 4, g, 2, h0, -89)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, ku, 3, ku, 1, qu, 5, ku, 5, g, 2, h0, -90)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, qu, Z0, Z0, ku, 2, qu, 4, ku, 1, g, 3, h0, -91)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, ku, 1, qu, 5, qu, 3, ku, 5, g, 2, h0, -92)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Z0, Z0, ku, kuB, 1, qu, 3, kuB, 4, g, 2, h0, -93)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, kuB, 3, kuB, 1, qu, 5, kuB, 5, g, 2, h0, -94)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, qu, Z0, Z0, kuB, 2, qu, 4, kuB, 1, g, 3, h0, -95)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, kuB, 1, qu, 5, qu, 3, kuB, 5, g, 2, h0, -96)));
      }
    }
    //now check which of those are really requested by the interfaced parameters
    for (vector<DiagPtr>::const_iterator di = allPossibleDiagrams.begin();
	 di != allPossibleDiagrams.end() ; di++) {
      if ( allowedDiagram(*di) ) {
	add(*di);
      }
    }
    return;
  }
}



void VBFNLOMEqq2hqqg::doinit(){ 
  VBFNLOMEVVJJNeutralBase::doinit();
  for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	q != theQuarkFlavours.end(); ++q ){
    if ( (**q).mass() != ZERO )
      Throw<InitException>() << "The matrix element '"
			     << name() << "' is only capable of "
			     << "producing massless quarks.";
    if ( (**q).id() < 1 || (**q).id() > 6) Throw<InitException>() << "Please insert  "
								  << "only quarks in the ME interface of" 
								  << name() <<".";
				       
  }
}

void VBFNLOMEqq2hqqg::doinitrun(){ 
  VBFNLOMEVVJJNeutralBase::doinitrun();
}

Selector<MEBase::DiagramIndex> 
VBFNLOMEqq2hqqg::diagrams(const DiagramVector & ) const {
  Selector<MEBase::DiagramIndex> sel;
  double w0 = sqr( sqr(generator()->maximumCMEnergy())/(2.*(meMomenta()[0]*meMomenta()[4])) );
  double w1 = sqr( sqr(generator()->maximumCMEnergy())/(2.*(meMomenta()[2]*meMomenta()[4])) );
  double w2 = sqr( sqr(generator()->maximumCMEnergy())/(2.*(meMomenta()[1]*meMomenta()[4])) );
  double w3 = sqr( sqr(generator()->maximumCMEnergy())/(2.*(meMomenta()[3]*meMomenta()[4])) );
  sel.insert(w0,0);
  sel.insert(w1,1);
  sel.insert(w2,2);
  sel.insert(w3,3);
  return sel;
}

Selector<const ColourLines *>
VBFNLOMEqq2hqqg::colourGeometries(tcDiagPtr diag) const {
  if (theInteger == 0){
    if (!theDiagramContainer){
      initDiagramContainers();
    }
    return theDiagramContainer->colourGeometries(diag);
  }

  bool useOldStyle = true;
  if (useOldStyle) {
    if (DecayChannel() != 4) {
      static const ColourLines cUL[4] = {ColourLines(" 1  8, -8  2  6,  5  7"),
					 ColourLines(" 1  6,  5  8, -8  4  7"),
					 ColourLines(" 1  5  8, -8  6,  4  7"),
					 ColourLines(" 1  6,  4  5  8, -8  7")};
      static const ColourLines cULbar[4] = {ColourLines(" 1  8, -8  2  6, -5 -7"),
					    ColourLines(" 1  6, -5 -8,  8 -4 -7"),
					    ColourLines(" 1  5  8, -8  6, -4 -7"),
					    ColourLines(" 1  6, -4 -5 -8,  8 -7")};
      static const ColourLines cUbarL[4] = {ColourLines("-1 -8,  8 -2 -6,  5  7"),
					    ColourLines("-1 -6,  5  8, -8  4  7"),
					    ColourLines("-1 -5 -8,  8 -6,  4  7"),
					    ColourLines("-1 -6,  4  5  8, -8  7")};
      static const ColourLines cUbarLbar[4] = {ColourLines("-1 -8,  8 -2 -6, -5 -7"),
					       ColourLines("-1 -6, -5 -8,  8 -4 -7"),
					       ColourLines("-1 -5 -8,  8 -6, -4 -7"),
					       ColourLines("-1 -6, -4 -5 -8,  8 -7")};   

      Selector<const ColourLines *> sel;
      if ( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() > 0)
	sel.insert(1.0, &cUL[abs(diag->id())-1]);
      else if ( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() < 0)
	sel.insert(1.0, &cULbar[abs(diag->id())-1]);
      else if ( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() > 0)
	sel.insert(1.0, &cUbarL[abs(diag->id())-1]);
      else if ( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() < 0)
	sel.insert(1.0, &cUbarLbar[abs(diag->id())-1]);
      return sel;
    }
    else {
      static const ColourLines cUL[4] = {ColourLines(" 1  8, -8  2  6,  5  7,  10 -11"),
					 ColourLines(" 1  6,  5  8, -8  4  7,  10 -11"),
					 ColourLines(" 1  5  8, -8  6,  4  7,  10 -11"),
					 ColourLines(" 1  6,  4  5  8, -8  7,  10 -11")};
      static const ColourLines cULbar[4] = {ColourLines(" 1  8, -8  2  6, -5 -7,  10 -11"),
					    ColourLines(" 1  6, -5 -8,  8 -4 -7,  10 -11"),
					    ColourLines(" 1  5  8, -8  6, -4 -7,  10 -11"),
					    ColourLines(" 1  6, -4 -5 -8,  8 -7,  10 -11")};
      static const ColourLines cUbarL[4] = {ColourLines("-1 -8,  8 -2 -6,  5  7,  10 -11"),
					    ColourLines("-1 -6,  5  8, -8  4  7,  10 -11"),
					    ColourLines("-1 -5 -8,  8 -6,  4  7,  10 -11"),
					    ColourLines("-1 -6,  4  5  8, -8  7,  10 -11")};
      static const ColourLines cUbarLbar[4] = {ColourLines("-1 -8,  8 -2 -6, -5 -7,  10 -11"),
					       ColourLines("-1 -6, -5 -8,  8 -4 -7,  10 -11"),
					       ColourLines("-1 -5 -8,  8 -6, -4 -7,  10 -11"),
					       ColourLines("-1 -6, -4 -5 -8,  8 -7,  10 -11")};

      Selector<const ColourLines *> sel;
      if ( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() > 0)
	sel.insert(1.0, &cUL[abs(diag->id())-1]);
      else if ( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() < 0)
	sel.insert(1.0, &cULbar[abs(diag->id())-1]);
      else if ( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() > 0)
	sel.insert(1.0, &cUbarL[abs(diag->id())-1]);
      else if ( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() < 0)
	sel.insert(1.0, &cUbarLbar[abs(diag->id())-1]);
      return sel;
    }
  }
  else{
    static const ColourLines diag1[1] = { ColourLines("-1 -6, 5 8, 7 -4 -8") }; 
    static const ColourLines diag2[1] = { ColourLines("-1 -6, 4 5 8, 7 -8") }; 
    static const ColourLines diag3[1] = { ColourLines("-1 -8, 5 7, -6 -2 8") }; 
    static const ColourLines diag4[1] = { ColourLines("-1 -5 -8, 4 7, -6 8") }; 
    static const ColourLines diag5[1] = { ColourLines("-1 -6, -5 -8, -7 4 8") }; 
    static const ColourLines diag6[1] = { ColourLines("-1 -6, -4 -5 -8, -7 8") }; 
    static const ColourLines diag7[1] = { ColourLines("-1 -8, -5 -7, -6 -2 8") }; 
    static const ColourLines diag8[1] = { ColourLines("-1 -5 -8, -4 -7, -6 8") }; 
    static const ColourLines diag9[1] = { ColourLines("1 6, -5 -8, -7 4 8") }; 
    static const ColourLines diag10[1] = { ColourLines("1 6, -4 -5 -8, -7 8") }; 
    static const ColourLines diag11[1] = { ColourLines("1 8, -5 -7, 6 2 -8") }; 
    static const ColourLines diag12[1] = { ColourLines("1 5 8, -4 -7, 6 -8") }; 
    static const ColourLines diag13[1] = { ColourLines("1 6, 5 8, 7 -4 -8") }; 
    static const ColourLines diag14[1] = { ColourLines("1 6, 4 5 8, 7 -8") }; 
    static const ColourLines diag15[1] = { ColourLines("1 8, 5 7, 6 2 -8") }; 
    static const ColourLines diag16[1] = { ColourLines("1 5 8, 4 7, 6 -8") }; 
    static const ColourLines diag17[1] = { ColourLines("-1 -6, -5 -8, -7 4 8") }; 
    static const ColourLines diag18[1] = { ColourLines("-1 -6, -4 -5 -8, -7 8") }; 
    static const ColourLines diag19[1] = { ColourLines("-1 -8, -5 -7, -6 -2 8") }; 
    static const ColourLines diag20[1] = { ColourLines("-1 -5 -8, -4 -7, -6 8") }; 
    static const ColourLines diag21[1] = { ColourLines("-1 -6, 5 8, 7 -4 -8") }; 
    static const ColourLines diag22[1] = { ColourLines("-1 -6, 4 5 8, 7 -8") }; 
    static const ColourLines diag23[1] = { ColourLines("-1 -8, 5 7, -6 -2 8") }; 
    static const ColourLines diag24[1] = { ColourLines("-1 -5 -8, 4 7, -6 8") }; 
    static const ColourLines diag25[1] = { ColourLines("1 6, 5 8, 7 -4 -8") }; 
    static const ColourLines diag26[1] = { ColourLines("1 6, 4 5 8, 7 -8") }; 
    static const ColourLines diag27[1] = { ColourLines("1 8, 5 7, 6 2 -8") }; 
    static const ColourLines diag28[1] = { ColourLines("1 5 8, 4 7, 6 -8") }; 
    static const ColourLines diag29[1] = { ColourLines("1 6, -5 -8, -7 4 8") }; 
    static const ColourLines diag30[1] = { ColourLines("1 6, -4 -5 -8, -7 8") }; 
    static const ColourLines diag31[1] = { ColourLines("1 8, -5 -7, 6 2 -8") }; 
    static const ColourLines diag32[1] = { ColourLines("1 5 8, -4 -7, 6 -8") }; 
    static const ColourLines diag33[1] = { ColourLines("-1 -6, 5 8, 7 -4 -8") }; 
    static const ColourLines diag34[1] = { ColourLines("-1 -6, 4 5 8, 7 -8") }; 
    static const ColourLines diag35[1] = { ColourLines("-1 -8, 5 7, -6 -2 8") }; 
    static const ColourLines diag36[1] = { ColourLines("-1 -5 -8, 4 7, -6 8") }; 
    static const ColourLines diag37[1] = { ColourLines("-1 -6, -5 -8, -7 4 8") }; 
    static const ColourLines diag38[1] = { ColourLines("-1 -6, -4 -5 -8, -7 8") }; 
    static const ColourLines diag39[1] = { ColourLines("-1 -8, -5 -7, -6 -2 8") }; 
    static const ColourLines diag40[1] = { ColourLines("-1 -5 -8, -4 -7, -6 8") }; 
    static const ColourLines diag41[1] = { ColourLines("-1 -6, 5 8, 7 -4 -8") }; 
    static const ColourLines diag42[1] = { ColourLines("-1 -6, 4 5 8, 7 -8") }; 
    static const ColourLines diag43[1] = { ColourLines("-1 -8, 5 7, -6 -2 8") }; 
    static const ColourLines diag44[1] = { ColourLines("-1 -5 -8, 4 7, -6 8") }; 
    static const ColourLines diag45[1] = { ColourLines("-1 -6, -5 -8, -7 4 8") }; 
    static const ColourLines diag46[1] = { ColourLines("-1 -6, -4 -5 -8, -7 8") }; 
    static const ColourLines diag47[1] = { ColourLines("-1 -8, -5 -7, -6 -2 8") }; 
    static const ColourLines diag48[1] = { ColourLines("-1 -5 -8, -4 -7, -6 8") }; 
    static const ColourLines diag49[1] = { ColourLines("1 6, 5 8, 7 -4 -8") }; 
    static const ColourLines diag50[1] = { ColourLines("1 6, 4 5 8, 7 -8") }; 
    static const ColourLines diag51[1] = { ColourLines("1 8, 5 7, 6 2 -8") }; 
    static const ColourLines diag52[1] = { ColourLines("1 5 8, 4 7, 6 -8") }; 
    static const ColourLines diag53[1] = { ColourLines("1 6, -5 -8, -7 4 8") }; 
    static const ColourLines diag54[1] = { ColourLines("1 6, -4 -5 -8, -7 8") }; 
    static const ColourLines diag55[1] = { ColourLines("1 8, -5 -7, 6 2 -8") }; 
    static const ColourLines diag56[1] = { ColourLines("1 5 8, -4 -7, 6 -8") }; 
    static const ColourLines diag57[1] = { ColourLines("1 6, 5 8, 7 -4 -8") }; 
    static const ColourLines diag58[1] = { ColourLines("1 6, 4 5 8, 7 -8") }; 
    static const ColourLines diag59[1] = { ColourLines("1 8, 5 7, 6 2 -8") }; 
    static const ColourLines diag60[1] = { ColourLines("1 5 8, 4 7, 6 -8") }; 
    static const ColourLines diag61[1] = { ColourLines("1 6, -5 -8, -7 4 8") }; 
    static const ColourLines diag62[1] = { ColourLines("1 6, -4 -5 -8, -7 8") }; 
    static const ColourLines diag63[1] = { ColourLines("1 8, -5 -7, 6 2 -8") }; 
    static const ColourLines diag64[1] = { ColourLines("1 5 8, -4 -7, 6 -8") }; 
    static const ColourLines diag65[1] = { ColourLines("-1 -6, 5 8, 7 -4 -8") }; 
    static const ColourLines diag66[1] = { ColourLines("-1 -6, 4 5 8, 7 -8") }; 
    static const ColourLines diag67[1] = { ColourLines("-1 -8, 5 7, -6 -2 8") }; 
    static const ColourLines diag68[1] = { ColourLines("-1 -5 -8, 4 7, -6 8") }; 
    static const ColourLines diag69[1] = { ColourLines("-1 -6, -5 -8, -7 4 8") }; 
    static const ColourLines diag70[1] = { ColourLines("-1 -6, -4 -5 -8, -7 8") }; 
    static const ColourLines diag71[1] = { ColourLines("-1 -8, -5 -7, -6 -2 8") }; 
    static const ColourLines diag72[1] = { ColourLines("-1 -5 -8, -4 -7, -6 8") }; 
    static const ColourLines diag73[1] = { ColourLines("-1 -6, 5 8, 7 -4 -8") }; 
    static const ColourLines diag74[1] = { ColourLines("-1 -6, 4 5 8, 7 -8") }; 
    static const ColourLines diag75[1] = { ColourLines("-1 -8, 5 7, -6 -2 8") }; 
    static const ColourLines diag76[1] = { ColourLines("-1 -5 -8, 4 7, -6 8") }; 
    static const ColourLines diag77[1] = { ColourLines("-1 -6, -5 -8, -7 4 8") }; 
    static const ColourLines diag78[1] = { ColourLines("-1 -6, -4 -5 -8, -7 8") }; 
    static const ColourLines diag79[1] = { ColourLines("-1 -8, -5 -7, -6 -2 8") }; 
    static const ColourLines diag80[1] = { ColourLines("-1 -5 -8, -4 -7, -6 8") }; 
    static const ColourLines diag81[1] = { ColourLines("1 6, 5 8, 7 -4 -8") }; 
    static const ColourLines diag82[1] = { ColourLines("1 6, 4 5 8, 7 -8") }; 
    static const ColourLines diag83[1] = { ColourLines("1 8, 5 7, 6 2 -8") }; 
    static const ColourLines diag84[1] = { ColourLines("1 5 8, 4 7, 6 -8") }; 
    static const ColourLines diag85[1] = { ColourLines("1 6, -5 -8, -7 4 8") }; 
    static const ColourLines diag86[1] = { ColourLines("1 6, -4 -5 -8, -7 8") }; 
    static const ColourLines diag87[1] = { ColourLines("1 8, -5 -7, 6 2 -8") }; 
    static const ColourLines diag88[1] = { ColourLines("1 5 8, -4 -7, 6 -8") }; 
    static const ColourLines diag89[1] = { ColourLines("1 6, 5 8, 7 -4 -8") }; 
    static const ColourLines diag90[1] = { ColourLines("1 6, 4 5 8, 7 -8") }; 
    static const ColourLines diag91[1] = { ColourLines("1 8, 5 7, 6 2 -8") }; 
    static const ColourLines diag92[1] = { ColourLines("1 5 8, 4 7, 6 -8") }; 
    static const ColourLines diag93[1] = { ColourLines("1 6, -5 -8, -7 4 8") }; 
    static const ColourLines diag94[1] = { ColourLines("1 6, -4 -5 -8, -7 8") }; 
    static const ColourLines diag95[1] = { ColourLines("1 8, -5 -7, 6 2 -8") }; 
    static const ColourLines diag96[1] = { ColourLines("1 5 8, -4 -7, 6 -8") }; 

    Selector <const ColourLines *> sel;

    if( diag->id() == -1 ) { sel.insert( 1.0,  &(diag1[0]) ); } 
    else if( diag->id() == -2 ) { sel.insert( 1.0,  &(diag2[0]) ); } 
    else if( diag->id() == -3 ) { sel.insert( 1.0,  &(diag3[0]) ); } 
    else if( diag->id() == -4 ) { sel.insert( 1.0,  &(diag4[0]) ); } 
    else if( diag->id() == -5 ) { sel.insert( 1.0,  &(diag5[0]) ); } 
    else if( diag->id() == -6 ) { sel.insert( 1.0,  &(diag6[0]) ); } 
    else if( diag->id() == -7 ) { sel.insert( 1.0,  &(diag7[0]) ); } 
    else if( diag->id() == -8 ) { sel.insert( 1.0,  &(diag8[0]) ); } 
    else if( diag->id() == -9 ) { sel.insert( 1.0,  &(diag9[0]) ); } 
    else if( diag->id() == -10 ) { sel.insert( 1.0,  &(diag10[0]) ); } 
    else if( diag->id() == -11 ) { sel.insert( 1.0,  &(diag11[0]) ); } 
    else if( diag->id() == -12 ) { sel.insert( 1.0,  &(diag12[0]) ); } 
    else if( diag->id() == -13 ) { sel.insert( 1.0,  &(diag13[0]) ); } 
    else if( diag->id() == -14 ) { sel.insert( 1.0,  &(diag14[0]) ); } 
    else if( diag->id() == -15 ) { sel.insert( 1.0,  &(diag15[0]) ); } 
    else if( diag->id() == -16 ) { sel.insert( 1.0,  &(diag16[0]) ); } 
    else if( diag->id() == -17 ) { sel.insert( 1.0,  &(diag17[0]) ); } 
    else if( diag->id() == -18 ) { sel.insert( 1.0,  &(diag18[0]) ); } 
    else if( diag->id() == -19 ) { sel.insert( 1.0,  &(diag19[0]) ); } 
    else if( diag->id() == -20 ) { sel.insert( 1.0,  &(diag20[0]) ); } 
    else if( diag->id() == -21 ) { sel.insert( 1.0,  &(diag21[0]) ); } 
    else if( diag->id() == -22 ) { sel.insert( 1.0,  &(diag22[0]) ); } 
    else if( diag->id() == -23 ) { sel.insert( 1.0,  &(diag23[0]) ); } 
    else if( diag->id() == -24 ) { sel.insert( 1.0,  &(diag24[0]) ); } 
    else if( diag->id() == -25 ) { sel.insert( 1.0,  &(diag25[0]) ); } 
    else if( diag->id() == -26 ) { sel.insert( 1.0,  &(diag26[0]) ); } 
    else if( diag->id() == -27 ) { sel.insert( 1.0,  &(diag27[0]) ); } 
    else if( diag->id() == -28 ) { sel.insert( 1.0,  &(diag28[0]) ); } 
    else if( diag->id() == -29 ) { sel.insert( 1.0,  &(diag29[0]) ); } 
    else if( diag->id() == -30 ) { sel.insert( 1.0,  &(diag30[0]) ); } 
    else if( diag->id() == -31 ) { sel.insert( 1.0,  &(diag31[0]) ); } 
    else if( diag->id() == -32 ) { sel.insert( 1.0,  &(diag32[0]) ); } 
    else if( diag->id() == -33 ) { sel.insert( 1.0,  &(diag33[0]) ); } 
    else if( diag->id() == -34 ) { sel.insert( 1.0,  &(diag34[0]) ); } 
    else if( diag->id() == -35 ) { sel.insert( 1.0,  &(diag35[0]) ); } 
    else if( diag->id() == -36 ) { sel.insert( 1.0,  &(diag36[0]) ); } 
    else if( diag->id() == -37 ) { sel.insert( 1.0,  &(diag37[0]) ); } 
    else if( diag->id() == -38 ) { sel.insert( 1.0,  &(diag38[0]) ); } 
    else if( diag->id() == -39 ) { sel.insert( 1.0,  &(diag39[0]) ); } 
    else if( diag->id() == -40 ) { sel.insert( 1.0,  &(diag40[0]) ); } 
    else if( diag->id() == -41 ) { sel.insert( 1.0,  &(diag41[0]) ); } 
    else if( diag->id() == -42 ) { sel.insert( 1.0,  &(diag42[0]) ); } 
    else if( diag->id() == -43 ) { sel.insert( 1.0,  &(diag43[0]) ); } 
    else if( diag->id() == -44 ) { sel.insert( 1.0,  &(diag44[0]) ); } 
    else if( diag->id() == -45 ) { sel.insert( 1.0,  &(diag45[0]) ); } 
    else if( diag->id() == -46 ) { sel.insert( 1.0,  &(diag46[0]) ); } 
    else if( diag->id() == -47 ) { sel.insert( 1.0,  &(diag47[0]) ); } 
    else if( diag->id() == -48 ) { sel.insert( 1.0,  &(diag48[0]) ); } 
    else if( diag->id() == -49 ) { sel.insert( 1.0,  &(diag49[0]) ); } 
    else if( diag->id() == -50 ) { sel.insert( 1.0,  &(diag50[0]) ); } 
    else if( diag->id() == -51 ) { sel.insert( 1.0,  &(diag51[0]) ); } 
    else if( diag->id() == -52 ) { sel.insert( 1.0,  &(diag52[0]) ); } 
    else if( diag->id() == -53 ) { sel.insert( 1.0,  &(diag53[0]) ); } 
    else if( diag->id() == -54 ) { sel.insert( 1.0,  &(diag54[0]) ); } 
    else if( diag->id() == -55 ) { sel.insert( 1.0,  &(diag55[0]) ); } 
    else if( diag->id() == -56 ) { sel.insert( 1.0,  &(diag56[0]) ); } 
    else if( diag->id() == -57 ) { sel.insert( 1.0,  &(diag57[0]) ); } 
    else if( diag->id() == -58 ) { sel.insert( 1.0,  &(diag58[0]) ); } 
    else if( diag->id() == -59 ) { sel.insert( 1.0,  &(diag59[0]) ); } 
    else if( diag->id() == -60 ) { sel.insert( 1.0,  &(diag60[0]) ); } 
    else if( diag->id() == -61 ) { sel.insert( 1.0,  &(diag61[0]) ); } 
    else if( diag->id() == -62 ) { sel.insert( 1.0,  &(diag62[0]) ); } 
    else if( diag->id() == -63 ) { sel.insert( 1.0,  &(diag63[0]) ); } 
    else if( diag->id() == -64 ) { sel.insert( 1.0,  &(diag64[0]) ); } 
    else if( diag->id() == -65 ) { sel.insert( 1.0,  &(diag65[0]) ); } 
    else if( diag->id() == -66 ) { sel.insert( 1.0,  &(diag66[0]) ); } 
    else if( diag->id() == -67 ) { sel.insert( 1.0,  &(diag67[0]) ); } 
    else if( diag->id() == -68 ) { sel.insert( 1.0,  &(diag68[0]) ); } 
    else if( diag->id() == -69 ) { sel.insert( 1.0,  &(diag69[0]) ); } 
    else if( diag->id() == -70 ) { sel.insert( 1.0,  &(diag70[0]) ); } 
    else if( diag->id() == -71 ) { sel.insert( 1.0,  &(diag71[0]) ); } 
    else if( diag->id() == -72 ) { sel.insert( 1.0,  &(diag72[0]) ); } 
    else if( diag->id() == -73 ) { sel.insert( 1.0,  &(diag73[0]) ); } 
    else if( diag->id() == -74 ) { sel.insert( 1.0,  &(diag74[0]) ); } 
    else if( diag->id() == -75 ) { sel.insert( 1.0,  &(diag75[0]) ); } 
    else if( diag->id() == -76 ) { sel.insert( 1.0,  &(diag76[0]) ); } 
    else if( diag->id() == -77 ) { sel.insert( 1.0,  &(diag77[0]) ); } 
    else if( diag->id() == -78 ) { sel.insert( 1.0,  &(diag78[0]) ); } 
    else if( diag->id() == -79 ) { sel.insert( 1.0,  &(diag79[0]) ); } 
    else if( diag->id() == -80 ) { sel.insert( 1.0,  &(diag80[0]) ); } 
    else if( diag->id() == -81 ) { sel.insert( 1.0,  &(diag81[0]) ); } 
    else if( diag->id() == -82 ) { sel.insert( 1.0,  &(diag82[0]) ); } 
    else if( diag->id() == -83 ) { sel.insert( 1.0,  &(diag83[0]) ); } 
    else if( diag->id() == -84 ) { sel.insert( 1.0,  &(diag84[0]) ); } 
    else if( diag->id() == -85 ) { sel.insert( 1.0,  &(diag85[0]) ); } 
    else if( diag->id() == -86 ) { sel.insert( 1.0,  &(diag86[0]) ); } 
    else if( diag->id() == -87 ) { sel.insert( 1.0,  &(diag87[0]) ); } 
    else if( diag->id() == -88 ) { sel.insert( 1.0,  &(diag88[0]) ); } 
    else if( diag->id() == -89 ) { sel.insert( 1.0,  &(diag89[0]) ); } 
    else if( diag->id() == -90 ) { sel.insert( 1.0,  &(diag90[0]) ); } 
    else if( diag->id() == -91 ) { sel.insert( 1.0,  &(diag91[0]) ); } 
    else if( diag->id() == -92 ) { sel.insert( 1.0,  &(diag92[0]) ); } 
    else if( diag->id() == -93 ) { sel.insert( 1.0,  &(diag93[0]) ); } 
    else if( diag->id() == -94 ) { sel.insert( 1.0,  &(diag94[0]) ); } 
    else if( diag->id() == -95 ) { sel.insert( 1.0,  &(diag95[0]) ); } 
    else if( diag->id() == -96 ) { sel.insert( 1.0,  &(diag96[0]) ); } 
    return sel; 
  }
}

ClassDescription<VBFNLOMEqq2hqqg> VBFNLOMEqq2hqqg::initVBFNLOMEqq2hqqg;

void VBFNLOMEqq2hqqg::persistentOutput(PersistentOStream & ) const {
}

void VBFNLOMEqq2hqqg::persistentInput(PersistentIStream & , int) {
}


void VBFNLOMEqq2hqqg::Init() {

  static ClassDocumentation<VBFNLOMEqq2hqqg> documentation
    ("VBFNLOMEqq2hqqg");

}
