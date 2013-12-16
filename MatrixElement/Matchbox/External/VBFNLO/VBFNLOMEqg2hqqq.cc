
#include "VBFNLOCommonBlocks.h"
#include "VBFNLOMEqg2hqqq.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Throw.h"
 
#include "Diagrams/DC_Hjjj_CC_gqB.h"
#include "Diagrams/DC_Hjjj_CC_gq.h"
#include "Diagrams/DC_Hjjj_CC_qBg.h"
#include "Diagrams/DC_Hjjj_CC_qg.h"
#include "Diagrams/DC_Hjjj_NC_gqB.h"
#include "Diagrams/DC_Hjjj_NC_gq.h"
#include "Diagrams/DC_Hjjj_NC_qBg.h"
#include "Diagrams/DC_Hjjj_NC_qg.h"

#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/DiagramDrawer.h"

using namespace Herwig;


VBFNLOMEqg2hqqq::VBFNLOMEqg2hqqq()
  : VBFNLOMEPP2hJetJetJet(), theWhichGluon(0){}

VBFNLOMEqg2hqqq::VBFNLOMEqg2hqqq(int current, bool incoming1, bool incoming2, int decay, bool nw, int gluon)
  : VBFNLOMEPP2hJetJetJet(current, incoming1, incoming2, decay, nw), theWhichGluon(gluon) {
}
VBFNLOMEqg2hqqq::~VBFNLOMEqg2hqqq(){}

IBPtr VBFNLOMEqg2hqqq::clone() const {
  return new_ptr(*this);
}

IBPtr VBFNLOMEqg2hqqq::fullclone() const {
  return new_ptr(*this);
}

void VBFNLOMEqg2hqqq::prepareMomenta(double (& pbar)[14][4],
					double (&qbar) [5],int gluonIndex) const{
  if (gluonIndex != theWhichGluon) throw Exception() << "The Method prepareMomenta() in matrix element "
						     << name() << " was called with gluonIndex="
						     << gluonIndex << " which is unequal theWhichGluon="
						     << theWhichGluon <<"." << Exception::eventerror;

  if (gluonIndex == 0) {
    L5MomToDouble( (lastMEMomenta()[4]), &pbar[0][0]);

    L5MomToDouble( (lastMEMomenta()[2]), &pbar[1][0]);

    L5MomToDouble( (mePartonData()[1]->id() > 0 ? 
		    lastMEMomenta()[1] : 
		    lastMEMomenta()[3]), &pbar[2][0]);

    L5MomToDouble( (mePartonData()[3]->id() > 0 ? 
		    lastMEMomenta()[3] : 
		    lastMEMomenta()[1]), &pbar[3][0]);
  }
  else if (gluonIndex == 1) {
    L5MomToDouble( (lastMEMomenta()[4]), &pbar[2][0]);

    L5MomToDouble( (lastMEMomenta()[3]), &pbar[3][0]);

    L5MomToDouble( (mePartonData()[0]->id() > 0 ? 
		    lastMEMomenta()[0]  : 
		    lastMEMomenta()[2]), &pbar[0][0]);
    
    L5MomToDouble( (mePartonData()[2]->id() > 0 ? 
		    lastMEMomenta()[2] : 
		    lastMEMomenta()[0]), &pbar[1][0]);
  }

  L5MomToDouble( (lastMEMomenta()[gluonIndex]), &qbar[0]);
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

void VBFNLOMEqg2hqqq::initDiagramContainers() const {
  if ( theCurrent == neutral ){
    if ( theIncoming1 && theWhichGluon == 1 ) theDiagramContainer = RCPtr<DC_Hjjj_NC_qg>::Create(DC_Hjjj_NC_qg(const_cast<VBFNLOMEqg2hqqq*>(this)));
    if ( !theIncoming1 && theWhichGluon == 1 ) theDiagramContainer = RCPtr<DC_Hjjj_NC_qBg>::Create(DC_Hjjj_NC_qBg(const_cast<VBFNLOMEqg2hqqq*>(this)));
    if ( theWhichGluon == 0 && theIncoming2 ) theDiagramContainer = RCPtr<DC_Hjjj_NC_gq>::Create(DC_Hjjj_NC_gq(const_cast<VBFNLOMEqg2hqqq*>(this)));
    if ( theWhichGluon == 0 && !theIncoming2 ) theDiagramContainer = RCPtr<DC_Hjjj_NC_gqB>::Create(DC_Hjjj_NC_gqB(const_cast<VBFNLOMEqg2hqqq*>(this)));
  }
  else if ( theCurrent == charged ){
    if ( theIncoming1 && theWhichGluon == 1 ) theDiagramContainer = RCPtr<DC_Hjjj_CC_qg>::Create(DC_Hjjj_CC_qg(const_cast<VBFNLOMEqg2hqqq*>(this)));
    if ( !theIncoming1 && theWhichGluon == 1 ) theDiagramContainer = RCPtr<DC_Hjjj_CC_qBg>::Create(DC_Hjjj_CC_qBg(const_cast<VBFNLOMEqg2hqqq*>(this)));
    if ( theWhichGluon == 0 && theIncoming2 ) theDiagramContainer = RCPtr<DC_Hjjj_CC_gq>::Create(DC_Hjjj_CC_gq(const_cast<VBFNLOMEqg2hqqq*>(this)));
    if ( theWhichGluon == 0 && !theIncoming2 ) theDiagramContainer = RCPtr<DC_Hjjj_CC_gqB>::Create(DC_Hjjj_CC_gqB(const_cast<VBFNLOMEqg2hqqq*>(this)));
  }
  else
    assert(false);
}

void VBFNLOMEqg2hqqq::getDiagrams() const {
  // cerr << "theCurrent = " << theCurrent
  //      << " theWhichGluon = " << theWhichGluon
  //      << " theIncoming1 = " << theIncoming1
  //      << " theIncoming2 = " << theIncoming2
  //      << "\n" << flush;
  // DiagramDrawer drawer;
  if (theInteger == 0) {
    if (!theDiagramContainer){
      initDiagramContainers();
    }
    allPossibleDiagrams = theDiagramContainer->getDiagrams();
    
    for (vector<DiagPtr>::const_iterator di = allPossibleDiagrams.begin();
	 di != allPossibleDiagrams.end() ; di++) {
      if ( allowedDiagram(*di,theWhichGluon) ) {
	add(*di);
      }
      // cerr << "diagram has id = " << (**di).id() << " and reads: "
      // 	   << (**di).partons()[0]->PDGName()
      // 	   << (**di).partons()[1]->PDGName()
      // 	   << (**di).partons()[2]->PDGName()
      // 	   << (**di).partons()[3]->PDGName()
      // 	   << (**di).partons()[4]->PDGName()
      // 	   << "\n" << flush;
    }
    return;
  }
  
  tcPDPtr Wplus = getParticleData(ParticleID::Wplus);
  tcPDPtr Wminus = getParticleData(ParticleID::Wminus);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  tcPDPtr h0 = getParticleData(ParticleID::h0);
  tcPDPtr g = getParticleData(ParticleID::g);

  bool useOldStyle = false;
  if (useOldStyle) {
    if (NDecayProducts() == 1) {   
      //MOD:
      if ( theCurrent == neutral ) {
	for ( PDVector::const_iterator q1 = theQuarkFlavours.begin();
	      q1 != theQuarkFlavours.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = theQuarkFlavours.begin();
		q2 != theQuarkFlavours.end(); ++q2 ) {
	    //no-emission line is fermionic
	    if (theWhichGluon == 0) {
	      if (requestedAsIncoming2(*q2)) {
		add(new_ptr((Tree2toNDiagram(5), g, *q1, Z0, Z0, *q2,  2, *q1,  4, *q2,  1, (**q1).CC(),  3, h0, -1)));
		add(new_ptr((Tree2toNDiagram(5), g, (**q1).CC(), Z0, Z0, *q2,  1, *q1,  4, *q2,  2, (**q1).CC(),  3, h0, -2)));
	      }
	    }
	    else {
	      if (requestedAsIncoming1(*q1)){
		add(new_ptr((Tree2toNDiagram(5), *q1, Z0, Z0, *q2, g,  1, *q1,  3, *q2,  4, (**q2).CC(),  2, h0, -3)));
		add(new_ptr((Tree2toNDiagram(5), *q1, Z0, Z0, (**q2).CC(), g,  1, *q1,  4, *q2,  3, (**q2).CC(),  2, h0, -4)));
	      }
	    }
	    //no-emission line is antifermionic
	    if (theWhichGluon == 0) {
	      if (requestedAsIncoming2((**q2).CC())){
		add(new_ptr((Tree2toNDiagram(5), g, *q1, Z0, Z0, (**q2).CC(),  2, *q1,  4, (**q2).CC(),  1, (**q1).CC(),  3, h0, -1)));
		add(new_ptr((Tree2toNDiagram(5), g, (**q1).CC(), Z0, Z0, (**q2).CC(),  1, *q1,  4, (**q2).CC(),  2, (**q1).CC(),   3, h0, -2)));
	      }
	    }
	    else {
	      if (requestedAsIncoming1((**q1).CC())){
		add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Z0, Z0, *q2, g,  1, (**q1).CC(),  3, *q2,  4, (**q2).CC(),  2, h0, -3)));
		add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Z0, Z0, (**q2).CC(), g,  1, (**q1).CC(),  4, *q2,  3, (**q2).CC(),  2, h0, -4)));
	      }
	    }
	  }
	}
      }
 
      else if ( theCurrent == charged ) {
	PDVector QuarksAndAntiQuarks;
	for ( PDVector::const_iterator q1 = theQuarkFlavours.begin();
	      q1 != theQuarkFlavours.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = theQuarkFlavours.begin();
		q2 != theQuarkFlavours.end(); ++q2 ) {

	    tcPDPtr q1prime = SU2Helper::SU2CC(*q1);
	    tcPDPtr q2prime = SU2Helper::SU2CC(*q2);
	    tcPDPtr qb1prime = SU2Helper::SU2CC((**q1).CC());
	    tcPDPtr qb2prime = SU2Helper::SU2CC((**q2).CC());

	    if (abs(( (-(**q2).charge() - (**q1).charge() + (*q1prime).charge() + (*q2prime).charge()))/eplus) < 0.00001){
	      //then no-emission line is lower line and fermionic
	      Charge transfer = ZERO;
	      if (theWhichGluon == 0){
		transfer = -(**q1).charge() + (*q1prime).charge();
		if (transfer < ZERO) {
		  if (requestedAsIncoming2(*q2)){
		    add(new_ptr((Tree2toNDiagram(5), g, *q1, Wplus, Wplus, *q2,  2, q1prime,  4, q2prime,  1, (**q1).CC(),  3, h0, -1)));
		    add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wplus, Wplus, *q2,  1, q1prime,  4, q2prime,  2, (**q1).CC(),  3, h0, -2)));
		  }
		}
		else if (transfer > ZERO) {
		  if (requestedAsIncoming2(*q2)){
		    add(new_ptr((Tree2toNDiagram(5), g, *q1, Wminus, Wminus, *q2,  2, q1prime,  4, q2prime,  1, (**q1).CC(),  3, h0, -1)));
		    add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wminus, Wminus, *q2,  1, q1prime,  4, q2prime,  2, (**q1).CC(),  3, h0, -2))); 
		  }
		}
	      }
	      else {
		//or no-emission line is upper line and fermionic
		transfer = -(**q2).charge() + (*q2prime).charge();
		if (transfer > ZERO) {
		  if (requestedAsIncoming1(*q1)){
		    add(new_ptr((Tree2toNDiagram(5), *q1, Wplus, Wplus, *q2, g,  1, q1prime,  3, q2prime,  4, (**q2).CC(),  2, h0, -3)));
		    add(new_ptr((Tree2toNDiagram(5), *q1, Wplus, Wplus, qb2prime, g,  1, q1prime,  4, q2prime,  3, (**q2).CC(),  2, h0, -4)));
		  }
		}
		else if (transfer < ZERO) {
		  if (requestedAsIncoming1(*q1)){
		    add(new_ptr((Tree2toNDiagram(5), *q1, Wminus, Wminus, *q2, g,  1, q1prime,  3, q2prime,  4, (**q2).CC(),  2, h0, -3)));
		    add(new_ptr((Tree2toNDiagram(5), *q1, Wminus, Wminus, qb2prime, g,  1, q1prime,  4, q2prime,  3, (**q2).CC(),  2, h0, -4)));
		  }
		}
	      }
	    }
	    else if (abs(( (-(**q2).CC()->charge() - (**q1).charge() + (*q1prime).charge() + (*qb2prime).charge()))/eplus) < 0.00001){
	      Charge transfer = ZERO;
	      if (theWhichGluon == 0) {
		//then no-emission line is lower line and antifermionic
		transfer = -(**q1).charge() + (*q1prime).charge();
		if (transfer < ZERO) {
		  if (requestedAsIncoming2((**q2).CC())){
		    add(new_ptr((Tree2toNDiagram(5), g, *q1, Wplus, Wplus, (**q2).CC(),  2, q1prime,  4, qb2prime,  1, (**q1).CC(),  3, h0, -1)));
		    add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wplus, Wplus, (**q2).CC(),  1, q1prime,  4, qb2prime,  2, (**q1).CC(),  3, h0, -2)));
		  }
		}
		else if (transfer > ZERO) {
		  if (requestedAsIncoming2((**q2).CC())){
		    add(new_ptr((Tree2toNDiagram(5), g, *q1, Wminus, Wminus, (**q2).CC(),  2, q1prime,  4, qb2prime,  1, (**q1).CC(),  3, h0, -1)));
		    add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wminus, Wminus, (**q2).CC(),  1, q1prime,  4, qb2prime,  2, (**q1).CC(),  3, h0, -2)));
		  }
		}
	      }
	      else {
		//or no-emission line is upper line and antifermionic
		transfer = -(**q2).charge() + (*q2prime).charge();
		if (transfer > ZERO) {
		  if (requestedAsIncoming1((**q1).CC())){
		    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wplus, Wplus, *q2, g,  1, qb1prime,  3, q2prime,  4, (**q2).CC(),  2, h0, -3)));
		    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wplus, Wplus, qb2prime, g,  1, qb1prime,  4, q2prime,  3, (**q2).CC(),  2, h0, -4)));
		  }
		}
		else if (transfer < ZERO) {
		  if (requestedAsIncoming1((**q1).CC())){
		    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wminus, Wminus, *q2, g,  1, qb1prime,  3, q2prime,  4, (**q2).CC(),  2, h0, -3)));
		    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wminus, Wminus, qb2prime, g,  1, qb1prime,  4, q2prime,  3, (**q2).CC(),  2, h0, -4)));
		  }
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
      if ( theCurrent == neutral ) {
	for ( PDVector::const_iterator q1 = theQuarkFlavours.begin();
	      q1 != theQuarkFlavours.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = theQuarkFlavours.begin();
		q2 != theQuarkFlavours.end(); ++q2 ) {
	    //no-emission line is fermionic
	    if (theWhichGluon == 0) {
	      if (requestedAsIncoming2(*q2)){
		add(new_ptr((Tree2toNDiagram(5), g, *q1, Z0, Z0, *q2,  2, *q1,  4, *q2,  1, (**q1).CC(),  3, h0,  9,dec1,  9,dec2, -1)));
		add(new_ptr((Tree2toNDiagram(5), g, (**q1).CC(), Z0, Z0, *q2,  1, *q1,  4, *q2,  2, (**q1).CC(),  3, h0,  9,dec1,  9,dec2, -2)));
	      }
	    }
	    else {
	      if (requestedAsIncoming1(*q1)){
		add(new_ptr((Tree2toNDiagram(5), *q1, Z0, Z0, *q2, g,  1, *q1,  3, *q2,  4, (**q2).CC(),  2, h0,  9,dec1,  9,dec2, -3)));
		add(new_ptr((Tree2toNDiagram(5), *q1, Z0, Z0, (**q2).CC(), g,  1, *q1,  4, *q2,  3, (**q2).CC(),  2, h0,  9,dec1,  9,dec2, -4)));
	      }
	    }
	    //no-emission line is antifermionic
	    if (theWhichGluon == 0) {
	      if (requestedAsIncoming2((**q2).CC())){
		add(new_ptr((Tree2toNDiagram(5), g, *q1, Z0, Z0, (**q2).CC(),  2, *q1,  4, (**q2).CC(),  1, (**q1).CC(),  3, h0,  9,dec1,  9,dec2, -1)));
		add(new_ptr((Tree2toNDiagram(5), g, (**q1).CC(), Z0, Z0, (**q2).CC(),  1, *q1,  4, (**q2).CC(),  2, (**q1).CC(),  3, h0,  9,dec1,  9,dec2, -2)));
	      }
	    }
	    else {
	      if (requestedAsIncoming1((**q1).CC())){
		add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Z0, Z0, *q2, g,  1, (**q1).CC(),  3, *q2,  4, (**q2).CC(),  2, h0,  9,dec1,  9,dec2, -3)));
		add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Z0, Z0, (**q2).CC(), g,  1, (**q1).CC(),  4, *q2,  3, (**q2).CC(),  2, h0,  9,dec1,  9,dec2, -4)));
	      }
	    }
	  }
	}
      }

      else if ( theCurrent == charged ) {
	PDVector QuarksAndAntiQuarks;
	for ( PDVector::const_iterator q1 = theQuarkFlavours.begin();
	      q1 != theQuarkFlavours.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = theQuarkFlavours.begin();
		q2 != theQuarkFlavours.end(); ++q2 ) {

	    tcPDPtr q1prime = SU2Helper::SU2CC(*q1);
	    tcPDPtr q2prime = SU2Helper::SU2CC(*q2);
	    tcPDPtr qb1prime = SU2Helper::SU2CC((**q1).CC());
	    tcPDPtr qb2prime = SU2Helper::SU2CC((**q2).CC());

	    if (abs(( (-(**q2).charge() - (**q1).charge() + (*q1prime).charge() + (*q2prime).charge()))/eplus) < 0.00001){
	      //then no-emission line is lower line and fermionic
	      Charge transfer = ZERO;
	      if (theWhichGluon == 0){
		transfer = -(**q1).charge() + (*q1prime).charge();
		if (transfer < ZERO) {
		  if (requestedAsIncoming2(*q2)){
		    add(new_ptr((Tree2toNDiagram(5), g, *q1, Wplus, Wplus, *q2,  2, q1prime,  4, q2prime,  1, (**q1).CC(),  3, h0,  9,dec1,  9,dec2, -1)));
		    add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wplus, Wplus, *q2,  1, q1prime,  4, q2prime,  2, (**q1).CC(),  3, h0,  9,dec1,  9,dec2, -2)));
		  }
		}
		else if (transfer > ZERO) {
		  if (requestedAsIncoming2(*q2)){
		    add(new_ptr((Tree2toNDiagram(5), g, *q1, Wminus, Wminus, *q2,  2, q1prime,  4, q2prime,  1, (**q1).CC(),  3, h0,  9,dec1,  9,dec2, -1)));
		    add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wminus, Wminus, *q2,  1, q1prime,  4, q2prime,  2, (**q1).CC(),  3, h0,  9,dec1,  9,dec2, -2)));
		  }
		}
	      }
	      else {
		//or no-emission line is upper line and fermionic
		transfer = -(**q2).charge() + (*q2prime).charge();
		if (transfer > ZERO) {
		  if (requestedAsIncoming1(*q1)){
		    add(new_ptr((Tree2toNDiagram(5), *q1, Wplus, Wplus, *q2, g,  1, q1prime,  3, q2prime,  4, (**q2).CC(),  2, h0,  9,dec1,  9,dec2, -3)));
		    add(new_ptr((Tree2toNDiagram(5), *q1, Wplus, Wplus, qb2prime, g,  1, q1prime,  4, q2prime,  3, (**q2).CC(),  2, h0,  9,dec1,  9,dec2, -4)));
		  }
		}
		else if (transfer < ZERO) {
		  if (requestedAsIncoming1(*q1)){
		    add(new_ptr((Tree2toNDiagram(5), *q1, Wminus, Wminus, *q2, g,  1, q1prime,  3, q2prime,  4, (**q2).CC(),  2, h0,  9,dec1,  9,dec2, -3)));
		    add(new_ptr((Tree2toNDiagram(5), *q1, Wminus, Wminus, qb2prime, g,  1, q1prime,  4, q2prime,  3, (**q2).CC(),  2, h0,  9,dec1,  9,dec2, -4)));
		  }
		}
	      }
	    }
	    else if (abs(( (-(**q2).CC()->charge() - (**q1).charge() + (*q1prime).charge() + (*qb2prime).charge()))/eplus) < 0.00001){
	      Charge transfer = ZERO;
	      if (theWhichGluon == 0) {
		//then no-emission line is lower line and antifermionic
		transfer = -(**q1).charge() + (*q1prime).charge();
		if (transfer < ZERO) {
		  if (requestedAsIncoming2((**q2).CC())){
		    add(new_ptr((Tree2toNDiagram(5), g, *q1, Wplus, Wplus, (**q2).CC(),  2, q1prime,  4, qb2prime,  1, (**q1).CC(),  3, h0,  9,dec1,  9,dec2, -1)));
		    add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wplus, Wplus, (**q2).CC(),  1, q1prime,  4, qb2prime,  2, (**q1).CC(),  1, q1prime,  3, h0,  9,dec1,  9,dec2, -2)));
		  }
		}
		else if (transfer > ZERO) {
		  if (requestedAsIncoming2((**q2).CC())){
		    add(new_ptr((Tree2toNDiagram(5), g, *q1, Wminus, Wminus, (**q2).CC(),  2, q1prime,  4, qb2prime,  1, (**q1).CC(),  3, h0,  9,dec1,  9,dec2, -1)));
		    add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wminus, Wminus, (**q2).CC(),  1, q1prime,  4, qb2prime,  2, (**q1).CC(),  3, h0,  9,dec1,  9,dec2, -2)));
		  }
		}
	      }
	      else {
		//or no-emission line is upper line and antifermionic
		transfer = -(**q2).charge() + (*q2prime).charge();
		if (transfer > ZERO) {
		  if (requestedAsIncoming1((**q1).CC())){
		    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wplus, Wplus, *q2, g,  1, qb1prime,  3, q2prime,  4, (**q2).CC(),  2, h0,  9,dec1,  9,dec2, -3)));
		    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wplus, Wplus, qb2prime, g,  1, qb1prime,  4, q2prime,  3, (**q2).CC(),  2, h0,  9,dec1,  9,dec2, -4)));
		  }
		}
		else if (transfer < ZERO) {
		  if (requestedAsIncoming1((**q1).CC())){
		    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wminus, Wminus, *q2, g,  1, qb1prime,  3, q2prime,  4, (**q2).CC(),  2, h0,  9,dec1,  9,dec2, -3)));
		    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wminus, Wminus, qb2prime, g,  1, qb1prime,  4, q2prime,  3, (**q2).CC(),  2, h0,  9,dec1,  9,dec2, -4)));
		  }
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
    else if (NDecayProducts() == 4){
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
	for ( PDVector::const_iterator q1 = theQuarkFlavours.begin();
	      q1 != theQuarkFlavours.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = theQuarkFlavours.begin();
		q2 != theQuarkFlavours.end(); ++q2 ) {
	    //no-emission line is fermionic
	    if (theWhichGluon == 0) {
	      if (requestedAsIncoming2(*q2)){
		add(new_ptr((Tree2toNDiagram(5), g, *q1, Z0, Z0, *q2,  2, *q1,  4, *q2,  1, (**q1).CC(),  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -1)));
		add(new_ptr((Tree2toNDiagram(5), g, (**q1).CC(), Z0, Z0, *q2,  1, *q1,  4, *q2,  2, (**q1).CC(),  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -2)));
	      }
	    }
	    else {
	      if (requestedAsIncoming1(*q1)){
		add(new_ptr((Tree2toNDiagram(5), *q1, Z0, Z0, *q2, g,  1, *q1,  3, *q2,  4, (**q2).CC(),  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -3)));
		add(new_ptr((Tree2toNDiagram(5), *q1, Z0, Z0, (**q2).CC(), g,  1, *q1,  4, *q2,  3, (**q2).CC(),  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -4)));
	      }
	    }
	    //no-emission line is antifermionic
	    if (theWhichGluon == 0) {
	      if (requestedAsIncoming2((**q2).CC())){
		add(new_ptr((Tree2toNDiagram(5), g, *q1, Z0, Z0, (**q2).CC(),  2, *q1,  4, (**q2).CC(),  1, (**q1).CC(),  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -1)));
		add(new_ptr((Tree2toNDiagram(5), g, (**q1).CC(), Z0, Z0, (**q2).CC(),  1, *q1,  4, (**q2).CC(),  2, (**q1).CC(),  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -2)));
	      }
	    }
	    else {
	      if (requestedAsIncoming1((**q1).CC())){
		add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Z0, Z0, *q2, g,  1, (**q1).CC(),  3, *q2,  4, (**q2).CC(),  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -3)));
		add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Z0, Z0, (**q2).CC(), g,  1, (**q1).CC(),  4, *q2,  3, (**q2).CC(),  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -4)));
	      }
	    }
	  }
	}
      }

      else if ( theCurrent == charged ) {
	PDVector QuarksAndAntiQuarks;
	for ( PDVector::const_iterator q1 = theQuarkFlavours.begin();
	      q1 != theQuarkFlavours.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = theQuarkFlavours.begin();
		q2 != theQuarkFlavours.end(); ++q2 ) {

	    tcPDPtr q1prime = SU2Helper::SU2CC(*q1);
	    tcPDPtr q2prime = SU2Helper::SU2CC(*q2);
	    tcPDPtr qb1prime = SU2Helper::SU2CC((**q1).CC());
	    tcPDPtr qb2prime = SU2Helper::SU2CC((**q2).CC());

	    if (abs(( (-(**q2).charge() - (**q1).charge() + (*q1prime).charge() + (*q2prime).charge()))/eplus) < 0.00001){
	      //then no-emission line is lower line and fermionic
	      Charge transfer = ZERO;
	      if (theWhichGluon == 0){
		transfer = -(**q1).charge() + (*q1prime).charge();
		if (transfer < ZERO) {
		  if (requestedAsIncoming2(*q2)){
		    add(new_ptr((Tree2toNDiagram(5), g, *q1, Wplus, Wplus, *q2,  2, q1prime,  4, q2prime,  1, (**q1).CC(),  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -1)));
		    add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wplus, Wplus, *q2,  1, q1prime,  4, q2prime,  2, (**q1).CC(),  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -2)));
		  }
		}
		else if (transfer > ZERO) {
		  if (requestedAsIncoming2(*q2)){
		    add(new_ptr((Tree2toNDiagram(5), g, *q1, Wminus, Wminus, *q2,  2, q1prime,  4, q2prime,  1, (**q1).CC(),  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -1)));
		    add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wminus, Wminus, *q2,  1, q1prime,  4, q2prime,  2, (**q1).CC(),  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -2)));
		  } 
		}
	      }
	      else {
		//or no-emission line is upper line and fermionic
		transfer = -(**q2).charge() + (*q2prime).charge();
		if (transfer > ZERO) {
		  if (requestedAsIncoming1(*q1)){
		    add(new_ptr((Tree2toNDiagram(5), *q1, Wplus, Wplus, *q2, g,  1, q1prime,  3, q2prime,  4, (**q2).CC(),  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -3)));
		    add(new_ptr((Tree2toNDiagram(5), *q1, Wplus, Wplus, qb2prime, g,  1, q1prime,  4, q2prime,  3, (**q2).CC(),  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -4)));
		  }
		}
		else if (transfer < ZERO) {
		  if (requestedAsIncoming1(*q1)){
		    add(new_ptr((Tree2toNDiagram(5), *q1, Wminus, Wminus, *q2, g,  1, q1prime,  3, q2prime,  4, (**q2).CC(),  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -3)));
		    add(new_ptr((Tree2toNDiagram(5), *q1, Wminus, Wminus, qb2prime, g,  1, q1prime,  4, q2prime,  3, (**q2).CC(),  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -4)));
		  }
		}
	      }
	    }
	    else if (abs(( (-(**q2).CC()->charge() - (**q1).charge() + (*q1prime).charge() + (*qb2prime).charge()))/eplus) < 0.00001){
	      Charge transfer = ZERO;
	      if (theWhichGluon == 0) {
		//then no-emission line is lower line and antifermionic
		transfer = -(**q1).charge() + (*q1prime).charge();
		if (transfer < ZERO) {
		  if (requestedAsIncoming2((**q2).CC())){
		    add(new_ptr((Tree2toNDiagram(5), g, *q1, Wplus, Wplus, (**q2).CC(),  2, q1prime,  4, qb2prime,  1, (**q1).CC(),  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -1)));
		    add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wplus, Wplus,  4, qb2prime,  2, (**q1).CC(), (**q2).CC(),  1, q1prime,  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -2)));
		  }
		}
		else if (transfer > ZERO) {
		  if (requestedAsIncoming2((**q2).CC())){
		    add(new_ptr((Tree2toNDiagram(5), g, *q1, Wminus, Wminus, (**q2).CC(),  2, q1prime,  4, qb2prime,  1, (**q1).CC(),  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -1)));
		    add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wminus, Wminus, (**q2).CC(),  1, q1prime,  4, qb2prime,  2, (**q1).CC(),  3, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -2)));
		  }
		}
	      }
	      else {
		//or no-emission line is upper line and antifermionic
		transfer = -(**q2).charge() + (*q2prime).charge();
		if (transfer > ZERO) {
		  if (requestedAsIncoming1((**q1).CC())){
		    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wplus, Wplus, *q2, g,  1, qb1prime,  3, q2prime,  4, (**q2).CC(),  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -3)));
		    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wplus, Wplus, qb2prime, g,  1, qb1prime,  4, q2prime,  3, (**q2).CC(),  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -4)));
		  }
		}
		else if (transfer < ZERO) {
		  if (requestedAsIncoming1((**q1).CC())){
		    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wminus, Wminus, *q2, g,  1, qb1prime,  3, q2prime,  4, (**q2).CC(),  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -3)));
		    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wminus, Wminus, qb2prime, g,  1, qb1prime,  4, q2prime,  3, (**q2).CC(),  2, h0,  9,VecBos1,  9,VecBos2,  10,dec1,  10,dec2,  11,dec3,  11,dec4, -4)));
		  }
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

  else { //!useOldStyle
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

	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, quBC, Wminus, Wminus, kd, 1, quC, 4, kdC, 2, quB, 3, h0, -1)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qu, Wminus, Wminus, kd, 2, quC, 4, kdC, 1, quB, 3, h0, -2)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qdBC, Wplus, Wplus, kdB, 1, qdC, 4, kdBC, 2, qdB, 3, h0, -3)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qd, Wplus, Wplus, kdB, 2, qdC, 4, kdBC, 1, qdB, 3, h0, -4)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qdBC, Wplus, Wplus, ku, 1, qdC, 4, kuC, 2, qdB, 3, h0, -5)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qd, Wplus, Wplus, ku, 2, qdC, 4, kuC, 1, qdB, 3, h0, -6)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, quBC, Wminus, Wminus, kuB, 1, quC, 4, kuBC, 2, quB, 3, h0, -7)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qu, Wminus, Wminus, kuB, 2, quC, 4, kuBC, 1, quB, 3, h0, -8)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, Wminus, Wminus, kdC, g, 1, qdBC, 4, kdC, 3, kdB, 2, h0, -9)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, Wminus, Wminus, kdB, g, 1, qdBC, 3, kdC, 4, kdB, 2, h0, -10)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Wplus, Wplus, kuC, g, 1, qdC, 4, kuC, 3, kuB, 2, h0, -11)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Wplus, Wplus, kuB, g, 1, qdC, 3, kuC, 4, kuB, 2, h0, -12)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, Wplus, Wplus, kuC, g, 1, quBC, 4, kuC, 3, kuB, 2, h0, -13)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, Wplus, Wplus, kuB, g, 1, quBC, 3, kuC, 4, kuB, 2, h0, -14)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Wminus, Wminus, kdC, g, 1, quC, 4, kdC, 3, kdB, 2, h0, -15)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Wminus, Wminus, kdB, g, 1, quC, 3, kdC, 4, kdB, 2, h0, -16)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qdB, Z0, Z0, kd, 1, qd, 4, kd, 2, qdB, 3, h0, -17)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qd, Z0, Z0, kd, 2, qd, 4, kd, 1, qdB, 3, h0, -18)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, quB, Z0, Z0, kd, 1, qu, 4, kd, 2, quB, 3, h0, -19)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qu, Z0, Z0, kd, 2, qu, 4, kd, 1, quB, 3, h0, -20)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qdB, Z0, Z0, kdB, 1, qd, 4, kdB, 2, qdB, 3, h0, -21)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qd, Z0, Z0, kdB, 2, qd, 4, kdB, 1, qdB, 3, h0, -22)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, quB, Z0, Z0, kdB, 1, qu, 4, kdB, 2, quB, 3, h0, -23)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qu, Z0, Z0, kdB, 2, qu, 4, kdB, 1, quB, 3, h0, -24)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qdB, Z0, Z0, ku, 1, qd, 4, ku, 2, qdB, 3, h0, -25)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qd, Z0, Z0, ku, 2, qd, 4, ku, 1, qdB, 3, h0, -26)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, quB, Z0, Z0, ku, 1, qu, 4, ku, 2, quB, 3, h0, -27)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qu, Z0, Z0, ku, 2, qu, 4, ku, 1, quB, 3, h0, -28)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qdB, Z0, Z0, kuB, 1, qd, 4, kuB, 2, qdB, 3, h0, -29)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qd, Z0, Z0, kuB, 2, qd, 4, kuB, 1, qdB, 3, h0, -30)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, quB, Z0, Z0, kuB, 1, qu, 4, kuB, 2, quB, 3, h0, -31)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qu, Z0, Z0, kuB, 2, qu, 4, kuB, 1, quB, 3, h0, -32)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, Z0, Z0, kd, g, 1, qdB, 4, kd, 3, kdB, 2, h0, -33)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, Z0, Z0, kdB, g, 1, qdB, 3, kd, 4, kdB, 2, h0, -34)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, Z0, Z0, ku, g, 1, qdB, 4, ku, 3, kuB, 2, h0, -35)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qdB, Z0, Z0, kuB, g, 1, qdB, 3, ku, 4, kuB, 2, h0, -36)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Z0, Z0, kd, g, 1, qd, 4, kd, 3, kdB, 2, h0, -37)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Z0, Z0, kdB, g, 1, qd, 3, kd, 4, kdB, 2, h0, -38)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Z0, Z0, ku, g, 1, qd, 4, ku, 3, kuB, 2, h0, -39)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Z0, Z0, kuB, g, 1, qd, 3, ku, 4, kuB, 2, h0, -40)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, Z0, Z0, kd, g, 1, quB, 4, kd, 3, kdB, 2, h0, -41)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, Z0, Z0, kdB, g, 1, quB, 3, kd, 4, kdB, 2, h0, -42)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, Z0, Z0, ku, g, 1, quB, 4, ku, 3, kuB, 2, h0, -43)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), quB, Z0, Z0, kuB, g, 1, quB, 3, ku, 4, kuB, 2, h0, -44)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Z0, Z0, kd, g, 1, qu, 4, kd, 3, kdB, 2, h0, -45)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Z0, Z0, kdB, g, 1, qu, 3, kd, 4, kdB, 2, h0, -46)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Z0, Z0, ku, g, 1, qu, 4, ku, 3, kuB, 2, h0, -47)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Z0, Z0, kuB, g, 1, qu, 3, ku, 4, kuB, 2, h0, -48)));
      }
    }

    //now check which of those are really requested by the interfaced parameters
    for (vector<DiagPtr>::const_iterator di = allPossibleDiagrams.begin();
	 di != allPossibleDiagrams.end() ; di++) {
      if ( allowedDiagram(*di,theWhichGluon) ) {
	add(*di);
      }
    }
  }
}

void VBFNLOMEqg2hqqq::doinit(){
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
  return;
}

void VBFNLOMEqg2hqqq::doinitrun(){ 
  VBFNLOMEVVJJNeutralBase::doinitrun();
}


Selector<MEBase::DiagramIndex> 
VBFNLOMEqg2hqqq::diagrams(const DiagramVector & diags) const {
  Selector<MEBase::DiagramIndex> sel;
  double w0 = 0, w1 = 0;
  if (mePartonData()[0]->id() == 21) {
    w0 = sqr( sqr(generator()->maximumCMEnergy())/(2.*(meMomenta()[0]*meMomenta()[4])) );
    w1 = sqr( sqr(generator()->maximumCMEnergy())/(2.*(meMomenta()[0]*meMomenta()[2])) );
  }
  else if (mePartonData()[1]->id() == 21){
    w0 = sqr( sqr(generator()->maximumCMEnergy())/(2.*(meMomenta()[1]*meMomenta()[4])) );
    w1 = sqr( sqr(generator()->maximumCMEnergy())/(2.*(meMomenta()[1]*meMomenta()[3])) );
  }
  else assert(false);
  assert(diags.size() == 2);
  sel.insert(w0,0);
  sel.insert(w1,1);
  return sel;
}

Selector<const ColourLines *>
VBFNLOMEqg2hqqq::colourGeometries(tcDiagPtr diag) const {
  if (theInteger == 0){
    if (!theDiagramContainer){
      initDiagramContainers();
    }
    return theDiagramContainer->colourGeometries(diag);
  }
  bool useOldStyle = false;
  if (useOldStyle){
    if (DecayChannel() != 4) {
      static const ColourLines cgF[4] =    {ColourLines("-1 -8,  1  2  6,  5  7"),
					    ColourLines("-1 -2 -8,  1  6,  5  7"),
					    ColourLines(" 1  6, -5 -8,  5  4  7"),
					    ColourLines(" 1  6, -5 -4 -8,  5  7")};
      static const ColourLines cgFbar[4] = {ColourLines("-1 -8,  1  2  6, -5 -7"),
					    ColourLines("-1 -2 -8,  1  6, -5 -7"),
					    ColourLines("-1 -6, -5 -8,  5  4  7"),
					    ColourLines("-1 -6, -5 -4 -8,  5  7")};
      Selector<const ColourLines *> sel;
      //!!!cleanup after checks: check of diagram ids in if conditions may be removed
      if      ( mePartonData()[0]->id() == 21 && mePartonData()[1]->id() > 0 && (abs(diag->id()) == 1 || abs(diag->id()) == 2)){
	sel.insert(1.0, &cgF[abs(diag->id())-1]);
      }
      else if ( mePartonData()[0]->id() == 21 && mePartonData()[1]->id() < 0 && (abs(diag->id()) == 1 || abs(diag->id()) == 2)){
	sel.insert(1.0, &cgFbar[abs(diag->id())-1]);
      }
      else if ( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() == 21 && (abs(diag->id()) == 3 || abs(diag->id()) == 4)){
	sel.insert(1.0, &cgF[abs(diag->id())-1]);
      }
      else if ( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() == 21 && (abs(diag->id()) == 3 || abs(diag->id()) == 4)){
	sel.insert(1.0, &cgFbar[abs(diag->id())-1]);
      }
      else cout << "something wrong in colour line selector"<<endl;
      //  cout<<"sel="<<sel[0]<<endl;
      return sel;
    }
    else{
      
      static const ColourLines cgF[4] =    {ColourLines("-1 -8,  1  2  6,  5  7,  10 -11"),
					    ColourLines("-1 -2 -8,  1  6,  5  7,  10 -11"),
					    ColourLines(" 1  6, -5 -8,  5  4  7,  10 -11"),
					    ColourLines(" 1  6, -5 -4 -8,  5  7,  10 -11")};
      static const ColourLines cgFbar[4] = {ColourLines("-1 -8,  1  2  6, -5 -7,  10 -11"),
					    ColourLines("-1 -2 -8,  1  6, -5 -7,  10 -11"),
					    ColourLines("-1 -6, -5 -8,  5  4  7,  10 -11"),
					    ColourLines("-1 -6, -5 -4 -8,  5  7,  10 -11")};
      Selector<const ColourLines *> sel;
      //!!!cleanup after checks: check of diagram ids in if conditions may be removed
      if      ( mePartonData()[0]->id() == 21 && mePartonData()[1]->id() > 0 && (abs(diag->id()) == 1 || abs(diag->id()) == 2)){
	sel.insert(1.0, &cgF[abs(diag->id())-1]);
      }
      else if ( mePartonData()[0]->id() == 21 && mePartonData()[1]->id() < 0 && (abs(diag->id()) == 1 || abs(diag->id()) == 2)){
	sel.insert(1.0, &cgFbar[abs(diag->id())-1]);
      }
      else if ( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() == 21 && (abs(diag->id()) == 3 || abs(diag->id()) == 4)){
	sel.insert(1.0, &cgF[abs(diag->id())-1]);
      }
      else if ( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() == 21 && (abs(diag->id()) == 3 || abs(diag->id()) == 4)){
	sel.insert(1.0, &cgFbar[abs(diag->id())-1]);
      }
      else cout << "something wrong in colour line selector"<<endl;
      return sel;
    }
  }
  else{
  static const ColourLines diag1[1] = { 
    ColourLines("1 6, -1 -2 -8, 5 7")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("1 2 6, -1 -8, 5 7")
  }; 
  static const ColourLines diag3[1] = { 
    ColourLines("1 6, -1 -2 -8, -5 -7")
  }; 
  static const ColourLines diag4[1] = { 
    ColourLines("1 2 6, -1 -8, -5 -7")
  }; 
  static const ColourLines diag5[1] = { 
    ColourLines("1 6, -1 -2 -8, 5 7")
  }; 
  static const ColourLines diag6[1] = { 
    ColourLines("1 2 6, -1 -8, 5 7")
  }; 
  static const ColourLines diag7[1] = { 
    ColourLines("1 6, -1 -2 -8, -5 -7")
  }; 
  static const ColourLines diag8[1] = { 
    ColourLines("1 2 6, -1 -8, -5 -7")
  }; 
  static const ColourLines diag9[1] = { 
    ColourLines("-1 -6, 5 7, -5 4 -8")
  }; 
  static const ColourLines diag10[1] = { 
    ColourLines("-1 -6, 5 -4 7, -5 -8")
  }; 
  static const ColourLines diag11[1] = { 
    ColourLines("1 6, 5 7, -5 4 -8")
  }; 
  static const ColourLines diag12[1] = { 
    ColourLines("1 6, 5 -4 7, -5 -8")
  }; 
  static const ColourLines diag13[1] = { 
    ColourLines("-1 -6, 5 7, -5 4 -8")
  }; 
  static const ColourLines diag14[1] = { 
    ColourLines("-1 -6, 5 -4 7, -5 -8")
  }; 
  static const ColourLines diag15[1] = { 
    ColourLines("1 6, 5 7, -5 4 -8")
  }; 
  static const ColourLines diag16[1] = { 
    ColourLines("1 6, 5 -4 7, -5 -8")
  }; 
  static const ColourLines diag17[1] = { 
    ColourLines("1 6, -1 -2 -8, 5 7")
  }; 
  static const ColourLines diag18[1] = { 
    ColourLines("1 2 6, -1 -8, 5 7")
  }; 
  static const ColourLines diag19[1] = { 
    ColourLines("1 6, -1 -2 -8, 5 7")
  }; 
  static const ColourLines diag20[1] = { 
    ColourLines("1 2 6, -1 -8, 5 7")
  }; 
  static const ColourLines diag21[1] = { 
    ColourLines("1 6, -1 -2 -8, -5 -7")
  }; 
  static const ColourLines diag22[1] = { 
    ColourLines("1 2 6, -1 -8, -5 -7")
  }; 
  static const ColourLines diag23[1] = { 
    ColourLines("1 6, -1 -2 -8, -5 -7")
  }; 
  static const ColourLines diag24[1] = { 
    ColourLines("1 2 6, -1 -8, -5 -7")
  }; 
  static const ColourLines diag25[1] = { 
    ColourLines("1 6, -1 -2 -8, 5 7")
  }; 
  static const ColourLines diag26[1] = { 
    ColourLines("1 2 6, -1 -8, 5 7")
  }; 
  static const ColourLines diag27[1] = { 
    ColourLines("1 6, -1 -2 -8, 5 7")
  }; 
  static const ColourLines diag28[1] = { 
    ColourLines("1 2 6, -1 -8, 5 7")
  }; 
  static const ColourLines diag29[1] = { 
    ColourLines("1 6, -1 -2 -8, -5 -7")
  }; 
  static const ColourLines diag30[1] = { 
    ColourLines("1 2 6, -1 -8, -5 -7")
  }; 
  static const ColourLines diag31[1] = { 
    ColourLines("1 6, -1 -2 -8, -5 -7")
  }; 
  static const ColourLines diag32[1] = { 
    ColourLines("1 2 6, -1 -8, -5 -7")
  }; 
  static const ColourLines diag33[1] = { 
    ColourLines("-1 -6, 5 7, -5 4 -8")
  }; 
  static const ColourLines diag34[1] = { 
    ColourLines("-1 -6, 5 -4 7, -5 -8")
  }; 
  static const ColourLines diag35[1] = { 
    ColourLines("-1 -6, 5 7, -5 4 -8")
  }; 
  static const ColourLines diag36[1] = { 
    ColourLines("-1 -6, 5 -4 7, -5 -8")
  }; 
  static const ColourLines diag37[1] = { 
    ColourLines("1 6, 5 7, -5 4 -8")
  }; 
  static const ColourLines diag38[1] = { 
    ColourLines("1 6, 5 -4 7, -5 -8")
  }; 
  static const ColourLines diag39[1] = { 
    ColourLines("1 6, 5 7, -5 4 -8")
  }; 
  static const ColourLines diag40[1] = { 
    ColourLines("1 6, 5 -4 7, -5 -8")
  }; 
  static const ColourLines diag41[1] = { 
    ColourLines("-1 -6, 5 7, -5 4 -8")
  }; 
  static const ColourLines diag42[1] = { 
    ColourLines("-1 -6, 5 -4 7, -5 -8")
  }; 
  static const ColourLines diag43[1] = { 
    ColourLines("-1 -6, 5 7, -5 4 -8")
  }; 
  static const ColourLines diag44[1] = { 
    ColourLines("-1 -6, 5 -4 7, -5 -8")
  }; 
  static const ColourLines diag45[1] = { 
    ColourLines("1 6, 5 7, -5 4 -8")
  }; 
  static const ColourLines diag46[1] = { 
    ColourLines("1 6, 5 -4 7, -5 -8")
  }; 
  static const ColourLines diag47[1] = { 
    ColourLines("1 6, 5 7, -5 4 -8")
  }; 
  static const ColourLines diag48[1] = { 
    ColourLines("1 6, 5 -4 7, -5 -8")
  }; 

  Selector <const ColourLines *> sel;

  if( diag->id() == -1 )  {
   sel.insert( 1.0,  &(diag1[0]) );
  } 
  else if( diag->id() == -2 )  {
   sel.insert( 1.0,  &(diag2[0]) );
  } 
  else if( diag->id() == -3 )  {
   sel.insert( 1.0,  &(diag3[0]) );
  } 
  else if( diag->id() == -4 )  {
   sel.insert( 1.0,  &(diag4[0]) );
  } 
  else if( diag->id() == -5 )  {
   sel.insert( 1.0,  &(diag5[0]) );
  } 
  else if( diag->id() == -6 )  {
   sel.insert( 1.0,  &(diag6[0]) );
  } 
  else if( diag->id() == -7 )  {
   sel.insert( 1.0,  &(diag7[0]) );
  } 
  else if( diag->id() == -8 )  {
   sel.insert( 1.0,  &(diag8[0]) );
  } 
  else if( diag->id() == -9 )  {
   sel.insert( 1.0,  &(diag9[0]) );
  } 
  else if( diag->id() == -10 )  {
   sel.insert( 1.0,  &(diag10[0]) );
  } 
  else if( diag->id() == -11 )  {
   sel.insert( 1.0,  &(diag11[0]) );
  } 
  else if( diag->id() == -12 )  {
   sel.insert( 1.0,  &(diag12[0]) );
  } 
  else if( diag->id() == -13 )  {
   sel.insert( 1.0,  &(diag13[0]) );
  } 
  else if( diag->id() == -14 )  {
   sel.insert( 1.0,  &(diag14[0]) );
  } 
  else if( diag->id() == -15 )  {
   sel.insert( 1.0,  &(diag15[0]) );
  } 
  else if( diag->id() == -16 )  {
   sel.insert( 1.0,  &(diag16[0]) );
  } 
  else if( diag->id() == -17 )  {
   sel.insert( 1.0,  &(diag17[0]) );
  } 
  else if( diag->id() == -18 )  {
   sel.insert( 1.0,  &(diag18[0]) );
  } 
  else if( diag->id() == -19 )  {
   sel.insert( 1.0,  &(diag19[0]) );
  } 
  else if( diag->id() == -20 )  {
   sel.insert( 1.0,  &(diag20[0]) );
  } 
  else if( diag->id() == -21 )  {
   sel.insert( 1.0,  &(diag21[0]) );
  } 
  else if( diag->id() == -22 )  {
   sel.insert( 1.0,  &(diag22[0]) );
  } 
  else if( diag->id() == -23 )  {
   sel.insert( 1.0,  &(diag23[0]) );
  } 
  else if( diag->id() == -24 )  {
   sel.insert( 1.0,  &(diag24[0]) );
  } 
  else if( diag->id() == -25 )  {
   sel.insert( 1.0,  &(diag25[0]) );
  } 
  else if( diag->id() == -26 )  {
   sel.insert( 1.0,  &(diag26[0]) );
  } 
  else if( diag->id() == -27 )  {
   sel.insert( 1.0,  &(diag27[0]) );
  } 
  else if( diag->id() == -28 )  {
   sel.insert( 1.0,  &(diag28[0]) );
  } 
  else if( diag->id() == -29 )  {
   sel.insert( 1.0,  &(diag29[0]) );
  } 
  else if( diag->id() == -30 )  {
   sel.insert( 1.0,  &(diag30[0]) );
  } 
  else if( diag->id() == -31 )  {
   sel.insert( 1.0,  &(diag31[0]) );
  } 
  else if( diag->id() == -32 )  {
   sel.insert( 1.0,  &(diag32[0]) );
  } 
  else if( diag->id() == -33 )  {
   sel.insert( 1.0,  &(diag33[0]) );
  } 
  else if( diag->id() == -34 )  {
   sel.insert( 1.0,  &(diag34[0]) );
  } 
  else if( diag->id() == -35 )  {
   sel.insert( 1.0,  &(diag35[0]) );
  } 
  else if( diag->id() == -36 )  {
   sel.insert( 1.0,  &(diag36[0]) );
  } 
  else if( diag->id() == -37 )  {
   sel.insert( 1.0,  &(diag37[0]) );
  } 
  else if( diag->id() == -38 )  {
   sel.insert( 1.0,  &(diag38[0]) );
  } 
  else if( diag->id() == -39 )  {
   sel.insert( 1.0,  &(diag39[0]) );
  } 
  else if( diag->id() == -40 )  {
   sel.insert( 1.0,  &(diag40[0]) );
  } 
  else if( diag->id() == -41 )  {
   sel.insert( 1.0,  &(diag41[0]) );
  } 
  else if( diag->id() == -42 )  {
   sel.insert( 1.0,  &(diag42[0]) );
  } 
  else if( diag->id() == -43 )  {
   sel.insert( 1.0,  &(diag43[0]) );
  } 
  else if( diag->id() == -44 )  {
   sel.insert( 1.0,  &(diag44[0]) );
  } 
  else if( diag->id() == -45 )  {
   sel.insert( 1.0,  &(diag45[0]) );
  } 
  else if( diag->id() == -46 )  {
   sel.insert( 1.0,  &(diag46[0]) );
  } 
  else if( diag->id() == -47 )  {
   sel.insert( 1.0,  &(diag47[0]) );
  } 
  else if( diag->id() == -48 )  {
   sel.insert( 1.0,  &(diag48[0]) );
  } 
  return sel;
  }
}


ClassDescription<VBFNLOMEqg2hqqq> VBFNLOMEqg2hqqq::initVBFNLOMEqg2hqqq;
// Definition of the static class description member.



void VBFNLOMEqg2hqqq::persistentOutput(PersistentOStream & os) const {
  os << theWhichGluon;
}

void VBFNLOMEqg2hqqq::persistentInput(PersistentIStream & is, int) {
  is >> theWhichGluon;
}


void VBFNLOMEqg2hqqq::Init() {

  static ClassDocumentation<VBFNLOMEqg2hqqq> documentation
    ("VBFNLOMEqg2hqqq");

  static Switch<VBFNLOMEqg2hqqq,int> interfaceWhichGluon
    ("WhichGluon",
     "Set the position of the incoming gluon.",
     &VBFNLOMEqg2hqqq::theWhichGluon, 0, false, false);  
  static SwitchOption interfaceWhichGluonFirst
    (interfaceWhichGluon,
     "First",
     "From first incoming hadron.",
     0);
  static SwitchOption interfaceWhichGluonSecond
    (interfaceWhichGluon,
     "Second",
     "From second incoming hadron.",
     1);

}
