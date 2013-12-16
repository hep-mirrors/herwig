
#include "VBFNLOCommonBlocks.h"
#include "VBFNLOMEqq2Zqqg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/Throw.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/DiagramDrawer.h"

#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"

using namespace Herwig;


VBFNLOMEqq2Zqqg::VBFNLOMEqq2Zqqg():VBFNLOMEPP2ZJetJetJet(){
}

VBFNLOMEqq2Zqqg::~VBFNLOMEqq2Zqqg(){}

IBPtr VBFNLOMEqq2Zqqg::clone() const {
  return new_ptr(*this);
}

IBPtr VBFNLOMEqq2Zqqg::fullclone() const {
  return new_ptr(*this);
}

void VBFNLOMEqq2Zqqg::prepareMomenta(double (& pbar)[14][4],
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

  L5MomToDouble( (lastMEMomenta()[5]), &pbar[4][0]);

  L5MomToDouble( (lastMEMomenta()[6]), &pbar[5][0]);

  return;
} 



void VBFNLOMEqq2Zqqg::getDiagrams() const {
  tcPDPtr Wplus = getParticleData(ParticleID::Wplus);
  tcPDPtr Wminus = getParticleData(ParticleID::Wminus);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  tcPDPtr g = getParticleData(ParticleID::g);

  PDVector QuarksAndAntiQuarks1,QuarksAndAntiQuarks2;
  for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	q != theQuarkFlavours.end(); ++q ){
    if (requestedAsIncoming1(*q)) QuarksAndAntiQuarks1.push_back(*q);
    if (requestedAsIncoming1((**q).CC())) QuarksAndAntiQuarks1.push_back((**q).CC());
    if (requestedAsIncoming2(*q)) QuarksAndAntiQuarks2.push_back(*q);
    if (requestedAsIncoming2((**q).CC())) QuarksAndAntiQuarks2.push_back((**q).CC());
  }
  // get the decay products
  tcPDPtr dec1, dec2;
  if (theDecayLeptons->id() > 0){
    dec2 = theDecayLeptons->CC();
    dec1 = theDecayLeptons;
  }
  else {
    dec2 = theDecayLeptons;
    dec1 = theDecayLeptons->CC();
  }
  if ( theCurrent == neutral ) {
    for ( PDVector::const_iterator q1 = QuarksAndAntiQuarks1.begin();
	  q1 != QuarksAndAntiQuarks1.end(); ++q1 ){
      for ( PDVector::const_iterator q2 = QuarksAndAntiQuarks2.begin();
	    q2 != QuarksAndAntiQuarks2.end(); ++q2 ) {
	add(new_ptr((Tree2toNDiagram(5), *q1, *q1, Z0, Z0, *q2, 2, *q1,  4, *q2,  1, g,  3, Z0,  9,dec1,  9,dec2,    -1)));
	add(new_ptr((Tree2toNDiagram(5), *q1, Z0, Z0, *q2, *q2, 1, *q1,  3, *q2,  4, g,  2, Z0,  9,dec1,  9,dec2,    -2)));
	add(new_ptr((Tree2toNDiagram(4), *q1, Z0, Z0, *q2,  1, *q1,  5, *q1,  3, *q2,  5, g,  2, Z0,  9,dec1,  9,dec2,  -3)));
	add(new_ptr((Tree2toNDiagram(4), *q1, Z0, Z0, *q2,  3, *q2,  1, *q1,  5, *q2,  5, g,  2, Z0,  9,dec1,  9,dec2,  -4)));
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
	      add(new_ptr((Tree2toNDiagram(5), *q1, *q1, Wplus, Wplus, *q2, 2, q1prime,  4, q2prime,  1, g,  3, Z0,  9,dec1,  9,dec2,    -1)));
	      add(new_ptr((Tree2toNDiagram(5), *q1, Wplus, Wplus, *q2, *q2, 1, q1prime,  3, q2prime,  4, g,  2, Z0,  9,dec1,  9,dec2,    -2)));
	      add(new_ptr((Tree2toNDiagram(4), *q1, Wplus, Wplus, *q2,  1, q1prime,  5, q1prime,  3, q2prime,  5, g,  2, Z0,  9,dec1,  9,dec2,  -3)));
	      add(new_ptr((Tree2toNDiagram(4), *q1, Wplus, Wplus, *q2,  3, q2prime,  1, q1prime,  5, q2prime,  5, g,  2, Z0,  9,dec1,  9,dec2,  -4)));
	    }
	    else if ( chargediff < ZERO ) {
	      add(new_ptr((Tree2toNDiagram(5), *q1, *q1, Wminus, Wminus, *q2, 2, q1prime,  4, q2prime,  1, g,  3, Z0,  9,dec1,  9,dec2,    -1)));
	      add(new_ptr((Tree2toNDiagram(5), *q1, Wminus, Wminus, *q2, *q2, 1, q1prime,  3, q2prime,  4, g,  2, Z0,  9,dec1,  9,dec2,    -2)));
	      add(new_ptr((Tree2toNDiagram(4), *q1, Wminus, Wminus, *q2,  1, q1prime,  5, q1prime,  3, q2prime,  5, g,  2, Z0,  9,dec1,  9,dec2,  -3)));
	      add(new_ptr((Tree2toNDiagram(4), *q1, Wminus, Wminus, *q2,  3, q2prime,  1, q1prime,  5, q2prime,  5, g,  2, Z0,  9,dec1,  9,dec2,  -4)));
	    }
	  }
	}
      }
    }
  }
  else {
    cout << "THIS SHOULD NEVER BE EXECUTED" << endl<<endl;
    throw ThePEG::Exception() << "Please insert a pointer to W or Z "
			      << "boson as ExchangedBoson "
			      << "in the ME interface of "
			      << name() << "."
			      << ThePEG::Exception::abortnow;
  }
}

void VBFNLOMEqq2Zqqg::doinit(){ 
  VBFNLOMEBase::doinit();
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


Selector<MEBase::DiagramIndex> 
VBFNLOMEqq2Zqqg::diagrams(const DiagramVector & ) const {
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
VBFNLOMEqq2Zqqg::colourGeometries(tcDiagPtr diag) const {
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


ClassDescription<VBFNLOMEqq2Zqqg> VBFNLOMEqq2Zqqg::initVBFNLOMEqq2Zqqg;
// Definition of the static class description member.

void VBFNLOMEqq2Zqqg::persistentOutput(PersistentOStream & os) const {
}

void VBFNLOMEqq2Zqqg::persistentInput(PersistentIStream & is, int) {
}

void VBFNLOMEqq2Zqqg::Init() {

  static ClassDocumentation<VBFNLOMEqq2Zqqg> documentation
    ("VBFNLOMEqq2Zqqg");

}
