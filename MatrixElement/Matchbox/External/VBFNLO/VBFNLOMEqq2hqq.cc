
#include "VBFNLOCommonBlocks.h"
#include "VBFNLOMEqq2hqq.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"

using namespace Herwig;


VBFNLOMEqq2hqq::VBFNLOMEqq2hqq()
  : VBFNLOMEPP2hJetJet(), theDecayChannel(0), theNarrowWidth(true) {
}

VBFNLOMEqq2hqq::~VBFNLOMEqq2hqq(){}

IBPtr VBFNLOMEqq2hqq::clone() const {
  return new_ptr(*this);
}

IBPtr VBFNLOMEqq2hqq::fullclone() const {
  return new_ptr(*this);
}




void VBFNLOMEqq2hqq::getDiagrams() const {
  bool useOldStyle = false;

  if (useOldStyle) {
    if(theQuarkFlavours.size() == 0) throw ThePEG::Exception() << "Please insert pointers to quark flavours "
							       << "in the ME interface."
							       << ThePEG::Exception::abortnow;

    tcPDPtr Wplus = getParticleData(ParticleID::Wplus);
    tcPDPtr Wminus = getParticleData(ParticleID::Wminus);
    tcPDPtr Z0 = getParticleData(ParticleID::Z0);
    tcPDPtr h0 = getParticleData(ParticleID::h0);


    PDVector QuarksAndAntiQuarks1,QuarksAndAntiQuarks2;
    for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	  q != theQuarkFlavours.end(); ++q ){
      if (requestedAsIncoming1(*q)) QuarksAndAntiQuarks1.push_back(*q);
      if (requestedAsIncoming1((**q).CC())) QuarksAndAntiQuarks1.push_back((**q).CC());
      if (requestedAsIncoming2(*q)) QuarksAndAntiQuarks2.push_back(*q);
      if (requestedAsIncoming2((**q).CC())) QuarksAndAntiQuarks2.push_back((**q).CC());
    }

    if (NDecayProducts() == 1) {
      if (theCurrent == neutral) {
	for ( PDVector::const_iterator q1 = QuarksAndAntiQuarks1.begin();
	      q1 != QuarksAndAntiQuarks1.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = QuarksAndAntiQuarks2.begin();
		q2 != QuarksAndAntiQuarks2.end(); ++q2 ) {
	    add(new_ptr((Tree2toNDiagram(4), *q1, Z0, Z0, *q2,   1, *q1,  3, *q2, 2, h0, -1)));
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
		  add(new_ptr((Tree2toNDiagram(4),*q1,Wplus,Wplus,*q2,   1,q1prime,  3,q2prime,  2,h0, -2)));
		}
		else if ( chargediff < ZERO ) {
		  add(new_ptr((Tree2toNDiagram(4),*q1,Wminus,Wminus,*q2,   1,q1prime,  3,q2prime,  2,h0, -2)));
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

      if (theDecayChannel == 1) {
	dec1 = getParticleData(ParticleID::gamma);
	dec2 = getParticleData(ParticleID::gamma);
      }
      else if (theDecayChannel == 2) {
	dec1 = getParticleData(ParticleID::muplus);
	dec2 = getParticleData(ParticleID::muminus);
      }
      else if (theDecayChannel == 3) {
	dec1 = getParticleData(ParticleID::tauplus);
	dec2 = getParticleData(ParticleID::tauminus);
      }
      else if (theDecayChannel == 4) {
	dec1 = getParticleData(ParticleID::b);
	dec2 = getParticleData(ParticleID::bbar);
      }

      if ( theCurrent == neutral ) {
	for ( PDVector::const_iterator q1 = QuarksAndAntiQuarks1.begin();
	      q1 != QuarksAndAntiQuarks1.end(); ++q1 ){
	  for ( PDVector::const_iterator q2 = QuarksAndAntiQuarks2.begin();
		q2 != QuarksAndAntiQuarks2.end(); ++q2 ) {
	    add(new_ptr((Tree2toNDiagram(4), *q1, Z0, Z0, *q2,   1, *q1,  3, *q2, 2, h0,  7,dec1,  7, dec2, -1)));
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
		  add(new_ptr((Tree2toNDiagram(4),*q1,Wplus,Wplus,*q2,   1,q1prime,  3,q2prime,  2,h0,  7,dec1,  7, dec2, -2)));
		}
		else if ( chargediff < ZERO ) {
		  add(new_ptr((Tree2toNDiagram(4),*q1,Wminus,Wminus,*q2,   1,q1prime,  3,q2prime,  2,h0,  7,dec1,  7, dec2, -2)));
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

      if (theDecayChannel == 5) {
	VecBos1 = getParticleData(ParticleID::Wplus);
	VecBos2 = getParticleData(ParticleID::Wminus);
	dec1 = getParticleData(ParticleID::nu_e);
	dec2 = getParticleData(ParticleID::eplus);
	dec3 = getParticleData(ParticleID::nu_ebar);
	dec4 = getParticleData(ParticleID::eminus);
      }
      else if (theDecayChannel == 6) {
	VecBos1 = getParticleData(ParticleID::Z0);
	VecBos2 = getParticleData(ParticleID::Z0);
	dec1 = getParticleData(ParticleID::eminus);
	dec2 = getParticleData(ParticleID::eplus);
	dec3 = getParticleData(ParticleID::eminus);
	dec4 = getParticleData(ParticleID::eplus);
      }
      else if (theDecayChannel == 7) {
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
	    add(new_ptr((Tree2toNDiagram(4), *q1, Z0, Z0, *q2,   1, *q1,  3, *q2, 2, h0,  7,VecBos1,  7,VecBos2,  8,dec1,  8,dec2,  9,dec3,  9,dec4, -1)));
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
		  add(new_ptr((Tree2toNDiagram(4),*q1,Wplus,Wplus,*q2,   1,q1prime,  3,q2prime,  2,h0,  7,VecBos1,  7,VecBos2,  8,dec1,  8,dec2,  9,dec3,  9,dec4, -2)));
		}
		else if ( chargediff < ZERO ) {
		  add(new_ptr((Tree2toNDiagram(4),*q1,Wminus,Wminus,*q2,   1,q1prime,  3,q2prime,  2,h0,  7,VecBos1,  7,VecBos2,  8,dec1,  8,dec2,  9,dec3,  9,dec4, -2)));
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
    else {
      throw ThePEG::Exception() << "The number of decay products isn't supported."
				<< "This may have something to do with an incorrect setting"
				<< "in the ME interface of the decay channel."
				<< ThePEG::Exception::abortnow;
    }
  }

  else { // ! useOldStyle

    PDPtr Wplus = getParticleData(ParticleID::Wplus); 
    PDPtr Wminus = getParticleData(ParticleID::Wminus); 
    PDPtr Z0 = getParticleData(ParticleID::Z0); 
    PDPtr h0 = getParticleData(ParticleID::h0); 
    
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

	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Wminus, Wminus, kd, 1, qdBC, 3, kdC, 2, h0, -1)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Wminus, Wminus, kuB, 1, qdBC, 3, kuBC, 2, h0, -2)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Wplus, Wplus, kdB, 1, qdC, 3, kdBC, 2, h0, -3)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Wplus, Wplus, ku, 1, qdC, 3, kuC, 2, h0, -4)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Wplus, Wplus, kdB, 1, quBC, 3, kdBC, 2, h0, -5)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Wplus, Wplus, ku, 1, quBC, 3, kuC, 2, h0, -6)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Wminus, Wminus, kd, 1, quC, 3, kdC, 2, h0, -7)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Wminus, Wminus, kuB, 1, quC, 3, kuBC, 2, h0, -8)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Z0, Z0, kd, 1, qdB, 3, kd, 2, h0, -9)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Z0, Z0, kdB, 1, qdB, 3, kdB, 2, h0, -10)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Z0, Z0, ku, 1, qdB, 3, ku, 2, h0, -11)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qdB, Z0, Z0, kuB, 1, qdB, 3, kuB, 2, h0, -12)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, kd, 1, qd, 3, kd, 2, h0, -13)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, kdB, 1, qd, 3, kdB, 2, h0, -14)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, ku, 1, qd, 3, ku, 2, h0, -15)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, kuB, 1, qd, 3, kuB, 2, h0, -16)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Z0, Z0, kd, 1, quB, 3, kd, 2, h0, -17)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Z0, Z0, kdB, 1, quB, 3, kdB, 2, h0, -18)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Z0, Z0, ku, 1, quB, 3, ku, 2, h0, -19)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), quB, Z0, Z0, kuB, 1, quB, 3, kuB, 2, h0, -20)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, kd, 1, qu, 3, kd, 2, h0, -21)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, kdB, 1, qu, 3, kdB, 2, h0, -22)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, ku, 1, qu, 3, ku, 2, h0, -23)));
	allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, kuB, 1, qu, 3, kuB, 2, h0, -24)));

      }
    }

    //now check which of those are really requested by the interfaced parameters
    int counter = 0;
    for (vector<DiagPtr>::const_iterator di = allPossibleDiagrams.begin();
	 di != allPossibleDiagrams.end() ; di++) {
      if ( allowedDiagram(*di) ) {
	add(*di);
	counter ++;
      }
    }
  }
}

Selector<MEBase::DiagramIndex> 
VBFNLOMEqq2hqq::diagrams(const DiagramVector & diags) const {
  Selector<MEBase::DiagramIndex> sel;
  sel.insert(1,0);
  return sel;
}

Selector<const ColourLines *>
VBFNLOMEqq2hqq::colourGeometries(tcDiagPtr diag) const {
  bool useOldStyle = false;
  if (useOldStyle) {
    if (theDecayChannel != 4) {
      static const ColourLines cUL      (" 1  5,  4  6");
      static const ColourLines cULbar   (" 1  5, -4 -6");
      static const ColourLines cUbarL   ("-1 -5,  4  6");
      static const ColourLines cUbarLbar("-1 -5, -4 -6");
      Selector<const ColourLines *> sel;
      if ( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() > 0)
	sel.insert(1.0, &cUL);
      else if ( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() < 0)
	sel.insert(1.0, &cULbar);
      else if ( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() > 0)
	sel.insert(1.0, &cUbarL);
      else if ( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() < 0)
	sel.insert(1.0, &cUbarLbar);
      return sel;
    }
    else {
      static const ColourLines cUL      (" 1  5,  4  6,  8 -9");
      static const ColourLines cULbar   (" 1  5, -4 -6,  8 -9");
      static const ColourLines cUbarL   ("-1 -5,  4  6,  8 -9");
      static const ColourLines cUbarLbar("-1 -5, -4 -6,  8 -9");
      Selector<const ColourLines *> sel;
      if ( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() > 0)
	sel.insert(1.0, &cUL);
      else if ( mePartonData()[0]->id() > 0 && mePartonData()[1]->id() < 0)
	sel.insert(1.0, &cULbar);
      else if ( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() > 0)
	sel.insert(1.0, &cUbarL);
      else if ( mePartonData()[0]->id() < 0 && mePartonData()[1]->id() < 0)
	sel.insert(1.0, &cUbarLbar);
      return sel;
    }
  }
  else {
    static const ColourLines diag1[1] = { ColourLines("-1 -5, 4 6") }; 
    static const ColourLines diag2[1] = { ColourLines("-1 -5, -4 -6") }; 
    static const ColourLines diag3[1] = { ColourLines("1 5, -4 -6") }; 
    static const ColourLines diag4[1] = { ColourLines("1 5, 4 6") }; 
    static const ColourLines diag5[1] = { ColourLines("-1 -5, -4 -6") }; 
    static const ColourLines diag6[1] = { ColourLines("-1 -5, 4 6") }; 
    static const ColourLines diag7[1] = { ColourLines("1 5, 4 6") }; 
    static const ColourLines diag8[1] = { ColourLines("1 5, -4 -6") }; 
    static const ColourLines diag9[1] = { ColourLines("-1 -5, 4 6") }; 
    static const ColourLines diag10[1] = { ColourLines("-1 -5, -4 -6") }; 
    static const ColourLines diag11[1] = { ColourLines("-1 -5, 4 6") }; 
    static const ColourLines diag12[1] = { ColourLines("-1 -5, -4 -6") }; 
    static const ColourLines diag13[1] = { ColourLines("1 5, 4 6") }; 
    static const ColourLines diag14[1] = { ColourLines("1 5, -4 -6") }; 
    static const ColourLines diag15[1] = { ColourLines("1 5, 4 6") }; 
    static const ColourLines diag16[1] = { ColourLines("1 5, -4 -6") }; 
    static const ColourLines diag17[1] = { ColourLines("-1 -5, 4 6") }; 
    static const ColourLines diag18[1] = { ColourLines("-1 -5, -4 -6") }; 
    static const ColourLines diag19[1] = { ColourLines("-1 -5, 4 6") }; 
    static const ColourLines diag20[1] = { ColourLines("-1 -5, -4 -6") }; 
    static const ColourLines diag21[1] = { ColourLines("1 5, 4 6") }; 
    static const ColourLines diag22[1] = { ColourLines("1 5, -4 -6") }; 
    static const ColourLines diag23[1] = { ColourLines("1 5, 4 6") }; 
    static const ColourLines diag24[1] = { ColourLines("1 5, -4 -6") }; 

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
    return sel;
  }
}

void VBFNLOMEqq2hqq::doinit(){
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

void VBFNLOMEqq2hqq::initPSGen() {

  if (NDecayProducts() != 4) {
    int bos=6;

    BLIPSIVNJ.RM2=BKOPOU.XM2[bos-1];
    BLIPSIVNJ.RMG=BKOPOU.XMG[bos-1];

    BLIPSIVNJ.RM2MIN=sqrt(BLIPSIVNJ.RM2)*( 1 - 150*BLIPSIVNJ.RMG/BLIPSIVNJ.RM2 );
    if (BLIPSIVNJ.RM2MIN > 0) BLIPSIVNJ.RM2MIN = BLIPSIVNJ.RM2MIN*BLIPSIVNJ.RM2MIN;
    else BLIPSIVNJ.RM2MIN = 0;

    BLIPSIVNJ.RM2MAX=BLIPSIVNJ.RM2*pow( 1 + 150*BLIPSIVNJ.RMG/BLIPSIVNJ.RM2, 2);
    BLIPSIVNJ.RM2MAX = min(BLIPSIVNJ.RM2MAX,pow(sqrt(BLIPSIVNJ.RM2) + 200, 2));

    BLIPSIVNJ.S=lastS()/GeV2;
    BLIPSIVNJ.M2MIN=0.001;

    //these should probably be taken from some cut class eventually
    BLIPSIVNJ.YJMIN[0]=0;
    BLIPSIVNJ.YJMIN[1]=0;

    BLIPSIVNJ.YJMAX[0]=5;
    BLIPSIVNJ.YJMAX[1]=5;

    BLIPSIVNJ.PTJMIN[0]=20;
    BLIPSIVNJ.PTJMIN[1]=20;

    BLIPSIVNJ.EJMIN[0]=0;
    BLIPSIVNJ.EJMIN[1]=0;

    BLIPSIVNJ.INFOJ[0]=-1;
    BLIPSIVNJ.INFOJ[1]=-1;
  }
  else {
    int bos[3];
    bos[0] = 6;
    if (theDecayChannel == 5) {
      bos[1] = 3;
      bos[2] = 4;
    }
    else {
      bos[1] = 2;
      bos[2] = 2;
    }

    BLIPSIVVNJ.S = lastS()/GeV2;
    BLIPSIVVNJ.M2MIN = 0.01;

    BLIPSIVVNJ.IWIDTH[0] = 0;
    BLIPSIVVNJ.IWIDTH[1] = 1;
    BLIPSIVVNJ.IWIDTH[2] = 1;

    for (int i = 0; i < 3; i++) {
      BLIPSIVVNJ.RM2[i] = BKOPOU.XM2[bos[i]-1];
      BLIPSIVVNJ.RMG[i] = BKOPOU.XMG[bos[i]-1];
    }

    BLIPSIVVNJ.RM2MIN[0] = sqrt(BLIPSIVVNJ.RM2[0])*( 1 - 150*BLIPSIVVNJ.RMG[0]/BLIPSIVVNJ.RM2[0] );
    if (BLIPSIVVNJ.RM2MIN[0] > 0) BLIPSIVVNJ.RM2MIN[0] = BLIPSIVVNJ.RM2MIN[0]*BLIPSIVVNJ.RM2MIN[0];
    else BLIPSIVVNJ.RM2MIN[0] = 0;
    BLIPSIVVNJ.RM2MAX[0] = BLIPSIVVNJ.RM2[0]*pow( 1 + 150*BLIPSIVVNJ.RMG[0]/BLIPSIVVNJ.RM2[0], 2);
    BLIPSIVVNJ.RM2MAX[0] = min(BLIPSIVVNJ.RM2MAX[0],BLIPSIVVNJ.S/4);

    for (int i = 1; i < 3; i++) {
      if (bos[i] == 2) {
	BLIPSIVVNJ.RM2MIN[i] = 225;
	BLIPSIVVNJ.RM2MAX[i] = BLIPSIVVNJ.RM2MAX[0];
      }
      else {
	BLIPSIVVNJ.RM2MIN[i] = 0.000001;
	BLIPSIVVNJ.RM2MAX[i] = BLIPSIVVNJ.S/4;
      }
    }

    //these should probably be taken from some cut class eventually
    BLIPSIVVNJ.YJMIN[0] = 0;
    BLIPSIVVNJ.YJMIN[1] = 0;

    BLIPSIVVNJ.YJMAX[0] = 5;
    BLIPSIVVNJ.YJMAX[1] = 5;

    BLIPSIVVNJ.PTJMIN[0] = 20;
    BLIPSIVVNJ.PTJMIN[1] = 20;

    BLIPSIVVNJ.EJMIN[0] = 0;
    BLIPSIVVNJ.EJMIN[1] = 0;

    BLIPSIVVNJ.INFOJ[0] = -1;
    BLIPSIVVNJ.INFOJ[1] = -1;
  }
  return;
}

bool VBFNLOMEqq2hqq::generateKinematics(const double * r){

  if ( phasespace() ) {
    initPSGen();
    return MatchboxMEBase::generateKinematics(r);
  }

  int Fn = Njets();
  double Frd[100], Frn, Fk1[4],Fk2[4],Fq[5],Fd[10][4],Fd1[4],Fd2[4],Fp[2][4],Fx1,Fx2,Fw;
  int Fnw = 0;

  if (theNarrowWidth) Fnw =1;

  for (int i=0; i<nDim()-1; i++){
    //nDim()-1 because one random number is needed for Frn
    Frd[i]=r[i];
  }
  Frn = r[nDim()-1];
  
  initPSGen();

  if (NDecayProducts() == 1) {
    LIPSN0(Fn,Frd,Frn,Fk1,Fk2,Fq,Fp,Fx1,Fx2,Fw);
  }
  else if (NDecayProducts() == 2) {
    LIPSN(Fn,Frd,Frn,Fk1,Fk2,Fq,Fd1,Fd2,Fp,Fx1,Fx2,Fw,Fnw);
  }
  else if (NDecayProducts() == 4) {
    LIPSNVV(Fn,Frd,Frn,Fk1,Fk2,Fd,Fp,Fx1,Fx2,Fw);
  }
  if(Fw > 0){
    lastMEMomenta()[0]=(Lorentz5Momentum(Fk1[1]*GeV,Fk1[2]*GeV,Fk1[3]*GeV,Fk1[0]*GeV));
    lastMEMomenta()[1]=(Lorentz5Momentum(Fk2[1]*GeV,Fk2[2]*GeV,Fk2[3]*GeV,Fk2[0]*GeV));
    for (int i=0; i<Fn; i++){
      lastMEMomenta()[2+i]=(Lorentz5Momentum(Fp[i][1]*GeV,Fp[i][2]*GeV,Fp[i][3]*GeV,Fp[i][0]*GeV));
    }

    if (NDecayProducts() == 1)
      lastMEMomenta()[4]=(Lorentz5Momentum(Fq[1]*GeV,Fq[2]*GeV,Fq[3]*GeV,Fq[0]*GeV));    
    else if (NDecayProducts() == 2) {
      lastMEMomenta()[4]=(Lorentz5Momentum(Fd1[1]*GeV,Fd1[2]*GeV,Fd1[3]*GeV,Fd1[0]*GeV));  
      lastMEMomenta()[5]=(Lorentz5Momentum(Fd2[1]*GeV,Fd2[2]*GeV,Fd2[3]*GeV,Fd2[0]*GeV));  
    }
    else if (NDecayProducts() == 4) {   
      lastMEMomenta()[4]=(Lorentz5Momentum(Fd[0][1]*GeV,Fd[0][2]*GeV,Fd[0][3]*GeV,Fd[0][0]*GeV));  
      lastMEMomenta()[5]=(Lorentz5Momentum(Fd[1][1]*GeV,Fd[1][2]*GeV,Fd[1][3]*GeV,Fd[1][0]*GeV));  
      lastMEMomenta()[6]=(Lorentz5Momentum(Fd[2][1]*GeV,Fd[2][2]*GeV,Fd[2][3]*GeV,Fd[2][0]*GeV));  
      lastMEMomenta()[7]=(Lorentz5Momentum(Fd[3][1]*GeV,Fd[3][2]*GeV,Fd[3][3]*GeV,Fd[3][0]*GeV));
    }
    //division by sHat takes place in LIPSN, but this is ok for 
    //hjj (2 to 3 process)
    //for processes with h decay, we need to divide by some more powers of sHat
    if (NDecayProducts() != 1){
      Energy2 thisSHat=(meMomenta()[0] + meMomenta()[1]).m2();
      Fw/=pow(thisSHat/GeV2,NDecayProducts()-1);
    }
    //revert the multiplication with (hbar*c)^2/2   
    //and divide by Feynman x-es
    jacobian(Fw*2/3.8937E11/Fx1/Fx2);
    setScale();
    logGenerateKinematics(r);
    return true; 
  } else { 
    jacobian(0.0);
    return false;
  }
}

int VBFNLOMEqq2hqq::NDecayProducts() const {
  if (theDecayChannel == 0) return 1;
  if (theDecayChannel < 5) return 2;
  if (theDecayChannel < 8) return 4;
  return 0;
}

double VBFNLOMEqq2hqq::BranchingRatio() const {
  if (theDecayChannel == 1) {
    return BRANCH.BHGAM;
  }
  if (theDecayChannel == 2) {
    return BRANCH.BHMU;
  }
  if (theDecayChannel == 3) {
    return BRANCH.BHTAU;
  }
  if (theDecayChannel == 4) {
    return BRANCH.BHBB;
  }
  return 1;
}

int VBFNLOMEqq2hqq::nDim() const {
  if (phasespace()) return phasespace()->nDim(3+NDecayProducts()-1);
  if (NDecayProducts() == 1) return 7;
  if (NDecayProducts() == 2) return 10;
  if (NDecayProducts() == 4) return 16;
  return 0;
}

ClassDescription<VBFNLOMEqq2hqq> VBFNLOMEqq2hqq::initVBFNLOMEqq2hqq;
// Definition of the static class description member.

void VBFNLOMEqq2hqq::persistentOutput(PersistentOStream & os) const {
  os  << theDecayChannel << theNarrowWidth;
}

void VBFNLOMEqq2hqq::persistentInput(PersistentIStream & is, int) {
  is  >> theDecayChannel >> theNarrowWidth;
}


void VBFNLOMEqq2hqq::Init() {

  static ClassDocumentation<VBFNLOMEqq2hqq> documentation
    ("VBFNLOMEqq2hqq");

  static Switch<VBFNLOMEqq2hqq,int> interfaceDecayChannel
    ("DecayChannel",
     "Choose the decay channel that is simulated by VBFNLO.",
     &VBFNLOMEqq2hqq::theDecayChannel, 0, true, false);
  static SwitchOption interfaceDecayChannelStable
    (interfaceDecayChannel,
     "Stable",
     "No Higgs decay is calculated within VBFNLO.",
     0);
  static SwitchOption interfaceDecayChannelAA
    (interfaceDecayChannel,
     "H -> A A",
     "Higgs decay into photons",
     1);
  static SwitchOption interfaceDecayChannelMu
    (interfaceDecayChannel,
     "H -> mu mu",
     "Higgs decay into muons",
     2);
  static SwitchOption interfaceDecayChannelTau
    (interfaceDecayChannel,
     "H -> tau tau",
     "Higgs decay into taus",
     3);
  static SwitchOption interfaceDecayChannelBBar
    (interfaceDecayChannel,
     "H -> b bbar",
     "Higgs decay into b anti-b",
     4);
  static SwitchOption interfaceDecayChannelWW
    (interfaceDecayChannel,
     "H -> W+ W-",
     "Higgs decay into W bosons",
     5);
  static SwitchOption interfaceDecayChannelZZ_ll
    (interfaceDecayChannel,
     "H -> Z Z -> l lbar",
     "Higgs decay into Z bosons into lepton antilepton",
     6);
  static SwitchOption interfaceDecayChannelZZ_lnu
    (interfaceDecayChannel,
     "H -> Z Z -> l nu",
     "Higgs decay into Z bosons into lepton neutrino",
     7);

  static Switch<VBFNLOMEqq2hqq,bool> interfaceNarrowWidth
    ("NarrowWidth",
     "Choose if VBFNLO simulates the Higgs decay in narrow width approximation",
     &VBFNLOMEqq2hqq::theNarrowWidth, true, true, false);
  static SwitchOption interfaceNarrowWidthTrue
    (interfaceNarrowWidth,
     "True",
     "Calculate with narrow width approximation",
     true);
  static SwitchOption interfaceNarrowWidthFalse
    (interfaceNarrowWidth,
     "False",
     "No narrow width approximation",
     false);

}
