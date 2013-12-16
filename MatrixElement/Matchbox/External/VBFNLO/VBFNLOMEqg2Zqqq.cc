
#include "VBFNLOCommonBlocks.h"
#include "VBFNLOMEqg2Zqqq.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Throw.h"

#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"

using namespace Herwig;


VBFNLOMEqg2Zqqq::VBFNLOMEqg2Zqqq()
  : VBFNLOMEPP2ZJetJetJet(), theWhichGluon(0){}

VBFNLOMEqg2Zqqq::~VBFNLOMEqg2Zqqq(){}

IBPtr VBFNLOMEqg2Zqqq::clone() const {
  return new_ptr(*this);
}

IBPtr VBFNLOMEqg2Zqqq::fullclone() const {
  return new_ptr(*this);
}

void VBFNLOMEqg2Zqqq::prepareMomenta(double (& pbar)[14][4],
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

  L5MomToDouble( (lastMEMomenta()[5]), &pbar[4][0]);

  L5MomToDouble( (lastMEMomenta()[6]), &pbar[5][0]);

  return;
} 



void VBFNLOMEqg2Zqqq::getDiagrams() const {

  tcPDPtr Wplus = getParticleData(ParticleID::Wplus);
  tcPDPtr Wminus = getParticleData(ParticleID::Wminus);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  tcPDPtr g = getParticleData(ParticleID::g);

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
    for ( PDVector::const_iterator q1 = theQuarkFlavours.begin();
	  q1 != theQuarkFlavours.end(); ++q1 ){
      for ( PDVector::const_iterator q2 = theQuarkFlavours.begin();
	    q2 != theQuarkFlavours.end(); ++q2 ) {
	//no-emission line is fermionic
	if (theWhichGluon == 0) {
	  if (requestedAsIncoming2(*q2)){
	    add(new_ptr((Tree2toNDiagram(5), g, *q1, Z0, Z0, *q2,  2, *q1,  4, *q2,  1, (**q1).CC(),  3, Z0,  9,dec1,  9,dec2, -1)));
	    add(new_ptr((Tree2toNDiagram(5), g, (**q1).CC(), Z0, Z0, *q2,  1, *q1,  4, *q2,  2, (**q1).CC(),  3, Z0,  9,dec1,  9,dec2, -2)));
	  }
	}
	else {
	  if (requestedAsIncoming1(*q1)){
	    add(new_ptr((Tree2toNDiagram(5), *q1, Z0, Z0, *q2, g,  1, *q1,  3, *q2,  4, (**q2).CC(),  2, Z0,  9,dec1,  9,dec2, -3)));
	    add(new_ptr((Tree2toNDiagram(5), *q1, Z0, Z0, (**q2).CC(), g,  1, *q1,  4, *q2,  3, (**q2).CC(),  2, Z0,  9,dec1,  9,dec2, -4)));
	  }
	}
	//no-emission line is antifermionic
	if (theWhichGluon == 0) {
	  if (requestedAsIncoming2((**q2).CC())){
	    add(new_ptr((Tree2toNDiagram(5), g, *q1, Z0, Z0, (**q2).CC(),  2, *q1,  4, (**q2).CC(),  1, (**q1).CC(),  3, Z0,  9,dec1,  9,dec2, -1)));
	    add(new_ptr((Tree2toNDiagram(5), g, (**q1).CC(), Z0, Z0, (**q2).CC(),  1, *q1,  4, (**q2).CC(),  2, (**q1).CC(),  3, Z0,  9,dec1,  9,dec2, -2)));
	  }
	}
	else {
	  if (requestedAsIncoming1((**q1).CC())){
	    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Z0, Z0, *q2, g,  1, (**q1).CC(),  3, *q2,  4, (**q2).CC(),  2, Z0,  9,dec1,  9,dec2, -3)));
	    add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Z0, Z0, (**q2).CC(), g,  1, (**q1).CC(),  4, *q2,  3, (**q2).CC(),  2, Z0,  9,dec1,  9,dec2, -4)));
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
		add(new_ptr((Tree2toNDiagram(5), g, *q1, Wplus, Wplus, *q2,  2, q1prime,  4, q2prime,  1, (**q1).CC(),  3, Z0,  9,dec1,  9,dec2, -1)));
		add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wplus, Wplus, *q2,  1, q1prime,  4, q2prime,  2, (**q1).CC(),  3, Z0,  9,dec1,  9,dec2, -2)));
	      }
	    }
	    else if (transfer > ZERO) {
	      if (requestedAsIncoming2(*q2)){
		add(new_ptr((Tree2toNDiagram(5), g, *q1, Wminus, Wminus, *q2,  2, q1prime,  4, q2prime,  1, (**q1).CC(),  3, Z0,  9,dec1,  9,dec2, -1)));
		add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wminus, Wminus, *q2,  1, q1prime,  4, q2prime,  2, (**q1).CC(),  3, Z0,  9,dec1,  9,dec2, -2)));
	      }
	    }
	  }
	  else {
	    //or no-emission line is upper line and fermionic
	    transfer = -(**q2).charge() + (*q2prime).charge();
	    if (transfer > ZERO) {
	      if (requestedAsIncoming1(*q1)){
		add(new_ptr((Tree2toNDiagram(5), *q1, Wplus, Wplus, *q2, g,  1, q1prime,  3, q2prime,  4, (**q2).CC(),  2, Z0,  9,dec1,  9,dec2, -3)));
		add(new_ptr((Tree2toNDiagram(5), *q1, Wplus, Wplus, qb2prime, g,  1, q1prime,  4, q2prime,  3, (**q2).CC(),  2, Z0,  9,dec1,  9,dec2, -4)));
	      }
	    }
	    else if (transfer < ZERO) {
	      if (requestedAsIncoming1(*q1)){
		add(new_ptr((Tree2toNDiagram(5), *q1, Wminus, Wminus, *q2, g,  1, q1prime,  3, q2prime,  4, (**q2).CC(),  2, Z0,  9,dec1,  9,dec2, -3)));
		add(new_ptr((Tree2toNDiagram(5), *q1, Wminus, Wminus, qb2prime, g,  1, q1prime,  4, q2prime,  3, (**q2).CC(),  2, Z0,  9,dec1,  9,dec2, -4)));
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
		add(new_ptr((Tree2toNDiagram(5), g, *q1, Wplus, Wplus, (**q2).CC(),  2, q1prime,  4, qb2prime,  1, (**q1).CC(),  3, Z0,  9,dec1,  9,dec2, -1)));
		add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wplus, Wplus, (**q2).CC(),  1, q1prime,  4, qb2prime,  2, (**q1).CC(),  3, Z0,  9,dec1,  9,dec2, -2)));
	      }
	    }
	    else if (transfer > ZERO) {
	      if (requestedAsIncoming2((**q2).CC())){
		add(new_ptr((Tree2toNDiagram(5), g, *q1, Wminus, Wminus, (**q2).CC(),  2, q1prime,  4, qb2prime,  1, (**q1).CC(),  3, Z0,  9,dec1,  9,dec2, -1)));
		add(new_ptr((Tree2toNDiagram(5), g, qb1prime, Wminus, Wminus, (**q2).CC(),  1, q1prime,  4, qb2prime,  2, (**q1).CC(),  3, Z0,  9,dec1,  9,dec2, -2)));
	      }
	    }
	  }
	  else {
	    //or no-emission line is upper line and antifermionic
	    transfer = -(**q2).charge() + (*q2prime).charge();
	    if (transfer > ZERO) {
	      if (requestedAsIncoming1((**q1).CC())){
		add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wplus, Wplus, *q2, g,  1, qb1prime,  3, q2prime,  4, (**q2).CC(),  2, Z0,  9,dec1,  9,dec2, -3)));
		add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wplus, Wplus, qb2prime, g,  1, qb1prime,  4, q2prime,  3, (**q2).CC(),  2, Z0,  9,dec1,  9,dec2, -4)));
	      }
	    }
	    else if (transfer < ZERO) {
	      if (requestedAsIncoming1((**q1).CC())){
		add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wminus, Wminus, *q2, g,  1, qb1prime,  3, q2prime,  4, (**q2).CC(),  2, Z0,  9,dec1,  9,dec2, -3)));
		add(new_ptr((Tree2toNDiagram(5), (**q1).CC(), Wminus, Wminus, qb2prime, g,  1, qb1prime,  4, q2prime,  3, (**q2).CC(),  2, Z0,  9,dec1,  9,dec2, -4)));
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


void VBFNLOMEqg2Zqqq::doinit(){
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
VBFNLOMEqg2Zqqq::diagrams(const DiagramVector & diags) const {
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
VBFNLOMEqg2Zqqq::colourGeometries(tcDiagPtr diag) const {


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
  return sel;
}


ClassDescription<VBFNLOMEqg2Zqqq> VBFNLOMEqg2Zqqq::initVBFNLOMEqg2Zqqq;
// Definition of the static class description member.

void VBFNLOMEqg2Zqqq::persistentOutput(PersistentOStream & os) const {
  os << theWhichGluon;
}

void VBFNLOMEqg2Zqqq::persistentInput(PersistentIStream & is, int) {
  is >> theWhichGluon;
}


void VBFNLOMEqg2Zqqq::Init() {

  static ClassDocumentation<VBFNLOMEqg2Zqqq> documentation
    ("VBFNLOMEqg2Zqqq");

  static Switch<VBFNLOMEqg2Zqqq,int> interfaceWhichGluon
    ("WhichGluon",
     "Set the position of the incoming gluon.",
     &VBFNLOMEqg2Zqqq::theWhichGluon, 0, false, false);  
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
