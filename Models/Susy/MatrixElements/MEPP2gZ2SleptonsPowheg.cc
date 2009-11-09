// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2gZ2SleptonsPowheg class.
//

#include "MEPP2gZ2SleptonsPowheg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig++/Models/Susy/SusyBase.h"

using namespace Herwig;

MEPP2gZ2SleptonsPowheg::MEPP2gZ2SleptonsPowheg() {
  massOption(true ,0);
  massOption(false,0);
}

void MEPP2gZ2SleptonsPowheg::getDiagrams() const {
}

unsigned int MEPP2gZ2SleptonsPowheg::orderInAlphaS() const {
  return 0;
}

unsigned int MEPP2gZ2SleptonsPowheg::orderInAlphaEW() const {
  return 2;
}

Selector<MEBase::DiagramIndex>
MEPP2gZ2SleptonsPowheg::diagrams(const DiagramVector & diags) const {
  // This example corresponds to the diagrams specified in the example
  // in the getDiagrams() function.

  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 ) sel.insert(1.0, i);
    else if ( diags[i]->id() == -2 )  sel.insert(1.0, i);
    else if ( diags[i]->id() == -3 )  sel.insert(1.0, i);
  // You probably do not want equal weights here...
  return sel;

  // If there is only one possible diagram you can override the
  // MEBase::diagram function instead.

}

Selector<const ColourLines *>
MEPP2gZ2SleptonsPowheg::colourGeometries(tcDiagPtr) const {
  static const ColourLines c1("1 -2");
  Selector<const ColourLines *> sel;
  sel.insert(1.0, &c1);
  return sel;
}


IBPtr MEPP2gZ2SleptonsPowheg::clone() const {
  return new_ptr(*this);
}

IBPtr MEPP2gZ2SleptonsPowheg::fullclone() const {
  return new_ptr(*this);
}

void MEPP2gZ2SleptonsPowheg::doinit() {
  NLODrellYanBase::doinit();
  Z0_    = getParticleData(ThePEG::ParticleID::Z0);
  gamma_ = getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  tcSusyBasePtr hwsm=ThePEG::dynamic_ptr_cast<tcSusyBasePtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    FFZVertex_ = hwsm->vertexFFZ();
    FFPVertex_ = hwsm->vertexFFP();
    WSSVertex_ = hwsm->vertexWSFSF();
  }
  else
    throw InitException() << "Must be the Herwig++ SusyBase class in "
			  << "MEPP2gZ2SleptonsPowheg::doinit" << Exception::abortnow;
}


void MEPP2gZ2SleptonsPowheg::persistentOutput(PersistentOStream & os) const {
  os << FFZVertex_ << FFPVertex_ << WSSVertex_ << Z0_ << gamma_;
}

void MEPP2gZ2SleptonsPowheg::persistentInput(PersistentIStream & is, int) {
  is >> FFZVertex_ >> FFPVertex_ >> WSSVertex_ >> Z0_ >> gamma_;
}

ClassDescription<MEPP2gZ2SleptonsPowheg> MEPP2gZ2SleptonsPowheg::initMEPP2gZ2SleptonsPowheg;
// Definition of the static class description member.

void MEPP2gZ2SleptonsPowheg::Init() {

  static ClassDocumentation<MEPP2gZ2SleptonsPowheg> documentation
    ("There is no documentation for the MEPP2gZ2SleptonsPowheg class");

}

NLODrellYanBase::Singular MEPP2gZ2SleptonsPowheg::virtualME() const {
  return NLODrellYanBase::Singular();
}

double MEPP2gZ2SleptonsPowheg::loME(const cPDVector & particles,
				    const vector<Lorentz5Momentum> & momenta) const {
  return 0.;
}

double MEPP2gZ2SleptonsPowheg::realME(const cPDVector & particles,
				      const vector<Lorentz5Momentum> & momenta) const {
  return 0.;
}
