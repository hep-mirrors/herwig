// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPWWHVertex class.
//

#include "LHTPWWHVertex.h"
#include "LHTPModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr LHTPWWHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPWWHVertex::fullclone() const {
  return new_ptr(*this);
}

void LHTPWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(coup_,GeV);
}

void LHTPWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coup_,GeV);
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPWWHVertex,VVSVertex>
describeHerwigLHTPWWHVertex("Herwig::LHTPWWHVertex", "HwLHTPModel.so");

void LHTPWWHVertex::Init() {

  static ClassDocumentation<LHTPWWHVertex> documentation
    ("The LHTPWWHVertex class implements the coupling of two electroweak"
     " gauge bosons to a Higgs boson in the Little Higgs Model with T-Parity"
     "including the additional heavy photon, Z and W bosons and the "
     "triplet Higgs bosons.");

}

LHTPWWHVertex::LHTPWWHVertex() : coupLast_(0.), q2Last_(0.*GeV2) {
  // order in the couplings
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void LHTPWWHVertex::doinit() {
  // W_L W_L H
  addToList(  24,  -24,    25);
  // Z_L Z_L H
  addToList(  23,   23,    25);
  // W_H W_H H
  addToList(  34,  -34,    25);
  // Z_H Z_H H
  addToList(  33,   33,    25);
  // A_H A_H H
  addToList(  32,   32,    25);
  // Z_H A_H H
  addToList(  33,   32,    25);
  // Z_L Z_H Phi0
  addToList(  23,   33,    35);
  // A_H Z_L Phi0
  addToList(  32,   23,    35);
  // W_L W_H PhiP
  addToList(  24,  -34,    36);
  addToList(  34,  -24,    36);
  // W_H Z_L Phi+/- 
  addToList(  34,   23,   -37);
  addToList( -34,   23,    37);
  // W_L A_H Phi+/-
  addToList(  24,   32,   -37);
  addToList( -24,   32,    37);
  // W_L Z_H Phi+/- 
  addToList(  24,   33,   -37);
  addToList( -24,   33,    37);
  // W_H A_L Phi+/- 
  addToList(  34,   22,   -37);
  addToList( -34,   22,    37);
  // W_L W_H Phi --/++
  addToList(  24,   34,   -38);
  addToList( -24,  -34,    38);
  // model
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHTPModel "
			  << " in LHTPWWHVertex::doinit()"
			  << Exception::runerror;
  // base class
  VVSVertex::doinit();
  // calculate the couplings for the different combinations of particles
  Energy fact = 0.5*model->vev()/model->sin2ThetaW();
  double sw(sqrt(model->sin2ThetaW())),cw(sqrt(1.-model->sin2ThetaW()));
  double vf(model->vev()/model->f());
  double r2(sqrt(2.));
  coup_.resize(14);
  // H 
  coup_[ 0] = fact        *(1.-sqr(vf)/3.);
  coup_[ 1] = fact/sqr(cw)*(1.-sqr(vf)/3.);
  coup_[ 2] =-fact;
  coup_[ 3] =-fact;
  coup_[ 4] =-fact*sqr(sw/cw);
  coup_[ 5] =-fact/cw*sw;
  // PhiP 
  coup_[ 6] = r2*fact*vf/3.;
  // Phi0
  coup_[ 7] =-fact*vf/r2/cw;
  coup_[ 8] = fact*vf/r2*sw/sqr(cw);
  // Phi+
  coup_[ 9] = fact*vf/6./cw*(1.+2.*sqr(sw));
  coup_[10] = fact*vf*sw/cw*0.5;
  coup_[11] = fact*vf*5./6.;
  coup_[12] =-fact*vf*sw/3.;
  // Phi++
  coup_[13] =-fact*vf;
}

void LHTPWWHVertex::setCoupling(Energy2 q2,tcPDPtr a,
					      tcPDPtr b, tcPDPtr c) {
  // first the overall normalisation
  if(q2!=q2Last_) {
    coupLast_ = sqr(electroMagneticCoupling(q2));
    q2Last_=q2;
  }
  long ih = abs(c->id());
  long ibos[2]={abs(a->id()),abs(b->id())};
  
  if(ih == 25) {
    if( ibos[0] == 24 && ibos[1] == 24) 
      norm(UnitRemoval::InvE *coupLast_*coup_[0]);
    else if( ibos[0] == 23 && ibos[1] == 23     ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[1]);
    else if( ibos[0] == 34 && ibos[1] == 34     ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[2]);
    else if( ibos[0] == 33 && ibos[1] == 33     ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[3]);
    else if( ibos[0] == 32 && ibos[1] == 32     ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[4]);
    else if((ibos[0] == 33 && ibos[1] == 32) ||
	    (ibos[0] == 32 && ibos[1] == 33)    ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[5]);
    else 
      assert(false);
  }
  else if(ih == 36) {
    if( a->id() == 34 || b->id() == 34) 
      norm( Complex(0.,1.)*UnitRemoval::InvE *coupLast_*coup_[6]);
    else
      norm(-Complex(0.,1.)*UnitRemoval::InvE *coupLast_*coup_[6]);
  }
  else if(ih == 35) {
    if((ibos[0] == 23 && ibos[1] == 33) ||
       (ibos[0] == 33 && ibos[1] == 23)    ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[7]);
    else if((ibos[0] == 23 && ibos[1] == 32) ||
	    (ibos[0] == 32 && ibos[1] == 23)    ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[8]);
    else 
      assert(false);
  }
  else if(ih == 37) {
    if((ibos[0] == 34 && ibos[1] == 23) ||
       (ibos[0] == 23 && ibos[1] == 34)    ) {
      norm(UnitRemoval::InvE *coupLast_*coup_[ 9]);
    }
    else if((ibos[0] == 24 && ibos[1] == 32) ||
	    (ibos[0] == 32 && ibos[1] == 24)    ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[10]);
    else if((ibos[0] == 24 && ibos[1] == 33) ||
	    (ibos[0] == 33 && ibos[1] == 24)    ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[11]);
    else if((ibos[0] == 34 && ibos[1] == 22) ||
	    (ibos[0] == 22 && ibos[1] == 34)    ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[12]);
    else
      assert(false);
  }
  else if(ih == 38) {
    norm(UnitRemoval::InvE *coupLast_*coup_[13]);
  }
  else
    assert(false);
}
