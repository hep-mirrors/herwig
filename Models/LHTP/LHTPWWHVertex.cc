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
  os << ounit(_coup,GeV);
}

void LHTPWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coup,GeV);
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

LHTPWWHVertex::LHTPWWHVertex() : _couplast(0.), _q2last(0.*GeV2) {
  // order in the couplings
  orderInGem(1);
  orderInGs(0);
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


  // Z_H Z_L H --- VIOLATES T-parity
//   addToList(  23,//    33,//     25);



  // Z_H A_H H
  addToList(  33,   32,    25);

  // W_L W_L Phi0 --- VIOLATES T-parity
//   addToList(  24,//   -24,//     35);
//   // Z_L Z_L Phi0
//   addToList(  23,//    23,//     35);


  // Z_L Z_H Phi0
  addToList(  23,   33,    35);

  // W_H W_H Phi0 -- VIOLATES T-parity
//   addToList(  34,//   -34,//     35);


  // A_H Z_L Phi0
  addToList(  32,   23,    35);

  // W_L Z_L Phi- --- VIOLATES T-parity 
//   addToList(  24,//    23,//    -37);
//   addToList( -24,//    23,//     37);

  // W_H Z_L Phi- 
  addToList(  34,   23,   -37);
  addToList( -34,   23,    37);
  // W_L A_H Phi-
  addToList(  24,   32,   -37);
  addToList( -24,   32,    37);

  // W_H A_H Phi- --- VIOLATES T-parity
//   addToList(  34,//    32,//    -37);
//   addToList( -34,//    32,//     37);

  // W_L Z_H Phi- 
  addToList(  24,   33,   -37);
  addToList( -24,   33,    37);

  // W_H Z_H Phi- --- VIOLATES T-parity
//   addToList(  34,//    33,//    -37);
//   addToList( -34,//    33,//     37);

  // W_H A_L Phi- 
  addToList(  34,   22,   -37);
  addToList( -34,   22,    37);
}

void LHTPWWHVertex::doinit() {
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
  double vf(sqr(model->vev()/model->f()));
  double r2(sqrt(2.));
  _coup.resize(18);
  _coup[ 0] = fact        *(1.-vf/3.);
  _coup[ 1] = fact/sqr(cw)*(1.-vf/3.);
  _coup[ 2] =-fact;
  _coup[ 3] =-fact;
  _coup[ 4] =-fact*sqr(sw/cw);
  _coup[ 5] =-fact/cw*sw;
  _coup[ 6] = fact*2.*r2;
  _coup[ 7] =-fact*2.*r2;
  _coup[ 8] = fact/sqr(cw)*4.*r2;
  _coup[ 9] =-fact*sqrt(vf)/r2/cw;
  _coup[10] = fact*sw/sqr(cw)*sqrt(vf)/r2/sw;
  _coup[11] =-fact/cw;
  _coup[12] = fact*sqrt(vf)/6./cw*(1.+2.*sqr(sw));
  _coup[13] =-fact*sqrt(vf)*sw/cw*0.5;
  _coup[14] =-fact*sw/cw*2.;
  _coup[15] = fact*sqrt(vf)*5./6.;
  _coup[16] =-fact*2.;
  _coup[17] =-fact*sqrt(vf)*sw/3.;
}

void LHTPWWHVertex::setCoupling(Energy2 q2,tcPDPtr a,
					      tcPDPtr b, tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = sqr(electroMagneticCoupling(q2));
    _q2last=q2;
  }
  long ih = abs(c->id());
  long ibos[2]={abs(a->id()),abs(b->id())};
  
  if(ih == 25) {
    if( ibos[0] == 24 && ibos[1] == 24) 
      norm(UnitRemoval::InvE *_couplast*_coup[0]);
    else if( ibos[0] == 23 && ibos[1] == 23     ) 
      norm(UnitRemoval::InvE *_couplast*_coup[1]);
    else if( ibos[0] == 34 && ibos[1] == 34     ) 
      norm(UnitRemoval::InvE *_couplast*_coup[2]);
    else if( ibos[0] == 33 && ibos[1] == 33     ) 
      norm(UnitRemoval::InvE *_couplast*_coup[3]);
    else if( ibos[0] == 32 && ibos[1] == 32     ) 
      norm(UnitRemoval::InvE *_couplast*_coup[4]);
    else if((ibos[0] == 33 && ibos[1] == 32) ||
	    (ibos[0] == 32 && ibos[1] == 33)    ) 
      norm(UnitRemoval::InvE *_couplast*_coup[5]);
    else 
      throw HelicityConsistencyError() << "LittleHiggsWWHVertex::setCoupling "
				       << "Invalid particles in WWH Vertex "
				       << a->PDGName() << " " << b->PDGName() << " "
				       << c->PDGName() 
				       << Exception::runerror;
  }
  else if(ih == 35) {
    //    if     ( ibos[0] == 24 && ibos[1] == 24     ) norm(UnitRemoval::InvE *_couplast*_coup[ 6]);
    //    if( ibos[0] == 34 && ibos[1] == 34     ) norm(UnitRemoval::InvE *_couplast*_coup[ 7]);
    //else if( ibos[0] == 23 && ibos[1] == 23     ) norm(UnitRemoval::InvE *_couplast*_coup[ 8]);
    if((ibos[0] == 23 && ibos[1] == 33) ||
       (ibos[0] == 33 && ibos[1] == 23)    ) 
      norm(UnitRemoval::InvE *_couplast*_coup[ 9]);
    else if((ibos[0] == 23 && ibos[1] == 32) ||
	    (ibos[0] == 32 && ibos[1] == 23)    ) 
      norm(UnitRemoval::InvE *_couplast*_coup[10]);
    else 
      throw HelicityConsistencyError() << "LittleHiggsWWHVertex::setCoupling "
				       << "Invalid particles in WWH Vertex " 
				       << a->PDGName() << " " << b->PDGName() << " "
				       << c->PDGName() 
				       << Exception::runerror;
  }
  else if(ih == 37) {
//     if     ((ibos[0] == 24 && ibos[1] == 23) ||
// 	    (ibos[0] == 23 && ibos[1] == 24)    ) norm(UnitRemoval::InvE *_couplast*_coup[11]);
    if((ibos[0] == 34 && ibos[1] == 23) ||
       (ibos[0] == 23 && ibos[1] == 34)    ) 
      norm(UnitRemoval::InvE *_couplast*_coup[12]);
    else if((ibos[0] == 24 && ibos[1] == 32) ||
	    (ibos[0] == 32 && ibos[1] == 24)    ) 
      norm(UnitRemoval::InvE *_couplast*_coup[13]);
    //     if((ibos[0] == 34 && ibos[1] == 32) ||
//        (ibos[0] == 32 && ibos[1] == 34)    ) norm(UnitRemoval::InvE *_couplast*_coup[14]);
    else if((ibos[0] == 24 && ibos[1] == 33) ||
       (ibos[0] == 33 && ibos[1] == 24)    ) 
      norm(UnitRemoval::InvE *_couplast*_coup[15]);
//     else if((ibos[0] == 34 && ibos[1] == 33) ||
// 	    (ibos[0] == 33 && ibos[1] == 34)    ) norm(UnitRemoval::InvE *_couplast*_coup[16]);
    else if((ibos[0] == 34 && ibos[1] == 22) ||
	    (ibos[0] == 22 && ibos[1] == 34)    ) 
      norm(UnitRemoval::InvE *_couplast*_coup[17]);
    else 
      throw HelicityConsistencyError() 
	<< "LittleHiggsWWHVertex::setCoupling "
	<< "Invalid particles in WWH Vertex " 
	<< a->PDGName() << " " << b->PDGName() << " " << c->PDGName() 
	<< Exception::runerror;
  }
}
