// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHWWHVertex class.
//

#include "LHWWHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LHWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coup,GeV);
}

void LHWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coup,GeV);
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHWWHVertex,VVSVertex>
describeHerwigLHWWHVertex("Herwig::LHWWHVertex", "HwLHModel.so");

void LHWWHVertex::Init() {

  static ClassDocumentation<LHWWHVertex> documentation
    ("The LHWWHVertex class implements the coupling of two electroweak"
     " gauge bosons to a Higgs boson in the Little Higgs Model including the "
     "additional heavy photon, Z and W bosons and the triplet Higgs bosons.");

}

LHWWHVertex::LHWWHVertex() 
  : _couplast(0.), _q2last(0.*GeV2) {
  // order in the couplings
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void LHWWHVertex::doinit() {
  // W_L W_L H
  addToList(  24,  -24,    25);
  // Z_L Z_L H
  addToList(  23,   23,    25);
  // W_L W_H H
  addToList(  24,  -34,    25);
  addToList(  34,  -24,    25);
  // Z_L A_H H
  addToList(  23,   32,    25);
  // W_H W_H H
  addToList(  34,  -34,    25);
  // Z_H Z_H H
  addToList(  33,   33,    25);
  // A_H A_H H
  addToList(  32,   32,    25);
  // Z_H Z_L H
  addToList(  23,   33,    25);
  // Z_H A_H H
  addToList(  33,   32,    25);
  // W_L W_L Phi0
  addToList(  24,  -24,    35);
  // W_L W_H Phi0
  addToList(  24,  -34,    35);
  addToList(  34,  -24,    35);
  // Z_L Z_L Phi0
  addToList(  23,   23,    35);
  // Z_L Z_H Phi0
  addToList(  23,   33,    35);
  // W_H W_H Phi0
  addToList(  34,  -34,    35);
  // Z_H Z_H Phi0
  addToList(  33,   33,    35);
  // A_H Z_H Phi0
  addToList(  32,   33,    35);
  // A_H Z_L Phi0
  addToList(  32,   23,    35);
  // A_H A_H Phi0
  addToList(  32,   32,    35);
  // W_L Z_L Phi-
  addToList(  24,   23,   -37);
  addToList( -24,   23,    37);
  // W_L A_H Phi-
  addToList(  24,   32,   -37);
  addToList( -24,   32,    37);
  // W_L Z_H Phi-
  addToList(  24,   33,   -37);
  addToList( -24,   33,    37);
  // W_H Z_L Phi-
  addToList(  34,   23,   -37);
  addToList( -34,   23,    37);
  // W_H A_H Phi-
  addToList(  34,   32,   -37);
  addToList( -34,   32,    37);
  // W_H Z_H Phi-
  addToList(  34,   33,   -37);
  addToList( -34,   33,    37);
  // W_L W_L Phi--
  addToList(  24,   24,   -38);
  addToList( -24,  -24,    38);
  // W_H W_H Phi--
  addToList(  34,   34,   -38);
  addToList( -34,  -34,    38);
  // W_L W_H Phi--
  addToList(  24,   34,   -38);
  addToList( -24,  -34,    38);
  // model
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) throw InitException() << "Must be using the LHModel "
				   << " in LHWWHVertex::doinit()"
				   << Exception::runerror;
  // base class
  VVSVertex::doinit();
  // calculate the couplings for the different combinations of particles
  double sw(sqrt(sin2ThetaW())),cw(sqrt(1.-sin2ThetaW()));
  Energy fact = getParticleData(ParticleID::Wplus)->mass()/sw;
  double vf(sqr(model->vev()/model->f()));
  double vr(model->vevPrime()/model->vev());
  double r2(sqrt(2.));
  double s (model->sinTheta()     ),c (model->cosTheta()     );
  double sp(model->sinThetaPrime()),cp(model->cosThetaPrime());
  double sPlus(model->sinTheta0());
  double s0(model->sinTheta0());
  _coup.resize(27);
  // couplings to SM higgs
  _coup[ 0] = fact        *(1.-vf/3.+0.5*vf* sqr(sqr(c)-sqr(s))
			    -0.5*sqr(s0)-2.*r2*s0*vr);
  _coup[ 1] = fact/sqr(cw)*(1.-vf/3.-0.5*vf*(sqr(sqr(c)-sqr(s))+5.*sqr(sqr(cp)-sqr(sp)))
			    -0.5*sqr(s0)+4.*r2*s0*vr);
  _coup[ 2] =-fact;
  _coup[ 3] =-fact;
  _coup[ 4] =-fact*sqr(sw/cw);
  _coup[ 5] =-fact*   0.5*(sqr(c)-sqr(s))/s/c;
  _coup[ 6] =-fact/cw*0.5*(sqr(c)-sqr(s))/s/c;
  _coup[ 7] =-fact/sqr(cw)*sw*0.5*(sqr(cp)-sqr(sp))/sp/cp;
  _coup[ 8] =-fact/cw*sw*0.5/s/c/sp/cp*(sqr(c*sp)+sqr(s*cp));
  _coup[ 9] =-fact*(s0-2.*r2*vr);
  _coup[10] = fact*(s0-2.*r2*vr);
  _coup[11] = fact*(s0-2.*r2*vr)*0.5*(sqr(c)-sqr(s))/s/c;
  _coup[12] =-fact/sqr(cw)*(s0-4.*r2*vr);
  _coup[13] = fact*(s0+sqr(sqr(c)-sqr(s))/sqr(s*c)*r2*vr);
  _coup[14] = fact/cw*0.5*(sqr(c)-sqr(s))/s/c*(s0-4.*r2*vr);
  _coup[15] = fact*sw/sqr(cw)*0.5*(sqr(cp)-sqr(sp))/sp/cp*(s0-4.*r2*vr);
  _coup[16] = fact*sw/cw*0.5/s/c/sp/cp*(s0*(sqr(c*sp)+sqr(s*cp))
					+2.*r2*(sqr(c)-sqr(s))*(sqr(cp)-sqr(sp))*vr);
  _coup[17] = fact*sqr(sw/cw)*(s0+r2*vr*sqr(sqr(cp)-sqr(sp))/sqr(sp*cp));
  _coup[18] =-2.*fact/cw*vr;
  _coup[19] = fact/cw*(sqr(c)-sqr(s))/s/c*vr;
  _coup[20] =-fact*sw/cw*0.5*(sqr(cp)-sqr(sp))/sp/cp*(sPlus-4.*vr);
  _coup[21] =-fact*sw/cw*(sqr(c*cp)+sqr(s*sp))/s/c/sp/cp*vr;
  _coup[22] = fact*(sqr(c)-sqr(s))/s/c*vr;
  _coup[23] =-fact*(pow(c,4)+pow(s,4))/sqr(s*c)*vr;
  _coup[24] = fact*4.*vr;
  _coup[25] = fact*2.*(pow(c,4)+pow(s,4))/sqr(s*c)*vr;
  _coup[26] =-fact*2.*vr*(sqr(c)-sqr(s))/s/c;
}

void LHWWHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = electroMagneticCoupling(q2);
    _q2last=q2;
  }
  int ih = abs(c->id());
  long int ibos[2]={abs(a->id()),abs(b->id())};
  if(ih==25) {
    if     ( ibos[0]==24&&ibos[1]==24     ) norm(UnitRemoval::InvE *_couplast*_coup[0]);
    else if( ibos[0]==23&&ibos[1]==23     ) norm(UnitRemoval::InvE *_couplast*_coup[1]);
    else if( ibos[0]==34&&ibos[1]==34     ) norm(UnitRemoval::InvE *_couplast*_coup[2]);
    else if( ibos[0]==33&&ibos[1]==33     ) norm(UnitRemoval::InvE *_couplast*_coup[3]);
    else if( ibos[0]==32&&ibos[1]==32     ) norm(UnitRemoval::InvE *_couplast*_coup[4]);
    else if((ibos[0]==24&&ibos[1]==34) ||
	    (ibos[0]==34&&ibos[1]==24)    ) norm(UnitRemoval::InvE *_couplast*_coup[5]);
    else if((ibos[0]==23&&ibos[1]==33) ||
	    (ibos[0]==33&&ibos[1]==23)    ) norm(UnitRemoval::InvE *_couplast*_coup[6]);
    else if((ibos[0]==23&&ibos[1]==32) ||
	    (ibos[0]==32&&ibos[1]==23)    ) norm(UnitRemoval::InvE *_couplast*_coup[7]);
    else if((ibos[0]==33&&ibos[1]==32) ||
	    (ibos[0]==32&&ibos[1]==33)    ) norm(UnitRemoval::InvE *_couplast*_coup[8]);
    else assert(false);
  }
  else if(ih==35) {
    if     ( ibos[0]==24&&ibos[1]==24     ) norm(UnitRemoval::InvE *_couplast*_coup[ 9]);
    else if( ibos[0]==34&&ibos[1]==34     ) norm(UnitRemoval::InvE *_couplast*_coup[10]);
    else if((ibos[0]==24&&ibos[1]==34) ||
	    (ibos[0]==34&&ibos[1]==24)    ) norm(UnitRemoval::InvE *_couplast*_coup[11]);
    else if( ibos[0]==23&&ibos[1]==23     ) norm(UnitRemoval::InvE *_couplast*_coup[12]);
    else if( ibos[0]==33&&ibos[1]==33     ) norm(UnitRemoval::InvE *_couplast*_coup[13]);
    else if((ibos[0]==23&&ibos[1]==33) ||
	    (ibos[0]==33&&ibos[1]==23)    ) norm(UnitRemoval::InvE *_couplast*_coup[14]);
    else if((ibos[0]==23&&ibos[1]==32) ||
	    (ibos[0]==32&&ibos[1]==23)    ) norm(UnitRemoval::InvE *_couplast*_coup[15]);
    else if((ibos[0]==33&&ibos[1]==32) ||
	    (ibos[0]==32&&ibos[1]==33)    ) norm(UnitRemoval::InvE *_couplast*_coup[16]);
    else if((ibos[0]==32&&ibos[1]==32)    ) norm(UnitRemoval::InvE *_couplast*_coup[17]);
    else assert(false);
  }
  else if(ih==37) {
    if     ((ibos[0]==24&&ibos[1]==23) ||
	    (ibos[0]==23&&ibos[1]==24)    ) norm(UnitRemoval::InvE *_couplast*_coup[18]);
    else if((ibos[0]==34&&ibos[1]==23) ||
	    (ibos[0]==23&&ibos[1]==34)    ) norm(UnitRemoval::InvE *_couplast*_coup[19]);
    else if((ibos[0]==24&&ibos[1]==32) ||
	    (ibos[0]==32&&ibos[1]==24)    ) norm(UnitRemoval::InvE *_couplast*_coup[20]);
    else if((ibos[0]==34&&ibos[1]==32) ||
	    (ibos[0]==32&&ibos[1]==34)    ) norm(UnitRemoval::InvE *_couplast*_coup[21]);
    else if((ibos[0]==24&&ibos[1]==33) ||
	    (ibos[0]==33&&ibos[1]==24)    ) norm(UnitRemoval::InvE *_couplast*_coup[22]);
    else if((ibos[0]==34&&ibos[1]==33) ||
	    (ibos[0]==33&&ibos[1]==34)    ) norm(UnitRemoval::InvE *_couplast*_coup[23]);
    else assert(false);
  }
  else if(ih==38) {
    if     ((ibos[0]==24&&ibos[1]==24)    ) norm(UnitRemoval::InvE *_couplast*_coup[24]);
    else if((ibos[0]==34&&ibos[1]==34)    ) norm(UnitRemoval::InvE *_couplast*_coup[24]);
    else if((ibos[0]==34&&ibos[1]==24) ||
	    (ibos[0]==24&&ibos[1]==34)    ) norm(UnitRemoval::InvE *_couplast*_coup[24]);
    else assert(false);
  }
  else assert(false);
}
