#include "VBFNLOCommonBlocks.h"
#include "VBFNLOAmplitudeqq2hqqg.h"
#include "VBFNLOMEqq2hqqg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/Throw.h"
 
#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"

using namespace Herwig;


VBFNLOAmplitudeqq2hqqg::VBFNLOAmplitudeqq2hqqg():VBFNLOAmplitudePP2hJetJetJet(){
}

VBFNLOAmplitudeqq2hqqg::~VBFNLOAmplitudeqq2hqqg(){}

IBPtr VBFNLOAmplitudeqq2hqqg::clone() const {
  return new_ptr(*this);
}

IBPtr VBFNLOAmplitudeqq2hqqg::fullclone() const {
  return new_ptr(*this);
}

bool VBFNLOAmplitudeqq2hqqg::allowedProcess(const PDVector& data) const {
  PDPtr Wplus = getParticleData(ParticleID::Wplus); 
  PDPtr Wminus = getParticleData(ParticleID::Wminus); 
  PDPtr Z0 = getParticleData(ParticleID::Z0); 

  if ( data[0]->id() == 21 || data[1]->id() == 21 ) return false;
  if ( theCurrent == neutral
       && ( data[0] != data[2]
	    || data[1] != data[3] ) ) return false;
  if ( theCurrent == charged
       &&  ( data[0] != SU2Helper::SU2CC(data[2])
	     || data[1] != SU2Helper::SU2CC(data[3]) ) ) return false;
  if ( ( data[0]->id() < 0 && theIncoming1)
       || (data[1]->id() < 0 && theIncoming2) ) return false;
  if ( ( data[0]->id() > 0 && !theIncoming1)
       || (data[1]->id() > 0 && !theIncoming2) ) return false;
  return true;
}

Ptr<MatchboxMEBase>::ptr VBFNLOAmplitudeqq2hqqg::makeME(const PDVector& ) const {
  return new_ptr(VBFNLOMEqq2hqqg(theCurrent, theIncoming1, theIncoming2, theDecayChannel, theNarrowWidth));
}

void VBFNLOAmplitudeqq2hqqg::prepareMomenta(double (& pbar)[14][4],
					double (&qbar) [5],int) const{
  L5MomToDouble( (mePartonData()[0]->id() > 0 ? 
  		  meMomenta()[0]  : 
  		  meMomenta()[2]), &pbar[0][0]);

  L5MomToDouble( (mePartonData()[2]->id() > 0 ? 
  		  meMomenta()[2] : 
  		  meMomenta()[0]), &pbar[1][0]);

  L5MomToDouble( (mePartonData()[1]->id() > 0 ? 
  		  meMomenta()[1] : 
  		  meMomenta()[3]), &pbar[2][0]);

  L5MomToDouble( (mePartonData()[3]->id() > 0 ? 
  		  meMomenta()[3] : 
  		  meMomenta()[1]), &pbar[3][0]);

  L5MomToDouble( (meMomenta()[4]), &qbar[0]);
  qbar[4]=0;

  pbar[4][0]=0;
  pbar[4][1]=0;
  pbar[4][2]=0;
  pbar[4][3]=0;

  if (NDecayProducts() == 1) {
    L5MomToDouble( (meMomenta()[5]), &pbar[5][0]);
  }
  else if (NDecayProducts() == 2) {
    L5MomToDouble( (meMomenta()[5]+meMomenta()[6]), &pbar[5][0]);
  }
  else if (NDecayProducts() == 4) {
    L5MomToDouble( (meMomenta()[5]+meMomenta()[6]+meMomenta()[7]+meMomenta()[8]), &pbar[5][0]);
  }


  return;
} 

void VBFNLOAmplitudeqq2hqqg::doinit(){ 
  VBFNLOAmplitudeVVJJNeutralBase::doinit();
}

void VBFNLOAmplitudeqq2hqqg::doinitrun(){ 
  VBFNLOAmplitudeVVJJNeutralBase::doinitrun();
}

ClassDescription<VBFNLOAmplitudeqq2hqqg> VBFNLOAmplitudeqq2hqqg::initVBFNLOAmplitudeqq2hqqg;

void VBFNLOAmplitudeqq2hqqg::persistentOutput(PersistentOStream & ) const {
}

void VBFNLOAmplitudeqq2hqqg::persistentInput(PersistentIStream & , int) {
}


void VBFNLOAmplitudeqq2hqqg::Init() {

  static ClassDocumentation<VBFNLOAmplitudeqq2hqqg> documentation
    ("VBFNLOAmplitudeqq2hqqg");

}

