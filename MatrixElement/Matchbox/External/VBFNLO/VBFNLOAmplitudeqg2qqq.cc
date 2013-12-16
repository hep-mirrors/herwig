#include "VBFNLOCommonBlocks.h"
#include "VBFNLOAmplitudeqg2hqqq.h"
#include "VBFNLOMEqg2hqqq.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Utilities/Throw.h"
 
#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"

using namespace Herwig;

VBFNLOAmplitudeqg2hqqq::VBFNLOAmplitudeqg2hqqq():VBFNLOAmplitudePP2hJetJetJet(){
}

VBFNLOAmplitudeqg2hqqq::~VBFNLOAmplitudeqg2hqqq(){}

IBPtr VBFNLOAmplitudeqg2hqqq::clone() const {
  return new_ptr(*this);
}

IBPtr VBFNLOAmplitudeqg2hqqq::fullclone() const {
  return new_ptr(*this);
}

Ptr<MatchboxMEBase>::ptr VBFNLOAmplitudeqg2hqqq::makeME(const PDVector& data) const {
  return new_ptr(VBFNLOMEqg2hqqq(theCurrent, theIncoming1, theIncoming2, theDecayChannel, theNarrowWidth, theWhichGluon));
}

void VBFNLOAmplitudeqg2hqqq::prepareMomenta(double (& pbar)[14][4],
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

void VBFNLOAmplitudeqg2hqqq::doinit(){ 
  VBFNLOAmplitudeVVJJNeutralBase::doinit();
}

void VBFNLOAmplitudeqg2hqqq::doinitrun(){ 
  VBFNLOAmplitudeVVJJNeutralBase::doinitrun();
}

ClassDescription<VBFNLOAmplitudeqg2hqqq> VBFNLOAmplitudeqg2hqqq::initVBFNLOAmplitudeqg2hqqq;


void VBFNLOAmplitudeqg2hqqq::persistentOutput(PersistentOStream & os) const {
  os << theWhichGluon;
}

void VBFNLOAmplitudeqg2hqqq::persistentInput(PersistentIStream & is, int) {
  is >> theWhichGluon;
}


void VBFNLOAmplitudeqg2hqqq::Init() {

  static ClassDocumentation<VBFNLOAmplitudeqg2hqqq> documentation
    ("VBFNLOAmplitudeqg2hqqq");

  static Switch<VBFNLOAmplitudeqg2hqqq,int> interfaceWhichGluon
    ("WhichGluon",
     "Set the position of the incoming gluon.",
     &VBFNLOAmplitudeqg2hqqq::theWhichGluon, 0, false, false);  
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
