#include "UA5.h"
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Handlers/CollisionHandler.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Handlers/DecayHandler.h>

using namespace std;
using namespace ThePEG;
using namespace Herwig;

// Default constructor
UA5Handler::UA5Handler() {}

// Copy constructor
UA5Handler::UA5Handler(const UA5Handler &h) :
  _globalParams(h._globalParams),
  _clusterFissioner(h._clusterFissioner),
  _clusterDecayer(h._clusterDecayer) {}

// Destructor, do nothing
UA5Handler::~UA5Handler() {}

// Saving things into run file
void UA5Handler::persistentOutput(PersistentOStream &os) const {
  os << _globalParams << _clusterFissioner << _clusterDecayer;
}

// Reading them back in, in the same order
void UA5Handler::persistentInput(PersistentIStream &is, int) {
  is >> _globalParams >> _clusterFissioner >> _clusterDecayer;
}

// We must define this static member for ThePEG
ClassDescription<UA5Handler> UA5Handler::initUA5Handler;

// This routine is for reading in from the in file
void UA5Handler::Init() {
  static ClassDocumentation<UA5Handler> documentation
    ("This is the simple UA5 model for the underlying event.");
  
  static Reference<UA5Handler,GlobalParameters>
    interfaceGlobalParameters("GlobalParameters",
			      "A reference to the GlobalParameters object",
			      &Herwig::UA5Handler::_globalParams,
			      false,false,true,false);

  static Reference<UA5Handler,ClusterFissioner>
    interfaceClusterFissioner("ClusterFissioner",
			      "A reference to the ClusterFissioner object",
			      &Herwig::UA5Handler::_clusterFissioner,
			      false,false,true,false);

  static Reference<UA5Handler,ClusterDecayer>
    interfaceClusterDecayer("ClusterDecayer",
			    "A reference to the ClusterDecayer object",
			    &Herwig::UA5Handler::_clusterDecayer,
			    false,false,true,false);
}

// This is the routine that is called to start the algorithm. 
void UA5Handler::handle(PartialCollisionHandler &ch, const tPVector &tagged,
			const Hint &hint) throw(Veto,Stop,Exception) {

  // A few notes:
  // to create a cluster: new_ptr(Cluster(PPtr,PPtr)) where the two
  // PPtr are the consituent partons. It can also be created with three
  // consituents new_ptr(Cluster(PPtr,PPtr,PPtr)); This produces a ClusterPtr
  // To decay a hadron into its products based on the branching ratio's use
  // the decayHadron routine provided in this class. 
}

ParticleVector UA5Handler::decayHadron(tPPtr &had) const 
  throw(Veto,Exception) {
  tDMPtr dm = had->data().selectMode(*had);
  if(!dm) throw DecHdlNoDecayMode(had->data());
  if(!dm->decayer()) throw DecHdlNoDecayer(had->data(), *dm);
  try {
    ParticleVector rval =dm->decayer()->decay(*dm,*had); 
    if(!rval.empty()) {
      had->decayMode(dm);
      had->scale(0.0*GeV2);
      return rval;
    }
  } catch(DecHdlChildFail) {
    throw;
  } catch(Veto) {}
  return ParticleVector();
}

// This is all just administrative functions for ThePEG structure
IBPtr UA5Handler::clone() const { return new_ptr(*this); }
  
IBPtr UA5Handler::fullclone() const { return clone(); }

void UA5Handler::doupdate() throw(UpdateException) {
  MultipleInteractionHandler::doupdate();
}

void UA5Handler::doinit() throw(InitException) {
  MultipleInteractionHandler::doinit();
}

void UA5Handler::dofinish() {
  MultipleInteractionHandler::dofinish();
}

void UA5Handler::rebind(const TranslationMap &trans) throw(RebindException) {
  MultipleInteractionHandler::rebind(trans);
}

IVector UA5Handler::getReferences() {
  IVector ref = MultipleInteractionHandler::getReferences();
  return ref;
}


 
