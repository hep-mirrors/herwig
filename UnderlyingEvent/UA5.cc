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

  /*****
   * A few notes:
   * to create a cluster: new_ptr(Cluster(PPtr,PPtr)) where the two
   * PPtr are the consituent partons. It can also be created with three
   * consituents new_ptr(Cluster(PPtr,PPtr,PPtr)); This produces a ClusterPtr
   * To decay a hadron into its products based on the branching ratio's use
   * the decayHadron routine provided in this class. 
   *
   * There is three ways to decay a cluster, they all return a PPair:
   * 1) _clusterFissioner->produceHadron(id1,id2,LorentzMomentum &a,
   *                                     LorentzPoint &b)
   *       This function takes the two ids given and finds the lightest hadron
   *       of those flavours. It creates the hadron and sets its momentum to
   *       a and its labVertex to b. This returns a PPair (pair of particles)
   *       where the first element is the hadron and the second is the quark
   *       drawn from the vacuum to create the hadron.
   * 2) _clusterFissioner->produceCluster(q,id,LorentzMomentum &a, 
   *                                      LorentzPoint &b, LorentzMomentum &c,
   *                                      LorentzMomentum &d)
   *       This function takes a consituent quark q and an id and produces
   *       a cluster from these. q is a PPtr and should be already in existence
   *       while id is a long indicating the flavour to draw from the vacuum.
   *       The cluster has momentum a and labVertex b. The consituent given
   *       by q is given momentum c (note that this doesn't change the original
   *       q's momentum, just sets it momentum of the copy inside the cluster)
   *       and the quark drawn from the vacuum is given momentum d. This 
   *       returns a PPair where the first element is the cluster and the 
   *       second is the quark drawn from the vacuum.
   * 3) _clusterDecayer->decayIntoTwoHadrons(ClusterPtr cluster)
   *       This function does what it says. It takes a cluster and returns
   *       two hadrons that this cluster decays into as given by the method
   *       of cluster decays (currently the new method created by Phil 
   *       Stephens). If the constituents are non-perturbative (not directly
   *       from the shower) then the clusters are decayed isotropically. If 
   *       Otherwise whichever cluster corresponds to the perturbative 
   *       constituent may have its momentum smeared in the direction of the
   *       constituents momenta.
   *
   * A few other ThePEG pointers:
   *  To find the colour partner to a particle use particle.colourNeighbour()
   *   or particle.antiColourNeighbour(). 
   *  To get the current step use ch->currentStep(), this allows access to the
   *   current set of particles. 
   *  To create a new step use ch->newStep();
   *  Steps are where changes to the EventRecord show up. If you don't use
   *   routines like addParticle, addIntermediate, addDecayProduct then the
   *   particles created won't appear in the EventRecord.
   *****/
}

ParticleVector UA5Handler::decayHadron(tPPtr &had) const 
  throw(Veto,Exception) {
  tDMPtr dm = had->data().selectMode(*had);
  if(!dm) throw DecHdlNoDecayMode(had->data());
  if(!dm->decayer()) throw DecHdlNoDecayer(had->data(), *dm);
  try {
    ParticleVector rval = dm->decayer()->decay(*dm,*had); 
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


 
