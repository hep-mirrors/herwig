#include "UA5.h"
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Handlers/CollisionHandler.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Handlers/DecayHandler.h>
#include "Herwig++/Hadronization/Cluster.h"
#include "Herwig++/Utilities/HwDebug.h"
using namespace std;
using namespace ThePEG;
using namespace Herwig;

// Default constructor
UA5Handler::UA5Handler() {}

// Copy constructor
UA5Handler::UA5Handler(const UA5Handler &h) :
  _globalParams(h._globalParams),
  _clusterFissioner(h._clusterFissioner),
  _clusterDecayer(h._clusterDecayer),
  _split(h._split),
  _N1(h._N1), _N2(h._N2), _N3(h._N3), _K1(h._K1), _K2(h._K2),
  _M1(h._M1), _M2(h._M2), _P1(h._P1), _P2(h._P2), _P3(h._P3),
  _probSoft(h._probSoft), _enhanceCM(h._enhanceCM)
{}

// Destructor, do nothing
UA5Handler::~UA5Handler() {}

// Saving things into run file
void UA5Handler::persistentOutput(PersistentOStream &os) const {
  os << _globalParams << _clusterFissioner << _clusterDecayer
     << _split
     << _N1 << _N2 << _N3 << _K1 << _K2 << _M1 << _M2 << _P1
     << _P2 << _P3 << _probSoft << _enhanceCM;
}

// Reading them back in, in the same order
void UA5Handler::persistentInput(PersistentIStream &is, int) {
  is >> _globalParams >> _clusterFissioner >> _clusterDecayer
     >> _split
     >> _N1 >> _N2 >> _N3 >> _K1 >> _K2 >> _M1 >> _M2 >> _P1 
     >> _P2 >> _P3 >> _probSoft >> _enhanceCM;
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

  static Reference<UA5Handler,PartonSplitter>
    interfaceSplitter("Splitter", "A reference to the PartonSplitter object",
		      &Herwig::UA5Handler::_split, false, false, true, false);

  static Parameter<UA5Handler,double>
    interfaceN1("N1", "Parameter N1 in the mean charge multiplicity",
                &Herwig::UA5Handler::_N1,9.11,0.,1000.);

  static Parameter<UA5Handler,double>
    interfaceN2("N2", "Parameter N2 in the mean charge multiplicity",
                &Herwig::UA5Handler::_N2,0.115,0.,1000.);

  static Parameter<UA5Handler,double>
    interfaceN3("N3", "Parameter N3 in the mean charge multiplicity",
                &Herwig::UA5Handler::_N3,-9.5,-1000.,1000.);

  static Parameter<UA5Handler,double>
    interfaceK1("K1", 
		"Parameter K1 used to generate the multiplicity distribution",
                &Herwig::UA5Handler::_K1,0.029,0.,1000.);

  static Parameter<UA5Handler,double>
    interfaceK2("K2", 
		"Parameter K2 used to generate the multiplicity distribution",
                &Herwig::UA5Handler::_K2,-0.104,-100.,100.);

  static Parameter<UA5Handler,double>
    interfaceM1("M1", "Parameter M1 used to generate soft cluster mass",
                &Herwig::UA5Handler::_M1,0.4,0.,100.);

  static Parameter<UA5Handler,double>
    interfaceM2("M2", "Parameter M2 used to generate soft cluster mass",
                &Herwig::UA5Handler::_M2,2.0,0.,100.);  

   static Parameter<UA5Handler,double>
    interfaceP1("P1", "Slope used to generate the pt of the u,d soft clusters",
                &Herwig::UA5Handler::_P1,5.2,0.,100.); 
   
   static Parameter<UA5Handler,double>
    interfaceP2("P2", "Slope used to generate the pt of the s,c soft clusters",
                &Herwig::UA5Handler::_P2,3.0,0.,100.);    

   static Parameter<UA5Handler,double>
    interfaceP3("P3", "Slope used to generate the pt of the qq soft clusters",
                &Herwig::UA5Handler::_P3,5.2,0.,100.); 

   static Parameter<UA5Handler,double>
    interfaceProbSoft("Prob", "Probability of underlying event",
		      &Herwig::UA5Handler::_probSoft,1.0,0.,1.); 

   static Parameter<UA5Handler,double>
    interfaceEnhance("Enhance", 
		     "Enhancement of the CM energy in the mult distribution",
		     &Herwig::UA5Handler::_enhanceCM,1.0,0.,10.);    
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
  cout << "In main code\n";
  // Find beam particles
  PPair beam = ch.currentCollision()->incoming();
  PVector::const_iterator it;
  PPtr rem1,rem2;
  for(it = beam.first->children().begin(); it != beam.first->children().end();
      it++)  if((*it)->children().size()==0) rem1 = *it;
  for(it = beam.second->children().begin(); it !=beam.second->children().end();
      it++)  if((*it)->children().size()==0) rem2 = *it;
  if(!rem1 && !rem2) return;

  // Now form the first two clusters, first split gluons
  ClusterVector clusters;
  StepPtr step = ch.newStep();
  tPVector tag;
  tag.push_back(rem1->colourNeighbour());
  tag.push_back(rem2->colourNeighbour());
  cout << "Calling split\n" << *step;
  _split->split(tag,step);

  ClusterPtr clu1 = new_ptr(Cluster(rem1, rem1->colourNeighbour()));
  ClusterPtr clu2 = new_ptr(Cluster(rem2, rem2->colourNeighbour()));
  Lorentz5Momentum cm = clu1->momentum() + clu2->momentum();
  cout << "rem1 := " << *rem1 << " and partner is " << *rem1->colourNeighbour()
       << endl
       << "rem2 := " << *rem2 << " and partner is " << *rem2->colourNeighbour()
       << endl;
  cout << "Mom: = " << cm << endl;
}

double UA5Handler::meanMultiplicity(Energy E) {
  return _N1*pow(E,_N2)+_N3;
}

double UA5Handler::negativeBinomial(int N, double mean, double ek) {
  if(N < 0) return 0.0;
  double r = mean/ek;
  double rval = pow(1.+r, -ek);
  r /= (1.+r);
  for(int i = 1; i<=N; i++) rval *=  r*(ek+double(i)-1.)/double(i);
  return rval;
}

int UA5Handler::multiplicity(Energy E) {
  double alogs = 2.*log(E);
  double rk = _K1*alogs+_K2;
  if(rk > 1000.) rk = 1000.;
  double ek = 1./rk;
  double mean = meanMultiplicity(E);
  if(mean < 1.) mean = 1.;
  vector<double> dist;
  double sum = 0.0;
  int imax = 0;
  int i;
  for(i = 0; i<500; i++) {
    int N = (i+1)*2;
    dist.push_back(negativeBinomial(N,mean,ek));
    if(dist[i] < 1e-7*sum) break;
    imax = i;
    sum += dist[i];
    dist[i] = sum;
  }
  for(i = 0; i<=imax; i++) dist[i] /= sum;
  double rn = rand();
  for(i = 0; i<=imax; i++) if(rn < dist[i]) break;
  return 2*(i+1);
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


 
