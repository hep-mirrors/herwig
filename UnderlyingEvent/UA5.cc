#include "UA5.h"
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/PDT/DecayMode.h>
// #include <ThePEG/Handlers/Handler.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Handlers/DecayHandler.h>
#include "Herwig++/Hadronization/Cluster.h"
#include "Herwig++/Utilities/HwDebug.h"
using namespace std;
using namespace ThePEG;
using namespace Herwig;

// Default constructor
UA5Handler::UA5Handler() : N1(9.11), N2(0.115), N3(-9.5),K1(0.029),K2(-0.104),
			   M1(0.4), M2(2.), P1(5.2),P2(3.0),P3(5.2),probSoft(1.0),
			   enhanceCM(1.) {}  
// Copy constructor
UA5Handler::UA5Handler(const UA5Handler &h) :
  MultipleInteractionHandler(h),
  globalParams(h.globalParams),
  clusterFissioner(h.clusterFissioner),
  clusterDecayer(h.clusterDecayer),
  split(h.split),
  N1(h.N1), N2(h.N2), N3(h.N3), K1(h.K1), K2(h.K2),
  M1(h.M1), M2(h.M2), P1(h.P1), P2(h.P2), P3(h.P3),
  probSoft(h.probSoft), enhanceCM(h.enhanceCM)
{}

// Destructor, do nothing
UA5Handler::~UA5Handler() {}

// Saving things into run file
void UA5Handler::persistentOutput(PersistentOStream &os) const {
  os << globalParams << clusterFissioner << clusterDecayer
     << split
     << N1 << N2 << N3 << K1 << K2 << M1 << M2 << P1
     << P2 << P3 << probSoft << enhanceCM;
}

// Reading them back in, in the same order
void UA5Handler::persistentInput(PersistentIStream &is, int) {
  is >> globalParams >> clusterFissioner >> clusterDecayer
     >> split
     >> N1 >> N2 >> N3 >> K1 >> K2 >> M1 >> M2 >> P1 
     >> P2 >> P3 >> probSoft >> enhanceCM;
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
			      &Herwig::UA5Handler::globalParams,
			      false,false,true,false);

  static Reference<UA5Handler,ClusterFissioner>
    interfaceClusterFissioner("ClusterFissioner",
			      "A reference to the ClusterFissioner object",
			      &Herwig::UA5Handler::clusterFissioner,
			      false,false,true,false);

  static Reference<UA5Handler,ClusterDecayer>
    interfaceClusterDecayer("ClusterDecayer",
			    "A reference to the ClusterDecayer object",
			    &Herwig::UA5Handler::clusterDecayer,
			    false,false,true,false);

  static Reference<UA5Handler,PartonSplitter>
    interfaceSplitter("Splitter", "A reference to the PartonSplitter object",
		      &Herwig::UA5Handler::split, false, false, true, false);

  static Parameter<UA5Handler,double>
    interfaceN1("N1", "Parameter N1 in the mean charge multiplicity",
                &Herwig::UA5Handler::N1,9.11,0.,1000.,false,false,false);

  static Parameter<UA5Handler,double>
    interfaceN2("N2", "Parameter N2 in the mean charge multiplicity",
                &Herwig::UA5Handler::N2,0.115,0.,1000.,false,false,false);

  static Parameter<UA5Handler,double>
    interfaceN3("N3", "Parameter N3 in the mean charge multiplicity",
                &Herwig::UA5Handler::N3,-9.5,-1000.,1000.,false,false,false);

  static Parameter<UA5Handler,double>
    interfaceK1("K1", 
		"Parameter K1 used to generate the multiplicity distribution",
                &Herwig::UA5Handler::K1,0.029,0.,1000.,false,false,false);

  static Parameter<UA5Handler,double>
    interfaceK2("K2", 
		"Parameter K2 used to generate the multiplicity distribution",
                &Herwig::UA5Handler::K2,-0.104,-100.,100.,false,false,false);

  static Parameter<UA5Handler,double>
    interfaceM1("M1", "Parameter M1 used to generate soft cluster mass",
                &Herwig::UA5Handler::M1,0.4,0.,100.,false,false,false);

  static Parameter<UA5Handler,double>
    interfaceM2("M2", "Parameter M2 used to generate soft cluster mass",
                &Herwig::UA5Handler::M2,2.0,0.,100.,false,false,false);  

   static Parameter<UA5Handler,double>
    interfaceP1("P1", "Slope used to generate the pt of the u,d soft clusters",
                &Herwig::UA5Handler::P1,5.2,0.,100.,false,false,false); 
   
   static Parameter<UA5Handler,double>
    interfaceP2("P2", "Slope used to generate the pt of the s,c soft clusters",
                &Herwig::UA5Handler::P2,3.0,0.,100.,false,false,false);    

   static Parameter<UA5Handler,double>
    interfaceP3("P3", "Slope used to generate the pt of the qq soft clusters",
                &Herwig::UA5Handler::P3,5.2,0.,100.,false,false,false); 

   static Parameter<UA5Handler,double>
    interfaceProbSoft("Prob", "Probability of underlying event",
		      &Herwig::UA5Handler::probSoft,1.0,0.,1.,false,false,false); 

   static Parameter<UA5Handler,double>
    interfaceEnhance("Enhance", 
		     "Enhancement of the CM energy in the mult distribution",
		     &Herwig::UA5Handler::enhanceCM,1.0,0.,10.,false,false,false);    
}

// This is the routine that is called to start the algorithm. 
void UA5Handler::handle(EventHandler &ch, const tPVector &tagged,
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
   * 1) clusterFissioner->produceHadron(id1,id2,LorentzMomentum &a,
   *                                     LorentzPoint &b)
   *       This function takes the two ids given and finds the lightest hadron
   *       of those flavours. It creates the hadron and sets its momentum to
   *       a and its labVertex to b. This returns a PPair (pair of particles)
   *       where the first element is the hadron and the second is the quark
   *       drawn from the vacuum to create the hadron.
   * 2) clusterFissioner->produceCluster(q,id,LorentzMomentum &a, 
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
   * 3) clusterDecayer->decayIntoTwoHadrons(ClusterPtr cluster)
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
  // Static data values, shouldn't need adjusting...~1 fm
  static double vclx = 4e-12;
  static double vcly = 4e-12;
  static double vclz = 4e-12;
  static double vclt = 4e-12;

  StepPtr step = ch.newStep();

  tPPtr clu[2];
  tPVector::const_iterator it;

  int i = 0;
  int ntry = 0;
  int netc = 0;
  bool multiplicityNeeded = true;
  bool newCluster = true;
  long id1,id2,id3;

  // Lets find the clusters
  for(it = tagged.begin(); it!=tagged.end(); it++) 
     if((*it)->id() == ExtraParticleID::Cluster) clu[i++] = step->copyParticle(*it);
  
  // and their cm
  Lorentz5Momentum cm = clu[0]->momentum() + clu[1]->momentum();

  // Now initialize charge information
  int theMult;
  tPVector tag, hadrons, stable;
  for(i = 0; i<2; i++) {
     ClusterPtr cluster = dynamic_ptr_cast<Cluster>(clu[i]);
     netc += (cluster->parents()[0]->data().charge() - cluster->particle(0)->data().charge() + cluster->particle(1)->data().charge());
  }
  if(netc == 0) id3 = rndbool() ? ParticleID::u : ParticleID::d;
  else if(netc == -1) id3 = ParticleID::u;
  else id3 = ParticleID::d; 
  int numberCluster = 0;
  double sumMasses = 0.;
  int modCharge = 0;
  int netCharge = 0.; 

  // Loop until we find a match to the charge multiplicity
  while(multiplicityNeeded && ntry < max_tries) {
     if(ntry % 50 == 0) theMult = multiplicity(sqrt(cm.m2()));
     ntry++;
     numberCluster = 0;
     tag.push_back(clu[0]);
     tag.push_back(clu[1]);
     sumMasses = 0.;
     modCharge = 0;
     netCharge = netc;
     newCluster = true;
     while(newCluster) {
        // Choose new constituents
        if(numberCluster < 2) {
           id1 = dynamic_ptr_cast<Cluster>(tag[numberCluster])->particle(0)->id();
           id2 = dynamic_ptr_cast<Cluster>(tag[numberCluster])->particle(1)->id();
        } else {
           id1 = -id2;
           if(numberCluster == 2) id1 = id3;
           id2 = rndbool() ? ParticleID::ubar ? ParticleID::dbar;
        }
        PPtr part1 = new_ptr(Particle(id1));
        PPtr part2 = new_ptr(Particle(id2));

        // Fortran code uses HWREXP(2/M2), this is P(x) = int_0^x dx' exp(- 2 x'/M2) -> x = log(1-R)*M2/2
        ClusterPtr cluster = new_ptr(Cluster(part1,part2));
        cluster->momentum()[3] = getParticleData(id1)->mass() + getParticleData(id2)->mass() + M1 + log(1.-rnd())*M2/2.;
        cluster->momentum()[4] = cluster->momentum()[3];
        cluster->momentum()[0] = cluster->momentum()[1] = cluster->momentum()[2] = 0.;
       
        // Now the gaussian distribution of the x,y,z components, and a time component given by
        // sqrt(vx^2+vy^2+vz^2) - vclt*log(r)
        LorentzPoint &point = cluster->vertex();
        point[0] = gaussDistribution(0.,vclx);
        point[1] = gaussDistribution(0.,vcly);
        point[2] = gaussDistribution(0.,vclz);
        point[3] = sqrt(sqr(point[0])+sqr(point[1])+sqr(point[2])) - vclt*log(rnd());
        
        // Now need to measure displacement relative to soft cm (TODO:)

        // Now we decay the clusters into hadrons
        PPair products = clusterDecayer->decayIntoTwoHadrons(cluster);
        if(numberCluster >= 2) tag.push_back(cluster);
        numberCluster++;
        int newHads;
        if(products.first == PPtr() || products.second == PPtr()) {
           products = clusterFissioner->produceHadron(id1,id2, cluster->momentum(), cluster->vertex());
           hadrons.push_back(products.first);
           newHads = 1;
        } else {
           hadrons.push_back(products.first);
           hadrons.push_back(products.second);
           newHads = 2;
        }
        sumMasses += cluster->momentum()[4];
      
        // TODO: Now we decays hadrons into stable particles. Count charged multiplicity
        
        // Now check which loop to do next
        if(netCharge < theMult) continue;
        else if(netCharge > theMul) {
           newCluster = false; 
           tag.clear();
           hadrons.clear();
           stable.clear(); 
        } else { newCluster = false; multiplicityNeeded = false; }
     }
     // Now do some checks on physical boundaries
     if(!multiplicityNeeded) { 
        if(sumMasses > cm) multiplicityNeeded = true;
        if(numberClusters == 0) multiplicityNeeded = true;
     } 
  }
  // TODO: Now generate momentum and add to event record
}

double UA5Handler::gaussDistribution(double mean, double stdev) {
   double x = rnd();
   x = sqrt(-2.*log(x));
   double c,s,cs;
   while(true) {
      c = 2.*rnd()-1.;
      s = 2.*rnd()-1.;
      cs = c*c+s*s;
      if(cs <= 1.) break;
   }
   double qt = x/cs;
   x = (c*c-s*s)*qt;
   return mean + stdev*x;
}

double UA5Handler::meanMultiplicity(Energy E) {
  return N1*pow(E,N2)+N3;
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
  double rk = K1*alogs+K2;
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
  double rn = rnd();
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


 
