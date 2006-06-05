#include "UA5.h"
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Repository/UseRandom.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Handlers/DecayHandler.h>
#include "Herwig++/Hadronization/Cluster.h"
#include "Herwig++/Utilities/HwDebug.h"

using namespace std;
using namespace ThePEG;
using namespace Herwig;

const int max_tries = 300;
const unsigned int max_had = 50;

// Default constructor
UA5Handler::UA5Handler() : N1(9.11), N2(0.115), N3(-9.5), K1(0.029), K2(-0.104),
			   M1(0.4*GeV), M2(2.*GeV), P1(5.2), P2(3.0), P3(5.2), probSoft(1.0),
			   enhanceCM(1.) {}  
// Copy constructor
UA5Handler::UA5Handler(const UA5Handler &h) :
  HadronizationHandler(h),
  globalParams(h.globalParams),
  clusterFissioner(h.clusterFissioner),
  clusterDecayer(h.clusterDecayer),
  split(h.split),
  decayer(h.decayer),
  N1(h.N1), N2(h.N2), N3(h.N3), K1(h.K1), K2(h.K2),
  M1(h.M1), M2(h.M2), P1(h.P1), P2(h.P2), P3(h.P3),
  probSoft(h.probSoft), enhanceCM(h.enhanceCM)
{}

// Destructor, do nothing
UA5Handler::~UA5Handler() {}

// Saving things into run file
void UA5Handler::persistentOutput(PersistentOStream &os) const {
  os << globalParams << clusterFissioner << clusterDecayer
     << split << decayer 
     << N1 << N2 << N3 << K1 << K2 << M1 << M2 << P1
     << P2 << P3 << probSoft << enhanceCM;
}

// Reading them back in, in the same order
void UA5Handler::persistentInput(PersistentIStream &is, int) {
  is >> globalParams >> clusterFissioner >> clusterDecayer
     >> split >> decayer
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

  static Reference<UA5Handler,HwDecayHandler> 
    interfacePartonicHadronizer("Decayer", "Pointer to the object which decays particles.",
                                &UA5Handler::decayer, false, false, true, false);

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
                &Herwig::UA5Handler::M1,0.4*GeV,0.,100.,false,false,false);

  static Parameter<UA5Handler,double>
    interfaceM2("M2", "Parameter M2 used to generate soft cluster mass",
                &Herwig::UA5Handler::M2,2.0*GeV,0.,100.,false,false,false);  

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
   // Constants that should not need changing. This corresponds to ~1fm
   static const double vclx = 4e-12; 
   static const double vcly = 4e-12; 
   static const double vclz = 4e-12; 
   static const double vclt = 4e-12; 

  // Now form the first two clusters, first split gluons
  ClusterVector clusters;
  StepPtr step = ch.newStep();

  tPPtr clu[2];
  Lorentz5Momentum cluP[2];
  tPVector::const_iterator it;

  int i = 0;
  int ntry = 0;
  int netc = 0;
  bool multiplicityNeeded = true;
  bool newCluster = true;
  long id1,id2,id3;
  id1 = id2 = id3 = 0;

  // Lets find the clusters, set the partons inside to be on shell and no momentum
  for(it = tagged.begin(); it!=tagged.end(); it++) {
     if((*it)->id() == ExtraParticleID::Cluster) {
        clu[i] = step->copyParticle(*it);
        cluP[i] = clu[i]->momentum();
        ClusterPtr c = dynamic_ptr_cast<ClusterPtr>(clu[i]);
        c->particle(0)->setMomentum(LorentzMomentum(0.,0.,0.,c->particle(0)->nominalMass()));
        c->particle(1)->setMomentum(LorentzMomentum(0.,0.,0.,c->particle(1)->nominalMass()));
        i++;
     }
  }
  if(i == 0) return;

  // and their cm
  Lorentz5Momentum cm = clu[0]->momentum() + clu[1]->momentum();
  double theCM = cm.mass();
  theCM *= enhanceCM;
  PPtr incomingHadron[2];
 
  // Now initialize charge information
  int theMult = 0;
  PVector tag;
  PPair incom = step->collision()->incoming();
  incomingHadron[0] = incom.first;
  incomingHadron[1] = incom.second;
  netc = 0;
  // Get total charge of collision
  //netc = (int)(incomingHadron[0]->data().charge() + incomingHadron[1]->data().charge());
  // Now subtract the charges that go to hard process
  //PPair incomPart = step->collision()->primarySubProcess()->incoming();
  //netc -= (int)(incomPart.first->data().charge() + incomPart.second->data().charge());

  if(netc == 0) id3 = rndbool() ? ParticleID::u : ParticleID::d;
  else if(netc == -1) id3 = ParticleID::u;
  else if(netc == 1) id3 = ParticleID::d; 
  else id3 = rndbool() ? ParticleID::u : ParticleID::d;

  //netc = 0;
  //id3 = rndbool() ? ParticleID::u : ParticleID::d;
  int numberCluster = 0;
  double sumMasses = 0.;
  int modCharge = 0;
  int numCharge = 0; 
  int newHads;
  pair<tPPtr,tPPtr> hads;
  int nppbar = 0;

  StepPtr newStep = ch.newStep();

  // Loop until we find a match to the charge multiplicity
  while(multiplicityNeeded && ntry < max_tries) {
     PPtr part1, part2;
     if(ntry % 50 == 0) { nppbar = multiplicity(theCM); }
     ntry++;
     numberCluster = 0;
     theMult = nppbar;
     sumMasses = 0.;
     modCharge = 0;

     // If we have restarted after failed attempt, clear step of old particles
     if(ntry > 1) { 
       ParticleSet parts = newStep->all();
       ParticleSet::iterator it;
       for(it = parts.begin(); it!=parts.end(); it++) { 
         if((*it)->birthStep()==newStep) newStep->removeParticle(*it);
       } 
     }
     tag.clear();

     numCharge = netc;
     newCluster = true;

     tag.push_back(clu[0]);
     tag.push_back(clu[1]);

     while(newCluster) {
        ClusterPtr cluster;
        // Choose new constituents
        if(numberCluster < 2) {
           cluster = dynamic_ptr_cast<ClusterPtr>(tag[numberCluster]);
           part1 = cluster->particle(0);
           part2 = cluster->particle(1);
           id1 = part1->id();
           id2 = part2->id();
        } else {
           id1 = -id2;
           if(numberCluster == 2) id1 = id3;
           id2 = rndbool() ? ParticleID::ubar : ParticleID::dbar;
           part1 = getParticle(id1);
           part2 = getParticle(id2); 
           cluster = new_ptr(Cluster(part1,part2));
           newStep->addIntermediate(cluster);
        }
        // Mass = Random number from dN/d(x**2)=exp(-b*x) with mean 1/M2
        double newMass = getParticleData(id1)->constituentMass() + getParticleData(id2)->constituentMass() 
                                 + M1 - log(rnd()*rnd())/M2;
        Lorentz5Momentum cp(0.,0.,0.,newMass,newMass);

        cluster->set5Momentum(cp);
       
        // Now the gaussian distribution of the x,y,z components, and a time component given by
        // sqrt(vx^2+vy^2+vz^2) - vclt*log(r)
        double vx = gaussDistribution(0.,vclx);
        double vy = gaussDistribution(0.,vcly);
        double vz = gaussDistribution(0.,vclz);
        LorentzPoint vert(vx,vy,vz, sqrt(sqr(vx)+sqr(vy)+sqr(vz)) - vclt*log(rnd()));
        cluster->setVertex(vert);
        
        // Now need to measure displacement relative to soft cm (TODO:)
	// Fortran code sets the vertex of the primary incoming particle to 0
        // then during boosting it adds the value of the vertex of particle 3 (???) to all the particles
        // for now I will not implement this. Perhaps Mike can decide if this is needed, as he added this
        // code to the fortran code on 07/03/05.

        // Now we decay the clusters into hadrons
        PPair products = clusterDecayer->decayIntoTwoHadrons(cluster);
        if(numberCluster >= 2) tag.push_back(cluster);
        addFission(products, newStep, cluster, hads, newHads, id1, id2);
        sumMasses += cluster->mass();
      
        // Now we decays hadrons into stable particles. Count charged multiplicity
        for(i = 0; i<newHads; i++) {
           tPPtr particle;
           if(i==0) particle = hads.first;
           else particle = hads.second;
           try { 
              if(!particle->data().stable()) decayer->performDecay(particle,*newStep);
           } catch(exception &e) {}
        }
        numCharge = 0;
        tPVector finalParts = newStep->getFinalState();
        countCharges(finalParts, numCharge, modCharge, newStep);

        // Count charges counts all charges, so only add modCharge for cluster 2
        if(numberCluster == 0) theMult = nppbar+netc;
        else if(numberCluster == 1) theMult += abs(modCharge);
        if(numberCluster == 1 && theMult < 0) theMult += 2;

        numberCluster++;

        // Now check which loop to do next
        if(numCharge < theMult) continue;
        else if(numCharge > theMult) { newCluster = false; }
        else { newCluster = false; multiplicityNeeded = false; }
     }
     // Now check the physical mass/energy boundary
     if(!multiplicityNeeded && sumMasses > theCM) { multiplicityNeeded = true; }
  }

  // Catch case of too many attempts
  if(ntry == max_tries) {
     if(ntry > 1) { 
       ParticleSet parts = newStep->all();
       ParticleSet::iterator it;
       for(it = parts.begin(); it!=parts.end(); it++) { 
         if((*it)->birthStep()==newStep) newStep->removeParticle(*it);
       } 
     }
     // Just hadronize and decay the two clusters
     for(int i =0; i<2; i++) {
        ClusterPtr cluster = dynamic_ptr_cast<ClusterPtr>(clu[i]);
        PPair products = clusterDecayer->decayIntoTwoHadrons(cluster); 
        addFission(products, newStep, cluster, hads, newHads, id1, id2);
        for(int i = 0; i<newHads; i++) {
           tPPtr particle;
           if(i==0) particle = hads.first;
           else particle = hads.second;
           try {
              if(!particle->data().stable()) decayer->performDecay(particle, *newStep);
           } catch(exception &e) {}
        }
     }
  }

  try {
     // Now generate momentum 
     generateMomentum(tag, theCM, cm);
  } catch (Veto &v) { throw(v); }
}

void UA5Handler::addFission(PPair &products, StepPtr &newStep, ClusterPtr &cluster, 
                            pair<tPPtr,tPPtr> &hads, int &newHads, long id1, long id2) {
  if(products.first == PPtr() || products.second == PPtr()) {
    Lorentz5Momentum mom = cluster->momentum();
    LorentzPoint vert = cluster->vertex();
    products = clusterFissioner->produceHadron(id1,id2, mom, vert);
    cluster->addChild(products.first);
    newStep->addDecayProduct(products.first);
    hads.first = products.first;
    newHads = 1;
  } else {
    cluster->addChild(products.first);
    cluster->addChild(products.second);
    newStep->addDecayProduct(products.first);
    newStep->addDecayProduct(products.second);
    hads.first = products.first;
    hads.second = products.second;
    newHads = 2;
  }
}

void UA5Handler::countCharges(tPVector &particles, int &numCharges, int &modCharge, StepPtr &newStep) {
   tPVector::iterator it;
   modCharge = 0;
   // Only count charges of particles produced in this step
   for(it = particles.begin(); it!=particles.end(); it++) {
      if((*it)->data().charged() && (*it)->birthStep() == newStep) {
         numCharges += abs((int)((*it)->data().charge()));
         modCharge += (int)((*it)->data().charge());
      }
   }
}

double UA5Handler::gaussDistribution(double mean, double stdev) {
   double x = rnd();
   x = sqrt(-2.*log(x));
   double y;
   randAzm(x,x,y);
   return mean + stdev*x;
}

// Generates a random azimuthal angle and creates a 2 vector of length x with angle phi
void UA5Handler::randAzm(double x, double &px, double &py) {
   double c,s,cs;
   while(true) {
      c = 2.*rnd()-1.;
      s = 2.*rnd()-1.;
      cs = c*c+s*s;
      if(cs <= 1.) break;
   }
   double qt = x/cs;
   px = (c*c-s*s)*qt;
   py = 2.*c*s*qt;
}

// This returns a random number with a flat distribution [-A,A] plus gaussian tail with stdev B
double UA5Handler::randUng(double A, double B) {
  double prun;
  if(A == 0.) prun = 0.;
  else prun = 1./(1.+B*1.2533/A);
  if(rnd() < prun) return 2.*(rnd()-0.5)*A;
  else {
    double temp = gaussDistribution(0.,B);
    if(temp < 0) return temp - abs(A);
    else return temp + abs(A);
  }
}

// This returns random number from dN/d(x**2)=exp(-B*TM) distribution, where
// TM = SQRT(X**2+AM0**2).  Uses Newton's method to solve F-R=0
double UA5Handler::randExt(double AM0, double B) {
  double r = rnd();
  // Starting value
  double am = AM0-log(r)/B;
  for(int i = 1; i<20; i++) {
    double a = exp(-B*(am-AM0))/(1.+B*AM0);
    double f = (1.+B*am)*a-r;
    double df = -B*B*am*a;
    double dam = -f/df;
    am += dam;
    if(am<AM0) am=AM0+.001;
    if(abs(dam) < .001) break;
  }
  return am;
}

// transforms B (given in rest from of A). Returns vector in lab frame
Lorentz5Momentum UA5Handler::transformToLab(Lorentz5Momentum &A, Lorentz5Momentum &B) {
   if(B.e() == B.mass()) return B;
   double pf4 = (B.px()*A.px() + B.py()*A.py() + B.pz()*A.pz() + B.e()*A.e())/A.mass();
   double fn = (pf4 + B.e())/(A.e() + A.mass());
   Lorentz5Momentum C;
   C.setMass(B.mass());
   C.setPx(B.px() + fn*A.px());
   C.setPy(B.py() + fn*A.py());
   C.setPz(B.pz() + fn*A.pz());
   C.setE(pf4);
   return C;
}
// This returns the mean multiplicity for the energy E amd the given parameters N1,N2,N3
double UA5Handler::meanMultiplicity(Energy E) {
  return N1*pow(E,N2)+N3;
}

// This returns the randomly generated value for the negative binomial
double UA5Handler::negativeBinomial(int N, double mean, double ek) {
  if(N < 0) return 0.0;
  double r = mean/ek;
  double rval = pow(1.+r, -ek);
  r /= (1.+r);
  for(int i = 1; i<=N; i++) rval *=  r*(ek+double(i)-1.)/double(i);
  return rval;
}

// Generates the multiplicity of the charged particles for the energy E
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
  double rn = rand();
  for(i = 0; i<=imax; i++) if(rn < dist[i]) break;
  return 2*(i+1);
}

LorentzRotation UA5Handler::rotate(LorentzMomentum &p) {
  LorentzRotation R;
  static double ptcut = 1e-20;
  //static double wn = 1.;

  double pt = sqr(p.px())+sqr(p.py());
  double pp = sqr(p.pz())+pt;
  double ct, cf;
  double phi, theta;
  if(pt <= pp*ptcut) {
     if(p.pz() > 0) theta = 0.;
     else theta = pi;
     phi = 0.;
  } else {
     pp = sqrt(pp);
     pt = sqrt(pt);
     ct = p.pz()/pp;
     cf = p.px()/pt;
     phi = -acos(cf);
     theta = acos(ct);
  }
  // Rotate first aronu the z axis to put p in the x-z plane
  R = R.rotateZ(phi);
  // Then rotate around the Y axis to put p on the z axis
  R = R.rotateY(theta);
  return R;
}

// Generate the momentum. This is based on the procedure of Jadach from Computer Physics Communications 9 (1975) 297-304
void UA5Handler::generateMomentum(PVector &clusters, double CME, Lorentz5Momentum cm) throw(Veto) {
  // begin with the cylindrical phase space generation described in the paper of Jadach
  generateCylindricalPS(clusters, CME);
  LorentzVector bmp = clusters[0]->momentum();
  bmp = bmp.boost(cm.findBoostToCM());
  // Rotation to put bmp on the z axis
  LorentzRotation R = rotate(bmp);
  // We need to use the inverse
  R = R.inverse();
  for(unsigned int i = 0; i<clusters.size(); i++) {
     // So we now take each cluster and transform it according to the rotation defined above
     Lorentz5Momentum origP = clusters[i]->momentum();
     clusters[i]->transform(R);
     // Then we transform back to the lab frame (away from cm frame)
     Lorentz5Momentum oldP = clusters[i]->momentum();
     Lorentz5Momentum newP;
     try {
        newP = transformToLab(cm, oldP);
        clusters[i]->deepBoost(newP.boostVector());
        clusters[i]->set5Momentum(newP);
     } catch(exception &e) {
        cout << "Caught an problem boosting. " << e.what() << endl;
        cout << "Must decide how to handle this..." << endl;
        cout << "Old p = " << oldP << endl << "New p = " << newP << endl;
        cout << "This is from cluster " << *clusters[i] << " and has z component > energy" << endl;
        cout << "Cluster 0 is = " << *clusters[0] << endl;
        cout << "Original momentum = " << origP << endl;
        cout << "The cm vector is = " << cm << endl;
        cout << "Mass error of original momentum is " << origP.massError() << endl;
        throw Veto();
     }
  }
}

void UA5Handler::generateCylindricalPS(PVector &clusters, double CME) {
  double alog = log(CME*CME);
  int ncl = clusters.size();
  double eps = 1e-10/double(ncl);
  double pt, pt2, px, py, sumx, sumy, sum1, sum2, sum3, sum4, sumpt, sumtm;
  double ximin, ximax;
  int its = 0;
  vector<double> xi(ncl);
  while(its < max_tries) {
    its++;
    sumx = sumy = 0.;
    for(int i = 0; i<ncl; i++) {
      // Decide between three options
      // If the q-qbar pair used to create the hadrons from the cluster is u or d, option1
      // If the q-qbar pair is a c or s, option 2
      // if the q-qbar pair is a diquark, option 3
      // if their was no q-qbar pair (for cluster->hadron) then option 1
      double slop;
      tPPtr hadron = clusters[i]->children()[0];
      long id1 = dynamic_ptr_cast<ClusterPtr>(clusters[i])->particle(0)->id();
      
      long ids[3];
      long hadId = hadron->id();
      ids[2] = (abs(hadId)/1000)%10;;
      ids[1] = (abs(hadId)/100)%10;
      ids[0] = (abs(hadId)/10)%10;
      if(ids[2] != 0 && hadId < 0) { ids[2] = -ids[2]; ids[1] = -ids[1]; ids[0] = -ids[0]; }
      else if(ids[2] == 0 && hadId < 0) ids[1] = -ids[1]; 
      else if(ids[2]) ids[0] = -ids[0];
      
      bool isAMeson = (ids[2] == 0);
      bool isSingle = (clusters[i]->children().size() == 1);
      if(isAMeson) {
         long newId = 0;
         if(ids[1] == id1) newId = abs(ids[0]);
         else if(ids[0] == id1) newId = abs(ids[1]);
         if(newId == ParticleID::s || id1 == ParticleID::c) slop = P2;
         else slop = P1; 
      } else if(isSingle) slop = P1;
      else slop = P3;

      // Now generate the pt given the parameter slop
      pt = randExt(clusters[i]->mass(), slop);
      Lorentz5Momentum mom = clusters[i]->momentum();
      double ptp = pt*pt - sqr(mom.mass());
      if(ptp <= 0) pt = - sqrt(-ptp);
      else pt = sqrt(ptp);
      randAzm(pt,px,py);
      mom.setPx(px);
      mom.setPy(py);
      sumx += px;
      sumy += py;
      clusters[i]->set5Momentum(mom);
    }
    sumx /= ncl;
    sumy /= ncl;
    sumpt = sumtm = 0.;
    for(int i = 0; i<ncl; i++) {
      Lorentz5Momentum mom = clusters[i]->momentum();
      mom.setPz(0.);
      mom.setX(mom.px()-sumx);
      mom.setY(mom.py()-sumy);
      pt2 = sqr(mom.px()) + sqr(mom.py());
      // Use the z component of the clusters momentum for temporary storage
      mom.setPz(sqrt(pt2+sqr(mom.mass())));
      sumtm += mom.pz();
      clusters[i]->set5Momentum(mom);
    }
    if(sumtm > CME) continue; 
    for(int i = 0; i<ncl; i++) xi[i] = randUng(0.6,1.0);
    sort(xi.begin(), xi.end());
    ximin = xi[0];
    ximax = xi[xi.size()-1]-ximin;
    xi[0] = 0.;
    for(int i = ncl-2; i>1; i--) xi[i+1] = (xi[i]-ximin)/ximax;
    xi[1] = 1.;
    double yy = log(CME*CME/(clusters[0]->momentum().pz()*clusters[1]->momentum().pz()));
    int j;
    for(j = 0; j<10; j++) {
      sum1 = sum2 = sum3 = sum4 = 0.;
      for(int i = 0; i<ncl; i++) {
        double tm = clusters[i]->momentum().pz();
        double ex = exp(yy*xi[i]);
        sum1 += tm*ex; sum2 += tm/ex;
        sum3 += (tm*ex)*xi[i]; sum4 += (tm/ex)*xi[i];
      }
      double fy = alog-log(sum1*sum2);
      double dd = (sum3*sum2 - sum1*sum4)/(sum1*sum2);
      double dyy = fy/dd;
      if(abs(dyy/yy) < eps) { yy += dyy; break; }
      yy += dyy;
    }
    if(j == 10) {
      if(its > 100) eps *= 10.; 
      continue;
    }
    double zz = log(CME/sum1);
    for(int i = 0; i<ncl; i++) {
      Lorentz5Momentum mom = clusters[i]->momentum();
      double tm = mom.pz();
      double E1 = exp(zz + yy*xi[i]);
      mom.setPz(tm/2.*(1./E1-E1));
      mom.setE(tm/2.*(1./E1+E1));
      clusters[i]->set5Momentum(mom);
    }
  } 
}

// This is all just administrative functions for ThePEG structure
IBPtr UA5Handler::clone() const { return new_ptr(*this); }
  
IBPtr UA5Handler::fullclone() const { return clone(); }

void UA5Handler::doupdate() throw(UpdateException) {
  HadronizationHandler::doupdate();
}

void UA5Handler::doinit() throw(InitException) {
  HadronizationHandler::doinit();
}

void UA5Handler::dofinish() {
  HadronizationHandler::dofinish();
}

void UA5Handler::rebind(const TranslationMap &trans) throw(RebindException) {
  HadronizationHandler::rebind(trans);
}

IVector UA5Handler::getReferences() {
  IVector ref = HadronizationHandler::getReferences();
  return ref;
}


 
