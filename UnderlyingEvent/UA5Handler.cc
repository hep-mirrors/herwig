#include "UA5Handler.h"
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Repository/UseRandom.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Handlers/DecayHandler.h>
#include "Herwig++/Hadronization/Cluster.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Timer.h"
#include <cassert>

using namespace std;
using namespace ThePEG;
using namespace Herwig;

// Default constructor
UA5Handler::UA5Handler() 
  : _n1(9.11), _n2(0.115), _n3(-9.5), _k1(0.029), _k2(-0.104),
    _m1(0.4*GeV), _m2(2./GeV), _p1(5.2/GeV), _p2(3.0/GeV), _p3(5.2/GeV),
    _probSoft(1.0), _enhanceCM(1.), _maxtries(300)
{}  

// Copy constructor
UA5Handler::UA5Handler(const UA5Handler &h) :
  HadronizationHandler(h),
  clusterFissioner(h.clusterFissioner),
  clusterDecayer(h.clusterDecayer),
  _n1(h._n1), _n2(h._n2), _n3(h._n3), _k1(h._k1), _k2(h._k2),
  _m1(h._m1), _m2(h._m2), _p1(h._p1), _p2(h._p2), _p3(h._p3),
  _probSoft(h._probSoft), _enhanceCM(h._enhanceCM),
  _maxtries(h._maxtries)
{}

// Saving things into run file
void UA5Handler::persistentOutput(PersistentOStream &os) const {
  os << clusterFissioner << clusterDecayer
     << _n1 << _n2 << _n3 << _k1 << _k2 << _m1 << _m2 << _p1
     << _p2 << _p3 << _probSoft << _enhanceCM << _maxtries;
}

// Reading them back in, in the same order
void UA5Handler::persistentInput(PersistentIStream &is, int) {
  is >> clusterFissioner >> clusterDecayer
     >> _n1 >> _n2 >> _n3 >> _k1 >> _k2 >> _m1 >> _m2 >> _p1
     >> _p2 >> _p3 >> _probSoft >> _enhanceCM >> _maxtries;
}

// We must define this static member for ThePEG
ClassDescription<UA5Handler> UA5Handler::initUA5Handler;

void UA5Handler::Init() {

  static ClassDocumentation<UA5Handler> documentation
    ("This is the simple UA5 model for the underlying event.",
     "The underlying event was simulated using the model of "
     "the UA5 collaboration, \\cite{Alner:1986is}.",
     "\\bibitem{Alner:1986is} G.~J.~Alner {\\it et al.}  "
     "UA5 Collaboration, Nucl.\\ Phys.\\ B{\bf 291} (1987) 445.");

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

  static Parameter<UA5Handler,double> interfaceN1
    ("N1",
     "The parameter n_1 in the mean charged multiplicity",
     &UA5Handler::_n1, 9.110, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<UA5Handler,double> interfaceN2
    ("N2",
     "The parameter n_2 in the mean charged multiplicity",
     &UA5Handler::_n2, 0.115, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<UA5Handler,double> interfaceN3
    ("N3",
     "The parameter n_3 in the mean charged multiplicity",
     &UA5Handler::_n3, -9.500, -100.0, 100.0,
     false, false, Interface::limited);

  static Parameter<UA5Handler,double> interfaceK1
    ("K1",
     "The parameter k_1 used to generate the multiplicity distribution",
     &UA5Handler::_k1, 0.029, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<UA5Handler,double> interfaceK2
    ("K2",
     "The parameter k_2 used to generate the multiplicity distribution",
     &UA5Handler::_k2, -0.104, -1.0, 1.0,
     false, false, Interface::limited);

  static Parameter<UA5Handler,Energy> interfaceM1
    ("M1",
     "The parameter m_1 used to generate the cluster mass distribution.",
     &UA5Handler::_m1, GeV, 0.4*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<UA5Handler,InvEnergy> interfaceM2
    ("M2",
     "The parameter m_2 used to generate the cluster mass distribution.",
     &UA5Handler::_m2, 1./GeV, 2.0/GeV, 0.1/GeV, 10.0/GeV,
     false, false, Interface::limited);

  static Parameter<UA5Handler,InvEnergy> interfaceP1
    ("P1",
     "Slope used to generate the pt of the u,d soft clusters",
     &UA5Handler::_p1, 1./GeV, 5.2/GeV, 0.1/GeV, 10.0/GeV,
     false, false, Interface::limited);

  static Parameter<UA5Handler,InvEnergy> interfaceP2
    ("P2",
     "Slope used to generate the pt of the s,c soft clusters",
     &UA5Handler::_p2, 1./GeV, 3.0/GeV, 0.1/GeV, 10.0/GeV,
     false, false, Interface::limited);

  static Parameter<UA5Handler,InvEnergy> interfaceP3
    ("P3",
     "Slope used to generate the pt of the qq soft clusters",
     &UA5Handler::_p3, 1./GeV, 5.2/GeV, 0.1/GeV, 10.0/GeV,
     false, false, Interface::limited);

  static Parameter<UA5Handler,double> interfaceProbSoft
    ("ProbSoft",
     "The probability of having a soft underlying event.",
     &UA5Handler::_probSoft, 1.0, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<UA5Handler,double> interfaceEnhanceCM
    ("EnhanceCM",
     "Enhancement of the CM energy in the mult distribution",
     &UA5Handler::_enhanceCM, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<UA5Handler,unsigned int> interfaceMaximumTries
    ("MaximumTries",
     "The maximum number of attempts to generate the clusters",
     &UA5Handler::_maxtries, 300, 100, 1000,
     false, false, Interface::limited);

}

void UA5Handler::insertParticle(PPtr particle,StepPtr step,bool all)
{ 
  if(all) step->addDecayProduct(particle);
  for(unsigned int ix=0;ix<particle->children().size();++ix)
    {insertParticle(particle->children()[ix],step,true);}
}

// This is all just administrative functions for ThePEG structure
IBPtr UA5Handler::clone() const { return new_ptr(*this); }
  
IBPtr UA5Handler::fullclone() const { return clone(); }

// Generates the multiplicity of the charged particles for the energy E
unsigned int UA5Handler::multiplicity(Energy E) {
  double alogs = 2.*log(E/GeV);
  double rk = _k1*alogs+_k2;
  if(rk > 1000.) rk = 1000.;
  double ek = 1./rk;
  double mean = meanMultiplicity(E);
  if(mean < 1.) mean = 1.;
  vector<double> dist;
  double sum = 0.0;
  unsigned int imax = 0;
  unsigned int i;
  for(i = 0; i<500; i++) {
    int N = (i+1)*2;
    dist.push_back(negativeBinomial(N,mean,ek));
    if(dist[i] < 1e-7*sum) break;
    imax = i;
    sum += dist[i];
    dist[i] = sum;
  }
  if (imax==0) dist[0]=1.;
  else if (imax==499) 
    throw Exception() << "Multiplicity too large in UA5Handler::multiplicity()" 
		      << Exception::eventerror;
  else {for(i = 0; i<=imax; i++) dist[i] /= sum;}
  double rn = rnd();
  for(i = 0; i<=imax; i++) if(rn < dist[i]) break;
  return 2*(i+1);
}

LorentzRotation UA5Handler::rotate(LorentzMomentum &p) {
  LorentzRotation R;
  static double ptcut = 1e-20;
  Energy pt = sqr(p.px())+sqr(p.py());
  Energy pp = sqr(p.pz())+pt;
  double phi, theta;
  if(pt <= pp*ptcut) {
     if(p.pz() > 0) theta = 0.;
     else theta = pi;
     phi = 0.;
  } else {
     pp = sqrt(pp);
     pt = sqrt(pt);
     double ct = p.pz()/pp;
     double cf = p.px()/pt;
     phi = -acos(cf);
     theta = acos(ct);
  }
  // Rotate first aronu the z axis to put p in the x-z plane
  R = R.rotateZ(phi);
  // Then rotate around the Y axis to put p on the z axis
  R = R.rotateY(theta);
  return R;
}

void UA5Handler::performDecay(PPtr parent,int & totalcharge,int & numbercharge)
{
  // for an already decayed particle apply to children
  if(!parent->children().empty())
    {
      for(unsigned int ix=0;ix<parent->children().size();++ix)
	performDecay(parent->children()[ix],totalcharge,numbercharge);
    }
  // for a stable particle just add the charge
  else if(parent->data().stable())
    {
      int charge=int(parent->data().charge());
      totalcharge  +=    charge ;
      numbercharge +=abs(charge);
    }
  // otherwise decay the particle
  else
    {
      unsigned int ntry = 0,maxtry=100;
      while (1)
	{
	  // veto if too many tries
	  if(++ntry>=maxtry) 
	    throw Exception() << "Too many attempts to decay " << parent->PDGName() 
			      << "in UA5Handler::performDecay(), probably needs a "
			      << "partonic decay of a heavy hadron which isn't"
			      << " implemented yet " << Exception::eventerror;
	  // select the decay mode
	  tDMPtr dm(parent->data().selectMode(*parent));
	  // throw event away if no decay mode
	  if ( !dm ) throw Exception() << "No decay mode for " << parent->PDGName()
				       << "in UA5Handler::performDecay()" 
				       << Exception::eventerror; 
	  // throw event away if no decayer
	  if( !dm->decayer() ) 
	    throw Exception() << "No decayer mode for " << parent->PDGName()
			      << "in UA5Handler::performDecay()" 
			      << Exception::eventerror;
	  try 
	    {
	      ParticleVector children = dm->decayer()->decay(*dm, *parent);
	      // decay generates decay products
	      if ( !children.empty() ) 
		{
		  // special for partonic decays
		  // see if colour particles produced
		  for(unsigned int ix=0;ix<children.size();++ix)
		    {if(children[ix]->coloured()){throw Veto();}}
		  // set up parent
		  parent->decayMode(dm);
		  // add children
		  for ( int i = 0, N = children.size(); i < N; ++i )
		    {
		      children[i]->setLabVertex(parent->labDecayVertex());
		      parent->addChild(children[i]);
		    }
		  parent->scale(0.0*GeV2);
		  // loop over the children and decay
		  for ( int i = 0, N = children.size(); i < N; ++i )
		    {performDecay(children[i],totalcharge,numbercharge);}
		  return;
		}
	    }
	  catch (Veto) {}
	}
    }
}

void UA5Handler::decayCluster(ClusterPtr cluster,bool single)
{
  PPair products = clusterDecayer->decayIntoTwoHadrons(cluster);
  if(products.first == PPtr() || products.second == PPtr()) 
    {
      if(!single) 
	throw Exception() << "Can't decay cluster in UA5Handler::decayCluster()"
			  << Exception::eventerror;
      // decay the cluster to one hadron
      Lorentz5Momentum mom = cluster->momentum();
      LorentzPoint vert = cluster->vertex();
      int id1=cluster->particle(0)->id();
      int id2=cluster->particle(1)->id();
      products=clusterFissioner->produceHadron(id1,id2, mom, vert);
      // put the cluster and the hadron on mass-shell
      Energy mass=products.first->nominalMass();
      Lorentz5Momentum newp(0.,0.,0.,mass,mass);
      cluster->set5Momentum(newp);
      products.first->set5Momentum(newp);
      cluster->addChild(products.first);
    }
  else
    {
      cluster->addChild(products.first);
      cluster->addChild(products.second);
    }
}

// This is the routine that is called to start the algorithm. 
void UA5Handler::handle(EventHandler &ch, const tPVector &tagged,
			const Hint &hint) throw(Veto,Stop,Exception) {
  useMe();
  Timer<10000> timer("UA5Handler::handle()");
  // create a new step for the products
  StepPtr newStep = ch.newStep();
  // Constants that should not need changing. This corresponds to ~1fm
  static const Length vclx = 4e-12; 
  static const Length vcly = 4e-12; 
  static const Length vclz = 4e-12; 
  static const Length vclt = 4e-12;
  // Find the first two clusters
  // Lets find the clusters, set the partons inside to be on shell and no momentum
  tClusterPtr clu[2];
  Lorentz5Momentum cluP[2];
  tPVector::const_iterator it;
  unsigned int i = 0;
  for(it = tagged.begin(); it!=tagged.end(); ++it) {
    if((*it)->id() != ExtraParticleID::Cluster) continue;
    clu[i] = dynamic_ptr_cast<ClusterPtr>(*it);
    cluP[i] = clu[i]->momentum();
    ++i;
    if(i>2) throw Exception() << "Must have at most two beam clusters in "
			      << "UA5Handler::handle " 
			      << Exception::eventerror;
  }
  if(i == 0) return;
  if(i!=2) throw Exception() << "Must have two and only two beam clusters in "
			     << "UA5Handler::handle " 
			     << Exception::eventerror;
  // if not generating the soft underlying event 
  // just hadronize and decay the two clusters
  if(rnd()>_probSoft)
    {
      for(int i =0; i<2; i++) {
	ClusterPtr cluster = clu[i];
	decayCluster(cluster,false);
	int totalcharge(0),numbercharge(0);
	for(unsigned int ix=0;ix<cluster->children().size();++ix)
	  {performDecay(cluster->children()[ix],totalcharge,numbercharge);}
	insertParticle(cluster,newStep,false);
      }
      // don't to the rest
      return;
    }
  // and their cm
  Lorentz5Momentum cm = clu[0]->momentum() + clu[1]->momentum();
  Energy theCM = cm.mass();
  // soft mass with optional enhancement for multiplicity only
  Energy softCM = _enhanceCM*theCM;
  int netc = 0;
  long id1(0),id2(0),id3= rndbool() ? ParticleID::u : ParticleID::d;
  // storage for the multiplicity
  int nppbar = 0;
  // Loop until we find a match to the charge multiplicity
  unsigned int ntry = 0;
  vector<ClusterPtr> clusters;
  bool multiplicityReached = false;
  while(!multiplicityReached && ntry < _maxtries) {
     PPtr part1, part2;
     // reset multiplicity every 50 tries
     if(ntry % 50 == 0) nppbar = multiplicity(softCM);
     ++ntry;
     unsigned int numberCluster = 0;
     int theMult = nppbar;
     Energy sumMasses = 0.;
     int modCharge = 0;
     // delete the particles from the previous attempt if needed
     if(ntry > 1) clusters.clear();
     int numCharge = netc;
     bool newCluster = true;
     while(newCluster) {
       ClusterPtr cluster;
       // Choose new constituents
       if(numberCluster < 2) {
	 part1 = clu[numberCluster]->particle(0);
	 part2 = clu[numberCluster]->particle(1);
	 id1 = part1->id();
	 id2 = part2->id();
	 cluster = new_ptr(Cluster(part1,part2));
       } else {
	 id1 = -id2;
	 if(numberCluster == 2) id1 = id3;
	 id2 = rndbool() ? ParticleID::ubar : ParticleID::dbar;
	 part1 = getParticle(id1);
	 part2 = getParticle(id2); 
	 cluster = new_ptr(Cluster(part1,part2));
       }
       // Mass = Random number from dN/d(x**2)=exp(-b*x) with mean 1/M2
       Energy newMass = 
	 getParticleData(id1)->constituentMass() + 
	 getParticleData(id2)->constituentMass() 
	 + _m1 - log(rnd()*rnd())/_m2;
       // set momentum of the cluster
       Lorentz5Momentum cp(0.,0.,0.,newMass,newMass);
       cluster->set5Momentum(cp);
       // Now the gaussian distribution of the x,y,z components,
       // and a time component given by
       // sqrt(vx^2+vy^2+vz^2) - vclt*log(r)
       Length vx = gaussDistribution(0.,vclx);
       Length vy = gaussDistribution(0.,vcly);
       Length vz = gaussDistribution(0.,vclz);
       LorentzPoint vert(vx,vy,vz, sqrt(sqr(vx)+sqr(vy)+sqr(vz)) - vclt*log(rnd()));
       cluster->setVertex(vert);
       // Now need to measure displacement relative to soft cm (TODO:)
       // Fortran code sets the vertex of the primary incoming particle to 0
       // then during boosting it adds the value of the vertex of particle 3 (???)
       // to all the particles
       // for now I will not implement this. Perhaps Mike can decide if this is 
       // needed, as he added this code to the fortran code on 07/03/05.

       // Add the cluster to the list
       clusters.push_back(cluster);
       // Now we decay the clusters into hadrons
       decayCluster(cluster,true);
       // sum of the cluster masses
       sumMasses += cluster->mass();
       // Now decay the hadrons into stable particles
       int totalcharge(0),numbercharge(0);
       for(unsigned int ix=0;ix<cluster->children().size();++ix)
	 {performDecay(cluster->children()[ix],totalcharge,numbercharge);}
       numCharge+=numbercharge;
       modCharge+=totalcharge;
       // Count charges counts all charges, so only add modCharge for cluster 2
       if(numberCluster == 0) theMult = nppbar+netc;
       else if(numberCluster == 1) theMult += abs(modCharge);
       if(numberCluster == 1 && theMult < 0) theMult += 2;
       numberCluster++;
       // Now check which loop to do next
       if(numCharge < theMult) continue;
       else if(numCharge > theMult) { newCluster = false; }
       else { newCluster = false; multiplicityReached = true; }
     }
     // Now check the physical mass/energy boundary
     if(multiplicityReached && (sumMasses > theCM || clusters.size()<2)) 
       { multiplicityReached = false; }
  }
  // Catch case of too many attempts
  if(ntry == _maxtries) {
    // Just hadronize and decay the two clusters
    for(int i =0; i<2; i++) {
      ClusterPtr cluster = clu[i];
      decayCluster(cluster,false);
      int totalcharge(0),numbercharge(0);
      for(unsigned int ix=0;ix<cluster->children().size();++ix)
	{performDecay(cluster->children()[ix],totalcharge,numbercharge);}
      insertParticle(cluster,newStep,false);
    }
    // don't to the rest
    return;
  }
  // Now generate momentum of the clusters
  try {generateMomentum(clusters, theCM, cm);} 
  catch (Veto &v) {throw(v);}
  // insert the particles into the event record
  for(unsigned int ix=0;ix<clusters.size();++ix)
    {
      clu[0]->addChild(clusters[ix]);
      clu[1]->addChild(clusters[ix]);
      newStep->addDecayProduct(clusters[ix]);
    }
  for(unsigned int ix=0;ix<clusters.size();++ix)
    {insertParticle(clusters[ix],newStep,false);}
}

// Generate the momentum. This is based on the procedure of
// Jadach from Computer Physics Communications 9 (1975) 297-304
void UA5Handler::generateMomentum(ClusterVector &clusters, double CME, Lorentz5Momentum cm) throw(Veto) {
  Timer<10001> timer("UA5Handler::generateMomentum()");
  // begin with the cylindrical phase space generation described in the paper of Jadach
  generateCylindricalPS(clusters, CME);
  LorentzVector bmp = clusters[0]->momentum();
  bmp = bmp.boost(cm.findBoostToCM());
  // Rotation to put bmp on the z axis
  LorentzRotation R = rotate(bmp);
  // We need to use the inverse
  R = R.inverse();
  Hep3Vector boostv=cm.boostVector();
  for(unsigned int i = 0; i<clusters.size(); i++) {
    // So we now take each cluster and transform it according to the
    // rotation defined above
    Lorentz5Momentum origP = clusters[i]->momentum();
    clusters[i]->transform(R);
    // Then we transform back to the lab frame (away from cm frame)
    Lorentz5Momentum oldP = clusters[i]->momentum();
    Lorentz5Momentum newP=oldP;
    try {
      newP.boost(boostv);
      clusters[i]->deepBoost(newP.boostVector());
      clusters[i]->set5Momentum(newP);
    } catch(exception &e) {
      cerr << "Caught an problem boosting. " << e.what() << endl;
      cerr << "Must decide how to handle this..." << endl;
      cerr << "Old p = " << oldP << endl << "New p = " << newP << endl;
      cerr << "This is from cluster " << *clusters[i] << " and has z component > energy" << endl;
      cerr << "Cluster 0 is = " << *clusters[0] << endl;
      cerr << "Original momentum = " << origP << endl;
      cerr << "The cm vector is = " << cm << endl;
      cerr << "Mass error of original momentum is " << origP.massError() << endl;
      throw Veto();
    }
  }
}

void UA5Handler::generateCylindricalPS(ClusterVector &clusters, double CME) {
  assert(clusters.size()>=2);
  double alog = log(CME*CME);
  unsigned int ncl = clusters.size();
  double eps = 1e-10/double(ncl);
  Energy pt, pt2, px, py, sumpt, sumtm;
  double sumx, sumy, sum1=0., sum2, sum3, sum4;
  double ximin, ximax;
  unsigned int its = 0;
  vector<double> xi(ncl);
  // calculate the slope parameters for the different clusters
  // outside loop to save time
  vector<InvEnergy> slope(ncl);
  vector<Lorentz5Momentum> mom(ncl);
  for(unsigned int ix=0;ix<ncl;++ix)
    {
      // Decide between three options
      // If the q-qbar pair used to create the hadrons from the cluster is u or d, option1
      // If the q-qbar pair is a c or s, option 2
      // if the q-qbar pair is a diquark, option 3
      // if their was no q-qbar pair (for cluster->hadron) then option 1
      long id1=clusters[ix]->particle(0)->id();
      long hadId = clusters[ix]->children()[0]->id();
      long ids[3]={(abs(hadId)/10)%10,(abs(hadId)/100)%10,(abs(hadId)/1000)%10};
      if(ids[2] != 0 && hadId < 0) { ids[2] = -ids[2]; ids[1] = -ids[1]; ids[0] = -ids[0]; }
      else if(ids[2] == 0 && hadId < 0) ids[1] = -ids[1]; 
      else if(ids[2]) ids[0] = -ids[0];
      if(ids[2] == 0) {
         long newId = 0;
         if(ids[1] == id1) newId = abs(ids[0]);
         else if(ids[0] == id1) newId = abs(ids[1]);
         if(newId == ParticleID::s || id1 == ParticleID::c) slope[ix] = _p2;
         else slope[ix] = _p1; 
      } else if(clusters[ix]->children().size() == 1) slope[ix] = _p1;
      else slope[ix] = _p3;
      mom[ix]=clusters[ix]->momentum();
    }
  // loop to generate momenta
  double yy(0.);
  while(its < _maxtries) {
    ++its;
    sumx = sumy = 0.;
    for(unsigned int i = 0; i<ncl; ++i) {
      // Now generate the pt given the parameter slope
      pt = randExt(clusters[i]->mass(), slope[i]);
      Energy2 ptp = pt*pt - sqr(mom[i].mass());
      if(ptp <= 0) pt = - sqrt(-ptp);
      else pt = sqrt(ptp);
      randAzm(pt,px,py);
      mom[i].setPx(px);
      mom[i].setPy(py);
      sumx += px;
      sumy += py;
    }
    sumx /= ncl;
    sumy /= ncl;
    sumpt = sumtm = 0.;
    for(unsigned int ix = 0; ix<ncl; ++ix) {
      mom[ix].setPz(0.);
      mom[ix].setX(mom[ix].px()-sumx);
      mom[ix].setY(mom[ix].py()-sumy);
      pt2 = mom[ix].perp2();
      // Use the z component of the clusters momentum for temporary storage
      mom[ix].setPz(sqrt(pt2+mom[ix].mass2()));
      sumtm += mom[ix].pz();
    }
    if(sumtm > CME) continue;
    for(unsigned int i = 0; i<ncl; i++) xi[i] = randUng(0.6,1.0);
    sort(xi.begin(), xi.end());
    ximin = xi[0];
    ximax = xi[xi.size()-1]-ximin;
    xi[0] = 0.;
    for(int i = ncl-2; i>1; i--) xi[i+1] = (xi[i]-ximin)/ximax;
    xi[1] = 1.;
    yy = log(CME*CME/(mom[0].pz()*mom[1].pz()));
    bool suceeded=false;
    for(unsigned int j = 0; j<10; j++) {
      sum1 = sum2 = sum3 = sum4 = 0.;
      for(unsigned int i = 0; i<ncl; i++) {
        double tm = mom[i].pz();
        double ex = exp(yy*xi[i]);
        sum1 += tm*ex;
	sum2 += tm/ex;
        sum3 += (tm*ex)*xi[i];
	sum4 += (tm/ex)*xi[i];      }
      double fy = alog-log(sum1*sum2);
      double dd = (sum3*sum2 - sum1*sum4)/(sum1*sum2);
      double dyy = fy/dd;
      if(abs(dyy/yy) < eps) { yy += dyy; suceeded=true; break;}
      yy += dyy;
    }
    if(suceeded) break;
    if(its > 100) eps *= 10.;
  }
  double zz = log(CME/sum1);
  for(unsigned int i = 0; i<ncl; i++) {
    double tm = mom[i].pz();
    double E1 = exp(zz + yy*xi[i]);
    mom[i].setPz(0.5*tm*(1./E1-E1));
    mom[i].setE( 0.5*tm*(1./E1+E1));
    clusters[i]->set5Momentum(mom[i]);
  }
}
