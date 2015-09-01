// -*- C++ -*-
//
// UA5Handler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#include <ThePEG/Repository/UseRandom.h>
#include "UA5Handler.h"
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Handlers/DecayHandler.h>
#include <ThePEG/Handlers/EventHandler.h>
#include "Herwig/Hadronization/Cluster.h"
#include "Herwig/Hadronization/ClusterFissioner.h"
#include "Herwig/Hadronization/ClusterDecayer.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Throw.h"
#include <cassert>

using namespace std;
using namespace ThePEG;
using namespace Herwig;

// Default constructor
UA5Handler::UA5Handler() 
  : _n1(9.11), _n2(0.115), _n3(-9.5), _k1(0.029), _k2(-0.104),
    _m1(0.4*GeV), _m2(2./GeV), _p1(5.2/GeV), _p2(3.0/GeV), _p3(5.2/GeV),
    _probSoft(1.0), _enhanceCM(1.), _maxtries(300), _needWarning(true)
{}  

// Saving things into run file
void UA5Handler::persistentOutput(PersistentOStream &os) const {
  os << clusterFissioner << clusterDecayer
     << _n1 << _n2 << _n3 << _k1 << _k2 
     << ounit(_m1,GeV) << ounit(_m2,InvGeV)
     << ounit(_p1,InvGeV) << ounit(_p2,InvGeV) << ounit(_p3,InvGeV) 
     << _probSoft << _enhanceCM << _maxtries << _needWarning;
}

// Reading them back in, in the same order
void UA5Handler::persistentInput(PersistentIStream &is, int) {
  is >> clusterFissioner >> clusterDecayer
     >> _n1 >> _n2 >> _n3 >> _k1 >> _k2 
     >> iunit(_m1,GeV) >> iunit(_m2,InvGeV) 
     >> iunit(_p1,InvGeV) >> iunit(_p2,InvGeV) >> iunit(_p3,InvGeV) 
     >> _probSoft >> _enhanceCM >> _maxtries >> _needWarning;
}

// We must define this static member for ThePEG
ClassDescription<UA5Handler> UA5Handler::initUA5Handler;

void UA5Handler::Init() {

  static ClassDocumentation<UA5Handler> documentation
    ("This is the simple UA5 model for the underlying event.",
     "The underlying event was simulated using the model of "
     "the UA5 collaboration, \\cite{Alner:1986is}.",
     "%\\cite{Alner:1986is}\n"
     "\\bibitem{Alner:1986is}\n"
     "  G.~J.~Alner {\\it et al.}  [UA5 Collaboration],\n"
     "  ``The UA5 High-Energy anti-p p Simulation Program,''\n"
     "  Nucl.\\ Phys.\\  B {\\bf 291}, 445 (1987).\n"
     "  %%CITATION = NUPHA,B291,445;%%\n"
     );

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
     &UA5Handler::_m1, GeV, 0.4*GeV, ZERO, 10.0*GeV,
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

  static Switch<UA5Handler,bool> interfaceWarning
    ("Warning",
     "Whether to issue a warning if UA5 and MPI are on at the same time.",
     &UA5Handler::_needWarning, true, false, false);
  static SwitchOption interfaceWarningYes
    (interfaceWarning,
     "Yes",
     "Warn if UA5 and MPI are on at the same time.",
     true);
  static SwitchOption interfaceWarningNo
    (interfaceWarning,
     "No",
     "Print no warnings.",
     false);

}
void UA5Handler::insertParticle(PPtr particle,StepPtr step,bool all) const
{ 
  if(all) step->addDecayProduct(particle);
  for(unsigned int ix=0; ix < particle->children().size(); ++ix) {
    insertParticle(particle->children()[ix],step,true);
  }
}

// This is all just administrative functions for ThePEG structure
IBPtr UA5Handler::clone() const { return new_ptr(*this); }
  
IBPtr UA5Handler::fullclone() const { return new_ptr(*this); }

// Generates the multiplicity of the charged particles for the energy E
unsigned int UA5Handler::multiplicity(Energy E) const {
  double alogs = 2.*log(E/GeV);
  double rk = _k1*alogs+_k2;
  if(rk > 1000.) rk = 1000.;
  double ek = 1./rk;
  double mean = meanMultiplicity(E);
  if(mean < 1.) mean = 1.;

  vector<double> dist;
  dist.reserve(500);
  double sum = 0.0;
  for(int i = 0; i<500; ++i) {
    int N = (i+1)*2;
    double negBin = negativeBinomial(N,mean,ek);
    if(negBin < 1e-7*sum) break;
    sum += negBin;
    dist.push_back(sum);
  }
  unsigned int imax = dist.size();
  if (imax==1) 
    dist[0]=1.;
  else if (imax==500) 
    throw Exception() << "Multiplicity too large in UA5Handler::multiplicity()" 
		      << Exception::eventerror;
  else {
    for(unsigned int i = 0; i<imax; ++i) 
      dist[i] /= sum;
  }
  double rn = rnd();
  unsigned int i=0;
  for(; i<imax; ++i) 
    if(rn < dist[i]) break;
  return 2*(i+1);
}

LorentzRotation UA5Handler::rotate(const LorentzMomentum &p) const {
  LorentzRotation R;
  static const double ptcut = 1e-20;
  Energy2 pt2 = sqr(p.x())+sqr(p.y());
  Energy2 pp2 = sqr(p.z())+pt2;
  double phi, theta;
  if(pt2 <= pp2*ptcut) {
     if(p.z() > ZERO) theta = 0.;
     else theta = Constants::pi;
     phi = 0.;
  } else {
     Energy pp = sqrt(pp2);
     Energy pt = sqrt(pt2);
     double ct = p.z()/pp;
     double cf = p.x()/pt;
     phi = -acos(cf);
     theta = acos(ct);
  }
  // Rotate first around the z axis to put p in the x-z plane
  // Then rotate around the Y axis to put p on the z axis
  R.rotateZ(phi).rotateY(theta);
  return R;
}

void UA5Handler::performDecay(PPtr parent,int & totalcharge,int & numbercharge) const
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
      int charge = parent->data().iCharge()/3;
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
		  parent->scale(ZERO);
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

void UA5Handler::decayCluster(ClusterPtr cluster,bool single) const
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
      tcPDPtr ptrQ = cluster->particle(0)->dataPtr();
      tPPtr newPtr = cluster->particle(1);
      products=clusterFissioner->produceHadron(ptrQ,newPtr,
					       mom, vert);
      // put the cluster and the hadron on mass-shell
      Energy mass=products.first->nominalMass();
      Lorentz5Momentum newp(ZERO,ZERO,ZERO,mass,mass);
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
			const Hint &) {
  // Warn if the event has multiple scatters. 
  // If so, UA5 often has been left on by accident.
  if (_needWarning 
      && ch.currentEvent()->primaryCollision()->subProcesses().size() > 1) {
    static const string message = "\n\n"
      "warning:\n"
      "  The use of UA5Handler for events with multiple hard subprocesses\n"
      "  is probably not intended as it applies two different models\n" 
      "  of the underlying event at the same time.\n" 
      "  UA5Handler can be disabled in the input files with\n"
      "  'set stdCluHadHandler:UnderlyingEventHandler NULL'\n";
    // here we should really ask the event handler 
    // for the name of the hadronization handler.
    cerr << message 
	 << "\n  This message can be disabled with\n  'set "
	 << fullName() << ":Warning No'\n\n";
    Throw<Exception>() << message << Exception::warning;
    _needWarning = false;
  }

  // create a new step for the products
  StepPtr newstep = newStep();
  // Constants that should not need changing.
  static const Length vclx = 4e-12*mm; 
  static const Length vcly = 4e-12*mm; 
  static const Length vclz = 4e-12*mm; 
  static const Length vclt = 4e-12*mm;
  // Find the first two clusters
  // Lets find the clusters, set the partons inside to be on shell and no momentum
  tClusterPtr clu[2];
  tPVector::const_iterator it;
  unsigned int i = 0;
  for(it = tagged.begin(); it!=tagged.end(); ++it) {
    if((*it)->id() != ParticleID::Cluster) continue;
    clu[i] = dynamic_ptr_cast<ClusterPtr>(*it);
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
	insertParticle(cluster,newstep,false);
      }
      // don't to the rest
      return;
    }
  useMe();
  // and their cm
  Lorentz5Momentum cm = clu[0]->momentum() + clu[1]->momentum();
  Energy theCM = cm.mass();
  // softCM = sqrt(S) with optional enhancement for multiplicity only 
  // (name of variable not very reasonable any more...)
  Lorentz5Momentum pcm=
    ch.currentEvent()->incoming().first ->momentum()+
    ch.currentEvent()->incoming().second->momentum();
  Energy softCM = _enhanceCM*pcm.m();
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
     Energy sumMasses = ZERO;
     // delete the particles from the previous attempt if needed
     if(ntry > 1) clusters.clear();
     int numCharge = 0;
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
       Lorentz5Momentum cp(ZERO,ZERO,ZERO,newMass,newMass);
       cluster->set5Momentum(cp);
       // Now the gaussian distribution of the x,y,z components,
       // and a time component given by
       // sqrt(vx^2+vy^2+vz^2) - vclt*log(r)
       Length vx = gaussDistribution(0.*mm,vclx);
       Length vy = gaussDistribution(0.*mm,vcly);
       Length vz = gaussDistribution(0.*mm,vclz);
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
       if(numberCluster == 0) theMult = nppbar+abs(totalcharge);
       if(numberCluster == 1) theMult += abs(totalcharge);
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
      insertParticle(cluster,newstep,false);
    }
    // don't to the rest
    return;
  }
  // Now generate momentum of the clusters
  generateMomentum(clu[0],clu[1],clusters, theCM, cm);
  // insert the particles into the event record
  for(unsigned int ix=0;ix<clusters.size();++ix) {
    clu[0]->addChild(clusters[ix]);
    clu[1]->addChild(clusters[ix]);
    newstep->addDecayProduct(clusters[ix]);
  }
  for(unsigned int ix=0;ix<clusters.size();++ix) {
    insertParticle(clusters[ix],newstep,false);
  }
}

// Generate the momentum. This is based on the procedure of
// Jadach from Computer Physics Communications 9 (1975) 297-304
void UA5Handler::generateMomentum(tClusterPtr clu1, tClusterPtr clu2,
				  const ClusterVector &clusters, 
				  Energy CME, const Lorentz5Momentum & cm) 
  const {
  // begin with the cylindrical phase space generation described in the paper of Jadach
  generateCylindricalPS(clusters, CME);
  // boost momentum of incoming cluster along z axis to cluster cmf frame
  if(clu2->momentum().z()>ZERO) swap(clu1,clu2);
  LorentzMomentum bmp = clu1->momentum();
  bmp = bmp.boost(cm.findBoostToCM());
  // Rotation to put bmp on the z axis
  LorentzRotation R = rotate(bmp);
  // We need to use the inverse
  R = R.inverse();
  Boost boostv=cm.boostVector();
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
      cerr << "Old p = " << oldP/GeV << endl << "New p = " << newP/GeV << endl;
      cerr << "This is from cluster " << *clusters[i] << " and has z component > energy" << endl;
      cerr << "Cluster 0 is = " << *clusters[0] << endl;
      cerr << "Original momentum = " << origP/GeV << endl;
      cerr << "The cm vector is = " << cm/GeV << endl;
      cerr << "Mass error of original momentum is " << origP.massError() << endl;
      throw Veto();
    }
  }
}

void UA5Handler::generateCylindricalPS(const ClusterVector &clusters, Energy CME) const {
  assert(clusters.size()>=2);
  double alog = log(CME*CME/GeV2);
  unsigned int ncl = clusters.size();
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
  // generate the momenta
  double eps = 1e-10/double(ncl);
  vector<double> xi(ncl);
  unsigned int its(0);
  Energy sum1(ZERO);
  double yy(0.);
  while(its < _maxtries) {
    ++its;
    Energy sumx = ZERO;
    Energy sumy = ZERO;
    for(unsigned int i = 0; i<ncl; ++i) {
      // Generate the pt given the parameter slope
      Energy   pt = randExt(clusters[i]->mass(), slope[i]);
      Energy2 ptp = pt*pt - sqr(mom[i].mass());
      if(ptp <= ZERO) pt = - sqrt(-ptp);
      else pt = sqrt(ptp);
      // randomize azimuth
      Energy px,py;
      randAzm(pt,px,py);
      // set transverse momentum
      mom[i].setX(px);
      mom[i].setY(py);
      sumx += px;
      sumy += py;
    }
    sumx /= ncl;
    sumy /= ncl;
    // find the sum of the transverse mass
    Energy sumtm=ZERO;
    for(unsigned int ix = 0; ix<ncl; ++ix) {
      mom[ix].setX(mom[ix].x()-sumx);
      mom[ix].setY(mom[ix].y()-sumy);
      Energy2 pt2 = mom[ix].perp2();
      // Use the z component of the clusters momentum for temporary storage
      mom[ix].setZ(sqrt(pt2+mom[ix].mass2()));
      sumtm += mom[ix].z();
    }
    // if transverse mass greater the CMS try again
    if(sumtm > CME) continue;
    for(unsigned int i = 0; i<ncl; i++) xi[i] = randUng(0.6,1.0);
    // sort into ascending order
    sort(xi.begin(), xi.end());
    double ximin = xi[0];
    double ximax = xi.back()-ximin;
    xi[0] = 0.;
    for(unsigned int i = ncl-2; i>=1; i--) xi[i+1] = (xi[i]-ximin)/ximax;
    xi[1] = 1.;
    yy= log(CME*CME/(mom[0].z()*mom[1].z()));
    bool suceeded=false;
    Energy sum2,sum3,sum4;
    for(unsigned int j = 0; j<10; j++) {
      sum1 = sum2 = sum3 = sum4 = ZERO;
      for(unsigned int i = 0; i<ncl; i++) {
        Energy tm = mom[i].z();
        double ex = exp(yy*xi[i]);
        sum1 += tm*ex;
	sum2 += tm/ex;
        sum3 += (tm*ex)*xi[i];
	sum4 += (tm/ex)*xi[i];  
      }
      double fy = alog-log(sum1*sum2/GeV2);
      double dd = (sum3*sum2 - sum1*sum4)/(sum1*sum2);
      double dyy = fy/dd;
      if(abs(dyy/yy) < eps) { 
	yy += dyy; 
	suceeded=true; 
	break;
      }
      yy += dyy;
    }
    if(suceeded) break;
    if(its > 100) eps *= 10.;
  }
  if(its==_maxtries) 
    throw Exception() << "Can't generate soft underlying event in "
		      << "UA5Handler::generateCylindricalPS"
		      << Exception::eventerror;
  double zz = log(CME/sum1);
  for(unsigned int i = 0; i<ncl; i++) {
    Energy tm = mom[i].z();
    double E1 = exp(zz + yy*xi[i]);
    mom[i].setZ(0.5*tm*(1./E1-E1));
    mom[i].setE( 0.5*tm*(1./E1+E1));
    clusters[i]->set5Momentum(mom[i]);
  }
}
