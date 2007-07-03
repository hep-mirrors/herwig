// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusterDecayer class.
//

#include "ClusterDecayer.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/Repository/EventGenerator.h>
#include "Herwig++/Utilities/Kinematics.h"
#include "CheckId.h"
#include "Herwig++/Utilities/Smearing.h"
#include "Cluster.h"

using namespace Herwig;

void ClusterDecayer::persistentOutput(PersistentOStream & os) const 
{
  os << _hadronsSelector << _ClDir1 << _ClDir2 
     << _ClSmr1 << _ClSmr2 << _onshell << _masstry;
}

void ClusterDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _hadronsSelector >> _ClDir1 >> _ClDir2
     >> _ClSmr1 >> _ClSmr2 >> _onshell >> _masstry;
}

ClassDescription<ClusterDecayer> ClusterDecayer::initClusterDecayer;
// Definition of the static class description member.


void ClusterDecayer::Init() {

  static ClassDocumentation<ClusterDecayer> documentation
    ("This class is responsible for the two-body decays of normal clusters");

  static Reference<ClusterDecayer,HadronSelector> 
    interfaceHadronSelector("HadronSelector", 
                             "A reference to the HadronSelector object", 
                             &Herwig::ClusterDecayer::_hadronsSelector,
			     false, false, true, false);
  static Parameter<ClusterDecayer,int> 
    interfaceClDir1 ("ClDir1", "cluster direction for non-b quarks",
                     &ClusterDecayer::_ClDir1, 0, 1 , 0 , 1,false,false,false);
  static Parameter<ClusterDecayer,int> 
    interfaceClDir2 ("ClDir2", "cluster direction for b quark",
                     &ClusterDecayer::_ClDir2, 0, 1 , 0 , 1,false,false,false);
  static Parameter<ClusterDecayer,double> 
    interfaceClSmr1 ("ClSmr1", "cluster direction Gaussian smearing for non-b quark",
                     &ClusterDecayer::_ClSmr1, 0, 0.0 , 0.0 , 2.0,false,false,false);
  static Parameter<ClusterDecayer,double> 
    interfaceClSmr2 ("ClSmr2", "cluster direction Gaussian smearing for b quark",
                     &ClusterDecayer::_ClSmr2, 0, 0.0 , 0.0 , 2.0,false,false,false);
  
  static Switch<ClusterDecayer,bool> interfaceOnShell
    ("OnShell",
     "Whether or not the hadrons produced should by on shell or generated using the"
     " mass generator.",
     &ClusterDecayer::_onshell, false, false, false);
  static SwitchOption interfaceOnShellOnShell
    (interfaceOnShell,
     "OnShell",
     "Produce the hadrons on shell",
     false);
  static SwitchOption interfaceOnShellOffShell
    (interfaceOnShell,
     "OffShell",
     "Generate the masses using the mass generator.",
     false);

  static Parameter<ClusterDecayer,unsigned int> interfaceMassTry
    ("MassTry",
     "The number attempts to generate the masses of the hadrons produced"
     " in the cluster decay.",
     &ClusterDecayer::_masstry, 20, 1, 50,
     false, false, Interface::limited);

}


void ClusterDecayer::decay(const StepPtr &pstep) 
  throw(Veto, Stop, Exception) {
  // Loop over all clusters, and if they are not too heavy (that is
  // intermediate clusters that have undergone to fission) or not 
  // too light (that is final clusters that have been already decayed 
  // into single hadron) then decay them into two hadrons.
  ClusterVector clusters; 
  for (ParticleSet::iterator it = pstep->particles().begin();
       it!= pstep->particles().end(); it++) { 
    if((*it)->id() == ExtraParticleID::Cluster) 
      clusters.push_back(dynamic_ptr_cast<ClusterPtr>(*it));
  }
  for (ClusterVector::const_iterator it = clusters.begin();
	 it != clusters.end(); ++it) {
    if ((*it)->isAvailable() && !(*it)->isStatusFinal() 
	&& (*it)->isReadyToDecay()) {   
      pair<PPtr,PPtr> prod = decayIntoTwoHadrons(*it);
      pstep->addDecayProduct(*it,prod.first);
      pstep->addDecayProduct(*it,prod.second);
    }
  }
}


pair<PPtr,PPtr> ClusterDecayer::decayIntoTwoHadrons(tClusterPtr ptr) 
  throw(Veto, Stop, Exception) {
  using Constants::pi;
  // To decay the cluster into two hadrons one must distinguish between
  // constituent quarks (or diquarks) that originate from perturbative
  // processes (hard process or parton shower) from those that are generated
  // by the non-perturbative gluon splitting or from fission of heavy clusters.
  // In the latter case the two body decay is assumed to be *isotropic*.
  // In the former case instead, if proper flag are activated, the two body 
  // decay is assumed to "remember" the direction of the constituents inside 
  // the cluster, in the cluster frame. The actual smearing of the hadron 
  // directions around the direction of the constituents, in the cluster 
  // frame, can be different between non-b hadrons and b-hadrons, but is given
  // by the same functional form:
  //          cosThetaSmearing = 1 + smearFactor * log( rnd() )
  // (repeated until cosThetaSmearing > -1)
  // where the factor smearFactor is different between b- and non-b hadrons.
  //
  // We need to require (at least at the moment, maybe in the future we 
  // could change it) that the cluster has exactly two components. 
  // If this is not the case, then send a warning because it is not suppose 
  // to happen, and then return.
  if ( ptr->numComponents() != 2 ) {
    generator()->logWarning( Exception("ClusterDecayer::decayIntoTwoHadrons "
				       "***Still cluster with not exactly 2 components*** ", 
				       Exception::warning) );
    return pair<PPtr,PPtr>();
  }
     
  // Extract the id and particle pointer of the two components of the cluster.
  long id1 = 0, id2 = 0;
  tPPtr ptr1 = tPPtr(), ptr2 = tPPtr();
  ptr1 = ptr->particle(0);
  id1 = ptr1->id();
  ptr2 = ptr->particle(1);
  id2 = ptr2->id();

  // Be careful that in these method we use "1" and "2" for the two hadrons
  // produced by the decay of the cluster, whereas "1" and "2" in the 
  // parameters _ClDir1, _ClDir2, _ClSmr1, _ClSmr2 have a different meaning:
  // "1" means non-b quark, and "2" means b quark.
  bool isHad1FlavB    = false;
  int cluDirHad1      = _ClDir1;
  double cluSmearHad1 = _ClSmr1;
  if ( CheckId::hasBeauty(id1) ) {
    isHad1FlavB  = true;
    cluDirHad1   = _ClDir2;
    cluSmearHad1 = _ClSmr2;
  } 
  bool isHad2FlavB    = false;
  int cluDirHad2      = _ClDir1;
  double cluSmearHad2 = _ClSmr1;
  if ( CheckId::hasBeauty(id2) ) {
    isHad2FlavB  = true;
    cluDirHad2   = _ClDir2;
    cluSmearHad2 = _ClSmr2;
  } 

  bool isOrigin1Perturbative = ptr->isPerturbative(0);
  bool isOrigin2Perturbative = ptr->isPerturbative(1);

  // We have to decide which, if any, of the two hadrons will have 
  // the momentum, in the cluster parent frame, smeared around the
  // direction of its constituent (for Had1 is the one pointed by
  // ptr1, and for Had2 is the one pointed by ptr2).
  // This happens only if the flag _ClDirX is 1 and the constituent is
  // perturbative (that is not coming from nonperturbative gluon splitting
  // or cluster fission). In the case that both the hadrons satisfy this
  // two requirements (of course only one must be treated, because the other
  // one will have the momentum automatically fixed by the momentum 
  // conservation) then more priority is given in the case of a b-hadron.
  // Finally, in the case that the two hadrons have same priority, then
  // we choose randomly, with equal probability, one of the two. 

  int priorityHad1 = 0;
  if ( cluDirHad1 == 1  &&  isOrigin1Perturbative ) {
    priorityHad1 = 1;
    if (isHad1FlavB) priorityHad1 = 2;
  }
  int priorityHad2 = 0;
  if ( cluDirHad2 == 1  &&  isOrigin2Perturbative ) {
    priorityHad2 = 1;
    if (isHad2FlavB) priorityHad2 = 2;
  }
  if ( priorityHad2  &&  priorityHad1 == priorityHad2  &&  UseRandom::rndbool() ) {
    priorityHad2 = 0;
  }

  Lorentz5Momentum pClu = ptr->momentum();
  Axis uSmear_v3;
  bool secondHad = false;
  if ( priorityHad1  ||  priorityHad2 ) { 

    double cluSmear;
    Lorentz5Momentum pQ;
    if ( priorityHad1 > priorityHad2 ) {
      pQ = ptr1->momentum();
      cluSmear = cluSmearHad1;
    } else {                                
      pQ = ptr2->momentum();
      cluSmear = cluSmearHad2;
      secondHad = true;
    }
    
    // To find the momenta of the two hadrons in the parent cluster frame
    // we proceed as follows. First, we determine the unit vector parallel
    // to the direction of the constituent in the cluster frame. Then we
    // have to smear such direction using the following prescription:
    //  --- in theta angle w.r.t. such direction (not the z-axis),
    //      the drawing of the cosine of such angle is done via:
    //                 1.0 + cluSmear*log( rnd() )
    //      (repeated until it gives a value above -1.0)
    //  --- in phi angle w.r.t. such direction, the drawing is simply flat.
    // Then, given the direction in the parent cluster frame of one of the
    // two hadrons, it is straighforward to get the momenta of both hadrons
    // (always in the same parent cluster frame).

    pQ.boost( -pClu.boostVector() );    // boost from Lab to Cluster frame 
    Axis u_v3 = unitVector(pQ.vect()); 
    uSmear_v3 = u_v3;
    if ( cluSmear > 0.001 ) {           // skip if cluSmear is too small
      double cosThetaQ = pQ.cosTheta();    
      double sinThetaQ = sqrt( 1.0 - cosThetaQ*cosThetaQ );
      double cosSmear;
      do {
	cosSmear = 1.0 + cluSmear*log( UseRandom::rnd() );            
      } while ( cosSmear < -1.0 );
      double sinSmear = sqrt( 1.0 - cosSmear*cosSmear );
      // Determine now the theta angle of the smeared direction 
      // w.r.t. to the usual x-y-z axes (rather than w.r.t. u_v3)
      double cosTheta = cosThetaQ*cosSmear - sinThetaQ*sinSmear;
      uSmear_v3.setTheta( acos( cosTheta ) );
      uSmear_v3.rotate( UseRandom::rnd( -pi , pi ) , u_v3 );  // smear in phi around u_v3
    }
  } 
  else {
    
    // Isotropic decay: flat in cosTheta and phi. 
    uSmear_v3 = Axis(1.0, 0.0, 0.0);  // just to set the rho to 1
    uSmear_v3.setTheta( acos( UseRandom::rnd( -1.0 , 1.0 ) ) );
    uSmear_v3.setPhi( UseRandom::rnd( -pi , pi ) );   
    
  }

  pair<tcPDPtr,tcPDPtr> dataPair = _hadronsSelector->chooseHadronPair(ptr->mass(),
								      id1,id2);
  if(dataPair.first  == tcPDPtr() || 
     dataPair.second == tcPDPtr()) return pair<PPtr,PPtr>();

  // Create the two hadron particle objects with the specified id.
  PPtr ptrHad1,ptrHad2;
  // produce the hadrons on mass shell
  if(_onshell) {
    ptrHad1 = dataPair.first ->produceParticle(dataPair.first ->mass());
    ptrHad2 = dataPair.second->produceParticle(dataPair.second->mass());
  }
  // produce the hadrons with mass given by the mass generator
  else {
    unsigned int ntry(0);
    do {
      ptrHad1 = dataPair.first ->produceParticle();
      ptrHad2 = dataPair.second->produceParticle();
      ++ntry;
    }
    while(ntry<_masstry&&ptrHad1->mass()+ptrHad2->mass()>ptr->mass());
    // if fails produce on shell and issue warning (should never happen??)
    if(ptrHad1->mass()+ptrHad2->mass()>ptr->mass()) {
      generator()->log() << "Failed to produce off-shell hadrons in "
			 << "ClusterDecayer::decayIntoTwoHadrons producing hadrons "
			 << "on-shell" << endl;
      ptrHad1 = dataPair.first ->produceParticle(dataPair.first ->mass());
      ptrHad2 = dataPair.second->produceParticle(dataPair.second->mass());
    }
  }
  
  if (!ptrHad1  ||  !ptrHad2) {
    ostringstream s;
    s << "ClusterDecayer::decayIntoTwoHadrons ***Cannot create the two hadrons***"
      << dataPair.first ->PDGName() << " and " 
      << dataPair.second->PDGName();
    cerr << s.str() << endl;
    generator()->logWarning( Exception(s.str(), Exception::warning) );
  } else {

    Lorentz5Momentum pHad1, pHad2;  // 5-momentum vectors that we have to determine
    if ( secondHad ) uSmear_v3 *= -1.0;

    if (pClu.m() < ptrHad1->mass()+ptrHad2->mass() ) {
      throw Exception() << "Impossible Kinematics in ClusterDecayer::decayIntoTwoHadrons()" 
			<< Exception::eventerror;
    }
    Kinematics::twoBodyDecay(pClu,ptrHad1->mass(),ptrHad2->mass(),uSmear_v3,
			     pHad1,pHad2);
    ptrHad1->set5Momentum(pHad1);
    ptrHad2->set5Momentum(pHad2);

    // Determine the positions of the two children clusters.
    LorentzPoint positionHad1 = LorentzPoint();
    LorentzPoint positionHad2 = LorentzPoint();
    calculatePositions(pClu, ptr->vertex(), pHad1, pHad2, positionHad1, positionHad2);
    ptrHad1->setLabVertex(positionHad1);
    ptrHad2->setLabVertex(positionHad2);
  }
  return pair<PPtr,PPtr>(ptrHad1,ptrHad2);
}


void ClusterDecayer::
calculatePositions(const Lorentz5Momentum &pClu, 
		   const LorentzPoint &positionClu, 
		   const Lorentz5Momentum &, 
		   const Lorentz5Momentum &, 
		   LorentzPoint &positionHad1, 
		   LorentzPoint &positionHad2 ) const 
{
  // First, determine the relative positions of the children hadrons
  // with respect to their parent cluster, in the cluster reference frame,
  // assuming gaussian smearing with width inversely proportional to the 
  // parent cluster mass.
  Length smearingWidth = hbarc / pClu.m();
  LorentzDistance distanceHad1, distanceHad2;
  for ( int i = 0; i < 2; i++ ) {   // i==0 is the first hadron; i==1 is the second one
    for ( int j = 0; j < 4; j++ ) {  // the four components of the LorentzDistance
      double delta;
      while ( ! Smearing::gaussianSmearing( 0.0, smearingWidth/mm, delta ) ) { }
      if ( i == 0 ) {
	// DGRELL fudge!
	distanceHad1 = LorentzDistance(0*mm,0*mm,delta*mm,0*mm);
      } else {
	distanceHad2 = LorentzDistance(0*mm,0*mm,delta*mm,0*mm);	
      }
    }
  }
  // Then, boost such relative positions of the children hadrons,
  // with respect to their parent cluster,
  // from the cluster reference frame to the Lab frame.
  distanceHad1.boost(pClu.boostVector()); 
  distanceHad2.boost(pClu.boostVector());  

  // Finally, determine the absolute positions of the children hadrons
  // in the Lab frame.
  positionHad1 = distanceHad1 + positionClu;
  positionHad2 = distanceHad2 + positionClu;

}

