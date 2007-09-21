// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KtTildeMeasure class.
//

#include "KtTildeMeasure.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/EventRecord/Particle.h"

#include "ThePEG/Utilities/Throw.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "KtTildeMeasure.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <cassert>

#include "Herwig++/Shower/CKKW/Clustering/EmitterSpectatorConfiguration.h"

// OS X workaround for isnan

extern "C" int isnan(double) throw();

using namespace Herwig;

KtTildeMeasure::~KtTildeMeasure() {}

void KtTildeMeasure::persistentOutput(PersistentOStream & ) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void KtTildeMeasure::persistentInput(PersistentIStream & , int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<KtTildeMeasure> KtTildeMeasure::initKtTildeMeasure;
// Definition of the static class description member.

void KtTildeMeasure::Init() {

  static ClassDocumentation<KtTildeMeasure> documentation
    ("This class implements the \\tilde{k}_\\perp "
     "jet measure, which generalizes the Durham measure "
     "to the Herwig++ evolution variable.");

}

// children, parents
Energy2 KtTildeMeasure::scale (const vector<tClusteringParticlePtr>& children,
			       const vector<ClusteringParticlePtr>&,
			       const tClusteringConfigurationPtr& config) {

  assert(children.size() == 3);

  if (children[0]->pData().partonId.state == ClusteringParticleState::initial ||
      children[1]->pData().partonId.state == ClusteringParticleState::initial ||
      children[2]->pData().partonId.state == ClusteringParticleState::initial)
    throw Exception() << "CKKW : KtTildeMeasure : Initial state not implemented yet." << Exception::runerror;

  tEmitterSpectatorConfigurationPtr conf =
    dynamic_ptr_cast<tEmitterSpectatorConfigurationPtr>(config);
  if(!conf) throw Exception() << "CKKW : KtTildeMeasure::scale : expecting EmitterSpectatorConfiguration"
			      << Exception::runerror;

  unsigned int spectator = conf->spectatorBeforeClusteringIndex();

  vector<tClusteringParticlePtr> emissions;

  for(unsigned int i=0; i<3; ++i)
    if (i != spectator) emissions.push_back(children[i]);

  double z = (sudakovBasis().second * emissions[0]->momentum())/
    (sudakovBasis().second * emissions[0]->momentum() + sudakovBasis().second * emissions[1]->momentum());

  Lorentz5Momentum q = emissions[0]->momentum() + emissions[1]->momentum();
  Energy2 q2 = q*q;

  Energy2 mi2 = sqr(getParticleData(conf->emission().first.partonId.PDGId)->mass());
  Energy2 mj2 = sqr(getParticleData(conf->emission().second.partonId.PDGId)->mass());
  Energy2 mij2 = sqr(getParticleData(conf->emitter().partonId.PDGId)->mass());

  Energy2 qt2 = (q2-mij2)/(z*(1.-z));

  Energy2 dm2 = mij2-(mi2+mj2);

  double zpm = zPM (mi2,mj2,mij2);

  Energy2 theScale;

  if (z > zpm)
    theScale = sqr(1.-z)*qt2 - mj2 + (1.-z)*dm2;
  else
    theScale = sqr(z)*qt2 - mi2 + z*dm2;

  return theScale;

}

double KtTildeMeasure::z (const vector<tClusteringParticlePtr>& children,
			  const vector<ClusteringParticlePtr>&,
			  const tClusteringConfigurationPtr& config) {

  assert(children.size() == 3);

  tEmitterSpectatorConfigurationPtr conf =
    dynamic_ptr_cast<tEmitterSpectatorConfigurationPtr>(config);
  if(!conf) throw Exception() << "CKKW : KtTildeMeasure::scale : expecting EmitterSpectatorConfiguration"
			      << Exception::runerror;

  unsigned int spectator = conf->spectatorBeforeClusteringIndex();

  vector<tClusteringParticlePtr> emissions;

  for(unsigned int i=0; i<3; ++i)
    if (i != spectator) emissions.push_back(children[i]);

  return (sudakovBasis().second * emissions[0]->momentum())/
    (sudakovBasis().second * emissions[0]->momentum() + sudakovBasis().second * emissions[1]->momentum());

}

Energy2 KtTildeMeasure::minResolvableScale (long, bool) {
  return resolution();
}

bool KtTildeMeasure::canHandle (const IdList& ids, bool) {

  assert(ids.size() == 3);
  return true;

}

double KtTildeMeasure::zPM (Energy2 mi2, Energy2 mj2, Energy2 mij2) {

  Energy2 dm2 = mij2-(mi2+mj2);

  if (mi2 == 0.*GeV2 && mj2 == 0.*GeV2 && mij2 == 0.*GeV2) return .5;
  if (mi2 == mj2) return .5;

  if(dm2 == 0.*GeV2) 
    return (mi2+resolution())/(mi2+resolution()+sqrt((mi2+resolution())*(mj2+resolution())));

  // for the general solution, disable unit checking

  double dmi2 = mi2/GeV2;
  double dmj2 = mj2/GeV2;
  double dmij2 = mij2/GeV2;
  double q02 = resolution()/GeV2;

  double general = -(-2*dmi2 + 3*dmij2 - 4*dmj2)/(6.*(dmi2 - dmij2 + dmj2)) - 
   (-pow(-2*dmi2 + 3*dmij2 - 4*dmj2,2) + 6*(dmi2 - dmij2 + dmj2)*(-dmi2 - dmij2 + dmj2 - 2*q02))/
    (3.*pow(2,2./3.)*(dmi2 - dmij2 + dmj2)*
      pow(-56*pow(dmi2,3) + 90*pow(dmi2,2)*dmij2 - 36*dmi2*pow(dmij2,2) - 
        48*pow(dmi2,2)*dmj2 + 36*pow(dmij2,2)*dmj2 + 48*dmi2*pow(dmj2,2) - 
        90*dmij2*pow(dmj2,2) + 56*pow(dmj2,3) - 36*pow(dmi2,2)*q02 + 36*dmi2*dmij2*q02 - 
        36*dmij2*dmj2*q02 + 36*pow(dmj2,2)*q02 + 
        sqrt(4*pow(-pow(-2*dmi2 + 3*dmij2 - 4*dmj2,2) + 
             6*(dmi2 - dmij2 + dmj2)*(-dmi2 - dmij2 + dmj2 - 2*q02),3) + 
          pow(-56*pow(dmi2,3) + 90*pow(dmi2,2)*dmij2 - 36*dmi2*pow(dmij2,2) - 
            48*pow(dmi2,2)*dmj2 + 36*pow(dmij2,2)*dmj2 + 48*dmi2*pow(dmj2,2) - 
            90*dmij2*pow(dmj2,2) + 56*pow(dmj2,3) - 36*pow(dmi2,2)*q02 + 
            36*dmi2*dmij2*q02 - 36*dmij2*dmj2*q02 + 36*pow(dmj2,2)*q02,2)),1./3.))
     + pow(-56*pow(dmi2,3) + 90*pow(dmi2,2)*dmij2 - 36*dmi2*pow(dmij2,2) - 
      48*pow(dmi2,2)*dmj2 + 36*pow(dmij2,2)*dmj2 + 48*dmi2*pow(dmj2,2) - 
      90*dmij2*pow(dmj2,2) + 56*pow(dmj2,3) - 36*pow(dmi2,2)*q02 + 36*dmi2*dmij2*q02 - 
      36*dmij2*dmj2*q02 + 36*pow(dmj2,2)*q02 + 
      sqrt(4*pow(-pow(-2*dmi2 + 3*dmij2 - 4*dmj2,2) + 
           6*(dmi2 - dmij2 + dmj2)*(-dmi2 - dmij2 + dmj2 - 2*q02),3) + 
        pow(-56*pow(dmi2,3) + 90*pow(dmi2,2)*dmij2 - 36*dmi2*pow(dmij2,2) - 
          48*pow(dmi2,2)*dmj2 + 36*pow(dmij2,2)*dmj2 + 48*dmi2*pow(dmj2,2) - 
          90*dmij2*pow(dmj2,2) + 56*pow(dmj2,3) - 36*pow(dmi2,2)*q02 + 36*dmi2*dmij2*q02 - 
          36*dmij2*dmj2*q02 + 36*pow(dmj2,2)*q02,2)),1./3.)/
    (6.*pow(2,1./3.)*(dmi2 - dmij2 + dmj2));

  if (general <= 0. || general >= 1. || isnan(general))
    throw Exception() << "CKKW : KtTildeMeasure::zPM : General solution gave unphysical result."
		      << Exception::eventerror;

  return general;

}

double KtTildeMeasure::showerJacobian (const IdList& ids, Energy2, double z, bool initial) {

  double j;

  if(!initial) {

    Energy2 mij2 = sqr(getParticleData(ids[0])->mass());
    Energy2 mi2 = sqr(getParticleData(ids[1])->mass());
    Energy2 mj2 = sqr(getParticleData(ids[2])->mass());

    double zpm = zPM (mi2,mj2,mij2);

    if (z>zpm)
      j= 1./sqr(1.-z);
    else
      j= 1./sqr(z);

  } else {

    // add initial state version here

    j = 0.;

  }

  return j;

}

pair<Energy2,double> KtTildeMeasure::invertClustering (const IdList& ids, Energy2 q2, double z, bool initial) {

  Energy2 qt2 = 0.*GeV2;

  if (!initial) {

    Energy2 mij2 = sqr(getParticleData(ids[0])->mass());
    Energy2 mi2 = sqr(getParticleData(ids[1])->mass());
    Energy2 mj2 = sqr(getParticleData(ids[2])->mass());

    Energy2 dm2 = mij2-(mi2+mj2);

    double zpm = zPM (mi2,mj2,mij2);
    
    if (z>zpm)
      qt2 = (q2+mj2-(1.-z)*dm2)/sqr(1.-z);
    else
      qt2 = (q2+mi2-z*dm2)/sqr(z);

  } else {

    // add initial state version here

  }

  return make_pair(qt2,z);

}


pair<double,double> KtTildeMeasure::zLimits (Energy2 q2, Energy2 Q2, const IdList& ids, bool initial) {

  double zmin = 0.;
  double zmax = 1.;

  if (!initial) {

    Energy2 mij2 = sqr(getParticleData(ids[0])->mass());
    Energy2 mi2 = sqr(getParticleData(ids[1])->mass());
    Energy2 mj2 = sqr(getParticleData(ids[2])->mass());
    
    Energy2 dm2 = mij2-(mi2+mj2);

    double zpm = zPM(mi2,mj2,mij2);

    Energy2 qt2maxi = (Q2+mi2-zpm*dm2)/sqr(zpm);
    Energy2 qt2maxj = (Q2+mj2-(1.-zpm)*dm2)/sqr(1.-zpm);
    
    double rhoi = dm2/(2.*qt2maxi);
    double rhoj = dm2/(2.*qt2maxj);

    zmin = sqrt(sqr(rhoi)+(mi2+q2)/qt2maxi)-rhoi;
    zmax = 1.- (sqrt(sqr(rhoj)+(mj2+q2)/qt2maxj)-rhoj);

  } else {

    Throw<InitException>()
      << "CKKW : KtTildeMeasure : Initial state not implemented yet.";

    // add initial state version here

  }

  return make_pair(zmin,zmax);

}

bool KtTildeMeasure::resolvable (const IdList& ids, Energy2 q2, Energy2 Q2, double z, bool) {
  pair<double,double> zlims = zLimits(q2,Q2,ids);
  return z > zlims.first && z < zlims.second && q2 > resolution();
}

bool KtTildeMeasure::resolvable (tcShowerParticlePtr, const Branching& br, bool initial) {

  // do not veto non-QCD branchings

  if (br.kinematics->splittingFn()->interactionType() != ShowerIndex::QCD)
    return false;

  if(!initial) {

    Energy2 mij2 = sqr(getParticleData(br.ids[0])->mass());
    Energy2 mi2 = sqr(getParticleData(br.ids[1])->mass());
    Energy2 mj2 = sqr(getParticleData(br.ids[2])->mass());
    
    Energy2 dm2 = mij2-(mi2+mj2);

    Energy2 qt2 = sqr(br.kinematics->scale());
    double z = br.kinematics->z();
    
    return qt2 > (resolution()+mi2-z*dm2)/sqr(z) && qt2 > (resolution()+mj2-(1.-z)*dm2)/sqr(1.-z);

  } else {

    throw Exception() << "CKKW : KtTildeMeasure : Initial state not implemented yet." << Exception::runerror;

    // add initial state version here

  }

}

Energy2 KtTildeMeasure::emitterMass (long id1, long id2, bool& is) {
  Energy2 m2;
  is = false;
  if ((id1+id2 == 0 && abs(id1) < 7) || (id1 == id2 && id1 == 21)) {
    m2 = 0.*GeV2;
    is = true;
  }
  if (abs(id1)<7 && id2 == 21) {
    m2 = sqr(getParticleData(id1)->mass());
    is = true;
  }
  if (abs(id2)<7 && id1 == 21) {
    m2 = sqr(getParticleData(id2)->mass());
    is = true;
  }
  return m2;
}


bool KtTildeMeasure::resolvable (const vector<pair<PPtr,bool> >& br) {

  assert(br.size() == 3);

  // abort on coloured initial state particles, not yet (fully) implemented
  for (vector<pair<PPtr,bool> >::const_iterator p = br.begin(); p != br.end(); ++p)
    if (p->second && p->first->coloured())
      throw Exception() << "CKKW : KtTildeMeasure : Initial state not implemented yet."
			<< Exception::runerror;

  // if not a qcd branching, we consider it in any case resolvable

  bool res = true;

  Energy2 mij2;
  Energy2 mi2;
  Energy2 mj2;
  Energy2 dm2;
  Energy2 qt2;

  Lorentz5Momentum n;
  Lorentz5Momentum q;
  double z; double zpm;
  bool isQCD;

  // (0 1) 2

  // fs clustering, only consider coloured spectator

  if (!br[0].second && !br[1].second && br[2].first->coloured()) {

    mij2 = emitterMass(br[0].first->id(),br[1].first->id(),isQCD);
    if (isQCD) {
      n = br[2].first->momentum(); n.setMass(0.*GeV); n.rescaleEnergy();
      q = br[0].first->momentum()+br[1].first->momentum();

      mi2 = sqr(getParticleData(br[0].first->id())->mass());
      mj2 = sqr(getParticleData(br[1].first->id())->mass());
      dm2 = mij2-(mi2+mj2);
      
      z = (n * br[0].first->momentum()) / (n*q);
      zpm = zPM(mi2,mj2,mij2);
      
      qt2 = (q*q-mij2)/(z*(1.-z));
      
      if (z > zpm)
	res &= sqr(1.-z)*qt2 - mj2 + (1.-z)*dm2 > resolution();
      else
	res &= sqr(z)*qt2 - mi2 + z*dm2 > resolution();

    }

  }

  else if (br[0].second && br[1].second) return true;

  else if ((br[0].second || br[1].second) && br[2].first->coloured()) {

    // add initial state here

  }

  // (0 2) 1

  // fs clustering

  if (!br[0].second && !br[2].second && br[1].first->coloured()) {

    mij2 = emitterMass(br[0].first->id(),br[2].first->id(),isQCD);
    if (isQCD) {
      n = br[1].first->momentum(); n.setMass(0.*GeV); n.rescaleEnergy();
      q = br[0].first->momentum()+br[2].first->momentum();
      
      mi2 = sqr(getParticleData(br[0].first->id())->mass());
      mj2 = sqr(getParticleData(br[2].first->id())->mass());
      dm2 = mij2-(mi2+mj2);
      
      z = (n * br[0].first->momentum()) / (n*q);
      zpm = zPM(mi2,mj2,mij2);
      
      qt2 = (q*q-mij2)/(z*(1.-z));

      if (z > zpm)
	res &= sqr(1.-z)*qt2 - mj2 + (1.-z)*dm2 > resolution();
      else
	res &= sqr(z)*qt2 - mi2 + z*dm2 > resolution();

    }

  }

  else if (br[0].second && br[2].second) return true;

  else if ((br[0].second || br[2].second) && br[1].first->coloured()) {

    // add initial state here

  }

  // (1 2) 0

  // fs clustering

  if (!br[1].second && !br[2].second && br[0].first->coloured()) {

    mij2 = emitterMass(br[1].first->id(),br[2].first->id(),isQCD);
    if (isQCD) {
      n = br[0].first->momentum(); n.setMass(0.*GeV); n.rescaleEnergy();
      q = br[1].first->momentum()+br[2].first->momentum();
      
      mi2 = sqr(getParticleData(br[1].first->id())->mass());
      mj2 = sqr(getParticleData(br[2].first->id())->mass());
      dm2 = mij2-(mi2+mj2);
      
      z = (n * br[1].first->momentum()) / (n*q);
      zpm = zPM(mi2,mj2,mij2);
      
      qt2 = (q*q-mij2)/(z*(1.-z));
      
      if (z > zpm)
	res &= sqr(1.-z)*qt2 - mj2 + (1.-z)*dm2 > resolution();
      else
	res &= sqr(z)*qt2 - mi2 + z*dm2 > resolution();

    }
  }

  else if (br[1].second && br[2].second) return true;

  else if ((br[1].second || br[2].second) && br[0].first->coloured()) {

    // add initial state here

  }

  return res;

}

bool KtTildeMeasure::hardScales (const vector <tClusteringParticlePtr>& hard) {

  map<tClusteringParticlePtr, tClusteringParticlePtr> partnerMap;

  for (vector <tClusteringParticlePtr>::const_iterator p = hard.begin();
       p != hard.end(); ++p) {

    // only consider coloured particles
    if ((**p).pData().colour != (**p).pData().antiColour) {
            
      // look up the partners

      // have we already got one?
      if (partnerMap.find(*p) == partnerMap.end()) {

	tClusteringParticlePtr cPartner;
	tClusteringParticlePtr acPartner;

	for (vector <tClusteringParticlePtr>::const_iterator q = hard.begin();
	     q != hard.end(); ++q)
	  if (*p != *q) {
	    // if p is an outgoing octet
	    if ((**p).pData().colour != 0 && (**p).pData().antiColour != 0 &&
		(**p).pData().colour != (**p).pData().antiColour &&
		(**p).pData().partonId.state == ClusteringParticleState::final) {
	      if ((**p).pData().colour == (**q).pData().antiColour &&
		  (**q).pData().partonId.state == ClusteringParticleState::final)
		cPartner = *q;
	      if ((**p).pData().antiColour == (**q).pData().colour &&
		  (**q).pData().partonId.state == ClusteringParticleState::final)
		acPartner = *q; 
	      if ((**p).pData().colour == (**q).pData().colour &&
		  (**q).pData().partonId.state == ClusteringParticleState::initial)
		cPartner = *q;
	      if ((**p).pData().antiColour == (**q).pData().antiColour &&
		  (**q).pData().partonId.state == ClusteringParticleState::initial)
		acPartner = *q; 
	    }

	    // if p is an incoming octet
	    if ((**p).pData().colour != 0 && (**p).pData().antiColour != 0 &&
		(**p).pData().colour != (**p).pData().antiColour &&
		(**p).pData().partonId.state == ClusteringParticleState::initial) {
	      if ((**p).pData().colour == (**q).pData().antiColour &&
		  (**q).pData().partonId.state == ClusteringParticleState::initial)
		cPartner = *q;
	      if ((**p).pData().antiColour == (**q).pData().colour &&
		  (**q).pData().partonId.state == ClusteringParticleState::initial)
		acPartner = *q; 
	      if ((**p).pData().colour == (**q).pData().colour &&
		  (**q).pData().partonId.state == ClusteringParticleState::final)
		cPartner = *q;
	      if ((**p).pData().antiColour == (**q).pData().antiColour &&
		  (**q).pData().partonId.state == ClusteringParticleState::final)
		acPartner = *q; 
	    }

	    // outgoing triplet
	    if ((**p).pData().colour != 0 && (**p).pData().antiColour == 0 &&
		(**p).pData().partonId.state == ClusteringParticleState::final) {
	      if ((**p).pData().colour == (**q).pData().antiColour &&
		  (**q).pData().partonId.state == ClusteringParticleState::final)
		cPartner = *q;
	      if ((**p).pData().colour == (**q).pData().colour &&
		  (**q).pData().partonId.state == ClusteringParticleState::initial)
		cPartner = *q;
	    }

	    // incoming triplet
	    if ((**p).pData().colour != 0 && (**p).pData().antiColour == 0 &&
		(**p).pData().partonId.state == ClusteringParticleState::initial) {
	      if ((**p).pData().colour == (**q).pData().antiColour &&
		  (**q).pData().partonId.state == ClusteringParticleState::initial)
		cPartner = *q;
	      if ((**p).pData().colour == (**q).pData().colour &&
		  (**q).pData().partonId.state == ClusteringParticleState::final)
		cPartner = *q;
	    }

	    // outgoing anti triplet
	    if ((**p).pData().colour == 0 && (**p).pData().antiColour != 0 &&
		(**p).pData().partonId.state == ClusteringParticleState::final) {
	      if ((**p).pData().antiColour == (**q).pData().colour &&
		  (**q).pData().partonId.state == ClusteringParticleState::final)
		acPartner = *q;
	      if ((**p).pData().antiColour == (**q).pData().antiColour &&
		  (**q).pData().partonId.state == ClusteringParticleState::initial)
		acPartner = *q;
	    }

	    // incoming anti triplet
	    if ((**p).pData().colour == 0 && (**p).pData().antiColour != 0 &&
		(**p).pData().partonId.state == ClusteringParticleState::initial) {
	      if ((**p).pData().antiColour == (**q).pData().colour &&
		  (**q).pData().partonId.state == ClusteringParticleState::initial)
		acPartner = *q;
	      if ((**p).pData().antiColour == (**q).pData().antiColour &&
		  (**q).pData().partonId.state == ClusteringParticleState::final)
		acPartner = *q;
	    }

	  }

	if (cPartner && acPartner) {
	  if (UseRandom::rndbool()) {
	    partnerMap.insert(make_pair(*p,cPartner));
	    partnerMap.insert(make_pair(cPartner,*p));
	  } else {
	    partnerMap.insert(make_pair(*p,acPartner));
	    partnerMap.insert(make_pair(acPartner,*p));
	  }
	}

	if (cPartner && !acPartner) {
	  partnerMap.insert(make_pair(*p,cPartner));
	  partnerMap.insert(make_pair(cPartner,*p));
	}

	if (!cPartner && acPartner) {
	  partnerMap.insert(make_pair(*p,acPartner));
	  partnerMap.insert(make_pair(acPartner,*p));
	}

      }

    }

  } // end of partner finding

  // now set the scales

  for (map<tClusteringParticlePtr,tClusteringParticlePtr>::iterator d = partnerMap.begin();
       d != partnerMap.end() ; ++d) {

#ifdef HERWIG_DEBUG_CKKW_EXTREME
    generator()->log() << "KtTildeMeasure::hardScales : considering " << endl;
    d->first->debugDump(generator()->log());
    d->second->debugDump(generator()->log());
#endif

    if (d->first->pData().partonId.state == ClusteringParticleState::final &&
	d->second->pData().partonId.state == ClusteringParticleState::final) {
    
      Lorentz5Momentum q = d->first->momentum() + d->second->momentum();

      Energy2 q2 = q*q;

      Energy2 mi2; mi2 = sqr(getParticleData(d->first->pData().partonId.PDGId)->mass());
      Energy2 mj2; mj2 = sqr(getParticleData(d->second->pData().partonId.PDGId)->mass());
      Energy2 mij2; mij2 = 0. * GeV2;

      double zpm = zPM(mi2,mj2,mij2);

#ifdef HERWIG_DEBUG_CKKW_EXTREME
      generator()->log() << "mi2/GeV2 = " << mi2/GeV2 << " mj2/GeV2 = " << mj2/GeV2
			 << " mij2 = " << mij2/GeV2 << endl
			 << "zpm = " << zpm << " q2/GeV2 = " << q2/GeV2 << endl;
#endif

      d->first->productionScale(zpm*q2/(1.-zpm) - mi2-zpm*(mi2+mj2));

    }

    if (d->first->pData().partonId.state == ClusteringParticleState::initial &&
	d->second->pData().partonId.state == ClusteringParticleState::final) {
      // add intial state scales here
    }

    if (d->first->pData().partonId.state == ClusteringParticleState::final &&
	d->second->pData().partonId.state == ClusteringParticleState::initial) {
      // add initial state scales here
    }

  }

  return true;

}
