// -*- C++ -*-
//
// ClusterFissioner.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// Thisk is the implementation of the non-inlined, non-templated member
// functions of the ClusterFissioner class.
//

#include "ClusterFissioner.h"
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Reference.h>
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Interface/Switch.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/PDT/EnumParticles.h>
#include "Herwig/Utilities/Kinematics.h"
#include "CheckId.h"
#include "Cluster.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include <ThePEG/Utilities/DescribeClass.h>

using namespace Herwig;

DescribeClass<ClusterFissioner,Interfaced>
describeClusterFissioner("Herwig::ClusterFissioner","");

ClusterFissioner::ClusterFissioner() :
  _clMaxLight(3.35*GeV),
  _clMaxBottom(3.35*GeV),
  _clMaxCharm(3.35*GeV),
  _clMaxExotic(3.35*GeV),
  _clPowLight(2.0),
  _clPowBottom(2.0),
  _clPowCharm(2.0),
  _clPowExotic(2.0),
  _pSplitLight(1.0),
  _pSplitBottom(1.0),
  _pSplitCharm(1.0),
  _pSplitExotic(1.0),
  _btClM(1.0*GeV),
  _iopRem(1),
  _kappa(1.0e15*GeV/meter)
{}

IBPtr ClusterFissioner::clone() const {
  return new_ptr(*this);
}

IBPtr ClusterFissioner::fullclone() const {
  return new_ptr(*this);
}

void ClusterFissioner::persistentOutput(PersistentOStream & os) const {
  os << _hadronsSelector << ounit(_clMaxLight,GeV)
     << ounit(_clMaxBottom,GeV) << ounit(_clMaxCharm,GeV)
     << ounit(_clMaxExotic,GeV) << _clPowLight << _clPowBottom
     << _clPowCharm << _clPowExotic << _pSplitLight 
     << _pSplitBottom << _pSplitCharm << _pSplitExotic << ounit(_btClM,GeV) 
     << _iopRem  << ounit(_kappa, GeV/meter);
}

void ClusterFissioner::persistentInput(PersistentIStream & is, int) {
  is >> _hadronsSelector >> iunit(_clMaxLight,GeV)
     >> iunit(_clMaxBottom,GeV) >> iunit(_clMaxCharm,GeV)
     >> iunit(_clMaxExotic,GeV) >> _clPowLight >> _clPowBottom
     >> _clPowCharm >> _clPowExotic >> _pSplitLight 
     >> _pSplitBottom >> _pSplitCharm >> _pSplitExotic
     >> iunit(_btClM,GeV) >> _iopRem
     >> iunit(_kappa, GeV/meter);
}

void ClusterFissioner::Init() {

  static ClassDocumentation<ClusterFissioner> documentation
    ("Class responsibles for chopping up the clusters");

  static Reference<ClusterFissioner,HadronSelector> 
    interfaceHadronSelector("HadronSelector", 
                             "A reference to the HadronSelector object", 
                             &Herwig::ClusterFissioner::_hadronsSelector,
			     false, false, true, false);
  
  // ClMax for light, Bottom, Charm and exotic (e.g. Susy) quarks
  static Parameter<ClusterFissioner,Energy>
    interfaceClMaxLight ("ClMaxLight","cluster max mass for light quarks (unit [GeV])",
                    &ClusterFissioner::_clMaxLight, GeV, 3.35*GeV, ZERO, 10.0*GeV,
		    false,false,false);
  static Parameter<ClusterFissioner,Energy>
    interfaceClMaxBottom ("ClMaxBottom","cluster max mass  for b quarks (unit [GeV])",
                    &ClusterFissioner::_clMaxBottom, GeV, 3.35*GeV, ZERO, 10.0*GeV,
		    false,false,false);
  static Parameter<ClusterFissioner,Energy>
    interfaceClMaxCharm ("ClMaxCharm","cluster max mass for c quarks  (unit [GeV])",
                    &ClusterFissioner::_clMaxCharm, GeV, 3.35*GeV, ZERO, 10.0*GeV,
		    false,false,false);
  static Parameter<ClusterFissioner,Energy>
    interfaceClMaxExotic ("ClMaxExotic","cluster max mass  for exotic quarks (unit [GeV])",
                    &ClusterFissioner::_clMaxExotic, GeV, 3.35*GeV, ZERO, 10.0*GeV,
		    false,false,false);
 
 // ClPow for light, Bottom, Charm and exotic (e.g. Susy) quarks
 static Parameter<ClusterFissioner,double>
    interfaceClPowLight ("ClPowLight","cluster mass exponent for light quarks",
                    &ClusterFissioner::_clPowLight, 0, 2.0, 0.0, 10.0,false,false,false);
 static Parameter<ClusterFissioner,double>
    interfaceClPowBottom ("ClPowBottom","cluster mass exponent for b quarks",
                    &ClusterFissioner::_clPowBottom, 0, 2.0, 0.0, 10.0,false,false,false);
 static Parameter<ClusterFissioner,double>
    interfaceClPowCharm ("ClPowCharm","cluster mass exponent for c quarks",
                    &ClusterFissioner::_clPowCharm, 0, 2.0, 0.0, 10.0,false,false,false);
 static Parameter<ClusterFissioner,double>
    interfaceClPowExotic ("ClPowExotic","cluster mass exponent for exotic quarks",
                    &ClusterFissioner::_clPowExotic, 0, 2.0, 0.0, 10.0,false,false,false);

 // PSplit for light, Bottom, Charm and exotic (e.g. Susy) quarks
  static Parameter<ClusterFissioner,double>
    interfacePSplitLight ("PSplitLight","cluster mass splitting param for light quarks",
                    &ClusterFissioner::_pSplitLight, 0, 1.0, 0.0, 10.0,false,false,false);
  static Parameter<ClusterFissioner,double>
    interfacePSplitBottom ("PSplitBottom","cluster mass splitting param for b quarks",
                    &ClusterFissioner::_pSplitBottom, 0, 1.0, 0.0, 10.0,false,false,false);
 static Parameter<ClusterFissioner,double>
    interfacePSplitCharm ("PSplitCharm","cluster mass splitting param for c quarks",
                    &ClusterFissioner::_pSplitCharm, 0, 1.0, 0.0, 10.0,false,false,false);
 static Parameter<ClusterFissioner,double>
    interfacePSplitExotic ("PSplitExotic","cluster mass splitting param for exotic quarks",
                    &ClusterFissioner::_pSplitExotic, 0, 1.0, 0.0, 10.0,false,false,false);


  static Switch<ClusterFissioner,int> interfaceRemnantOption
    ("RemnantOption",
     "Option for the treatment of remnant clusters",
     &ClusterFissioner::_iopRem, 1, false, false);
  static SwitchOption interfaceRemnantOptionSoft
    (interfaceRemnantOption,
     "Soft",
     "Both clusters produced in the fission of the beam cluster"
     " are treated as soft clusters.",
     0);
  static SwitchOption interfaceRemnantOptionHard
    (interfaceRemnantOption,
     "Hard",
     "Only the cluster containing the remnant is treated as a soft cluster.",
     1);
  static SwitchOption interfaceRemnantOptionVeryHard
    (interfaceRemnantOption,
     "VeryHard",
     "Even remnant clusters are treated as hard, i.e. all clusters the same",
     2);

  static Parameter<ClusterFissioner,Energy> interfaceBTCLM
    ("SoftClusterFactor",
     "Parameter for the mass spectrum of remnant clusters",
     &ClusterFissioner::_btClM, GeV, 1.*GeV, 0.1*GeV, 10.0*GeV,
     false, false, Interface::limited);

  
  static Parameter<ClusterFissioner,Tension> interfaceStringTension
    ("StringTension",
     "String tension used in vertex displacement calculation",
     &ClusterFissioner::_kappa, GeV/meter, 
     1.0e15*GeV/meter, ZERO, ZERO,
     false, false, Interface::lowerlim);

}

tPVector ClusterFissioner::fission(ClusterVector & clusters, bool softUEisOn) {
  // return if no clusters
  if (clusters.empty()) return tPVector();

  /*****************
   * Loop over the (input) collection of cluster pointers, and store in 
   * the vector  splitClusters  all the clusters that need to be split
   * (these are beam clusters, if soft underlying event is off, and 
   *  heavy non-beam clusters).
   ********************/
  
  stack<ClusterPtr> splitClusters; 
  for(ClusterVector::iterator it = clusters.begin() ; 
      it != clusters.end() ; ++it) {
    
    /**************
     * Skip 3-component clusters that have been redefined (as 2-component 
     * clusters) or not available clusters. The latter check is indeed 
     * redundant now, but it is used for possible future extensions in which, 
     * for some reasons, some of the clusters found by ClusterFinder are tagged
     * straight away as not available.
     **************/
    if((*it)->isRedefined() || !(*it)->isAvailable()) continue;
    // if the cluster is a beam cluster add it to the vector of clusters
    // to be split or if it is heavy
    if((*it)->isBeamCluster() || isHeavy(*it)) splitClusters.push(*it);
  }    
  tPVector finalhadrons;
  cut(splitClusters, clusters, finalhadrons, softUEisOn);
  return finalhadrons;
}

void ClusterFissioner::cut(stack<ClusterPtr> & clusterStack, 
   			   ClusterVector &clusters, tPVector & finalhadrons, 
			   bool softUEisOn) {
  /**************************************************
   * This method does the splitting of the cluster pointed by  cluPtr
   * and "recursively" by all of its cluster children, if heavy. All of these
   * new children clusters are added (indeed the pointers to them) to the
   * collection of cluster pointers  collecCluPtr. The method works as follows.
   * Initially the vector vecCluPtr contains just the input pointer to the
   * cluster to be split. Then it will be filled "recursively" by all
   * of the cluster's children that are heavy enough to require, in their turn,
   * to be split. In each loop, the last element of the vector vecCluPtr is 
   * considered (only once because it is then removed from the vector).
   * This approach is conceptually recursive, but avoid the overhead of
   * a concrete recursive function. Furthermore it requires minimal changes
   * in the case that the fission of an heavy cluster could produce more
   * than two cluster children as assumed now. 
   *
   * Draw the masses: for normal, non-beam clusters a power-like mass dist
   * is used, whereas for beam clusters a fast-decreasing exponential mass 
   * dist is used instead (to avoid many iterative splitting which could 
   * produce an unphysical large transverse energy from a supposed soft beam
   * remnant process).
   ****************************************/
  // Here we recursively loop over clusters in the stack and cut them
  while (!clusterStack.empty()) {
    // take the last element of the vector
    ClusterPtr iCluster = clusterStack.top(); clusterStack.pop();   
    // split it
    cutType ct = iCluster->numComponents() == 2 ?
      cutTwo(iCluster, finalhadrons, softUEisOn) : 
      cutThree(iCluster, finalhadrons, softUEisOn);

    // There are cases when we don't want to split, even if it fails mass test
    if(!ct.first.first || !ct.second.first) {
      // if an unsplit beam cluster leave if for the underlying event
      if(iCluster->isBeamCluster() && softUEisOn)
	iCluster->isAvailable(false);
      continue;
    }
    // check if clusters
    ClusterPtr one = dynamic_ptr_cast<ClusterPtr>(ct.first.first);
    ClusterPtr two = dynamic_ptr_cast<ClusterPtr>(ct.second.first);
    // is a beam cluster must be split into two clusters
    if(iCluster->isBeamCluster() && (!one||!two) && softUEisOn) {
      iCluster->isAvailable(false);
      continue;
    }

    // There should always be a intermediate quark(s) from the splitting
    assert(ct.first.second && ct.second.second);
    /// \todo sort out motherless quark pairs here. Watch out for 'quark in final state' errors
    iCluster->addChild(ct.first.first);
    //    iCluster->addChild(ct.first.second);
    //    ct.first.second->addChild(ct.first.first);
 
    iCluster->addChild(ct.second.first);
    //    iCluster->addChild(ct.second.second);
    //    ct.second.second->addChild(ct.second.first);
 
    // Sometimes the clusters decay C -> H + C' rather then C -> C' + C''
    if(one) {
      clusters.push_back(one);
      if(one->isBeamCluster() && softUEisOn)
	one->isAvailable(false);
      if(isHeavy(one) && one->isAvailable())
	clusterStack.push(one);
    }
    if(two) {
      clusters.push_back(two);
      if(two->isBeamCluster() && softUEisOn)
	two->isAvailable(false);
      if(isHeavy(two) && two->isAvailable())
	clusterStack.push(two);
    }
  }
}

namespace {

  /**
   *  Check if can't make a hadron from the partons
   */
  bool cantMakeHadron(tcPPtr p1, tcPPtr p2) {
    return ! CheckId::canBeHadron(p1->dataPtr(), p2->dataPtr());
  }
  
  /**
   *  Check if can't make a diquark from the partons
   */
  bool cantMakeDiQuark(tcPPtr p1, tcPPtr p2) {
    long id1 = p1->id(), id2 = p2->id();
    return ! (QuarkMatcher::Check(id1) && QuarkMatcher::Check(id2) && id1*id2>0);
  }
}

ClusterFissioner::cutType
ClusterFissioner::cutTwo(ClusterPtr & cluster, tPVector & finalhadrons,
			 bool softUEisOn) {
  // need to make sure only 2-cpt clusters get here
  assert(cluster->numComponents() == 2);
  tPPtr ptrQ1 = cluster->particle(0);
  tPPtr ptrQ2 = cluster->particle(1);
  Energy Mc = cluster->mass();
  assert(ptrQ1);
  assert(ptrQ2);
  
  // And check if those particles are from a beam remnant
  bool rem1 = cluster->isBeamRemnant(0);
  bool rem2 = cluster->isBeamRemnant(1);
  // workout which distribution to use
  bool soft1(false),soft2(false);
  switch (_iopRem) {
  case 0:
    soft1 = rem1 || rem2;
    soft2 = rem2 || rem1;
    break;
  case 1:
    soft1 = rem1;
    soft2 = rem2;
    break;
  }
  // Initialization for the exponential ("soft") mass distribution.
  static const int max_loop = 1000;
  int counter = 0;
  Energy Mc1 = ZERO, Mc2 = ZERO,m1=ZERO,m2=ZERO,m=ZERO;
  tcPDPtr toHadron1, toHadron2;
  PPtr newPtr1 = PPtr ();
  PPtr newPtr2 = PPtr ();
  bool succeeded = false;
  do
    {
      succeeded = false;
      ++counter;
      
      drawNewFlavour(newPtr1,newPtr2); 
      // check for right ordering
      assert (ptrQ2);
      assert (newPtr2);
      assert (ptrQ2->dataPtr());
      assert (newPtr2->dataPtr());
      if(cantMakeHadron(ptrQ1, newPtr1) || cantMakeHadron(ptrQ2, newPtr2)) {
	swap(newPtr1, newPtr2);
	// check again
	if(cantMakeHadron(ptrQ1, newPtr1) || cantMakeHadron(ptrQ2, newPtr2)) {
	  throw Exception() 
	    << "ClusterFissioner cannot split the cluster ("
	    << ptrQ1->PDGName() << ' ' << ptrQ2->PDGName()
	    << ") into hadrons.\n" << Exception::runerror;
	}
      }
      // Check that new clusters can produce particles and there is enough
      // phase space to choose the drawn flavour
      m1 = ptrQ1->data().constituentMass();
      m2 = ptrQ2->data().constituentMass();
      m  = newPtr1->data().constituentMass();
      // Do not split in the case there is no phase space available
      if(Mc <  m1+m + m2+m) continue;
      // power for splitting
      double exp1=_pSplitLight;
      double exp2=_pSplitLight;

      if     (CheckId::isExotic(ptrQ1->dataPtr())) exp1 = _pSplitExotic;
      else if(CheckId::hasBottom(ptrQ1->dataPtr()))exp1 = _pSplitBottom;
      else if(CheckId::hasCharm(ptrQ1->dataPtr())) exp1 = _pSplitCharm;

      if     (CheckId::isExotic(ptrQ2->dataPtr()))  exp2 = _pSplitExotic;
      else if(CheckId::hasBottom(ptrQ2->dataPtr())) exp2 = _pSplitBottom;
      else if(CheckId::hasCharm(ptrQ2->dataPtr()))  exp2 = _pSplitCharm;

      // If, during the drawing of candidate masses, too many attempts fail 
      // (because the phase space available is tiny)
      /// \todo run separate loop here?
      Mc1 = drawChildMass(Mc,m1,m2,m,exp1,soft1);
      Mc2 = drawChildMass(Mc,m2,m1,m,exp2,soft2);
      if(Mc1 < m1+m || Mc2 < m+m2 || Mc1+Mc2 > Mc) continue;
      /**************************
       * New (not present in Fortran Herwig):
       * check whether the fragment masses  Mc1  and  Mc2  are above the 
       * threshold for the production of the lightest pair of hadrons with the 
       * right flavours. If not, then set by hand the mass to the lightest 
       * single hadron with the right flavours, in order to solve correctly
       * the kinematics, and (later in this method) create directly such hadron
       * and add it to the children hadrons of the cluster that undergoes the
       * fission (i.e. the one pointed by iCluPtr). Notice that in this special
       * case, the heavy cluster that undergoes the fission has one single 
       * cluster child and one single hadron child. We prefer this approach,
       * rather than to create a light cluster, with the mass set equal to
       * the lightest hadron, and let then the class LightClusterDecayer to do 
       * the job to decay it to that single hadron, for two reasons: 
       * First, because the sum of the masses of the two constituents can be, 
       * in this case, greater than the mass of that hadron, hence it would
       * be impossible to solve the kinematics for such two components, and
       * therefore we would have a cluster whose components are undefined.
       * Second, the algorithm is faster, because it avoids the reshuffling
       * procedure that would be necessary if we used LightClusterDecayer
       * to decay the light cluster to the lightest hadron.   
       ****************************/
      toHadron1 = _hadronsSelector->chooseSingleHadron(ptrQ1->dataPtr(), newPtr1->dataPtr(),Mc1);
      if(toHadron1) Mc1 = toHadron1->mass();
      toHadron2 = _hadronsSelector->chooseSingleHadron(ptrQ2->dataPtr(), newPtr2->dataPtr(),Mc2);
      if(toHadron2) Mc2 = toHadron2->mass();
      // if a beam cluster not allowed to decay to hadrons
      if(cluster->isBeamCluster() && (toHadron1||toHadron2) && softUEisOn)
	continue;
      // Check if the decay kinematics is still possible: if not then 
      // force the one-hadron decay for the other cluster as well.
      if(Mc1 + Mc2  >  Mc) {
	if(!toHadron1) {
	  toHadron1 = _hadronsSelector->chooseSingleHadron(ptrQ1->dataPtr(), newPtr1->dataPtr(),Mc-Mc2);
	  if(toHadron1) Mc1 = toHadron1->mass();
	} 
	else if(!toHadron2) {
	  toHadron2 = _hadronsSelector->chooseSingleHadron(ptrQ2->dataPtr(), newPtr2->dataPtr(),Mc-Mc1);
	  if(toHadron2) Mc2 = toHadron2->mass();
	}
      }
      succeeded = (Mc >= Mc1+Mc2);
    }
  while (!succeeded && counter < max_loop);

  if(counter >= max_loop) {
    static const PPtr null = PPtr();
    return cutType(PPair(null,null),PPair(null,null));
  }

  // Determined the (5-components) momenta (all in the LAB frame)
  Lorentz5Momentum pClu = cluster->momentum(); // known
  Lorentz5Momentum p0Q1 = ptrQ1->momentum(); // known (mom Q1 before fission)
  Lorentz5Momentum pClu1, pClu2, pQ1, pQone, pQtwo, pQ2; //unknown
  pClu1.setMass(Mc1);
  pClu2.setMass(Mc2);
  pQ1.setMass(m1);
  pQ2.setMass(m2);
  pQone.setMass(m); 
  pQtwo.setMass(m);
  calculateKinematics(pClu,p0Q1,toHadron1,toHadron2,
		      pClu1,pClu2,pQ1,pQone,pQtwo,pQ2);                // out

  /******************
   * The previous methods have determined the kinematics and positions
   * of C -> C1 + C2. 
   * In the case that one of the two product is light, that means either
   * decayOneHadronClu1 or decayOneHadronClu2 is true, then the momenta
   * of the components of that light product have not been determined,
   * and a (light) cluster will not be created: the heavy father cluster
   * decays, in this case, into a single (not-light) cluster and a
   * single hadron. In the other, "normal", cases the father cluster
   * decays into two clusters, each of which has well defined components.
   * Notice that, in the case of components which point to particles, the
   * momenta of the components is properly set to the new values, whereas
   * we do not change the momenta of the pointed particles, because we 
   * want to keep all of the information (that is the new momentum of a
   * component after the splitting, which is contained in the _momentum
   * member of the Component class, and the (old) momentum of that component
   * before the splitting, which is contained in the momentum of the
   * pointed particle). Please not make confusion of this only apparent
   * inconsistency!
   ********************/
  LorentzPoint posC,pos1,pos2;
  posC = cluster->vertex();
  calculatePositions(pClu, posC, pClu1, pClu2, pos1, pos2);
  cutType rval;
  if(toHadron1) {
    rval.first = produceHadron(toHadron1, newPtr1, pClu1, pos1);
    finalhadrons.push_back(rval.first.first);
  }
  else {
    rval.first = produceCluster(ptrQ1, newPtr1, pClu1, pos1, pQ1, pQone, rem1);
  }
  if(toHadron2) {
    rval.second = produceHadron(toHadron2, newPtr2, pClu2, pos2);
    finalhadrons.push_back(rval.second.first);
  }
  else {
    rval.second = produceCluster(ptrQ2, newPtr2, pClu2, pos2, pQ2, pQtwo, rem2);
  }
  return rval;
}


ClusterFissioner::cutType
ClusterFissioner::cutThree(ClusterPtr & cluster, tPVector & finalhadrons,
			   bool softUEisOn) {
  // need to make sure only 3-cpt clusters get here
  assert(cluster->numComponents() == 3);
  // extract quarks
  tPPtr ptrQ[3] = {cluster->particle(0),cluster->particle(1),cluster->particle(2)};
  assert( ptrQ[0] && ptrQ[1] && ptrQ[2] );
  // find maximum mass pair
  Energy mmax(ZERO);
  Lorentz5Momentum pDiQuark;
  int iq1(-1),iq2(-1);
  Lorentz5Momentum psum;
  for(int q1=0;q1<3;++q1) {
    psum+= ptrQ[q1]->momentum();
    for(int q2=q1+1;q2<3;++q2) {
      Lorentz5Momentum ptest = ptrQ[q1]->momentum()+ptrQ[q2]->momentum();
      ptest.rescaleMass();
      Energy mass = ptest.m();
      if(mass>mmax) {
	mmax = mass;
	pDiQuark = ptest;
	iq1  = q1;
	iq2  = q2;
      }
    }
  }
  // and the spectators
  int iother(-1);
  for(int ix=0;ix<3;++ix) if(ix!=iq1&&ix!=iq2) iother=ix;
  assert(iq1>=0&&iq2>=0&&iother>=0);

  // And check if those particles are from a beam remnant
  bool rem1 = cluster->isBeamRemnant(iq1);
  bool rem2 = cluster->isBeamRemnant(iq2);
  // workout which distribution to use
  bool soft1(false),soft2(false);
  switch (_iopRem) {
  case 0:
    soft1 = rem1 || rem2;
    soft2 = rem2 || rem1;
    break;
  case 1:
    soft1 = rem1;
    soft2 = rem2;
    break;
  }
  // Initialization for the exponential ("soft") mass distribution.
  static const int max_loop = 1000;
  int counter = 0;
  Energy Mc1 = ZERO, Mc2 = ZERO, m1=ZERO, m2=ZERO, m=ZERO;
  tcPDPtr toHadron;
  bool toDiQuark(false);
  PPtr newPtr1 = PPtr(),newPtr2 = PPtr();
  PDPtr diquark;
  bool succeeded = false;
  do {
    succeeded = false;
    ++counter;
    drawNewFlavour(newPtr1,newPtr2);
    // randomly pick which will be (anti)diquark and which a mesonic cluster
    if(UseRandom::rndbool()) {
      swap(iq1,iq2);
      swap(rem1,rem2);
    }
    // check first order
    if(cantMakeHadron(ptrQ[iq1], newPtr1) || cantMakeDiQuark(ptrQ[iq2], newPtr2)) {
      swap(newPtr1,newPtr2);
    }
    // check again
    if(cantMakeHadron(ptrQ[iq1], newPtr1) || cantMakeDiQuark(ptrQ[iq2], newPtr2)) {
      throw Exception()
	<< "ClusterFissioner cannot split the cluster ("
	<< ptrQ[iq1]->PDGName() << ' ' << ptrQ[iq2]->PDGName()
	<< ") into a hadron and diquark.\n" << Exception::runerror;
    }
    // Check that new clusters can produce particles and there is enough
    // phase space to choose the drawn flavour
    m1 = ptrQ[iq1]->data().constituentMass();
    m2 = ptrQ[iq2]->data().constituentMass();
    m  = newPtr1->data().constituentMass();
    // Do not split in the case there is no phase space available
    if(mmax <  m1+m + m2+m) continue;
    // power for splitting
    double exp1(_pSplitLight),exp2(_pSplitLight);
    if     (CheckId::isExotic (ptrQ[iq1]->dataPtr())) exp1 = _pSplitExotic;
    else if(CheckId::hasBottom(ptrQ[iq1]->dataPtr())) exp1 = _pSplitBottom;
    else if(CheckId::hasCharm (ptrQ[iq1]->dataPtr())) exp1 = _pSplitCharm;

    if     (CheckId::isExotic (ptrQ[iq2]->dataPtr())) exp2 = _pSplitExotic;
    else if(CheckId::hasBottom(ptrQ[iq2]->dataPtr())) exp2 = _pSplitBottom;
    else if(CheckId::hasCharm (ptrQ[iq2]->dataPtr())) exp2 = _pSplitCharm;

    // If, during the drawing of candidate masses, too many attempts fail
    // (because the phase space available is tiny)
    /// \todo run separate loop here?
    Mc1 = drawChildMass(mmax,m1,m2,m,exp1,soft1);
    Mc2 = drawChildMass(mmax,m2,m1,m,exp2,soft2);
    if(Mc1 < m1+m || Mc2 < m+m2 || Mc1+Mc2 > mmax) continue;
    // check if need to force meson clster to hadron
    toHadron = _hadronsSelector->chooseSingleHadron(ptrQ[iq1]->dataPtr(), newPtr1->dataPtr(),Mc1);
    if(toHadron) Mc1 = toHadron->mass();
    // check if need to force diquark cluster to be on-shell
    toDiQuark = false;
    diquark = CheckId::makeDiquark(ptrQ[iq2]->dataPtr(), newPtr2->dataPtr());
    if(Mc2 < diquark->constituentMass()) {
      Mc2 = diquark->constituentMass();
      toDiQuark = true;
    }
    // if a beam cluster not allowed to decay to hadrons
    if(cluster->isBeamCluster() && toHadron && softUEisOn)
      continue;
    // Check if the decay kinematics is still possible: if not then
    // force the one-hadron decay for the other cluster as well.
    if(Mc1 + Mc2  >  mmax) {
      if(!toHadron) {
	toHadron = _hadronsSelector->chooseSingleHadron(ptrQ[iq1]->dataPtr(), newPtr1->dataPtr(),mmax-Mc2);
	if(toHadron) Mc1 = toHadron->mass();
      }
      else if(!toDiQuark) {
	Mc2 = _hadronsSelector->massLightestHadron(ptrQ[iq2]->dataPtr(), newPtr2->dataPtr());
	toDiQuark = true;
      }
    }
    succeeded = (mmax >= Mc1+Mc2);
  }
  while (!succeeded && counter < max_loop);
  // check no of tries
  if(counter >= max_loop) return cutType();

  // Determine the (5-components) momenta (all in the LAB frame)
  Lorentz5Momentum p0Q1 = ptrQ[iq1]->momentum();
  // to be determined
  Lorentz5Momentum pClu1(Mc1), pClu2(Mc2), pQ1(m1), pQone(m), pQtwo(m), pQ2(m2);
  calculateKinematics(pDiQuark,p0Q1,toHadron,toDiQuark,
		      pClu1,pClu2,pQ1,pQone,pQtwo,pQ2);
  // positions of the new clusters
  LorentzPoint pos1,pos2;
  Lorentz5Momentum pBaryon = pClu2+ptrQ[iother]->momentum();
  calculatePositions(cluster->momentum(), cluster->vertex(), pClu1, pBaryon, pos1, pos2);
  // first the mesonic cluster/meson
  cutType rval;
   if(toHadron) {
     rval.first = produceHadron(toHadron, newPtr1, pClu1, pos1);
     finalhadrons.push_back(rval.first.first);
   }
   else {
     rval.first = produceCluster(ptrQ[iq1], newPtr1, pClu1, pos1, pQ1, pQone, rem1);
   }
   if(toDiQuark) {
     rem2 |= cluster->isBeamRemnant(iother);
     PPtr newDiQuark = diquark->produceParticle(pClu2);
     rval.second = produceCluster(newDiQuark, ptrQ[iother], pBaryon, pos2, pClu2,
      				  ptrQ[iother]->momentum(), rem2);
   }
   else {
     rval.second = produceCluster(ptrQ[iq2], newPtr2, pBaryon, pos2, pQ2, pQtwo, rem2,
				  ptrQ[iother],cluster->isBeamRemnant(iother));
   }
   cluster->isAvailable(false);
   return rval;
}

ClusterFissioner::PPair 
ClusterFissioner::produceHadron(tcPDPtr hadron, tPPtr newPtr, const Lorentz5Momentum &a,
				const LorentzPoint &b) const {
  PPair rval;
  if(hadron->coloured()) {
    rval.first = (_hadronsSelector->lightestHadron(hadron,newPtr->dataPtr()))->produceParticle();
  }
  else
    rval.first = hadron->produceParticle();
  rval.second = newPtr;
  rval.first->set5Momentum(a);
  rval.first->setVertex(b);
  return rval;
}

ClusterFissioner::PPair ClusterFissioner::produceCluster(tPPtr ptrQ, tPPtr newPtr,
							 const Lorentz5Momentum & a,
				                         const LorentzPoint & b,
							 const Lorentz5Momentum & c,
				                         const Lorentz5Momentum & d,
							 bool isRem,
							 tPPtr spect, bool remSpect) const {
  PPair rval;
  rval.second = newPtr;
  ClusterPtr cluster = !spect ? new_ptr(Cluster(ptrQ,rval.second)) : new_ptr(Cluster(ptrQ,rval.second,spect));
  rval.first = cluster;
  cluster->set5Momentum(a);
  cluster->setVertex(b);
  assert(cluster->particle(0)->id() == ptrQ->id());
  cluster->particle(0)->set5Momentum(c);
  cluster->particle(1)->set5Momentum(d);
  cluster->setBeamRemnant(0,isRem);
  if(remSpect) cluster->setBeamRemnant(2,remSpect);
  return rval;
}

void ClusterFissioner::drawNewFlavour(PPtr& newPtrPos,PPtr& newPtrNeg) const {

  // Flavour is assumed to be only  u, d, s,  with weights
  // (which are not normalized probabilities) given
  // by the same weights as used in HadronsSelector for
  // the decay of clusters into two hadrons. 
  double prob_d = _hadronsSelector->pwtDquark();
  double prob_u = _hadronsSelector->pwtUquark();
  double prob_s = _hadronsSelector->pwtSquark();
  int choice = UseRandom::rnd3(prob_u, prob_d, prob_s);
  long idNew = 0;
  switch (choice) {
  case 0: idNew = ThePEG::ParticleID::u; break;  
  case 1: idNew = ThePEG::ParticleID::d; break;
  case 2: idNew = ThePEG::ParticleID::s; break;
  }
  newPtrPos = getParticle(idNew);
  newPtrNeg = getParticle(-idNew);
  assert (newPtrPos);
  assert(newPtrNeg);
  assert (newPtrPos->dataPtr());
  assert(newPtrNeg->dataPtr());

}

Energy ClusterFissioner::drawChildMass(const Energy M, const Energy m1,
				       const Energy m2, const Energy m,
				       const double expt, const bool soft) const {

  /***************************
   * This method, given in input the cluster mass Mclu of an heavy cluster C,
   * made of consituents of masses m1 and m2, draws the masses Mclu1 and Mclu2
   * of, respectively, the children cluster C1, made of constituent masses m1 
   * and m, and cluster C2, of mass Mclu2 and made of constituent masses m2 
   * and m. The mass is extracted from one of the two following mass 
   * distributions:
   *   --- power-like ("normal" distribution)
   *                        d(Prob) / d(M^exponent) = const
   *       where the exponent can be different from the two children C1 (exp1)
   *       and C2 (exponent2).
   *   --- exponential ("soft" distribution) 
   *                        d(Prob) / d(M^2) = exp(-b*M) 
   *       where b = 2.0 / average.
   * Such distributions are limited below by the masses of
   * the constituents quarks, and above from the mass of decaying cluster C. 
   * The choice of which of the two mass distributions to use for each of the
   * two cluster children is dictated by  iRemnant  (see below).
   * If the number of attempts to extract a pair of mass values that are 
   * kinematically acceptable is above some fixed number (max_loop, see below)
   * the method gives up and returns false; otherwise, when it succeeds, it
   * returns true. 
   *
   * These distributions have been modified from HERWIG:   
   * Before these were:
   *      Mclu1 = m1 + (Mclu - m1 - m2)*pow( rnd(), 1.0/exponent1 );
   * The new one coded here is a more efficient version, same density 
   * but taking into account 'in phase space from' beforehand      
   ***************************/
  // hard cluster
  if(!soft) {
    return pow(UseRandom::rnd(pow((M-m1-m2-m)*UnitRemoval::InvE, expt), 
			      pow(m*UnitRemoval::InvE, expt)), 1./expt
	       )*UnitRemoval::E + m1;
  }
  // Otherwise it uses a soft mass distribution
  else { 
    static const InvEnergy b = 2.0 / _btClM;
    Energy max = M-m1-m2-2.0*m;
    double rmin = b*max;
    rmin = ( rmin < 50 ) ? exp(-rmin) : 0.;
    double r1;
    do {
      r1 = UseRandom::rnd(rmin, 1.0) * UseRandom::rnd(rmin, 1.0);
    }
    while (r1 < rmin);
    return m1 + m - log(r1)/b;
  }
}


void ClusterFissioner::calculateKinematics(const Lorentz5Momentum & pClu,
					   const Lorentz5Momentum & p0Q1,
					   const bool toHadron1,
					   const bool toHadron2,
					   Lorentz5Momentum & pClu1, 
					   Lorentz5Momentum & pClu2, 
					   Lorentz5Momentum & pQ1,
					   Lorentz5Momentum & pQbar,
					   Lorentz5Momentum & pQ,     
					   Lorentz5Momentum & pQ2bar) const {

  /******************
   * This method solves the kinematics of the two body cluster decay:
   *    C (Q1 Q2bar)  --->  C1 (Q1 Qbar)  +  C2 (Q Q2bar)
   * In input we receive the momentum of C, pClu, and the momentum 
   * of the quark Q1 (constituent of C), p0Q1, both in the LAB frame.
   * Furthermore, two boolean variables inform whether the two fission
   * products (C1, C2) decay immediately into a single hadron (in which
   * case the cluster itself is identify with that hadron) and we do 
   * not have to solve the kinematics of the components (Q1,Qbar) for
   * C1 and (Q,Q2bar) for C2.
   * The output is given by the following momenta (all 5-components, 
   * and all in the LAB frame):
   *   pClu1 , pClu2   respectively of   C1 , C2  
   *   pQ1 , pQbar     respectively of   Q1 , Qbar  in  C1
   *   pQ  , pQ2bar    respectively of   Q  , Q2    in  C2
   * The assumption, suggested from the string model, is that, in C frame,
   * C1 and its constituents Q1 and Qbar are collinear, and collinear to
   * the direction of Q1 in C (that is before cluster decay); similarly,
   * (always in the C frame) C2 and its constituents Q and Q2bar are
   * collinear (and therefore anti-collinear with C1,Q1,Qbar).
   * The solution is then obtained by using Lorentz boosts, as follows.
   * The kinematics of C1 and C2 is solved in their parent C frame,
   * and then boosted back in the LAB. The kinematics of Q1 and Qbar
   * is solved in their parent C1 frame and then boosted back in the LAB; 
   * similarly, the kinematics of Q and Q2bar is solved in their parent 
   * C2 frame and then boosted back in the LAB. In each of the three
   * "two-body decay"-like cases, we use the fact that the direction 
   * of the motion of the decay products is known in the rest frame of 
   * their parent. This is obvious for the first case in which the 
   * parent rest frame is C; but it is also true in the other two cases
   * where the rest frames are C1 and C2. This is because C1 and C2 
   * are boosted w.r.t. C in the same direction where their components, 
   * respectively (Q1,Qbar) and (Q,Q2bar) move in C1 and C2 rest frame
   * respectively.      
   * Of course, although the notation used assumed that C = (Q1 Q2bar)
   * where Q1 is a quark and Q2bar an antiquark, indeed everything remain
   * unchanged also in all following cases:
   *  Q1 quark, Q2bar antiquark;           --> Q quark;
   *  Q1 antiquark , Q2bar quark;          --> Q antiquark;  
   *  Q1 quark, Q2bar diquark;             --> Q quark
   *  Q1 antiquark, Q2bar anti-diquark;    --> Q antiquark
   *  Q1 diquark, Q2bar quark              --> Q antiquark
   *  Q1 anti-diquark, Q2bar antiquark;    --> Q quark
   **************************/

  // Calculate the unit three-vector, in the C frame, along which
  // all of the constituents and children clusters move.
  Lorentz5Momentum u(p0Q1);
  u.boost( -pClu.boostVector() );        // boost from LAB to C
  // the unit three-vector is then  u.vect().unit()
						       
  // Calculate the momenta of C1 and C2 in the (parent) C frame first, 
  // where the direction of C1 is u.vect().unit(), and then boost back in the 
  // LAB frame.

  if (pClu.m() < pClu1.mass() + pClu2.mass() ) {
    throw Exception() << "Impossible Kinematics in ClusterFissioner::calculateKinematics() (A)" 
		      << Exception::eventerror;
  }
  Kinematics::twoBodyDecay(pClu, pClu1.mass(), pClu2.mass(), 
			   u.vect().unit(), pClu1, pClu2);

  // In the case that cluster1 does not decay immediately into a single hadron,
  // calculate the momenta of Q1 (as constituent of C1) and Qbar in the
  // (parent) C1 frame first, where the direction of Q1 is u.vect().unit(), 
  // and then boost back in the LAB frame. 
  if(!toHadron1) {
    if (pClu1.m() < pQ1.mass() + pQbar.mass() ) {
      throw Exception() << "Impossible Kinematics in ClusterFissioner::calculateKinematics() (B)" 
			<< Exception::eventerror;
    }
    Kinematics::twoBodyDecay(pClu1, pQ1.mass(), pQbar.mass(), 
			     u.vect().unit(), pQ1, pQbar);
  }

  // In the case that cluster2 does not decay immediately into a single hadron,
  // Calculate the momenta of Q and Q2bar (as constituent of C2) in the
  // (parent) C2 frame first, where the direction of Q is u.vect().unit(), 
  // and then boost back in the LAB frame. 
  if(!toHadron2) {
    if (pClu2.m() < pQ.mass() + pQ2bar.mass() ) {
      throw Exception() << "Impossible Kinematics in ClusterFissioner::calculateKinematics() (C)" 
			<< Exception::eventerror;
    }
    Kinematics::twoBodyDecay(pClu2, pQ.mass(), pQ2bar.mass(), 
			     u.vect().unit(), pQ, pQ2bar);
  }
}


void ClusterFissioner::calculatePositions(const Lorentz5Momentum & pClu,
					  const LorentzPoint & positionClu,
					  const Lorentz5Momentum & pClu1,
					  const Lorentz5Momentum & pClu2,
					  LorentzPoint & positionClu1,
					  LorentzPoint & positionClu2) const {
  // Determine positions of cluster children.
  // See Marc Smith's thesis, page 127, formulas (4.122) and (4.123).
  Energy Mclu  = pClu.m();
  Energy Mclu1 = pClu1.m();
  Energy Mclu2 = pClu2.m();

  // Calculate the unit three-vector, in the C frame, along which 
  // children clusters move.
  Lorentz5Momentum u(pClu1);
  u.boost( -pClu.boostVector() );        // boost from LAB to C frame
  // the unit three-vector is then  u.vect().unit()

  Energy pstarChild = Kinematics::pstarTwoBodyDecay(Mclu,Mclu1,Mclu2);

  // First, determine the relative positions of the children clusters
  // in the parent cluster reference frame.
  Length x1 = ( 0.25*Mclu + 0.5*( pstarChild + (sqr(Mclu2) - sqr(Mclu1))/(2.0*Mclu)))/_kappa;
  Length t1 = Mclu/_kappa - x1; 
  LorentzDistance distanceClu1( x1 * u.vect().unit(), t1 );

  Length x2 = (-0.25*Mclu + 0.5*(-pstarChild + (sqr(Mclu2) - sqr(Mclu1))/(2.0*Mclu)))/_kappa;
  Length t2 = Mclu/_kappa + x2;
  LorentzDistance distanceClu2( x2 * u.vect().unit(), t2 );

  // Then, transform such relative positions from the parent cluster
  // reference frame to the Lab frame.
  distanceClu1.boost( pClu.boostVector() );
  distanceClu2.boost( pClu.boostVector() );

  // Finally, determine the absolute positions in the Lab frame.
  positionClu1 = positionClu + distanceClu1;
  positionClu2 = positionClu + distanceClu2;

}

bool ClusterFissioner::isHeavy(tcClusterPtr clu) {
  // default
  double clpow = _clPowLight;
  Energy clmax = _clMaxLight;
  // particle data for constituents
  tcPDPtr cptr[3]={tcPDPtr(),tcPDPtr(),tcPDPtr()};
  for(int ix=0;ix<min(clu->numComponents(),3);++ix) {
    cptr[ix]=clu->particle(ix)->dataPtr();
  }
  // different parameters for exotic, bottom and charm clusters
  if(CheckId::isExotic(cptr[0],cptr[1],cptr[1])) {
    clpow = _clPowExotic;
    clmax = _clMaxExotic;
  }
  else if(CheckId::hasBottom(cptr[0],cptr[1],cptr[1])) {
    clpow = _clPowBottom;
    clmax = _clMaxBottom;
  }
  else if(CheckId::hasCharm(cptr[0],cptr[1],cptr[1])) {
    clpow = _clPowCharm;
    clmax = _clMaxCharm;
  }
  bool aboveCutoff = (
		      pow(clu->mass()*UnitRemoval::InvE , clpow) 
		      > 
		      pow(clmax*UnitRemoval::InvE, clpow) 
		      + pow(clu->sumConstituentMasses()*UnitRemoval::InvE, clpow)
		      );
  // required test for SUSY clusters, since aboveCutoff alone
  // cannot guarantee (Mc > m1 + m2 + 2*m) in cut()
  static const Energy minmass 
    = getParticleData(ParticleID::d)->constituentMass();
  bool canSplitMinimally 
    = clu->mass() > clu->sumConstituentMasses() + 2.0 * minmass;
  return aboveCutoff && canSplitMinimally;
}
