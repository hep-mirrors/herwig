// -*- C++ -*-
//
// ClusterFissioner.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
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
#include "Cluster.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include <ThePEG/Utilities/DescribeClass.h>
#include "ThePEG/Interface/ParMap.h"
#include "Herwig/Utilities/AlphaS.h"

using namespace Herwig;

DescribeClass<ClusterFissioner,Interfaced>
describeClusterFissioner("Herwig::ClusterFissioner","Herwig.so");

ClusterFissioner::ClusterFissioner() :
  _clMaxLight(3.35*GeV),
  _clMaxDiquark(3.35*GeV),
  _clMaxExotic(3.35*GeV),
  _clPowLight(2.0),
  _clPowDiquark(2.0),
  _clPowExotic(2.0),
  _pSplitLight(1.0),
  _pSplitExotic(1.0),
  _phaseSpaceWeights(0),
  _dim(4),
  _fissionCluster(0),
  _kinematicThresholdChoice(0),
  _pwtDIquark(0.0),
  _diquarkClusterFission(0),
  _btClM(1.0*GeV),
  _iopRem(1),
  _kappa(1.0e15*GeV/meter),
  _enhanceSProb(0),
  _m0Fission(2.*GeV),
  _massMeasure(0),
  _probPowFactor(4.0),
  _probShift(0.0),
  _kinThresholdShift(1.0*sqr(GeV)),
  _strictDiquarkKinematics(0),
  _covariantBoost(false),
	_hadronizingStrangeDiquarks(2),
	_writeOut(0)
{
}
ClusterFissioner::~ClusterFissioner(){
}
IBPtr ClusterFissioner::clone() const {
  return new_ptr(*this);
}

IBPtr ClusterFissioner::fullclone() const {
  return new_ptr(*this);
}

void ClusterFissioner::persistentOutput(PersistentOStream & os) const {
  os << ounit(_clMaxLight,GeV) << ounit(_clMaxHeavy,GeV) << ounit(_clMaxDiquark,GeV) << ounit(_clMaxExotic,GeV)
		 << _clPowLight << _clPowHeavy << _clPowDiquark << _clPowExotic
		 << _pSplitLight << _pSplitHeavy << _pSplitExotic
     << _fissionCluster << _fissionPwt
	 << _pwtDIquark
	 << _diquarkClusterFission
     << ounit(_btClM,GeV)
     << _iopRem  << ounit(_kappa, GeV/meter)
     << _enhanceSProb << ounit(_m0Fission,GeV) << _massMeasure
		 << _dim << _phaseSpaceWeights
     << _hadronSpectrum << _kinematicThresholdChoice
     << _probPowFactor << _probShift << ounit(_kinThresholdShift,sqr(GeV))
	 << _strictDiquarkKinematics
	 << _covariantBoost
	 << _hadronizingStrangeDiquarks
	 << _writeOut
	 ;
}

void ClusterFissioner::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_clMaxLight,GeV) >> iunit(_clMaxHeavy,GeV) >> iunit(_clMaxDiquark,GeV) >> iunit(_clMaxExotic,GeV)
		 >> _clPowLight >> _clPowHeavy >> _clPowDiquark >> _clPowExotic
		 >> _pSplitLight >> _pSplitHeavy >> _pSplitExotic
     >> _fissionCluster >> _fissionPwt
	 >> _pwtDIquark
	 >> _diquarkClusterFission
   >> iunit(_btClM,GeV)
	 	 >> _iopRem >> iunit(_kappa, GeV/meter)
     >> _enhanceSProb >> iunit(_m0Fission,GeV) >> _massMeasure
		 >> _dim >> _phaseSpaceWeights
     >> _hadronSpectrum >> _kinematicThresholdChoice
     >> _probPowFactor >> _probShift >> iunit(_kinThresholdShift,sqr(GeV))
	 >> _strictDiquarkKinematics
	 >> _covariantBoost
	 >> _hadronizingStrangeDiquarks
	 >> _writeOut
	 ;
}

void ClusterFissioner::doinit() {
  Interfaced::doinit();
	if (_writeOut){
		std::ofstream out("data_CluFis.dat", std::ios::out);
		out.close();
	}
  for ( const long& id : spectrum()->heavyHadronizingQuarks() ) {
    if ( _pSplitHeavy.find(id) == _pSplitHeavy.end() ||
	 _clPowHeavy.find(id) == _clPowHeavy.end() ||
	 _clMaxHeavy.find(id) == _clMaxHeavy.end() ){
			std::cout << "id = "<<id << std::endl;
      throw InitException() << "not all parameters have been set for heavy quark cluster fission";
  }
  }
  // for default Pwts not needed to initialize
  if (_fissionCluster==0) return;  
  for ( const long& id : spectrum()->lightHadronizingQuarks() ) {
    if ( _fissionPwt.find(id) == _fissionPwt.end() )
		  // check that all relevant weights are set
      throw InitException() << "fission weights for light quarks have not been set";
  }
  double pwtDquark=_fissionPwt.find(ParticleID::d)->second;
  double pwtUquark=_fissionPwt.find(ParticleID::u)->second;
  double pwtSquark=_fissionPwt.find(ParticleID::s)->second;
  // TODO better solution for this magic  number alternative
  _fissionPwt[1103] =       _pwtDIquark * pwtDquark * pwtDquark;
  _fissionPwt[2101] = 0.5 * _pwtDIquark * pwtUquark * pwtDquark;
  _fissionPwt[2203] =       _pwtDIquark * pwtUquark * pwtUquark;
	if (_hadronizingStrangeDiquarks>0) {
		_fissionPwt[3101] = 0.5 * _pwtDIquark * pwtSquark * pwtDquark;
		_fissionPwt[3201] = 0.5 * _pwtDIquark * pwtSquark * pwtUquark;
		if (_hadronizingStrangeDiquarks==2) {
			_fissionPwt[3303] =       _pwtDIquark* pwtSquark * pwtSquark;
		}
	}
}

void ClusterFissioner::Init() {

  static ClassDocumentation<ClusterFissioner> documentation
    ("Class responsibles for chopping up the clusters");

  static Reference<ClusterFissioner,HadronSpectrum> interfaceHadronSpectrum
    ("HadronSpectrum",
     "Set the Hadron spectrum for this cluster fissioner.",
     &ClusterFissioner::_hadronSpectrum, false, false, true, false);

  // ClMax for light, Bottom, Charm and exotic (e.g. Susy) quarks
  static Parameter<ClusterFissioner,Energy>
    interfaceClMaxLight ("ClMaxLight","cluster max mass for light quarks (unit [GeV])",
                    &ClusterFissioner::_clMaxLight, GeV, 3.35*GeV, ZERO, 100.0*GeV,
		    false,false,false);  
  static Parameter<ClusterFissioner,Energy>
    interfaceClMaxDiquark ("ClMaxDiquark","cluster max mass for light hadronizing diquarks (unit [GeV])",
                    &ClusterFissioner::_clMaxDiquark, GeV, 3.35*GeV, ZERO, 100.0*GeV,
		    false,false,false);  
  static ParMap<ClusterFissioner,Energy> interfaceClMaxHeavy
    ("ClMaxHeavy",
     "ClMax for heavy quarks",
     &ClusterFissioner::_clMaxHeavy, GeV, -1, 3.35*GeV, ZERO, 100.0*GeV,
     false, false, Interface::upperlim);
  static Parameter<ClusterFissioner,Energy>
    interfaceClMaxExotic ("ClMaxExotic","cluster max mass  for exotic quarks (unit [GeV])",
                    &ClusterFissioner::_clMaxExotic, GeV, 3.35*GeV, ZERO, 100.0*GeV,
		    false,false,false);

 // ClPow for light, Bottom, Charm and exotic (e.g. Susy) quarks
 static Parameter<ClusterFissioner,double>
    interfaceClPowLight ("ClPowLight","cluster mass exponent for light quarks",
                    &ClusterFissioner::_clPowLight, 0, 2.0, 0.0, 10.0,false,false,false); 
  static ParMap<ClusterFissioner,double> interfaceClPowHeavy
    ("ClPowHeavy",
     "ClPow for heavy quarks",
     &ClusterFissioner::_clPowHeavy, -1, 1.0, 0.0, 10.0,
     false, false, Interface::upperlim); 
 static Parameter<ClusterFissioner,double>
    interfaceClPowDiquark ("ClPowDiquark","cluster mass exponent for light hadronizing diquarks",
                    &ClusterFissioner::_clPowDiquark, 0, 2.0, 0.0, 10.0,false,false,false); 
 static Parameter<ClusterFissioner,double>
    interfaceClPowExotic ("ClPowExotic","cluster mass exponent for exotic quarks",
                    &ClusterFissioner::_clPowExotic, 0, 2.0, 0.0, 10.0,false,false,false);

 // PSplit for light, Bottom, Charm and exotic (e.g. Susy) quarks
  static Parameter<ClusterFissioner,double>
    interfacePSplitLight ("PSplitLight","cluster mass splitting param for light quarks",
                    &ClusterFissioner::_pSplitLight, 0, 1.0, 0.0, 10.0,false,false,false);
  static ParMap<ClusterFissioner,double> interfacePSplitHeavy
    ("PSplitHeavy",
     "PSplit for heavy quarks",
     &ClusterFissioner::_pSplitHeavy, -1, 1.0, 0.0, 10.0,
     false, false, Interface::upperlim);
 static Parameter<ClusterFissioner,double>
    interfacePSplitExotic ("PSplitExotic","cluster mass splitting param for exotic quarks",
                    &ClusterFissioner::_pSplitExotic, 0, 1.0, 0.0, 10.0,false,false,false);


  static Switch<ClusterFissioner,int> interfaceFission
    ("Fission",
     "Option for different Fission options",
     &ClusterFissioner::_fissionCluster, 1, false, false);
  static SwitchOption interfaceFissionDefault
    (interfaceFission,
     "Default",
     "Normal cluster fission which depends on the hadron spectrum class.",
     0);
  static SwitchOption interfaceFissionNew
    (interfaceFission,
     "New",
     "Alternative cluster fission which does not depend on the hadron spectrum class",
     1);
  static SwitchOption interfaceFissionNewDiquarkSuppression
    (interfaceFission,
     "NewDiquarkSuppression",
     "Alternative cluster fission which does not depend on the hadron spectrum class"
		 " and includes a suppression of AlphaS^2(Mc) for Diquark Production during "
		 "Cluster Fission",
     -1);


  static Switch<ClusterFissioner,int> interfaceDiquarkClusterFission
    ("DiquarkClusterFission",
     "Allow clusters to fission to 1 or 2 diquark Clusters or Turn off diquark fission completely",
     &ClusterFissioner::_diquarkClusterFission, 0, false, false);
  static SwitchOption interfaceDiquarkClusterFissionAll
	(interfaceDiquarkClusterFission,
     "All",
     "Allow diquark clusters and baryon clusters to fission to new diquark Clusters",
     2);
  static SwitchOption interfaceDiquarkClusterFissionOnlyBaryonClusters
	(interfaceDiquarkClusterFission,
     "OnlyBaryonClusters",
     "Allow only baryon clusters to fission to new diquark Clusters",
     1);
  static SwitchOption interfaceDiquarkClusterFissionNo
    (interfaceDiquarkClusterFission,
     "No",
     "Don't allow clusters to fission to new diquark Clusters",
     0);
  static SwitchOption interfaceDiquarkClusterFissionOff
    (interfaceDiquarkClusterFission,
     "Off",
     "Don't allow clusters fission to draw diquarks ",
     -1);

  static ParMap<ClusterFissioner,double> interfaceFissionPwt
    ("FissionPwt",
     "The weights for quarks in the fission process.",
     &ClusterFissioner::_fissionPwt, -1, 1.0, 0.0, 10.0,
     false, false, Interface::upperlim);

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

  static Switch<ClusterFissioner,int> interfaceEnhanceSProb
    ("EnhanceSProb",
     "Option for enhancing strangeness",
     &ClusterFissioner::_enhanceSProb, 0, false, false);
  static SwitchOption interfaceEnhanceSProbNo
    (interfaceEnhanceSProb,
     "No",
     "No strangeness enhancement.",
     0);
  static SwitchOption interfaceEnhanceSProbScaled
    (interfaceEnhanceSProb,
     "Scaled",
     "Scaled strangeness enhancement",
     1);
  static SwitchOption interfaceEnhanceSProbExponential
    (interfaceEnhanceSProb,
     "Exponential",
     "Exponential strangeness enhancement",
     2);

   static Switch<ClusterFissioner,int> interfaceMassMeasure
     ("MassMeasure",
      "Option to use different mass measures",
      &ClusterFissioner::_massMeasure,0,false,false);
   static SwitchOption interfaceMassMeasureMass
     (interfaceMassMeasure,
      "Mass",
      "Mass Measure",
      0);
   static SwitchOption interfaceMassMeasureLambda
     (interfaceMassMeasure,
      "Lambda",
      "Lambda Measure",
      1);

  static Parameter<ClusterFissioner,Energy> interfaceFissionMassScale
    ("FissionMassScale",
     "Cluster fission mass scale",
     &ClusterFissioner::_m0Fission, GeV, 2.0*GeV, 0.1*GeV, 50.*GeV,
     false, false, Interface::limited);

  static Parameter<ClusterFissioner,double> interfaceProbPowFactor
     ("ProbabilityPowerFactor",
      "Power factor in ClusterFissioner bell probablity function",
      &ClusterFissioner::_probPowFactor, 2.0, 0.001, 20.0,
      false, false, Interface::limited);

  static Parameter<ClusterFissioner,double> interfaceProbShift
     ("ProbabilityShift",
      "Shifts from the center in ClausterFissioner bell probablity function",
      &ClusterFissioner::_probShift, 0.0, -10.0, 10.0,
      false, false, Interface::limited);

  static Parameter<ClusterFissioner,Energy2> interfaceKineticThresholdShift
     ("KineticThresholdShift",
      "Shifts from the kinetic threshold in ClusterFissioner",
      &ClusterFissioner::_kinThresholdShift, sqr(GeV), 0.*sqr(GeV), -10.0*sqr(GeV), 10.0*sqr(GeV),
      false, false, Interface::limited);

  static Switch<ClusterFissioner,int> interfaceKinematicThreshold
    ("KinematicThreshold",
     "Option for using static or dynamic kinematic thresholds in cluster splittings",
     &ClusterFissioner::_kinematicThresholdChoice, 0, false, false);
  static SwitchOption interfaceKinematicThresholdStatic
    (interfaceKinematicThreshold,
     "Static",
     "Set static kinematic thresholds for cluster splittings.",
     0);
  static SwitchOption interfaceKinematicThresholdDynamic
    (interfaceKinematicThreshold,
     "Dynamic",
     "Set dynamic kinematic thresholds for cluster splittings.",
     1);

  static Switch<ClusterFissioner,bool> interfaceCovariantBoost
    ("CovariantBoost",
     "Use single Covariant Boost for Cluster Fission",
     &ClusterFissioner::_covariantBoost, false, false, false);
  static SwitchOption interfaceCovariantBoostYes
    (interfaceCovariantBoost,
     "Yes",
     "Use Covariant boost",
     true);
  static SwitchOption interfaceCovariantBoostNo
    (interfaceCovariantBoost,
     "No",
     "Do NOT use Covariant boost",
     false);
  

  static Switch<ClusterFissioner,int> interfaceStrictDiquarkKinematics
    ("StrictDiquarkKinematics",
     "Option for selecting different selection criterions of diquarks for ClusterFission",
     &ClusterFissioner::_strictDiquarkKinematics, 0, false, false);
  static SwitchOption interfaceStrictDiquarkKinematicsLoose
    (interfaceStrictDiquarkKinematics,
     "Loose",
     "No kinematic threshold for diquark selection except for Mass bigger than 2 baryons",
     0);
  static SwitchOption interfaceStrictDiquarkKinematicsStrict
    (interfaceStrictDiquarkKinematics,
     "Strict",
     "Resulting clusters are at least as heavy as 2 lightest baryons",
     1);

  static Parameter<ClusterFissioner,double> interfacePwtDIquark
     ("PwtDIquark",
      "specific probability for choosing a d diquark",
      &ClusterFissioner::_pwtDIquark, 0.0, 0.0, 10.0,
      false, false, Interface::limited);

  static Switch<ClusterFissioner,int> interfacePhaseSpaceWeights
    ("PhaseSpaceWeights",
     "Include phase space weights.",
     &ClusterFissioner::_phaseSpaceWeights, 0, false, false);
  static SwitchOption interfacePhaseSpaceWeightsNo
    (interfacePhaseSpaceWeights,
     "No",
     "Do not include the effect of cluster phase space",
     0);
  static SwitchOption interfacePhaseSpaceWeightsYes
    (interfacePhaseSpaceWeights,
     "Yes",
     "Do include the effect of cluster fission phase space "
		 "related to constituent masses."
		 "Note: Need static Threshold choice",
     1);
  static SwitchOption interfacePhaseSpaceWeightsUseHadronMasses
    (interfacePhaseSpaceWeights,
     "UseHadronMasses",
     "Do include the effect of cluster fission phase space "
		 "related to hadron masses."
		 "Note: Need static Threshold choice",
     2);
  static SwitchOption interfacePhaseSpaceWeightsNoConstituentMasses
    (interfacePhaseSpaceWeights,
     "NoConstituentMasses",
     "Do not include the effect of cluster fission phase space "
		 "related to constituent masses."
		 "Note: Need static Threshold choice",
     3);

  static Parameter<ClusterFissioner,double>
    interfaceDim ("Dimension","Dimension in which phase space weights are calculated",
		  &ClusterFissioner::_dim, 0, 4.0, 0.0, 10.0,false,false,false);

	// Allowing for strange diquarks in the ClusterFission
  static Switch<ClusterFissioner,unsigned int> interfaceHadronizingStrangeDiquarks
    ("HadronizingStrangeDiquarks",
     "Option for adding strange diquarks to Cluster Fission (if Fission = New or Hybrid is enabled)",
     &ClusterFissioner::_hadronizingStrangeDiquarks, 0, false, false);
  static SwitchOption interfaceHadronizingStrangeDiquarksNo
    (interfaceHadronizingStrangeDiquarks,
     "No",
     "No strangeness containing diquarks during Cluster Fission",
     0);
  static SwitchOption interfaceHadronizingStrangeDiquarksOnlySingleStrange
    (interfaceHadronizingStrangeDiquarks,
     "OnlySingleStrange",
     "Only one strangeness containing diquarks during Cluster Fission i.e. su,sd",
     1);
  static SwitchOption interfaceHadronizingStrangeDiquarksAll
    (interfaceHadronizingStrangeDiquarks,
     "All",
     "All strangeness containing diquarks during Cluster Fission  i.e. su,sd,ss",
     2);

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

    // Sometimes the clusters decay C -> H + C' or C -> H + H' rather then C -> C' + C''
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
  Lorentz5Momentum pClu1, pClu2, pQ1, pQone, pQtwo, pQ2;
  do
    {
      succeeded = false;
      ++counter;
      // get a flavour for the qqbar pair
      drawNewFlavour(newPtr1,newPtr2,cluster);
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

      pQ1.setMass(m1);
      pQone.setMass(m);
      pQtwo.setMass(m);
      pQ2.setMass(m2);

      double weightMasses = drawNewMasses(Mc, soft1, soft2, pClu1, pClu2,
					      ptrQ1, pQ1, newPtr1, pQone,
					      newPtr2, pQtwo, ptrQ2, pQ2);
			if (weightMasses==0.0)
				continue;

      // derive the masses of the children
      Mc1 = pClu1.mass();
      Mc2 = pClu2.mass();
      // static kinematic threshold
      if(_kinematicThresholdChoice == 0) {
        if (Mc1 < m1+m || Mc2 < m+m2 || Mc1+Mc2 > Mc) continue;
				if (_phaseSpaceWeights==2 && 
						(   Mc1 < spectrum()->massLightestHadronPair(ptrQ1->dataPtr(),newPtr1->dataPtr())
						 || Mc2 < spectrum()->massLightestHadronPair(ptrQ2->dataPtr(),newPtr2->dataPtr()) ))
						continue;
      // dynamic kinematic threshold
      }
      else if(_kinematicThresholdChoice == 1) {
        bool C1 = ( sqr(Mc1) )/( sqr(m1) + sqr(m) + _kinThresholdShift ) < 1.0 ? true : false;
        bool C2 = ( sqr(Mc2) )/( sqr(m2) + sqr(m) + _kinThresholdShift ) < 1.0 ? true : false;
        bool C3 = ( sqr(Mc1) + sqr(Mc2) )/( sqr(Mc) ) > 1.0 ? true : false;

        if( C1 || C2 || C3 ) continue;
      }
      if ( _phaseSpaceWeights && phaseSpaceVeto(Mc,Mc1,Mc2,m,m1,m2, ptrQ1, ptrQ2, newPtr1, 0.0) ) {
				// reduce counter as it regards only the mass sampling
				counter--;
	  continue;
      }

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
      // override chosen masses if needed
      toHadron1 = _hadronSpectrum->chooseSingleHadron(ptrQ1->dataPtr(), newPtr1->dataPtr(),Mc1);
      if(toHadron1) { Mc1 = toHadron1->mass(); pClu1.setMass(Mc1); }
      toHadron2 = _hadronSpectrum->chooseSingleHadron(ptrQ2->dataPtr(), newPtr2->dataPtr(),Mc2);
      if(toHadron2) { Mc2 = toHadron2->mass(); pClu2.setMass(Mc2); }
      // if a beam cluster not allowed to decay to hadrons
      if(cluster->isBeamCluster() && (toHadron1||toHadron2) && softUEisOn)
	continue;
      // Check if the decay kinematics is still possible: if not then
      // force the one-hadron decay for the other cluster as well.
      if(Mc1 + Mc2  >  Mc) {
	if(!toHadron1) {
	  toHadron1 = _hadronSpectrum->chooseSingleHadron(ptrQ1->dataPtr(), newPtr1->dataPtr(),Mc-Mc2);
	  if(toHadron1) { Mc1 = toHadron1->mass(); pClu1.setMass(Mc1); }
	}
	else if(!toHadron2) {
	  toHadron2 = _hadronSpectrum->chooseSingleHadron(ptrQ2->dataPtr(), newPtr2->dataPtr(),Mc-Mc1);
	  if(toHadron2) { Mc2 = toHadron2->mass(); pClu2.setMass(Mc2); }
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
  calculateKinematics(pClu,p0Q1,toHadron1,toHadron2,
		      pClu1,pClu2,pQ1,pQone,pQtwo,pQ2);

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
  Lorentz5Momentum pClu1, pClu2, pQ1, pQone, pQtwo, pQ2;
  do {
    succeeded = false;
    ++counter;

    // get a flavour for the qqbar pair
    drawNewFlavour(newPtr1,newPtr2,cluster);
    
    // randomly pick which will be (anti)diquark and which a mesonic cluster
    if(UseRandom::rndbool()) {
      swap(iq1,iq2);
      swap(rem1,rem2);
    }
    // check first order
    if(cantMakeHadron(ptrQ[iq1], newPtr1) || !spectrum()->canMakeDiQuark(ptrQ[iq2], newPtr2)) {
      swap(newPtr1,newPtr2);
    }
    // check again
    if(cantMakeHadron(ptrQ[iq1], newPtr1) || !spectrum()->canMakeDiQuark(ptrQ[iq2], newPtr2)) {
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

    pQ1.setMass(m1);
    pQone.setMass(m);
    pQtwo.setMass(m);
    pQ2.setMass(m2);

    double weightMasses = drawNewMasses(mmax, soft1, soft2, pClu1, pClu2,
					    ptrQ[iq1], pQ1, newPtr1, pQone,
					    newPtr2, pQtwo, ptrQ[iq1], pQ2);

		if (weightMasses == 0.0) continue;

    Mc1 = pClu1.mass();
		Mc2 = pClu2.mass();

    if(Mc1 < m1+m || Mc2 < m+m2 || Mc1+Mc2 > mmax) continue;

    if ( _phaseSpaceWeights && phaseSpaceVeto(mmax,Mc1,Mc2,m,m1,m2) ) {
				// reduce counter as it regards only the mass sampling
				counter--;
	continue;
    }
    
    // check if need to force meson clster to hadron
    toHadron = _hadronSpectrum->chooseSingleHadron(ptrQ[iq1]->dataPtr(), newPtr1->dataPtr(),Mc1);
    if(toHadron) { Mc1 = toHadron->mass(); pClu1.setMass(Mc1); }
    // check if need to force diquark cluster to be on-shell
    toDiQuark = false;
    diquark = spectrum()->makeDiquark(ptrQ[iq2]->dataPtr(), newPtr2->dataPtr());
    if(Mc2 < diquark->constituentMass()) {
      Mc2 = diquark->constituentMass(); pClu2.setMass(Mc2);
      toDiQuark = true;
    }
    // if a beam cluster not allowed to decay to hadrons
    if(cluster->isBeamCluster() && toHadron && softUEisOn)
      continue;
    // Check if the decay kinematics is still possible: if not then
    // force the one-hadron decay for the other cluster as well.
    if(Mc1 + Mc2  >  mmax) {
      if(!toHadron) {
	toHadron = _hadronSpectrum->chooseSingleHadron(ptrQ[iq1]->dataPtr(), newPtr1->dataPtr(),mmax-Mc2);
	if(toHadron) { Mc1 = toHadron->mass(); pClu1.setMass(Mc1); }
      }
      else if(!toDiQuark) {
	Mc2 = _hadronSpectrum->massLightestHadron(ptrQ[iq2]->dataPtr(), newPtr2->dataPtr()); pClu2.setMass(Mc2);
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
    rval.first = (_hadronSpectrum->lightestHadron(hadron,newPtr->dataPtr()))->produceParticle();
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

/**
 * Calculate the phase space weight for  M1*M2*(2 body PhaseSpace) ignore constituent masses
 */
double ClusterFissioner::weightFlatPhaseSpaceNoConstituentMasses(const Energy Mc, const Energy Mc1, const Energy Mc2) const {
	double M_temp = Mc/GeV;
	double M1_temp = Mc1/GeV;
	double M2_temp = Mc2/GeV;
	if (sqr(M_temp)<sqr(M1_temp+M2_temp)) {
		// This should be checked before
		throw Exception()
			<< "ERROR in ClusterFissioner::weightFlatPhaseSpaceNoConstituentMasses\n"
			<< "ClusterFissioner has not checked Masses properly\n"
			<< "Mc  = " << M_temp << "\n"
			<< "Mc1 = " << M1_temp << "\n"
			<< "Mc2 = " << M2_temp << "\n"
			<< Exception::warning;
		return 0.0;
	}
	double lam = Kinematics::kaellen(M_temp,  M1_temp, M2_temp);
	double ratio;
	//  new weight with the Jacobi factor M1*M2 of the Mass integration
	double PSweight = M1_temp*M2_temp*pow(sqrt(lam),_dim-3.);
	// overestimate only possible for dim>=3.0
	assert(_dim>=3.0);
	//  new improved overestimate with the Jacobi factor M1*M2 of the Mass integration
	double overEstimate = pow(6.0*sqrt(3.0), 3.0 - _dim)*pow(M_temp, 2.*(_dim-2.));
	ratio = PSweight/overEstimate;
	if (!(ratio>=0) || !(ratio<=1)) {
		throw Exception()
			<< "ERROR in ClusterFissioner::weightFlatPhaseSpaceNoConstituentMasses\n"
			<< "ratio = " <<ratio
			<<" M "<<M_temp
			<<" M1 "<<M1_temp
			<<" M2 "<<M2_temp <<"\t"<<_dim<<"\t" << lam 
			<<"\t"<< overEstimate<<"\n\n"
			<< Exception::runerror;
	}
	return ratio;
}
/**
 * Calculate the phase space weight for  M1*M2*(2 body PhaseSpace)^3
 */
double ClusterFissioner::weightPhaseSpaceConstituentMasses(const Energy Mc, const Energy Mc1, const Energy Mc2,
		const Energy m, const Energy m1, const Energy m2, const double power) const {
	double M_temp = Mc/GeV;
	double M1_temp = Mc1/GeV;
	double M2_temp = Mc2/GeV;
	double m_temp = m/GeV;
	double m1_temp = m1/GeV;
	double m2_temp = m2/GeV;
	if (sqr(M_temp)<sqr(M1_temp+M2_temp)
			|| sqr(M1_temp)<sqr(m1_temp+m_temp)
			|| sqr(M2_temp)<sqr(m2_temp+m_temp)
			) {
		// This should be checked before
		throw Exception()
			<< "ERROR in ClusterFissioner::weightPhaseSpaceConstituentMasses\n"
			<< "ClusterFissioner has not checked Masses properly\n"
			<< "Mc  = " << M_temp << "\n"
			<< "Mc1 = " << M1_temp << "\n"
			<< "Mc2 = " << M2_temp << "\n"
			<< "m1  = " << m1_temp << "\n"
			<< "m2  = " << m2_temp << "\n"
			<< "m   = " << m_temp << "\n"
			<< Exception::warning;
		return 0.0;
	}
	double lam1 = Kinematics::kaellen(M1_temp, m1_temp, m_temp);
	double lam2 = Kinematics::kaellen(M2_temp, m2_temp, m_temp);
	double lam3 = Kinematics::kaellen(M_temp,  M1_temp, M2_temp);
	double ratio;
	//  new weight with the Jacobi factor M1*M2 of the Mass integration
	double PSweight = pow(lam1*lam2*lam3,(_dim-3.)/2.0)*pow(M1_temp*M2_temp,3.-_dim);
	// overestimate only possible for dim>=3.0
	assert(_dim>=3.0);
	//  new improved overestimate with the Jacobi factor M1*M2 of the Mass integration
	double overEstimate = pow(6.0*sqrt(3.0), 3.0 - _dim)*pow(M_temp, 4.*_dim-12.);
	ratio = PSweight/overEstimate;
	if (!(ratio>=0)) std::cout << "ratio = " <<ratio<<" M "<<M_temp<<" M1 "<<M1_temp<<" M2 "<<M2_temp<<" m1 "<<m1_temp<<" m2 "<<m2_temp<<" m "<<m_temp<<"\t"<<_dim<<"\t" << lam1<<"\t"<< lam2<<"\t" << lam3 <<"\t"<< overEstimate<<"\n\n";
	if (!(ratio>=0) || !(ratio<=1)) {
		throw Exception()
			<< "ERROR in ClusterFissioner::weightPhaseSpaceConstituentMasses\n"
			<< "ratio = " <<ratio
			<<" M "<<M_temp
			<<" M1 "<<M1_temp
			<<" M2 "<<M2_temp
			<<" m1 "<<m1_temp
			<<" m2 "<<m2_temp
			<<" m "<<m_temp <<"\t"<<_dim
			<<"\t" << lam1<<"\t"<< lam2<<"\t" << lam3
			<<"\t"<< overEstimate<<"\n\n"
			<< Exception::runerror;
	}
	// multiply by overestimate of power of matrix element to modulate the phase space with (M1*M2)^power
	if (power) {
		double powerLawOver = power<0 ? pow(Mc1*Mc2/((m1+m)*(m2+m)),power):pow(Mc1*Mc2/((Mc-(m1+m))*(Mc-(m2+m))),power);
		ratio*=powerLawOver;
	}
	return ratio;
}
/**
 * Calculate the phase space weight for  M1*M2*(2 body PhaseSpace)^3
 * using Hadron Masses
 */
double ClusterFissioner::weightFlatPhaseSpaceHadronMasses(const Energy Mc, const Energy Mc1, const Energy Mc2, tcPPtr pQ, tcPPtr pQ1, tcPPtr pQ2) const {
	auto LHP1 = spectrum()->lightestHadronPair(pQ1->dataPtr(),pQ->dataPtr());
	auto LHP2 = spectrum()->lightestHadronPair(pQ2->dataPtr(),pQ->dataPtr());
	if (sqr(Mc1)<sqr(LHP1.first->mass()+LHP1.second->mass()))
		return true;
	if (sqr(Mc2)<sqr(LHP2.first->mass()+LHP2.second->mass()))
		return true;
	double lam1 = sqrt(Kinematics::kaellen(Mc1/GeV, LHP1.first->mass()/GeV, LHP1.second->mass()/GeV));
	double lam2 = sqrt(Kinematics::kaellen(Mc2/GeV, LHP2.first->mass()/GeV, LHP2.second->mass()/GeV));
	double lam3 = sqrt(Kinematics::kaellen(Mc/GeV,  Mc1/GeV, Mc2/GeV));
	double ratio;
	//  new weight with the Jacobi factor M1*M2 of the Mass integration
	double PSweight = pow(lam1*lam2*lam3,_dim-3.)*pow(Mc1*Mc2/GeV2,3.-_dim);
	// overestimate only possible for dim>=3.0
	assert(_dim>=3.0);
	//  new improved overestimate with the Jacobi factor M1*M2 of the Mass integration
	double overEstimate = pow(6.0*sqrt(3.0), 3.0 - _dim)*pow(Mc/GeV, 4.*_dim-12.);
	ratio = PSweight/overEstimate;
	if (!(ratio>=0) || !(ratio<=1)) {
		throw Exception()
			<< "ERROR in ClusterFissioner::weightFlatPhaseSpaceHadronMasses\n"
			<< "ratio = " <<ratio
			<<" M "<<Mc/GeV
			<<" M1 "<<Mc1/GeV
			<<" M2 "<<Mc2/GeV <<"\t"<<_dim<<"\t" << lam1<<"\t" << lam2  <<"\t" << lam3 
			<<"\t"<< overEstimate<<"\n\n"
			<< Exception::runerror;
	}
	return ratio;
}

/**
 * Veto for the phase space weight
 * returns true if proposed Masses are rejected
 * 					else returns false
 */
bool ClusterFissioner::phaseSpaceVeto(const Energy Mc, const Energy Mc1, const Energy Mc2,
		const Energy m, const Energy m1, const Energy m2, tcPPtr pQ1, tcPPtr pQ2, tcPPtr pQ, const double power) const {
	switch (_phaseSpaceWeights)
	{
		case 1:
			return phaseSpaceVetoConstituentMasses(Mc, Mc1, Mc2, m, m1, m2, power);
		case 2:
			return phaseSpaceVetoHadronPairs(Mc, Mc1, Mc2, pQ, pQ1, pQ2);
		case 3:
			return phaseSpaceVetoNoConstituentMasses(Mc, Mc1, Mc2);
		default:
			assert(false);
	}
}

/**
 * Veto for the phase space weight
 * returns true if proposed Masses are rejected
 * 					else returns false
 */
bool ClusterFissioner::phaseSpaceVetoConstituentMasses(const Energy Mc, const Energy Mc1, const Energy Mc2,
		const Energy m, const Energy m1, const Energy m2, const double power) const {
	return (UseRandom::rnd()>weightPhaseSpaceConstituentMasses(Mc, Mc1, Mc2, m, m1, m2, power));
}
bool ClusterFissioner::phaseSpaceVetoNoConstituentMasses(const Energy Mc, const Energy Mc1, const Energy Mc2) const {
	return (UseRandom::rnd()>weightFlatPhaseSpaceNoConstituentMasses(Mc, Mc1, Mc2));
}

bool ClusterFissioner::phaseSpaceVetoHadronPairs(const Energy Mc, const Energy Mc1, const Energy Mc2, tcPPtr pQ, tcPPtr pQ1, tcPPtr pQ2) const {
	return (UseRandom::rnd()>weightFlatPhaseSpaceHadronMasses(Mc, Mc1, Mc2, pQ, pQ1, pQ2));
}

/**
 * Calculate the masses and possibly kinematics of the cluster
 * fission at hand; if calculateKineamtics is perfomring non-trivial
 * steps kinematics claulcated here will be overriden. Currentl;y resorts to the default
 */
double ClusterFissioner::drawNewMasses(const Energy Mc, const bool soft1, const bool soft2,
		Lorentz5Momentum& pClu1, Lorentz5Momentum& pClu2,
		tcPPtr ptrQ1, const Lorentz5Momentum& pQ1,
		tcPPtr, const Lorentz5Momentum& pQone,
		tcPPtr, const Lorentz5Momentum& pQtwo,
		tcPPtr ptrQ2,  const Lorentz5Momentum& pQ2) const {
	// power for splitting
	double exp1 = !spectrum()->isExotic(ptrQ1->dataPtr()) ? _pSplitLight : _pSplitExotic;
	double exp2 = !spectrum()->isExotic(ptrQ2->dataPtr()) ? _pSplitLight : _pSplitExotic;
	for ( const long& id : spectrum()->heavyHadronizingQuarks() ) {
		assert(_pSplitHeavy.find(id) != _pSplitHeavy.end());
		if ( spectrum()->hasHeavy(id,ptrQ1->dataPtr()) ) exp1 = _pSplitHeavy.find(id)->second;
		if ( spectrum()->hasHeavy(id,ptrQ2->dataPtr()) ) exp2 = _pSplitHeavy.find(id)->second;
	}

	Energy M1 = drawChildMass(Mc,pQ1.mass(),pQ2.mass(),pQone.mass(),exp1,soft1);
	Energy M2 = drawChildMass(Mc,pQ2.mass(),pQ1.mass(),pQtwo.mass(),exp2,soft2);

	pClu1.setMass(M1);
	pClu2.setMass(M2);

	return 1.0; // succeeds
}


void ClusterFissioner::drawNewFlavourDiquarks(PPtr& newPtrPos,PPtr& newPtrNeg,
		const ClusterPtr & clu) const {

	// Flavour is assumed to be only  u, d, s,  with weights
	// (which are not normalized probabilities) given
	// by the same weights as used in HadronsSelector for
	// the decay of clusters into two hadrons.

	unsigned hasDiquarks=0;
	assert(clu->numComponents()==2);
	tcPDPtr pD1=clu->particle(0)->dataPtr();
	tcPDPtr pD2=clu->particle(1)->dataPtr();
	bool isDiq1=DiquarkMatcher::Check(pD1->id());
	if (isDiq1) hasDiquarks++;
	bool isDiq2=DiquarkMatcher::Check(pD2->id());
	if (isDiq2) hasDiquarks++;
	assert(hasDiquarks<=2);
	Energy Mc=(clu->momentum().mass());
  Energy minMass;
  double weight;
  Selector<long> choice;
  // adding quark-antiquark pairs to the selection list
  for ( const long& id : spectrum()->lightHadronizingQuarks() ) {
		minMass=spectrum()->massLightestHadronPair(pD1,pD2);
	  if (_fissionCluster==0) choice.insert(_hadronSpectrum->pwtQuark(id),id);
		else if (abs(_fissionCluster)==1) choice.insert(_fissionPwt.find(id)->second,id);
	  else assert(false);
  }
  // adding diquark-antidiquark pairs to the selection list
  switch (hasDiquarks)
  {
  	case 0:
		for ( const long& id : spectrum()->lightHadronizingDiquarks() ) {
			if (_strictDiquarkKinematics) {
				tPDPtr cand = getParticleData(id);
				Energy mH1=spectrum()->massLightestHadron(pD2,cand);
				Energy mH2=spectrum()->massLightestHadron(cand,pD1);
				minMass = mH1 + mH2;
			}
			else {
				minMass = spectrum()->massLightestBaryonPair(pD1,pD2);
			}
			if (Mc < minMass) continue;
			if (_fissionCluster==0) weight = _hadronSpectrum->pwtQuark(id);
			else if (abs(_fissionCluster)==1) weight = _fissionPwt.find(id)->second;
			else assert(false); 
			if (_fissionCluster==-1)
				weight*=sqr(Herwig::Math::alphaS(Mc, 0.25*GeV,3, 2));
			choice.insert(weight,id);
		}
  		break;
  	case 1:
		if (_diquarkClusterFission<1) break;
		for ( const long& id : spectrum()->lightHadronizingDiquarks() ) {
			tPDPtr diq = getParticleData(id);
			if (isDiq1)
				minMass = spectrum()->massLightestHadron(pD2,diq)
						+ spectrum()->massLightestBaryonPair(diq,pD1);
			else
				minMass = spectrum()->massLightestHadron(pD1,diq)
						+ spectrum()->massLightestBaryonPair(diq,pD2);
			if (Mc < minMass) continue;
			if (_fissionCluster==0) weight = _hadronSpectrum->pwtQuark(id);
			else if (abs(_fissionCluster)==1) weight = _fissionPwt.find(id)->second;
			else assert(false); 
			if (_fissionCluster==-1)
				weight*=sqr(Herwig::Math::alphaS(Mc, 0.25*GeV,3, 2));
			choice.insert(weight,id);
		}
		break;
  	case 2:
		if (_diquarkClusterFission<2) break;
		for ( const long& id : spectrum()->lightHadronizingDiquarks() ) {
			tPDPtr diq = getParticleData(id);
			if (Mc < spectrum()->massLightestBaryonPair(pD1,pD2)) {
				throw Exception() << "Found Diquark Cluster:\n" << *clu << "\nwith  MassCluster = "
					<< Mc/GeV <<" GeV MassLightestBaryonPair = "
					<< spectrum()->massLightestBaryonPair(pD1,pD2)/GeV
					<< " GeV cannot decay" << Exception::eventerror;
			}
			minMass = spectrum()->massLightestBaryonPair(pD1,diq)
					+ spectrum()->massLightestBaryonPair(diq,pD2);
			if (Mc < minMass) continue;
			if (_fissionCluster==0) weight = _hadronSpectrum->pwtQuark(id);
			else if (abs(_fissionCluster)==1) weight = _fissionPwt.find(id)->second;
			else assert(false); 
			if (_fissionCluster==-1)
				weight*=sqr(Herwig::Math::alphaS(Mc, 0.25*GeV,3, 2));
			choice.insert(weight,id);
		}
  		break;
  	default:
  		assert(false);
  }
  assert(choice.size()>0);
  long idNew = choice.select(UseRandom::rnd());
  newPtrPos = getParticle(idNew);
  newPtrNeg = getParticle(-idNew);
  assert(newPtrPos);
  assert(newPtrNeg);
  assert(newPtrPos->dataPtr());
  assert(newPtrNeg->dataPtr());

}

void ClusterFissioner::drawNewFlavourQuarks(PPtr& newPtrPos,PPtr& newPtrNeg) const {

  // Flavour is assumed to be only  u, d, s,  with weights
  // (which are not normalized probabilities) given
  // by the same weights as used in HadronsSelector for
  // the decay of clusters into two hadrons.


  Selector<long> choice;
  switch(abs(_fissionCluster)){
  case 0:
    for ( const long& id : spectrum()->lightHadronizingQuarks() )
      choice.insert(_hadronSpectrum->pwtQuark(id),id);
    break;
  case 1:
    for ( const long& id : spectrum()->lightHadronizingQuarks() )
      choice.insert(_fissionPwt.find(id)->second,id);
    break;
  default :
    assert(false);
  }
  long idNew = choice.select(UseRandom::rnd());
  newPtrPos = getParticle(idNew);
  newPtrNeg = getParticle(-idNew);
  assert (newPtrPos);
  assert(newPtrNeg);
  assert (newPtrPos->dataPtr());
  assert(newPtrNeg->dataPtr());

}

void ClusterFissioner::drawNewFlavourEnhanced(PPtr& newPtrPos,PPtr& newPtrNeg,
                                              Energy2 mass2) const {

  if ( spectrum()->gluonId() != ParticleID::g )
    throw Exception() << "strange enhancement only working with Standard Model hadronization"
		      << Exception::runerror;
  
  // Flavour is assumed to be only  u, d, s,  with weights
  // (which are not normalized probabilities) given
  // by the same weights as used in HadronsSelector for
  // the decay of clusters into two hadrons.

    double prob_d = 0.;
    double prob_u = 0.;
    double prob_s = 0.;
    double scale = abs(double(sqr(_m0Fission)/mass2));
    // Choose which splitting weights you wish to use
switch(abs(_fissionCluster)){
  // 0: ClusterFissioner and ClusterDecayer use the same weights
  case 0:
    prob_d = _hadronSpectrum->pwtQuark(ParticleID::d);
     prob_u = _hadronSpectrum->pwtQuark(ParticleID::u);
     /* Strangeness enhancement:
        Case 1: probability scaling
        Case 2: Exponential scaling
     */
     if (_enhanceSProb == 1)
        prob_s = (_maxScale < scale) ? 0. : pow(_hadronSpectrum->pwtQuark(ParticleID::s),scale);
     else if (_enhanceSProb == 2)
        prob_s = (_maxScale < scale) ? 0. : exp(-scale);
    break;
    /* 1: ClusterFissioner uses its own unique set of weights,
       i.e. decoupled from ClusterDecayer */
  case 1:
    prob_d = _fissionPwt.find(ParticleID::d)->second;
    prob_u = _fissionPwt.find(ParticleID::u)->second;
     if (_enhanceSProb == 1)
       prob_s = (_maxScale < scale) ? 0. : pow(_fissionPwt.find(ParticleID::s)->second,scale);
     else if (_enhanceSProb == 2)
        prob_s = (_maxScale < scale) ? 0. : exp(-scale);
    break;
  default:
    assert(false);
  }

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


Energy2 ClusterFissioner::clustermass(const ClusterPtr & cluster) const {
  Lorentz5Momentum pIn = cluster->momentum();
  Energy2 endpointmass2 = sqr(cluster->particle(0)->mass() +
  cluster->particle(1)->mass());
  Energy2 singletm2 = pIn.m2();
  // Return either the cluster mass, or the lambda measure
  return (_massMeasure == 0) ? singletm2 : singletm2 - endpointmass2;
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

  Energy2 mag2 = u.vect().mag2();
  InvEnergy fact = mag2>ZERO ? 1./sqrt(mag2) : 1./GeV;

  Length x1 = ( 0.25*Mclu + 0.5*( pstarChild + (sqr(Mclu2) - sqr(Mclu1))/(2.0*Mclu)))/_kappa;
  Length t1 = Mclu/_kappa - x1;
  LorentzDistance distanceClu1( x1 * fact * u.vect(), t1 );

  Length x2 = (-0.25*Mclu + 0.5*(-pstarChild + (sqr(Mclu2) - sqr(Mclu1))/(2.0*Mclu)))/_kappa;
  Length t2 = Mclu/_kappa + x2;
  LorentzDistance distanceClu2( x2 * fact * u.vect(), t2 );

  // Then, transform such relative positions from the parent cluster
  // reference frame to the Lab frame.
  distanceClu1.boost( pClu.boostVector() );
  distanceClu2.boost( pClu.boostVector() );

  // Finally, determine the absolute positions in the Lab frame.
  positionClu1 = positionClu + distanceClu1;
  positionClu2 = positionClu + distanceClu2;

}

bool ClusterFissioner::ProbabilityFunction(double scale, double threshold) {
  double cut = UseRandom::rnd(0.0,1.0);
  return 1./(1.+pow(abs((threshold-_probShift)/scale),_probPowFactor)) > cut ? true : false;
}
bool ClusterFissioner::ProbabilityFunctionPower(double Mass, double threshold) {
  double cut = UseRandom::rnd(0.0,1.0);
	if ((Mass-threshold)<=0)
		return false;
	return 1.0/(1.0 + _probPowFactor*pow(1.0/(Mass-threshold),_clPowLight)) > cut ? true : false;
}


bool ClusterFissioner::isHeavy(tcClusterPtr clu) {
  // particle data for constituents
  tcPDPtr cptr[3]={tcPDPtr(),tcPDPtr(),tcPDPtr()};
	bool hasDiquark=0;
  for(size_t ix=0;ix<min(clu->numComponents(),3);++ix) {
    cptr[ix]=clu->particle(ix)->dataPtr();
		// Assuming diquark masses are ordered with larger id corresponding to larger masses
		if (DiquarkMatcher::Check(*(cptr[ix]))) {
			hasDiquark=true;
			break;
		}
  }
  // different parameters for exotic, bottom and charm clusters
  double clpow = !spectrum()->isExotic(cptr[0],cptr[1],cptr[1]) ? _clPowLight : _clPowExotic;
  Energy clmax = !spectrum()->isExotic(cptr[0],cptr[1],cptr[1]) ? _clMaxLight : _clMaxExotic;
	// if no heavy quark is found in the cluster, but diquarks are present use 
	// different ClMax and ClPow
	if ( hasDiquark) {
		clpow = _clPowDiquark;
		clmax = _clMaxDiquark;
	}

  for ( const long& id : spectrum()->heavyHadronizingQuarks() ) {
		if ( spectrum()->hasHeavy(id,cptr[0],cptr[1],cptr[1])) {
      clpow = _clPowHeavy[id];
      clmax = _clMaxHeavy[id];
    }
  }
   // required test for SUSY clusters, since aboveCutoff alone
  // cannot guarantee (Mc > m1 + m2 + 2*m) in cut()
  static const Energy minmass
    = getParticleData(ParticleID::d)->constituentMass();
  bool aboveCutoff = false, canSplitMinimally = false;
  // static kinematic threshold
  if(_kinematicThresholdChoice == 0) {
    aboveCutoff = (
    	      pow(clu->mass()*UnitRemoval::InvE , clpow)
    	      >
    	      pow(clmax*UnitRemoval::InvE, clpow)
    	      + pow(clu->sumConstituentMasses()*UnitRemoval::InvE, clpow)
    	      );

    canSplitMinimally = clu->mass() > clu->sumConstituentMasses() + 2.0 * minmass;
  }
  // dynamic kinematic threshold
  else if(_kinematicThresholdChoice == 1) {
    //some smooth probablity function to create a dynamic thershold
    double scale     = pow(clu->mass()/GeV , clpow);
    double threshold = pow(clmax/GeV, clpow)
                     + pow(clu->sumConstituentMasses()/GeV, clpow);
		aboveCutoff = ProbabilityFunction(scale,threshold);

    scale     = clu->mass()/GeV;
    threshold = clu->sumConstituentMasses()/GeV + 2.0 * minmass/GeV;

		canSplitMinimally = ProbabilityFunction(scale,threshold);
	}
	// probablistic kinematic threshold
	else if(_kinematicThresholdChoice == 2) {
		// Consistent power law for CF probability
		double Mass     = clu->mass()/GeV;
		double threshold = clu->sumConstituentMasses()/GeV + 2.0 * minmass/GeV;
		aboveCutoff = ProbabilityFunctionPower(Mass,threshold + clmax/GeV);

		canSplitMinimally = Mass - threshold>ZERO;
  }

  return aboveCutoff && canSplitMinimally;
}
