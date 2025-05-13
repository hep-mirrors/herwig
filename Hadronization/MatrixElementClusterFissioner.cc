// -*- C++ -*-
//
// MatrixElementClusterFissioner.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// Thisk is the implementation of the non-inlined, non-templated member
// functions of the MatrixElementClusterFissioner class.
//

#include "MatrixElementClusterFissioner.h"
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

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <cassert>
#include <vector>
using namespace Herwig;

DescribeClass<MatrixElementClusterFissioner,ClusterFissioner>
describeMatrixElementClusterFissioner("Herwig::MatrixElementClusterFissioner","Herwig.so");

// Initialize counters
unsigned int MatrixElementClusterFissioner::_counterMaxLoopViolations=0;
unsigned int MatrixElementClusterFissioner::_counterPaccGreater1=0;
unsigned int MatrixElementClusterFissioner::_counterFissionMatrixElement=0;
MatrixElementClusterFissioner::MatrixElementClusterFissioner() :
  _allowHadronFinalStates(2),
  _massSampler(0),
  _phaseSpaceSamplerCluster(0),
  _phaseSpaceSamplerConstituent(0),
  _matrixElement(0),
	_fissionApproach(0),
	_powerLawPower(0.0),
	_maxLoopFissionMatrixElement(5000000),
	_safetyFactorMatrixElement(10.0),
	_epsilonResolution(1.0),
	_failModeMaxLoopMatrixElement(0),
	_writeOut(0)
{}

MatrixElementClusterFissioner::~MatrixElementClusterFissioner(){
	if (_counterMaxLoopViolations){
		std::cout << "WARNING: There have been "
			<< _counterMaxLoopViolations
			<< "/" << _counterFissionMatrixElement
			<< " Maximum tries violations during the Simulation! ("
			<< 100.0*double(_counterMaxLoopViolations)/double(_counterFissionMatrixElement)
			<< " %)\n"
			<< "You may want to reduce MatrixElementClusterFissioner:SafetyFactorOverEst (=>potentially Pacc>=1) "
			<< "and/or increase MatrixElementClusterFissioner:MaxLoopMatrixElement (=>slower performance)\n";
	}
	if (_counterPaccGreater1){
		std::cout << "WARNING: There have been "
			<< _counterPaccGreater1
			<< "/" << _counterFissionMatrixElement
			<< " Pacc>=1 exceptions during the Simulation! ("
			<< 100.0*double(_counterPaccGreater1)/double(_counterFissionMatrixElement)
			<< " %)\n"
			<< "You may want to increase MatrixElementClusterFissioner:SafetyFactorOverEst to reduce this "
			<< "\nNOTE: look at the log file and look for the larges Pacc accepted and multiply "
			<< "\n      MatrixElementClusterFissioner:SafetyFactorOverEst by this factor";
	}
}
IBPtr MatrixElementClusterFissioner::clone() const {
  return new_ptr(*this);
}

IBPtr MatrixElementClusterFissioner::fullclone() const {
  return new_ptr(*this);
}

void MatrixElementClusterFissioner::persistentOutput(PersistentOStream & os) const {
  os << _allowHadronFinalStates
	 << _massSampler
	 << _phaseSpaceSamplerCluster
	 << _phaseSpaceSamplerConstituent
	 << _matrixElement
	 << _fissionApproach
	 << _powerLawPower
	 << _maxLoopFissionMatrixElement
	 << _safetyFactorMatrixElement
	 << _epsilonResolution
	 << _failModeMaxLoopMatrixElement
	 << _writeOut
	 ;
}
void MatrixElementClusterFissioner::persistentInput(PersistentIStream & is, int) {
  is >> _allowHadronFinalStates
	 >> _massSampler
	 >> _phaseSpaceSamplerCluster
	 >> _phaseSpaceSamplerConstituent
	 >> _matrixElement
	 >> _fissionApproach
	 >> _powerLawPower
	 >> _maxLoopFissionMatrixElement
	 >> _safetyFactorMatrixElement
	 >> _epsilonResolution
	 >> _failModeMaxLoopMatrixElement
	 >> _writeOut
	 ;
}
/*
namespace{
	void printV(Lorentz5Momentum p) {
		std::cout << "("<<p.e()/GeV<<"|"<<p.vect().x()/GeV<<","<<p.vect().y()/GeV<<","<<p.vect().z()/GeV<<") Mass = "<<p.mass()/GeV<<" m = "<<p.m()/GeV<<"\n";
	}
}
*/
void MatrixElementClusterFissioner::doinit() {
	ClusterFissioner::doinit();
	// TODO: Some User warnings/errors but not complete list
	if (_matrixElement!=0 && _fissionApproach==0) generator()->logWarning( Exception(
				"For non-trivial MatrixElement you need to enable FissionApproach=New or Hybrid\n",
				Exception::warning));
}

void MatrixElementClusterFissioner::Init() {
	// TODO test for what can be copied
	// std::cout << "sizeof Cluster "  << sizeof(Cluster) << std::endl;
	// std::cout << "sizeof Axis "  << sizeof(Axis) << std::endl;
	// std::cout << "sizeof Axis & "  << sizeof(Axis&) << std::endl;
	// std::cout << "sizeof Axis * "  << sizeof(Axis*) << std::endl;
	// std::cout << "sizeof energy "  << sizeof(Energy) << std::endl;
	// std::cout << "sizeof energy2 "  << sizeof(Energy2) << std::endl;
	// std::cout << "sizeof bool "  << sizeof(bool) << std::endl;
	// std::cout << "sizeof energy &  "  << sizeof(Energy&) << std::endl;
	// std::cout << "sizeof energy *  "  << sizeof(Energy*) << std::endl;

  static ClassDocumentation<MatrixElementClusterFissioner> documentation
    ("Class responsibles for chopping up the clusters");

  static Switch<MatrixElementClusterFissioner,int> interfaceFissionApproach
    ("FissionApproach",
     "Option for different Cluster Fission approaches",
     &MatrixElementClusterFissioner::_fissionApproach, 0, false, false);
  static SwitchOption interfaceFissionApproachDefault
    (interfaceFissionApproach,
     "Default",
     "Default Herwig-7.3.0 cluster fission without restructuring",
     0);
  static SwitchOption interfaceFissionApproachNew
    (interfaceFissionApproach,
     "New",
     "New cluster fission which allows to choose MassSampler"
		 ", PhaseSpaceSampler and MatrixElement",
     1);
  static SwitchOption interfaceFissionApproachHybrid
    (interfaceFissionApproach,
     "Hybrid",
     "New cluster fission which allows to choose MassSampler"
		 ", PhaseSpaceSampler and MatrixElement, but uses Default"
		 " Approach for BeamClusters",
     2);
  static SwitchOption interfaceFissionApproachPreservePert
		(interfaceFissionApproach,
     "PreservePert",
     "New cluster fission which allows to choose MassSampler"
		 ", PhaseSpaceSampler and MatrixElement, but uses Default"
		 " Approach for BeamClusters",
     3);
  static SwitchOption interfaceFissionApproachPreserveNonPert
		(interfaceFissionApproach,
     "PreserveNonPert",
     "New cluster fission which allows to choose MassSampler"
		 ", PhaseSpaceSampler and MatrixElement, but uses Default"
		 " Approach for BeamClusters",
     4);
  static SwitchOption interfaceFissionApproachPreserveFirstPert
		(interfaceFissionApproach,
     "PreserveFirstPert",
     "New cluster fission which allows to choose MassSampler"
		 ", PhaseSpaceSampler and MatrixElement, but uses Default"
		 " Approach for BeamClusters",
     5);

  // Switch C->H1,C2 C->H1,H2 on or off
  static Switch<MatrixElementClusterFissioner,int> interfaceAllowHadronFinalStates
    ("AllowHadronFinalStates",
     "Option for allowing hadron final states of cluster fission",
     &MatrixElementClusterFissioner::_allowHadronFinalStates, 0, false, false);
  static SwitchOption interfaceAllowHadronFinalStatesNone
    (interfaceAllowHadronFinalStates,
     "None",
     "Option for disabling hadron final states of cluster fission",
     0);
  static SwitchOption interfaceAllowHadronFinalStatesSemiHadronicOnly
    (interfaceAllowHadronFinalStates,
     "SemiHadronicOnly",
     "Option for allowing hadron final states of cluster fission of type C->H1,C2 "
		 "(NOT YET USABLE WITH MatrixElement only use option None)",
     1);
  static SwitchOption interfaceAllowHadronFinalStatesAll
    (interfaceAllowHadronFinalStates,
     "All",
     "Option for allowing hadron final states of cluster fission "
		 "of type C->H1,C2 or C->H1,H2 "
		 "(NOT YET USABLE WITH MatrixElement only use option None)",
     2);

  // Mass Sampler Switch
  static Switch<MatrixElementClusterFissioner,int> interfaceMassSampler
    ("MassSampler",
     "Option for different Mass sampling options",
     &MatrixElementClusterFissioner::_massSampler, 0, false, false);
  static SwitchOption interfaceMassSamplerDefault
    (interfaceMassSampler,
     "Default",
     "Choose H7.2.3 default mass sampling using PSplitX",
     0);
  static SwitchOption interfaceMassSamplerUniform
    (interfaceMassSampler,
     "Uniform",
     "Choose Uniform Mass sampling in M1,M2 space",
     1);
  static SwitchOption interfaceMassSamplerFlatPhaseSpace
    (interfaceMassSampler,
     "FlatPhaseSpace",
     "Choose Flat Phase Space sampling of Mass in M1,M2 space (Recommended)",
     2);
  static SwitchOption interfaceMassSamplerPhaseSpacePowerLaw
    (interfaceMassSampler,
     "PhaseSpacePowerLaw",
     "Choose Phase Space times a power law sampling of Mass in M1,M2 space.",
     3);

  static Parameter<MatrixElementClusterFissioner,double>
    interfaceMassPowerLaw("MassPowerLaw",
				"power of mass power law modulation of FlatPhaseSpace",
                    &MatrixElementClusterFissioner::_powerLawPower, 0.0, -10.0, 0.0, 10.0,
		    false,false,false);

  // Phase Space Sampler Switch
  static Switch<MatrixElementClusterFissioner,int> interfacePhaseSpaceSamplerCluster
    ("PhaseSpaceSamplerCluster",
     "Option for different phase space sampling options",
     &MatrixElementClusterFissioner::_phaseSpaceSamplerCluster, 0, false, false);
  static SwitchOption interfacePhaseSpaceSamplerClusterAligned
    (interfacePhaseSpaceSamplerCluster,
     "Aligned",
     "Draw the momentum of child clusters to be aligned"
		 " the relative momentum of the mother cluster",
     0);
  static SwitchOption interfacePhaseSpaceSamplerClusterIsotropic
    (interfacePhaseSpaceSamplerCluster,
     "Isotropic",
     "Draw the momentum of child clusters to be isotropic"
		 " the relative momentum of the mother cluster",
     1);
  static SwitchOption interfacePhaseSpaceSamplerClusterTchannel
    (interfacePhaseSpaceSamplerCluster,
     "Tchannel",
     "Draw the momentum of child clusters to be t-channel like"
		 " the relative momentum of the mother cluster"
		 " i.e. cosTheta ~ 1/(1+A-cosTheta)^2",
     2);

  static Switch<MatrixElementClusterFissioner,int> interfacePhaseSpaceSamplerConstituents
    ("PhaseSpaceSamplerConstituents",
     "Option for different phase space sampling options",
     &MatrixElementClusterFissioner::_phaseSpaceSamplerConstituent, 0, false, false);
  static SwitchOption interfacePhaseSpaceSamplerConstituentsAligned
    (interfacePhaseSpaceSamplerConstituents,
     "Aligned",
     "Draw the momentum of the constituents of the child clusters "
		 "to be aligned relative to the direction of the child cluster",
     0);
  static SwitchOption interfacePhaseSpaceSamplerConstituentsIsotropic
    (interfacePhaseSpaceSamplerConstituents,
     "Isotropic",
     "Draw the momentum of the constituents of the child clusters "
		 "to be isotropic relative to the direction of the child cluster",
     1);
  static SwitchOption interfacePhaseSpaceSamplerConstituentsTchannel
    (interfacePhaseSpaceSamplerConstituents,
     "Exponential",
     "Draw the momentum of the constituents of the child clusters "
		 "to be exponential relative to the direction of the child cluster"
		 " i.e. cosTheta ~ exp(lam*(cosTheta-1)",
     2);



  // Matrix Element Choice Switch
  static Switch<MatrixElementClusterFissioner,int> interfaceMatrixElement
    ("MatrixElement",
     "Option for different Matrix Element options",
     &MatrixElementClusterFissioner::_matrixElement, 0, false, false);
  static SwitchOption interfaceMatrixElementDefault
    (interfaceMatrixElement,
     "Default",
     "Choose trivial matrix element i.e. whatever comes from the mass and "
		 "phase space sampling",
     0);
  static SwitchOption interfaceMatrixElementSoftQQbarFinalFinal
    (interfaceMatrixElement,
     "SoftQQbarFinalFinal",
		 "Choose Soft p1,p2->q1,q2,g*->q1,q2,q,qbar matrix element"
		 "NOTE: Here the matrix element depends on qi.q(bar)",
     1);
  static SwitchOption interfaceMatrixElementSoftQQbarInitialFinal
    (interfaceMatrixElement,
     "SoftQQbarInitialFinal",
		 "Choose Soft p1,p2->q1,q2,g*->q1,q2,q,qbar matrix element "
		 "NOTE: Here the matrix element depends on pi.q(bar)",
     2);
  static SwitchOption interfaceMatrixElementTest
    (interfaceMatrixElement,
     "Test",
		 "Choose Soft p1,p2->q1,q2,g*->q1,q2,q,qbar matrix element "
		 "NOTE: Here the matrix element depends on pi.q(bar)",
     4);
  static SwitchOption interfaceMatrixElementSoftQQbarInitialFinalTchannel
    (interfaceMatrixElement,
     "SoftQQbarInitialFinalTchannel",
		 "Choose Soft p1,p2->q1,q2,g*->q1,q2,q,qbar matrix element "
		 "NOTE: Here the matrix element depends on pi.q(bar)",
     5);
  static SwitchOption interfaceMatrixElementSoftQQbarInitialFinalTchannelOnlyMassDep
    (interfaceMatrixElement,
     "SoftQQbarInitialFinalTchannelOnlyMassDep",
		 "Choose Soft p1,p2->q1,q2,g*->q1,q2,q,qbar matrix element, but ignoring the dependence"
		 "on the angles NOTE: Here the matrix element depends on pi.q(bar)",
     6);
	// Technical Max loop parameter for New ClusterFission approach and Matrix Element 
	// Rejection sampling
  static Parameter<MatrixElementClusterFissioner,int> interfaceMaxLoopMatrixElement
    ("MaxLoopMatrixElement",
     "Technical Parameter for how many tries are allowed to sample the "
		 "Cluster Fission matrix element before reverting to fissioning "
		 "using the default Fission Aproach",
     &MatrixElementClusterFissioner::_maxLoopFissionMatrixElement, 5000000, 100, 1e8,
		 false, false, Interface::limited);
  static Switch<MatrixElementClusterFissioner,int> interfaceFailModeMaxLoopMatrixElement
    ("FailModeMaxLoopMatrixElement",
     "How to fail if we try more often than MaxLoopMatrixElement",
     &MatrixElementClusterFissioner::_failModeMaxLoopMatrixElement, 0, false, false);
  static SwitchOption interfaceDiquarkClusterFissionOldFission
	(interfaceFailModeMaxLoopMatrixElement,
     "OldFission",
     "Use old Cluster Fission Aproach if we try more often than MaxLoopMatrixElement",
     0);
  static SwitchOption interfaceDiquarkClusterFissionRejectEvent
	(interfaceFailModeMaxLoopMatrixElement,
     "RejectEvent",
     "Reject Event if we try more often than MaxLoopMatrixElement",
     1);

  static Switch<MatrixElementClusterFissioner,int> interfaceDiquarkClusterFission
    ("DiquarkClusterFission",
     "Allow clusters to fission to 1 or 2 diquark Clusters or Turn off diquark fission completely",
     &MatrixElementClusterFissioner::_diquarkClusterFission, 0, false, false);
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

  // Matrix Element Choice Switch
  static Parameter<MatrixElementClusterFissioner,double>
    interfaceSafetyFactorOverEst("SafetyFactorOverEst","Safety factor with which the numerical overestimate is calculated",
		  &MatrixElementClusterFissioner::_safetyFactorMatrixElement, 0, 10.0, 1.0, 1000.0,false,false,false);
  // Matrix Element IR cutoff
  static Parameter<MatrixElementClusterFissioner,double>
		interfaceMatrixElementIRCutOff("MatrixElementResolutionCutOff",
				"Resolution CutOff for MatrixElement t-channel with the Coloumb divergence. "
				"for MatrixElement SoftQQbarInitialFinalTchannel e.g. such that 1/(t^2)->1/((t-_epsilonResolution*mGluonConst^2)^2)",
		  &MatrixElementClusterFissioner::_epsilonResolution, 0, 0.01, 1e-6, 1000.0,false,false,false);
	
  // Matrix Element Choice Switch
  static Switch<MatrixElementClusterFissioner,int> interfaceWriteOut
    ("WriteOut",
     "Option for different Matrix Element options",
     &MatrixElementClusterFissioner::_writeOut, 0, false, false);
  static SwitchOption interfaceWriteOutNo
    (interfaceWriteOut,
     "No",
     "Choose trivial matrix element i.e. whatever comes from the mass and "
	 "phase space sampling",
     0);
  static SwitchOption interfaceWriteOutYes
    (interfaceWriteOut,
     "Yes",
	 "Choose Soft q1,q2->q1,q2,g*->q1,q2,q,qbar matrix element",
     1);
}
tPVector MatrixElementClusterFissioner::fission(ClusterVector & clusters, bool softUEisOn) {
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

void MatrixElementClusterFissioner::cut(stack<ClusterPtr> & clusterStack,
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
MatrixElementClusterFissioner::cutType
MatrixElementClusterFissioner::cutTwo(ClusterPtr & cluster, tPVector & finalhadrons,
			 bool softUEisOn) {
	switch (_fissionApproach)
	{
		case 0:
			return ClusterFissioner::cutTwo(cluster, finalhadrons, softUEisOn);
			break;
		case 1:
			return cutTwoMatrixElement(cluster, finalhadrons, softUEisOn);
			break;
		case 2:
			if (cluster->isBeamCluster()) {
				return ClusterFissioner::cutTwo(cluster, finalhadrons, softUEisOn);
			}
			else {
				return cutTwoMatrixElement(cluster, finalhadrons, softUEisOn);
			}
			break;
		case 3:
			{
				bool isPert=false;
				for (unsigned int i = 0; i < cluster->numComponents(); i++)
				{
					if (cluster->isPerturbative(i))
						isPert=true;
				}
				if (isPert)
					return ClusterFissioner::cutTwo(cluster, finalhadrons, softUEisOn);
				else
					return cutTwoMatrixElement(cluster, finalhadrons, softUEisOn);
			break;
			}
		case 4:
			{
				bool isPert=false;
				for (unsigned int i = 0; i < cluster->numComponents(); i++)
				{
					if (cluster->isPerturbative(i))
						isPert=true;
				}
				if (!isPert)
					return ClusterFissioner::cutTwo(cluster, finalhadrons, softUEisOn);
				else
					return cutTwoMatrixElement(cluster, finalhadrons, softUEisOn);
			break;
			}
		case 5:
			{
				bool isPert=true;
				for (unsigned int i = 0; i < cluster->numComponents(); i++)
				{
					if (!cluster->isPerturbative(i)) {
						isPert=false;
						break;
					}
				}
				if (isPert)
					return ClusterFissioner::cutTwo(cluster, finalhadrons, softUEisOn);
				else
					return cutTwoMatrixElement(cluster, finalhadrons, softUEisOn);
			break;
			}
		default:
			assert(false);
	}
}
namespace {
int areNotSame(const Lorentz5Momentum & p1,const Lorentz5Momentum & p2){
	double tol=1e-5;
	if (abs(p1.vect().x()-p2.vect().x())>tol*GeV){
		std::cout << "disagreeing x " <<std::setprecision(-log10(tol)+2) << abs(p1.vect().x()-p2.vect().x())/GeV<< std::endl;
		return 1;
	}
	if (abs(p1.vect().y()-p2.vect().y())>tol*GeV){
		std::cout << "disagreeing y " <<std::setprecision(-log10(tol)+2)<< abs(p1.vect().y()-p2.vect().y())/GeV<< std::endl;
		return 2;
	}
	if (abs(p1.vect().z()-p2.vect().z())>tol*GeV){
		std::cout << "disagreeing z " <<std::setprecision(-log10(tol)+2)<< abs(p1.vect().z()-p2.vect().z())/GeV<< std::endl;
		return 3;
	}
	if (abs(p1.e()-p2.e())>tol*GeV){
		std::cout << "disagreeing e " <<std::setprecision(-log10(tol)+2)<< abs(p1.e()-p2.e())/GeV<< std::endl;
		return 4;
	}
	if (abs(p1.m()-p2.m())>tol*GeV){
		std::cout << "disagreeing m " <<std::setprecision(-log10(tol)+2)<< abs(p1.m()-p2.m())/GeV<< std::endl;
		return 4;
	}
	if (abs(p1.m()-p1.mass())>tol*GeV){
		std::cout << "disagreeing p1 m - mass " <<std::setprecision(-log10(tol)+2)<< abs(p1.m()-p1.mass())/GeV<< std::endl;
		return 4;
	}
	if (abs(p2.m()-p2.mass())>tol*GeV){
		std::cout << "disagreeing p2 m - mass " <<std::setprecision(-log10(tol)+2)<< abs(p2.m()-p2.mass())/GeV<< std::endl;
		return 4;
	}
	return 0;
}
}
MatrixElementClusterFissioner::cutType
MatrixElementClusterFissioner::cutTwoMatrixElement(ClusterPtr & cluster, tPVector & finalhadrons,
			 bool softUEisOn) {
  // need to make sure only 2-cpt clusters get here
  assert(cluster->numComponents() == 2);
  tPPtr ptrQ1 = cluster->particle(0);
  tPPtr ptrQ2 = cluster->particle(1);
	Energy Mc = cluster->mass();
	// TODO BEGIN Changed comp to default
	if ( Mc < spectrum()->massLightestHadronPair(ptrQ1->dataPtr(),ptrQ2->dataPtr())) {
		static const PPtr null = PPtr();
		return cutType(PPair(null,null),PPair(null,null));
	}
	// TODO END Changed comp to default
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
  int counter_MEtry = 0;
  Energy Mc1 = ZERO;
	Energy Mc2 = ZERO;
	Energy m  = ZERO;
	Energy m1 = ptrQ1->data().constituentMass();
	Energy m2 = ptrQ2->data().constituentMass();
	Energy mMin = getParticleData(ParticleID::d)->constituentMass();
	// Minimal threshold for non-zero Mass PhaseSpace
	if ( Mc < (m1 + m2 + 2*mMin )) {
		static const PPtr null = PPtr();
		return cutType(PPair(null,null),PPair(null,null));
	}
  tcPDPtr toHadron1, toHadron2;
  PPtr newPtr1 = PPtr ();
  PPtr newPtr2 = PPtr ();
  Lorentz5Momentum pClu1, pClu2, pQ1, pQone, pQtwo, pQ2;
	Lorentz5Momentum pClu = cluster->momentum(); // known
	const Lorentz5Momentum p0Q1 = ptrQ1->momentum(); // known (mom Q1 before fission)
	const Lorentz5Momentum p0Q2 = ptrQ2->momentum(); // known (mom Q2 before fission)
	// where to return to in case of rejected sample
	enum returnTo {
		FlavourSampling,
		MassSampling,
		PhaseSpaceSampling,
		MatrixElementSampling,
		Done
	};
	// start with FlavourSampling
	returnTo retTo=FlavourSampling;
	int letFissionToXHadrons = _allowHadronFinalStates;
	// if a beam cluster not allowed to decay to hadrons
	if (cluster->isBeamCluster() && softUEisOn) letFissionToXHadrons = 0;
	// ### Flavour, Mass, PhaseSpace and MatrixElement Sampling loop until accepted: ###
	bool escape = false;
	double weightMasses,weightPhaseSpaceAngles;
	// TODO make this better
	// TODO make this independent of explicit u/d quark
	Energy mLightestQuark = getParticleData(ThePEG::ParticleID::u)->constituentMass();
	if (mLightestQuark>getParticleData(ThePEG::ParticleID::d)->constituentMass())
		mLightestQuark = getParticleData(ThePEG::ParticleID::d)->constituentMass();
	double SQME,SQMEoverEstimate;
	double weightFlatPS;
  do
  {
	  switch (retTo)
	  {
			case FlavourSampling:
			{
				// ## Flavour sampling and kinematic constraints ##
				drawNewFlavour(newPtr1,newPtr2,cluster);
				// get a flavour for the qqbar pair
				// check for right ordering
				// careful if DiquarkClusters can exist
				bool Q1diq = DiquarkMatcher::Check(ptrQ1->id());
				bool Q2diq = DiquarkMatcher::Check(ptrQ2->id());
				bool newQ1diq = DiquarkMatcher::Check(newPtr1->id());
				bool newQ2diq = DiquarkMatcher::Check(newPtr2->id());
				bool diqClu1 = Q1diq && newQ1diq;
				bool diqClu2 = Q2diq && newQ2diq;
				// DEBUG only:
				// std::cout << "Make Clusters: ( " << ptrQ1->PDGName() << " " << newPtr1->PDGName() << " ), ( "
				// << ptrQ2->PDGName() << " " << newPtr2->PDGName() << " )\n";
				// check if Hadron formation is possible
				if (!( diqClu1 || diqClu2 )
						&& (cantMakeHadron(ptrQ1, newPtr1) || cantMakeHadron(ptrQ2, newPtr2))) {
					swap(newPtr1, newPtr2);
					// check again
					if(cantMakeHadron(ptrQ1, newPtr1) || cantMakeHadron(ptrQ2, newPtr2)) {
						throw Exception()
							<< "MatrixElementClusterFissioner cannot split the cluster ("
							<< ptrQ1->PDGName() << ' ' << ptrQ2->PDGName()
							<< ") into hadrons.\n"
							<< "drawn Flavour: "<< newPtr1->PDGName()<<"\n"<< Exception::runerror;
					}
				}
				else if ( diqClu1 || diqClu2 ){
					bool swapped=false;
					if ( !diqClu1 && cantMakeHadron(ptrQ1,newPtr1) ) {
						swap(newPtr1, newPtr2);
						swapped=true;
					}
					if ( !diqClu2 && cantMakeHadron(ptrQ2,newPtr2) ) {
						assert(!swapped);
						swap(newPtr1, newPtr2);
					}
					if ( diqClu2 && diqClu1 && ptrQ1->id()*newPtr1->id()>0 ) {
						assert(!swapped);
						swap(newPtr1, newPtr2);
					}
					if (!diqClu1) assert(!cantMakeHadron(ptrQ1,newPtr1));
					if (!diqClu2) assert(!cantMakeHadron(ptrQ2,newPtr2));
				}
				// Check that new clusters can produce particles and there is enough
				// phase space to choose the drawn flavour
				m  = newPtr1->data().constituentMass();
				// Do not split in the case there is no phase space available + permille security
				double eps=0.1;
				if(Mc <  (1+eps)*(m1 + m + m2 + m))  {
					retTo = FlavourSampling;
					// escape if no flavour phase space possibile without fission
					if (fabs((m - mMin)/GeV) < 1e-3) {
						escape = true;
						retTo = Done;
					}
					continue;
				}
				pQ1.setMass(m1);
				pQone.setMass(m);
				pQtwo.setMass(m);
				pQ2.setMass(m2);
				// Determined the (5-components) momenta (all in the LAB frame)
				// p0Q1.setMass(ptrQ1->mass()); // known (mom Q1 before fission)
				// p0Q1.rescaleEnergy(); // TODO check if needed
				// p0Q2.setMass(ptrQ2->mass()); // known (mom Q2 before fission)
				// p0Q2.rescaleEnergy();// TODO check if needed
				pClu.rescaleMass();

				Energy MLHP1 = spectrum()->hadronPairThreshold(ptrQ1->dataPtr(),newPtr1->dataPtr());
				Energy MLHP2 = spectrum()->hadronPairThreshold(ptrQ2->dataPtr(),newPtr2->dataPtr());

				Energy MLH1 = spectrum()->lightestHadron(ptrQ1->dataPtr(),newPtr1->dataPtr())->mass();
				Energy MLH2 = spectrum()->lightestHadron(ptrQ2->dataPtr(),newPtr2->dataPtr())->mass();

				bool canBeSingleHadron1 = (m1 + m) < MLHP1;
				bool canBeSingleHadron2 = (m2 + m) < MLHP2;

				Energy mThresh1 = (m1 + m);
				Energy mThresh2 = (m2 + m);

				if (canBeSingleHadron1) mThresh1 = MLHP1;
				if (canBeSingleHadron2) mThresh2 = MLHP2;

				switch (letFissionToXHadrons)
				{
					case 0:
						{
							// Option None: only C->C1,C2 allowed
							// check if mass phase space is non-zero 
							// resample or escape if only allowed mass phase space is for C->H1,H2 or C->H1,C2
							if ( Mc < (mThresh1 + mThresh2)) {
								// escape if not even the lightest flavour phase space is possibile
								if ( fabs((m - mLightestQuark)/GeV) < 1e-3	) {
									escape = true;
									retTo = Done;
									continue;
								}
								else {
									retTo = FlavourSampling;
									continue;
								}
							}
							break;
						}
					case 1:
						{
							// Option SemiHadronicOnly: C->H,C allowed
							// NOTE: TODO implement matrix element for this
							// resample or escape if only allowed mass phase space is for C->H1,H2
							// First case is for ensuring the enough mass to  be available and second one rejects disjoint mass regions
							if ( ( (canBeSingleHadron1 && canBeSingleHadron2)
									    &&  Mc < (mThresh1 + mThresh2) )
									|| 
									 ( (canBeSingleHadron1 || canBeSingleHadron2)
								      && (canBeSingleHadron1 ? Mc-(m2+m) < MLH1:false ||  canBeSingleHadron2 ? Mc-(m1+m) < MLH2:false) )
								 ){
								// escape if not even the lightest flavour phase space is possibile
								if (fabs((m - mLightestQuark)/GeV) < 1e-3 ) {
									escape = true;
									retTo = Done;
									continue;
								}
								else {
									retTo = FlavourSampling;
									continue;
								}
							}
							break;
						}
					case 2:
						{
							// Option All: C->H,C and C->H,H allowed
							// NOTE: TODO implement matrix element for this
							// Mass Phase space for all option can always be found if cluster massive enough to go
							// to the lightest 2 hadrons
							break;
						}
					default:
						assert(false);
				}
				// Total overestimate of MatrixElement independent 
				// of the PhaseSpace point and independent of M1,M2
				SQMEoverEstimate = calculateSQME_OverEstimate(Mc,m1,m2,m); 
				// Note: want to fallthrough (in C++17 one could uncomment
				// 			 the line below to show that this is intentional)
				[[fallthrough]];
			}
			/*
			 * MassSampling choices:
			 * 	- Default (default)
			 * 	- Uniform
			 * 	- FlatPhaseSpace
			 * 	- SoftMEPheno
			 * 	*/
			case MassSampling:
			{
				weightMasses = drawNewMasses(Mc, soft1, soft2, pClu1, pClu2,
						ptrQ1, pQ1, newPtr1, pQone,
						newPtr2, pQtwo, ptrQ2, pQ2);

				// TODO IF C->C1,C2 (and not C->C,H or H1,H2) masses sampled and is in PhaseSpace must push through
				// because otherwise no matrix element
				if(weightMasses==0.0) {
					// TODO check which option is better
					retTo = FlavourSampling;
					// retTo = MassSampling;
					continue;
				}

				// derive the masses of the children
				Mc1 = pClu1.mass();
				Mc2 = pClu2.mass();
				// static kinematic threshold
				if(_kinematicThresholdChoice != 1 )
				{
					// NOTE: that the squared thresholds below are 
					// numerically more stringent which is needed
					// kaellen function calculation
					if(sqr(Mc1) < sqr(m1+m) || sqr(Mc2) < sqr(m+m2) || sqr(Mc1+Mc2) > sqr(Mc)){
						// TODO check which option is better
						// retTo = FlavourSampling;
						retTo = MassSampling;
						continue;
					}
					// dynamic kinematic threshold
				} 
				else if(_kinematicThresholdChoice == 1) {
					bool C1 = ( sqr(Mc1) )/( sqr(m1) + sqr(m) + _kinThresholdShift ) < 1.0 ? true : false;
					bool C2 = ( sqr(Mc2) )/( sqr(m2) + sqr(m) + _kinThresholdShift ) < 1.0 ? true : false;
					bool C3 = ( sqr(Mc1) + sqr(Mc2) )/( sqr(Mc) ) > 1.0 ? true : false;

					if( C1 || C2 || C3 ) {
						// TODO check which option is better
						// retTo = FlavourSampling;
						retTo = MassSampling;
						continue;
					}
				}
				else 
					assert(false);
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
				toHadron1 = spectrum()->chooseSingleHadron(ptrQ1->dataPtr(), newPtr1->dataPtr(),Mc1);
				if ( letFissionToXHadrons == 0 && toHadron1 ) {
					// reject mass samples which would force C->C,H or C->H1,H2 fission if desired 
					//
					// Check if Mc1max < MLHP1, in which case we might need to choose a different flavour
					// else resampling the masses should be sufficient
					Energy MLHP1 = spectrum()->massLightestHadronPair(ptrQ1->dataPtr(), newPtr1->dataPtr());
					Energy MLHP2 = spectrum()->massLightestHadronPair(ptrQ2->dataPtr(), newPtr2->dataPtr());
					Energy Mc1max = (m2+m) > MLHP2 ? (Mc-(m2+m)):(Mc-MLHP2);
					// for avoiding inf loops set min threshold
					if ( Mc1max - MLHP1 < 1e-2*GeV ) {
						retTo = FlavourSampling;
					}
					else {
						retTo = MassSampling;
					}
					continue;
				}
				if(toHadron1) {
					Mc1 = toHadron1->mass();
					pClu1.setMass(Mc1);
				}
				toHadron2 = spectrum()->chooseSingleHadron(ptrQ2->dataPtr(), newPtr2->dataPtr(),Mc2);
				if ( letFissionToXHadrons == 0 && toHadron2 ) {
					// reject mass samples which would force C->C,H or C->H1,H2 fission if desired 
					//
					// Check if Mc2max < MLHP2, in which case we might need to choose a different flavour
					// else resampling the masses should be sufficient
					Energy MLHP1 = spectrum()->massLightestHadronPair(ptrQ1->dataPtr(), newPtr1->dataPtr());
					Energy MLHP2 = spectrum()->massLightestHadronPair(ptrQ2->dataPtr(), newPtr2->dataPtr());
					Energy Mc2max = (m1+m) > MLHP1 ? (Mc-(m1+m)):(Mc-MLHP1);
					// for avoiding inf loops set min threshold
					if ( Mc2max - MLHP2 < 1e-2*GeV ) {
						retTo = FlavourSampling;
					}
					else {
						retTo = MassSampling;
					}
					continue;
				}
				if(toHadron2) {
					Mc2 = toHadron2->mass();
					pClu2.setMass(Mc2);
				}
				if (letFissionToXHadrons == 1 && (toHadron1 && toHadron2) ) {
					// reject mass samples which would force C->H1,H2 fission if desired 
					// TODO check which option is better
					// retTo = FlavourSampling;
					retTo = MassSampling;
					continue;
				}
				// Check if the decay kinematics is still possible: if not then
				// force the one-hadron decay for the other cluster as well.
				// NOTE: that the squared thresholds below are 
				// numerically more stringent which is needed
				// kaellen function calculation
				if(sqr(Mc1 + Mc2)  >  sqr(Mc)) {
					// reject if we would need to create a C->H,H process if letFissionToXHadrons<2
					if (letFissionToXHadrons < 2) {
						// TODO check which option is better
						// retTo = FlavourSampling;
						retTo = MassSampling;
						continue;
					}

					// TODO forbid other cluster!!!! to be also at hadron
					if(!toHadron1) {
						toHadron1 = spectrum()->chooseSingleHadron(ptrQ1->dataPtr(), newPtr1->dataPtr(),Mc-Mc2);
						// toHadron1 = spectrum()->chooseSingleHadron(ptrQ1->dataPtr(), newPtr1->dataPtr(),ZERO);
						if(toHadron1) {
							Mc1 = toHadron1->mass();
							pClu1.setMass(Mc1);
						}
					}
					else if(!toHadron2) {
						toHadron2 = spectrum()->chooseSingleHadron(ptrQ2->dataPtr(), newPtr2->dataPtr(),Mc-Mc1);
						// toHadron2 = spectrum()->chooseSingleHadron(ptrQ2->dataPtr(), newPtr2->dataPtr(),ZERO);
						if(toHadron2) {
							Mc2 = toHadron2->mass();
							pClu2.setMass(Mc2);
						}
					}
				}
				// NOTE: that the squared thresholds below are 
				// numerically more stringent which is needed
				// kaellen function calculation
				if (sqr(Mc) <= sqr(Mc1+Mc2)){
					// escape if already lightest quark drawn
					if (fabs((m - mLightestQuark)/GeV) < 1e-3 ) {
						// escape = true;
						retTo = Done;
					}
					else {
						// Try again with lighter quark
						retTo = FlavourSampling;
					}
					continue;
				}
				// weight(M1,M2) for M1*M2*(Two body PhaseSpace)**3 should be in [0,1]
				weightFlatPS = weightFlatPhaseSpace(Mc, Mc1, Mc2, m, m1, m2, ptrQ1, ptrQ2, newPtr1);
				// Note: want to fallthrough (in C++17 one could uncomment
				// 			 the line below to show that this is intentional)
				[[fallthrough]];
			}
			/*
			 * PhaseSpaceSampling choices:
			 * 	- FullyAligned (default)
			 * 	- AlignedIsotropic
			 * 	- FullyIsotropic
			 * 	*/
			case PhaseSpaceSampling:
			{
				// ### Sample the Phase Space with respect to Matrix Element: ###
				// TODO insert here PhaseSpace sampler
				weightPhaseSpaceAngles = drawKinematics(pClu,p0Q1,p0Q2,toHadron1,toHadron2,
						pClu1,pClu2,pQ1,pQone,pQtwo,pQ2);
				if(weightPhaseSpaceAngles==0.0) {
					// TODO check which option is better
					retTo = MassSampling;
					continue;
				}
				// Activate only if needed
				if (0)
				{ //sanity
					if (areNotSame(pClu1+pClu2,pClu)) {
						std::cout << "PCLU" << std::endl;
					}
					if (areNotSame(pClu1,pQ1+pQone)) {
						std::cout << "PCLU1" << std::endl;
					}
					if (areNotSame(pClu2,pQ2+pQtwo)) {
						std::cout << "PCLU2" << std::endl;
					}
					if (areNotSame(pClu,pQ1+pQone+pQ2+pQtwo)) {
						std::cout << "PTOT" << std::endl;
					}
				}

				// Should be precise i.e. no rejection expected
				// Note: want to fallthrough (in C++17 one could uncomment
				// 			 the line below to show that this is intentional)
				[[fallthrough]];
			}
			/*
			 * MatrixElementSampling choices:
			 * 	- Default (default)
			 * 	- SoftQQbar
			 * 	*/
			case MatrixElementSampling:
			{
				counter_MEtry++;
				// TODO maybe bridge this to work more neatly
				// Ignore matrix element for C->C,H or C->H1,H2 fission
				if (toHadron1 || toHadron2) {
					retTo = Done;
					break;
				}
				// Actual MatrixElement evaluated at sampled PhaseSpace point
				SQME = calculateSQME(p0Q1,p0Q2,pQ1,pQone,pQ2,pQtwo);
				// weight for MatrixElement*PhaseSpace must be in [0:1]
				double weightSQME = SQME/SQMEoverEstimate;
				assert(weightSQME>0.0);
				if (weightFlatPS<0 || weightFlatPS>1.0){
					throw Exception() << "weightFlatPS = "<< weightFlatPS << " > 1 or negative in MatrixElementClusterFissioner::cutTwo"
						<< "Mc  = " << Mc/GeV
						<< "Mc1 = " << Mc1/GeV
						<< "Mc2 = " << Mc2/GeV
						<< "m1  = " << m1/GeV
						<< "m2  = " << m2/GeV
						<< "m   = " << m/GeV
						<< "SQME   = " << SQME
						<< "SQME_OE   = " << SQMEoverEstimate
						<< "PS   = " << weightFlatPS
						<< Exception::runerror;
				}
				// current phase space point is distributed according to weightSamp
				double Pacc = weightFlatPS * weightSQME;
				double SamplingWeight = weightMasses * weightPhaseSpaceAngles;
				if (!(SamplingWeight >= 0 ) || std::isnan(SamplingWeight) || std::isinf(SamplingWeight)){
					throw Exception() << "SamplingWeight = "<< SamplingWeight << " in MatrixElementClusterFissioner::cutTwo"
						<< "Mc  = " << Mc/GeV
						<< "Mc1 = " << Mc1/GeV
						<< "Mc2 = " << Mc2/GeV
						<< "m1  = " << m1/GeV
						<< "m2  = " << m2/GeV
						<< "m   = " << m/GeV
						<< "weightMasses   = " << weightMasses
						<< "weightPhaseSpaceAngles   = " << weightPhaseSpaceAngles
						<< "PS   = " << weightFlatPS
						<< Exception::runerror;
				}
				Pacc/=SamplingWeight;
				if (!(Pacc >= 0 ) || std::isnan(Pacc) || std::isinf(Pacc)){
					throw Exception() << "Pacc = "<< Pacc << " < 0 in MatrixElementClusterFissioner::cutTwo"
						<< "Mc  = " << Mc/GeV
						<< "Mc1 = " << Mc1/GeV
						<< "Mc2 = " << Mc2/GeV
						<< "m1  = " << m1/GeV
						<< "m2  = " << m2/GeV
						<< "m   = " << m/GeV
						<< "SQME   = " << SQME
						<< "SQME_OE   = " << SQMEoverEstimate
						<< "PS   = " << weightFlatPS
						<< Exception::runerror;
				}
				assert(Pacc >= 0.0);
				if (_matrixElement!=0) {
					if ( Pacc > 1.0){
						countPaccGreater1();
						throw Exception() << "Pacc = "<< Pacc << " > 1 in MatrixElementClusterFissioner::cutTwo"
							<< "Mc  = " << Mc/GeV
							<< "Mc1 = " << Mc1/GeV
							<< "Mc2 = " << Mc2/GeV
							<< "m1  = " << m1/GeV
							<< "m2  = " << m2/GeV
							<< "m   = " << m/GeV
							<< Exception::warning;
					}
				}
				static int first=_writeOut;
				if (_matrixElement==0  || UseRandom::rnd()<Pacc) { //Always accept a sample if trivial matrix element
					if (_writeOut) {
						// std::cout << "\nAccept Pacc = "<<Pacc<<"\n";
						std::ofstream out("data_CluFis.dat", std::ios::app | std::ios::out);
						Lorentz5Momentum p0Q1tmp(p0Q1);
						Lorentz5Momentum pClu1tmp(pClu1);
						Lorentz5Momentum pClu2tmp(pClu2);
						p0Q1tmp.boost(-pClu.boostVector());
						pClu1tmp.boost(-pClu.boostVector());
						double cosThetaP1C1=p0Q1tmp.vect().cosTheta(pClu1tmp.vect());
						out << Pacc << "\t"
							<< Mc/GeV  << "\t"
							<< pClu1.mass()/GeV << "\t"
							<< pClu2.mass()/GeV << "\t"
							<< pQone*pQtwo/(pQone.mass()*pQtwo.mass()) << "\t"
							<< pQ1*pQ2/(pQ1.mass()*pQ2.mass()) << "\t"
							<< pQ1*pQtwo/(pQ1.mass()*pQtwo.mass()) << "\t"
							<< pQ2*pQone/(pQ2.mass()*pQone.mass()) << "\t"
							<< p0Q1*(pQone+pQtwo)/(p0Q1.mass()*(pQone.mass()+pQtwo.mass())) << "\t"
							<< p0Q2*(pQone+pQtwo)/(p0Q2.mass()*(pQone.mass()+pQtwo.mass())) << "\t"
							<< m/GeV << "\t"
							<< cosThetaP1C1
							<< "\n";
						out.close();
						Energy MbinStart=2.0*GeV;
						Energy MbinEnd=91.0*GeV;
						Energy dMbin=1.0*GeV;
						Energy MbinLow;
						Energy MbinHig;
						if (   fabs((m-getParticleData(ParticleID::d)->constituentMass())/GeV)<1e-5
							&& fabs((m1-getParticleData(ParticleID::d)->constituentMass())/GeV)<1e-5
							&& fabs((m2-getParticleData(ParticleID::d)->constituentMass())/GeV)<1e-5) {
							int cnt = 0;
							for (Energy MbinIt=MbinStart; MbinIt < MbinEnd; MbinIt+=dMbin)
							{
								if (first){
									first=0;
									std::cout << "\nFirst\n" << std::endl;
									int ctr=0;
									for (Energy MbinIt=MbinStart; MbinIt < MbinEnd; MbinIt+=dMbin)
									{
										std::string name = "WriteOut/data_CluFis_BinM_"+std::to_string(ctr)+".dat";
										std::ofstream out2(name,std::ios::out);
										out2 << "# Binned from "<<MbinIt/GeV << "\tto\t" << (MbinIt+dMbin)/GeV<<"\n";
										out2 << "# m=m1=m2=0.325\n";
										out2 << "# (1) Pacc\t"
											<< "(2) Mc/GeV\t"
											<< "(3) M1/GeV\t"
											<< "(4) M2/GeV\t"
											<< "(5) q.qbar/(mq^2)\t"
											<< "(6) q1.q2/(m1*m2)\t"
											<< "(7) q1.qbar/(m1*mq)\t"
											<< "(8) q2.q/(m2*mq)\t"
											<< "(9) p1.(q+qbar)/(m1*2*mq)\t"
											<< "(10) p2.(q+qbar)/(m2*2*mq)\t"
											<< "(11) m/GeV\t"
											<< "(12) cos(theta_p1_Q1)\n";
										out2.close();
										ctr++;
									}
									// first=0;
								}
								MbinLow  = MbinIt;
								MbinHig  = MbinLow+dMbin;
								if (Mc>MbinLow && Mc<MbinHig) {
									std::string name = "WriteOut/data_CluFis_BinM_"+std::to_string(cnt)+".dat";
									std::ofstream out2(name, std::ios::app );
									Lorentz5Momentum p0Q1tmp(p0Q1);
									Lorentz5Momentum pClu1tmp(pClu1);
									Lorentz5Momentum pClu2tmp(pClu2);
									p0Q1tmp.boost(-pClu.boostVector());
									pClu1tmp.boost(-pClu.boostVector());
									double cosThetaP1C1=p0Q1tmp.vect().cosTheta(pClu1tmp.vect());
									out2 << Pacc << "\t"
										<< Mc/GeV  << "\t"
										<< pClu1.mass()/GeV << "\t"
										<< pClu2.mass()/GeV << "\t"
										<< pQone*pQtwo/(pQone.mass()*pQtwo.mass()) << "\t"
										<< pQ1*pQ2/(pQ1.mass()*pQ2.mass()) << "\t"
										<< pQ1*pQtwo/(pQ1.mass()*pQtwo.mass()) << "\t"
										<< pQ2*pQone/(pQ2.mass()*pQone.mass()) << "\t"
										<< p0Q1*(pQone+pQtwo)/(p0Q1.mass()*(pQone.mass()+pQtwo.mass())) << "\t"
										<< p0Q2*(pQone+pQtwo)/(p0Q2.mass()*(pQone.mass()+pQtwo.mass())) << "\t"
										<< m/GeV << "\t"
										<< cosThetaP1C1
										<< "\n";
									out2.close();
									break;
								}
								// if (Mc<MbinIt) break;
								cnt++;
							}
						}
					}
					retTo=Done;
					break;
				}
				retTo = MassSampling;
				continue;
			}
			default:
			{
				assert(false);
			}
	  }
  }
  while (retTo!=Done && !escape && counter_MEtry < _maxLoopFissionMatrixElement);

  if(escape) {
		// happens if we get at too light cluster to begin with
		static const PPtr null = PPtr();
		return cutType(PPair(null,null),PPair(null,null));
  }

  if(counter_MEtry >= _maxLoopFissionMatrixElement) {
		countMaxLoopViolations();
		// happens if we get too massive clusters where the matrix element
		// is very inefficiently sampled
		std::stringstream warning;
		warning
			<< "Matrix Element rejection sampling tried more than "
			<< counter_MEtry
			<< " times.\nMcluster = " << Mc/GeV
			<< "GeV\nm1 = " << m1/GeV
			<< "GeV\nm2 = " << m2/GeV
			<< "GeV\nm  = " << m/GeV
			<< "GeV\n"
			<< "isBeamCluster = " << cluster->isBeamCluster()
			<<"Using default as MatrixElementClusterFissioner::cutTwo as a fallback.\n"
			<< "***This Exception should not happen too often!*** ";
		generator()->logWarning( Exception(warning.str(),Exception::warning));
		switch (_failModeMaxLoopMatrixElement)
		{
			case 0:
				// OldFission
				return ClusterFissioner::cutTwo(cluster,finalhadrons,softUEisOn);
				break;
			case 1:
				// RejectEvent
				throw Exception() << warning.str()
					<< Exception::eventerror;
				break;
			default:
				assert(false);
		}
  }
	assert(abs(Mc1-pClu1.m())<1e-2*GeV);
	assert(abs(Mc2-pClu2.m())<1e-2*GeV);
	assert(abs(Mc1-pClu1.mass())<1e-2*GeV);
	assert(abs(Mc2-pClu2.mass())<1e-2*GeV);
	{
		Lorentz5Momentum pp1(p0Q1);
		Lorentz5Momentum pclu1(pClu1);
		pclu1.boost(-pClu.boostVector());
		pp1.boost(-pClu.boostVector());
		// std::ofstream out("testCF.dat", std::ios::app);
		// out << pclu1.vect().cosTheta(pp1.vect())<<"\n";
		// out.close();
 }
	// Count successfull ClusterFission
	countFissionMatrixElement();
  // ==> full sample generated
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

/**
 * Calculate the masses and possibly kinematics of the cluster
 * fission at hand; if claculateKineamtics is perfomring non-trivial
 * steps kinematics calulcated here will be overriden. Currently resorts to the default
 * @return the potentially non-trivial distribution weight=f(M1,M2)
 *         On Failure we return 0
 */ 
double MatrixElementClusterFissioner::drawNewMasses(const Energy Mc, const bool soft1, const bool soft2,
		Lorentz5Momentum& pClu1, Lorentz5Momentum& pClu2,
		tcPPtr ptrQ1,    const Lorentz5Momentum& pQ1, 
		tcPPtr newPtr1, const Lorentz5Momentum& pQone,
		tcPPtr, const Lorentz5Momentum& pQtwo,
		tcPPtr ptrQ2,    const Lorentz5Momentum& pQ2) const {
	// TODO	add precise weightMS that could be used used for improving the rejection Sampling
	switch (_massSampler)
	{
		case 0:
			return ClusterFissioner::drawNewMasses(Mc, soft1, soft2, pClu1, pClu2, ptrQ1, pQ1, tcPPtr(), pQone, tcPPtr(),pQtwo, ptrQ2, pQ2);
			break;
		case 1: 
			return drawNewMassesUniform(Mc, pClu1, pClu2, pQ1, pQone, pQ2);
			break;
		case 2:
			return drawNewMassesPhaseSpace(Mc, pClu1, pClu2, pQ1, pQone, pQ2, ptrQ1, newPtr1, ptrQ2);
			break;
		case 3:
			return drawNewMassesPhaseSpacePowerLaw(Mc, pClu1, pClu2, pQ1, pQone, pQ2, ptrQ1, newPtr1, ptrQ2);
			break;
		default:
			assert(false);
	}
	return 0;// failure
}


/**
 * Sample the masses for flat phase space
 * */
double MatrixElementClusterFissioner::drawNewMassesUniform(const Energy Mc, Lorentz5Momentum& pClu1, Lorentz5Momentum& pClu2,
		const Lorentz5Momentum& pQ1, 
		const Lorentz5Momentum& pQ,
		const Lorentz5Momentum& pQ2) const {

	Energy M1,M2;
	const Energy m1 = pQ1.mass();
	const Energy m2 = pQ2.mass();
	const Energy m  = pQ.mass();
	const Energy M1min = m1 + m;
	const Energy M2min = m2 + m;
	const Energy M1max = Mc - M2min;
	const Energy M2max = Mc - M1min;

	assert(M1max-M1min>ZERO);
	assert(M2max-M2min>ZERO);

	double r1;
	double r2;

	int counter = 0;
	const int max_counter = 100;

	while (counter < max_counter) {
		r1 = UseRandom::rnd();
		r2 = UseRandom::rnd();

		M1 = (M1max-M1min)*r1 + M1min;
		M2 = (M2max-M2min)*r2 + M2min;

		counter++;
		if ( Mc > M1 + M2) break;
	}

	if (counter==max_counter
			|| Mc < M1 + M2
			|| M1 <= M1min
			|| M2 <= M2min ) return 0.0; // failure

	pClu1.setMass(M1);
	pClu2.setMass(M2);

	return 1.0; // succeeds
}
/**
 * Sample the masses for flat phase space
 * */
double MatrixElementClusterFissioner::drawNewMassesPhaseSpacePowerLaw(const Energy Mc,
		Lorentz5Momentum& pClu1, Lorentz5Momentum& pClu2,
		const Lorentz5Momentum& pQ1, 
		const Lorentz5Momentum& pQ,
		const Lorentz5Momentum& pQ2,
		tcPPtr ptrQ1, tcPPtr ptrQ2, tcPPtr ptrQ) const {
	Energy M1,M2,MuS;
	const Energy m1 = pQ1.mass();
	const Energy m2 = pQ2.mass();
	const Energy m  = pQ.mass();
	const Energy M1min = m1 + m;
	const Energy M2min = m2 + m;
	const Energy M1max = Mc - M2min;
	const Energy M2max = Mc - M1min;

	assert(M1max-M1min>ZERO);
	assert(M2max-M2min>ZERO);

	double r1;
	double r2;

	int counter = 0;
	const int max_counter = 200;

	while (counter < max_counter) {
		r1 = UseRandom::rnd();
		r2 = UseRandom::rnd();

		M1 = (M1max-M1min)*r1 + M1min;
		M2 = (M2max-M2min)*r2 + M2min;

		counter++;
		if (sqr(M1+M2)>sqr(Mc))
			continue;

		// if (!phaseSpaceVeto(Mc,M1,M2,m,m1,m2) ) break; // For FlatPhaseSpace sampling vetoing
		if (!phaseSpaceVeto(Mc,M1,M2,m,m1,m2, ptrQ1, ptrQ2, ptrQ,_powerLawPower) ) break; // For FlatPhaseSpace sampling vetoing
	}
	if (counter==max_counter) return 0.0; // failure

	pClu1.setMass(M1);
	pClu2.setMass(M2);

	return weightPhaseSpaceConstituentMasses(Mc, M1, M2, m, m1, m2, _powerLawPower); // succeeds return weight
}


/**
 * Sample the masses for flat phase space
 * */
double MatrixElementClusterFissioner::drawNewMassesPhaseSpace(const Energy Mc,
		Lorentz5Momentum& pClu1, Lorentz5Momentum& pClu2,
		const Lorentz5Momentum& pQ1, 
		const Lorentz5Momentum& pQ,
		const Lorentz5Momentum& pQ2,
		tcPPtr ptrQ1, tcPPtr ptrQ2, tcPPtr ptrQ ) const {
	Energy M1,M2,MuS;
	const Energy m1 = pQ1.mass();
	const Energy m2 = pQ2.mass();
	const Energy m  = pQ.mass();
	const Energy M1min = m1 + m;
	const Energy M2min = m2 + m;
	const Energy M1max = Mc - M2min;
	const Energy M2max = Mc - M1min;

	assert(M1max-M1min>ZERO);
	assert(M2max-M2min>ZERO);

	double r1;
	double r2;

	int counter = 0;
	const int max_counter = 200;

	while (counter < max_counter) {
		r1 = UseRandom::rnd();
		r2 = UseRandom::rnd();

		M1 = (M1max-M1min)*r1 + M1min;
		M2 = (M2max-M2min)*r2 + M2min;

		counter++;
		if (sqr(M1+M2)>sqr(Mc))
			continue;

		// if (!phaseSpaceVeto(Mc,M1,M2,m,m1,m2) ) break; // For FlatPhaseSpace sampling vetoing
		if (!phaseSpaceVeto(Mc,M1,M2,m,m1,m2, ptrQ1, ptrQ2, ptrQ) ) break; // For FlatPhaseSpace sampling vetoing
	}
	if (counter==max_counter) return 0.0; // failure

	pClu1.setMass(M1);
	pClu2.setMass(M2);

	return weightFlatPhaseSpace(Mc, M1, M2, m, m1, m2, ptrQ1, ptrQ2, ptrQ); // succeeds return weight
}
std::pair<Axis,double> MatrixElementClusterFissioner::sampleDirectionConstituents(
		const Lorentz5Momentum & pClu, const Energy Mcluster)  const {
		switch (_phaseSpaceSamplerConstituent)
		{
			case 0:
				{
					// Aligned
					return std::make_pair(pClu.vect().unit(),1.0);
				}
			case 1:
				{
					// Isotropic
					return std::make_pair(sampleDirectionIsotropic(),1.0);
				}
			case 2:
			{
				// Exponential
				// New
				// double C_est_Exponential=1.18;
				// double power_est_Exponential=0.65;
				// Old 
				double C_est_Exponential=0.63;
				double power_est_Exponential=0.89;
				double lambda=C_est_Exponential*pow(Mcluster/GeV,power_est_Exponential);
				return sampleDirectionExponential(pClu.vect().unit(), lambda);
			}
			default:
				assert(false);
				break;
		}
	return std::make_pair(Axis(),0.0); // Failure
}

std::pair<Axis,double> MatrixElementClusterFissioner::sampleDirectionCluster(
		const Lorentz5Momentum & pQ,
		const Lorentz5Momentum & pClu)  const {
	switch (_phaseSpaceSamplerCluster)
	{
		case 0:
			{
				// Aligned
				Axis dir;
				if (_covariantBoost)
					// in Covariant Boost the positive z-Axis is defined as the direction of
					// the pQ vector in the Cluster rest frame
					dir = Axis(0,0,1);
				else 
					dir = sampleDirectionAligned(pQ, pClu);
				return std::make_pair(dir,1.0);
			}
		case 1:
			{
				// Isotropic
				return std::make_pair(sampleDirectionIsotropic(),1.0);
			}
		case 2:
			{
				// Tchannel
				Energy M=pClu.mass();
				// New
				// double C_est_Tchannel=3.66;
				// double power_est_Tchannel=-1.95;
				// Old
				double C_est_Tchannel=9.78;
				double power_est_Tchannel=-2.5;
				double A = C_est_Tchannel*pow(M/GeV,power_est_Tchannel);
				double Ainv = 1.0/(1.0+A);
				return sampleDirectionTchannel(sampleDirectionAligned(pQ,pClu),Ainv);
			}
		default:
			assert(false);  
	}
	return std::make_pair(Axis(),0.0); // Failure
}

std::pair<Axis,double> MatrixElementClusterFissioner::sampleDirectionExponential(const Axis & dirQ, const double lambda) const {
	Axis FinalDir = dirQ;
	double cosTheta = Kinematics::sampleCosExp(lambda);
	double phi;
	// std::cout << "cosThetaSamp = "<< cosTheta <<"\tlambda = " <<lambda << std::endl;
	// If no change in angle keep the direction fixed
	if (fabs(cosTheta-1.0)>1e-14) {
		// rotate to sampled angles
		FinalDir.rotate(acos(cosTheta),dirQ.orthogonal());
		phi  = UseRandom::rnd(-Constants::pi,Constants::pi);
		FinalDir.rotate(phi,dirQ);
	}
	// std::cout << "cosThetaTrue = "<< FinalDir.cosTheta(dirQ) <<"\tlambda = " <<lambda << std::endl;
	double weight = exp(lambda*(cosTheta-1.0));
	if (std::isnan(weight) ||std::isinf(weight) )
		std::cout << "lambda = " << lambda << "\t cos " << cosTheta<< std::endl;
	assert(!std::isnan(weight));
	assert(!std::isinf(weight));
	if (fabs(FinalDir.cosTheta(dirQ)-cosTheta)>1e-8) std::cout << "cosThetaTrue = "<< FinalDir.cosTheta(dirQ) <<"\tlambda = " <<lambda << std::endl;
	return std::make_pair(FinalDir,weight);
}
std::pair<Axis,double> MatrixElementClusterFissioner::sampleDirectionTchannel(const Axis & dirQ, const double Ainv) const {
	double cosTheta = Kinematics::sampleCosTchannel(Ainv);
	double phi;
	Axis FinalDir = dirQ;
	// If no change in angle keep the direction fixed
	if (fabs(cosTheta-1.0)>1e-14) {
		// rotate to sampled angles
		FinalDir.rotate(acos(cosTheta),dirQ.orthogonal());
		phi  = UseRandom::rnd(-Constants::pi,Constants::pi);
		FinalDir.rotate(phi,dirQ);
	}
	double weight = 0.0;
	if (1.0-Ainv>1e-10)
		weight=1.0/pow(1.0/Ainv-cosTheta,2.0);
	else
		weight=1.0/pow((1.0-cosTheta)+(1.0-Ainv),2.0);
	if (std::isnan(weight) ||std::isinf(weight) )
		std::cout << "Ainv = " << Ainv << "\t cos " << cosTheta<< std::endl;
	assert(!std::isnan(weight));
	assert(!std::isinf(weight));
	if (fabs(FinalDir.cosTheta(dirQ)-cosTheta)>1e-8) std::cout << "cosThetaTrue = "<< FinalDir.cosTheta(dirQ) <<"\tAinv = " <<Ainv << std::endl;
	return std::make_pair(FinalDir,weight);
}

Axis MatrixElementClusterFissioner::sampleDirectionAligned(const Lorentz5Momentum & pQ, const Lorentz5Momentum & pClu) const {
	Lorentz5Momentum pQinCOM(pQ);
	pQinCOM.setMass(pQ.m());
	pQinCOM.boost( -pClu.boostVector() );        // boost from LAB to C
	return pQinCOM.vect().unit();
}

Axis MatrixElementClusterFissioner::sampleDirectionAlignedSmeared(const Lorentz5Momentum & pQ, const Lorentz5Momentum & pClu) const {
	Lorentz5Momentum pQinCOM(pQ);
	pQinCOM.setMass(pQ.m());
	pQinCOM.boost( -pClu.boostVector() );        // boost from LAB to C
	if (pQinCOM.vect().mag2()<=ZERO){
		return sampleDirectionAlignedSmeared(pClu-pQ,pClu);
	}
	const Axis dir = pQinCOM.vect().unit();
	double cluSmear = 0.5;
	double cosTheta;
	do {
		cosTheta = 1.0 + cluSmear*log( UseRandom::rnd() );
	} 
	while (fabs(cosTheta)>1.0 || std::isnan(cosTheta) || std::isinf(cosTheta));

	if (!(dir.mag2()>ZERO)){
		std::cout << "\nDRI = 0"<< dir.mag2() << std::endl;
	}
	Axis dirSmeared = dir;
	if (fabs(cosTheta-1.0)>1e-10) {
		// rotate to sampled angles
		dirSmeared.rotate(acos(cosTheta),dir.orthogonal());
		double phi  = UseRandom::rnd(-M_PI,M_PI);
		dirSmeared.rotate(phi,dir);
    if (!(dirSmeared.mag2()>ZERO)){
      std::cout << "\nDIR SMR = 0" <<cosTheta<< "\t" <<acos(cosTheta)<< "\t"<<  phi<< "\t"<< dirSmeared.mag2() << std::endl;
    }
	}
	if (!(dirSmeared.mag2()>ZERO)){
		std::cout << "\nDIR SMR = 0" <<cosTheta<< "\t" <<acos(cosTheta)<< "\t"<< dirSmeared.mag2() << std::endl;
		std::cout << dir.orthogonal() << std::endl;
	}
	return dirSmeared;
}

Axis MatrixElementClusterFissioner::sampleDirectionIsotropic() const {
  double cosTheta = -1 + 2.0 * UseRandom::rnd();
  double sinTheta = sqrt(1.0-cosTheta*cosTheta);
  double Phi = 2.0 * Constants::pi * UseRandom::rnd();
  Axis uClusterUniform(cos(Phi)*sinTheta, sin(Phi)*sinTheta, cosTheta);
  return uClusterUniform.unit();
}

Axis MatrixElementClusterFissioner::sampleDirectionSemiUniform(const Lorentz5Momentum & pQ, const Lorentz5Momentum & pClu) const {
  Axis dir = sampleDirectionAligned(pQ,pClu);
  Axis res;
  do {
	  res=sampleDirectionIsotropic();
  }
  while (dir*res<0);
  return res;
}
namespace {
	double SoftFunction(
			const Lorentz5Momentum & pi,
			const Lorentz5Momentum & pj,
			const Lorentz5Momentum & q,
			const Lorentz5Momentum & qbar) {
		Energy2 mq2 = q*q;
		Energy2 qqbar = q*qbar;
		Energy2 piq = pi*q;
		Energy2 piqbar = pi*qbar;
		Energy2 pjq = pj*q;
		Energy2 pjqbar = pj*qbar;
		double numerator   = -(pi*pj)*(qqbar+mq2)/sqr(GeV2);
		double denominator = sqr(qqbar+mq2)*(piq+piqbar)*(pjq+pjqbar)/sqr(GeV2*GeV2);
		int cataniScheme = 1;
		// Different expressions which should yield similar results
		switch (cataniScheme) {
			case 0:
				numerator += ((piq)*(pjqbar) + (pjq)*(piqbar))/sqr(GeV2);
				break;
			case 1:
				numerator += - (0.5*(piq-piqbar)*(pjq-pjqbar))/sqr(GeV2);
				break;
			default:
				assert(false);
		}
		double Iij = numerator/denominator;
		return Iij;
	}
}
/* SQME for p1,p2->C1(q1,q),C2(q2,qbar) 
 * Note that:
 *      p0Q1  -> p1
 *      p0Q2  -> p2
 *      pQ1   -> q1
 *      pQone -> q
 *      pQ2   -> q2
 *      pQone -> qbar
 * */
double MatrixElementClusterFissioner::calculateSQME(
		const Lorentz5Momentum & p1,
		const Lorentz5Momentum & p2,
		const Lorentz5Momentum & q1,
		const Lorentz5Momentum & q,
		const Lorentz5Momentum & q2,
		const Lorentz5Momentum & qbar) const {
	double SQME;
	switch (_matrixElement)
	{
		case 0:
				SQME = 1.0;
				break;
		case 1:
			{
				/*
				// Energy2 p1p2 = p1*p2;
				Energy2 q1q2    = q1 * q2;
				Energy2 q1q     = q1 * q ;
				Energy2 q2qbar  = q2 * qbar;
				Energy2 q2q     = q2 * q ;
				Energy2 q1qbar  = q1 * qbar;
				Energy2 qqbar   = q  * qbar;
				Energy2 mq2 = q.mass2();
				double Numerator = q1q2 * (qqbar + mq2)/sqr(GeV2);
				Numerator += 0.5 * (q1q - q1qbar)*(q2q - q2qbar)/sqr(GeV2);
				double Denominator = sqr(qqbar + mq2)*(q1q + q1qbar)*(q2q + q2qbar)/sqr(sqr(GeV2));
				SQME = Numerator/Denominator;
				*/
				double I11 = SoftFunction(q1,q1,q,qbar);
				double I22 = SoftFunction(q2,q2,q,qbar);
				double I12 = SoftFunction(q1,q2,q,qbar);
				SQME = (I11+I22-2.0*I12);
				if (SQME<0.0) {
					std::cout << "I11 = " << I11<< std::endl;
					std::cout << "I22 = " << I22<< std::endl;
					std::cout << "2*I12 = " << I12<< std::endl;
					std::cout << "soft = " << I11+I22-2.0*I12<< std::endl;
					std::cout << "M  = " << (p1+p2).m()/GeV<< std::endl;
					std::cout << "M1 = " << (q1+q).m()/GeV<< std::endl;
					std::cout << "M2 = " << (q2+qbar).m()/GeV<< std::endl;
					std::cout << "m1 = " << (q1).m()/GeV<< std::endl;
					std::cout << "m2 = " << (q2).m()/GeV<< std::endl;
					std::cout << "m  = " << (q).m()/GeV<< std::endl;
				}
				break;
			}
		case 2:
			{
				/*
				Energy2 p1p2    = p1 * p2;
				Energy2 p1q     = p1 * q ;
				Energy2 p2qbar  = p2 * qbar;
				Energy2 p2q     = p2 * q ;
				Energy2 p1qbar  = p1 * qbar;
				Energy2 qqbar   = q  * qbar;
				Energy2 mq2 = q.mass2();
				double Numerator = p1p2 * (qqbar + mq2)/sqr(GeV2);
				Numerator += 0.5 * (p1q - p1qbar)*(p2q - p2qbar)/sqr(GeV2);
				double Denominator = sqr(qqbar + mq2)*(p1q + p1qbar)*(p2q + p2qbar)/sqr(sqr(GeV2));
				SQME = Numerator/Denominator;
				*/
				double I11 = SoftFunction(p1,p1,q,qbar);
				double I22 = SoftFunction(p2,p2,q,qbar);
				double I12 = SoftFunction(p1,p2,q,qbar);
				SQME = (I11+I22-2.0*I12);
				if (SQME<0.0) {
					std::cout << "I11 = " << I11<< std::endl;
					std::cout << "I22 = " << I22<< std::endl;
					std::cout << "2*I12 = " << I12<< std::endl;
					std::cout << "soft = " << I11+I22-2.0*I12<< std::endl;
					std::cout << "M  = " << (p1+p2).m()/GeV<< std::endl;
					std::cout << "M1 = " << (q1+q).m()/GeV<< std::endl;
					std::cout << "M2 = " << (q2+qbar).m()/GeV<< std::endl;
					std::cout << "m1 = " << (q1).m()/GeV<< std::endl;
					std::cout << "m2 = " << (q2).m()/GeV<< std::endl;
					std::cout << "m  = " << (q).m()/GeV<< std::endl;
				}
				break;
			}
		case 4:
			{
				SQME = sqr(sqr(GeV2))/(sqr(q*qbar-q*q)*(p1*(q1+q)-p1*p1-sqrt((p1*p1)*(q*q)))*(p2*(q2+qbar)-p2*p2-sqrt((p2*p2)*(q*q))));
				break;
			}
		case 5:
			{
				/*
				Energy2 p1p2    = p1 * p2;
				Energy2 p1q     = p1 * q ;
				Energy2 p2qbar  = p2 * qbar;
				Energy2 p2q     = p2 * q ;
				Energy2 p1qbar  = p1 * qbar;
				Energy2 qqbar   = q  * qbar;
				Energy2 mq2 = q.mass2();
				double Numerator = p1p2 * (qqbar + mq2)/sqr(GeV2);
				Numerator += 0.5 * (p1q - p1qbar)*(p2q - p2qbar)/sqr(GeV2);
				double Denominator = sqr(qqbar + mq2)*(p1q + p1qbar)*(p2q + p2qbar)/sqr(sqr(GeV2));
				// add test Tchannel Matrix Element interference with gluon mass regulated
				// Denominator *= sqr(sqr(p1-q1)-_epsilonResolution*sqrt((p1*p1)*(p2*p2)))/sqr(GeV2);
				*/
				static Energy mg=getParticleData(spectrum()->gluonId())->constituentMass();
				double Denominator = (sqr(p1-q1)-_epsilonResolution*mg*mg)/(GeV2);
				Denominator *= (sqr(p2-q2)-_epsilonResolution*mg*mg)/(GeV2);
				double Numerator = ((q1*q2)*(p1*p2)+(p2*q1)*(p1*q2))/sqr(GeV2);
				// double I11 = SoftFunction(p1,p1,q,qbar);
				// double I22 = SoftFunction(p2,p2,q,qbar);
				// double I12 = SoftFunction(p1,p2,q,qbar);
				double I11 = SoftFunction(q1,q1,q,qbar);
				double I22 = SoftFunction(q2,q2,q,qbar);
				double I12 = SoftFunction(q1,q2,q,qbar);
				SQME = Numerator*(I11+I22-2.0*I12)/Denominator;
				if (SQME<0.0) {
					std::cout << "I11 = " << I11<< std::endl;
					std::cout << "I22 = " << I22<< std::endl;
					std::cout << "2*I12 = " << I12<< std::endl;
					std::cout << "soft = " << I11+I22-2.0*I12<< std::endl;
					std::cout << "Numerator = " << Numerator<< std::endl;
					std::cout << "Denominator = " << Denominator<< std::endl;
					std::cout << "M  = " << (p1+p2).m()/GeV<< std::endl;
					std::cout << "M1 = " << (q1+q).m()/GeV<< std::endl;
					std::cout << "M2 = " << (q2+qbar).m()/GeV<< std::endl;
					std::cout << "m1 = " << (q1).m()/GeV<< std::endl;
					std::cout << "m2 = " << (q2).m()/GeV<< std::endl;
					std::cout << "m  = " << (q).m()/GeV<< std::endl;
				}
				break;
			}
		case 6:
			{
				Energy M = (p1+p2).m();
				Energy M1 = (q1+q).m();
				Energy M2 = (q2+qbar).m();
				Energy m1 = q1.mass();
				Energy m2 = q2.mass();
				Energy m  = q.mass();
				Energy Pcom = Kinematics::pstarTwoBodyDecay(M,m1,m2);
				Axis z(0,0,1);
				Lorentz5Momentum p1_ali = Lorentz5Momentum(m1, Pcom*z);
				Lorentz5Momentum p2_ali = Lorentz5Momentum(m2,-Pcom*z);
				Lorentz5Momentum Q1_ali = Lorentz5Momentum(M1,GeV*Axis(0,0,0));
				Lorentz5Momentum Q2_ali = Lorentz5Momentum(M2,GeV*Axis(0,0,0));
				Lorentz5Momentum q1_ali = Lorentz5Momentum(m1,GeV*Axis(0,0,0));
				Lorentz5Momentum q2_ali = Lorentz5Momentum(m2,GeV*Axis(0,0,0));
				Lorentz5Momentum q_ali = Lorentz5Momentum(m,GeV*Axis(0,0,0));
				Lorentz5Momentum qbar_ali = Lorentz5Momentum(m,GeV*Axis(0,0,0));
				Lorentz5Momentum pClutmp(p1_ali+p2_ali);
				Kinematics::twoBodyDecay(pClutmp, M1, M2, z, Q1_ali, Q2_ali);
				Axis direction1 = Q1_ali.vect().unit();
				// Need to boost constituents first into the pClu rest frame
				//  boost from Cluster1 rest frame to Cluster COM Frame
				Kinematics::twoBodyDecay(Q1_ali, m1, m, direction1, q1_ali, q_ali);
				Axis direction2 = Q2_ali.vect().unit();
				//  boost from Cluster2 rest frame to Cluster COM Frame 
				Kinematics::twoBodyDecay(Q2_ali, m2, m, direction2, q2_ali, qbar_ali);

				if (areNotSame(p1_ali+p2_ali, q1_ali+q2_ali+q_ali+qbar_ali))
						std::cout << "ERROR TOT " << std::endl;
				if (areNotSame(Q1_ali, q1_ali+q_ali))
						std::cout << "ERROR Q1 " << std::endl;
				if (areNotSame(Q2_ali, q2_ali+qbar_ali))
						std::cout << "ERROR Q2 " << std::endl;
				static Energy mg=getParticleData(spectrum()->gluonId())->constituentMass();
				double Denominator = (sqr(p1_ali-q1_ali)-_epsilonResolution*mg*mg)/(GeV2);
				Denominator *= (sqr(p2_ali-q2_ali)-_epsilonResolution*mg*mg)/(GeV2);
				double Numerator = ((q1_ali*q2_ali)*(p1_ali*p2_ali)+(p2_ali*q1_ali)*(p1_ali*q2_ali))/sqr(GeV2);
				// double I11 = SoftFunction(p1_ali,p1_ali,q_ali,qbar_ali);
				// double I22 = SoftFunction(p2_ali,p2_ali,q_ali,qbar_ali);
				// double I12 = SoftFunction(p1_ali,p2_ali,q_ali,qbar_ali);
				double I11 = SoftFunction(q1_ali,q1_ali,q_ali,qbar_ali);
				double I22 = SoftFunction(q2_ali,q2_ali,q_ali,qbar_ali);
				double I12 = SoftFunction(q1_ali,q2_ali,q_ali,qbar_ali);
				SQME = Numerator*(I11+I22-2.0*I12)/Denominator;
				if (SQME<0.0) {
					std::cout << "I11 = " << I11<< std::endl;
					std::cout << "I22 = " << I22<< std::endl;
					std::cout << "2*I12 = " << I12<< std::endl;
					std::cout << "soft = " << I11+I22-2.0*I12<< std::endl;
					std::cout << "Numerator = " << Numerator<< std::endl;
					std::cout << "Denominator = " << Denominator<< std::endl;
					std::cout << "M  = " << (p1_ali+p2_ali).m()/GeV<< std::endl;
					std::cout << "M1 = " << (q1_ali+q_ali).m()/GeV<< std::endl;
					std::cout << "M2 = " << (q2_ali+qbar_ali).m()/GeV<< std::endl;
					std::cout << "m1 = " << (q1_ali).m()/GeV<< std::endl;
					std::cout << "m2 = " << (q2_ali).m()/GeV<< std::endl;
					std::cout << "m  = " << (q_ali).m()/GeV<< std::endl;
				}
				break;
			}
		default:
			assert(false);
	}
	if (SQME < 0) throw	Exception()
		<< "Squared Matrix Element = "<< SQME <<" < 0 in MatrixElementClusterFissioner::calculateSQME() "
		<< Exception::runerror;
	return SQME;
}
/* Overestimate for SQME for p1,p2->C1(q1,q),C2(q2,qbar) 
 * Note that:
 *      p0Q1  -> p1   where Mass -> m1
 *      p0Q2  -> p2   where Mass -> m2
 *      pQ1   -> q1   where Mass -> m1
 *      pQone -> q    where Mass -> mq
 *      pQ2   -> q2   where Mass -> m2
 *      pQone -> qbar where Mass -> mq
 * */
double MatrixElementClusterFissioner::calculateSQME_OverEstimate(
		const Energy& Mc,
		const Energy& m1,
		const Energy& m2,
		const Energy& mq
		) const {
	double SQME_OverEstimate;
	switch (_matrixElement)
	{
		case 0:
				SQME_OverEstimate = 1.0;
				break;
		case 1:
			{
				// Fit factor for guess of best overestimate
				double A = 0.25;
				SQME_OverEstimate = _safetyFactorMatrixElement*A*pow(mq/GeV,-4);
				break;
			}
		case 2:
			{
				// Fit factor for guess of best overestimate
				double A = 0.25;
				SQME_OverEstimate = _safetyFactorMatrixElement*A*pow(mq/GeV,-4);
				break;
			}
		case 4:
			{
				// Fit factor for guess of best overestimate
				SQME_OverEstimate = _safetyFactorMatrixElement*sqr(sqr(GeV2))/(sqr(mq*mq)*m1*m2*(m1+mq)*(m2+mq));
				break;
			}
		case 5:
		case 6:
			{
				// Fit factor for guess of best overestimate
				// Energy MassMax = 8.550986578849006037e+00*GeV; // mass max of python
				// double wTotMax = 5.280507222727160297e+03; // at MassMax
				// double argExp = 55.811; // from fit of python package
				// double powLog = 0.08788; // from fit of python package

				// Old
				static Energy MassMax = 5.498037126562503119e+01*GeV; // mass max of python
				static const double wTotMax = 1.570045806367627112e+06; // at MassMax
				static const double argExp = 16.12; // from fit of python package
				static const double powLog = 0.3122; // from fit of python package
				// New fits from python package (most effective for all masses m_ud)
				// static Energy MassMax = 1.072117239679688225e+02*GeV; // mass max of python
				// static const double wTotMax = 6.572201964467836914e+11; // at MassMax
				// static const double argExp = 1.216961110145260250e+01; // from fit of python package
				// static const double powLog = 6.437455646293682721e-01; // from fit of python package
				// Energy2 p1p2=0.5*(Mc*Mc-m1*m1-m2*m2);
				// SQME_OverEstimate = _safetyFactorMatrixElement*A*pow(mq/GeV,-4)*(2*sqr(p1p2))/sqr(_epsilonResolution*(m1*m2));
				Energy Mmin = (m1+m2+2*mq);
				// Energy MdblPrec = Mmin*exp(pow(pow(log(MassMax/Mmin),powLog)+600.0/argExp,1.0/powLog));
				// if (Mc >MdblPrec)
					// std::cout << "errroRRRRR R MC = " << Mc/GeV << "\t Mcmax "<< MdblPrec/GeV<< std::endl;
				if (MassMax<Mmin) MassMax=Mmin;
				double Overestimate = wTotMax*exp(argExp*(pow(log(Mc/Mmin),powLog)-pow(log(MassMax/Mmin),powLog)));
				SQME_OverEstimate = _safetyFactorMatrixElement*Overestimate;//*A*pow(mq/GeV,-4)*(2*sqr(p1p2))/sqr(_epsilonResolution*(m1*m2));
				if (SQME_OverEstimate==0 || std::isinf(SQME_OverEstimate)|| std::isnan(SQME_OverEstimate)) {
					std::cout << "SQME_OverEstimate is " << SQME_OverEstimate << std::endl;
					std::cout << "Overestimate is " << Overestimate << std::endl;
					std::cout << "MC " << Mc/GeV << std::endl;
					std::cout << "m " << mq/GeV << std::endl;
					std::cout << "m1 " << m1/GeV << std::endl;
					std::cout << "m2 " << m2/GeV << std::endl;
				}
				break;
			}
		default:
			assert(false);
	}
	return SQME_OverEstimate;
}

double MatrixElementClusterFissioner::drawKinematics(
             const Lorentz5Momentum & pClu,
					   const Lorentz5Momentum & p0Q1,
					   const Lorentz5Momentum & p0Q2,
					   const bool toHadron1,
					   const bool toHadron2,
					   Lorentz5Momentum & pClu1,
					   Lorentz5Momentum & pClu2,
					   Lorentz5Momentum & pQ1,
					   Lorentz5Momentum & pQbar,
					   Lorentz5Momentum & pQ,
					   Lorentz5Momentum & pQ2bar) const {
	double weightTotal=0.0;
  if (pClu.mass() < pClu1.mass() + pClu2.mass()
		  || pClu1.mass()<ZERO
		  || pClu2.mass()<ZERO  ) {
    throw Exception() << "Impossible Kinematics in MatrixElementClusterFissioner::drawKinematics() (A)\n"
					<< "Mc  = "<< pClu.mass()/GeV <<" GeV\n"
					<< "Mc1 = "<< pClu1.mass()/GeV <<" GeV\n"
					<< "Mc2 = "<< pClu2.mass()/GeV <<" GeV\n"
		      << Exception::eventerror;
  }

	// Sample direction of the daughter clusters
	std::pair<Axis,double> directionCluster = sampleDirectionCluster(p0Q1, pClu);
	Axis DirToClu = directionCluster.first;
	weightTotal   = directionCluster.second;
	if (0 && _writeOut) {
		ofstream out("WriteOut/test_Cluster.dat",std::ios::app);
		Lorentz5Momentum pQinCOM(p0Q1);
		pQinCOM.setMass(p0Q1.m());
		pQinCOM.boost( -pClu.boostVector() );        // boost from LAB to C
		out << DirToClu.cosTheta(pQinCOM.vect()) << "\n";
	}
	if (_covariantBoost) {
		const Energy M  = pClu.mass();
		const Energy M1 = pClu1.mass();
		const Energy M2 = pClu2.mass();
		const Energy PcomClu=Kinematics::pstarTwoBodyDecay(M,M1,M2);
		Momentum3 pClu1sampVect( PcomClu*DirToClu);
		Momentum3 pClu2sampVect(-PcomClu*DirToClu);
		pClu1.setMass(M1);
		pClu1.setVect(pClu1sampVect);
		pClu1.rescaleEnergy();
		pClu2.setMass(M2);
		pClu2.setVect(pClu2sampVect);
		pClu2.rescaleEnergy();
	}
	else {
		Lorentz5Momentum pClutmp(pClu);
		pClutmp.boost(-pClu.boostVector());
		Kinematics::twoBodyDecay(pClutmp, pClu1.mass(), pClu2.mass(),DirToClu, pClu1, pClu2);
	}
	// In the case that cluster1 does not decay immediately into a single hadron,
	// calculate the momenta of Q1 (as constituent of C1) and Qbar in the
	// (parent) C1 frame first, where the direction of Q1 is u.vect().unit(),
	// and then boost back in the LAB frame.
	if(!toHadron1) {
		if (pClu1.m() < pQ1.mass() + pQbar.mass() ) {
			throw Exception() << "Impossible Kinematics in MatrixElementClusterFissioner::drawKinematics() (B)"
				<< Exception::eventerror;
		}
		std::pair<Axis,double> direction1 = sampleDirectionConstituents(pClu1,pClu.mass());
		// Need to boost constituents first into the pClu rest frame
		//  boost from Cluster1 rest frame to Cluster COM Frame
		// Kinematics::twoBodyDecay(pClu1, pQ1.mass(), pQbar.mass(), DirClu1, pQ1, pQbar);
		Kinematics::twoBodyDecay(pClu1, pQ1.mass(), pQbar.mass(), direction1.first, pQ1, pQbar);
		weightTotal*=direction1.second;
		if (0 && _writeOut) {
			ofstream out("WriteOut/test_Constituents1.dat",std::ios::app);
			Lorentz5Momentum pQ1inCOM(pQ1);
			pQ1inCOM.setMass(pQ1.m());
			pQ1inCOM.boost( -pClu1.boostVector() );        // boost from LAB to C
			out << pQ1inCOM.vect().cosTheta(pClu1.vect()) << "\n";
		}
	}

	// In the case that cluster2 does not decay immediately into a single hadron,
	// Calculate the momenta of Q and Q2bar (as constituent of C2) in the
	// (parent) C2 frame first, where the direction of Q is u.vect().unit(),
	// and then boost back in the LAB frame.
	if(!toHadron2) {
		if (pClu2.m() < pQ.mass() + pQ2bar.mass() ) {
			throw Exception() << "Impossible Kinematics in MatrixElementClusterFissioner::drawKinematics() (C)"
				<< Exception::eventerror;
		}
		std::pair<Axis,double> direction2 = sampleDirectionConstituents(pClu2,pClu.mass());
		//  boost from Cluster2 rest frame to Cluster COM Frame 
		Kinematics::twoBodyDecay(pClu2, pQ2bar.mass(), pQ.mass(), direction2.first, pQ2bar, pQ);
		weightTotal*=direction2.second;
		if (0 && _writeOut) {
			ofstream out("WriteOut/test_Constituents2.dat",std::ios::app);
			Lorentz5Momentum pQ2inCOM(pQ2bar);
			pQ2inCOM.setMass(pQ2bar.m());
			pQ2inCOM.boost( -pClu2.boostVector() );        // boost from LAB to C
			out << pQ2inCOM.vect().cosTheta(pClu2.vect()) << "\n";
		}
	}
	// Boost all momenta from the Cluster COM frame to the Lab frame
	if (_covariantBoost) {
			std::vector<Lorentz5Momentum *> momenta;
			momenta.push_back(&pClu1);
			momenta.push_back(&pClu2);
			if (!toHadron1) {
				momenta.push_back(&pQ1);
				momenta.push_back(&pQbar);
			}
			if (!toHadron2) {
				momenta.push_back(&pQ);
				momenta.push_back(&pQ2bar);
			}
			Kinematics::BoostIntoTwoParticleFrame(pClu.mass(),p0Q1, p0Q2, momenta);
	}
	else {
		pClu1.boost(pClu.boostVector());
		pClu2.boost(pClu.boostVector());
		pQ1.boost(pClu.boostVector());
		pQ2bar.boost(pClu.boostVector());
		pQ.boost(pClu.boostVector());
		pQbar.boost(pClu.boostVector());
	}
	return weightTotal; // success
}

