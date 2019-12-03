// -*- C++ -*-
//
// MatchboxPhasespace.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxPhasespace_H
#define HERWIG_MatchboxPhasespace_H
//
// This is the declaration of the MatchboxPhasespace class.
//

#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/MatrixElement/Matchbox/Utility/LastMatchboxXCombInfo.h"
#include "Herwig/MatrixElement/Matchbox/Utility/ProcessData.fh"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.fh"
#include "Herwig/MatrixElement/Matchbox/Phasespace/PhasespaceCouplings.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Wrap around a vector of random numbers to behave as a stream
 * of those.
 */
struct StreamingRnd {

  /**
   * The random numbers
   */
  const double* numbers;

  /**
   * The number of random numbers available.
   */
  size_t nRnd;

  /**
   * Default constructor.
   */
  StreamingRnd()
    : numbers(0), nRnd(0) {}

  /**
   * Construct from random numbers.
   */
  explicit StreamingRnd(const double* newNumbers,
			size_t n)
    : numbers(newNumbers), nRnd(n) {}

  /**
   * Return next random number
   */
  inline double operator()() {
    assert(numbers && nRnd > 0);
    const double ret = numbers[0];
    ++numbers; --nRnd;
    return ret;
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxPhasespace defines an abstract interface to a phase
 * space generator.
 *
 */
class MatchboxPhasespace: 
    public HandlerBase, 
    public LastXCombInfo<StandardXComb>,
    public LastMatchboxXCombInfo {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxPhasespace();

  /**
   * The destructor.
   */
  virtual ~MatchboxPhasespace();
  //@}

public:

  /**
   * Set the XComb object steering the Born matrix
   * element this class represents virtual corrections to.
   */
  virtual void setXComb(tStdXCombPtr xc) { 
    theLastXComb = xc;
    lastMatchboxXComb(xc);
  }

  /**
   * Return the factory object
   */
  Ptr<MatchboxFactory>::tcptr factory() const;

  /**
   * Return the process data object
   */
  Ptr<ProcessData>::tptr processData() const;

  /**
   * Generate a phase space point and return its weight.
   */
  virtual double generateKinematics(const double* r,
				    vector<Lorentz5Momentum>& momenta);

  /**
   * Generate a phase space point and return its weight.
   */
  virtual double generateTwoToNKinematics(const double*,
					  vector<Lorentz5Momentum>& momenta) = 0;

  /**
   * Generate a 2 -> 1 phase space point and return its weight.
   */
  virtual double generateTwoToOneKinematics(const double*,
					    vector<Lorentz5Momentum>& momenta);

  /**
   * Return the number of random numbers required to produce a given
   * multiplicity final state.
   */
  virtual int nDim(const cPDVector&) const;

  /**
   * Return the number of random numbers required to produce a given
   * multiplicity final state.
   */
  virtual int nDimPhasespace(int nFinal) const = 0;

  /**
   * Return true, if this phasespace generator will generate incoming
   * partons itself.
   */
  virtual bool haveX1X2() const { return false; }

  /**
   * Return true, if this phase space generator expects
   * the incoming partons in their center-of-mass system
   */
  virtual bool wantCMS() const { return true; }

  /**
   * True, if mass generators should be used instead of fixed masses
   */
  bool useMassGenerators() const { return theUseMassGenerators; }

  /**
   * Fill a diagram selector for the last phase space point.
   */
  virtual Selector<MEBase::DiagramIndex> selectDiagrams(const MEBase::DiagramVector&) const;

  /**
   * Return the momentum and weight appropriate to the given timelike
   * branch of the diagram.
   */
  pair<double,Lorentz5Momentum> timeLikeWeight(const Tree2toNDiagram& diag,
					       int branch, double flatCut) const;

  /**
   * Return the weight appropriate to the given spacelike branch of
   * the diagram.
   */
  double spaceLikeWeight(const Tree2toNDiagram& diag,
			 const Lorentz5Momentum& incoming,
			 int branch, double flatCut) const;

  /**
   * Return the weight appropriate to the given diagram.
   */
  double diagramWeight(const Tree2toNDiagram& diag) const {
    assert( !diagramWeights().empty() );
    return diagramWeights().find(diag.id())->second;
  }

  /**
   * Fill the diagram weights.
   */
  void fillDiagramWeights(double flatCut = 0.0);

  /**
   * Clear the diagram weights.
   */
  void clearDiagramWeights() {
    diagramWeights().clear();
  }

  /**
   * Clone this phase space generator.
   */
  Ptr<MatchboxPhasespace>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<MatchboxPhasespace>::ptr>(clone());
  }

  /**
   * Clone the dependencies, using a given prefix.
   */
  virtual void cloneDependencies(const std::string& prefix = "");

public:

  /**
   * Return true, if this phase space generator is invertible
   */
  virtual bool isInvertible() const { return false; }

  /**
   * Invert the given phase space point to the random numbers which
   * would have generated it.
   */
  virtual double invertKinematics(const vector<Lorentz5Momentum>& momenta,
				  double* r) const;

  /**
   * Invert the given phase space point to the random numbers which
   * would have generated it.
   */
  virtual double invertTwoToNKinematics(const vector<Lorentz5Momentum>&,
					double*) const {
    return 0.;
  }

  /**
   * Invert the given 2 -> 1 phase space point to the random numbers which
   * would have generated it.
   */
  virtual double invertTwoToOneKinematics(const vector<Lorentz5Momentum>&, double*) const;

public:

  /**
   * Limit phasespace generation to a given collinear or soft limit.
   */
  void singularLimit(size_t i, size_t j) {
    if ( i > j )
      swap(i,j);
    singularLimits().insert(make_pair(i,j));
  }

  /**
   * Return the last matched singular limit.
   */
  const pair<size_t,size_t>& lastSingularIndices() const {
    assert(lastSingularLimit() != singularLimits().end());
    return *lastSingularLimit();
  }

  /**
   * Return true, if constraints on phasespace generation have been met.
   */
  bool matchConstraints(const vector<Lorentz5Momentum>& momenta);

protected:

  /**
   * Set a coupling for the given vertex; the convention is that all
   * legs are outgoing, and all possible crossings will be taken care
   * of. If not set, coupling weights default to one.
   */
  void setCoupling(long a, long b, long c,
		   double coupling, bool includeCrossings = true);

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * A cutoff below which a region is considered singular.
   */
  Energy singularCutoff;

  /**
   * True, if mass generators should be used instead of fixed masses
   */
  bool theUseMassGenerators;

  /**
   * Couplings to be used in diagram weighting
   */
  Ptr<PhasespaceCouplings>::ptr theCouplings;

  /**
   * Interface function to setcoupling
   */
  string doSetCoupling(string);

  /**
   * Interface function to setcoupling
   */
  string doSetPhysicalCoupling(string);

  /**
   * The first id in a range of id's meant to denote fictitious
   * 'ghost' particles to be used by the diagram generator
   * in loop induced processes.
   */
  int theLoopParticleIdMin;

  /**
   * The last id in a range of id's meant to denote fictitious
   * 'ghost' particles to be used by the diagram generator
   * in loop induced processes.
   */
  int theLoopParticleIdMax;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxPhasespace & operator=(const MatchboxPhasespace &) = delete;

};

}

#endif /* HERWIG_MatchboxPhasespace_H */
