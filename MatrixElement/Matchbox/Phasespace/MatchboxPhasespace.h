// -*- C++ -*-
//
// MatchboxPhasespace.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
    public HandlerBase, public LastXCombInfo<StandardXComb> {

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
   * Prepare a phase space generator for the given xcomb object.
   */
  virtual void prepare(tStdXCombPtr, bool verbose = false) = 0;

  /**
   * Generate a phase space point and return its weight.
   */
  virtual double generateKinematics(const double*,
				    vector<Lorentz5Momentum>& momenta) = 0;

  /**
   * Return the number of random numbers required to produce a given
   * multiplicity final state.
   */
  virtual int nDim(int nFinal) const = 0;

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
    assert( !diagramWeights.empty() );
    return diagramWeights.find(diag.id())->second;
  }

  /**
   * Fill the diagram weights.
   */
  void fillDiagramWeights(double flatCut = 0.0);

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

  /**
   * Dump xcomb hierarchies.
   */
  void dumpInfo(const string& prefix = "") const;

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
   * The diagram weights indexed by diagram id.
   */
  map<int,double> diagramWeights;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxPhasespace & operator=(const MatchboxPhasespace &);

};

}

#endif /* HERWIG_MatchboxPhasespace_H */
