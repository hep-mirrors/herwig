// -*- C++ -*-
//
// PowhegInclusiveME.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PowhegInclusiveME_H
#define HERWIG_PowhegInclusiveME_H
//
// This is the declaration of the PowhegInclusiveME class.
//

#include "ThePEG/MatrixElement/MEGroup.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxNLOME.h"
#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"
#include "Herwig++/MatrixElement/Matchbox/Base/SubtractedME.h"
#include "Herwig++/MatrixElement/Matchbox/Powheg/PowhegSplittingKernel.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief PowhegInclusiveME represents a BBar function.
 *
 */
class PowhegInclusiveME: public MEGroup {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  PowhegInclusiveME();

  /**
   * The destructor.
   */
  virtual ~PowhegInclusiveME();
  //@}

public:

  /**
   * Setup from NLO ME and Subtracted ME's
   */
  void setup(Ptr<MatchboxNLOME>::ptr,
	     const vector<Ptr<SubtractedME>::ptr>&,
	     bool bornScreening);

  /**
   * Return the dipoles used
   */
  vector<Ptr<SubtractionDipole>::ptr> dipoles() const;

  /**
   * Return the splitting kernels
   */
  const vector<Ptr<PowhegSplittingKernel>::ptr>& splittingKernels() const {
    return theSplittingKernels;
  }

  /**
   * Access the splitting kernels
   */
  vector<Ptr<PowhegSplittingKernel>::ptr>& splittingKernels() {
    return theSplittingKernels;
  }

  /**
   * Return true, if SubProcessGroups should be
   * setup from this MEGroup. If not, a single SubProcess
   * is constructed from the data provided by the
   * head matrix element.
   */
  virtual bool subProcessGroups() const { return false; }

public:

  /**
   * Return true, if the same additional random numbers
   * should be presented to any of the dependent
   * matrix elements.
   */
  virtual bool uniformAdditional() const { return true; }

  /**
   * Given a process from the head matrix element,
   * return a list of diagrams which should be considered for
   * the given dependent matrix element.
   */
  virtual MEBase::DiagramVector dependentDiagrams(const cPDVector& proc,
						  tMEPtr depME) const;

  /**
   * Set the XComb object to be used in the next call to
   * generateKinematics() and dSigHatDR().
   */
  virtual void setXComb(tStdXCombPtr);

  /**
   * Return true for MC summation of dependent
   * matrix elements, if feasible.
   */
  virtual bool mcSumDependent() const { return theMCSum; }

  /**
   * Switch on MC summation of dependent
   * matrix elements, if feasible.
   */
  virtual void doMCSum() { theMCSum = true; }

public:

  /**
   * Dump the setup to an ostream
   */
  void print(ostream&) const;

  /**
   * Collect information on the last evaluated phasespace
   * point for verification or debugging purposes. This
   * only called, if the StdXCombGroup did accumulate
   * a non-zero cross section from this ME group.
   */
  virtual void lastEventStatistics();

  /**
   * Print debug information on the last event
   */
  void printLastEvent(ostream&) const;

  /**
   * Switch on verbosity for this subtracted ME
   */
  void beVerbose() { theVerbose = true; }

  /**
   * Switch off verbosity for this subtracted ME
   */
  void beQuiet() { theVerbose = false; }

  /**
   * Return true, if verbose
   */
  bool verbose() const { return theVerbose; }

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

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}


private:

  /**
   * The splitting kernels which need to be considered
   */
  vector<Ptr<PowhegSplittingKernel>::ptr> theSplittingKernels;

  /**
   * Map dependent me's to splitting kernels to
   * properly setup the xcombs
   */
  map<Ptr<MEBase>::ptr,Ptr<PowhegSplittingKernel>::ptr> theKernelMap;

  /**
   * Switch to print full information on the
   * last phase space point.
   */
  bool theVerbose;

  /**
   * Wether or not the real contributions should be MC summed
   */
  bool theMCSum;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PowhegInclusiveME & operator=(const PowhegInclusiveME &);

};

}

#endif /* HERWIG_PowhegInclusiveME_H */
