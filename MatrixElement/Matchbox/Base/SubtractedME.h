// -*- C++ -*-
//
// SubtractedME.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SubtractedME_H
#define HERWIG_SubtractedME_H
//
// This is the declaration of the SubtractedME class.
//

#include "ThePEG/MatrixElement/MEGroup.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig++/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief SubtractedME represents a subtracted real emission matrix element.
 *
 * @see \ref SubtractedMEInterfaces "The interfaces"
 * defined for SubtractedME.
 */
class SubtractedME: public MEGroup {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SubtractedME();

  /**
   * The destructor.
   */
  virtual ~SubtractedME();
  //@}

public:

  /** @name Phasespace and subprocess information */
  //@{
  /**
   * Set the XComb object to be used in the next call to
   * generateKinematics() and dSigHatDR().
   */
  virtual void setXComb(tStdXCombPtr);

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
   * Return true, if SubProcessGroups should be
   * setup from this MEGroup. If not, a single SubProcess
   * is constructed from the data provided by the
   * head matrix element.
   */
  virtual bool subProcessGroups() const { return theSubProcessGroups; }

  /**
   * Switch on or off producing subprocess groups.
   */
  void setSubProcessGroups(bool on = true) { theSubProcessGroups = on; }

  //@}

  /** @name Matrix element and dipole information */
  //@{

  /**
   * Return the subtraction dipoles.
   */
  vector<Ptr<SubtractionDipole>::ptr> dipoles();

  /**
   * Return the underlying born matrix elements.
   */
  const vector<Ptr<MatchboxMEBase>::ptr>& borns() const { return theBorns; }

  /**
   * Access the underlying born matrix elements.
   */
  vector<Ptr<MatchboxMEBase>::ptr>& borns() { return theBorns; }

  /**
   * Build up dipoles needed.
   */
  void getDipoles();

  /**
   * Return all dipoles matching the given Born process
   */
  vector<Ptr<SubtractionDipole>::ptr> splitDipoles(const cPDVector&);

  /**
   * Return the subtraction dipoles.
   */
  const vector<Ptr<SubtractionDipole>::ptr>& allDipoles() const { return theDipoles; }

  /**
   * Access the subtraction dipoles.
   */
  vector<Ptr<SubtractionDipole>::ptr>& allDipoles() { return theDipoles; }

  //@}

  /** @name Veto scale settings */
  //@{
  /**
   * Set veto scales on the particles at the given
   * SubProcess which has been generated using this
   * matrix element.
   */
  virtual void setVetoScales(tSubProPtr) const;

  /**
   * Return true, if veto scales should be set
   * for the real emission
   */
  bool vetoScales() const { return theVetoScales; }

  /**
   * Switch on setting veto scales
   */
  void doVetoScales() { theVetoScales = true; }

  /**
   * Switch off setting veto scales
   */
  void noVetoScales() { theVetoScales = true; }
  //@}

  /** @name Diagnostic information */
  //@{

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
   * Check the subtraction for the last event
   */
  void lastEventSubtraction();

  /**
   * Set a file to print subtraction check to
   */
  void subtractionData(string n) { theSubtractionData = n; }

  /**
   * Return true, if verbose
   */
  bool verbose() const { return theVerbose; }

  /**
   * Switch on or off verbosity for this subtracted ME
   */
  void setVerbose(bool on = true) { theVerbose = on; }

  /**
   * Dump xcomb hierarchies.
   */
  void dumpInfo(const string& prefix = "") const;

  //@}

  /** @name Setup of Subtracted ME objects */
  //@{

  /**
   * Return true if this object needs to be initialized before all
   * other objects (except those for which this function also returns
   * true).  This default version always returns false, but subclasses
   * may override it to return true.
   */
  virtual bool preInitialize() const { return true; }

  //@}

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
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

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
   * The dipoles to be considered; the dipoles generated
   * can be accessed throught the dependent() matrxi element
   * vector, provided the head() is a MatchboxMEBase object.
   */
  vector<Ptr<SubtractionDipole>::ptr> theDipoles;

  /**
   * The underlying Born matrix elements to be considered
   */
  vector<Ptr<MatchboxMEBase>::ptr> theBorns;

  /**
   * File name to dump subtraction check to
   */
  string theSubtractionData;

  /**
   * Simple envelope histogram to keep track of subtraction
   */
  struct SubtractionHistogram {

    /**
     * The lower bound
     */
    double lower;

    /**
     * The bins, indexed by upper bound.
     */
    map<double,pair<double,double> > bins;

    /**
     * Constructor
     */
    SubtractionHistogram(double low = 0.01, 
			 double up = 100., 
			 unsigned int nbins = 100);

    /**
     * Book an event.
     */
    void book(double inv, double ratio) {
      map<double,pair<double,double> >::iterator b =
	bins.upper_bound(inv);
      if ( b == bins.end() ) return;
      b->second.first = min(b->second.first,abs(ratio));
      b->second.second = max(b->second.second,abs(ratio));
    }

    /**
     * Write to file given name and invariant.
     */
    void dump(const std::string& prefix, 
	      const cPDVector& proc,
	      int i, int j) const;

  };

  /**
   * Define the key for the collinear subtraction data.
   */
  typedef pair<cPDVector,pair<size_t, size_t> > CollinearSubtractionIndex;

  /**
   * subtraction data for collinear limits.
   */
  map<CollinearSubtractionIndex,SubtractionHistogram> collinearHistograms;

  /**
   * Define the key for the soft subtraction data.
   */
  typedef pair<cPDVector,size_t> SoftSubtractionIndex;

  /**
   * subtraction data for soft limits.
   */
  map<SoftSubtractionIndex,SubtractionHistogram> softHistograms;

  /**
   * Switch to print full information on the
   * last phase space point.
   */
  bool theVerbose;

  /**
   * True, if SubProcessGroups should be
   * setup from this MEGroup. If not, a single SubProcess
   * is constructed from the data provided by the
   * head matrix element.
   */
  bool theSubProcessGroups;

  /**
   * True, if veto scales should be set
   * for the real emission
   */
  bool theVetoScales;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SubtractedME & operator=(const SubtractedME &);

};

}

#endif /* HERWIG_SubtractedME_H */
