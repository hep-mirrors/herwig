// -*- C++ -*-
//
// SubtractedME.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SubtractedME_H
#define HERWIG_SubtractedME_H
//
// This is the declaration of the SubtractedME class.
//

#include "ThePEG/MatrixElement/MEGroup.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"
#include "Herwig/MatrixElement/Matchbox/Utility/LastMatchboxXCombInfo.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.fh"

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
class SubtractedME: 
    public MEGroup, 
    public LastMatchboxXCombInfo {

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

  /**
   * Return the factory which produced this matrix element
   */
  Ptr<MatchboxFactory>::tcptr factory() const;

  /**
   * Set the factory which produced this matrix element
   */
  void factory(Ptr<MatchboxFactory>::tcptr f);

  /** @name Phasespace and subprocess information */
  //@{

  /**
   * For the given event generation setup return an xcomb object
   * appropriate to this matrix element.
   */
  virtual StdXCombPtr makeXComb(Energy newMaxEnergy, const cPDPair & inc,
				tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
				tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
				const PBPair & newPartonBins, tCutsPtr newCuts,
				const DiagramVector & newDiagrams, bool mir,
				const PartonPairVec& allPBins,
				tStdXCombPtr newHead = tStdXCombPtr(),
				tMEPtr newME = tMEPtr());

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
   * Return true, if the XComb steering this matrix element
   * should keep track of the random numbers used to generate
   * the last phase space point
   */
  virtual bool keepRandomNumbers() const { return true; }

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
  virtual bool subProcessGroups() const;

  /**
   * Return true, if one of the dependent subprocesses should be
   * constructed in place of the one driven by the head matrix element
   * or a full subprocess group.
   */
  virtual bool selectDependentSubProcess() const { return false; }

  /**
   * Fill the projectors object of xcombs to choose subprocesses
   * different than the one currently integrated.
   */
  virtual void fillProjectors();

  /**
   * Return true, if projectors will be used
   */
  virtual bool willProject() const { 
    return virtualShowerSubtraction() || loopSimSubtraction();
  }

  /**
   * Return true, if this MEGroup will reweight the contributing cross
   * sections.
   */
  virtual bool groupReweighted() const { 
    return showerApproximation();
  }

  /**
   * Reweight the head cross section
   */
  virtual double reweightHead(const vector<tStdXCombPtr>&);

  /**
   * Reweight the dependent cross section
   */
  virtual double reweightDependent(tStdXCombPtr, const vector<tStdXCombPtr>&);

  /**
   * Switch on or off that scales should be calculated from real emission kinematics
   */
  void doRealEmissionScales();

  //@}

  /** @name Methods relevant to matching */
  //@{

  /**
   * Inform this matrix element that a new phase space
   * point is about to be generated, so all caches should
   * be flushed.
   */
  virtual void flushCaches() { 
    MEGroup::flushCaches();
    if ( showerApproximation() )
      showerApproximation()->resetBelowCutoff();
  }

  /**
   * Return the shower approximation.
   */
  Ptr<ShowerApproximation>::tptr showerApproximation() const;

  /**
   * Indicate that the shower real emission contribution should be subtracted.
   */
  void doRealShowerSubtraction();

  /**
   * Return true, if the shower real emission contribution should be subtracted.
   */
  bool realShowerSubtraction() const { return theRealShowerSubtraction; }

  /**
   * Indicate that the shower virtual contribution should be subtracted.
   */
  void doVirtualShowerSubtraction();

  /**
   * Return true, if the shower virtual contribution should be subtracted.
   */
  bool virtualShowerSubtraction() const { return theVirtualShowerSubtraction; }

  /**
   * Indicate that the loopsim matched virtual contribution should be subtracted.
   */
  void doLoopSimSubtraction();

  /**
   * Return true, if the loopsim matched virtual contribution should be subtracted.
   */
  bool loopSimSubtraction() const { return theLoopSimSubtraction; }

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
  const vector<Ptr<MatchboxMEBase>::ptr>& borns() const;

  /**
   * Access the underlying born matrix elements,
   * overriding the ones contained in the factory object.
   */
  void setBorns(const vector<Ptr<MatchboxMEBase>::ptr>& newBorns) { theBorns = newBorns; }

  /**
   * Build up dipoles needed.
   */
  void getDipoles();

  /**
   * Clone all dipoles.
   */
  void cloneDipoles(const string& prefix = "");

  /**
   * Clone the real emission matrix element.
   */
  void cloneRealME(const string& prefix = "");

  /**
   * Clone all dependencies.
   */
  void cloneDependencies(const string& prefix = "") {
    cloneDipoles(prefix);
    cloneRealME(prefix);
  }

  /**
   * Return all dipoles matching the given Born process
   */
  vector<Ptr<SubtractionDipole>::ptr> splitDipoles(const cPDVector&);

  //@}

  /** @name Veto scale settings */
  //@{
  /**
   * Set veto scales on the particles at the given
   * SubProcess which has been generated using this
   * matrix element.
   */
  virtual void setVetoScales(tSubProPtr) const;
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
   * Return true, if verbose
   */
  bool verbose() const;

  /**
   * Return true, if verbose
   */
  bool initVerbose() const;

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
    SubtractionHistogram(double low = 0.001, 
			 double up = 10., 
			 unsigned int nbins = 100);

    /**
     * Book an event.
     */
    void book(double inv, double diff) {
      map<double,pair<double,double> >::iterator b =
	bins.upper_bound(inv);
      if ( b == bins.end() ) return;
      b->second.first = min(b->second.first,diff);
      b->second.second = max(b->second.second,diff);
    }

    /**
     * Write to file given name and invariant.
     */
    void dump(const std::string& prefix, 
        const int& plottype,
        const bool& scatterplot,
	      const cPDVector& proc,
	      int i, int j) const;

    /**
     * Write to persistent ostream
     */
    void persistentOutput(PersistentOStream&) const;

    /**
     * Read from persistent istream
     */
    void persistentInput(PersistentIStream&);

  };

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

  //@}

private:

  /**
   * The factory which produced this matrix element
   */
  Ptr<MatchboxFactory>::tcptr theFactory;

  /**
   * The underlying born matrix elements, overriding the ones
   * contained in the factory object.
   */
  vector<Ptr<MatchboxMEBase>::ptr> theBorns;

  /**
   * Pointer to the head real emission ME casted to a MatchboxMEBase
   * object.
   */
  Ptr<MatchboxMEBase>::ptr theReal;

  /**
   * Define the key for the collinear subtraction data.
   */
  typedef pair<cPDVector,pair<size_t, size_t> > CollinearSubtractionIndex;

  /**
   * subtraction data for collinear limits.
   */
  map<CollinearSubtractionIndex,SubtractionHistogram> collinearHistograms;

  /**
   * names of files to which subtraction data is written for all phase space points individually
   */
  map<CollinearSubtractionIndex,string> fnamesCollinearSubtraction;

  /**
   * Define the key for the soft subtraction data.
   */
  typedef pair<cPDVector,size_t> SoftSubtractionIndex;

  /**
   * subtraction data for soft limits.
   */
  map<SoftSubtractionIndex,SubtractionHistogram> softHistograms;

  /**
   * names of files to which subtraction data is written for all phase space points individually
   */
  map<SoftSubtractionIndex,string> fnamesSoftSubtraction;

  /**
   * True, if the shower real emission contribution should be subtracted.
   */
  bool theRealShowerSubtraction;

  /**
   * True, if the shower virtual contribution should be subtracted.
   */
  bool theVirtualShowerSubtraction;

  /**
   * True, if the loopsim matched virtual contribution should be subtracted.
   */
  bool theLoopSimSubtraction;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SubtractedME & operator=(const SubtractedME &);

};

inline PersistentOStream& operator<<(PersistentOStream& os,
				     const SubtractedME::SubtractionHistogram& h) {
  h.persistentOutput(os);
  return os;
}

inline PersistentIStream& operator>>(PersistentIStream& is,
				     SubtractedME::SubtractionHistogram& h) {
  h.persistentInput(is);
  return is;
}

}

#endif /* HERWIG_SubtractedME_H */
