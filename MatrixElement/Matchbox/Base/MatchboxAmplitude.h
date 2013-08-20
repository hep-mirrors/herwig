// -*- C++ -*-
//
// MatchboxAmplitude.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxAmplitude_H
#define HERWIG_MatchboxAmplitude_H
//
// This is the declaration of the MatchboxAmplitude class.
//

#include "ThePEG/MatrixElement/Amplitude.h"
#include "ThePEG/Handlers/LastXCombInfo.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/ColourBasis.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/LastMatchboxXCombInfo.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/MatchboxXComb.h"
#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.fh"
#include "Herwig++/MatrixElement/Matchbox/MatchboxFactory.fh"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Process information with coupling order
 */
struct Process {

  PDVector legs;
  unsigned int orderInAlphaS;
  unsigned int orderInAlphaEW;

  Process()
    : orderInAlphaS(0), orderInAlphaEW(0) {}

  Process(const PDVector& p,
	  unsigned int oas,
	  unsigned int oae)
    : legs(p), orderInAlphaS(oas), orderInAlphaEW(oae) {}

  bool operator==(const Process& other) const {
    return
      legs == other.legs &&
      orderInAlphaS == other.orderInAlphaS &&
      orderInAlphaEW == other.orderInAlphaEW;
  }

  bool operator<(const Process& other) const {
    if ( orderInAlphaS != other.orderInAlphaS )
      return orderInAlphaS < other.orderInAlphaS;
    if ( orderInAlphaEW != other.orderInAlphaEW )
      return orderInAlphaEW < other.orderInAlphaEW;
    return legs < other.legs;
  }

  void persistentOutput(PersistentOStream & os) const {
    os << legs << orderInAlphaS << orderInAlphaEW;
  }

  void persistentInput(PersistentIStream & is) {
    is >> legs >> orderInAlphaS >> orderInAlphaEW;
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Enumerate the type of calculation required
 */
namespace ProcessType {

  enum Types {

    treeME2 = 0,
    colourCorrelatedME2,
    spinColourCorrelatedME2,
    oneLoopInterference

  };

}


/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxAmplitude is the base class for amplitude
 * implementations inside Matchbox.
 *
 * @see \ref MatchboxAmplitudeInterfaces "The interfaces"
 * defined for MatchboxAmplitude.
 */
class MatchboxAmplitude: 
    public Amplitude, 
    public LastXCombInfo<StandardXComb>, 
    public LastMatchboxXCombInfo {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxAmplitude();

  /**
   * The destructor.
   */
  virtual ~MatchboxAmplitude();
  //@}

public:

  /**
   * Return the amplitude. Needs to be implemented from
   * ThePEG::Amplitude but is actually ill-defined, as colours of the
   * external particles are not specified. To this extent, this
   * implementation just asserts.
   */
  virtual Complex value(const tcPDVector & particles,
			const vector<Lorentz5Momentum> & momenta, 
			const vector<int> & helicities);

  /** @name Subprocess information */
  //@{

  /**
   * Return true, if this amplitude can handle the given process.
   */
  virtual bool canHandle(const PDVector& p,
			 Ptr<MatchboxFactory>::tptr) const { return canHandle(p); }

  /**
   * Return true, if this amplitude can handle the given process.
   */
  virtual bool canHandle(const PDVector&) const { return false; }

  /**
   * Return the number of random numbers required to evaluate this
   * amplitude at a fixed phase space point.
   */
  virtual int nDimAdditional() const { return 0; }

  /**
   * Return a ME instance appropriate for this amplitude and the given
   * subprocesses
   */
  virtual Ptr<MatchboxMEBase>::ptr makeME(const PDVector&) const;

  /**
   * Set the (tree-level) order in \f$g_S\f$ in which this matrix
   * element should be evaluated.
   */
  virtual void orderInGs(unsigned int) {}

  /**
   * Return the (tree-level) order in \f$g_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGs() const = 0;

  /**
   * Set the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element should be evaluated.
   */
  virtual void orderInGem(unsigned int) {}

  /**
   * Return the (tree-level) order in \f$g_{EM}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInGem() const = 0;

  /**
   * Return the Herwig++ StandardModel object
   */
  Ptr<StandardModel>::tcptr standardModel() { 
    if ( !hwStandardModel() )
      hwStandardModel(dynamic_ptr_cast<Ptr<StandardModel>::tcptr>(HandlerBase::standardModel()));
    return hwStandardModel();
  }

  /**
   * Tell whether the outgoing partons should be sorted when determining
   * allowed subprocesses. Otherwise, all permutations are counted as
   * separate subprocesses.
   */
  virtual bool sortOutgoing() { return true; }

  /**
   * Return true, if this amplitude already includes averaging over
   * incoming parton's quantum numbers.
   */
  virtual bool hasInitialAverage() const { return false; }

  /**
   * Return true, if this amplitude already includes symmetry factors
   * for identical outgoing particles.
   */
  virtual bool hasFinalStateSymmetry() const { return false; }

  /**
   * Return true, if this amplitude is handled by a BLHA one-loop provider
   */
  virtual bool isOLPTree() const { return false; }

  /**
   * Return true, if this amplitude is handled by a BLHA one-loop provider
   */
  virtual bool isOLPLoop() const { return false; }

  /**
   * Write the order file header
   */
  virtual void olpOrderFileHeader(ostream&) const;

  /**
   * Write the order file process list
   */
  virtual void olpOrderFileProcessGroup(ostream&,
					const string&,
					const set<Process>&) const;

  /**
   * Write the order file process list
   */
  virtual void olpOrderFileProcesses(ostream&,
				     const map<pair<Process,int>,int>& procs) const;

  /**
   * Start the one loop provider, if appropriate, giving order and
   * contract files
   */
  virtual void signOLP(const string&, const string&) { }

  /**
   * Start the one loop provider, if appropriate
   */
  virtual void startOLP(const string&, int& status) { status = -1; }

  /**
   * Start the one loop provider, if appropriate. This default
   * implementation writes an BLHA 2.0 order file and starts the OLP
   */
  virtual bool startOLP(const map<pair<Process,int>,int>& procs);


  //@}

  /** @name Colour basis. */
  //@{

  /**
   * Return the colour basis.
   */
  Ptr<ColourBasis>::tptr colourBasis() const { return theColourBasis; }

  /**
   * Return true, if this amplitude will not require colour correlations.
   */
  virtual bool noCorrelations() const { return !haveOneLoop(); }  

  /**
   * Return true, if the colour basis is capable of assigning colour
   * flows.
   */
  virtual bool haveColourFlows() const { 
    return colourBasis() ? colourBasis()->haveColourFlows() : false;
  }

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *> colourGeometries(tcDiagPtr diag) const;

  /**
   * Return an ordering identifier for the current subprocess and
   * colour absis tensor index.
   */
  const string& colourOrderingString(size_t id) const;

  /**
   * Return an ordering identifier for the current subprocess and
   * colour absis tensor index.
   */
  const vector<vector<size_t> >& colourOrdering(size_t id) const;

  //@}

  /** @name Phasespace point, crossing and helicities */
  //@{

  /**
   * Set the xcomb object.
   */
  virtual void setXComb(tStdXCombPtr xc);

  /**
   * Return the momentum as crossed appropriate for this amplitude.
   */
  Lorentz5Momentum amplitudeMomentum(int) const;

  /**
   * Perform a normal ordering of external legs and fill the
   * crossing information as. This default implementation sorts
   * lexicographically in (abs(colour)/spin/abs(charge)), putting pairs
   * of particles/anti-particles where possible.
   */
  virtual void fillCrossingMap(size_t shift = 0);

  /**
   * Generate the helicity combinations.
   */
  virtual set<vector<int> > generateHelicities() const;

  //@}

  /** @name Tree-level amplitudes */
  //@{

  /**
   * Calculate the tree level amplitudes for the phasespace point
   * stored in lastXComb.
   */
  virtual void prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr);

  /**
   * Return the matrix element squared.
   */
  virtual double me2() const;

  /**
   * Return the colour correlated matrix element.
   */
  virtual double colourCorrelatedME2(pair<int,int> ij) const;

  /**
   * Return the large-N colour correlated matrix element.
   */
  virtual double largeNColourCorrelatedME2(pair<int,int> ij,
					   Ptr<ColourBasis>::tptr largeNBasis) const;

  /**
   * Return a positive helicity polarization vector for a gluon of
   * momentum p (with reference vector n) to be used when evaluating
   * spin correlations.
   */
  virtual LorentzVector<Complex> plusPolarization(const Lorentz5Momentum& p,
						  const Lorentz5Momentum& n,
						  int id = -1) const;

  /**
   * Return the colour and spin correlated matrix element.
   */
  virtual double spinColourCorrelatedME2(pair<int,int> emitterSpectator,
					 const SpinCorrelationTensor& c) const;


  /**
   * Return true, if tree-level contributions will be evaluated at amplitude level.
   */
  virtual bool treeAmplitudes() const { return true; }

  /**
   * Evaluate the amplitude for the given colour tensor id and
   * helicity assignment
   */
  virtual Complex evaluate(size_t, const vector<int>&, Complex&) { return 0.; }

  //@}

  /** @name One-loop amplitudes */
  //@{

  /**
   * Return true, if this amplitude is capable of calculating one-loop
   * (QCD) corrections.
   */
  virtual bool haveOneLoop() const { return false; }

  /**
   * Return true, if this amplitude only provides
   * one-loop (QCD) corrections.
   */
  virtual bool onlyOneLoop() const { return false; }

  /**
   * Return true, if one-loop contributions will be evaluated at amplitude level.
   */
  virtual bool oneLoopAmplitudes() const { return true; }

  /**
   * Return true, if one loop corrections have been calculated in
   * dimensional reduction. Otherwise conventional dimensional
   * regularization is assumed. Note that renormalization is always
   * assumed to be MSbar.
   */
  virtual bool isDR() const { return false; }

  /**
   * Return true, if one loop corrections are given in the conventions
   * of the integrated dipoles.
   */
  virtual bool isCS() const { return false; }

  /**
   * Return true, if one loop corrections are given in the conventions
   * of BDK.
   */
  virtual bool isBDK() const { return false; }

  /**
   * Return true, if one loop corrections are given in the conventions
   * of everything expanded.
   */
  virtual bool isExpanded() const { return false; }

  /**
   * Return the value of the dimensional regularization
   * parameter. Note that renormalization scale dependence is fully
   * restored in DipoleIOperator.
   */
  virtual Energy2 mu2() const { return 0.*GeV2; }

  /**
   * If defined, return the coefficient of the pole in epsilon^2
   */
  virtual double oneLoopDoublePole() const { return 0.; }

  /**
   * If defined, return the coefficient of the pole in epsilon
   */
  virtual double oneLoopSinglePole() const { return 0.; }

  /**
   * Calculate the one-loop amplitudes for the phasespace point
   * stored in lastXComb, if provided.
   */
  virtual void prepareOneLoopAmplitudes(Ptr<MatchboxMEBase>::tcptr);

  /**
   * Return the one-loop/tree interference.
   */
  virtual double oneLoopInterference() const;

  /**
   * Evaluate the amplitude for the given colour tensor id and
   * helicity assignment
   */
  virtual Complex evaluateOneLoop(size_t, const vector<int>&) { return 0.; }

  //@}

  /** @name Caching and helpers to setup amplitude objects. */
  //@{

  /**
   * Flush all cashes.
   */
  virtual void flushCaches() {}

  /**
   * Clone this amplitude.
   */
  Ptr<MatchboxAmplitude>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<MatchboxAmplitude>::ptr>(clone());
  }

  /**
   * Clone the dependencies, using a given prefix.
   */
  virtual void cloneDependencies(const std::string& prefix = "");

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

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * Recursively generate helicities
   */
  void doGenerateHelicities(set<vector<int> >& res,
			    vector<int>& current,
			    size_t pos) const;

  /**
   * The colour basis implementation to be used.
   */
  Ptr<ColourBasis>::ptr theColourBasis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxAmplitude & operator=(const MatchboxAmplitude &);

};

inline PersistentOStream& operator<<(PersistentOStream& os,
				     const Process& h) {
  h.persistentOutput(os);
  return os;
}

inline PersistentIStream& operator>>(PersistentIStream& is,
				     Process& h) {
  h.persistentInput(is);
  return is;
}

}

#endif /* HERWIG_MatchboxAmplitude_H */
