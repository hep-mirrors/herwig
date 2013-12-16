// -*- C++ -*-
//
// HEJMEBase.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_HEJMEBase_H
#define Herwig_HEJMEBase_H
//
// This is the declaration of the HEJMEBase class.
//

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"
#include "HEJFactory.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Keys for XComb meta information
 */
struct HEJMetaKeys {

  enum Keys {
    Jets
  };

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Convert CLHEP momenta
 */
struct CLHEPConverter {

  Lorentz5Momentum operator()(const CLHEP::HepLorentzVector& p) const {
    Lorentz5Momentum ret(p.x()*GeV,p.y()*GeV,p.z()*GeV,p.t()*GeV);
    ret.rescaleMass();
    return ret;
  }

  void operator()(const CLHEP::HepLorentzVector& p, Lorentz5Momentum& q) const {
    q = Lorentz5Momentum(p.x()*GeV,p.y()*GeV,p.z()*GeV,p.t()*GeV);
    q.rescaleMass();
  }

  CLHEP::HepLorentzVector operator()(const Lorentz5Momentum& p) const {
    CLHEP::HepLorentzVector ret(p.x()/GeV,p.y()/GeV,p.z()/GeV,p.t()/GeV);
    return ret;
  }

  void operator()(const Lorentz5Momentum& p, CLHEP::HepLorentzVector& q) const {
    q = CLHEP::HepLorentzVector(p.x()/GeV,p.y()/GeV,p.z()/GeV,p.t()/GeV);
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief HEJMEBase is the base class for HEJ matrix elements.
 *
 * @see \ref HEJMEBaseInterfaces "The interfaces"
 * defined for HEJMEBase.
 */
class HEJMEBase: public Herwig::MatchboxMEBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  HEJMEBase();

  /**
   * The destructor.
   */
  virtual ~HEJMEBase();
  //@}

  /**
   * Clone this matrix element.
   */
  Ptr<HEJMEBase>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<HEJMEBase>::ptr>(clone());
  }

public:

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   */
  virtual double me2() const;

public:

  /**
   * Get the flavour of the first incoming parton
   */
  cPDPtr firstIncoming() const { return theFirstIncoming; }

  /**
   * Get the flavour of the second incoming parton
   */
  cPDPtr secondIncoming() const { return theSecondIncoming; }

  /**
   * Get the number of additional gluons between the outermost partons.
   */
  unsigned int nGluons() const { return theNGluons; }

  /**
   * Get the jets object pointer
   */
  CMultijet* jets() const { return theJets; }

  /**
   * Set the flavour of the first incoming parton
   */
  void firstIncoming(cPDPtr d) { theFirstIncoming = d; }

  /**
   * Set the flavour of the second incoming parton
   */
  void secondIncoming(cPDPtr d) { theSecondIncoming = d; }

  /**
   * Set the number of additional gluons between the outermost partons.
   */
  void nGluons(unsigned int n) { theNGluons = n; }

  /**
   * Set the jets object pointer
   */
  void jets(CMultijet* j) { theJets = j; }

  /**
   * Set the factory pointer
   */
  void factory(Ptr<HEJFactory>::ptr);

public:

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;

  /**
   * With the information previously supplied with the
   * setKinematics(...) method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector &) const;
  using MatchboxMEBase::diagrams;

  /**
   * Return the number of light flavours, this matrix
   * element is calculated for.
   */
  virtual unsigned int nLight() const { return 4; }

  /**
   * Return the renormalization scale for the last generated phasespace point.
   */
  virtual Energy2 factorizationScale() const;

  /**
   * Return the renormalization scale for the last generated phasespace point. 
   */
  virtual Energy2 renormalizationScale() const;

public:

  /**
   * Set the XComb object to be used in the next call to
   * generateKinematics() and dSigHatDR().
   */
  virtual void setXComb(tStdXCombPtr xc);

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
   * The flavour of the first incoming parton
   */
  cPDPtr theFirstIncoming;

  /**
   * The flavour of the second incoming parton
   */
  cPDPtr theSecondIncoming;

  /**
   * The number of additional gluons between the outermost partons.
   */
  unsigned int theNGluons;

  /**
   * A pointer to the HEJFactory object which created this matrix element.
   */
  Ptr<HEJFactory>::ptr theFactory;

  /**
   * A pointer to the CMultijet object
   */
  CMultijet* theJets;

  /**
   * The current colour lines used
   */
  mutable ColourLines currentColourLines;

  /**
   * First incoming momentum
   */
  mutable CLHEP::HepLorentzVector pa;

  /**
   * Second incoming momentum
   */
  mutable CLHEP::HepLorentzVector pb;

  /**
   * Outgoing momenta
   */
  mutable vector<CLHEP::HepLorentzVector> pOut;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HEJMEBase & operator=(const HEJMEBase &);

};

}

#endif /* Herwig_HEJMEBase_H */
