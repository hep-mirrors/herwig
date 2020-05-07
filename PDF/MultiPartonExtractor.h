// -*- C++ -*-
#ifndef Herwig_MultiPartonExtractor_H
#define Herwig_MultiPartonExtractor_H
//
// This is the declaration of the MultiPartonExtractor class.
//

#include "ThePEG/PDF/PartonExtractor.h"
#include <deque>

namespace Herwig {

using namespace ThePEG;

/**
 * The MultiPartonExtractor class inherits from the PartonExtractor of ThePEG
 * but allows more control over the PDFs used in the case that there are multiple
 * stages of parton extraction
 *
 * @see \ref MultiPartonExtractorInterfaces "The interfaces"
 * defined for MultiPartonExtractor.
 */
class MultiPartonExtractor: public PartonExtractor {

public:

  /**
   * The default constructor.
   */
  MultiPartonExtractor() {};

  /**
   * Return a vector of possible pairs of parton bins which can be
   * produced within a given maximum total particle-particle
   * invariant mass squared, \a maxEnergy sBin.
   */
  virtual PartonPairVec getPartons(Energy maxEnergy, const cPDPair &,
				   const Cuts &) const;

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

  /**
   * Add parton bins to pbins for the given incoming particle and the
   * specified cuts.
   */
  virtual void addPartons(tPBPtr incoming ,const PDFCuts & cuts,
			  std::deque<tcPDFPtr> pdf ,PartonVector & pbins) const;
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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MultiPartonExtractor & operator=(const MultiPartonExtractor &) = delete;


  /**
   *  PDFBase object to override  first PDF
   */
  vector<PDFPtr> firstPDF_;

  /**
   *  PDFBase object to override second PDF
   */
  vector<PDFPtr> secondPDF_;
  
};

}

#endif /* Herwig_MultiPartonExtractor_H */
