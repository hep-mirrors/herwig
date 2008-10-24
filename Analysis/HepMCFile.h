// -*- C++ -*-
//
// HepMCFile.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_HepMCFile_H
#define THEPEG_HepMCFile_H
//
// This is the declaration of the HepMCFile class.
//
#include <iostream>
#include <fstream>
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "HepMC/IO_BaseClass.h"

namespace Herwig {
  using namespace ThePEG;

/**
 * Here is the documentation of the HepMCFile class.
 *
 * @see \ref HepMCFileInterfaces "The interfaces"
 * defined for HepMCFile.
 */
class HepMCFile: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  HepMCFile();

  /**
   * The copy constructor.
   */
  HepMCFile(const HepMCFile &);
  //@}

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);
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

protected:

  /** @name Standard Interfaced functions. */
  //@{
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<HepMCFile> initHepMCFile;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HepMCFile & operator=(const HepMCFile &);

private:

  /**
   * Last event that should be written out as HepMC format
   */
  long _eventNumber;

  /**
   *  The HepMC format
   */
  int _format;

  /**
   *  The HepMC filename
   */
  string _filename;

  /**
   *  The HepMC I/O handler
   */
  HepMC::IO_BaseClass *_hepmcio;

  /**
   *  The HepMC dump file
   */
  ofstream _hepmcdump;

  /**
   * Selector for the choice of units
   */
  int _unitchoice;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HepMCFile. */
template <>
struct BaseClassTrait<Herwig::HepMCFile,1> {
  /** Typedef of the first base class of HepMCFile. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HepMCFile class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HepMCFile>
  : public ClassTraitsBase<Herwig::HepMCFile> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HepMCFile"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the HepMCFile class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so"; }
};

/** @endcond */

}

#endif /* THEPEG_HepMCFile_H */
