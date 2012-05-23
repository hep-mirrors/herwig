// -*- C++ -*-
//
// BasicLesHouchesFileReader.h is a part of Herwig++ - A multi-purpose
// Monte Carlo event generator.
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BasicLesHouchesFileReader class.
//
#ifndef HERWIG_BasicLesHouchesFileReader_H
#define HERWIG_BasicLesHouchesFileReader_H
// This is the declaration of the BasicLesHouchesFileReader class.

#include "ThePEG/LesHouches/LesHouchesReader.h"
#include "BasicLesHouchesFileReader.fh"
#include "ThePEG/PDT/Decayer.h"
#include "ThePEG/Utilities/CFileLineReader.h"
#include <stdio.h>

namespace Herwig {

using namespace ThePEG;

/**
 * BasicLesHouchesFileReader derives from the LesHouchesReader base class
 * to be used for objects which read event files from matrix element
 * generators. It extends LesHouchesReader by defining a file handle to be
 * read from, which is opened and closed by the open() and close()
 * functions. Note that the file handle is a standard C filehandle and
 * not a C++ stream. This is because there is no standard way in C++
 * to connect a pipe to a stream for reading eg. gzipped files. This
 * class is able to read plain event files conforming to the Les
 * Houches Event File accord.
 *
 * @see \ref BasicLesHouchesFileReaderInterfaces "The interfaces"
 * defined for BasicLesHouchesFileReader.
 * @see LesHouchesReader
 * @see BasicLesHouchesFileReader
 */
class BasicLesHouchesFileReader: public LesHouchesReader {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  BasicLesHouchesFileReader() : neve(0), ieve(0) {}

  /**
   * Copy-constructor. Note that a file which is opened in the object
   * copied from will have to be reopened in this.
   */
  BasicLesHouchesFileReader(const BasicLesHouchesFileReader &);

  /**
   * Destructor.
   */
  virtual ~BasicLesHouchesFileReader();
  //@}

public:

  /** @name Virtual functions specified by the LesHouchesReader base class. */
  //@{
  /**
   * Initialize. This function is called by the LesHouchesEventHandler
   * to which this object is assigned.
   */
  virtual void initialize(LesHouchesEventHandler & eh);

  /**
   * Calls readEvent() or uncacheEvent() to read information into the
   * LesHouches common block variables. This function is called by the
   * LesHouchesEventHandler if this reader has been selectod to
   * produce an event.
   *
   * @return the weight asociated with this event. If negative weights
   * are allowed it should be between -1 and 1, otherwise between 0
   * and 1. If outside these limits the previously estimated maximum
   * is violated. Note that the estimated maximum then should be
   * updated from the outside.
   */
  virtual double getEvent();

  /**
   * Calls doReadEvent() and performs pre-defined reweightings. A
   * sub-class overrides this function it must make sure that the
   * corresponding reweightings are done.
   */
  virtual bool readEvent();

  /**
   * Skip \a n events. Used by LesHouchesEventHandler to make sure
   * that a file is scanned an even number of times in case the events
   * are not ramdomly distributed in the file.
   */
  virtual void skip(long n);

  /**
   * Scan the file or stream to obtain information about cross section
   * weights and particles etc. This function should fill the
   * variables corresponding to the /HEPRUP/ common block. The
   * function returns the number of events scanned.
   */
  virtual long scan();

  /**
   * Open a file with events. Derived classes should overwrite it and
   * first calling it before reading in the run information into the
   * corresponding protected variables.
   */
  virtual void open();

  /**
   * Close the file from which events have been read.
   */
  virtual void close();

  /**
   * Read the next event from the file or stream into the
   * corresponding protected variables. Return false if there is no
   * more events or if this was not a LHF event file.
   */
  virtual bool doReadEvent();
  //@}

  /**
   * Return the name of the file from where to read events.
   */
  string filename() const { return theFileName; }

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
   * Standard Init function used to initialize the interfaces.
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

  /** @name Standard (and non-standard) Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Return true if this object needs to be initialized before all
   * other objects because it needs to extract PDFs from the event file.
   */
  virtual bool preInitialize() const;
  //@

protected:

  /**
   * The wrapper around the C FILE stream from which to read
   */
  CFileLineReader cfile;

protected:

  /**
   * The number of events in this file.
   */
  long neve;

  /**
   * The current event number.
   */
  long ieve;

  /**
   * If the file is a standard Les Houches formatted file (LHF) this
   * is its version number. If empty, this is not a Les Houches
   * formatted file
   */
  string LHFVersion;

  /**
   * If LHF. All lines (since the last open() or readEvent()) outside
   * the header, init and event tags.
   */
  string outsideBlock;

  /**
   * If LHF. All lines from the header block.
   */
  string headerBlock;

  /**
   * If LHF. Additional comments found in the init block.
   */
  string initComments;

  /**
   * If LHF. Map of attributes (name-value pairs) found in the init
   * tag.
   */
  map<string,string> initAttributes;

  /**
   * If LHF. Additional comments found with the last read event.
   */
  string eventComments;

  /**
   * If LHF. Map of attributes (name-value pairs) found in the last
   * event tag.
   */
  map<string,string> eventAttributes;

private:

  /**
   * The name of the file from where to read events.
   */
  string theFileName;

  /**
   * Determines whether events in the LH file are or are not read
   * more than once in order to generate the requested number of
   * events.
   */
  bool overSampling_;

private:

  /**
   * Describe an abstract base class with persistent data.
   */
  static ClassDescription<BasicLesHouchesFileReader> initBasicLesHouchesFileReader;

  /**
   * Private and non-existent assignment operator.
   */
  BasicLesHouchesFileReader & operator=(const BasicLesHouchesFileReader &);

public:

  /** @cond EXCEPTIONCLASSES */
  /** Exception class used by BasicLesHouchesFileReader if reading the file
   *  fails. */
  class LesHouchesFileError: public Exception {};
  /** @endcond */

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of BasicLesHouchesFileReader.
 */
template <>
struct BaseClassTrait<Herwig::BasicLesHouchesFileReader,1> {
  /** Typedef of the base class of BasicLesHouchesFileReader. */
  typedef LesHouchesReader NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * BasicLesHouchesFileReader class and the shared object where it is
 * defined.
 */
template <>
struct ClassTraits<Herwig::BasicLesHouchesFileReader>
  : public ClassTraitsBase<Herwig::BasicLesHouchesFileReader> {
  /**
   * Return the class name.
   */
  static string className() { return "Herwig::BasicLesHouchesFileReader"; }
  /**
   * Return the name of the shared library to be loaded to get access
   * to the BasicLesHouchesFileReader class and every other class it uses
   * (except the base class).
   */
  static string library() { return "BasicLesHouchesFileReader.so"; }

};

/** @endcond */

}

#endif /* HERWIG_BasicLesHouchesFileReader_H */
