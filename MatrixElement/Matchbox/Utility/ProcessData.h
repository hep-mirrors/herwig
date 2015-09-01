// -*- C++ -*-
//
// ProcessData.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_ProcessData_H
#define Herwig_ProcessData_H
//
// This is the declaration of the ProcessData class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "Herwig/MatrixElement/Matchbox/Utility/ProcessData.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief Provide storage for process data
 *
 * @see \ref ProcessDataInterfaces "The interfaces"
 * defined for ProcessData.
 */
class ProcessData: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ProcessData();

  /**
   * The destructor.
   */
  virtual ~ProcessData();
  //@}

public:

  /**
   * Access diagrams contributing to a given subprocess
   */
  map<PDVector,vector<Ptr<Tree2toNDiagram>::ptr> >& diagramMap() { return theDiagramMap; }

  /**
   * Return diagrams contributing to a given subprocess
   */
  const map<PDVector,vector<Ptr<Tree2toNDiagram>::ptr> >& diagramMap() const { return theDiagramMap; }

  /**
   * Access colour flows contributing to a given diagram
   */
  map<Ptr<Tree2toNDiagram>::tcptr,vector<ColourLines*> >& colourFlowMap() { return theColourFlowMap; }

  /**
   * Return colour flows contributing to a given diagram
   */
  const map<Ptr<Tree2toNDiagram>::tcptr,vector<ColourLines*> >& colourFlowMap() const { return theColourFlowMap; }

  /**
   * Access the mass generators to be used for the given particles
   */
  map<cPDPtr,tGenericMassGeneratorPtr>& massGenerators() { return theMassGenerators; }

  /**
   * Return the mass generators to be used for the given particles
   */
  const map<cPDPtr,tGenericMassGeneratorPtr>& massGenerators() const { return theMassGenerators; }

  /**
   * Return the mass generator (or NULL)
   */
  tGenericMassGeneratorPtr massGenerator(cPDPtr);

  /**
   * Fill the mass generators for the given process
   */
  void fillMassGenerators(const PDVector&);

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

private:

  /**
   * The diagrams contributing to a subprocess
   */
  map<PDVector,vector<Ptr<Tree2toNDiagram>::ptr> > theDiagramMap;

  /**
   * Colour flows contributing to a given diagram
   */
  map<Ptr<Tree2toNDiagram>::tcptr,vector<ColourLines*> > theColourFlowMap;

  /**
   * The mass generators to be used for the given particles
   */
  map<cPDPtr,tGenericMassGeneratorPtr> theMassGenerators;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ProcessData & operator=(const ProcessData &);

};

}

#endif /* Herwig_ProcessData_H */
