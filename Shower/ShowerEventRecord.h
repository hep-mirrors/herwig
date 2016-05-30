// -*- C++ -*-
#ifndef Herwig_ShowerEventRecord_H
#define Herwig_ShowerEventRecord_H
//
// This is the declaration of the ShowerEventRecord class.
//

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDF/PDF.h"
#include "Herwig/MatrixElement/Matchbox/Matching/ShowerApproximation.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the ShowerEventRecord class.
 */
class ShowerEventRecord: public Base {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ShowerEventRecord();

  /**
   * The destructor.
   */
  virtual ~ShowerEventRecord();
  //@}

public:

  /**
   * Return the subprocess currently showered
   */
  tSubProPtr subProcess() const { return subProcess_; }

  /**
   * Return the XComb describing the hard process.
   */
  tStdXCombPtr xcombPtr() const { return XComb_; }

  /**
   * Set the XComb describing the hard process.
   */
  void xcombPtr(tStdXCombPtr in) { XComb_ = in; }

  /**
   * Return the XComb describing the hard process.
   */
  const StandardXComb& xcomb() const { return *XComb_; }

public:

  /**
   * Return the incoming partons at the current 
   * stage of the evolution.
   */
  PPair& incoming() { return incoming_; }

  /**
   * Return the incoming partons at the current 
   * stage of the evolution.
   */
  const PPair& incoming() const { return incoming_; }

  /**
   * Return the outgoing partons at the current
   * stage of the evolution.
   */
  PList& outgoing() { return outgoing_; }

  /**
   * Return the outgoing partons at the current
   * stage of the evolution.
   */
  const PList& outgoing() const { return outgoing_; }

  /**
   * Return the intermediate particles at the current
   * stage of the evolution.
   */
  PList& intermediates() { return intermediates_; }

  /**
   * Return the intermediate particles at the current
   * stage of the evolution.
   */
  const PList& intermediates() const { return intermediates_; }

  /**
   * Return the momentum fractions.
   */
  const pair<double,double>& fractions() const { return fractions_; }

  /**
   * Return the momentum fractions.
   */
  pair<double,double>& fractions() { return fractions_; }

  /**
   * Return the PDFs
   */
  const pair<PDF,PDF>& pdfs() const { return PDFs_; }

public:

  /** @name MC@NLO diagnostics */
  //@{
  /**
   * True, if Matchbox MC@NLO S-event
   */
  bool isMCatNLOSEvent() const { return isMCatNLOSEvent_; }

  /**
   * True, if matchbox MC@NLO H-event
   */
  bool isMCatNLOHEvent() const { return isMCatNLOHEvent_; }

  /**
   * True, if Matchbox MC@NLO S-event
   */
  bool isPowhegSEvent() const { return isPowhegSEvent_; }

  /**
   * True, if matchbox MC@NLO H-event
   */
  bool isPowhegHEvent() const { return isPowhegHEvent_; }

  /**
   * True, if Matchbox MC@NLO S-event
   */
  void isMCatNLOSEvent(bool in) { isMCatNLOSEvent_ = in; }

  /**
   * True, if matchbox MC@NLO H-event
   */
  void isMCatNLOHEvent(bool in) { isMCatNLOHEvent_ = in; }

  /**
   *  Access to the shower approximation
   */
  Ptr<ShowerApproximation>::tptr showerApproximation() {
    return showerApproximation_;
  }

  /**
   *  True if truncated shower required
   */
  bool truncatedShower() const { return truncatedShower_; }

  /**
   *  True if truncated shower required
   */
  void truncatedShower(bool in) { truncatedShower_=in; }
  //@}

public:

  /**
   * Clear the event record: Give up ownership
   * on any object involved in the evolution.
   */
  virtual void clear();

protected:

  /**
   *  Set the subprocess
   */
  void subProcess(tSubProPtr in) { subProcess_ = in; }

  /**
   * Return the PDFs
   */
  pair<PDF,PDF>& pdfs() { return PDFs_; }

public:

  /**
   *  Identify the type of event
   */
  void identifyEventType();

protected:

  /**
   *  Map of colour lines used to reset colours when inserted into the event
   */
  map<ColinePtr,ColinePtr> & colourLines() {return colourLines_;}

  /**
   * Isolate the colour of the process from the rest of the event.
   * Called in the constructor
   */
  void colourIsolate(const vector<PPtr> & original, const vector<PPtr> & copy);

  /**
   * Update the colour information of a particle prior to insertion into the
   * event record.
   */
  void updateColour(PPtr particle, bool recursive);

  /**
   *  Map the colours for a given particle
   */
  void mapColour(PPtr original, PPtr copy);

private:
 
  /**
   *  Isolate a specific colour line
   */
  void isolateLine(vector<PPair>::const_iterator cit,
		   vector<PPair> & particles,
		   tcColinePtr oldline,
		   tColinePtr newline);

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerEventRecord & operator=(const ShowerEventRecord &);


private:

  /**
   * The subprocess currently showered.
   */
  SubProPtr subProcess_;

  /**
   * Pointer to the XComb which generated the hard process.
   */
  StdXCombPtr XComb_;

  /**
   * The incoming partons at the current
   * stage of the evolution.
   */
  PPair incoming_;

  /**
   * The outgoing partons at the current stage of the evolution.
   */
  PList outgoing_;

  /**
   * The intermediate particles at the current
   * stage of the evolution.
   */
  PList intermediates_;

  /**
   * The PDFs to be considered.
   */
  pair<PDF,PDF> PDFs_;

  /**
   * Momentum fractions of the incoming partons.
   */
  pair<double,double> fractions_;

private:
  
  /**
   *  Type of event
   */
  //@{
  /**
   * True, if Matchbox MC@NLO S-event
   */
  bool isMCatNLOSEvent_;

  /**
   * True, if matchbox MC@NLO H-event
   */
  bool isMCatNLOHEvent_;

  /**
   * True, if Matchbox Powheg S-event
   */
  bool isPowhegSEvent_;

  /**
   * True, if matchbox Powheg H-event
   */
  bool isPowhegHEvent_;

  /**
   * The shower approximation to provide the hard scale profile
   */
  Ptr<ShowerApproximation>::tptr showerApproximation_;

  /**
   *  Map of colour lines used to reset colours when inserted into the event
   */
  map<ColinePtr,ColinePtr> colourLines_;

  /**
   *  Whether or ont needs truncated shower
   */
  bool truncatedShower_;
  //@}

};

}

#endif /* Herwig_ShowerEventRecord_H */
