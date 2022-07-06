// -*- C++ -*-
//
// LesHouches.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_LesHouches_H
#define THEPEG_LesHouches_H
//
// This is the declaration of the LesHouches class.
//

#include "ThePEG/Config/ThePEG.h"

namespace ThePEG {

/**
 * The HEPRUP class is a simple container corresponding to the Les
 * Houches accord (hep-ph/0109068) common block with the same
 * name. The members are named in the same way as in the common
 * block. However, fortran arrays are represented by vectors, except
 * for the arrays of length two which are represented by pair objects.
 */
class HEPRUP {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  HEPRUP() : IDWTUP(0), NPRUP(0) {}
  //@}

public:

  /**
   * Set the NPRUP variable, corresponding to the number of
   * sub-processes, to \a nrup, and resize all relevant vectors
   * accordingly.
   */
  void resize(int nrup) {
    NPRUP = nrup;
    resize();
  }

  /**
   * Assuming the NPRUP variable, corresponding to the number of
   * sub-processes, is correctly set, resize the relevant vectors
   * accordingly.
   */
  void resize() {
    XSECUP.resize(NPRUP);
    XERRUP.resize(NPRUP);
    XMAXUP.resize(NPRUP);
    LPRUP.resize(NPRUP);
  }

  /**
   * PDG id's of beam particles. (first/second is in +/-z direction).
   */
  pair<long,long> IDBMUP;

  /**
   * Energy of beam particles given in GeV.
   */
  pair<double,double> EBMUP;

  /**
   * The author group for the PDF used for the beams according to the
   * PDFLib specification.
   */
  pair<int,int> PDFGUP;

  /**
   * The id number the PDF used for the beams according to the
   * PDFLib specification.
   */
  pair<int,int> PDFSUP;

  /**
   * Master switch indicating how the ME generator envisages the
   * events weights should be interpreted according to the Les Houches
   * accord.
   */
  int IDWTUP;

  /**
   * The number of different subprocesses in this file (should
   * typically be just one)
   */
  int NPRUP;

  /**
   * The cross sections for the different subprocesses in pb.
   */
  vector<double> XSECUP;

  /**
   * The statistical error in the cross sections for the different
   * subprocesses in pb.
   */
  vector<double> XERRUP;

  /**
   * The maximum event weights (in XWGTUP) for different subprocesses.
   */
  vector<double> XMAXUP;

  /**
   * The subprocess code for the different subprocesses.
   */
  vector<int> LPRUP;

};


/**
 * The HEPEUP class is a simple container corresponding to the Les
 * Houches accord (hep-ph/0109068) common block with the same
 * name. The members are named in the same way as in the common
 * block. However, fortran arrays are represented by vectors, except
 * for the arrays of length two which are represented by pair objects.
 */
class HEPEUP {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  HEPEUP()
    : NUP(0), IDPRUP(0), XWGTUP(0.0), XPDWUP(0.0, 0.0),
      SCALUP(0.0), AQEDUP(0.0), AQCDUP(0.0) {}
  //@}

public:

  /**
   * Set the NUP variable, corresponding to the number of particles in
   * the current event, to \a nup, and resize all relevant vectors
   * accordingly.
   */
  void resize(int nup) {
    NUP = nup;
    resize();
  }

  /**
   * Assuming the NUP variable, corresponding to the number of
   * particles in the current event, is correctly set, resize the
   * relevant vectors accordingly.
   */
  void resize() {
    IDUP.resize(NUP);
    ISTUP.resize(NUP);
    MOTHUP.resize(NUP);
    ICOLUP.resize(NUP);
    PUP.resize(NUP);
    VTIMUP.resize(NUP);
    SPINUP.resize(NUP);
  }

  /**
   * The number of particle entries in the current event.
   */
  int NUP;

  /**
   * The subprocess code for this event (as given in LPRUP).
   */
  int IDPRUP;

  /**
   * The weight for this event.
   */
  double XWGTUP;

  /**
   * The PDF weights for the two incoming partons. Note that this
   * variable is not present in the current LesHouches accord
   * (hep-ph/0109068), hopefully it will be present in a future
   * accord.
   */
  pair<double,double> XPDWUP;

  /**
   * The scale in GeV used in the calculation of the PDF's in this
   * event.
   */
  double SCALUP;

  /**
   * The value of the QED coupling used in this event.
   */
  double AQEDUP;

  /**
   * The value of the QCD coupling used in this event.
   */
  double AQCDUP;

  /**
   * The PDG id's for the particle entries in this event.
   */
  vector<long> IDUP;

  /**
   * The status codes for the particle entries in this event.
   */
  vector<int> ISTUP;

  /**
   * Indices for the first and last mother for the particle entries in
   * this event.
   */
  vector< pair<int,int> > MOTHUP;

  /**
   * The colour-line indices (first(second) is (anti)colour) for the
   * particle entries in this event.
   */
  vector< pair<int,int> > ICOLUP;

  /**
   * Lab frame momentum (Px, Py, Pz, E and M in GeV) for the particle
   * entries in this event.
   */
  vector< array<double,5> > PUP;

  /**
   * Invariant lifetime (c*tau, distance from production to decay im
   * mm) for the particle entries in this event.
   */
  vector<double> VTIMUP;

  /**
   * Spin info for the particle entries in this event given as the
   * cosine of the angle between the spin vector of a particle and the
   * 3-momentum of the decaying particle, specified in the lab frame.
   */
  vector<double> SPINUP;

};

}

#endif /* THEPEG_LesHouches_H */
