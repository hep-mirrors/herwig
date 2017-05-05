// -*- C++ -*-
//
// PhasespaceHelpers.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PhasespaceHelpers_H
#define HERWIG_PhasespaceHelpers_H

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/PDT/ParticleData.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Phasespace/MatchboxPhasespace.h"

namespace Herwig {

using namespace ThePEG;

  namespace PhasespaceHelpers {

    /**
     * \ingroup Matchbox
     * \author Simon Platzer
     * \brief General information for phase space generation
     */
    struct PhasespaceInfo {

      /**
       * The center of mass energy squared.
       */
      Energy2 sHat;

      /**
       * The center of mass energy.
       */
      Energy sqrtSHat;

      /**
       * The phase space weight.
       */
      double weight;

      /**
       * The random number generator.
       */
      StreamingRnd rnd;

      /**
       * Parameter steering from which on propagator virtualities are
       * sampled flat.
       */
      double x0;

      /**
       * Parameter steering at which virtuality singularities of
       * propagators are actually cut off.
       */
      double xc;

      /**
       * Parameter steering from which on propagator virtualities are
       * sampled flat.
       */
      Energy M0;

      /**
       * Parameter steering at which virtuality singularities of
       * propagators are actually cut off.
       */
      Energy Mc;

      /**
       * Generate a mass for the given particle type and mass range.
       */
      Energy generateMass(tcPDPtr, const pair<Energy,Energy>&);

      /**
       * Calculate a transverse momentum for the given momenta,
       * invariant pt and azimuth.
       */
      Lorentz5Momentum generateKt(const Lorentz5Momentum& p1,
				  const Lorentz5Momentum& p2,
				  Energy pt);

    };

    /**
     * \ingroup Matchbox
     * \author Simon Platzer, Ken Arnold
     * \brief A phase space tree.
     */
    struct PhasespaceTree {

      /**
       * Default constructor
       */
      PhasespaceTree()
	: massRange(ZERO,ZERO), externalId(-1), 
	  spacelike(false) {}

      /**
       * The particle running along this line.
       */
      tcPDPtr data;

      /**
       * The allowed mass range for this line.
       */
      pair<Energy,Energy> massRange;

      /**
       * The momentum running along this line.
       */
      Lorentz5Momentum momentum;

      /**
       * The external leg id of this line, if external.
       */
      int externalId;

      /**
       * The children lines; if empty this is an external line.
       */
      vector<PhasespaceTree> children;

      /**
       * External lines originating from this line.
       */
      set<int> leafs;

      /**
       * Wether or not this is a spacelike line.
       */
      bool spacelike;

      /**
       * Wether or not this is a mirrored channel.
       */
      bool doMirror;
      
      /**
       * Setup from diagram at given position.
       */
      void setup(const Tree2toNDiagram&, int pos = 0);

      /**
       * Setup mirror from diagram at given position.
       */
      void setupMirrored(const Tree2toNDiagram& diag, int pos);

      /**
       * Initialize using masses as given by mass() members of the
       * final state momenta
       */
      void init(const vector<Lorentz5Momentum>&);

      /**
       * Generate kinematics for the children
       */
      void generateKinematics(PhasespaceInfo&,
			      vector<Lorentz5Momentum>&);

      /**
       * Write phase space tree to ostream
       */
      void put(PersistentOStream&) const;

      /**
       * Read phase space tree from istream
       */
      void get(PersistentIStream&);

      /**
       * Print tree, only for debugging purposes
       */
      void print(int in = 0);

    };

  }

}

#endif // HERWIG_PhasespaceHelpers_H
