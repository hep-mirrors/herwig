// -*- C++ -*-
//
// Herwig.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration, 2015 Marco A. Harrendorf
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef SRC_HERWIG_H
#define SRC_HERWIG_H

namespace Herwig {

class HerwigUI;

/**
 * A very high-level API for interacting with Herwig's different run modes.
 *
 * It's not very convenient (yet), since you'll have to provide your own 
 * HerwigUI-derived object with some fairly obscure settings.
 *
 * Much more fine-grained control is available through ThePEG::Repository.
 */
namespace API {

/**
 * Herwig read mode. Prepares a generator .run file.
 */
void read(const HerwigUI &);


/**
 * Herwig build mode
 *
 *
 */
void build(const HerwigUI &);


/**
 * Herwig integrate mode
 *
 *
 */
void integrate(const HerwigUI &);


/**
 * Herwig run mode
 *
 *
 */
void run(const HerwigUI &);


/**
 * Herwig init mode. Creates a new default repository. 
 *
 * Usually, the default repo created during the Herwig installation
 * is fine and external users will not need this mode.
 */
void init(const HerwigUI &);

}

}

#endif
