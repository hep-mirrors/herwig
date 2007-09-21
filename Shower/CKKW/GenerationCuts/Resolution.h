// -*- C++ -*-
#ifndef HERWIG_Resolution_H
#define HERWIG_Resolution_H

#include "CutsInterface.h"

/**
 * Set the resolution parameter (in GeV^2)
 */
const double resolution = sqr(10.);

/**
 * Define the resolution definition to be
 * used.
 */
#define useResolution_KtTilde

/**
 * Return true, if the event passes the cut.
 * passFlag > 0 indicates, that the event did pass
 * the cut. passFlag < 0 indicates, that something
 * has gone wrong. passFlag = 0 indicates, that event
 * did not pass the cut.
 */
extern "C" void checkCut_ (int * passFlag);

#endif
