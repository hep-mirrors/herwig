// -*- C++ -*-
//
// config.h is part of ExSample -- A Library for Sampling Sudakov-Type Distributions
//
// Copyright (C) 2008-2011 Simon Platzer -- simon.plaetzer@desy.de
//
// ExSample is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef EXSAMPLE_config_h_included
#define EXSAMPLE_config_h_included

#include <vector>
#include <map>
#include <set>
#include <string>

#include <cmath>
#include <cassert>
#include <climits>

#include <algorithm>
#include <numeric>
#include <limits>

#include <boost/utility.hpp>
#include <boost/scoped_array.hpp>
#include <boost/scoped_ptr.hpp>

#define EXSAMPLE_has_ThePEG

#ifdef EXSAMPLE_has_ThePEG

#include "ThePEG/Persistency/PersistentOStream.h"

#endif // EXSAMPLE_has_ThePEG

namespace exsample {

  static const unsigned long parameter_hash_bits = 512;

}

#endif // EXSAMPLE_config_h_included
