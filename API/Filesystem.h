// -*- C++ -*-
//
// Filesystem.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2017-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include <string>

namespace Herwig {

namespace filesystem {

	/// Check if file/dir location exists
	bool exists(const std::string & path);

	/// Check if location is a directory
	bool is_directory(const std::string & path);

	/// Make a directory
	bool create_directory(std::string path);

}

}
