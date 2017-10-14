// -*- C++ -*-
//
// Filesystem.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2017-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "Filesystem.h"
#include <cstdio>

#include <config.h>
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#else
#error "Need sys/stat.h"
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#else
#error "Need unistd.h"
#endif

namespace Herwig {

namespace filesystem {

	/// Check if file/dir location exists
	bool exists(const std::string & path) {
		struct stat sb;
		return stat(path.c_str(), &sb) == 0;
	}

	/// Check if location is a directory
	bool is_directory(const std::string & path) {
		struct stat sb;
		return stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode);
	}

	/// Make a directory
	bool create_directory(const std::string & path) {
		const mode_t mode = 0755;
		if (mkdir(path.c_str(), mode) == 0) {
			return true;
		}
		else {
			perror("Herwig::filesystem::mkdir()");
			return false;
		}

	}

}

}
