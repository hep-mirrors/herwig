// -*- C++ -*-
//
// Filesystem.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2017-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "Filesystem.h"
#include <cstdio>
#include <stack>

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

#include <iostream>

namespace {
	std::string & strip_last_dir(std::string & s) {
		size_t last = s.rfind('/');
		if ( last != std::string::npos )
			s.assign(s,0,last);
		else
			s.assign("");
		return s;
	}
}

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
	bool create_directory(std::string path) {
		static const mode_t mode = 0755;
		std::stack<std::string> dirs;
		dirs.push(path);
		while ( strip_last_dir(path) != "" ) {
			dirs.push(path);
		}

		while ( !dirs.empty() ) {
			const auto & top = dirs.top();

			if ( is_directory(top) ) {
				dirs.pop();
				continue;
			}
			
			if ( exists(top) ) {
				std::cerr << "Path " << top << " exists, but isn't a directory" << '\n';
				return false;
			}

			if ( mkdir(top.c_str(), mode) ) {
				static const std::string msg = "Herwig::filesystem::mkdir("+path+")";
				perror(msg.c_str());
				return false;
			}

			dirs.pop();
		}

		return true;
	}

}

}

