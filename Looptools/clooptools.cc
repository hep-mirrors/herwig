/* -*- C++ -*-
  clooptools.cc
  the C++ file with the definitions for fortran IO redirection
  Output redirected to log file. 2007-07-18 dgrell
*/
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <cstdio>
#include <cassert>
#include <string>

#ifdef HAVE_UNISTD_H
# include "ThePEG/Repository/CurrentGenerator.h"
#endif

extern "C" {

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

  void ffini_();
  void ffexi_();
}


namespace {

#ifdef HAVE_UNISTD_H
  int start_redirection(std::string logfilename) {
    if ( ThePEG::CurrentGenerator::current().useStdOut() ) return -1;
    // redirect C stdout --- unix specific solution,
    // see C FAQ: http://c-faq.com/stdio/undofreopen.html
    int    fd;
    fflush(stdout);
    fd = dup(fileno(stdout));
    freopen(logfilename.c_str(), "a", stdout);
    return fd;
  }
  
  void stop_redirection(int fd) {
    if ( ThePEG::CurrentGenerator::current().useStdOut() ) return;
    fflush(stdout);
    close(fileno(stdout));
    dup2(fd, fileno(stdout));
    close(fd);
    clearerr(stdout);
  }
#else
  int start_redirection(std::string) {
    return -1;
  }
  
  void stop_redirection(int) {}
#endif

} // namespace

namespace Herwig {
  namespace Looptools {

    static int initcount = 0;

    void ffini(std::string logfilename) {
      assert( initcount >= 0 );
      if ( initcount == 0 ) {
	int rd = start_redirection(logfilename);
	ffini_();
	stop_redirection(rd);
      }
      ++initcount;
    }

    void ffexi(std::string logfilename) {
      assert( initcount > 0 );
      --initcount;
      if ( initcount == 0 ) {
	int rd = start_redirection(logfilename);
	ffexi_();
	stop_redirection(rd);
      }
    }

  } // namespace Looptools
} // namespace Herwig
