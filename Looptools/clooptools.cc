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

extern "C" {

#ifdef HAVE_UNISTD_H
# include <unistd.h>
#endif

  void ffini_();
  void ffexi_();
}


namespace {
  struct RedirectionInfo {
    RedirectionInfo(int fdin) 
      : fd(fdin) {}
    int fd;
  };

#ifdef HAVE_UNISTD_H
  RedirectionInfo start_redirection(std::string logfilename) {
    // redirect C stdout --- unix specific solution,
    // see C FAQ: http://c-faq.com/stdio/undofreopen.html
    int    fd;
    fflush(stdout);
    fd = dup(fileno(stdout));
    freopen(logfilename.c_str(), "a", stdout);
    return RedirectionInfo(fd);
  }
  
  void stop_redirection(RedirectionInfo rdinfo) {
    fflush(stdout);
    close(fileno(stdout));
    dup2(rdinfo.fd, fileno(stdout));
    close(rdinfo.fd);
    clearerr(stdout);
  }
#else
  RedirectionInfo start_redirection(std::string) {
    return RedirectionInfo();
  }
  
  void stop_redirection(RedirectionInfo) {}
#endif

} // namespace

namespace Herwig {
  namespace Looptools {

    static int initcount = 0;

    void ffini(std::string logfilename) {
      assert( initcount >= 0 );
      if ( initcount == 0 ) {
	RedirectionInfo rd = start_redirection(logfilename);
	ffini_();
	stop_redirection(rd);
      }
      ++initcount;
    }

    void ffexi(std::string logfilename) {
      assert( initcount > 0 );
      --initcount;
      if ( initcount == 0 ) {
	RedirectionInfo rd = start_redirection(logfilename);
	ffexi_();
	stop_redirection(rd);
      }
    }

  } // namespace Looptools
} // namespace Herwig
