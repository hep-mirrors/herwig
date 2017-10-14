// -*- C++ -*-
//
// Progress.h is a part of Herwig - A multi-purpose Monte Carlo event generator
//

//  boost progress.hpp header file  ------------------------------------------//

//  Copyright Beman Dawes 1994-99.  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/timer for documentation.

//  Revision History

//  2017-10-14 Modified for Herwig (D. Grellscheid)

//   1 Dec 01  Add leading progress display strings (suggested by Toon Knapen)
//  20 May 01  Introduce several static_casts<> to eliminate warning messages
//             (Fixed by Beman, reported by Herve Bronnimann)
//  12 Jan 01  Change to inline implementation to allow use without library
//             builds. See docs for more rationale. (Beman Dawes) 
//  22 Jul 99  Name changed to .hpp
//  16 Jul 99  Second beta
//   6 Jul 99  Initial boost version

#ifndef HERWIG_PROGRESS_H
#define HERWIG_PROGRESS_H

#include <iosfwd>           // for ostream, cout, etc
#include <string>             // for string

namespace Herwig {

//  progress_display  --------------------------------------------------------//

//  progress_display displays an appropriate indication of 
//  progress at an appropriate place in an appropriate form.

// NOTE: (Jan 12, 2001) Tried to change unsigned long to boost::uintmax_t, but
// found some compilers couldn't handle the required conversion to double.
// Reverted to unsigned long until the compilers catch up. 

class progress_display
{
 public:
  explicit progress_display( unsigned long expected_count,
                             std::ostream & os,
                             const std::string & s1 = "\n", //leading strings
                             const std::string & s2 = "",
                             const std::string & s3 = "" );
  progress_display( const progress_display& ) = delete;
  progress_display& operator=( const progress_display& ) = delete;
  void           restart( unsigned long expected_count );

  unsigned long  operator+=( unsigned long increment );

  unsigned long  operator++()           { return operator+=( 1 ); }
  unsigned long  count() const          { return _count; }
  unsigned long  expected_count() const { return _expected_count; }

  private:
  std::ostream &     m_os;  // may not be present in all imps
  const std::string  m_s1;  // string is more general, safer than 
  const std::string  m_s2;  //  const char *, and efficiency or size are
  const std::string  m_s3;  //  not issues

  unsigned long _count, _expected_count, _next_tic_count;
  unsigned int  _tic;
  void display_tic();
};

} // namespace Herwig

#endif  // HERWIG_PROGRESS_H
