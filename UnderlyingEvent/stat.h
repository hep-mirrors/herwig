// -*- C++ -*-
#ifndef HERWIG_Stat_H
#define HERWIG_Stat_H

namespace Herwig{
struct Stat {

  Stat() : attempted(0), accepted(0), sumw(0.0), maxXSec(CrossSection()),
       totsum(0.0) {}
  Stat(long att, long acc, double w, CrossSection x, double sumw)
    : attempted(att), accepted(acc), sumw(w), maxXSec(x),
       totsum(sumw) {}

  inline CrossSection xSec() const {
    return maxXSec*sumw/totsum;
  }

  long attempted;
  long accepted;
  double sumw;
  CrossSection maxXSec;
  double totsum;

  const Stat & operator+= (const Stat & s) {
    attempted += s.attempted;
    accepted += s.accepted;
    sumw += s.sumw;
    totsum = max(totsum, s.totsum);
    if ( totsum > 0.0 )
      maxXSec = max(maxXSec, s.maxXSec);
    else
      maxXSec += s.maxXSec;
    return *this;
  }
};
}
#endif
