// -*- C++ -*-
#ifndef LWH_AIHistogram1D_H
#define LWH_AIHistogram1D_H

#ifndef LWH_USING_AIDA

/** @cond DONT_DOCUMENT_STRIPPED_DOWN_AIDA_INTERFACES */

#include "AIBaseHistogram.h"

namespace AIDA {

class IAnnotation;


class IHistogram : virtual public IBaseHistogram {

public:

  virtual ~IHistogram() {}

  virtual int allEntries() const = 0;
  virtual int extraEntries() const = 0;
  virtual double equivalentBinEntries() const = 0;
  virtual double sumBinHeights() const = 0;
  virtual double sumAllBinHeights() const = 0;
  virtual double sumExtraBinHeights() const = 0;
  virtual double minBinHeight() const = 0;
  virtual double maxBinHeight() const = 0;

};

class IAxis;

class IHistogram1D: virtual public IHistogram {

public:

  virtual ~IHistogram1D() {}

  virtual bool fill(double x, double weight = 1.) = 0;
  virtual double binMean(int index) const = 0;
  virtual int binEntries(int index) const = 0;
  virtual double binHeight(int index) const = 0;
  virtual double binError(int index) const = 0;
  virtual double mean() const = 0;
  virtual double rms() const = 0;
  virtual const IAxis & axis() const = 0;
  virtual int coordToIndex(double coord) const = 0;
  virtual bool add(const IHistogram1D & hist) = 0;
  virtual bool scale(double scaleFactor) = 0;
  virtual void setAttempts(long attempts) = 0;
  virtual void finishEvent() = 0;

};

}

/** @endcond */

#else
#include "AIDA/IHistogram1D.h"
#endif

#endif /* LWH_AIHistogram1D_H */
