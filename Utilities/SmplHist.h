/*
 *    --- patched version!!!, see file 'libg-diffs' ---
 */
// This may look like C code, but it is really -*- C++ -*-

#ifndef SampleHistogram_h
#ifdef __GNUG__
#pragma interface
#endif
#define SampleHistogram_h 1

#include <iostream>
#include <fstream>
#include <vector>
#include "SmplStat.h"

extern const int SampleHistogramMinimum;
extern const int SampleHistogramMaximum;
using std::vector;



/** \ingroup Utilities
 *
 *  Copyright (C) 1988 Free Software Foundation
 *  written by Dirk Grunwald (grunwald@cs.uiuc.edu)
 *
 *  This file is part of the GNU C++ Library.  This library is free
 *  software; you can redistribute it and/or modify it under the terms of
 *  the GNU Library General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your
 *  option) any later version.  This library is distributed in the hope
 *  that it will be useful, but WITHOUT ANY WARRANTY; without even the
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *  PURPOSE.  See the GNU Library General Public License for more details.
 *  You should have received a copy of the GNU Library General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */
class SampleHistogram : public SampleStatistic {

protected:

  /**
   *  Number of bins
   */
  short howManyBuckets;

  /**
   *  Number of entries in a bin
   */
  vector<int> bucketCount;
  
  /**
   *  Limits for the bins
   */
  vector<double> bucketLimit;

public:
    
  /**
   * Constructor with uniform bins
   * @param low The lower limit
   * @param h The upper limits
   * @param  bucketWidth The width of the bins
   */
  SampleHistogram(double low=0., double h=0., double bucketWidth = -1.0);

  /**
   * Constructor with non-uniform bins
   * @param low limits of the bins
   * @param size The number of bins
   */
  SampleHistogram(double loVals[], int size);

  /**
   * Destructor 
   */
  ~SampleHistogram();
  
  /**
   *  Reset the entries
   */
  void reset();

  /**
   *  Add a point
   */
  void operator+=(double);
  
  int similarSamples(double);

  /**
   *  Number of bins
   */  
  int buckets();  

  /**
   *  Lower limit of bin
   * @param i The bin
   */
  double bucketThreshold(int i);

  /**
   *  Number of entries in a given bin 
   */
  int inBucket(int i);
  void printBuckets(std::ostream&);
  void printGnuplot(char* name);
  void printMoments(char*,double,double,double,double);

};


inline int SampleHistogram:: buckets() { return(howManyBuckets); }


inline double SampleHistogram:: bucketThreshold(int i) {
    if (i < 0 || i >= howManyBuckets)
        error("invalid bucket access");
    return(bucketLimit[i]);
}


inline int SampleHistogram:: inBucket(int i) {
    if (i < 0 || i >= howManyBuckets)
        error("invalid bucket access");
    return(bucketCount[i]);
}

#endif
