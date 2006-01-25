// shamelessly copied from libg++ and modified for our own needs
// SG 08/2002

/*
 *    --- patched version!!!, see file 'libg-diffs' ---
 */
// This may look like C code, but it is really -*- C++ -*-
/* 
Copyright (C) 1988 Free Software Foundation
    written by Dirk Grunwald (grunwald@cs.uiuc.edu)

This file is part of the GNU C++ Library.  This library is free
software; you can redistribute it and/or modify it under the terms of
the GNU Library General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.  This library is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU Library General Public License for more details.
You should have received a copy of the GNU Library General Public
License along with this library; if not, write to the Free Software
Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include <iostream>
#include <limits>
#include <cmath>
#include "SmplHist.h"

namespace {
  const double myHUGE_VAL = std::numeric_limits<double>::max();
}

using std::ostream;
using std::cout;   
using std::cerr;   
using std::endl;   
using std::ofstream;

const int SampleHistogramMinimum = -2;
const int SampleHistogramMaximum = -1;

SampleHistogram::SampleHistogram(double low, double high, double width)
{
    if (high < low) {
	double t = high;
	high = low;
	low = t;
    }

    if (width < 0) {
      if (width > -1) {
	cerr << "Herwig++/Utilities/SampleHistogram::SampleHistogram(): zero bins!" << endl;
      }
      int bins = -int(width);
      width = (high - low)/(bins);
      //      howManyBuckets = int((high - low) / width) + 2;            
      howManyBuckets = bins+2;            
    } else {     
      howManyBuckets = int((high - low) / width) + 2;
    }
    
    bucketCount.resize(howManyBuckets);
    bucketLimit.resize(howManyBuckets);
    double lim = low;
    for (int i = 0; i < howManyBuckets-1; i++) {
	bucketCount[i] = 0;
	bucketLimit[i] = lim;
	lim += width;
    }
    bucketLimit[howManyBuckets-1] = myHUGE_VAL;	/* from math.h */
}

SampleHistogram::SampleHistogram(double loVals[], int size)
{
  howManyBuckets = size+1;
  bucketCount.resize(howManyBuckets);
  bucketLimit.resize(howManyBuckets);
  for (int i = 0; i < howManyBuckets-1; i++) {
    bucketCount[i] = 0;
    bucketLimit[i] = loVals[i];
  }
  bucketLimit[howManyBuckets-1] = myHUGE_VAL;	/* from math.h */
}

SampleHistogram::~SampleHistogram() {}

void
SampleHistogram::operator+=(double value)
{
  if(isnan(value)) return;
    int i;
    for (i = 0; i < howManyBuckets; i++) {
	if (value < bucketLimit[i]) break;
    }
    bucketCount[i]++;
    this->SampleStatistic::operator+=(value);
}

int
SampleHistogram::similarSamples(double d)
{
    int i;
    for (i = 0; i < howManyBuckets; i++) {
	if (d < bucketLimit[i]) return(bucketCount[i]);
    }
    return(0);
}

void
SampleHistogram::printBuckets(ostream& s)
{
    for(int i = 0; i < howManyBuckets; i++) {
	if (bucketLimit[i] >= myHUGE_VAL) {
	    s << "< max : " << bucketCount[i] << "\n";
	} else {
	    s << "< " << bucketLimit[i] << " : " << bucketCount[i] << "\n";
	}
    }
}


// void SampleHistogram::printGnuplot(char* name)
// {
//   ofstream out(name);
//   if (!out) {
//     cerr << "SampleHistoGram::printGnuplot: ERROR! Can't open file" << endl;
//   }

//   time_t now_t;
//   now_t = time(0);
//   out << "# created " << ctime(&now_t)
//       << "# by SampleHistogram::printGnuplot" 
//       << " (simply GNUPLOT plot with histeps)" << endl
//       << "# " << this->samples() << " entries, mean +/- sigma = " 
//       << this->mean() << " +/- " << this->stdDev() << endl
//       << "# 1:xmid 2:entr 3:entr n1 4:estd err " 
//       << "5:err n1 6:xlow 7:xhi 8:entr/tot"
//       << endl;
  
//   double delta;
//   double norm = 0.0; 
//   for(int i = 0; i < howManyBuckets-2; i++) 
//     norm += double(bucketCount[i+1])*(bucketLimit[i+1] - bucketLimit[i]);

//   for(int i = 0; i < howManyBuckets-2; i++) {
//     delta = (bucketLimit[i+1] - bucketLimit[i])/2.;
//     out << bucketLimit[i] + delta << "\t" 
// 	<< bucketCount[i+1] << "\t"
// 	<< bucketCount[i+1]/norm << "\t"
// 	<< (bucketCount[i+1] == 0 ? 0.0 : sqrt(double(bucketCount[i+1])))
// 	<< "\t"
// 	<< (bucketCount[i+1] == 0 ? 0.0 : sqrt(double(bucketCount[i+1]))/norm) 
// 	<< "\t"
// 	<< bucketLimit[i] << "\t" 
// 	<< bucketLimit[i] + 2.*delta << "\t"
// 	<< (double) bucketCount[i+1]/(this->samples()) << "\t"
// 	<< "\n";
//   }
//   out.close();
//}

void SampleHistogram::printGnuplot(char* name)
{
  ofstream out(name);
  if (!out) {
    cerr << "SampleHistoGram::printGnuplot: ERROR! Can't open file" << endl;
  }

  time_t now_t;
  now_t = time(0);
  out << "# created " << ctime(&now_t)
      << "# by SampleHistogram::printGnuplot (simply GNUPLOT plot with histeps)" << endl
      << "# " << this->samples() << " entries, mean +/- sigma = " 
      << this->mean() << " +/- " << this->stdDev() << endl
      << "# 1:xmid 2:entr 3:entr n1 4:estd err 5:err n1 6:xlow 7:xhi 8:entr/tot"
      << endl;
  double delta;
  for(int i = 0; i < howManyBuckets-1; i++) {
    delta = (bucketLimit[i+1] - bucketLimit[i])/2.;
    
    out << bucketLimit[i] + delta << "\t" 
	<< bucketCount[i+1] << "\t"
	<< double(bucketCount[i+1]/(2.*delta*this->samples())) << "\t"
	<< (bucketCount[i+1] == 0 ? 0.0 : 
	    sqrt(double(bucketCount[i+1]))) 
	<< "\t"
	<< (bucketCount[i+1] == 0 ? 0.0 : 
	    sqrt(double(bucketCount[i+1]))/(2.*delta*(this->samples()))) 
	<< "\t"
	<< bucketLimit[i] << "\t" 
	<< bucketLimit[i] + 2.*delta << "\t"
	<< double(bucketCount[i+1]/(this->samples())) << "\t"
	<< "\n";
  }
  out.close();
}

void SampleHistogram::printMoments(char* name, double Nmax, double dN, 
				   double x0, double x1) {
  ofstream out(name);
  if (!out) {
    cerr << "SampleHistoGram::printMoments: ERROR! Can't open file" << endl;
  }

  time_t now_t;
  now_t = time(0);
  out << "# created " << ctime(&now_t)
      << "# by SampleHistogram::printMoments(..., "
      << Nmax << ", " << dN << ")" << endl
      << "# " << this->samples() << " entries, mean +/- sigma = " 
      << this->mean() << " +/- " << this->stdDev() << endl;

  double x0N, x1N, delta, hi;
  for (double N=dN; N < Nmax; N += dN) {
    double fN = 0.0;
    for(int i = 0; i < howManyBuckets-1; i++) {
      x0N = pow(bucketLimit[i], N);
      x1N = pow(bucketLimit[i+1], N);
      delta = (bucketLimit[i+1] - bucketLimit[i]);
      if (delta > 0 && this->samples() > 0 
	  && bucketLimit[i] >= x0 && bucketLimit[i] <= x1
	  && bucketLimit[i+1] >= x0 && bucketLimit[i+1] <= x1) {
	hi = double(bucketCount[i+1]/(delta*(this->samples())));
	fN += hi*(x1N-x0N)/N;
      }
    }
    out << N 
	<< "  " 
	<< fN << endl;
  }
  out.close();
}


void
SampleHistogram::reset()
{
    this->SampleStatistic::reset();
    if (howManyBuckets > 0) {
	for (register int i = 0; i < howManyBuckets; i++) {
	    bucketCount[i] = 0;
	}
    }
}

