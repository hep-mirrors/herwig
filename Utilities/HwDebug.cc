//========================================================================
// file HwDebug.cc
//========================================================================

#include "HwDebug.h"

using namespace Herwig;

int HwDebug::level = 0;

SampleHistogram HwDebug::lambda1Histo(0., 1., 1./100.);
SampleHistogram HwDebug::lambda2Histo(0., 1., 1./100.);
SampleHistogram HwDebug::lambda3Histo(0., 1., 1./100.);
SampleHistogram HwDebug::CparameterHisto(0., 1., 1./100.);
SampleHistogram HwDebug::DparameterHisto(0., 1., 1./100.);
SampleHistogram HwDebug::multiplicityHisto(0., 100., 1.);
