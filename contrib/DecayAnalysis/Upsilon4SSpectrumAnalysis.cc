// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Upsilon4SSpectrumAnalysis class.
//

#include "Upsilon4SSpectrumAnalysis.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Herwig;

IBPtr Upsilon4SSpectrumAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr Upsilon4SSpectrumAnalysis::fullclone() const {
  return new_ptr(*this);
}

void Upsilon4SSpectrumAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  set<tPPtr> particles;
  event->select(inserter(particles),ThePEG::AllSelector());
  set<tPPtr> finalState;
  event->selectFinalState(inserter(finalState));
  for(set<tPPtr>::const_iterator it=finalState.begin();it!=finalState.end();++it) {
    double p = (**it).momentum().vect().mag()/GeV;
    switch (abs((**it).id())) {
    case ParticleID::piplus:
      if(abs((**it).parents()[0]->id())!=ParticleID::Lambda0&&
	 abs((**it).parents()[0]->id())!=ParticleID::K_S0)
	if(p>=0.075&&p<=2.8) *_spectrumpipNo += p;
      if(p>=0.075&&p<=2.8) *_spectrumpipAll += p;
      break;
    case ParticleID::Kplus:
      if(p>=0.175&&p<=2.8) *_spectrumKpA +=p;
      if(p>=0.2  &&p<=2.8) *_spectrumKpB +=p;
      break;
    case ParticleID::pplus:
      if(abs((**it).parents()[0]->id())!=ParticleID::Lambda0)
	if(p>=0.3&&p<=2.2)  *_spectrumprotonNo += p;
      if(p>=0.3&&p<=2.2)  *_spectrumprotonAll += p;
    case ParticleID::K_L0:
      if(p>=0.2&&p<=2.7) *_spectrumK0 += p;
    }
  }
  for(set<tPPtr>::const_iterator it=particles.begin();it!=particles.end();++it) {
    if(!(**it).children().empty()&&(**it).children()[0]->id()==(**it).id())
      continue;
    double p = (**it).momentum().vect().mag()/GeV;
    switch (abs((**it).id())) {
    case ParticleID::pi0:
      if(p<=3.) *_spectrumpi0 += p;
//     case ParticleID::K_S0:
//       if(p>=0.2&&p<=2.7) *_spectrumK0 += p;
    }
  }
}

NoPIOClassDescription<Upsilon4SSpectrumAnalysis> 
Upsilon4SSpectrumAnalysis::initUpsilon4SSpectrumAnalysis;
// Definition of the static class description member.

void Upsilon4SSpectrumAnalysis::Init() {

  static ClassDocumentation<Upsilon4SSpectrumAnalysis> documentation
    ("There is no documentation for the Upsilon4SSpectrumAnalysis class");

}

void Upsilon4SSpectrumAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + 
    string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  // output chi sq to log file for analysis
  double chisq=0.,minfrac=0.05;
  unsigned int ndegrees;
  _spectrumpipNo->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square  = " << chisq << " for " << ndegrees 
		     << " degrees of freedom for ARGUS pi+ distribution "
		     << "excluding K_S0 and Lambda decay products\n";
  _spectrumpipAll->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square  = " << chisq << " for " << ndegrees 
		     << " degrees of freedom for ARGUS pi+ distribution "
		     << "including K_S0 and Lambda decay products\n";
  _spectrumpi0->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square  = " << chisq << " for " << ndegrees 
		     << " degrees of freedom for Belle pi0 distribution\n";
  _spectrumKpA->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square  = " << chisq << " for " << ndegrees 
		     << " degrees of freedom for ARGUS K+ distribution\n";
  _spectrumKpB->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square  = " << chisq << " for " << ndegrees 
		     << " degrees of freedom for ARGUS K+ decay in flight"
		     << "distribution\n";
  _spectrumprotonNo->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square  = " << chisq << " for " << ndegrees 
		     << " degrees of freedom for ARGUS proton distribution "
		     << "excluding Lambda decay products\n";
  _spectrumprotonAll->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square  = " << chisq << " for " << ndegrees 
		     << " degrees of freedom for ARGUS proton distribution "
		     << "including Lambda decay products\n";
  _spectrumK0->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square  = " << chisq << " for " << ndegrees 
		     << " degrees of freedom for ARGUS K0 distribution\n";
  // output the histograms
  using namespace HistogramOptions;
  _spectrumpipNo->topdrawOutput(output,Frame|Errorbars,"RED",
				"P2+3 excluding K2030S1 and L decay products (ARGUS)",
				"GXMX            X XX X     F               ",
				"1/SdS/dp/GeV2-13","  G G       X  X",
				"p/GeV","     ");
  _spectrumpipAll->topdrawOutput(output,Frame|Errorbars,"RED",
				 "P2+3 including K2030S1 and L decay products (ARGUS)",
				 "GX X            X XX X     F               ",
				 "1/SdS/dp/GeV2-13","  G G       X  X",
				 "p/GeV","     ");
  _spectrumpi0->topdrawOutput(output,Frame|Errorbars,"RED",
			      "P203 (Belle)",
			      "GX X        ",
			      "1/SdS/dp/GeV2-13","  G G       X  X",
			      "p/GeV","     ");
  _spectrumKpA->topdrawOutput(output,Frame|Errorbars,"RED",
			      "K2+3 (ARGUS)",
			      " XMX",
			      "1/SdS/dp/GeV2-13","  G G       X  X",
			      "p/GeV","     ");
  _spectrumKpB->topdrawOutput(output,Frame|Errorbars,"RED",
			      "K2+3 decay in flight (ARGUS)",
			      " XMX",
			      "1/SdS/dp/GeV2-13","  G G       X  X",
			      "p/GeV","     ");
  _spectrumprotonNo->topdrawOutput(output,Frame|Errorbars,"RED",
				   "proton excluding L decay products (ARGUS)",
				   "                 F                       ",
				   "1/SdS/dp/GeV2-13","  G G       X  X",
				   "p/GeV","     ");
  _spectrumprotonAll->topdrawOutput(output,Frame|Errorbars,"RED",
				    "proton including L decay products (ARGUS)",
				    "                 F                       ",
				    "1/SdS/dp/GeV2-13","  G G       X  X",
				    "p/GeV","     ");
  _spectrumK0->topdrawOutput(output,Frame|Errorbars,"RED",
			     "K203 (ARGUS)",
			     " X X        ",
			     "1/SdS/dp/GeV2-13","  G G       X  X",
			     "p/GeV","     ");
}

void Upsilon4SSpectrumAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  // ARGUS pi+ spectrum excluding K^0_S and Lambda decays
  double x1[]={0.075,0.100,0.125,0.150,0.175,0.200,0.225,0.250,0.275,0.300,
	       0.325,0.350,0.375,0.400,0.425,0.450,0.475,0.500,0.525,0.550,
	       0.575,0.600,0.625,0.650,0.675,0.700,0.725,0.750,0.775,0.800,
	       0.900,1.000,1.100,1.200,1.300,1.400,1.500,1.600,1.700,1.800,
	       1.900,2.000,2.100,2.200,2.300,2.400,2.500,2.600,2.700,2.800};
  double y1[]={ 8.104, 9.439, 9.967,11.238,11.217,12.208,12.153,12.148,12.194,11.811,
	       11.650,11.010,10.085,10.160, 9.378, 8.957, 8.198, 7.801, 6.806, 6.175,
		6.222, 5.380, 5.309, 5.096, 4.193, 3.882, 3.939, 3.785, 3.370, 2.732,
		1.861, 1.160, 0.811, 0.622, 0.433, 0.454, 0.208, 0.141, 0.175, 0.102,
		0.090, 0.123, 0.066, 0.140, 0.010, 0.021, 0.026,-0.015,-0.011};
  double st1[]={0.255,0.273,0.283,0.310,0.335,0.351,0.365,0.369,0.367,0.363,
		0.358,0.351,0.344,0.335,0.325,0.317,0.307,0.301,0.293,0.284,
		0.273,0.267,0.258,0.251,0.247,0.238,0.232,0.224,0.222,0.112,
		0.105,0.105,0.110,0.115,0.110,0.096,0.080,0.071,0.064,0.059,
		0.053,0.049,0.044,0.038,0.035,0.031,0.027,0.025,0.024};
  double sy1[]={0.332,0.239,0.385,0.376,0.307,0.337,0.292,0.293,0.366,0.248,
		0.325,0.286,0.232,0.188,0.182,0.231,0.145,0.262,0.103,0.056,
		0.151,0.073,0.091,0.105,0.129,0.215,0.134,0.148,0.116,0.098,
		0.072,0.089,0.089,0.136,0.064,0.075,0.072,0.058,0.043,0.028,
		0.029,0.029,0.033,0.017,0.027,0.034,0.014,0.025,0.016};
  double total=0.0;
  for(unsigned int ix=0;ix<49;++ix) {
    total+=y1[ix]*(x1[ix+1]-x1[ix]);
    st1[ix] = sqrt(sqr(st1[ix])+sqr(sy1[ix]));
  }
  vector<double> bins,data,error;
  bins  = vector<double>( x1, x1+50);
  data  = vector<double>( y1, y1+49);
  error = vector<double>(st1,st1+49);
  for(unsigned int ix=0;ix<data.size();++ix) {
    data[ix]/=total;
    error[ix]/=total;
  }
  _spectrumpipNo = new_ptr(Histogram(bins,data,error));
  // ARGUS pi+ spectrum including K^0_s and Lambda decays
  double y2[]={ 9.356,11.121,12.036,13.575,13.772,14.786,14.514,14.268,14.179,13.658,
		13.367,12.649,11.562,11.514,10.650,10.112, 9.271, 8.762, 7.645, 6.957,
		6.924, 6.021, 5.877, 5.636, 4.671, 4.315, 4.341, 4.157, 3.689, 2.983,
		2.042, 1.280, 0.891, 0.676, 0.473, 0.477, 0.223, 0.150, 0.181, 0.105,
		0.091, 0.125, 0.067, 0.140, 0.010, 0.021, 0.026,-0.015,-0.011};
  double st2[]={0.255,0.273,0.283,0.310,0.335,0.351,0.365,0.369,0.367,0.363,
		0.358,0.351,0.344,0.335,0.325,0.317,0.307,0.301,0.293,0.284,
		0.273,0.267,0.258,0.251,0.247,0.238,0.232,0.224,0.222,0.112,
		0.105,0.105,0.110,0.115,0.110,0.096,0.080,0.071,0.064,0.059,
		0.053,0.049,0.044,0.038,0.035,0.031,0.027,0.025,0.024};
  double sy2[]={0.383,0.281,0.465,0.455,0.377,0.408,0.348,0.344,0.426,0.286,
		0.373,0.328,0.266,0.213,0.207,0.261,0.164,0.294,0.116,0.063,
		0.168,0.082,0.100,0.116,0.144,0.239,0.147,0.162,0.127,0.107,
		0.079,0.099,0.098,0.148,0.070,0.079,0.077,0.061,0.045,0.029,
		0.030,0.030,0.033,0.017,0.027,0.034,0.014,0.025,0.016};
  total=0.0;
  for(unsigned int ix=0;ix<49;++ix) {
    total+=y2[ix]*(x1[ix+1]-x1[ix]);
    st2[ix] = sqrt(sqr(st2[ix])+sqr(sy2[ix]));
  }
  bins  = vector<double>( x1, x1+50);
  data  = vector<double>( y2, y2+49);
  error = vector<double>(st2,st2+49);
  for(unsigned int ix=0;ix<data.size();++ix) {
    data[ix]/=total;
    error[ix]/=total;
  }
  _spectrumpipAll = new_ptr(Histogram(bins,data,error));
  // ARGUS K+- spectrum
  double x3[]={0.175,0.200,0.225,0.250,0.275,0.300,0.325,0.350,0.375,0.400,
	       0.425,0.450,0.475,0.500,0.525,0.550,0.575,0.600,0.625,0.650,
	       0.675,0.700,0.725,0.750,0.775,0.800,0.900,1.000,1.100,1.200,
	       1.300,1.400,1.500,1.600,1.700,1.800,1.900,2.000,2.100,2.200,
	       2.300,2.400,2.500,2.600,2.700,2.800};
  double y3[]={ 0.678, 1.032, 1.240, 1.386, 1.423, 1.501, 1.611, 1.893, 1.834, 2.083,
		2.009, 2.000, 2.130, 1.952, 1.887, 1.822, 1.809, 1.904, 1.732, 1.571,
		1.521, 1.405, 1.392, 1.379, 1.200, 0.917, 0.957, 0.659, 0.403, 0.572,
		0.535, 0.208, 0.312, 0.160, 0.124, 0.013, 0.065,-0.062,-0.026, 0.024,
	       -0.021, 0.014,-0.032, 0.040, 0.046};
  double st3[]={0.081,0.086,0.087,0.087,0.093,0.092,0.094,0.098,0.100,0.099,
		0.103,0.103,0.104,0.106,0.107,0.108,0.108,0.109,0.111,0.113,
		0.112,0.113,0.114,0.112,0.118,0.060,0.066,0.079,0.094,0.103,
		0.101,0.084,0.067,0.056,0.049,0.044,0.039,0.039,0.037,0.035,
		0.033,0.036,0.033,0.034,0.037};
  double sy3[]={0.097,0.095,0.108,0.071,0.126,0.071,0.096,0.084,0.057,0.075,
		0.068,0.056,0.084,0.094,0.086,0.086,0.072,0.071,0.086,0.063,
		0.060,0.133,0.091,0.156,0.115,0.051,0.054,0.092,0.063,0.131,
		0.029,0.074,0.042,0.019,0.024,0.030,0.011,0.021,0.010,0.027,
		0.023,0.023,0.023,0.011,0.022};
  total=0.0;
  for(unsigned int ix=0;ix<45;++ix) {
    total+=y3[ix]*(x3[ix+1]-x3[ix]);
    st3[ix] = sqrt(sqr(st3[ix])+sqr(sy3[ix]));
  }
  bins  = vector<double>( x3, x3+46);
  data  = vector<double>( y3, y3+45);
  error = vector<double>(st3,st3+45);
  for(unsigned int ix=0;ix<data.size();++ix) {
    data[ix]/=total;
    error[ix]/=total;
  }
  _spectrumKpA = new_ptr(Histogram(bins,data,error));
  // kaon decay in flight
  double x4[]={0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900,1.000,1.100,
	       1.200,1.300,1.400,1.500,1.600,1.700,1.800,1.900,2.000,2.100,
	       2.200,2.300,2.400,2.500,2.600,2.700,2.800};
  double y4[]={ 1.21, 1.71, 1.75, 2.01, 1.45, 1.20, 0.98, 0.97, 0.64, 0.13,
		0.56, 0.23, 0.09, 0.15, 0.09,-0.11, 0.01,-0.07,-0.03, 0.19,
		0.27, 0.15,-0.11,-0.06, 0.04, 0.02};
   double st4[]={0.09,0.12,0.15,0.18,0.18,0.18,0.18,0.20,0.17,0.17,
		 0.16,0.17,0.13,0.15,0.13,0.17,0.15,0.09,0.09,0.19,
		 0.29,0.19,0.16,0.14,0.09,0.04};
   double sy4[]={0.13,0.17,0.17,0.20,0.14,0.12,0.09,0.09,0.06,0.06,
		 0.05,0.03,0.02,0.02,0.01,0.03,0.02,0.01,0.01,0.02,
		 0.04,0.01,0.01,0.01,0.00,0.02};
  total=0.0;
  for(unsigned int ix=0;ix<26;++ix) {
    total+=y4[ix]*(x4[ix+1]-x4[ix]);
    st4[ix] = sqrt(sqr(st4[ix])+sqr(sy4[ix]));
  }
  bins  = vector<double>( x4, x4+27);
  data  = vector<double>( y4, y4+26);
  error = vector<double>(st4,st4+26);
  for(unsigned int ix=0;ix<data.size();++ix) {
    data[ix]/=total;
    error[ix]/=total;
  }
  _spectrumKpB = new_ptr(Histogram(bins,data,error));
  // proton excluding lambda products
  double x5[]={0.300,0.400,0.500,0.600,0.700,0.800,0.900,1.000,1.100,1.200,
	       1.300,1.400,1.500,1.600,1.700,1.800,1.900,2.000,2.100,2.200};
  double y5[]={ 0.117, 0.115, 0.201, 0.120, 0.125, 0.088, 0.112, 0.054, 0.060,-0.012,
		0.010,-0.021,-0.021,-0.011,-0.016, 0.034,-0.020, 0.031,-0.003};
 double st5[]={0.025,0.032,0.033,0.035,0.034,0.033,0.031,0.031,0.030,0.022,
	       0.022,0.023,0.025,0.027,0.029,0.032,0.033,0.038,0.038};
 double sy5[]={0.012,0.002,0.005,0.008,0.010,0.011,0.001,0.009,0.003,0.019,
	       0.029,0.013,0.014,0.013,0.014,0.016,0.027,0.015,0.025};
  total=0.0;
  for(unsigned int ix=0;ix<19;++ix) {
    total+=y5[ix]*(x5[ix+1]-x5[ix]);
    st5[ix] = sqrt(sqr(st5[ix])+sqr(sy5[ix]));
  }
  bins  = vector<double>( x5, x5+20);
  data  = vector<double>( y5, y5+19);
  error = vector<double>(st5,st5+19);
  for(unsigned int ix=0;ix<data.size();++ix) {
    data[ix]/=total;
    error[ix]/=total;
  }
  _spectrumprotonNo = new_ptr(Histogram(bins,data,error));
 // proton including lambda products
  double y6[]={ 0.172, 0.188, 0.274, 0.190, 0.175, 0.124, 0.136, 0.077, 0.072,-0.004,
		0.014,-0.018,-0.020,-0.010,-0.016, 0.034,-0.020, 0.031,-0.003};
  double st6[]={0.025,0.032,0.033,0.035,0.034,0.033,0.031,0.031,0.030,0.022,
		0.022,0.023,0.025,0.027,0.029,0.032,0.033,0.038,0.038};
  double sy6[]={0.017,0.003,0.007,0.013,0.013,0.016,0.001,0.012,0.003,0.019,
		0.040,0.013,0.014,0.013,0.014,0.016,0.027,0.015,0.025};
  total=0.0;
  for(unsigned int ix=0;ix<19;++ix) {
    total+=y6[ix]*(x5[ix+1]-x5[ix]);
    st6[ix] = sqrt(sqr(st6[ix])+sqr(sy6[ix]));
  }
  bins  = vector<double>( x5, x5+20);
  data  = vector<double>( y6, y6+19);
  error = vector<double>(st6,st6+19);
  for(unsigned int ix=0;ix<data.size();++ix) {
    data[ix]/=total;
    error[ix]/=total;
  }
 _spectrumprotonAll = new_ptr(Histogram(bins,data,error));
 // neutral kaons
 double x7[]={0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,
	      1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.7};
 double y7[]={26.1E-03,37.2E-03,41.7E-03,39.0E-03,33.3E-03,
	      28.8E-03,23.4E-03,17.9E-03,13.3E-03,10.4E-03,
	       9.4E-03, 6.9E-03, 4.6E-03, 2.9E-03, 3.2E-03,
	       1.6E-03, 1.4E-03, 1.0E-03, 2.0E-03, 0.7E-03};
  double st7[]={1.6E-03,1.8E-03,1.7E-03,1.6E-03,1.4E-03,
		1.3E-03,1.2E-03,1.1E-03,0.9E-03,0.9E-03,
		0.7E-03,0.7E-03,0.6E-03,0.6E-03,0.5E-03,
		0.5E-03,0.5E-03,0.4E-03,0.6E-03,1.0E-03};
 double sy7[]={1.8E-03,2.5E-03,2.8E-03,2.6E-03,2.2E-03,
	       1.9E-03,1.6E-03,1.2E-03,0.9E-03,0.7E-03,
	       0.7E-03,0.5E-03,0.4E-03,0.3E-03,0.3E-03,
	       0.2E-03,0.2E-03,0.1E-03,0.3E-03,0.6E-03};
  total=0.0;
  for(unsigned int ix=0;ix<20;++ix) {
    y7[ix]/=(x7[ix+1]-x7[ix]);
    total+=y7[ix]*(x7[ix+1]-x7[ix]);
    st7[ix] = sqrt(sqr(st7[ix])+sqr(sy7[ix]));
    st7[ix]/=(x7[ix+1]-x7[ix]);
  }
  bins  = vector<double>( x7, x7+21);
  data  = vector<double>( y7, y7+20);
  error = vector<double>(st7,st7+20);
  for(unsigned int ix=0;ix<data.size();++ix) {
    data[ix]/=total;
    error[ix]/=total;
  }
  _spectrumK0 = new_ptr(Histogram(bins,data,error));
  // pi0 spectrum
 double x8[]={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
	      1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,
	      2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0};
 double y8[]={4.483, 8.893, 8.106, 6.635, 5.245, 3.749, 2.780, 1.942, 1.370, 1.017,
	      0.701, 0.535, 0.402, 0.294, 0.223, 0.152, 0.129, 0.101, 0.084, 0.061,
	      0.059, 0.044, 0.027,-0.003, 0.002, 0.005,-0.003, 0.002,-0.003,-0.001};
 double sy8[]={0.143,0.297,0.182,0.117,0.083,0.062,0.048,0.039,0.031,0.026,
	       0.018,0.015,0.013,0.012,0.011,0.009,0.008,0.007,0.006,0.006,
	       0.005,0.005,0.005,0.005,0.004,0.004,0.004,0.004,0.004,0.004};
 double st8[]={0.644,1.182,1.012,0.836,0.699,0.506,0.377,0.290,0.209,0.161,
	       0.125,0.095,0.076,0.059,0.049,0.039,0.031,0.026,0.022,0.018,
	       0.017,0.014,0.011,0.011,0.009,0.007,0.007,0.006,0.005,0.005};
  total=0.0;
  for(unsigned int ix=0;ix<30;++ix) {
    total+=y8[ix]*(x8[ix+1]-x8[ix]);
    st8[ix] = sqrt(sqr(st8[ix])+sqr(sy8[ix]));
  }
  bins  = vector<double>( x8, x8+31);
  data  = vector<double>( y8, y8+30);
  error  = vector<double>(st8,st8+30);
  for(unsigned int ix=0;ix<data.size();++ix) {
    data[ix]/=total;
    error[ix]/=total;
  }
 _spectrumpi0 = new_ptr(Histogram(bins,data,error));
}




 


