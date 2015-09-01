// -*- C++ -*-
//
// LEPEventShapes.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LEPEventShapes class.
//

#include "LEPEventShapes.h"
#include "EventShapes.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;

void LEPEventShapes::analyze(tEventPtr event, long ieve,
			     int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  if ( loop > 0 || state != 0 || !event ) return;
  // get the final-state particles
  tPVector hadrons=event->getFinalState();
  // event shapes
}

LorentzRotation LEPEventShapes::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void LEPEventShapes::analyze(const tPVector & ) {
  double eventweight = generator()->currentEvent()->weight();
  _omthr ->addWeighted( 1.-_shapes->thrust() ,eventweight);
  _maj ->addWeighted( _shapes->thrustMajor() ,eventweight);
  _min ->addWeighted( _shapes->thrustMinor() ,eventweight);
  _obl ->addWeighted( _shapes->oblateness() ,eventweight); 
  _c ->addWeighted( _shapes->CParameter() ,eventweight); 
  _d ->addWeighted( _shapes->DParameter() ,eventweight); 
  _sph ->addWeighted( _shapes->sphericity() ,eventweight);
  _apl ->addWeighted( _shapes->aplanarity() ,eventweight);
  _pla ->addWeighted( _shapes->planarity() ,eventweight); 
  _mhi ->addWeighted( _shapes->Mhigh2() ,eventweight);
  _mlo ->addWeighted( _shapes->Mlow2() ,eventweight); 
  _mdiff ->addWeighted( _shapes->Mdiff2() ,eventweight); 
  _bmax ->addWeighted( _shapes->Bmax() ,eventweight); 
  _bmin ->addWeighted( _shapes->Bmin() ,eventweight); 
  _bsum ->addWeighted( _shapes->Bsum() ,eventweight); 
  _bdiff ->addWeighted( _shapes->Bdiff() ,eventweight); 
}

void LEPEventShapes::persistentOutput(PersistentOStream & os) const {
  os << _shapes;
}

void LEPEventShapes::persistentInput(PersistentIStream & is, int) {
  is >> _shapes;
}

ClassDescription<LEPEventShapes> LEPEventShapes::initLEPEventShapes;
// Definition of the static class description member.

void LEPEventShapes::Init() {

  static ClassDocumentation<LEPEventShapes> documentation
    ("The LEPEventShapes class compares event shapes at the Z mass"
     "with experimental results",
     "The LEP EventShapes analysis uses data from \\cite{Pfeifenschneider:1999rz,Abreu:1996na}.",
     "%\\cite{Pfeifenschneider:1999rz}\n"
     "\\bibitem{Pfeifenschneider:1999rz}\n"
     "  P.~Pfeifenschneider {\\it et al.}  [JADE collaboration and OPAL\n"
     "                  Collaboration],\n"
     "   ``QCD analyses and determinations of alpha(s) in e+ e- annihilation at\n"
     "  %energies between 35-GeV and 189-GeV,''\n"
     "  Eur.\\ Phys.\\ J.\\  C {\\bf 17}, 19 (2000)\n"
     "  [arXiv:hep-ex/0001055].\n"
     "  %%CITATION = EPHJA,C17,19;%%\n"
     "%\\cite{Abreu:1996na}\n"
     "\\bibitem{Abreu:1996na}\n"
     "  P.~Abreu {\\it et al.}  [DELPHI Collaboration],\n"
     "   ``Tuning and test of fragmentation models based on identified particles  and\n"
     "  %precision event shape data,''\n"
     "  Z.\\ Phys.\\  C {\\bf 73}, 11 (1996).\n"
     "  %%CITATION = ZEPYA,C73,11;%%\n"
     );

  static Reference<LEPEventShapes,EventShapes> interfaceEventShapes
    ("EventShapes",
     "Pointer to the object which calculates the event shapes",
     &LEPEventShapes::_shapes, false, false, true, false, false);

}

void LEPEventShapes::dofinish() {
  useMe();
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + 
    string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  _omthr->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"1-T compared to DELPHI data",
			"                           ",
			"1/SdS/d(1-T)",
			"  G G       ",
			"1-T",
			"   ");
  _maj->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "Thrust Major compared to DELPHI data",
		      "                                    ",
		      "1/SdS/dMajor",
		      "  G G       ",
		      "Major",
		      "     ");
  _min->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "Thrust Minor compared to DELPHI data",
		      "                           ",
		      "1/SdS/dMinor",
		      "  G G       ",
		      "Minor",
		      "     ");
  _obl->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "Oblateness compared to DELPHI data",
		      "                                  ",
		      "1/SdS/dO",
		      "  G G   ",
		      "O",
		      " ");
  _sph->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "Sphericity compared to DELPHI data",
		      "                                  ",
		      "1/SdS/dS",
		      "  G G   ",
		      "S",
		      " ");
  _apl->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "Aplanarity compared to DELPHI data",
		      "                                  ",
		      "1/SdS/dA",
		      "  G G   ",
		      "A",
		      " ");
  _pla->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "Planarity compared to DELPHI data",
		      "                                 ",
		      "1/SdS/dP",
		      "  G G   ",
		      "P",
		      " ");
  _c->topdrawOutput(output,Frame|Errorbars|Ylog,
		    "RED",
		    "C parameter compared to DELPHI data",
		    "                                   ",
		    "1/SdS/dC",
		    "  G G   ",
		    "C",
		    " ");
  _d->topdrawOutput(output,Frame|Errorbars|Ylog,
		    "RED",
		    "D parameter compared to DELPHI data",
		    "                                   ",
		    "1/SdS/dD",
		    "  G G   ",
		    "D",
		    " ");
  _mhi->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "High hemisphere mass compared to DELPHI data",
		      "                                            ",
		      "1/SdS/dM0high1",
		      "  G G   X    X",
		      "M0high1",
		      " X    X");
  _mlo->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "Low hemisphere mass compared to DELPHI data",
		      "                                           ",
		      "1/SdS/dM0low1",
		      "  G G   X   X",
		      "M0low1",
		      " X   X");
  _mdiff->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "Difference in hemisphere masses compared to DELPHI data",
		      "                                                       ",
		      "1/SdS/dM0diff1",
		      "  G G   X    X",
		      "M0diff1",
		      " X    X");
  _bmax->topdrawOutput(output,Frame|Errorbars|Ylog,
		       "RED",
		       "Wide jet broadening measure compared to DELPHI data",
		       "                                                   ",
		       "1/SdS/dB0max1",
		       "  G G   X   X",
		       "B0max1",
		       " X   X");
  _bmin->topdrawOutput(output,Frame|Errorbars|Ylog,
		       "RED",
		       "Narrow jet broadening measure compared to DELPHI data",
		       "                                                     ",
		       "1/SdS/dB0min1",
		       "  G G   X   X",
		       "B0min1",
		       " X   X");
  _bsum->topdrawOutput(output,Frame|Errorbars|Ylog,
		       "RED",
		       "Sum of jet broadening measures compared to DELPHI data",
		       "                                                      ",
		       "1/SdS/dB0sum1",
		       "  G G   X   X",
		       "B0sum1",
		       " X   X");
  _bdiff->topdrawOutput(output,Frame|Errorbars|Ylog,
		       "RED",
		       "Difference of jet broadenings measure compared to DELPHI data",
		       "                                                             ",
		       "1/SdS/dB0diff1",
		       "  G G   X    X",
		       "B0diff1",
		       " X    X");
  // chi squareds
  double chisq=0.,minfrac=0.05;
  unsigned int ndegrees;
  _omthr->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI thrust distribution\n";
  _maj->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI major distribution\n";
  _min->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI minor distribution\n";
  _obl->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI oblateness distribution\n";
  _sph->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI sphericity distribution\n";
  _apl->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI aplanarity distribution\n";
  _pla->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI planarity distribution\n";
  _c->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI C distribution\n";
  _d->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI D distribution\n";
  _mhi->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI m_high distribution\n";
  _mlo->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI m_low distribution\n";
  _mdiff->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI m_diff distribution\n";
  _bmax->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI B_max distribution\n";
  _bmin->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI B_min distribution\n";
  _bsum->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI B_sum distribution\n";
  _bdiff->chiSquared(chisq,ndegrees,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << ndegrees 
			  << " degrees of freedom for DELPHI B_diff distribution\n";
}

void LEPEventShapes::doinitrun() {
  AnalysisHandler::doinitrun();
  vector<double> bins,data,error;
  // 1-T
  double vals1[] = {0.000, 0.010, 0.020, 0.030, 0.040, 
		    0.050, 0.060, 0.070, 0.080, 0.090, 
		    0.100, 0.120, 0.140, 0.160, 0.180, 
		    0.200, 0.250, 0.300, 0.350, 0.400, 
		    0.500};
  double data1[]= { 1.030,10.951,17.645,14.192,10.009,
		    7.572, 5.760, 4.619, 3.792, 3.176, 
		    2.456, 1.825, 1.401, 1.074, 0.8262,
		    0.5525,0.3030,0.1312,0.0238,0.0007};
  double error1stat[]={0.019  ,0.051  ,0.066  ,0.061  ,0.050  ,
		       0.044  ,0.038  ,0.034  ,0.031  ,0.028  ,
		       0.018  ,0.015  ,0.013  ,0.011  ,0.0100 ,
		       0.0051 ,0.0038 ,0.0025 ,0.0012 ,0.0002 };
  double error1syst[]={0.076  , 0.527	, 0.547	, 0.292	, 0.152	,
		       0.101	, 0.076	, 0.062	, 0.051	, 0.042	,
		       0.032	, 0.022	, 0.016	, 0.011	, 0.0083,
		       0.0065	, 0.0058, 0.0044, 0.0014, 0.0001};
  double error1[20];
  for(unsigned int ix=0;ix<20;++ix){error1[ix]=sqrt(sqr(error1stat[ix])+
						    sqr(error1syst[ix]));}
  bins  = vector<double>(vals1 ,vals1 +21);
  data  = vector<double>(data1 ,data1 +20);
  error = vector<double>(error1,error1+20);
  _omthr= new_ptr(Histogram(bins,data,error));
  // major
  double vals2[] = {0.000, 0.020, 0.040, 0.050, 0.060, 
		    0.070, 0.080, 0.100, 0.120, 0.140, 
		    0.160, 0.200, 0.240, 0.280, 0.320, 
		    0.360, 0.400, 0.440, 0.480, 0.520, 
		    0.560, 0.600, 0.640};
  double data2[]={0.00040  ,0.0590   ,0.642    ,2.178    ,4.303 ,
		  5.849    ,6.889    ,6.342    ,4.890    ,3.900 ,
		  2.960    ,2.124    ,1.5562   ,1.1807   ,0.8693,
		  0.6493   ,0.4820   ,0.3493   ,0.2497   ,0.1489,
		  0.0714   ,0.0203};
  double error2stat[]={0.00090  ,0.0030   ,0.013    ,0.024    ,0.034    ,
		       0.039    ,0.030    ,0.028    ,0.024    ,0.021    ,
		       0.013    ,0.011    ,0.0095   ,0.0083   ,0.0071   ,
		       0.0061   ,0.0052   ,0.0044   ,0.0037   ,0.0028   ,
		       0.0019   ,0.0010};
  double error2syst[]={0.00005 ,0.0058  ,0.028	 ,0.086	 ,0.155	 ,
		       0.192   ,0.194	,0.143	 ,0.085	 ,0.050	 ,
		       0.030   ,0.021	,0.0156  ,0.0118 ,0.0087 ,
		       0.0065  ,0.0048  ,0.0055  ,0.0065 ,0.0058 ,
		       0.0038  ,0.0014};
  double error2[22];
  for(unsigned int ix=0;ix<22;++ix){error2[ix]=sqrt(sqr(error2stat[ix])+
						    sqr(error2syst[ix]));}
  bins  = vector<double>(vals2 ,vals2 +23);
  data  = vector<double>(data2 ,data2 +22);
  error = vector<double>(error2,error2+22);
  _maj= new_ptr(Histogram(bins,data,error));
  // minor
  double vals3[] = {0.000, 0.020, 0.040, 0.050, 0.060, 
		    0.070, 0.080, 0.100, 0.120, 0.140, 
		    0.160, 0.200, 0.240, 0.280, 0.320, 
		    0.400}; 
  double data3[]={ 0.0156  , 1.236   , 5.706   , 9.714   ,12.015   ,
		   12.437   ,10.404   , 6.918   , 4.250   , 2.517   ,
		   1.2561  , 0.4895  , 0.2112  , 0.0879  , 0.0250  };
  double error3stat[]={0.0017  ,0.013   ,0.037   ,0.048   ,0.054   ,
		     0.055   ,0.036   ,0.029   ,0.023   ,0.017   ,
		     0.0086  ,0.0054  ,0.0036  ,0.0023  ,0.0009};
  double error3syst[]={0.0036,0.066 ,0.073 ,0.125 ,0.155 ,
		     0.161 ,0.136 ,0.092 ,0.058 ,0.035 ,
		     0.0187,0.0080,0.0039,0.0018,0.0006};
  double error3[15];
  for(unsigned int ix=0;ix<15;++ix){error3[ix]=sqrt(sqr(error3stat[ix])+
						    sqr(error3syst[ix]));}
  bins  = vector<double>(vals3 ,vals3 +16);
  data  = vector<double>(data3 ,data3 +15);
  error = vector<double>(error3,error3+15);
  _min= new_ptr(Histogram(bins,data,error));
  // oblateness
  double vals4[] = {0.000, 0.020, 0.040, 0.060, 0.080, 
		    0.100, 0.120, 0.140, 0.160, 0.200, 
		    0.240, 0.280, 0.320, 0.360, 0.400, 
		    0.440, 0.520};
   double data4[]={ 9.357   ,11.508   , 7.215   , 4.736   , 3.477   ,
		    2.696   , 2.106   , 1.690   , 1.2648  , 0.8403  ,
		    0.5674  , 0.3842  , 0.2573  , 0.1594  , 0.0836  ,
		    0.0221  };
   double error4stat[]={0.036   ,0.038   ,0.029   ,0.023   ,0.020   ,
			0.018   ,0.016   ,0.014   ,0.0085  ,0.0069  ,
			0.0056  ,0.0046  ,0.0037  ,0.0029  ,0.0020  ,
			0.0007};
   double error4syst[]={0.178 ,0.140 ,0.072 ,0.047 ,0.035 ,
			0.027 ,0.021 ,0.017 ,0.0126,0.0087,
			0.0065,0.0050,0.0043,0.0037,0.0030,
			0.0015};
  double error4[16];
  for(unsigned int ix=0;ix<16;++ix){error4[ix]=sqrt(sqr(error4stat[ix])+
						    sqr(error4syst[ix]));}
   bins  = vector<double>(vals4 ,vals4 +17);
   data  = vector<double>(data4 ,data4 +16);
   error = vector<double>(error4,error4+16);
   _obl= new_ptr(Histogram(bins,data,error));
  // sphericity
  double vals5[] = {0.000, 0.010, 0.020, 0.030, 0.040, 
		    0.050, 0.060, 0.080, 0.100, 0.120, 
		    0.160, 0.200, 0.250, 0.300, 0.350, 
		    0.400, 0.500, 0.600, 0.700, 0.850};
  double data5[]={16.198    ,20.008    ,12.896    , 8.237    , 5.885    ,
		  4.458    , 3.272    , 2.290    , 1.699    , 1.2018   ,
		  0.7988   , 0.5610   , 0.3926   , 0.2810   , 0.2099   ,
		  0.1441   , 0.0842   , 0.04160  , 0.00758 };
  double error5stat[]={0.067    ,0.072    ,0.056    ,0.043    ,0.037    ,
		       0.032    ,0.019    ,0.016    ,0.014    ,0.0082   ,
		       0.0067   ,0.0050   ,0.0042   ,0.0035   ,0.0030   ,
		       0.0018   ,0.0013   ,0.00092  ,0.00032};
  double error5syst[]={0.208  ,0.246  ,0.153  ,0.094  ,0.065  ,
		       0.048  ,0.034  ,0.023  ,0.017  ,0.0120 ,
		       0.0080 ,0.0063 ,0.0051 ,0.0043 ,0.0037 ,
		       0.0032 ,0.0023 ,0.00129,0.00024};
  double error5[19];
  for(unsigned int ix=0;ix<19;++ix){error5[ix]=sqrt(sqr(error5stat[ix])+
						    sqr(error5syst[ix]));}
  bins  = vector<double>(vals5 ,vals5 +20);
  data  = vector<double>(data5 ,data5 +19);
  error = vector<double>(error5,error5+19);
  _sph=new_ptr(Histogram(bins,data,error));
  // aplanarity
  double vals6[] = {0.000, 0.005, 0.010, 0.015, 0.020, 
		    0.030, 0.040, 0.060, 0.080, 0.100, 
		    0.120, 0.140, 0.160, 0.200, 0.250, 
		    0.300};
  double data6[]={75.10     ,55.31     ,26.03     ,13.927    , 6.768    ,
		  3.014    , 1.281    , 0.5181   , 0.2619   , 0.1461   ,
		  0.0758   , 0.0467   , 0.0234   , 0.00884  , 0.00310  };
  double error6stat[]={0.19       ,0.17       ,0.11       ,0.079    ,0.038    ,
		       0.025      ,0.012    ,0.0075   ,0.0054   ,0.0041   ,
		       0.0029   ,0.0023   ,0.0011   ,0.00061  ,0.00040 };
  double error6syst[]={0.75	,0.55	,0.28	,0.176  ,0.098  ,
		       0.056  ,0.035  ,0.0188 ,0.0118 ,0.0079 ,
		       0.0043 ,0.0027 ,0.0014 ,0.00052,0.00018};
  double error6[15];
  for(unsigned int ix=0;ix<15;++ix){error6[ix]=sqrt(sqr(error6stat[ix])+
						    sqr(error6syst[ix]));}
  bins  = vector<double>(vals6 ,vals6 +16);
  data  = vector<double>(data6 ,data6 +15);
  error = vector<double>(error6,error6+15);
  _apl= new_ptr(Histogram(bins,data,error));
  // planarity
  double vals7[] = {0.000, 0.005, 0.010, 0.015, 0.020, 
		    0.025, 0.030, 0.035, 0.040, 0.050, 
		    0.060, 0.080, 0.100, 0.120, 0.160, 
		    0.200, 0.250, 0.300, 0.350, 0.400, 
		    0.500};
  double data7[]={68.69    ,31.66    ,17.091   ,11.370   , 8.417   ,
		  6.578   , 5.479   , 4.493   , 3.610   , 2.749   ,
		  1.987   , 1.362   , 1.008   , 0.6676  , 0.4248  ,
		  0.2692  , 0.1742  , 0.1042  , 0.0566  , 0.0145  };
  double error7stat[]={0.19    ,0.12    ,0.088   ,0.072   ,0.062   ,
		       0.055   ,0.050   ,0.045   ,0.029   ,0.025   ,
		       0.015   ,0.012   ,0.011   ,0.0061  ,0.0048  ,
		       0.0034  ,0.0028  ,0.0021  ,0.0015  ,0.0006};
  double error7syst[]={0.74  ,0.35  ,0.188 ,0.127 ,0.095 ,
		       0.075 ,0.063 ,0.052 ,0.042 ,0.033 ,
		       0.024 ,0.017 ,0.013 ,0.0093,0.0063,
		       0.0042 ,0.0029,0.0019,0.0011,0.0003};
  double error7[20];
  for(unsigned int ix=0;ix<20;++ix){error7[ix]=sqrt(sqr(error7stat[ix])+
						    sqr(error7syst[ix]));}
  bins  = vector<double>(vals7 ,vals7 +21);
  data  = vector<double>(data7 ,data7 +20);
  error = vector<double>(error7,error7+20);
  _pla= new_ptr(Histogram(bins,data,error));
  // C
  double vals8[] = {0.000, 0.040, 0.080, 0.120, 0.160, 
		    0.200, 0.240, 0.280, 0.320, 0.360, 
		    0.400, 0.440, 0.480, 0.520, 0.560, 
		    0.600, 0.640, 0.680, 0.720, 0.760, 
		    0.800, 0.840, 0.880, 0.920};
  double data8[]={0.0881  ,1.5383  ,3.909   ,3.833   ,2.835   ,
		  2.164   ,1.716   ,1.3860  ,1.1623  ,0.9720  ,
		  0.8349  ,0.7161  ,0.6205  ,0.5441  ,0.4844  ,
		  0.4209  ,0.3699  ,0.3286  ,0.2813  ,0.2178  ,
		  0.1287  ,0.0542  ,0.0212  };
  double error8stat[]={0.0030   ,0.0100   ,0.016    ,0.016    ,0.013    ,
		       0.012    ,0.010    ,0.0092   ,0.0084   ,0.0077   ,
		       0.0072   ,0.0066   ,0.0061   ,0.0057   ,0.0054   ,
		       0.0050   ,0.0046   ,0.0044   ,0.0040   ,0.0033   ,
		       0.0026   ,0.0016   ,0.0009};
  double error8syst[]={0.0067,0.0831,0.142 ,0.088 ,0.040 ,
		       0.022 ,0.017 ,0.0139,0.0116,0.0097,
		       0.0083,0.0072,0.0062,0.0054,0.0050,
		       0.0063,0.0079,0.0099,0.0129,0.0151,
		       0.0130,0.0076,0.0040};
  double error8[23];
  for(unsigned int ix=0;ix<23;++ix){error8[ix]=sqrt(sqr(error8stat[ix])+
						    sqr(error8syst[ix]));}
  bins  = vector<double>(vals8 ,vals8 +24);
  data  = vector<double>(data8 ,data8 +23);
  error = vector<double>(error8,error8+23);
  _c= new_ptr(Histogram(bins,data,error));
  // D
  double vals9[] = {0.000, 0.008, 0.016, 0.030, 0.044, 
		    0.066, 0.088, 0.112, 0.136, 0.162, 
		    0.188, 0.218, 0.248, 0.284, 0.320, 
		    0.360, 0.400, 0.450, 0.500, 0.560, 
		    0.620, 0.710, 0.800};
  double data9[]={22.228   ,22.766   ,12.107   , 6.879   , 4.284  ,
		  2.727   , 1.909   , 1.415   , 1.051   , 0.7977  , 
		  0.6155  , 0.4566  , 0.3341  , 0.2452  , 0.1774  ,
		  0.1234  , 0.0902  , 0.0603  , 0.0368  , 0.0222  ,
		  0.0128  , 0.0052  };
  double error9stat[]={0.082   ,0.085   ,0.047   ,0.035   ,0.022   ,
		       0.018   ,0.014   ,0.012   ,0.010   ,0.0089  ,
		       0.0073  ,0.0063  ,0.0049  ,0.0042  ,0.0033  ,
		       0.0028  ,0.0021  ,0.0017  ,0.0012  ,0.0009  ,
		       0.0006  ,0.0004};
  double error9syst[]={0.868 ,0.440 ,0.150 ,0.079 ,0.053 ,
		       0.036 ,0.028 ,0.022 ,0.018 ,0.0145,
		       0.0117   ,0.0089,0.0065,0.0049,0.0037  ,
		       0.0028,0.0023  ,0.0018,0.0013,0.0009,
		       0.0006,0.0003};
  double error9[22];
  for(unsigned int ix=0;ix<22;++ix){error9[ix]=sqrt(sqr(error9stat[ix])+
						    sqr(error9syst[ix]));}
  bins  = vector<double>(vals9 ,vals9 +23);
  data  = vector<double>(data9 ,data9 +22);
  error = vector<double>(error9,error9+22);
  _d= new_ptr(Histogram(bins,data,error));
  // M high
  double vals10[] = {0.000, 0.010, 0.020, 0.030, 0.040, 
		     0.050, 0.060, 0.080, 0.100, 0.120, 
		     0.140, 0.160, 0.200, 0.250, 0.300, 
		     0.350, 0.400};
  double data10[]={ 1.994   ,18.580   ,20.678   ,13.377   , 8.965   ,
		    6.558   , 4.515   , 2.914   , 1.991   , 1.406   ,
		    1.010   , 0.6319  , 0.3085  , 0.1115  , 0.0184  ,
		    0.0008   };
  double error10stat[]={0.027   ,0.065   ,0.076   ,0.060   ,0.049   ,
			0.041   ,0.024   ,0.019   ,0.016   ,0.013   ,
			0.011   ,0.0063  ,0.0039  ,0.0022  ,0.0008  ,
			0.0002 };
  double error10syst[]={0.166 ,0.709 ,0.729 ,0.412 ,0.239 ,
			0.151 ,0.082 ,0.037 ,0.020 ,0.014 ,
			0.010 ,0.0063,0.0051,0.0039,0.0012,
			0.0001};
  double error10[16];
  for(unsigned int ix=0;ix<16;++ix){error10[ix]=sqrt(sqr(error10stat[ix])+
						     sqr(error10syst[ix]));}
  bins  = vector<double>(vals10 ,vals10 +17);
  data  = vector<double>(data10 ,data10 +16);
  error = vector<double>(error10,error10+16);
  _mhi= new_ptr(Histogram(bins,data,error));
  // M low
  double vals11[] = {0.000, 0.010, 0.020, 0.030, 0.040, 
		     0.050, 0.060, 0.080, 0.100, 0.120};
  double data11[]={23.414   ,39.12    ,18.080   , 7.704   , 3.922   ,
		   2.128   , 1.013   , 0.3748  , 0.1412 };
  double error11stat[]={0.074   ,0.11      ,0.081    ,0.052   ,0.036   ,
			0.026   ,0.013  ,0.0079    ,0.0050};
  double error11syst[]={1.595 ,2.65  ,1.215 ,0.514 ,0.260 ,
			0.140 ,0.066 ,0.0241,0.0089};
  double error11[9];
  for(unsigned int ix=0;ix<9;++ix){error11[ix]=sqrt(sqr(error11stat[ix])+
						    sqr(error11syst[ix]));}
  bins  = vector<double>(vals11 ,vals11 +10);
  data  = vector<double>(data11 ,data11 + 9);
  error = vector<double>(error11,error11+ 9);
  _mlo= new_ptr(Histogram(bins,data,error));
  // M diff
  double vals12[] = {0.000, 0.010, 0.020, 0.030, 0.040, 
		     0.060, 0.080, 0.120, 0.160, 0.200, 
		     0.250, 0.300, 0.350, 0.400};
  double data12[]={35.393   ,20.745   ,11.426   , 7.170   , 4.344   ,
		   2.605   , 1.4238  , 0.7061  , 0.3831  , 0.1836  ,
		   0.0579  , 0.0075  , 0.0003};
  double error12stat[]={0.092   ,0.071   ,0.052   ,0.041   ,0.023   ,
			0.017   ,0.0092  ,0.0064  ,0.0046  ,0.0028  ,
			0.0015  ,0.0006  ,0.0002};
  double error12syst[]={0.354 ,0.207 ,0.114 ,0.072 ,0.043 ,
			0.026 ,0.0142,0.0071,0.0044,0.0032,
			0.0018,0.0006,0.0001};
  double error12[13];
  for(unsigned int ix=0;ix<13;++ix){error12[ix]=sqrt(sqr(error12stat[ix])+
						     sqr(error12syst[ix]));}
  bins  = vector<double>(vals12 ,vals12 +14);
  data  = vector<double>(data12 ,data12 +13);
  error = vector<double>(error12,error12+13);
  _mdiff= new_ptr(Histogram(bins,data,error));
  // Bmax
  double vals13[] = {0.010, 0.020, 0.030, 0.040, 0.050, 
		     0.060, 0.070, 0.080, 0.100, 0.120, 
		     0.140, 0.170, 0.200, 0.240, 0.280, 
		     0.320};
  double data13[]={0.6707  , 7.538   ,14.690   ,13.942   ,11.298   ,
		   9.065   , 7.387   , 5.445   , 3.796   , 2.670   ,
		   1.756   , 1.0580  , 0.5288  , 0.1460  , 0.0029  };
  double error13stat[]={0.0096  ,0.038   ,0.058   ,0.057   ,0.053   ,
			0.048   ,0.043   ,0.026   ,0.022   ,0.018   ,
			0.012   ,0.0092  ,0.0056  ,0.0028  ,0.0004};
  double error13syst[]={0.1077,0.809 ,0.745 ,0.592 ,0.379 ,
			0.266 ,0.222 ,0.176 ,0.127 ,0.087 ,
			0.051 ,0.0218,0.0053,0.0071,0.0003};
  double error13[15];
  for(unsigned int ix=0;ix<15;++ix){error13[ix]=sqrt(sqr(error13stat[ix])+
						     sqr(error13syst[ix]));}
  bins  = vector<double>(vals13 ,vals13 +16);
  data  = vector<double>(data13 ,data13 +15);
  error = vector<double>(error13,error13+15);
  _bmax= new_ptr(Histogram(bins,data,error));
  // Bmin
  double vals14[] = {0.000, 0.010, 0.020, 0.030, 0.040, 
		     0.050, 0.060, 0.080, 0.100, 0.120, 
		     0.150, 0.180}; 
  double data14[]={0.645   ,11.169   ,28.908   ,25.972   ,14.119   ,
		   7.500   , 3.405   , 1.320   , 0.5448  , 0.1916  ,
		   0.0366};
  double error14stat[]={0.010   ,0.045   ,0.082   ,0.083   ,0.061   ,
			0.044   ,0.021   ,0.013   ,0.0082  ,0.0040  ,
			0.0017};
  double error14syst[]={0.096 ,1.006 ,1.823 ,1.478 ,0.860 ,
			0.494 ,0.233 ,0.089 ,0.0328,0.0104,
			0.0034};
  double error14[11];
  for(unsigned int ix=0;ix<11;++ix){error14[ix]=sqrt(sqr(error14stat[ix])+
						     sqr(error14syst[ix]));}
  bins  = vector<double>(vals14 ,vals14 +12);
  data  = vector<double>(data14 ,data14 +11);
  error = vector<double>(error14,error14+11);
  _bmin= new_ptr(Histogram(bins,data,error));
  // Bsum
  double vals15[] = {0.020, 0.030, 0.040, 0.050, 0.060, 
		     0.070, 0.080, 0.090, 0.100, 0.110, 
		     0.130, 0.150, 0.170, 0.190, 0.210, 
		     0.240, 0.270, 0.300, 0.330, 0.360};
  double data15[]={0.2030  ,1.628   ,4.999   ,8.190   ,9.887   ,
		   9.883   ,9.007   ,7.746   ,6.714   ,5.393   ,
		   3.998   ,2.980   ,2.294   ,1.747   ,1.242   ,
		   0.8125  ,0.4974  ,0.2285  ,0.0732  };
  double error15stat[]={0.0055  ,0.015   ,0.031   ,0.041   ,0.047   ,
			0.049   ,0.047   ,0.044   ,0.041   ,0.026   ,
			0.023   ,0.019   ,0.017   ,0.015   ,0.010   ,
			0.0080  ,0.0062  ,0.0041  ,0.0024};
  double error15syst[]={0.0383,0.183 ,0.463 ,0.644 ,0.661 ,
			0.564 ,0.443 ,0.332 ,0.255 ,0.180 ,
			0.125 ,0.098 ,0.085 ,0.075 ,0.063 ,
			0.0469,0.0296,0.0119,0.0007};
  double error15[19];
  for(unsigned int ix=0;ix<19;++ix){error15[ix]=sqrt(sqr(error15stat[ix])+
						     sqr(error15syst[ix]));}
  bins  = vector<double>(vals15 ,vals15 +20);
  data  = vector<double>(data15 ,data15 +19);
  error = vector<double>(error15,error15+19);
  _bsum= new_ptr(Histogram(bins,data,error));
  // Bdiff
  double vals16[] = {0.000, 0.010, 0.020, 0.030, 0.040, 
		     0.050, 0.060, 0.070, 0.080, 0.090, 
		     0.100, 0.120, 0.140, 0.160, 0.180, 
		     0.200, 0.240, 0.280}; 
  double data16[]={26.630   ,18.684   ,12.343   , 8.819   , 6.688   ,
		   5.111   , 4.071   , 3.271   , 2.681   , 2.233   ,
		   1.647   , 1.111   , 0.7618  , 0.5138  , 0.3167  ,
		   0.1265  , 0.0117};
  double error16stat[]={0.081   ,0.066     ,0.054   ,0.046   ,0.040 ,
			0.035   ,0.031     ,0.028   ,0.025   ,0.023 ,
			0.014   ,0.011   ,0.0095  ,0.0078  ,0.0062  ,
			0.0026  ,0.0008};
  double error16syst[]={0.459 ,0.292 ,0.186 ,0.134 ,0.106 ,
			0.084 ,0.068 ,0.054 ,0.043 ,0.035 ,
			0.026 ,0.019 ,0.0144,0.0119,0.0098,
			0.0056,0.0008};
  double error16[17];
  for(unsigned int ix=0;ix<17;++ix){error16[ix]=sqrt(sqr(error16stat[ix])+
						     sqr(error16syst[ix]));}
  bins  = vector<double>(vals16 ,vals16 +18);
  data  = vector<double>(data16 ,data16 +17);
  error = vector<double>(error16,error16+17);
  _bdiff= new_ptr(Histogram(bins,data,error));
}
