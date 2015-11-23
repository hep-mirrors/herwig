// -*- C++ -*-
//
// SingleParticleAnalysis.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SingleParticleAnalysis class.
//

#include "SingleParticleAnalysis.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void SingleParticleAnalysis::analyze(const tPVector & particles) {
  useMe();
  int Ncharged = 0; 
  for(unsigned int ix=0;ix<particles.size();++ix) {
    const Lorentz5Momentum p = particles[ix]->momentum();
    if(!particles[ix]->data().charged()) continue;
    Ncharged++;
    *_yT += abs(_shapes->yT(p));
    *_yS += abs(_shapes->yS(p));
    *_ptinT += _shapes->ptInT(p)/GeV;
    *_ptoutT += _shapes->ptOutT(p)/GeV;
    *_ptinS += _shapes->ptInS(p)/GeV;
    *_ptoutS += _shapes->ptOutS(p)/GeV;
  }
  // make sure the right bin is booked by subtracting a small amount. 
  *_nch += Ncharged-0.00001;
}

void SingleParticleAnalysis::persistentOutput(PersistentOStream & os) const {
  os << _shapes;
}

void SingleParticleAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> _shapes;
}

ClassDescription<SingleParticleAnalysis> 
SingleParticleAnalysis::initSingleParticleAnalysis;
// Definition of the static class description member.

void SingleParticleAnalysis::Init() {

  static ClassDocumentation<SingleParticleAnalysis> documentation
    ("LEP SingleParticle analysis class",
     "The LEP SingleParticle analysis uses data from \\cite{Abreu:1996na}.",
     "%\\cite{Abreu:1996na}\n"
     "\\bibitem{Abreu:1996na}\n"
     "  P.~Abreu {\\it et al.}  [DELPHI Collaboration],\n"
     "   ``Tuning and test of fragmentation models based on identified particles  and\n"
     "  %precision event shape data,''\n"
     "  Z.\\ Phys.\\  C {\\bf 73}, 11 (1996).\n"
     "  %%CITATION = ZEPYA,C73,11;%%\n"
     );

  static Reference<SingleParticleAnalysis,EventShapes> interfaceEventShapes
    ("EventShapes",
     "Pointer to the object which calculates the event shapes",
     &SingleParticleAnalysis::_shapes, false, false, true, false, false);
}


void SingleParticleAnalysis::dofinish() {
  useMe();
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  // normalize the data
  _yT->normaliseToData();
  _yS->normaliseToData();
  _ptinT->normaliseToData();
  _ptoutT->normaliseToData();
  _ptinS->normaliseToData();
  _ptoutS->normaliseToData();
  _nch->prefactor(2.);
  // chisq
  double chisq,minfrac=0.05;
  unsigned int npoint;
  _yT->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for DELPHI rapidity wrt "
		     << "thrust axis distribution or " 
		     << chisq/npoint << " per degree of freedom \n";  
  _yS->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for DELPHI rapidity wrt "
		     << "sphericity axis distribution or " 
		     << chisq/npoint << " per degree of freedom \n";  
  _ptinT->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for DELPHI pT in the plane wrt "
		     << "thrust axis distribution " 
		     << chisq/npoint << " per degree of freedom \n";  
  _ptoutT->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for DELPHI pT out of the plane wrt "
		     << "thrust axis distribution " 
		     << chisq/npoint << " per degree of freedom \n";  
  _ptinS->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for DELPHI pT in the plane wrt "
		     << "sphericity axis distribution " 
		     << chisq/npoint << " per degree of freedom \n";  
  _ptoutS->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for DELPHI pT out of the wrt "
		     << "sphericity axis  plane distribution " 
		     << chisq/npoint << " per degree of freedom \n";  
  _nch->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL  number of charged "
		     << "particles distribution or " 
		     << chisq/npoint << " per degree of freedom \n";  
  // output the plots
  using namespace HistogramOptions;
  _yT->topdrawOutput(output,Frame|Errorbars|Ylog,
		     "RED",
		     "Radidity with respect to the thrust axis compared to DELPHI data",
		     "                                                                ",
		     "1/NdN/dy0T1",
		     "        X X",
		     "y0T1",
		     " X X");
  _yS->topdrawOutput(output,Frame|Errorbars|Ylog,
		     "RED",
		     "Radidity with respect to the sphericity axis"
		     " compared to DELPHI data",
		     "                                            "
		     "                        ",
		     "1/NdN/dy0S1",
		     "        X X",
		     "y0S1",
		     " X X");
  _ptinT->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"p0T1 in the plane with respect to the thrust axis"
			" compared to DELPHI data",
			" X X                                             "
			"                        ",
			"1/NdN/dp0T,in1",
			"        X    X",
			"p0T,in1",
			" X    X");
  _ptoutT->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"p0T1 out of the plane with respect to the thrust axis"
			" compared to DELPHI data",
			" X X                                                 "
			"                        ",
			"1/NdN/dp0T,out1",
			"        X     X",
			"p0T,out1",
			" X     X");
  _ptinS->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"p0T1 in the plane with respect to the sphericity axis"
			" compared to DELPHI data",
			" X X                                                 "
			"                        ",
			"1/NdN/dp0T,in1",
			"        X    X",
			"p0T,in1",
			" X    X");
  _ptoutS->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"p0T1 out of the plane with respect to the sphericity axis"
			" compared to DELPHI data",
			" X X                                                     ",
			"1/NdN/dp0T,out1",
			"        X     X",
			"p0T,out1",
			" X     X");
  _nch->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "Multiplcity of charged particles compared to OPAL data",
		      "                                                      ",
		      "1/SdS/dN0charged1",
		      "  G G   X       X",
		      "N0charged1",
		      " X       X");
}

void SingleParticleAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  vector<double> bins,data,error;
  // rapidity with respect to thrust axis
  double vals19[] = {0.000, 0.250, 0.500, 0.750, 1.000, 
		     1.250, 1.500, 1.750, 2.000, 2.250, 
		     2.500, 2.750, 3.000, 3.250, 3.500, 
		     3.750, 4.000, 4.250, 4.500, 5.000, 
		     5.500, 6.000}; 
  double data19[]={5.9517,6.4291,6.6831,6.7763,6.7650,
		   6.7230,6.6085,6.4346,6.1697,5.7692,
		   5.1450,4.3511,3.4481,2.5852,1.7999,
		   1.1669,0.7054,0.3997,0.15673,0.03374,
		   0.00502};
  double error19stat[]={0.0095,0.0095,0.0094,0.0089 ,0.0085,
			0.0083,0.0081,0.0080,0.0078 ,0.0076,
			0.0072,0.0066,0.0059,0.0051 ,0.0043,
			0.0035,0.0027,0.0020,0.00089,0.00041,
			0.00016};
  double error19syst[]={0.5628 ,0.4417 ,0.3319 ,0.2429 ,0.1755 ,
			0.1277 ,0.0939 ,0.0710 ,0.0617 ,0.0577 ,
			0.0514 ,0.0435 ,0.0345 ,0.0259 ,0.0180 ,
			0.0117 ,0.0071 ,0.0041 ,0.00177,0.00043,
			0.00007};
  double error19[21];
  for(unsigned int ix=0;ix<21;++ix){error19[ix]=sqrt(sqr(error19stat[ix])+
						     sqr(error19syst[ix]));}
  bins  = vector<double>(vals19 ,vals19 +22);
  data  = vector<double>(data19 ,data19 +21);
  error = vector<double>(error19,error19+21);
  _yT= new_ptr(Histogram(bins,data,error));
  // rapidity with respect to the sphericity axis
  double vals20[] = {0.000, 0.250, 0.500, 0.750, 1.000, 
		     1.250, 1.500, 1.750, 2.000, 2.250, 
		     2.500, 2.750, 3.000, 3.250, 3.500, 
		     3.750, 4.000, 4.250, 4.500, 5.000, 
		     5.500, 6.000};
  double data20[]={6.5680 ,6.5901 ,6.6094 ,6.6152 ,6.5917 ,
		   6.5817 ,6.5221 ,6.4097 ,6.1741 ,5.7542 ,
		   5.1066 ,4.2721 ,3.3718 ,2.5185 ,1.7588 ,
		   1.1589 ,0.7327 ,0.4402 ,0.1952 ,0.05574,
		   0.01306};
  double error20stat[]={0.0097   ,  0.0096   ,  0.0093   ,  0.0088   ,  0.0084   ,
			0.0082   ,  0.0081   ,  0.0080   ,  0.0079   ,  0.0076   ,
			0.0072   ,  0.0066   ,  0.0059   ,  0.0051   ,  0.0042   ,
			0.0034   ,  0.0027   ,  0.0021   ,  0.0010   ,  0.00050  ,
			0.00024};
  double error20syst[]={0.5323 ,0.4246 ,0.3329 ,0.2554 ,0.1908 ,
			0.1393 ,0.0983 ,0.0673 ,0.0617 ,0.0575 ,
			0.0511 ,0.0427 ,0.0337 ,0.0252 ,0.0176 ,
			0.0130 ,0.0105 ,0.0078 ,0.0046 ,0.00180,
			0.00055};
  double error20[21];
  for(unsigned int ix=0;ix<21;++ix){error20[ix]=sqrt(sqr(error20stat[ix])+
						     sqr(error20syst[ix]));}
  bins  = vector<double>(vals20 ,vals20 +22);
  data  = vector<double>(data20 ,data20 +21);
  error = vector<double>(error20,error20+21);
  _yS= new_ptr(Histogram(bins,data,error));
  // pt_in with respect to thrust axis
  double vals21[] = {0.000,  0.100,  0.200,  0.300,  0.400,  
		     0.500,  0.600,  0.700,  0.800,  1.000,  
		     1.200,  1.400,  1.600,  1.800,  2.000,  
		     2.500,  3.000,  3.500,  4.000,  5.000,  
		     6.000,  7.000,  8.000, 10.000, 12.000, 
		     14.000}; 
  double data21[]={46.663    ,39.823    ,29.351    ,21.034    ,15.156    ,
		   11.149    , 8.348    , 6.430    , 4.5131   , 2.9522   ,
		   2.0401   , 1.4597   , 1.0796   , 0.8155   , 0.5326   ,
		   0.2988   , 0.18067  , 0.11471  , 0.06305  , 0.03040  ,
		   0.01501  , 0.00858  , 0.00376  , 0.00123  , 0.00044};
  double error21stat[]={0.037    ,0.033    ,0.028    ,0.024    ,0.020    ,
			0.017    ,0.015    ,0.013    ,0.0076   ,0.0062   ,
			0.0052   ,0.0044   ,0.0038   ,0.0033   ,0.0017   ,
			0.0013   ,0.00099  ,0.00079  ,0.00042  ,0.00029  ,
			0.00021  ,0.00016  ,0.00008  ,0.00004  ,0.00003};
  double error21syst[]={1.758  ,1.092  ,0.608  ,0.350  ,0.219  ,
			0.150  ,0.111  ,0.087  ,0.0624 ,0.0420 ,
			0.0299 ,0.0222 ,0.0171 ,0.0136 ,0.0095 ,
			0.0057 ,0.00383,0.00273,0.00171,0.00095,
			0.00054,0.00035,0.00017,0.00006,0.00002};
  double error21[25];
  for(unsigned int ix=0;ix<25;++ix){error21[ix]=sqrt(sqr(error21stat[ix])+
						     sqr(error21syst[ix]));}
  bins  = vector<double>(vals21 ,vals21 +26);
  data  = vector<double>(data21 ,data21 +25);
  error = vector<double>(error21,error21+25);
  _ptinT= new_ptr(Histogram(bins,data,error));
  // pt_out  with respect to thrust axis
  double vals22[] = {0.000, 0.100, 0.200, 0.300, 0.400, 
		     0.500, 0.600, 0.700, 0.800, 1.000, 
		     1.200, 1.400, 1.600, 1.800, 2.000, 
		     2.500, 3.000, 3.500}; 
  double data22[]={66.160    ,49.794    ,33.544    ,21.407    ,13.466    ,
		   8.527    , 5.448    , 3.5845   , 2.0309   , 0.9959   ,
		   0.5288   , 0.2987   , 0.1755   , 0.1086   , 0.05266  ,
		   0.01885  , 0.00814};
  double error22stat[]={0.043    ,0.037    ,0.030    ,0.024    ,0.019    ,
			0.015    ,0.012    ,0.0098   ,0.0052   ,0.0037   ,
			0.0028   ,0.0021   ,0.0016   ,0.0013   ,0.00058  ,
			0.00035  ,0.00023};
  double error22syst[]={1.822  ,1.149  ,0.678  ,0.397  ,0.239  ,
			0.150  ,0.097  ,0.0658 ,0.0398 ,0.0216 ,
			0.0127 ,0.0079 ,0.0051 ,0.0034 ,0.00189,
			0.00080,0.00040};
  double error22[17];
  for(unsigned int ix=0;ix<17;++ix){error22[ix]=sqrt(sqr(error22stat[ix])+
						     sqr(error22syst[ix]));}
  bins  = vector<double>(vals22 ,vals22 +18);
  data  = vector<double>(data22 ,data22 +17);
  error = vector<double>(error22,error22+17);
  _ptoutT= new_ptr(Histogram(bins,data,error));
  // pt_in with respect to sphericity axis 
  double vals23[] = {0.000,  0.100,  0.200,  0.300,  0.400,  
		     0.500,  0.600,  0.700,  0.800,  1.000,  
		     1.200,  1.400,  1.600,  1.800,  2.000,  
		     2.500,  3.000,  3.500,  4.000,  5.000,  
		     6.000,  7.000,  8.000, 10.000, 12.000,  
		     14.000};
  double data23[]={49.206    ,38.461    ,28.203    ,20.391    ,14.926    ,
		   11.133    , 8.458    , 6.548    , 4.6706   , 3.0684   ,
		   2.1299   , 1.5201   , 1.1143   , 0.8398   , 0.5334   ,
		   0.2968   , 0.17343  , 0.10741  , 0.05615  , 0.02473  ,
		   0.01157  , 0.00561  , 0.00204  , 0.00049  , 0.00012 };
  double error23stat[]={0.038    ,0.033    ,0.027    ,0.023    ,0.020    ,
			0.017    ,0.015    ,0.013    ,0.0078   ,0.0064   ,
			0.0053   ,0.0045   ,0.0039   ,0.0034   ,0.0017   ,
			0.0013   ,0.00098  ,0.00078  ,0.00040  ,0.00027  ,
			0.00019  ,0.00013  ,0.00006  ,0.00003  ,0.00001};
  double error23syst[]={1.672  ,0.984  ,0.571  ,0.349  ,0.233  ,
			0.168  ,0.129  ,0.102  ,0.0747 ,0.0504 ,
			0.0359 ,0.0264 ,0.0201 ,0.0159 ,0.0107 ,
			0.0065 ,0.00418,0.00292,0.00176,0.00089,
			0.00048,0.00026,0.00010,0.00002,0.00001};
  double error23[25];
  for(unsigned int ix=0;ix<25;++ix){error23[ix]=sqrt(sqr(error23stat[ix])+
						     sqr(error23syst[ix]));}
  bins  = vector<double>(vals23 ,vals23 +26);
  data  = vector<double>(data23 ,data23 +25);
  error = vector<double>(error23,error23+25);
  _ptinS= new_ptr(Histogram(bins,data,error));
  // pt_out with respect to sphericity axis 
  double vals24[] = {0.000, 0.100, 0.200, 0.300, 0.400, 
		     0.500, 0.600, 0.700, 0.800, 1.000, 
		     1.200, 1.400, 1.600, 1.800, 2.000, 
		     2.500, 3.000, 3.500};
  double data24[]={66.825    ,50.556    ,34.241    ,21.708    ,13.481    ,
		   8.314    , 5.180    , 3.2986   , 1.7559   , 0.8187   ,
		   0.4064   , 0.2175   , 0.1232   , 0.0712   , 0.03217  ,
		   0.01112  , 0.00387};
  double error24stat[]={0.043    ,0.037    ,0.030    ,0.024    ,0.019    ,
			0.015    ,0.012    ,0.0094   ,0.0049   ,0.0034   ,
			0.0024   ,0.0018   ,0.0014   ,0.0011   ,0.00047  ,
			0.00029  ,0.00017};
  double error24syst[]={1.506  ,1.102  ,0.726  ,0.451  ,0.277  ,
			0.170  ,0.106  ,0.0679 ,0.0370 ,0.0181 ,
			0.0096 ,0.0055 ,0.0034 ,0.0022 ,0.00115,
			0.00050,0.00021};
  double error24[17];
  for(unsigned int ix=0;ix<17;++ix){error24[ix]=sqrt(sqr(error24stat[ix])+
						     sqr(error24syst[ix]));}
  bins  = vector<double>(vals24 ,vals24 +18);
  data  = vector<double>(data24 ,data24 +17);
  error = vector<double>(error24,error24+17);
  _ptoutS= new_ptr(Histogram(bins,data,error));
  // number of charged particles
  double vals25[]={01., 3., 5., 7., 9.,
		   11.,13.,15.,17.,19.,
		   21.,23.,25.,27.,29.,
		   31.,33.,35.,37.,39.,
		   41.,43.,45.,47.,49.,
		   51.,53.,55.};
  double data25[]={0.0010 ,0.016  ,0.16   ,0.68   , 2.08  ,
		   4.69  , 8.00  ,10.79  ,12.61  ,12.85  ,
		   11.83  , 9.99  , 7.85  , 5.95, 4.35  ,
		   2.97  , 2.02  , 1.29  , 0.81  , 0.47  ,
		   0.26  , 0.17  , 0.089 , 0.042 , 0.025 ,
		   0.011 , 0.004};
  double error25stat[]={0.0010  ,0.020   ,0.03   ,0.05   , 0.08  ,
			0.12  , 0.16  , 0.18  , 0.19  , 0.20  ,
			0.19  , 0.17  , 0.15  , 0.13  , 0.11  ,
			0.09  , 0.08  , 0.06  , 0.05  , 0.04  ,
			0.03  , 0.02  , 0.016 , 0.011 , 0.009 ,
			0.007 , 0.004};
  double error25syst[]={0.0,0.0,0.10, 0.18,  0.19,
			0.23,  0.19,  0.40,  0.24,  0.34,
			0.20,  0.35,  0.14,  0.17,  0.17,
			0.09,  0.09,  0.11,  0.06,  0.05,
			0.04,  0.05, 0.038, 0.020,  0.015,
			0.007,  0.004};
  double error25[27];
  for(unsigned int ix=0;ix<27;++ix)
    {
      error25[ix]=sqrt(sqr(error25stat[ix])+
		       sqr(error25syst[ix]))/100.;
      data25[ix]/=100.;
    }
  bins  = vector<double>(vals25 ,vals25 +28);
  data  = vector<double>(data25 ,data25 +27);
  error = vector<double>(error25,error25+27);
  _nch=new_ptr(Histogram(bins,data,error));
}
