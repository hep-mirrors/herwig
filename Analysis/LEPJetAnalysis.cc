// -*- C++ -*-
//
// LEPJetAnalysis.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LEPJetAnalysis class.
//

#include "LEPJetAnalysis.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig/Utilities/HerwigStrategy.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace Herwig;

void LEPJetAnalysis::analyze(tEventPtr event, long, int, int ) {
  ++_nevent;
  tPVector particles;
  event->selectFinalState(back_inserter(particles));

  //  copy fastjet particles from event record.  Templated fastjet
  //  method might leave units ambigouos.  Loop with integer index
  //  allows backtracing ThePEG particles if needed.
  vector<fastjet::PseudoJet> fastjet_particles;

  for (unsigned int j=0; j<particles.size(); j++) {
    fastjet::PseudoJet p(particles[j]->momentum().x()/GeV, 
			 particles[j]->momentum().y()/GeV, 
			 particles[j]->momentum().z()/GeV, 
			 particles[j]->momentum().e()/GeV);
    p.set_user_index(j);
    fastjet_particles.push_back(p);
  }
  
  fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;
  fastjet::Strategy strategy = fastjet::Best;
  fastjet::JetDefinition jet_def(fastjet::ee_kt_algorithm, 
				 recomb_scheme, strategy);
  fastjet::ClusterSequence cs(fastjet_particles, jet_def);
  
  // ynm distributions
  *_y23 += cs.exclusive_ymerge(2); 
  *_y34 += cs.exclusive_ymerge(3); 
  *_y45 += cs.exclusive_ymerge(4); 
  *_y56 += cs.exclusive_ymerge(5); 

  const int entr_fr = 37;
  const int ddentr = 16;
  for(int j = 0; j<entr_fr; j++) {
    int Nfound = cs.n_exclusive_jets_ycut(_yc_frac[j]);
    _njet[j] += double(Nfound);
    switch (Nfound) {	 
    case 0: break; 
    case 1:  _frac1[j]++; break;
    case 2:  _frac2[j]++; break; 
    case 3:  _frac3[j]++; break; 
    case 4:  _frac4[j]++; break; 
    case 5:  _frac5[j]++; break; 
    default: _frac6[j]++;
    }
  }
  // DnD-rates 
  for(int j = 0; j<ddentr; j++) {
    int Nfound = cs.n_exclusive_jets_ycut(_d2dbins[j]);
    if (Nfound == 2) _d2dN2[j]++;
  }
  for(int j = 0; j<ddentr; j++) {
    int Nfound = cs.n_exclusive_jets_ycut(_d3dbins[j]);
    if (Nfound == 2) _d3dN2[j]++;
    else if (Nfound == 3) _d3dN3[j]++;
  }
  for(int j = 0; j<ddentr; j++) {
    int Nfound = cs.n_exclusive_jets_ycut(_d4dbins[j]);
    if (Nfound == 2) _d4dN2[j]++;
    else if (Nfound == 3) _d4dN3[j]++;
    else if (Nfound == 4) _d4dN4[j]++;
  }
}

NoPIOClassDescription<LEPJetAnalysis> LEPJetAnalysis::initLEPJetAnalysis;
// Definition of the static class description member.

void LEPJetAnalysis::Init() {

  static ClassDocumentation<LEPJetAnalysis> documentation
    ("LEP Jet data analysis",
     "The LEP Jet analysis uses data from \\cite{Pfeifenschneider:1999rz,Abreu:1996na}. ",
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

}

void LEPJetAnalysis::dofinish() {
  useMe();
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") 
    + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;
  double chisq,minfrac=0.05;
  unsigned int npoint;
  _y23->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL y_23 distribution or " 
		     << chisq/npoint << " per degree of freedom \n";
  _y23->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "y0231 compared to OPAL data",
		      " X  X                      ",
		      "1/NdN/dy0231",
		      "        X  X",
		       "y0231",
		       " X  X");
  _y34->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL y_34 distribution or " 
		     << chisq/npoint << " per degree of freedom \n";
  _y34->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "y0341 compared to OPAL data",
		      " X  X                      ",
		      "1/NdN/dy0341",
		      "        X  X",
		       "y0341",
		       " X  X");
  _y45->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL y_45 distribution or " 
		     << chisq/npoint << " per degree of freedom \n";
  _y45->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "y0451 compared to OPAL data",
		      " X  X                      ",
		      "1/NdN/dy0451",
		      "        X  X",
		       "y0451",
		       " X  X");
  _y56->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL y_56 distribution or " 
		     << chisq/npoint << " per degree of freedom \n";
  _y56->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "y0561 compared to OPAL data",
		      " X  X                      ",
		      "1/NdN/dy0561",
		      "        X  X",
		       "y0561",
		       " X  X");
  // data for jet fractions
  double jet2data[]={0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0    ,
		     0.0   ,  0.0   ,  0.0   ,  0.02  ,  0.06   ,
		     0.17  ,  0.43  ,  1.07  ,  2.35  ,  4.69   ,
		     8.36  , 13.46  , 19.77  , 26.96  , 34.42   ,
		     41.62  , 48.29  , 54.25  , 59.88  , 65.05  ,
		     70.02  , 74.71  , 79.04  , 83.01  , 86.68  ,
		     89.89  , 92.73  , 95.08  , 97.02  , 98.51  , 99.52  , 99.99};
  double jet2stat[]={0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,
		     0.0   ,0.0   ,0.0   ,0.00  ,0.00  ,
		     0.01  ,0.01  ,0.02  ,0.03  ,0.04  ,
		     0.05  ,0.07  ,0.08  ,0.08  ,0.09  ,
		     0.09  ,0.09  ,0.09  ,0.09  ,0.08  ,
		     0.08  ,0.08  ,0.07  ,0.07  ,0.06  ,
		     0.05  ,0.05  ,0.04  ,0.03  ,0.02  ,0.01  ,0.00};
  double jet2syst[]={0.0, 0.0, 0.0, 0.0, 0.0,
		     0.0, 0.0, 0.0,0.01,0.03,
		     0.06,0.13,0.27,0.51,0.88,
		     1.41,1.95,2.39,2.66,2.71,
		     2.63,2.46,2.26,2.08,1.87,
		     1.67,1.49,1.31,1.10,0.88,
		     0.69,0.50,0.31,0.17,0.09,0.04,0.04};
  double jet3data[]={0.0   , 0.0   , 0.0   , 0.0   , 0.0   ,
		     0.01  , 0.03  , 0.09  , 0.24  , 0.67  , 
		     1.63  , 3.57  , 6.82  ,11.57  ,17.39  ,
		     23.77  ,29.92  ,34.68  ,37.87  ,39.20  ,
		     38.91  ,37.48  ,35.47  ,32.83  ,29.92  ,
		     26.68  ,23.25  ,19.79  ,16.38  ,13.06  ,
		     10.03  , 7.26  , 4.92  , 2.98  , 1.48  , 0.47  , 0.0 };
  double jet3stat[]={0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,
		     0.00  ,0.00  ,0.01  ,0.01  ,0.01  ,
		     0.02  ,0.03  ,0.05  ,0.06  ,0.07  ,
		     0.08  ,0.08  ,0.08  ,0.08  ,0.08  ,
		     0.08  ,0.08  ,0.08  ,0.07  ,0.07  ,
		     0.07  ,0.07  ,0.06  ,0.06  ,0.05  ,
		     0.05  ,0.04  ,0.03  ,0.03  ,0.02  ,0.01  ,0.0};
  double jet3syst[]={0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,
		     0.01,0.01,0.03,0.08,0.21,
		     0.44,0.78,1.23,1.68,1.97,
		     2.02,1.73,1.14,0.64,0.59,
		     0.88,1.12,1.25,1.35,1.37,
		     1.32,1.24,1.13,0.98,0.81,
		     0.63,0.46,0.27,0.14,0.06,0.03,0.0};
  double jet4data[]={ 0.0   ,  0.0 ,  0.02  ,  0.03  ,  0.04  ,
		      0.09  ,  0.23  ,  0.60  ,  1.55  ,  3.45  ,
		      6.89  , 11.80  , 17.67  , 23.30  , 27.89  ,
		      30.19  , 29.93  , 27.81  , 24.05  , 19.72  ,
		      15.64  , 12.09  ,  9.13  ,  6.74  ,  4.79  ,
		      3.24  ,  2.06  ,  1.22  ,  0.66  ,  0.30  ,
		      0.12  ,  0.03  ,  0.0   ,  0.0  ,  0.0   ,  0.0  ,  0.0};
  double jet4stat[]={0.0   ,0.0   ,0.00  ,0.00  ,0.00  ,
		     0.01  ,0.01  ,0.01  ,0.02  ,0.03  ,
		     0.05  ,0.06  ,0.07  ,0.08  ,0.08  ,
		     0.08  ,0.08  ,0.07  ,0.07  ,0.06  ,
		     0.06  ,0.05  ,0.05  ,0.04  ,0.04  ,
		     0.03  ,0.02  ,0.02  ,0.01  ,0.01  ,
		     0.01  ,0.00  ,0.0   ,0.0   ,0.0   ,0.0   ,0.0};
  double jet4syst[]={0.0 ,0.0 ,0.02,0.01,0.01,
		     0.02,0.05,0.14,0.34,0.72,
		     1.30,1.87,2.09,1.89,1.35,
		     0.90,0.91,1.15,1.36,1.39,
		     1.23,1.01,0.83,0.63,0.46,
		     0.32,0.23,0.16,0.10,0.06,
		     0.03,0.02,0.0 ,0.0 ,0.0 ,0.0 ,0.0};
  double jet5data[]={0.02  ,  0.03  ,  0.04  ,  0.07  ,  0.17  ,
		     0.44  ,  1.14  ,  2.68  ,  5.52  ,  9.82  ,
		     15.10  , 20.28  , 24.06  , 25.57  , 24.61 ,
		     21.66  , 17.47  , 12.91  ,  8.83  ,  5.67 ,
		     3.43  ,  2.02  ,  1.14  ,  0.60  ,  0.30  ,
		     0.13  ,  0.05  ,  0.01  ,  0.0   ,  0.0   ,
		     0.0   ,  0.0   ,  0.0  ,  0.0    ,  0.0   ,  0.0  ,  0.0};
  double jet5stat[]={0.00  ,0.00  ,0.00  ,0.00  ,0.01  ,
		     0.01  ,0.02  ,0.03  ,0.04  ,0.06  ,
		     0.07  ,0.07  ,0.08  ,0.07  ,0.07  ,
		     0.07  ,0.06  ,0.05  ,0.04  ,0.04  ,
		     0.03  ,0.02  ,0.02  ,0.01  ,0.01  ,
		     0.01  ,0.00  ,0.00  ,0.0   ,0.0   ,
		     0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,0.0   ,0.0};
  double jet5syst[]={0.01,0.01,0.01,0.01,0.02,
		     0.06,0.20,0.52,1.04,1.61,
		     1.92,1.78,1.22,0.78,0.85,
		     1.13,1.30,1.26,1.04,0.73,
		     0.46,0.27,0.15,0.08,0.04,
		     0.02,0.01,0.00,0.0 ,0.0 ,
		     0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0};
  double jet6data[]={99.98  , 99.97  , 99.95  , 99.90  , 99.78  ,
		     99.45  , 98.58  , 96.63  , 92.76  , 86.25  ,
		     76.64  , 64.51  , 51.02  , 37.79  , 25.95  ,
		     16.48  ,  9.62  ,  5.15  ,  2.53  ,  1.14  ,
		     0.48  ,  0.20  ,  0.07  ,  0.03  ,  0.01   ,
		     0.0   ,  0.0  ,  0.0    ,  0.0   ,  0.0    ,
		     0.0   ,  0.0  ,  0.0    ,  0.0   ,  0.0    ,  0.0  ,  0.0};
  double jet6stat[]={0.00 ,0.00 ,0.00 ,0.01 ,0.01 ,
		     0.01 ,0.02 ,0.03 ,0.04 ,0.05 ,
		     0.06 ,0.07 ,0.07 ,0.07 ,0.07 ,
		     0.06 ,0.04 ,0.03 ,0.02 ,0.02 ,
		     0.01 ,0.01 ,0.00 ,0.00 ,0.00 ,
		     0.0  ,0.0  ,0.0  ,0.0  ,0.0  ,
		     0.0  ,0.0  ,0.0  ,0.0  ,0.0  ,0.0  ,0.0};
  double jet6syst[]={ 0.1, 0.1, 0.1, 0.2, 0.3,
		      0.4, 0.6, 1.0, 1.6, 2.4,
		      3.3, 4.0, 4.2, 3.8, 3.2,
		      2.4, 1.7, 1.0, 0.6, 0.3,
		      0.1, 0.0, 0.0, 0.0, 0.0,
		      0.0, 0.0, 0.0, 0.0, 0.0,
		      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double jet2error[37],jet3error[37],jet4error[37],jet5error[37],jet6error[38];
  for(unsigned int ix=0;ix<37;++ix) {
    jet2data[ix] /=100.;
    jet3data[ix] /=100.;
    jet4data[ix] /=100.;
    jet5data[ix] /=100.;
    jet6data[ix] /=100.;
    jet2error[ix] = sqrt(sqr(jet2stat[ix])+sqr(jet2syst[ix]))/100.;
    jet3error[ix] = sqrt(sqr(jet3stat[ix])+sqr(jet3syst[ix]))/100.;
    jet4error[ix] = sqrt(sqr(jet4stat[ix])+sqr(jet4syst[ix]))/100.;
    jet5error[ix] = sqrt(sqr(jet5stat[ix])+sqr(jet5syst[ix]))/100.;
    jet6error[ix] = sqrt(sqr(jet6stat[ix])+sqr(jet6syst[ix]))/100.;
  }
  for(unsigned int ix=2;ix<7;++ix) {
    vector<double> data,error;
    vector<int> obs;
    switch(ix) {
    case 2: 
      data  = vector<double>(jet2data ,jet2data+37 );
      error = vector<double>(jet2error,jet2error+37);
      obs   = _frac2;
      break;
    case 3: 
      data  = vector<double>(jet3data ,jet3data+37 );
      error = vector<double>(jet3error,jet3error+37);
      obs   = _frac3;
      break;
    case 4: 
      data  = vector<double>(jet4data ,jet4data+37 );
      error = vector<double>(jet4error,jet4error+37);
      obs   = _frac4;
      break;
    case 5: 
      data  = vector<double>(jet5data ,jet5data+37 );
      error = vector<double>(jet5error,jet5error+37);
      obs   = _frac5;
      break;
    case 6: 
      data  = vector<double>(jet6data ,jet6data+37 );
      error = vector<double>(jet6error,jet6error+37);
      obs   = _frac6;
      break;
    }
    output << "NEW FRAME\n";
    output << "SET WINDOW X 1.6 8 Y 3.5 9\n"; 
    output << "SET FONT DUPLEX\n";
    output << "TITLE TOP \"R0" << ix << "1 vs y0cut1\"\n";
    output << "CASE      \" X X     X   X\"\n";
    output << "TITLE LEFT \"R0" << ix << "1\"\n";
    output << "CASE       \" X X\"\n";
    if (HerwigStrategy::version != "") {
      output << "TITLE RIGHT \"" << HerwigStrategy::version << "\"\n";
      output << "CASE        \"\"\n";
    }
    output << "SET AXIS BOTTOM OFF\n";
    output << "SET ORDER X Y DY\n";
    output << "SET LIMITS X " << _yc_frac[0] << " " << _yc_frac[36] << "\n";
    output << "SET SCALE X LOG\n";
    for(unsigned int iy=0;iy<37;++iy) {
      output << _yc_frac[iy] << "\t" << double(obs[iy])/double(_nevent) << "\n";
    }
    output << "JOIN RED\n";
    for(unsigned int iy=0;iy<37;++iy) {
      output << _yc_frac[iy] << "\t" << data[iy] << "\t" << error[iy] << "\n";
    }
    output << "PLOT " << endl;
    output << "SET WINDOW X 1.6 8 Y 2.5 3.5\n";
    output << "SET LIMITS X " << _yc_frac[0] << " " << _yc_frac[36] << "\n";
    output << "SET SCALE X LOG\n";
    double ymax=0.;
    for(unsigned int iy=0;iy<37;++iy) {
      double y = data[iy]>0. ? error[iy]/data[iy] : 1.;
      if(y>ymax) ymax=y;
      output << _yc_frac[iy] << "\t" << y << "\n";
    }
    for(int iy=36;iy>=0;--iy) {
      double y = data[iy]>0. ? error[iy]/data[iy] : 1.;
      output << _yc_frac[iy] << "\t" << -y << "\n";
    }
    output << "set limits y " << -ymax << " " << ymax << "\n";
    output << "set fill full\n";
    output << "join yellow fill yellow\n";
    for(unsigned int iy=0;iy<37;++iy) {
      output << _yc_frac[iy] << "\t" 
	     << (double(obs[iy])/double(_nevent)-data[iy])/data[iy] << "\n";
    }
    output << "JOIN\n";
    output << "SET WINDOW X 1.6 8 Y 1.6 2.5\n";
    output << "SET LIMITS X " << _yc_frac[0] << " " << _yc_frac[36] << "\n";
    output << "SET SCALE X LOG\n";
    output << "SET AXIS BOTTOM ON\n";
    output << "TITLE BOTTOM \"y0" << ix << "1\"\n";
    output << "CASE         \" X X\n";
    ymax =0.;
    double ymin=0.,chisq=0.;
    int npoint=0;
    for(unsigned int iy=0;iy<37;++iy) {
      double point = data[iy]>0.&&error[iy]>0. ? 
	(double(obs[iy])/double(_nevent)-data[iy])/
	sqrt(sqr(error[iy])+double(obs[iy])/sqr(double(_nevent))) : 0.;
      if(point!=0.) ++npoint;
      if(point<ymin) ymin=point;
      if(point>ymax) ymax=point;
      output << _yc_frac[iy] << "\t" << point << "\n";
      if(point!=0.) {
	if(error[iy]>0.05*data[ix]) chisq+=sqr(point);
	else chisq+=sqr(double(obs[iy])/double(_nevent)-data[iy])/
	  (sqr(0.05*data[ix])+double(obs[iy])/sqr(double(_nevent)));
      }
    }
    output << "set limits y " << ymin << " " << ymax << "\n";
    output << "JOIN" << endl;
    generator()->log() << "Chi Square = " << chisq << " for " << npoint
		       << " degrees of freedom for OPAL R_" << ix 
		       << " distribution\n";
  }
  // n jet distributions
  double njetdata[]={15.662  , 14.500  , 13.377  , 12.296  , 11.264  ,
		     10.286  ,  9.363  ,  8.501  ,  7.701  ,  6.964  ,
		     6.294  ,  5.688  ,  5.142  ,  4.652  ,  4.219  ,
		     3.837  ,  3.510  ,  3.233  ,  3.008  ,  2.829  ,
		     2.688  ,  2.576  ,  2.483  ,  2.405  ,  2.336  ,
		     2.275  ,  2.222  ,  2.177  ,  2.136  ,  2.102  ,
		     2.073  ,  2.049  ,  2.030  ,  2.015  ,  2.005  ,  2.000};
  double njetstat[]={0.006 ,0.006 ,0.006 ,0.006 ,0.005 ,
		     0.005 ,0.005 ,0.005 ,0.004 ,0.004 ,
		     0.004 ,0.004 ,0.003 ,0.003 ,0.003 ,
		     0.003 ,0.003 ,0.003 ,0.003 ,0.003 ,
		     0.002 ,0.002 ,0.002 ,0.002 ,0.002 ,
		     0.002 ,0.002 ,0.002 ,0.002 ,0.002 ,
		     0.002 ,0.002 ,0.002 ,0.002 ,0.002 ,0.002 };
  double njetsyst[]={0.46, 0.42, 0.38, 0.35, 0.33,
		     0.32, 0.31, 0.30, 0.28, 0.26, 
		     0.23, 0.20, 0.17, 0.14, 0.12, 
		     0.10, 0.09, 0.07, 0.06, 0.04, 
		     0.03, 0.03, 0.02, 0.02, 0.02, 
		     0.01, 0.01, 0.01, 0.01, 0.00, 
		     0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
  double njeterror[36];
  for(unsigned int ix=0;ix<36;++ix)
    njeterror[ix]=sqrt(sqr(njetstat[ix])+sqr(njetsyst[ix]));
  output << "NEW FRAME\n";
  output << "SET WINDOW X 1.6 8 Y 3.5 9\n"; 
  output << "SET FONT DUPLEX\n";
  output << "TITLE TOP \" N0jets1\"\n";
  output << "CASE      \"  X    X\"\n";
  output << "TITLE LEFT \" <N0jets1>\"\n";
  output << "CASE       \"   X    X \"\n";
  if (HerwigStrategy::version != "") {
    output << "TITLE RIGHT \"" << HerwigStrategy::version << "\"\n";
    output << "CASE        \"\"\n";
  }
  output << "SET AXIS BOTTOM OFF\n";
  output << "SET ORDER X Y DY\n";
  output << "SET LIMITS X " << _yc_frac[1] << " " << _yc_frac[36] << "\n";
  output << "SET SCALE X LOG Y LOG\n";
  output << "SET LIMITS Y 2 16\n";
  for(unsigned int ix=1;ix<37;++ix) {
    output << _yc_frac[ix] << "\t" << _njet[ix].mean() << "\n";
  }
  output << "JOIN RED\n";
  for(unsigned int ix=0;ix<36;++ix) {
    output << _yc_frac[ix+1] << "\t" << njetdata[ix] << "\t" << njeterror[ix] << "\n";
  }
  output << "PLOT " << endl;
  output << "SET WINDOW X 1.6 8 Y 2.5 3.5\n";
  output << "SET LIMITS X " << _yc_frac[1] << " " << _yc_frac[36] << "\n";
  output << "SET SCALE X LOG Y LIN\n";
  double ymax=0.;
  for(unsigned int ix=0;ix<36;++ix) {
    double y = njetdata[ix]>0. ? njeterror[ix]/njetdata[ix] : 1.;
    if(y>ymax) ymax=y;
    output << _yc_frac[ix] << "\t" << y << "\n";
  }
  for(int ix=35;ix>=0;--ix) {
    double y = njetdata[ix]>0. ? njeterror[ix]/njetdata[ix] : 1.;
    output << _yc_frac[ix] << "\t" << -y << "\n";
  }
  output << "set limits y " << -ymax << " " << ymax << "\n";
  output << "set fill full\n";
  output << "join yellow fill yellow\n";
  for(unsigned int ix=1;ix<37;++ix) {
    output << _yc_frac[ix] << "\t" 
	   << (_njet[ix].mean()-njetdata[ix-1])/njetdata[ix-1] << "\n";
  }
  output << "JOIN\n";
  output << "SET WINDOW X 1.6 8 Y 1.6 2.5\n";
  output << "SET SCALE X LOG\n";
  output << "SET LIMITS X " << _yc_frac[1] << " " << _yc_frac[36] << "\n";
  output << "SET AXIS BOTTOM ON\n";
  output << "TITLE BOTTOM \"y0cut1\"\n";
  output << "CASE         \" X   X\"\n";
  ymax =0.;
  double ymin=0.;
  chisq=0.;
  npoint=0;
  for(unsigned int ix=1;ix<37;++ix) {
    double point = njetdata[ix-1]>0.&&njeterror[ix-1]>0. ? 
      (_njet[ix].mean()-njetdata[ix-1])/
      sqrt(sqr(njeterror[ix-1])+_njet[ix].mean_var()) : 0.;
    if(point!=0.) ++npoint;
    if(point<ymin) ymin=point;
    if(point>ymax) ymax=point;
    output << _yc_frac[ix] << "\t" << point << "\n";
    if(point!=0.) {
      if(njeterror[ix-1]>0.05*njetdata[ix-1]) chisq+=sqr(point);
      else chisq+=sqr(_njet[ix].mean()-njetdata[ix-1])/
	(sqr(0.05*njetdata[ix-1])+_njet[ix].mean_var());
    }
  }
  output << "set limits y " << ymin << " " << ymax << "\n";
  output << "JOIN" << endl;
  generator()->log() << "Chi Square = " << chisq << " for " << npoint
		     << " degrees of freedom for OPAL N_jet "
		     << " distribution\n";
  const int ddentr = 16;
  // D_2
  vector<double> d2,d2error;
  for(int j = 0; j<ddentr-1; j++) {
    double dy = _d2dbins[j+1]-_d2dbins[j];
    d2.push_back((_d2dN2[j+1]-_d2dN2[j])/(dy*double(_nevent)));
    d2error.push_back(sqrt(_d2dN2[j+1]+_d2dN2[j])/(dy*double(_nevent)));
  }
  // D_3
  vector<double> d3,d3error;
  for(int j = 0; j<ddentr-1; j++) {
    double dy = _d3dbins[j+1]-_d3dbins[j];
    d3.push_back((_d3dN2[j+1]-_d3dN2[j]+_d3dN3[j+1]-_d3dN3[j])/(dy*double(_nevent)));
    d3error.push_back(sqrt(_d3dN2[j+1]+_d3dN2[j]+_d3dN3[j+1]+_d3dN3[j])/
		      (dy*double(_nevent)));
  }
  // D_4
  vector<double> d4,d4error;
  for(int j = 0; j<ddentr-1; j++) {
    double dy = _d4dbins[j+1]-_d4dbins[j];
    d4.push_back((_d4dN2[j+1]-_d4dN2[j]+_d4dN3[j+1]-_d4dN3[j]
		  +_d4dN4[j+1]-_d4dN4[j])/(dy*double(_nevent)));
    d4error.push_back(sqrt(_d4dN2[j+1]+_d4dN2[j]+_d4dN3[j+1]+_d4dN3[j]
			   +_d4dN4[j+1]+_d4dN4[j])/(dy*double(_nevent)));
  }
  double d2data[]={63.81    ,11.589   , 6.040   , 3.823   , 2.697   ,
		   2.041   , 1.410   , 0.952   , 0.6708  , 0.4831  ,
		   0.3660  , 0.2830  , 0.2190  , 0.1357  , 0.0605};
  double d2stat[]={0.12    ,0.052   ,0.038   ,0.030   ,0.026   ,
		   0.022   ,0.013   ,0.011   ,0.0090  ,0.0077  ,
		   0.0068  ,0.0061  ,0.0054  ,0.0028  ,0.0022};
  double d2syst[]={0.64  ,0.116 ,0.060 ,0.038 ,0.027 ,
		   0.020 ,0.014 ,0.011 ,0.0102,0.0098,
		   0.0098,0.0102,0.0109,0.0092,0.0055};
  double d2total[15];
  double d3data[]={333.97    , 85.51    , 29.49    , 14.78    ,  7.583   ,
		   3.707   ,  2.119   ,  1.348   ,  0.906   ,  0.646   ,
		   0.366   ,  0.1965  ,  0.0877  ,  0.0367  ,  0.0073};
  double d3stat[]={0.64      ,0.31      ,0.18      ,0.13      ,0.067     ,
		   0.047     ,0.037     ,0.030     ,0.024     ,0.021     ,
		   0.011     ,0.0076  ,0.0050  ,0.0024  ,0.0010};
  double d3syst[]={3.34  ,0.86  ,0.29  ,0.15  ,0.076 ,
		   0.037 ,0.042 ,0.042 ,0.038 ,0.034 ,
		   0.025 ,0.0177,0.0100,0.0054,0.0014};
  double d3total[15];
  double d4data[]={711.0     ,184.32    , 53.35    , 21.18    , 10.15    ,
		   5.67    ,  2.768   ,  1.203   ,  0.617   ,  0.357   ,
		   0.215   ,  0.144   ,  0.0745  ,  0.0434  ,  0.0133};
  double d4stat[]={1.3     ,0.66    ,0.35    ,0.22    ,0.15    ,
		   0.11    ,0.057   ,0.039   ,0.028   ,0.022   ,
		   0.017   ,0.015   ,0.0103  ,0.0052  ,0.0028};
  double d4syst[]={7.1	,1.84	,1.06	,0.53	,0.30	,
		   0.20	,0.114	,0.060	,0.036	,0.024	,
		   0.016,0.011	,0.0064,0.0042,0.0014};
  double d4total[15];
  for(unsigned int ix=0;ix<15;++ix) {
    d2total[ix]=sqrt(sqr(d2stat[ix])+sqr(d2syst[ix]));
    d3total[ix]=sqrt(sqr(d3stat[ix])+sqr(d3syst[ix]));
    d4total[ix]=sqrt(sqr(d4stat[ix])+sqr(d4syst[ix]));
  }
  for(unsigned int ix=2;ix<5;++ix) {
    vector<double> data,error,obs,obserr,bins;
    switch(ix) {
    case 2: 
      data   = vector<double>(d2data ,d2data+15 );
      error  = vector<double>(d2total,d2total+15);
      obs    = d2;
      obserr = d2error;
      bins   = _d2dbins;
      break;
    case 3: 
      data   = vector<double>(d3data ,d3data+15 );
      error  = vector<double>(d3total,d3total+15);
      obs    = d3;
      obserr = d3error;
      bins   = _d3dbins;
      break;
    case 4: 
      data   = vector<double>(d4data ,d4data+15 );
      error  = vector<double>(d4total,d4total+15);
      obs    = d4;
      obserr = d4error;
      bins   = _d4dbins;
      break;
    }
    output << "NEW FRAME\n";
    output << "SET WINDOW X 1.6 8 Y 3.5 9\n"; 
    output << "SET FONT DUPLEX\n";
    output << "TITLE TOP \"D0" << ix << "1 vs y0cut1\"\n";
    output << "CASE      \" X X     X   X\"\n";
    output << "TITLE LEFT \"D0" << ix << "1\"\n";
    output << "CASE       \" X X\"\n";
    if (HerwigStrategy::version != "") {
      output << "TITLE RIGHT \"" << HerwigStrategy::version << "\"\n";
      output << "CASE        \"\"\n";
    }
    output << "SET AXIS BOTTOM OFF\n";
    output << "SET ORDER X Y DY\n";
    output << "SET LIMITS X " << bins[0] << " " << bins[ddentr-1] << "\n";
    output << "SET SCALE Y LOG\n";
    for(int iy=0;iy<ddentr-1;++iy) {
      output << 0.5*(bins[iy]+bins[iy+1]) << "\t" << obs[iy] << "\n";
    }
    output << "JOIN RED\n";
    for(int iy=0;iy<ddentr-1;++iy) {
      output << 0.5*(bins[iy]+bins[iy+1]) << "\t" 
	     << data[iy] << "\t" << error[iy] << "\n";
    }
    output << "PLOT " << endl;
    output << "SET WINDOW X 1.6 8 Y 2.5 3.5\n";
    output << "SET LIMITS X " << bins[0] << " " << bins[15] << "\n";
    output << "SET SCALE Y LIN\n";
    double ymax=0.;
    for(unsigned int iy=0;iy<ddentr-1;++iy) {
      double y = data[iy]>0. ? error[iy]/data[iy] : 1.;
      if(y>ymax) ymax=y;
      output << 0.5*(bins[iy]+bins[iy+1]) << "\t" << y << "\n";
    }
    for(int iy=ddentr-2;iy>=0;--iy) {
      double y = data[iy]>0. ? error[iy]/data[iy] : 1.;
      output << 0.5*(bins[iy]+bins[iy+1]) << "\t" << -y << "\n";
    }
    output << "set limits y " << -ymax << " " << ymax << "\n";
    output << "set fill full\n";
    output << "join yellow fill yellow\n";
    for(unsigned int iy=0;iy<ddentr-1;++iy) {
      output << 0.5*(bins[iy]+bins[iy+1]) << "\t" 
	     << (obs[iy]-data[iy])/data[iy] << "\n";
    }
    output << "JOIN\n";
    output << "SET WINDOW X 1.6 8 Y 1.6 2.5\n";
    output << "SET LIMITS X " << bins[0] << " " << bins[15] << "\n";;
    output << "SET AXIS BOTTOM ON\n";
    output << "TITLE BOTTOM \"y0" << ix << "1\"\n";
    output << "CASE         \" X X\n";
    ymax =0.;
    double ymin=0.,chisq=0.;
    int npoint=0;
    for(unsigned int iy=0;iy<ddentr-1;++iy) {
      double point = data[iy]>0.&&error[iy]>0. ? 
	(double(obs[iy])-data[iy])/
	sqrt(sqr(error[iy])+sqr(obserr[iy])) : 0.;
      if(point!=0.) ++npoint;
      if(point<ymin) ymin=point;
      if(point>ymax) ymax=point;
      output << 0.5*(bins[iy]+bins[iy+1]) << "\t" << point << "\n";
      if(point!=0.) {
	if(error[iy]>0.05*data[ix]) chisq+=sqr(point);
	else chisq+=sqr(obs[iy]-data[iy])/
	  (sqr(0.05*data[ix])+sqr(obserr[iy]));
      }
    }
    output << "set limits y " << ymin << " " << ymax << "\n";
    output << "JOIN" << endl;
    generator()->log() << "Chi Square = " << chisq << " for " << npoint
		       << " degrees of freedom for DELPHI D_" << ix 
		       << " distribution\n";
  }
}

void LEPJetAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  _nevent=0;
  // y distributions
  double diffbins[] = {1.00E-05, 1.33E-05, 1.78E-05, 2.37E-05, 3.16E-05, 
		       4.22E-05, 5.62E-05, 7.50E-05, 1.00E-04, 1.33E-04, 
		       1.78E-04, 2.37E-04, 3.16E-04, 4.22E-04, 5.62E-04, 
		       7.50E-04, 1.00E-03, 1.33E-03, 1.78E-03, 2.37E-03, 
		       3.16E-03, 4.22E-03, 5.62E-03, 7.50E-03, 1.00E-02, 
		       1.33E-02, 1.78E-02, 2.37E-02, 3.16E-02, 4.22E-02, 
		       5.62E-02, 7.50E-02, 1.00E-01, 1.33E-01, 1.78E-01, 
		       2.37E-01, 3.16E-01};
  double y23data[]= {0.0    ,  0.0    ,  0.0    ,  0.0    ,  2.672  ,
		     2.984  ,  2.482  ,  5.199  , 11.11   , 24.82   , 
		     43.76   , 80.73   ,121.9    ,166.2    ,195.5    ,
		     203.8    ,189.7    ,162.3    ,126.6    , 91.87   , 
		     63.59   , 42.77   , 30.13   , 20.75   , 14.93   , 
		     10.56   ,  7.318  ,  5.045  ,  3.483  ,  2.289  ,  
		     1.521  ,  0.942  ,  0.583  ,  0.336  ,  0.170  ,  0.058};
  double y23stat[]={0.0    ,0.0    ,0.0    ,0.0    ,2.248  ,
		    1.821  ,0.662  ,0.867  ,1.05   ,1.45   ,
		    1.56   ,1.93   ,2.1    ,2.1    ,2.0    ,
		    1.7    ,1.4    ,1.1    ,0.8    ,0.56   ,
		    0.40   ,0.28   ,0.20   ,0.15   ,0.11   ,
		    0.08   ,0.057  ,0.041  ,0.030  ,0.020  ,
		    0.015  ,0.010  ,0.007  ,0.004  ,0.003  ,0.002};
  double y23syst[]={0.0  ,0.0  ,0.0  ,0.0  ,2.296,
		    2.024,1.502,2.085, 4.56, 6.57,
		    10.18,15.48,19.2 ,23.3 ,25.6 ,
		    21.5 ,14.3 , 8.6 , 5.2 , 4.43,
		    3.77, 2.50, 1.38, 0.92, 0.65, 
		    0.44,0.321,0.270,0.209,0.146,
		    0.104,0.077,0.046,0.022,0.009,0.011};
  vector<double> bins(diffbins,diffbins+37);
  vector<double> data(y23data,y23data+36);
  vector<double> error(36);
  for(unsigned int ix=0;ix<36;++ix) error[ix]=sqrt(sqr(y23stat[ix])+sqr(y23syst[ix]));
  _y23 = new_ptr(Histogram(bins,data,error));
  double y34data[]={0.0   ,  0.0   ,  0.0   ,  5.910 ,  7.934 ,
		    14.89  , 34.42  , 64.74  ,138.3   ,242.9   ,
		    369.6   ,492.9   ,571.7   ,580.0   ,535.1   ,
		    448.7   ,333.2   ,234.0   ,149.3   , 88.43  ,
		    50.04  , 28.38  , 15.95  ,  9.092 ,  5.181 ,
		    2.848 ,  1.479 ,  0.727 ,  0.338 ,  0.134 , 
		    0.045 ,  0.011 ,  0.0   ,  0.0   ,  0.0   ,  0.0};
  double y34stat[]={0.0    , 0.0    , 0.0    , 1.837  , 1.577  ,
		    1.77   , 2.40   , 2.79   , 3.8    , 4.4    ,
		    4.7    , 4.7    , 4.3    , 3.7    , 3.0    ,
		    2.3    , 1.7    , 1.2    , 0.8    , 0.52   ,
		    0.34   , 0.22   , 0.15   , 0.096  , 0.064  ,
		    0.042  , 0.026  , 0.016  , 0.010  , 0.005  ,
		    0.003  , 0.001  , 0.0    , 0.0    , 0.0    , 0.0};
  double y34syst[]={0.0   ,0.0   ,0.0   ,3.084 ,2.557 ,
		    4.84 , 8.66 ,20.29 ,38.3  ,50.2  ,
		    60.3  ,66.6  ,59.8  ,45.1  ,30.5  ,
		    18.1  ,11.6  ,10.0  , 9.0  ,7.43  ,
		    5.09  ,2.94  ,1.75  ,1.052 ,0.580 ,
		    0.313 ,0.166 ,0.098 ,0.053 ,0.026 ,
		    0.014 ,0.007 ,0.0   ,0.0   ,0.0   ,0.0}; 
  data=vector<double>(y34data,y34data+36);
  for(unsigned int ix=0;ix<36;++ix) error[ix]=sqrt(sqr(y34stat[ix])+sqr(y34syst[ix]));   
  _y34 = new_ptr(Histogram(bins,data,error));    
  double y45data[]={17.50   ,  18.42   ,  19.98   ,  28.13   ,  56.68   ,
		    113.1    , 233.0    , 446.7    , 711.5    ,1018.6    ,
		    1197.7    ,1236.0    ,1108.2    , 907.8    , 658.8    ,
		    439.2    , 269.8    , 150.1    ,  76.52   ,  36.67   ,
		    16.23   ,   7.154  ,   3.150  ,   1.267  ,   0.537  , 
		    0.179  ,   0.064  ,   0.014  ,   0.003  ,   0.0    ,
		    0.0    ,   0.0    ,   0.0    ,   0.0    ,   0.0    ,   0.0};
  double y45stat[]={8.19     ,5.57          ,3.46     ,3.06     ,3.75        ,
		    4.6      ,6.1      ,7.7      ,8.6      ,9.1      ,
		    8.4      ,7.2      ,5.7      ,4.3      ,3.1      ,
		    2.1      ,1.4      ,0.9      ,0.55     ,0.33     ,
		    0.19     ,0.112    ,0.067    ,0.037    ,0.023    ,
		    0.011    ,0.006    ,0.003    ,0.001    ,0.0    ,
		    0.0      ,0.0      ,0.0      ,0.0      ,0.0      ,0.0};
  double y45syst[]={10.58   ,10.58   , 7.51   , 8.89   ,10.54   ,
		    19.5   , 48.8   , 90.5   ,138.9   ,169.5   ,
		    156.2   ,106.5   , 52.8   , 28.5   , 25.0   ,
		    23.3   , 20.0   , 15.7   ,10.60   , 6.24   ,
		    3.15   ,1.43    ,0.63    ,0.28    ,0.11    ,
		    0.04    ,0.01    ,0.00    ,0.00    ,0.0     ,
		    0.0     ,0.0     ,0.0     ,0.0     ,0.0     ,0.0}; 
  data=vector<double>(y45data,y45data+36);
  for(unsigned int ix=0;ix<36;++ix) error[ix]=sqrt(sqr(y45stat[ix])+sqr(y45syst[ix]));   
  _y45 = new_ptr(Histogram(bins,data,error)); 
  double y56data[]={42.43   ,  44.80   ,  78.16   , 155.9    , 309.6    ,
		    619.9    ,1055.3    ,1586.5    ,2003.2    ,2204.6    ,
		    2072.6    ,1715.4    ,1255.9    , 841.2    , 503.6    ,
		    273.2    , 133.7    ,  58.54   ,  23.36   ,   8.336  , 
		    2.735  ,   0.874  ,   0.253  ,   0.084  ,   0.012  , 
		    0.003  ,   0.0         ,   0.0    ,   0.0    ,   0.0    ,
		    0.0    ,   0.0    ,   0.0    ,   0.0    ,   0.0    ,   0.0};
  double y56stat[]={7.35     , 5.01     , 5.56     , 7.0      , 8.9      ,
		    11.7      ,13.8      ,15.0      ,14.5      ,13.0      ,
		    10.5      , 8.0      , 5.7      , 4.0      , 2.6      ,
		    1.6      , 1.0      , 0.55     , 0.30     , 0.157    ,
		    0.080    , 0.042    , 0.020    , 0.012    , 0.003    ,
		    0.002    , 0.0    , 0.0    , 0.0    , 0.0    ,
		    0.0      , 0.0    , 0.0    , 0.0    , 0.0    , 0.0};
  double y56syst[]={14.5,12.8,16.0, 20., 44.,
		    103.,195.,280.,307.,256.,
		    154., 84., 61., 51., 44.,
		    33., 21.,12.1, 6.2,2.86,
		    1.24,0.49,0.19,0.07,0.02,
		    0.01,0.0 ,0.0 ,0.0 ,0.0 ,
		    0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 }; 
  data=vector<double>(y56data,y56data+36);
  for(unsigned int ix=0;ix<36;++ix) error[ix]=sqrt(sqr(y56stat[ix])+sqr(y56syst[ix]));   
  _y56 = new_ptr(Histogram(bins,data,error));  
  double yc_frac[] = {1.00E-05, 1.33E-05, 1.78E-05, 2.37E-05, 3.16E-05, 
		      4.22E-05, 5.62E-05, 7.50E-05, 1.00E-04, 1.33E-04, 
		      1.78E-04, 2.37E-04, 3.16E-04, 4.22E-04, 5.62E-04, 
		      7.50E-04, 1.00E-03, 1.33E-03, 1.78E-03, 2.37E-03, 
		      3.16E-03, 4.22E-03, 5.62E-03, 7.50E-03, 1.00E-02, 
		      1.33E-02, 1.78E-02, 2.37E-02, 3.16E-02, 4.22E-02, 
		      5.62E-02, 7.50E-02, 1.00E-01, 1.33E-01, 1.78E-01, 
		      2.37E-01, 3.16E-01};
  _yc_frac=vector<double>(yc_frac,yc_frac+37);
  const int entr_fr = 37;
  _frac1.resize(entr_fr,0);
  _frac2.resize(entr_fr,0);
  _frac3.resize(entr_fr,0);
  _frac4.resize(entr_fr,0);
  _frac5.resize(entr_fr,0);
  _frac6.resize(entr_fr,0);
  _njet.resize(entr_fr);
  _njet = vector<Statistic>(entr_fr);
  const int ddentr = 16;
  double d2dbins[] = {0.000, 0.010, 0.020, 0.030, 0.040, 
		      0.050, 0.060, 0.080, 0.100, 0.120, 
		      0.140, 0.160, 0.180, 0.200, 0.250, 
		      0.300};
  double d3dbins[] = {0.000, 0.002, 0.004, 0.006, 0.008, 
		      0.012, 0.016, 0.020, 0.024, 0.028, 
		      0.032, 0.040, 0.050, 0.060, 0.080, 
		      0.100};
  double d4dbins[] = {0.000, 0.001, 0.002, 0.003, 0.004, 
		      0.005, 0.006, 0.008, 0.010, 0.012, 
		      0.014, 0.016, 0.018, 0.020, 0.025, 
		      0.030};
  _d2dbins=vector<double>(d2dbins,d2dbins+ddentr);
  _d3dbins=vector<double>(d3dbins,d3dbins+ddentr);
  _d4dbins=vector<double>(d4dbins,d4dbins+ddentr);
  _d2dN2.resize(ddentr);
  _d3dN2.resize(ddentr);
  _d3dN3.resize(ddentr);
  _d4dN2.resize(ddentr);
  _d4dN3.resize(ddentr);
  _d4dN4.resize(ddentr);
}
