// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TauTo3MesonAnalysis class.
//

#include "TauTo3MesonAnalysis.h"
#include "Herwig/Utilities/Histogram.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Herwig;

void TauTo3MesonAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  for(unsigned int ix=0;ix<4;++ix) {
    _m3pippimpim   .push_back(new_ptr(Histogram(0.,1.8,200)));
    _m3pi0pi0pim   .push_back(new_ptr(Histogram(0.,1.8,200)));
    _m3kmpimkp     .push_back(new_ptr(Histogram(0.,1.8,200)));
    _m3k0pimk0     .push_back(new_ptr(Histogram(0.,1.8,200)));
    _m3kmpi0k0     .push_back(new_ptr(Histogram(0.,1.8,200)));
    _m3pi0pi0km    .push_back(new_ptr(Histogram(0.,1.8,200)));
    _m3kmpimpip    .push_back(new_ptr(Histogram(0.,1.8,200)));
    _m3pimk0pi0    .push_back(new_ptr(Histogram(0.,1.8,200)));
    _m3pimpi0eta   .push_back(new_ptr(Histogram(0.,1.8,200)));
    _m3pimpi0gamma .push_back(new_ptr(Histogram(0.,1.8,200)));
    _m3kspimks     .push_back(new_ptr(Histogram(0.,1.8,200)));
    _m3klpimkl     .push_back(new_ptr(Histogram(0.,1.8,200)));
    _m3kspimkl     .push_back(new_ptr(Histogram(0.,1.8,200)));
  }
}

void TauTo3MesonAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  using namespace HistogramOptions;

  _m3pippimpim[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P2+3P2-3P2-3 mass in T2-3RN0T1P2+3P2-3P2-3",
			    "GX XGX XGX X         GX XWGXGXGX XGX XGX X",
			    "1/SdS/dm0P2+3P2-3P2-31/GeV2-13",
			    "  G G   XGX XGX XGX XX    X  X",
			    "m0P2+3P2-3P2-31/GeV",
			    " XGX XGX XGX XX    ");
  _m3pippimpim[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P2-3P2-3 mass in T2-3RN0T1P2+3P2-3P2-3",
			    "GX XGX X         GX XWGXGXGX XGX XGX X",
			    "1/SdS/dm0P2-3P2-31/GeV2-13",
			    "  G G   XGX XGX XX    X  X",
			    "m0P2-3P2-31/GeV",
			    " XGX XGX XX    ");
  _m3pippimpim[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P2+3P2-3 mass in T2-3RN0T1P2+3P2-3P2-3",
			    "GX XGX X         GX XWGXGXGX XGX XGX X",
			    "1/SdS/dm0P2+3P2-31/GeV2-13",
			    "  G G   XGX XGX XX    X  X",
			    "m0P2+3P2-31/GeV",
			    " XGX XGX XX    ");
  _m3pi0pi0pim[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P2-3P203P203 mass in T2-3RN0T1P2-3P203P203",
			    "GX XGX XGX X         GX XWGXGXGX XGX XGX X",
			    "1/SdS/dm0P2-3P203P2031/GeV2-13",
			    "  G G   XGX XGX XGX XX    X  X",
			    "m0P2-3P203P2031/GeV",
			    " XGX XGX XGX XX    ");
  _m3pi0pi0pim[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P203P203 mass in T2-3RN0T1P2-3P203P203",
			    "GX XGX X         GX XWGXGXGX XGX XGX X",
			    "1/SdS/dm0P203P2031/GeV2-13",
			    "  G G   XGX XGX XX    X  X",
			    "m0P203P2031/GeV",
			    " XGX XGX XX    ");
  _m3pi0pi0pim[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P2-3P203 mass in T2-3RN0T1P2-3P203P203",
			    "GX XGX X         GX XWGXGXGX XGX XGX X",
			    "1/SdS/dm0P2-3P2031/GeV2-13",
			    "  G G   XGX XGX XX    X  X",
			    "m0P2-3P2031/GeV",
			    " XGX XGX XX    ");
  _m3kmpimkp[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "K2-3P2-3K2+3 mass in T2-3RN0T1K2-3P2-3K2+3",
			    " X XGX X X X         GX XWGXGX X XGX X X X",
			    "1/SdS/dm0K2-3P2-3K2+31/GeV2-13",
			    "  G G   X X XGX X X XX    X  X",
			    "m0K2-3P2-3K2+31/GeV",
			    " X X XGX X X XX    ");
  _m3kmpimkp[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "K2-3P2-3 mass in T2-3RN0T1K2-3P2-3K2+3",
			    " X XGX X         GX XWGXGX X XGX X X X",
			    "1/SdS/dm0K2-3P2-31/GeV2-13",
			    "  G G   X X XGX XX    X  X",
			    "m0K2-3P2-31/GeV",
			    " X X XGX XX    ");
  _m3kmpimkp[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "K2-3K2+3 mass in T2-3RN0T1K2-3P2-3K2+3",
			    " X X X X         GX XWGXGX X XGX X X X",
			    "1/SdS/dm0K2-3K2+31/GeV2-13",
			    "  G G   X X X X XX    X  X",
			    "m0K2-3K2+31/GeV",
			    " X X X X XX    ");
  _m3kmpimkp[3]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P2-3K2+3 mass in T2-3RN0T1K2-3P2-3K2+3",
			    "GX X X X         GX XWGXGX X XGX X X X",
			    "1/SdS/dm0P2-3K2+31/GeV2-13",
			    "  G G   XGX X X XX    X  X",
			    "m0P2-3K2+31/GeV",
			    " XGX X X XX    ");
  _m3k0pimk0[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K203P2-3K0O203 mass in T2-3RN0T1K203P2-3K0O203",
			       " X XGX X UDX X         GX XWGXGX X XGX X UDX X",
			       "1/SdS/dm0K203P2-3K0O2031/GeV2-13",
			       "  G G   X X XGX X UDX XX    X  X",
			       "m0K203P2-3K0O2031/GeV",
			       " X X XGX X UDX XX    ");
  _m3k0pimk0[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K203P2-3 mass in T2-3RN0T1K203P2-3K0O2033",
			       " X XGX X         GX XWGXGX X XGX X UDX XX",
			       "1/SdS/dm0K203P2-31/GeV2-13",
			       "  G G   X X XGX XX    X  X",
			       "m0K203P2-31/GeV",
			       " X X XGX XX    ");
  _m3k0pimk0[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K203K0O203 mass in T2-3RN0T1K203P2-3K0O203",
			       " X X UDX X         GX XWGXGX X XGX X UDX X",
			       "1/SdS/dm0K203K0O2031/GeV2-13",
			       "  G G   X X X UDX XX    X  X",
			       "m0K203K0O2031/GeV",
			       " X X X UDX XX    ");
  _m3k0pimk0[3]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "P2-3K0O203 mass in T2-3RN0T1K203P2-3K0O203",
			       "GX X UDX X         GX XWGXGX X XGX X UDX X",
			       "1/SdS/dm0P2-3K0O2031/GeV2-13",
			       "  G G   XGX X UDX XX    X  X",
			       "m0P2-3K0O2031/GeV",
			       " XGX X UDX XX    ");
  _m3kmpi0k0[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K2-3P203K203 mass in T2-3RN0T1K2-3P203K203",
			       " X XGX X X X         GX XWGXGX X XGX X X X",
			       "1/SdS/dm0K2-3P203K2031/GeV2-13",
			       "  G G   X X XGX X X XX    X  X",
			       "m0K2-3P203K2031/GeV",
			       " X X XGX X X XX    ");
  _m3kmpi0k0[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K2-3P203 mass in T2-3RN0T1K2-3P203K203",
			       " X XGX X         GX XWGXGX X XGX X X X",
			       "1/SdS/dm0K2-3P2031/GeV2-13",
			       "  G G   X X XGX XX    X  X",
			       "m0K2-3P2031/GeV",
			       " X X XGX XX    ");
  _m3kmpi0k0[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K2-3K203 mass in T2-3RN0T1K2-3P203K203",
			       " X X X X         GX XWGXGX X XGX X X X",
			       "1/SdS/dm0K2-3K2031/GeV2-13",
			       "  G G   X X X X XX    X  X",
			       "m0K2-3K2031/GeV",
			       " X X X X XX    ");
  _m3kmpi0k0[3]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "P203K203 mass in T2-3RN0T1K2-3P203K203",
			       "GX X X X         GX XWGXGX X XGX X X X",
			       "1/SdS/dm0P203K2031/GeV2-13",
			       "  G G   XGX X X XX    X  X",
			       "m0P203K2031/GeV",
			       " XGX X X XX    ");
  _m3pi0pi0km[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",
				"P203P203K2-3 mass in T2-3RN0T1P203P203K2-3",
				"GX XGX X X X         GX XWGXGXGX XGX X X X",
				"1/SdS/dm0P203P203K2-31/GeV2-13",
				"  G G   XGX XGX X X XX    X  X",
				"m0P203P203K2-31/GeV",
				" XGX XGX X X XX    ");
  _m3pi0pi0km[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",
				"P203P203 mass in T2-3RN0T1P203P203K2-3",
				"GX XGX X         GX XWGXGXGX XGX X X X",
				"1/SdS/dm0P203P2031/GeV2-13",
				"  G G   XGX XGX XX    X  X",
				"m0P203P2031/GeV",
				" XGX XGX XX    ");
  _m3pi0pi0km[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",
				"P203K2-3 mass in T2-3RN0T1P203P203K2-3",
				"GX X X X         GX XWGXGXGX XGX X X X",
				"1/SdS/dm0P203K2-31/GeV2-13",
				"  G G   XGX X X XX    X  X",
				"m0P203K2-31/GeV",
				" XGX X X XX    ");
  _m3kmpimpip[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",
				"K2-3P2-3P2+3 mass in T2-3RN0T1K2-3P2-3P2+3",
				" X XGX XGX X         GX XWGXGX X XGX XGX X",
				"1/SdS/dm0K2-3P2-3P2+31/GeV2-13",
				"  G G   X X XGX XGX XX    X  X",
				"m0K2-3P2-3P2+31/GeV",
				" X X XGX XGX XX    ");
  _m3kmpimpip[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",
				"K2-3P2-3 mass in T2-3RN0T1K2-3P2-3P2+3",
				" X XGX X         GX XWGXGX X XGX XGX X",
				"1/SdS/dm0K2-3P2-31/GeV2-13",
				"  G G   X X XGX XX    X  X",
				"m0K2-3P2-31/GeV",
				" X X XGX XX    ");
  _m3kmpimpip[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",
				"K2-3P2+3 mass in T2-3RN0T1K2-3P2-3P2+3",
				" X XGX X         GX XWGXGX X XGX XGX X",
				"1/SdS/dm0P2-3P2+31/GeV2-13",
				"  G G   XGX XGX XX    X  X",
				"m0K2-3P2+31/GeV",
				" X X XGX XX    ");
  _m3kmpimpip[3]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",
				"P2-3P2+3 mass in T2-3RN0T1K2-3P2-3P2+3",
				"GX XGX X         GX XWGXGX X XGX XGX X",
				"1/SdS/dm0P2-3P2+31/GeV2-13",
				"  G G   XGX XGX XX    X  X",
				"m0P2-3P2+31/GeV",
				" XGX XGX XX    ");
  _m3pimk0pi0[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",
				"P2-3P203K0O203 mass in T2-3RN0T1P2-3P203K0O203",
				"GX XGX X UDX X         GX XWGXGXGX XGX X UDX X",
				"1/SdS/dm0P2-3P203K0O2031/GeV2-13",
				"  G G   XGX XGX X UDX XX    X  X",
				"m0P2-3P203K0O2031/GeV",
				" XGX XGX X UDX XX    ");
  _m3pimk0pi0[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",
				"P2-3K0O203 mass in T2-3RN0T1P2-3P203K0O203",
				"GX X UDX X         GX XWGXGXGX XGX X UDX X",
				"1/SdS/dm0P2-3K0O2031/GeV2-13",
				"  G G   XGX X UDX XX    X  X",
				"m0P2-3K0O2031/GeV",
				" XGX X UDX XX    ");
  _m3pimk0pi0[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",
				"P2-3P203 mass in T2-3RN0T1P2-3P203K0O203",
				"GX XGX X         GX XWGXGXGX XGX X UDX X",
				"1/SdS/dm0P2-3P2031/GeV2-13",
				"  G G   XGX XGX XX    X  X",
				"m0P2-3P2031/GeV",
				" XGX XGX XX    ");
  _m3pimk0pi0[3]->topdrawOutput(output,Frame|Errorbars|Ylog,
				"RED",
				"P203K0O203 mass in T2-3RN0T1P2-3P203K0O203",
				"GX X UDX X         GX XWGXGXGX XGX X UDX X",
				"1/SdS/dm0P203K0O2031/GeV2-13",
				"  G G   XGX X UDX XX    X  X",
				"m0P203K0O2031/GeV",
				" XGX X UDX XX    ");
  _m3pimpi0eta[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P2-3P203H mass in T2-3RN0T1P2-3P203H",
			    "GX XGX XG         GX XWGXGXGX XGX XG",
			    "1/SdS/dm0P2-3P203H1/GeV2-13",
			    "  G G   XGX XGX XGX    X  X",
			    "m0P2-3P203H1/GeV",
			    " XGX XGX XGX    ");
  _m3pimpi0eta[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P2-3P203 mass in T2-3RN0T1P2-3P203H",
			    "GX XGX X         GX XWGXGXGX XGX XG",
			    "1/SdS/dm0P2-3P2031/GeV2-13",
			    "  G G   XGX XGX XX    X  X",
			    "m0P2-3P2031/GeV",
			    " XGX XGX XX    ");
  _m3pimpi0eta[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P2-3H mass in T2-3RN0T1P2-3P203H",
			    "GX XG         GX XWGXGXGX XGX XG",
			    "1/SdS/dm0P2-3H1/GeV2-13",
			    "  G G   XGX XGX    X  X",
			    "m0P2-3H1/GeV",
			    " XGX XGX    ");
  _m3pimpi0eta[3]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P203H mass in T2-3RN0T1P2-3P203H",
			    "GX XG         GX XWGXGXGX XGX XG",
			    "1/SdS/dm0P203H1/GeV2-13",
			    "  G G   XGX XGX    X  X",
			    "m0P203H1/GeV",
			    " XGX XGX    ");
  _m3pimpi0gamma[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P2-3P203G mass in T2-3RN0T1P2-3P203G",
			    "GX XGX XG         GX XWGXGXGX XGX XG",
			    "1/SdS/dm0P2-3P203G1/GeV2-13",
			    "  G G   XGX XGX XGX    X  X",
			    "m0P2-3P203G1/GeV",
			    " XGX XGX XGX    ");
  _m3pimpi0gamma[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P2-3P203 mass in T2-3RN0T1P2-3P203G",
			    "GX XGX X         GX XWGXGXGX XGX XG",
			    "1/SdS/dm0P2-3P2031/GeV2-13",
			    "  G G   XGX XGX XX    X  X",
			    "m0P2-3P2031/GeV",
			    " XGX XGX XX    ");
  _m3pimpi0gamma[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P2-3G mass in T2-3RN0T1P2-3P203G",
			    "GX XG         GX XWGXGXGX XGX XG",
			    "1/SdS/dm0P2-3G1/GeV2-13",
			    "  G G   XGX XGX    X  X",
			    "m0P2-3G1/GeV",
			    " XGX XGX    ");
  _m3pimpi0gamma[3]->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "P203G mass in T2-3RN0T1P2-3P203G",
			    "GX XG         GX XWGXGXGX XGX XG",
			    "1/SdS/dm0P203G1/GeV2-13",
			    "  G G   XGX XGX    X  X",
			    "m0P203G1/GeV",
			    " XGX XGX    ");
  _m3kspimks[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K2030S1P2-3K2030S1 mass in T2-3RN0T1K2030S1P2-3K2030S1",
			       " X XX XGX X X XX X         GX XWGXGX X XX XGX X X XX X",
			       "1/SdS/dm0K2030S1P2-3K2030S11/GeV2-13",
			       "  G G   X X XX XGX X X XX XX    X  X",
			       "m0K2030S1P2-3K2030S11/GeV",
			       " X X XX XGX X X XX XX    ");
  _m3kspimks[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K2030S1P2-3 mass in T2-3RN0T1K2030S1P2-3K2030S1",
			       " X XX XGX X         GX XWGXGX X XX XGX X X XX X",
			       "1/SdS/dm0K2030S1P2-31/GeV2-13",
			       "  G G   X X XX XGX XX    X  X",
			       "m0K2030S1P2-31/GeV",
			       " X X XX XGX XX    ");
  _m3kspimks[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K2030S1K2030S1 mass in T2-3RN0T1K2030S1P2-3K2030S1",
			       " X XX X X XX X         GX XWGXGX X XX XGX X X XX X",
			       "1/SdS/dm0K2030S1K2030S11/GeV2-13",
			       "  G G   X X XX X X XX XX    X  X",
			       "m0K2030S1K2030S11/GeV",
			       " X X XX X X XX XX    ");
  _m3klpimkl[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K2030L1P2-3K2030L1 mass in T2-3RN0T1K2030L1P2-3K2030L1",
			       " X XX XGX X X XX X         GX XWGXGX X XX XGX X X XX X",
			       "1/SdS/dm0K2030L1P2-3K2030L11/GeV2-13",
			       "  G G   X X XX XGX X X XX XX    X  X",
			       "m0K2030L1P2-3K2030L11/GeV",
			       " X X XX XGX X X XX XX    ");
  _m3klpimkl[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K2030L1P2-3 mass in T2-3RN0T1K2030L1P2-3K2030L1",
			       " X XX XGX X         GX XWGXGX X XX XGX X X XX X",
			       "1/SdS/dm0K2030L1P2-31/GeV2-13",
			       "  G G   X X XX XGX XX    X  X",
			       "m0K2030L1P2-31/GeV",
			       " X X XX XGX XX    ");
  _m3klpimkl[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K2030L1K2030L1 mass in T2-3RN0T1K2030L1P2-3K2030L1",
			       " X XX X X XX X         GX XWGXGX X XX XGX X X XX X",
			       "1/SdS/dm0K2030L1K2030L11/GeV2-13",
			       "  G G   X X XX X X XX XX    X  X",
			       "m0K2030L1K2030L11/GeV",
			       " X X XX X X XX XX    ");
  _m3kspimkl[0]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K2030S1P2-3K0L1 mass in T2-3RN0T1K2030S1P2-3K0L1",
			       " X XX XGX X X X         GX XWGXGX X XX XGX X X X",
			       "1/SdS/dm0K2030S1P2-3K0L11/GeV2-13",
			       "  G G   X X XX XGX X X XX    X  X",
			       "m0K2030S1P2-3K0L11/GeV",
			       " X X XX XGX X X XX    ");
  _m3kspimkl[1]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K2030S1P2-3 mass in T2-3RN0T1K2030S1P2-3K0L1",
			       " X XX XGX X         GX XWGXGX X XX XGX X X X",
			       "1/SdS/dm0K2030S1P2-31/GeV2-13",
			       "  G G   X X XX XGX XX    X  X",
			       "m0K2030S1P2-31/GeV",
			       " X X XX XGX XX    ");
  _m3kspimkl[2]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "K2030S1K0L1 mass in T2-3RN0T1K2030S1K0L1",
			       " X XX X X X         GX XWGXGX X XX X X X",
			       "1/SdS/dm0K2030S1K0L11/GeV2-13",
			       "  G G   X X XX X X XX    X  X",
			       "m0K2030S1K0L11/GeV",
			       " X X XX X X XX    ");
  _m3kspimkl[3]->topdrawOutput(output,Frame|Errorbars|Ylog,
			       "RED",
			       "1P2-3K0L1 mass in T2-3RN0T1K2030S1P2-3K0L1",
			       "XGX X X X         GX XWGXGX X XX XGX X X X",
			       "1/SdS/dm0P2-3K0L11/GeV2-13",
			       "  G G   XGX X X XX    X  X",
			       "m0P2-3K0L11/GeV",
			       " XGX X X XX    ");
}



void TauTo3MesonAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  tPVector hadrons=event->getFinalState();
  map<tPPtr,ParticleVector> taus;
  for(unsigned int ix=0;ix<hadrons.size();++ix) {
    PPtr mother=hadrons[ix];
    do {
      if(!mother->parents().empty()) mother=mother->parents()[0];
      else                           mother=tPPtr();
    }
    while(mother&&abs(mother->id())!=ParticleID::tauminus);
    if(mother&&abs(mother->id())==ParticleID::tauminus) {
      if(taus.find(mother)==taus.end()) {
	taus.insert(make_pair(mother,ParticleVector()));
      }
      taus[mother].push_back(hadrons[ix]);
    }
  }
  map<tPPtr,ParticleVector>::const_iterator tit;
  for(tit=taus.begin();tit!=taus.end();++tit) {
    if(tit->second.size()!=4) continue;
    ParticleVector decay=tit->second;
    int tsign=tit->first->id()/abs(tit->first->id());
    vector<Lorentz5Momentum> ppi0,ppip,ppim,pkp,pkm,pk0,pk0bar,peta,pgamma,pkl,pks;
    Lorentz5Momentum ptotal;
    for(unsigned int ix=0;ix<decay.size();++ix) {
       long id = decay[ix]->id()*tsign;
       if(abs(id)!=ParticleID::nu_tau) ptotal+=decay[ix]->momentum();
       if(id==ParticleID::piplus)            ppip  .push_back(decay[ix]->momentum());
       else if(id==ParticleID::piminus)      ppim  .push_back(decay[ix]->momentum());
       else if(abs(id)==ParticleID::pi0)     ppi0  .push_back(decay[ix]->momentum());
       else if(id==ParticleID::Kplus)        pkp   .push_back(decay[ix]->momentum());
       else if(id==ParticleID::Kminus)       pkm   .push_back(decay[ix]->momentum());
       else if(id==ParticleID::K0)           pk0   .push_back(decay[ix]->momentum());
       else if(id==ParticleID::Kbar0)        pk0bar.push_back(decay[ix]->momentum());
       else if(id==ParticleID::eta)          peta  .push_back(decay[ix]->momentum());
       else if(id==ParticleID::gamma)        pgamma.push_back(decay[ix]->momentum());
       else if(abs(id)==ParticleID::K_S0)    pks   .push_back(decay[ix]->momentum());
       else if(abs(id)==ParticleID::K_L0)    pkl   .push_back(decay[ix]->momentum());
    }
    if(ppim.size()==2&&ppip.size()==1) {
      *_m3pippimpim[0] += ptotal.m()/GeV;
      *_m3pippimpim[1] += (ppim[0]+ppim[1]).m()/GeV;
      *_m3pippimpim[2] += (ppim[0]+ppip[0]).m()/GeV;
      *_m3pippimpim[2] += (ppim[1]+ppip[0]).m()/GeV;
    }
    else if(ppim.size()==1&&ppi0.size()==2) {
      *_m3pi0pi0pim[0] += ptotal.m()/GeV;
      *_m3pi0pi0pim[1] += (ppi0[0]+ppi0[1]).m()/GeV;
      *_m3pi0pi0pim[2] += (ppi0[0]+ppim[0]).m()/GeV;
      *_m3pi0pi0pim[2] += (ppi0[1]+ppim[0]).m()/GeV;
    }
    else if(pkm.size()==1&&pkp.size()==1&&ppim.size()==1) {
      *_m3kmpimkp[0] += ptotal.m()/GeV;
      *_m3kmpimkp[1] += (pkm[0]+ppim[0]).m()/GeV;
      *_m3kmpimkp[2] += (pkm[0]+ pkp[0]).m()/GeV;
      *_m3kmpimkp[3] += (pkp[0]+ppim[0]).m()/GeV;
    }
    else if(pk0.size()==1&&pk0bar.size()==1&&ppim.size()==1) {
      *_m3k0pimk0[0] += ptotal.m()/GeV;
      *_m3k0pimk0[1] += (pk0[0]   +ppim[0]  ).m()/GeV;
      *_m3k0pimk0[2] += (pk0[0]   +pk0bar[0]).m()/GeV;
      *_m3k0pimk0[3] += (pk0bar[0]+ppim[0]  ).m()/GeV;
    }
    else if(pk0.size()==1&&pkm.size()==1&&ppi0.size()==1) {
      *_m3kmpi0k0[0] += ptotal.m()/GeV;
      *_m3kmpi0k0[1] += (pkm[0]+ppi0[0]).m()/GeV;
      *_m3kmpi0k0[2] += (pkm[0]+pk0[0] ).m()/GeV;
      *_m3kmpi0k0[3] += (pk0[0]+ppi0[0]).m()/GeV;
    }
    else if(ppi0.size()==2&&pkm.size()==1) {
      *_m3pi0pi0km[0] += ptotal.m()/GeV;
      *_m3pi0pi0km[1] += (ppi0[0]+ppi0[1]).m()/GeV;
      *_m3pi0pi0km[2] += (ppi0[0]+pkm[0] ).m()/GeV;
      *_m3pi0pi0km[3] += (ppi0[1]+pkm[0] ).m()/GeV;
    }
    else if(pkm.size()==1&&ppim.size()==1&&ppip.size()==1) {
      *_m3kmpimpip[0] += ptotal.m()/GeV;
      *_m3kmpimpip[1] += (pkm[0]+ppim[0]).m()/GeV;
      *_m3kmpimpip[2] += (pkm[0]+ppip[0] ).m()/GeV;
      *_m3kmpimpip[3] += (ppip[0]+ppim[0] ).m()/GeV;
    }
    else if(ppim.size()==1&&pk0bar.size()==1&&ppi0.size()==1) {
      *_m3pimk0pi0[0] += ptotal.m()/GeV;
      *_m3pimk0pi0[1] += (ppim[0]+pk0bar[0]).m()/GeV;
      *_m3pimk0pi0[2] += (ppim[0]+ppi0[0]  ).m()/GeV;
      *_m3pimk0pi0[3] += (pk0bar[0]+ppi0[0]).m()/GeV;
    }
    else if(ppim.size()==1&&ppi0.size()==1&&peta.size()==1) {
      *_m3pimpi0eta[0] += ptotal.m()/GeV;
      *_m3pimpi0eta[1] += (ppim[0]+ppi0[0]).m()/GeV;
      *_m3pimpi0eta[2] += (ppim[0]+peta[0]).m()/GeV;
      *_m3pimpi0eta[3] += (ppi0[0]+peta[0]).m()/GeV;
    }
    else if(ppim.size()==1&&ppi0.size()==1&&pgamma.size()==1) {
      *_m3pimpi0gamma[0] += ptotal.m()/GeV;
      *_m3pimpi0gamma[1] += (ppim[0]+ppi0[0]).m()/GeV;
      *_m3pimpi0gamma[2] += (ppim[0]+pgamma[0]).m()/GeV;
      *_m3pimpi0gamma[3] += (ppi0[0]+pgamma[0]).m()/GeV;
    }
    else if(pks.size()==2&&ppim.size()==1) {
      *_m3kspimks[0] += ptotal.m()/GeV;
      *_m3kspimks[1] += (ppim[0]+pks[0]).m()/GeV;
      *_m3kspimks[1] += (ppim[0]+pks[1]).m()/GeV;
      *_m3kspimks[2] += (pks [0]+pks[1]).m()/GeV;
    }
    else if(pkl.size()==2&&ppim.size()==1) {
      *_m3klpimkl[0] += ptotal.m()/GeV;
      *_m3klpimkl[1] += (ppim[0]+pkl[0]).m()/GeV;
      *_m3klpimkl[1] += (ppim[0]+pkl[1]).m()/GeV;
      *_m3klpimkl[2] += (pkl [0]+pkl[1]).m()/GeV;
    }
    else if(pks.size()==1&&pkl.size()==1&&ppim.size()==1) {
      *_m3kspimkl[0] += ptotal.m()/GeV;
      *_m3kspimkl[1] += (ppim[0]+pks[0]).m()/GeV;
      *_m3kspimkl[2] += (pks[0] +pkl[0]).m()/GeV;
      *_m3kspimkl[3] += (ppim[0]+pkl[0]).m()/GeV;
    }
  }
}

NoPIOClassDescription<TauTo3MesonAnalysis> TauTo3MesonAnalysis::initTauTo3MesonAnalysis;
// Definition of the static class description member.

void TauTo3MesonAnalysis::Init() {

  static ClassDocumentation<TauTo3MesonAnalysis> documentation
    ("The TauTo2MesonAnalysis class plots the mass distributions of tau decays to"
     "three mesons.");

}

