// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#include "GnuplotOutput.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GnuplotOutput.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Analysis2;

GnuplotOutput::~GnuplotOutput() {}

void GnuplotOutput::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _subproMult << _ratio << _chi2 << _gnuplot;
}

void GnuplotOutput::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _subproMult >> _ratio >> _chi2 >> _gnuplot;
}

ClassDescription<GnuplotOutput> GnuplotOutput::initGnuplotOutput;
// Definition of the static class description member.

void GnuplotOutput::Init() {

  static ClassDocumentation<GnuplotOutput> documentation
    ("Histogram2 output to gnuplot");


  static Parameter<GnuplotOutput,string> interfaceMultiplicityRange
    ("MultiplicityRange",
     "The range of hard jet multiplicities to plot seperatly",
     &GnuplotOutput::_subproMult, "",
     false, false);


  static Switch<GnuplotOutput,bool> interfacePlotRatio
    ("Ratio",
     "Wether or not to include a ratio plot",
     &GnuplotOutput::_ratio, true, false, false);
  static SwitchOption interfacePlotRatioOn
    (interfacePlotRatio,
     "On",
     "Include a ratio plot, if available",
     true);
  static SwitchOption interfacePlotRatioOff
    (interfacePlotRatio,
     "Off",
     "Do not include a ratio plot",
     false);


  static Switch<GnuplotOutput,bool> interfaceChi2Plot
    ("Chi2",
     "Wether or not include a chi^2 plot",
     &GnuplotOutput::_chi2, true, false, false);
  static SwitchOption interfaceChi2PlotOn
    (interfaceChi2Plot,
     "On",
     "Include a chi^2 plot, if available",
     true);
  static SwitchOption interfaceChi2PlotOff
    (interfaceChi2Plot,
     "Off",
     "Do not include a chi^2 plot",
     false);


  static Parameter<GnuplotOutput,string> interfaceGnuplot
    ("Gnuplot",
     "Location of the gnuplot program, e.g. /usr/bin/gnuplot",
     &GnuplotOutput::_gnuplot, "gnuplot",
     false, false);

}

void GnuplotOutput::put (Histogram2Ptr histo, const Histogram2Options& out, const string& datachannel) {

  if (!histo->haveChannel("MC")) return;

  // initialize to the histogram's name
  Histogram2Output::initialize(histo->name() + ".dat");

  // open a gnuplot file
  ofstream gpfile ((prefix() + histo->name() + ".gp").c_str());

  // get the multiplicity range
  istringstream mrange_in (_subproMult);
  pair<unsigned int, unsigned int> mrange;
  mrange_in >> mrange.first >> mrange.second;

  bool havedata = (datachannel != "");
  if (!histo->haveChannel(datachannel)) havedata = false;

  bool errorbars = (out.plotFlags & HistogramOutput::Errorbars) == HistogramOutput::Errorbars;

  bool xlog = (out.plotFlags & HistogramOutput::Xlog) == HistogramOutput::Xlog;
  bool ylog = (out.plotFlags & HistogramOutput::Ylog) == HistogramOutput::Ylog;

  string ratiochannel = "";
  string chi2channel = "";

  bool ratio = _ratio;
  bool chi2 = _chi2;

  if (!havedata) {
    ratio = false; chi2 = false;
  } else {
    // figure out the ratio and chi2 channel names
    vector<string> allchannels = histo->channels();
    for (vector<string>::iterator ch = allchannels.begin();
	 ch != allchannels.end(); ++ch) {
      if ((*ch).find("delta") != string::npos && ratiochannel == "")
	ratiochannel = *ch;
      if ((*ch).find("chi2") != string::npos && chi2channel == "")
	chi2channel = *ch;
    }
    if (ratiochannel == "") { ratio = false; chi2 = false; }
    if (chi2channel == "") { chi2 = false; }
  }

  bool wantmult = (mrange.first != mrange.second);

  gpfile << "set terminal epslatex color;" << endl
         << "set output '" << histo->name() << ".eps';" << endl
	 << "set size .6,1.1;" << endl
	 << "set title '" << out.title << "';" << endl
	 << "set lmargin .2;" << endl;

  if (ratio || chi2) gpfile << "set multiplot;" << endl;

  if (ratio && chi2) gpfile << "set origin 0,.35;" << endl
			    << "set size .6,.75;" << endl;

  if (ratio && !chi2) gpfile << "set origin 0,.3;" << endl
			     << "set size .6,.8;" << endl;

  if (xlog) gpfile << "set log x;" << endl;
  if (ylog) gpfile << "set log y;" << endl;

  gpfile << "set key samplen .7;" << endl
	 << "set label 1 '{\\definecolor{Gray}{gray}{.5}\\textcolor{Gray}{\\sffamily\\bfseries "
	 << "MC" << "}}' at graph .14,.94 back;" << endl
	 << "set label 2 '" << out.ylabel << "' at graph -.2,.8 rotate by 90 left front;" << endl;
  gpfile << "set label 3 '{\\definecolor{Gray}{gray}{.5}\\textcolor{Gray}{\\tiny\\sffamily "
	 << "Analysis2 powered" << "}}' rotate by 90 at graph -.2,.4 left front;";

  gpfile << "set bars small;" << endl;

  if (ratio || chi2) gpfile << "set format x \"\";" << endl;

  if (!ratio && !chi2) gpfile << "set xlabel '" << out.xlabel << "';" << endl;

  // main plot

  histo->output(currentOStream(),"MC");
  currentOStream() << endl << endl;

  if (mctitle() != "")
    gpfile << "plot '" << histo->name() << ".dat' index 0 u (($1+$2)/2):3:($2-$1) "
	   <<"title '" << mctitle()
	   << "' w boxes lt rgbcolor \"#EFEF00\" fillstyle solid 0.5 noborder\\" << endl;
  else
    gpfile << "plot '" << histo->name() << ".dat' index 0 u (($1+$2)/2):3:($2-$1) "
	   <<"notitle w boxes lt rgbcolor \"#EFEF00\" fillstyle solid 0.5 noborder\\" << endl;

  if (errorbars)
    gpfile << ", '" << histo->name()  << ".dat' index 0 u (($1+$2)/2):3:($3-$4):($3+$4) "
	   << "notitle w errorbars pt -1 lt -1 lc rgbcolor \"#A0A0A0\" lw 1\\" << endl;

  unsigned int currentindex = 0;

  if (wantmult) {
    vector<pair<string,string> > putthem;
    for (unsigned int i = mrange.first; i < mrange.second; ++i) {
      ostringstream multname;
      multname << "MC-" << i;
      ostringstream key;
      key << "$+" << (i-mrange.first) << "\\ {\\rm jet";
      if (i-mrange.first != 1) key << "s";
      key << "}$";
      if (histo->haveChannel(multname.str())) {
	putthem.push_back(make_pair(multname.str(),key.str()));
      } else {
      }
    }


    // check for higher multiplicities and add them together

    ostringstream checkchannel;
    checkchannel << "MC-" << (mrange.second + 1);

    bool havehigher = histo->haveChannel(checkchannel.str());

    if (!havehigher) {
      ostringstream multname; multname << "MC-" << mrange.second;
      ostringstream key; key << "$+" << (mrange.second-mrange.first) << "\\ {\\rm jets}$";
      if (histo->haveChannel(multname.str()))
	putthem.push_back(make_pair(multname.str(),key.str()));
    } else {

      for(unsigned int i = mrange.second; i < 11; ++i) {
	ostringstream multname; multname << "MC-" << i;
	if (i==mrange.second) {
	  histo->insertChannel("MC-largerthanmax",histo->channel(multname.str()));
	} else {
	  if (histo->haveChannel(multname.str()))
	    histo->channel("MC-largerthanmax") += histo->channel(multname.str());
	}
      }

      ostringstream key; key << "$+ \\ge" << (mrange.second-mrange.first) << "\\ {\\rm jets}$";
      putthem.push_back(make_pair("MC-largerthanmax",key.str()));

    }

    unsigned int count = 1;
    for (vector<pair<string,string> >::iterator ch = putthem.begin(); ch != putthem.end(); ++ch) {
      histo->output(currentOStream(),ch->first);
      currentOStream() << endl << endl;
      currentindex += 1;
      count += 1;
      gpfile << ", '" << histo->name()  << ".dat' index " << currentindex
	     << " u (($1+$2)/2):3 title '" << ch->second << "' w histeps lt " << count
	     << " lc rgbcolor \"#000000\"\\" << endl;
    }
  }

  int dataindex = -1;

  if (havedata) {
    currentindex += 1;
    dataindex = currentindex;
    histo->output(currentOStream(),datachannel);
    currentOStream() << endl << endl;
    if (out.datatitle != "")
      gpfile << ", '" << histo->name() << ".dat' index " << currentindex
	     << " u (($1+$2)/2):3:1:2:($3-$4):($3+$4) w xyerrorbars title '"
	     << out.datatitle << "' lt -1 lw 2 pt 7\\" << endl;
    else
      gpfile << ", '" << histo->name() << ".dat' index " << currentindex
	     << " u (($1+$2)/2):3:1:2:($3-$4):($3+$4) w xyerrorbars notitle"
	     << " lt -1 lw 2 pt 7\\" << endl;
  }

  gpfile << ";" << endl;

  // ratio plot

  if (ratio) {

    bool havechi2 = (chi2channel != "");
    double thechi2 = 0.;
    if (havechi2) {
      thechi2 = histo->channel(chi2channel).average(histo->binning()).first;
    }

    currentindex += 1;
    histo->output(currentOStream(),ratiochannel);
    currentOStream() << endl << endl;

    if (chi2) {
      gpfile << "set origin 0,.2;" << endl
	     << "set size .6,.246;" << endl
	     << "set yrange [-.2:.2];" << endl
	     << "set ytics -.1,.1,.1;" << endl;
    } else {
      gpfile << "set origin 0,.0;" << endl
	     << "set size .6,.397;" << endl
	     << "set yrange [-.3:.3];" << endl
	     << "set ytics -.2,.1,.2;" << endl;
    }

    if (ylog) gpfile << "unset log y;" << endl;

    gpfile << "unset title;" << endl;

    if (chi2) {
      gpfile << "set format x \"\";" << endl;
    } else {
      gpfile << "set format x;" << endl;
      gpfile << "set xlabel '" << out.xlabel << "';" << endl;
    }

    gpfile << "unset label 1;" << endl
	   << "unset label 2;" << endl
	   << "unset label 3;" << endl;;

    if (havechi2) {
      gpfile << "set label 1 '{\\tiny $\\chi^2/{\\rm DOF}=" << thechi2
	     << "$}' at graph .2,.87 left front;" << endl;
    }

    gpfile << "plot '" << histo->name() << ".dat' index " << dataindex
	   << " u (($1+$2)/2):($4/$3):($2-$1) notitle w boxes lt rgbcolor \"#00EF40\" fillstyle solid 0.5 noborder, \\"
	   << endl
	   << "'" << histo->name() << ".dat' index " << dataindex
	   << " u (($1+$2)/2):(-$4/$3):($2-$1) notitle w boxes lt rgbcolor \"#00EF40\" fillstyle solid 0.5 noborder,\\"
	   << endl
	   << "'" << histo->name() << ".dat' index 0 u (($1+$2)/2):(0):1:2 notitle w xerrorbars"
	   << " pt -1 lt -1 lc rgbcolor \"#A0A0A0\" lw 1\\" << endl;

    if (errorbars) {
      gpfile << ", '" << histo->name() << ".dat' index 0 u (($1+$2)/2):($4/$3) notitle"
	     << " w histeps lt 3 lc rgbcolor \"#000000\", \\" << endl
	     << "'" << histo->name() << ".dat' index 0 u (($1+$2)/2):(-$4/$3) notitle"
	     << " w histeps lt 3 lc rgbcolor \"#000000\"\\" << endl;
    }

    gpfile << ", '" << histo->name() << ".dat' index " << currentindex
	   << " u (($1+$2)/2):3 notitle w histeps lt -1 lc rgbcolor \"#000000\";"
	   << endl;

  }

  // chi2 plot

  if (chi2) {

    currentindex += 1;
    histo->output(currentOStream(),chi2channel);
    currentOStream() << endl << endl;

    gpfile << "set origin 0,0;" << endl
	   << "set size .6,.296;" << endl
	   << "set format x;" << endl
	   << "unset label 1;" << endl
	   << "set xlabel '" << out.xlabel << "';" << endl
	   << "set yrange [0:10];" << endl
	   << "set ytics 0,4,8;" << endl
	   << "plot '" << histo->name() << ".dat' index " << currentindex 
	   << " u (($1+$2)/2):3:($2-$1) notitle w boxes lt rgbcolor \"#A0A0A0\" fillstyle solid 0.8 noborder;"
	   << endl;

  }

  gpfile << "reset;" << endl;

  if (ratio || chi2)
    gpfile << "unset multiplot;" << endl
	   << "reset;" << endl;

  _allHistos += 1;

  _TeXfile << "\\input{" << histo->name() << "}";

  if (_allHistos % 2 == 0)
    _TeXfile << endl << endl;
  else
    _TeXfile << "\\hfill";

  _makefile << "\t" << _gnuplot << " " << histo->name() << ".gp" << endl;

}

