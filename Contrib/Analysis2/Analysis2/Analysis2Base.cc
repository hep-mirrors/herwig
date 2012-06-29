// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#include "Analysis2Base.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#include "ThePEG/Repository/Repository.h"

#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Analysis2Base.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/SubProcess.h"

#include "ThePEG/Utilities/StringUtils.h"

#include <cassert>

using namespace Analysis2;

Analysis2Base::~Analysis2Base() {}

void Analysis2Base::analyze(tEventPtr event, long ieve, int loop, int state) {

  _lastEvent = event;
  if (_minMult == 0) 
    _minMult = _lastEvent->primarySubProcess()->outgoing().size();
  else if (_lastEvent->primarySubProcess()->outgoing().size() < _minMult)
    _minMult = _lastEvent->primarySubProcess()->outgoing().size();
  if (_maxMult == 0)
    _maxMult = _lastEvent->primarySubProcess()->outgoing().size();
  else if (_lastEvent->primarySubProcess()->outgoing().size() > _maxMult)
    _maxMult = _lastEvent->primarySubProcess()->outgoing().size();

  _eventExtractor->use(event);

  AnalysisHandler::analyze(event, ieve, loop, state);

}

LorentzRotation Analysis2Base::transform(tEventPtr) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void Analysis2Base::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.

}

void Analysis2Base::analyze(tPPtr) {}

void Analysis2Base::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _output << _histograms << _outputOptions << _datachannels << _normalisation
     << _bookPerSubprocess << _lastEvent << _minMult << _maxMult
     << _parallel << _jetFinder << _eventExtractor;
}

void Analysis2Base::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _output >> _histograms >> _outputOptions >> _datachannels >> _normalisation
     >> _bookPerSubprocess >> _lastEvent >> _minMult >> _maxMult
     >> _parallel >> _jetFinder >> _eventExtractor;
}

AbstractClassDescription<Analysis2Base> Analysis2Base::initAnalysis2Base;
// Definition of the static class description member.

void Analysis2Base::Init() {

  static ClassDocumentation<Analysis2Base> documentation
    ("Analysis2Base is the base class for all analysis handlers"
     " used by Analysis2");


  static Reference<Analysis2Base,Histogram2Output> interfaceHistogramOutput
    ("HistogramOutput",
     "Set the HistogramOutput object",
     &Analysis2Base::_output, false, false, true, false, false);

  static Reference<Analysis2Base,JetFinder> interfaceJetFinder
    ("JetFinder",
     "Set the JetFinder object",
     &Analysis2Base::_jetFinder, false, false, true, true, false);

  static Reference<Analysis2Base,EventExtractor> interfaceEventExtractor
    ("EventExtractor",
     "Set the EventExtractor object",
     &Analysis2Base::_eventExtractor, false, false, true, true, false);

  static Switch<Analysis2Base,bool> interfaceBookPerSubprocess
    ("BookPerSubprocess",
     "Book events per subprocess multiplicity",
     &Analysis2Base::_bookPerSubprocess, false, false, false);
  static SwitchOption interfaceBookPerSubprocessOn
    (interfaceBookPerSubprocess,
     "Yes",
     "Book per subprocess multiplicity",
     true);
  static SwitchOption interfaceBookPerSubprocessOff
    (interfaceBookPerSubprocess,
     "No",
     "Do not book per subprocess multiplicity",
     false);


  static Switch<Analysis2Base,bool> interfaceParallel
    ("Parallel",
     "Wether or not his analysis is carried out as part of a parallel run.",
     &Analysis2Base::_parallel, false, false, false);
  static SwitchOption interfaceParallelOn
    (interfaceParallel,
     "Yes",
     "This analysis is carried out on part of a parallel run.",
     true);
  static SwitchOption interfaceParallelOff
    (interfaceParallel,
     "No",
     "This analysis is carried out on a serial run",
     false);


  static Command<Analysis2Base> interfaceStartCombine
    ("StartCombine",
     "Start combining parallel runs",
     &Analysis2Base::startCombine, false);

  static Command<Analysis2Base> interfaceCombine
    ("Combine",
     "Combine parallel runs",
     &Analysis2Base::combine, false);

  static Command<Analysis2Base> interfaceFinishCombine
    ("FinishCombine",
     "Finish combining parallel runs",
     &Analysis2Base::finishCombine, false);

}

void Analysis2Base::insert (const string& name,
			    const string& instring,
			    Histogram2Options& options) {
  string theinstring = StringUtils::stripws(instring);
  string error = "Analysis2Base::insert : warning : "
    + instring + " is not a valid options tag.\nNo events for the observable "
    + name + " will be booked.\n";
  if (instring.find("<options") == string::npos) {
    Repository::clog() << error;
    return;
  }
  map<string,string> attributes = StringUtils::xmlAttributes("options",theinstring);

  if (attributes.find("title") != attributes.end())
    options.title = attributes.find("title")->second;
  else
    options.title = "";

  if (attributes.find("xlabel") != attributes.end())
    options.xlabel = attributes.find("xlabel")->second;
  else
    options.xlabel = "";

  if (attributes.find("ylabel") != attributes.end())
    options.ylabel = attributes.find("ylabel")->second;
  else
    options.ylabel = "";

  string norm = "";
  if (attributes.find("norm") != attributes.end())
    norm = attributes.find("norm")->second;

  int normMode = NormaliseToXSec;

  if (norm == "none")
    normMode = NoNormalisation;
  else if (norm == "unity")
    normMode = NormaliseToUnity;
  else if (norm == "xsec")
    normMode = NormaliseToXSec;
  else if (norm == "data")
    normMode = NormaliseToData;
  else {
    Repository::clog() << error;
    return;
  }

  string dataname = "";
  string datafile = "";

  string blow = "";
  string bhigh = "";
  string nbins = "";

  if (attributes.find("datafile") != attributes.end())
    datafile = attributes.find("datafile")->second;
  if (attributes.find("dataname") != attributes.end()) {
    dataname = attributes.find("dataname")->second;
  }
  if (attributes.find("lowerbound") != attributes.end())
    blow = attributes.find("lowerbound")->second;
  if (attributes.find("upperbound") != attributes.end())
    bhigh = attributes.find("upperbound")->second;
  if (attributes.find("numbins") != attributes.end())
    nbins = attributes.find("numbins")->second;

  bool fromdata = true;
  if (blow != "" && bhigh != "" && nbins != "") {
    fromdata=false;
    if (datafile != "" || dataname != "") {
      Repository::clog() << error;
      return;
    }
  }

  if (fromdata && datafile == "") {
    Repository::clog() << error;
    return;
  }

  bool log = false;

  if (attributes.find("log") != attributes.end())
    if (attributes.find("log")->second == "true")
      log = true;


  if (!fromdata) {

    istringstream inblow (blow);
    istringstream inbhigh (bhigh);
    istringstream innbins (nbins);

    double binlower; inblow >> binlower;
    double binhigher; inbhigh >> binhigher;
    unsigned int numbins; innbins >> numbins;

    if (binhigher <= binlower || numbins <= 0) {
      Repository::clog() << error;
      return;
    }
    if (!log)
      insertObservable(name,new_ptr(Histogram2(binlower,binhigher,numbins,"MC")),options);
    else {
      if (binlower == 0) {
	Repository::clog() << error;
	return;
      }
      double c = log10(binhigher/binlower)/numbins;
      vector<pair<double,double> > binning;
      double bi = binlower;
      for(unsigned int i = 0; i< numbins; ++i) {
	binning.push_back(make_pair(bi*pow(10.,c*i),bi*pow(10.,c*(i+1))));
      }
      insertObservable(name,new_ptr(Histogram2(binning,"MC")),options);
    }

  } else {

    options.datatitle = dataname;
    insertObservable(name,dataname,datafile,options);

  }

  _normalisation.insert(make_pair(name,normMode));

}

void Analysis2Base::book (double evt, const string& name, double weight) {
  map<string,Histogram2Ptr>::iterator h = _histograms.find(name);
  if (h != _histograms.end()) {
    assert(h->second);
    assert(lastEvent());
    if (!_bookPerSubprocess) {
      h->second->book("MC",evt,weight);
    } else {
      // figure out the final state multiplicity
      unsigned int mult = lastEvent()->primarySubProcess()->outgoing().size();
      ostringstream chName; chName << "MC-" << mult;
      if (!h->second->haveChannel(chName.str())) {
	// insert the channel
	h->second->insertChannel(chName.str(),HistogramChannel(h->second->numberOfBins()));
      }
      h->second->book(chName.str(),evt,weight);
    }
  }
}

void Analysis2Base::finish (const string& name,
			    int norm,
			    bool combined) {

  Histogram2Ptr theHisto;
  if (_histograms.find(name) != _histograms.end()) {
    theHisto = _histograms.find(name)->second;
  } else {
    generator()->log() << "Analysis2Base::finish : could not finish "
		       << name << " : No such observable." << endl;
    return;
  }

  string data = "";

  if (_datachannels.find(name) != _datachannels.end())
    data = _datachannels.find(name)->second;

  vector<string> allchannels = theHisto->channels();
  for (vector<string>::iterator c = allchannels.begin(); c != allchannels.end(); ++c)
    if (*c != data)
      theHisto->finish(*c);

  // return here, if part of a parallel run
  if (_parallel) return;

  int normMode = NormaliseToXSec;

  if (norm == FromMap) {
    if (_normalisation.find(name) != _normalisation.end())
      normMode = _normalisation.find(name)->second;
  } else {
    normMode = norm;
  }

  if (!combined) {
    theHisto->xSec(generator()->currentEventHandler()->histogramScale());
  }

  if (_outputOptions[name].differential)
    for (vector<string>::iterator c = allchannels.begin(); c != allchannels.end(); ++c)
      if (*c != data)
	theHisto->differential(*c);

  // normalize

  if (!_bookPerSubprocess) {

    if (theHisto->haveChannel("MC") && normMode != NoNormalisation) {

      if (normMode == NormaliseToData && !theHisto->haveChannel(data))
	normMode = NormaliseToXSec;
      switch (normMode) {
      case NormaliseToUnity:
	theHisto->normalise("MC");
	break;
      case NormaliseToXSec:
	theHisto->normaliseToCrossSection("MC");
	break;
      case NormaliseToData:
	theHisto->normalise("MC",data);
	break;
      }

    }

  } else {

    // booking per subprocess

    // sum up the channels
    
    HistogramChannel sum (theHisto->numberOfBins());
    
    for (unsigned int i = minMult(); i<= maxMult(); ++i) {
      ostringstream chName; chName << "MC-" << i;    
      sum += theHisto->channel(chName.str());
    }
    
    theHisto->removeChannel("MC");
    theHisto->insertChannel("MC", sum);
    
    // normalize
    
    // booking per subprocess needs to know the normalisation factor

    if (normMode != NoNormalisation) {
    
      pair<double,double> fromIntegral;
      pair<double,double> toIntegral;
      
      fromIntegral = theHisto->integrate("MC");
      if (normMode == NormaliseToData && !theHisto->haveChannel(data))
	normMode = NormaliseToXSec;
      switch (normMode) {
      case NormaliseToUnity:
	toIntegral = make_pair(1.,0.);
	break;
      case NormaliseToXSec:
	toIntegral = make_pair(generator()->currentEventHandler()->histogramScale()/nanobarn,0.);
	break;
      case NormaliseToData:
	toIntegral = theHisto->integrate(data);
	break;
      }
      
      theHisto->channel("MC") /= fromIntegral;
      theHisto->channel("MC") *= toIntegral;
      
      for (unsigned int i = minMult(); i<= maxMult(); ++i) {
	ostringstream chName; chName << "MC-" << i;    
	theHisto->channel(chName.str()) /= fromIntegral;
	theHisto->channel(chName.str()) *= toIntegral;
      }

    }

  }

  // compute statistical tests, if data present

  if (theHisto->haveChannel("MC") && theHisto->haveChannel(data)) {

    // ratio
    theHisto->insertChannel("delta (MC,"+data+")",theHisto->channel("MC").delta(theHisto->channel(data)));
    
    // chi2
    theHisto->insertChannel("chi2 (MC,"+data+")",theHisto->channel("MC").chi2(theHisto->channel(data)));
    if (!combined)
      generator()->log() << name << " : chi2/DOF (MC,"+data+") = "
			 << theHisto->channel("chi2 (MC,"+data+")").average(theHisto->binning()).first
			 << endl;
    else
      Repository::clog() << name << " : chi2/DOF (MC,"+data+") = "
			 << theHisto->channel("chi2 (MC,"+data+")").average(theHisto->binning()).first
			 << endl;
  }

  // output to plot file

  if (_outputOptions.find(name) == _outputOptions.end())
    output()->put(theHisto,Histogram2Options(),data);
  else
    output()->put(theHisto,_outputOptions.find(name)->second,data);

}

void Analysis2Base::combineObservable (const string& prefix,
				       const string& name,
				       unsigned int numRuns,
				       int norm) {

  int normMode = NormaliseToXSec;

  if (norm == -1) {
    if (_normalisation.find(name) != _normalisation.end())
      normMode = _normalisation.find(name)->second;
  }

  string dataChannel = "";

  if(_datachannels.find(name) != _datachannels.end())
    dataChannel = _datachannels.find(name)->second;

  Histogram2Ptr sum = new_ptr(Histogram2());
  sum->combine(prefix,name,numRuns,dataChannel,"MC");

  insertObservable(name,sum);
  finish(name,normMode,true);

  _histograms.find(name)->second->store(name);

}

string Analysis2Base::combine (string argstring) {
  string usage = "Combine [run prefix] [number of runs] [observable name] \\\n"
                 "        {normalization mode : none, unity, xsec, data}";
  StringUtils::StringVector args = StringUtils::split(argstring," ");
  if (args.size() < 3) return usage;
  string prefix = args[0];
  Repository::clog() << prefix << endl;
  istringstream nrunin (args[1]);
  unsigned int nruns; nrunin >> nruns;
  string name = args[2];
  int normMode = FromMap;
  if (args.size() == 4) {
    if (args[3] == "none")
      normMode = NoNormalisation;
    else if (args[3] == "unity")
      normMode = NormaliseToUnity;
    else if (args[3] == "xsec")
      normMode = NormaliseToXSec;
    else if (args[3] == "data")
      normMode = NormaliseToData;
    else normMode = FromMap;
  }

  if (name == "*") {
    for (map<string,Histogram2Ptr>::iterator obs = _histograms.begin();
	 obs != _histograms.end(); ++obs)
      combineObservable(prefix,obs->first,nruns,normMode);
  }
  else
    combineObservable(prefix,name,nruns,normMode);

  return "";
  
}
