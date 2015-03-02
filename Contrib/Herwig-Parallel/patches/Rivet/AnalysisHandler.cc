// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/ParticleName.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/Beam.hh"
#include "LWH/AIManagedObject.h"

using namespace AIDA;

namespace Rivet {


  AnalysisHandler::AnalysisHandler(const string& runname)
    : _runname(runname), _numEvents(0),
      _sumOfWeights(0.0), _xs(-1.0),
      _initialised(false), _ignoreBeams(false)
  {
    _theAnalysisFactory.reset( createAnalysisFactory() );
    _setupFactories();
  }


  AnalysisHandler::AnalysisHandler(const string& basefilename,
                                   const string& runname, HistoFormat storetype)
    : _runname(runname), _numEvents(0),
      _sumOfWeights(0.0), _xs(-1.0),
      _initialised(false), _ignoreBeams(false)
  {
    cerr << "AnalysisHandler(basefilename, runname, format) constructor is deprecated: "
         << "please migrate your code to use the one-arg constructor" << endl;
    _theAnalysisFactory.reset( createAnalysisFactory() );
    _setupFactories(basefilename, storetype);
  }


  AnalysisHandler::~AnalysisHandler()
  {  }


  Log& AnalysisHandler::getLog() const {
    return Log::getLog("Rivet.Analysis.Handler");
  }


  void AnalysisHandler::init(const GenEvent& ge) {
    assert(!_initialised);
    setRunBeams(Rivet::beams(ge));
    MSG_DEBUG("Initialising the analysis handler");
    _numEvents = 0;
    _sumOfWeights = 0.0;

    // Check that analyses are beam-compatible, and remove those that aren't
    const size_t num_anas_requested = analysisNames().size();
    vector<string> anamestodelete;
    foreach (const AnaHandle a, _analyses) {
      if ((!a->isCompatible(beams())) && (!_ignoreBeams)) {
        //MSG_DEBUG(a->name() << " requires beams " << a->requiredBeams() << " @ " << a->requiredEnergies() << " GeV");
        anamestodelete.push_back(a->name());
      }
    }
    foreach (const string& aname, anamestodelete) {
      MSG_WARNING("Analysis '" << aname << "' is incompatible with the provided beams: removing");
      removeAnalysis(aname);
    }
    if (num_anas_requested > 0 && analysisNames().size() == 0) {
      cerr << "All analyses were incompatible with the first event's beams\n"
           << "Exiting, since this probably wasn't intentional!" << endl;
      exit(1);
    }

    // Warn if any analysis' status is not unblemished
    foreach (const AnaHandle a, analyses()) {
      if (toUpper(a->status()) == "PRELIMINARY") {
        MSG_WARNING("Analysis '" << a->name() << "' is preliminary: be careful, it may change and/or be renamed!");
      } else if (toUpper(a->status()) == "OBSOLETE") {
        MSG_WARNING("Analysis '" << a->name() << "' is obsolete: please update!");
      } else if (toUpper(a->status()).find("UNVALIDATED") != string::npos) {
        MSG_WARNING("Analysis '" << a->name() << "' is unvalidated: be careful, it may be broken!");
      }
    }

    // Initialize the remaining analyses
    foreach (AnaHandle a, _analyses) {
      MSG_DEBUG("Initialising analysis: " << a->name());
      try {
        // Allow projection registration in the init phase onwards
        a->_allowProjReg = true;
        a->init();
        //MSG_DEBUG("Checking consistency of analysis: " << a->name());
        //a->checkConsistency();
      } catch (const Error& err) {
        cerr << "Error in " << a->name() << "::init method: " << err.what() << endl;
        exit(1);
      }
      MSG_DEBUG("Done initialising analysis: " << a->name());
    }
    _initialised = true;
    MSG_DEBUG("Analysis handler initialised");
  }

  
  void AnalysisHandler::finishEvent(){
    foreach (AnaHandle a, _analyses){
      a->finishEvent();
    }
  }

  void AnalysisHandler::analyze(const GenEvent& ge) {
    // Call init with event as template if not already initialised
    if (!_initialised) {
      init(ge);
    }
    // Proceed with event analysis
    assert(_initialised);
    // Ensure that beam details match those from first event
    const PdgIdPair beams = Rivet::beamIds(ge);
    const double sqrts = Rivet::sqrtS(ge);
    if (!compatible(beams, _beams) || !fuzzyEquals(sqrts, sqrtS())) {
      cerr     << "Event beams mismatch: "
               << toBeamsString(beams) << " @ " << sqrts/GeV << " GeV" << " vs. first beams "
               << this->beams() << " @ " << this->sqrtS()/GeV << " GeV" << endl;
      exit(1);
    }


    Event event(ge);
    _numEvents++;
    // Weights
    const double weight = event.weight();
    _sumOfWeights += weight;
    MSG_DEBUG("Event #" << _numEvents << " weight = " << weight);
    #ifdef HEPMC_HAS_CROSS_SECTION
    if (ge.cross_section()) {
      const double xs = ge.cross_section()->cross_section();
      setCrossSection(xs);
    }
    #endif
    foreach (AnaHandle a, _analyses) {
      //MSG_DEBUG("About to run analysis " << a->name());
      try {
        a->analyze(event);
      } catch (const Error& err) {
        cerr     << "Error in " << a->name() << "::analyze method: "
                 << err.what() << endl;
        exit(1);
      }
      //MSG_DEBUG("Finished running analysis " << a->name());
    }
  }

  void AnalysisHandler::setAttempts(long attempts) {
    foreach (AnaHandle a, _analyses) {
      if(!a->histos.empty())
        for(vector<AIDA::IHistogram1D *>::iterator it= a->histos.begin();it!= a->histos.end();it++)
          if((*it)) (*it)->setAttempts(attempts);
    }
  }

  void AnalysisHandler::finalize() {
    assert(_initialised);
    MSG_INFO("Finalising analyses");
    foreach (AnaHandle a, _analyses) {
      try {
        a->finalize();
      } catch (const Error& err) {
        cerr     << "Error in " << a->name() << "::finalize method: "
                 << err.what() << endl;
        exit(1);
      }
    }

    // Print out number of events processed
    MSG_INFO("Processed " << _numEvents << " event" << (_numEvents == 1 ? "" : "s"));

    // Change AIDA histos into data point sets
    MSG_DEBUG("Converting histograms to scatter plots");
    assert(_theTree != 0);
    _normalizeTree(tree());

    // Delete analyses
    MSG_DEBUG("Deleting analyses");
    _analyses.clear();

    // Print out MCnet boilerplate
    cout << endl;
    cout << "The MCnet usage guidelines apply to Rivet: see http://www.montecarlonet.org/GUIDELINES" << endl;
    cout << "Please acknowledge plots made with Rivet analyses, and cite arXiv:1003.0694 (http://arxiv.org/abs/1003.0694)" << endl;
  }


  AnalysisHandler& AnalysisHandler::addAnalysis(const string& analysisname) {
    // Check for a duplicate analysis
    /// @todo Might we want to be able to run an analysis twice, with different params?
    ///       Requires avoiding histo tree clashes, i.e. storing the histos on the analysis objects.
    foreach (const AnaHandle& a, _analyses) {
      if (a->name() == analysisname) {
        MSG_WARNING("Analysis '" << analysisname << "' already registered: skipping duplicate");
        return *this;
      }
    }
    AnaHandle analysis( AnalysisLoader::getAnalysis(analysisname) );
    if (analysis.get() != 0) { // < Check for null analysis.
      MSG_DEBUG("Adding analysis '" << analysisname << "'");
      analysis->_analysishandler = this;
      _analyses.insert(analysis);
    } else {
      MSG_WARNING("Analysis '" << analysisname << "' not found.");
    }
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeAnalysis(const string& analysisname) {
    shared_ptr<Analysis> toremove;
    foreach (const AnaHandle a, _analyses) {
      if (a->name() == analysisname) {
        toremove = a;
        break;
      }
    }
    if (toremove.get() != 0) {
      MSG_DEBUG("Removing analysis '" << analysisname << "'");
      _analyses.erase(toremove);
    }
    return *this;
  }


  void AnalysisHandler::_setupFactories(const string& basefilename, HistoFormat storetype) {
    string filename(basefilename), storetypestr("");
    if (storetype == AIDAML) {
      if (!endsWith(filename, ".aida")) filename += ".aida";
      storetypestr = "xml";
    } else if (storetype == FLAT) {
      if (!endsWith(filename, ".data")) filename += ".data";
      storetypestr = "flat";
    } else if (storetype == ROOT) {
      if (!endsWith(filename, ".root")) filename += ".root";
      storetypestr = "root";
    }
    _theTreeFactory = _theAnalysisFactory->createTreeFactory();
    _theTree = _theTreeFactory->create(filename, storetypestr, false, true);
    _theHistogramFactory = _theAnalysisFactory->createHistogramFactory(tree());
    _theDataPointSetFactory = _theAnalysisFactory->createDataPointSetFactory(tree());
  }


  void AnalysisHandler::_setupFactories() {
    _theTreeFactory = _theAnalysisFactory->createTreeFactory();
    _theTree = _theTreeFactory->create();
    _theHistogramFactory = _theAnalysisFactory->createHistogramFactory(tree());
    _theDataPointSetFactory = _theAnalysisFactory->createDataPointSetFactory(tree());
  }


  void AnalysisHandler::commitData() {
    tree().commit();
  }


  void AnalysisHandler::writeData(const string& filename) {
    tree().commit(filename);
  }


  void AnalysisHandler::_normalizeTree(ITree& tree) {
    const vector<string> paths = tree.listObjectNames("/", true); // args set recursive listing
    MSG_TRACE("Number of objects in AIDA tree = " << paths.size());
    const string tmpdir = "/RivetNormalizeTmp";
    tree.mkdir(tmpdir);
    foreach (const string& path, paths) {

      IManagedObject* hobj = tree.find(path);
      if (hobj) {

        // Try to cast to specific histogram types
        const IProfile1D* prof = dynamic_cast<IProfile1D*>(hobj);
        const IHistogram1D* histo = (prof) ? 0 : dynamic_cast<IHistogram1D*>(hobj);
        const IHistogram2D* histo2 = (prof || histo) ? 0 : dynamic_cast<IHistogram2D*>(hobj);
        if (!(histo || histo2 || prof)) {
          MSG_TRACE("Could not find the type of histo for " << path << ": it's probably already a DPS");
          continue;
        }

        // AIDA path mangling
        const size_t lastslash = path.find_last_of("/");
        const string basename = path.substr(lastslash+1, path.length() - (lastslash+1));
        const string tmppath = tmpdir + "/" + basename;

        // If it's a normal histo:
        tree.mv(path, tmpdir);
        if (histo) {
          MSG_TRACE("Converting histo " << path << " to DPS");
          IHistogram1D* tmphisto = dynamic_cast<IHistogram1D*>(tree.find(tmppath));
          if (tmphisto) datapointsetFactory().create(path, *tmphisto);
        }
        // If it's a 2D histo:
        else if (histo2) {
          MSG_TRACE("Converting 2D histo " << path << " to DPS");
          IHistogram2D* tmphisto2 = dynamic_cast<IHistogram2D*>(tree.find(tmppath));
          if (tmphisto2) datapointsetFactory().create(path, *tmphisto2);
        }
        // If it's a profile histo:
        else if (prof) {
          MSG_TRACE("Converting profile histo " << path << " to DPS");
          IProfile1D* tmpprof = dynamic_cast<IProfile1D*>(tree.find(tmppath));
          if (tmpprof) datapointsetFactory().create(path, *tmpprof);
        }
        tree.rm(tmppath);

      }

    }
    tree.rmdir(tmpdir);
  }


  string AnalysisHandler::runName() const { return _runname; }
  size_t AnalysisHandler::numEvents() const { return _numEvents; }
  double AnalysisHandler::sumOfWeights() const { return _sumOfWeights; }


  void AnalysisHandler::setSumOfWeights(const double& sum) {
    _sumOfWeights=sum;
  }


  std::vector<std::string> AnalysisHandler::analysisNames() const {
    std::vector<std::string> rtn;
    foreach (AnaHandle a, _analyses) {
      rtn.push_back(a->name());
    }
    return rtn;
  }


  AnalysisHandler& AnalysisHandler::addAnalyses(const std::vector<std::string>& analysisnames) {
    foreach (const string& aname, analysisnames) {
      //MSG_DEBUG("Adding analysis '" << aname << "'");
      addAnalysis(aname);
    }
    return *this;
  }


  AnalysisHandler& AnalysisHandler::removeAnalyses(const std::vector<std::string>& analysisnames) {
    foreach (const string& aname, analysisnames) {
      removeAnalysis(aname);
    }
    return *this;
  }



  AIDA::IAnalysisFactory& AnalysisHandler::analysisFactory() {
    return *_theAnalysisFactory;
  }


  AIDA::ITree& AnalysisHandler::tree() {
    return *_theTree;
  }


  AIDA::IHistogramFactory& AnalysisHandler::histogramFactory() {
    return *_theHistogramFactory;
  }


  AIDA::IDataPointSetFactory& AnalysisHandler::datapointsetFactory() {
    return *_theDataPointSetFactory;
  }


  bool AnalysisHandler::needCrossSection() const {
    bool rtn = false;
    foreach (const AnaHandle a, _analyses) {
      if (!rtn) rtn = a->needsCrossSection();
      if (rtn) break;
    }
    return rtn;
  }


  AnalysisHandler& AnalysisHandler::setCrossSection(double xs) {
    _xs = xs;
    foreach (AnaHandle a, _analyses) {
      a->setCrossSection(xs);
    }
    return *this;
  }


  bool AnalysisHandler::hasCrossSection() const {
    return (!std::isnan(crossSection()));
  }


  AnalysisHandler& AnalysisHandler::addAnalysis(Analysis* analysis) {
    analysis->_analysishandler = this;
    _analyses.insert(AnaHandle(analysis));
    return *this;
  }


  PdgIdPair AnalysisHandler::beamIds() const {
    return Rivet::beamIds(beams());
  }


  double AnalysisHandler::sqrtS() const {
    return Rivet::sqrtS(beams());
  }

  void AnalysisHandler::setIgnoreBeams(bool ignore) {
    _ignoreBeams=ignore;
  }


}
