// -*- C++ -*-
#include "Rivet/Rivet.hh"
#include "Rivet/RivetAIDA.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "LWH/AIManagedObject.h"
using namespace AIDA;

namespace Rivet {


  Analysis::Analysis(const string& name)
    : _crossSection(-1.0),
      _gotCrossSection(false),
      _analysishandler(NULL),
      _madeHistoDir(false)
  {
    ProjectionApplier::_allowProjReg = false;
    _defaultname = name;

    AnalysisInfo* ai = AnalysisInfo::make(name);
    assert(ai != 0);
    _info.reset(ai);
    assert(_info.get() != 0);
  }


  IAnalysisFactory& Analysis::analysisFactory() {
    return handler().analysisFactory();
  }


  ITree& Analysis::tree() {
    return handler().tree();
  }


  IHistogramFactory& Analysis::histogramFactory() {
    return handler().histogramFactory();
  }


  IDataPointSetFactory& Analysis::datapointsetFactory() {
    return handler().datapointsetFactory();
  }


  double Analysis::sqrtS() const {
    return handler().sqrtS();
  }

  const ParticlePair& Analysis::beams() const {
    return handler().beams();
  }

  const PdgIdPair Analysis::beamIds() const {
    return handler().beamIds();
  }


  const string Analysis::histoDir() const {
    /// @todo This doesn't change: calc and cache at first use!
    string path = "/" + name();
    if (handler().runName().length() > 0) {
      path = "/" + handler().runName() + path;
    }
    while (find_first(path, "//")) {
      replace_all(path, "//", "/");
    }
    return path;
  }


  const string Analysis::histoPath(const string& hname) const {
    const string path = histoDir() + "/" + hname;
    return path;
  }


  const string Analysis::histoPath(size_t datasetId, size_t xAxisId, size_t yAxisId) const {
    return histoDir() + "/" + makeAxisCode(datasetId, xAxisId, yAxisId);
  }


  const string Analysis::makeAxisCode(size_t datasetId, size_t xAxisId, size_t yAxisId) const {
    stringstream axisCode;
    axisCode << "d";
    if (datasetId < 10) axisCode << 0;
    axisCode << datasetId;
    axisCode << "-x";
    if (xAxisId < 10) axisCode << 0;
    axisCode << xAxisId;
    axisCode << "-y";
    if (yAxisId < 10) axisCode << 0;
    axisCode << yAxisId;
    return axisCode.str();
  }


  Log& Analysis::getLog() const {
    string logname = "Rivet.Analysis." + name();
    return Log::getLog(logname);
  }


  size_t Analysis::numEvents() const {
    return handler().numEvents();
  }


  double Analysis::sumOfWeights() const {
    return handler().sumOfWeights();
  }


  ///////////////////////////////////////////


  bool Analysis::isCompatible(const ParticlePair& beams) const {
    return isCompatible(beams.first.pdgId(),  beams.second.pdgId(),
                        beams.first.energy(), beams.second.energy());
  }


  bool Analysis::isCompatible(PdgId beam1, PdgId beam2, double e1, double e2) const {
    PdgIdPair beams(beam1, beam2);
    pair<double,double> energies(e1, e2);
    return isCompatible(beams, energies);
  }


  bool Analysis::isCompatible(const PdgIdPair& beams, const pair<double,double>& energies) const {
    // First check the beam IDs
    bool beamIdsOk = false;
    foreach (const PdgIdPair& bp, requiredBeams()) {
      if (compatible(beams, bp)) {
        beamIdsOk =  true;
        break;
      }
    }
    if (!beamIdsOk) return false;

    // Next check that the energies are compatible (within 1%, to give a bit of UI forgiveness)
    bool beamEnergiesOk = requiredEnergies().size()>0 ? false : true;
    typedef pair<double,double> DoublePair;
    foreach (const DoublePair& ep, requiredEnergies()) {
      if ((fuzzyEquals(ep.first, energies.first, 0.01) && fuzzyEquals(ep.second, energies.second, 0.01)) ||
          (fuzzyEquals(ep.first, energies.second, 0.01) && fuzzyEquals(ep.second, energies.first, 0.01))) {
        beamEnergiesOk =  true;
        break;
      }
    }
    return beamEnergiesOk;

    /// @todo Need to also check internal consistency of the analysis'
    /// beam requirements with those of the projections it uses.
  }


  ///////////////////////////////////////////


  Analysis& Analysis::setCrossSection(double xs) {
    _crossSection = xs;
    _gotCrossSection = true;
    return *this;
  }

  double Analysis::crossSection() const {
    if (!_gotCrossSection || std::isnan(_crossSection)) {
      string errMsg = "You did not set the cross section for the analysis " + name();
      throw Error(errMsg);
    }
    return _crossSection;
  }

  double Analysis::crossSectionPerEvent() const {
    const double sumW = sumOfWeights();
    assert(sumW != 0.0);
    return _crossSection / sumW;
  }



  ////////////////////////////////////////////////////////////
  // Histogramming


  void Analysis::_cacheBinEdges() const {
    _cacheXAxisData();
    if (_histBinEdges.empty()) {
      MSG_TRACE("Getting histo bin edges from AIDA for paper " << name());
      _histBinEdges = getBinEdges(_dpsData);
    }
  }


  void Analysis::_cacheXAxisData() const {
    if (_dpsData.empty()) {
      MSG_TRACE("Getting DPS x-axis data from AIDA for paper " << name());
      _dpsData = getDPSXValsErrs(name());
    }
  }


  const BinEdges& Analysis::binEdges(const string& hname) const {
    _cacheBinEdges();
    MSG_TRACE("Using histo bin edges for " << name() << ":" << hname);
    const BinEdges& edges = _histBinEdges.find(hname)->second;
    if (getLog().isActive(Log::TRACE)) {
      stringstream edges_ss;
      foreach (const double be, edges) {
        edges_ss << " " << be;
      }
      MSG_TRACE("Edges:" << edges_ss.str());
    }
    return edges;
  }


  const BinEdges& Analysis::binEdges(size_t datasetId, size_t xAxisId, size_t yAxisId) const {
    const string hname = makeAxisCode(datasetId, xAxisId, yAxisId);
    return binEdges(hname);
  }


  BinEdges Analysis::logBinEdges(size_t nbins, double lower, double upper) {
    return logspace(nbins, lower, upper);
  }


  IHistogram1D* Analysis::bookHistogram1D(size_t datasetId, size_t xAxisId,
                                          size_t yAxisId, const string& title,
                                          const string& xtitle, const string& ytitle)
  {
    const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
    return bookHistogram1D(axisCode, title, xtitle, ytitle);
  }


  IHistogram1D* Analysis::bookHistogram1D(const string& hname, const string& title,
                                          const string& xtitle, const string& ytitle)
  {
    // Get the bin edges (only read the AIDA file once)
    const BinEdges edges = binEdges(hname);
    _makeHistoDir();
    const string path = histoPath(hname);
    if (path.find(" ") != string::npos) {
      throw Error("Histogram path '" + path + "' is invalid: spaces are not permitted in paths");
    }
    IHistogram1D* hist = histogramFactory().createHistogram1D(path, title, edges);
    MSG_TRACE("Made histogram " << hname <<  " for " << name());
    hist->setXTitle(xtitle);
    hist->setYTitle(ytitle);
    
    histos.push_back(hist);
    return hist;
  }


  IHistogram1D* Analysis::bookHistogram1D(const string& hname,
                                          size_t nbins, double lower, double upper,
                                          const string& title,
                                          const string& xtitle, const string& ytitle) {
    _makeHistoDir();
    const string path = histoPath(hname);
    if (path.find(" ") != string::npos) {
      throw Error("Histogram path '" + path + "' is invalid: spaces are not permitted in paths");
    }
    IHistogram1D* hist = histogramFactory().createHistogram1D(path, title, nbins, lower, upper);
    MSG_TRACE("Made histogram " << hname <<  " for " << name());
    hist->setXTitle(xtitle);
    hist->setYTitle(ytitle);
    histos.push_back(hist);
    return hist;
  }


  IHistogram1D* Analysis::bookHistogram1D(const string& hname,
                                          const vector<double>& binedges,
                                          const string& title,
                                          const string& xtitle, const string& ytitle) {
    _makeHistoDir();
    const string path = histoPath(hname);
    if (path.find(" ") != string::npos) {
      throw Error("Histogram path '" + path + "' is invalid: spaces are not permitted in paths");
    }
    IHistogram1D* hist = histogramFactory().createHistogram1D(path, title, binedges);
    MSG_TRACE("Made histogram " << hname <<  " for " << name());
    hist->setXTitle(xtitle);
    hist->setYTitle(ytitle);
    histos.push_back(hist);
    return hist;
  }

  IHistogram2D*
  Analysis::bookHistogram2D(const string& hname,
                            size_t nxbins, double xlower, double xupper,
                            size_t nybins, double ylower, double yupper,
                            const string& title, const string& xtitle,
                            const string& ytitle, const string& ztitle) {
    _makeHistoDir();
    const string path = histoPath(hname);
    if (path.find(" ") != string::npos) {
      throw Error("Histogram path '" + path + "' is invalid: spaces are not permitted in paths");
    }
    IHistogram2D* hist =
      histogramFactory().createHistogram2D(path, title, nxbins, xlower, xupper,
                                           nybins, ylower, yupper);
    MSG_TRACE("Made 2D histogram " << hname <<  " for " << name());
    hist->setXTitle(xtitle);
    hist->setYTitle(ytitle);
    hist->setZTitle(ztitle);

    return hist;
  }


  IHistogram2D*
  Analysis::bookHistogram2D(const string& hname,
                            const vector<double>& xbinedges,
                            const vector<double>& ybinedges,
                            const string& title, const string& xtitle,
                            const string& ytitle, const string& ztitle) {
    _makeHistoDir();
    const string path = histoPath(hname);
    if (path.find(" ") != string::npos) {
      throw Error("Histogram path '" + path + "' is invalid: spaces are not permitted in paths");
    }
    IHistogram2D* hist =
      histogramFactory().createHistogram2D(path, title, xbinedges, ybinedges);
    MSG_TRACE("Made 2D histogram " << hname <<  " for " << name());
    hist->setXTitle(xtitle);
    hist->setYTitle(ytitle);
    hist->setZTitle(ztitle);
    return hist;
  }


  /////////////////


  IProfile1D* Analysis::bookProfile1D(size_t datasetId, size_t xAxisId,
                                      size_t yAxisId, const string& title,
                                      const string& xtitle, const string& ytitle) {
    const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
    return bookProfile1D(axisCode, title, xtitle, ytitle);
  }


  IProfile1D* Analysis::bookProfile1D(const string& hname, const string& title,
                                      const string& xtitle, const string& ytitle)
  {
    // Get the bin edges (only read the AIDA file once)
    const BinEdges edges = binEdges(hname);
    _makeHistoDir();
    const string path = histoPath(hname);
    if (path.find(" ") != string::npos) {
      throw Error("Histogram path '" + path + "' is invalid: spaces are not permitted in paths");
    }
    IProfile1D* prof = histogramFactory().createProfile1D(path, title, edges);
    MSG_TRACE("Made profile histogram " << hname <<  " for " << name());
    prof->setXTitle(xtitle);
    prof->setYTitle(ytitle);
    return prof;
  }


  IProfile1D* Analysis::bookProfile1D(const string& hname,
                                      size_t nbins, double lower, double upper,
                                      const string& title,
                                      const string& xtitle, const string& ytitle) {
    _makeHistoDir();
    const string path = histoPath(hname);
    if (path.find(" ") != string::npos) {
      throw Error("Histogram path '" + path + "' is invalid: spaces are not permitted in paths");
    }
    IProfile1D* prof = histogramFactory().createProfile1D(path, title, nbins, lower, upper);
    MSG_TRACE("Made profile histogram " << hname <<  " for " << name());
    prof->setXTitle(xtitle);
    prof->setYTitle(ytitle);
    return prof;
  }


  IProfile1D* Analysis::bookProfile1D(const string& hname,
                                      const vector<double>& binedges,
                                      const string& title,
                                      const string& xtitle, const string& ytitle) {
    _makeHistoDir();
    const string path = histoPath(hname);
    if (path.find(" ") != string::npos) {
      throw Error("Histogram path '" + path + "' is invalid: spaces are not permitted in paths");
    }
    IProfile1D* prof = histogramFactory().createProfile1D(path, title, binedges);
    MSG_TRACE("Made profile histogram " << hname <<  " for " << name());
    prof->setXTitle(xtitle);
    prof->setYTitle(ytitle);
    return prof;
  }


  ///////////////////



  IDataPointSet* Analysis::bookDataPointSet(const string& hname, const string& title,
                                            const string& xtitle, const string& ytitle) {
    _makeHistoDir();
    const string path = histoPath(hname);
    if (path.find(" ") != string::npos) {
      throw Error("Histogram path '" + path + "' is invalid: spaces are not permitted in paths");
    }
    IDataPointSet* dps = datapointsetFactory().create(path, title, 2);
    MSG_TRACE("Made data point set " << hname <<  " for " << name());
    dps->setXTitle(xtitle);
    dps->setYTitle(ytitle);
    return dps;
  }


  IDataPointSet* Analysis::bookDataPointSet(const string& hname,
                                            size_t npts, double lower, double upper,
                                            const string& title,
                                            const string& xtitle, const string& ytitle) {
    IDataPointSet* dps = bookDataPointSet(hname, title, xtitle, ytitle);
    for (size_t pt = 0; pt < npts; ++pt) {
      const double binwidth = (upper-lower)/npts;
      const double bincentre = lower + (pt + 0.5) * binwidth;
      dps->addPoint();
      IMeasurement* meas = dps->point(pt)->coordinate(0);
      meas->setValue(bincentre);
      meas->setErrorPlus(binwidth/2.0);
      meas->setErrorMinus(binwidth/2.0);
    }
    return dps;
  }


  IDataPointSet* Analysis::bookDataPointSet(size_t datasetId, size_t xAxisId,
                                            size_t yAxisId, const string& title,
                                            const string& xtitle, const string& ytitle) {
    // Get the bin edges (only read the AIDA file once)
    _cacheXAxisData();
    // Build the axis code
    const string axisCode = makeAxisCode(datasetId, xAxisId, yAxisId);
    //const map<string, vector<DPSXPoint> > xpoints = getDPSXValsErrs(papername);
    MSG_TRACE("Using DPS x-positions for " << name() << ":" << axisCode);
    IDataPointSet* dps = bookDataPointSet(axisCode, title, xtitle, ytitle);
    const vector<DPSXPoint> xpts = _dpsData.find(axisCode)->second;
    for (size_t pt = 0; pt < xpts.size(); ++pt) {
      dps->addPoint();
      IMeasurement* meas = dps->point(pt)->coordinate(0);
      meas->setValue(xpts[pt].val);
      meas->setErrorPlus(xpts[pt].errplus);
      meas->setErrorMinus(xpts[pt].errminus);
    }
    MSG_TRACE("Made DPS " << axisCode <<  " for " << name());
    return dps;
  }


  ////////////////////


  void Analysis::_makeHistoDir() {
    if (!_madeHistoDir) {
      if (! name().empty()) {
        // vector<string> dirs;
        // split(dirs, histoDir(), "/");
        // string pathpart;
        // foreach (const string& d, dirs) {
        //tree().mkdir();
        //}
        tree().mkdirs(histoDir());
      }
      _madeHistoDir = true;
    }
  }


  void Analysis::normalize(AIDA::IHistogram1D*& histo, double norm, bool includeoverflows) {
    if (!histo) {
      MSG_ERROR("Failed to normalize histo=NULL in analysis "
                << name() << " (norm=" << norm << ")");
      return;
    }
    const string hpath = tree().findPath(dynamic_cast<const AIDA::IManagedObject&>(*histo));
    MSG_TRACE("Normalizing histo " << hpath << " to " << norm);

    // Get integral
    double oldintg = 0.0;
    int nBins = histo->axis().bins();
    for (int iBin = 0; iBin != nBins; ++iBin) {
      // Leaving out factor of binWidth because AIDA's "height" already includes a width factor.
      oldintg += histo->binHeight(iBin); // * histo->axis().binWidth(iBin);
    }
    if (includeoverflows) {
      // Include overflow bins in the integral
      oldintg += histo->binHeight(AIDA::IAxis::UNDERFLOW_BIN);
      oldintg += histo->binHeight(AIDA::IAxis::OVERFLOW_BIN);
    }
    // Sanity check
    if (oldintg == 0.0) {
      MSG_WARNING("Histo " << hpath << " has null integral during normalization");
      return;
    }

    // Scale by the normalisation factor.
    scale(histo, norm/oldintg);
  }

  void Analysis::finishEvent() {
    if(!histos.empty())
      for(vector<AIDA::IHistogram1D *>::iterator it= histos.begin();it!= histos.end();it++)
        if((*it)) (*it)->finishEvent();
  }
  
  void Analysis::scale(AIDA::IHistogram1D*& histo, double scale) {
    if (!histo) {
      MSG_ERROR("Failed to scale histo=NULL in analysis "
                << name() << " (scale=" << scale << ")");
      return;
    }
    const string hpath = tree().findPath(dynamic_cast<const AIDA::IManagedObject&>(*histo));
    MSG_TRACE("Scaling histo " << hpath);

    vector<double> x, y, ex, ey;
    for (size_t i = 0, N = histo->axis().bins(); i < N; ++i) {
      x.push_back(0.5 * (histo->axis().binLowerEdge(i) + histo->axis().binUpperEdge(i)));
      ex.push_back(histo->axis().binWidth(i)*0.5);

      // "Bin height" is a misnomer in the AIDA spec: width is neglected.
      // We'd like to do this: y.push_back(histo->binHeight(i) * scale);
      y.push_back(histo->binHeight(i)*scale/histo->axis().binWidth(i));

      // "Bin error" is a misnomer in the AIDA spec: width is neglected.
      // We'd like to do this: ey.push_back(histo->binError(i) * scale);
      ey.push_back(histo->binError(i)*scale/histo->axis().binWidth(i));
    }

    string title = histo->title();
    string xtitle = histo->xtitle();
    string ytitle = histo->ytitle();

    tree().mkdir("/tmpnormalize");
    tree().mv(hpath, "/tmpnormalize");

    if (hpath.find(" ") != string::npos) {
      throw Error("Histogram path '" + hpath + "' is invalid: spaces are not permitted in paths");
    }
    AIDA::IDataPointSet* dps = datapointsetFactory().createXY(hpath, title, x, y, ex, ey);
    dps->setXTitle(xtitle);
    dps->setYTitle(ytitle);

    tree().rm(tree().findPath(dynamic_cast<AIDA::IManagedObject&>(*histo)));
    tree().rmdir("/tmpnormalize");

    // Set histo pointer to null - it can no longer be used.
    histo = 0;
  }


  void Analysis::normalize(AIDA::IHistogram2D*& histo, double norm) {
    if (!histo) {
      MSG_ERROR("Failed to normalize histo=NULL in analysis "
                << name() << " (norm=" << norm << ")");
      return;
    }
    const string hpath = tree().findPath(dynamic_cast<const AIDA::IManagedObject&>(*histo));
    MSG_TRACE("Normalizing histo " << hpath << " to " << norm);

    double oldintg = 0.0;
    int nxBins = histo->xAxis().bins();
    int nyBins = histo->yAxis().bins();
    for (int ixBin = 0; ixBin != nxBins; ++ixBin)
      for (int iyBin = 0; iyBin != nyBins; ++iyBin) {
      // Leaving out factor of binWidth because AIDA's "height"
      // already includes a width factor.
        oldintg += histo->binHeight(ixBin, iyBin); // * histo->axis().binWidth(iBin);
    }
    if (oldintg == 0.0) {
      MSG_WARNING("Histo " << hpath << " has null integral during normalization");
      return;
    }

    // Scale by the normalisation factor.
    scale(histo, norm/oldintg);
  }


  void Analysis::scale(AIDA::IHistogram2D*& histo, double scale) {
    if (!histo) {
      MSG_ERROR("Failed to scale histo=NULL in analysis "
                << name() << " (scale=" << scale << ")");
      return;
    }
    const string hpath =
      tree().findPath(dynamic_cast<const AIDA::IManagedObject&>(*histo));
    MSG_TRACE("Scaling histo " << hpath);

    vector<double> x, y, z, ex, ey, ez;
    for (size_t ix = 0, Nx = histo->xAxis().bins(); ix < Nx; ++ix)
      for (size_t iy = 0, Ny = histo->yAxis().bins(); iy < Ny; ++iy) {
        x.push_back(0.5 * (histo->xAxis().binLowerEdge(ix) +
                           histo->xAxis().binUpperEdge(ix)));
        ex.push_back(histo->xAxis().binWidth(ix)*0.5);
        y.push_back(0.5 * (histo->yAxis().binLowerEdge(iy) +
                           histo->yAxis().binUpperEdge(iy)));
        ey.push_back(histo->yAxis().binWidth(iy)*0.5);

        // "Bin height" is a misnomer in the AIDA spec: width is neglected.
        // We'd like to do this: y.push_back(histo->binHeight(i) * scale);
        z.push_back(histo->binHeight(ix, iy)*scale/
                    (histo->xAxis().binWidth(ix)*histo->yAxis().binWidth(iy)));
        // "Bin error" is a misnomer in the AIDA spec: width is neglected.
        // We'd like to do this: ey.push_back(histo->binError(i) * scale);
        ez.push_back(histo->binError(ix, iy)*scale/
                     (histo->xAxis().binWidth(ix)*histo->yAxis().binWidth(iy)));
    }

    string title = histo->title();
    string xtitle = histo->xtitle();
    string ytitle = histo->ytitle();
    string ztitle = histo->ztitle();

    tree().mkdir("/tmpnormalize");
    tree().mv(hpath, "/tmpnormalize");

    if (hpath.find(" ") != string::npos) {
      throw Error("Histogram path '" + hpath + "' is invalid: spaces are not permitted in paths");
    }
    AIDA::IDataPointSet* dps =
      datapointsetFactory().createXYZ(hpath, title, x, y, z, ex, ey, ez);
    dps->setXTitle(xtitle);
    dps->setYTitle(ytitle);
    dps->setZTitle(ztitle);

    tree().rm(tree().findPath(dynamic_cast<AIDA::IManagedObject&>(*histo)));
    tree().rmdir("/tmpnormalize");

    // Set histo pointer to null - it can no longer be used.
    histo = 0;
  }


}
