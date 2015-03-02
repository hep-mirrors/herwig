// -*- C++ -*-
#ifndef RIVET_Analysis_HH
#define RIVET_Analysis_HH

#include "Rivet/Rivet.hh"
#include "Rivet/Analysis.fhh"
#include "Rivet/AnalysisInfo.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/ProjectionApplier.hh"
#include "Rivet/ProjectionHandler.hh"
#include "Rivet/Constraints.hh"
#include "Rivet/AnalysisHandler.fhh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/Logging.fhh"
#include "Rivet/RivetAIDA.fhh"


/// @def vetoEvent
/// Preprocessor define for vetoing events, including the log message and return.
#define vetoEvent \
  do { MSG_DEBUG("Vetoing event on line " << __LINE__ << " of " << __FILE__); return; } while(0)

/// @def DECLARE_RIVET_PLUGIN
/// Preprocessor define to prettify the global-object plugin hook mechanism.
#define DECLARE_RIVET_PLUGIN(clsname) Rivet::AnalysisBuilder<clsname> plugin_ ## clsname



namespace Rivet {


  /// @brief This is the base class of all analysis classes in Rivet.
  ///
  /// There are
  /// three virtual functions which should be implemented in base classes:
  ///
  /// void init() is called by Rivet before a run is started. Here the
  /// analysis class should book necessary histograms. The needed
  /// projections should probably rather be constructed in the
  /// constructor.
  ///
  /// void analyze(const Event&) is called once for each event. Here the
  /// analysis class should apply the necessary Projections and fill the
  /// histograms.
  ///
  /// void finalize() is called after a run is finished. Here the analysis
  /// class should do whatever manipulations are necessary on the
  /// histograms. Writing the histograms to a file is, however, done by
  /// the Rivet class.
  class Analysis : public ProjectionApplier {

    /// The AnalysisHandler is a friend.
    friend class AnalysisHandler;


  public:

    /// @name Standard constructors and destructors.
    //@{

    // /// The default constructor.
    // Analysis();

    /// Constructor
    Analysis(const std::string& name);

    /// The destructor.
    virtual ~Analysis() {}

    //@}


  public:

    /// @name Main analysis methods
    //@{

    /// Initialize this analysis object. A concrete class should here
    /// book all necessary histograms. An overridden function must make
    /// sure it first calls the base class function.
    virtual void init() { }

    /// Analyze one event. A concrete class should here apply the
    /// necessary projections on the \a event and fill the relevant
    /// histograms. An overridden function must make sure it first calls
    /// the base class function.
    virtual void analyze(const Event& event) = 0;

    /// Finalize this analysis object. A concrete class should here make
    /// all necessary operations on the histograms. Writing the
    /// histograms to a file is, however, done by the Rivet class. An
    /// overridden function must make sure it first calls the base class
    /// function.
    virtual void finalize() { }

    //@}


  public:

    /// @name Metadata
    /// Metadata is used for querying from the command line and also for
    /// building web pages and the analysis pages in the Rivet manual.
    //@{

    /// Get the actual AnalysisInfo object in which all this metadata is stored.
    const AnalysisInfo& info() const {
      assert(_info.get() != 0 && "No AnalysisInfo object :O");
      return *_info;
    }

    /// @brief Get the name of the analysis.
    ///
    /// By default this is computed by combining the results of the experiment,
    /// year and Spires ID metadata methods and you should only override it if
    /// there's a good reason why those won't work.
    virtual std::string name() const {
      return (info().name().empty()) ? _defaultname : info().name();
    }

    /// Get the Inspire ID code for this analysis.
    virtual std::string inspireId() const {
      return info().inspireId();
    }

    /// Get the SPIRES ID code for this analysis (~deprecated).
    virtual std::string spiresId() const {
      return info().spiresId();
    }

    /// @brief Names & emails of paper/analysis authors.
    ///
    /// Names and email of authors in 'NAME \<EMAIL\>' format. The first
    /// name in the list should be the primary contact person.
    virtual std::vector<std::string> authors() const {
      return info().authors();
    }

    /// @brief Get a short description of the analysis.
    ///
    /// Short (one sentence) description used as an index entry.
    /// Use @a description() to provide full descriptive paragraphs
    /// of analysis details.
    virtual std::string summary() const {
      return info().summary();
    }

    /// @brief Get a full description of the analysis.
    ///
    /// Full textual description of this analysis, what it is useful for,
    /// what experimental techniques are applied, etc. Should be treated
    /// as a chunk of restructuredText (http://docutils.sourceforge.net/rst.html),
    /// with equations to be rendered as LaTeX with amsmath operators.
    virtual std::string description() const {
      return info().description();
    }

    /// @brief Information about the events needed as input for this analysis.
    ///
    /// Event types, energies, kinematic cuts, particles to be considered
    /// stable, etc. etc. Should be treated as a restructuredText bullet list
    /// (http://docutils.sourceforge.net/rst.html)
    virtual std::string runInfo() const {
      return info().runInfo();
    }

    /// Experiment which performed and published this analysis.
    virtual std::string experiment() const {
      return info().experiment();
    }

    /// Collider on which the experiment ran.
    virtual std::string collider() const {
      return info().collider();
    }

    /// When the original experimental analysis was published.
    virtual std::string year() const {
      return info().year();
    }

    /// Journal, and preprint references.
    virtual std::vector<std::string> references() const {
      return info().references();
    }

    /// BibTeX citation key for this article.
    virtual std::string bibKey() const {
      return info().bibKey();
    }

    /// BibTeX citation entry for this article.
    virtual std::string bibTeX() const {
      return info().bibTeX();
    }

    /// Whether this analysis is trusted (in any way!)
    virtual std::string status() const {
      return (info().status().empty()) ? "UNVALIDATED" : info().status();
    }

    /// Any work to be done on this analysis.
    virtual std::vector<std::string> todos() const {
      return info().todos();
    }


    /// Return the allowed pairs of incoming beams required by this analysis.
    virtual const std::vector<PdgIdPair>& requiredBeams() const {
      return info().beams();
    }
    /// Declare the allowed pairs of incoming beams required by this analysis.
    virtual Analysis& setRequiredBeams(const std::vector<PdgIdPair>& requiredBeams) {
      info().setBeams(requiredBeams);
      return *this;
    }


    /// Sets of valid beam energy pairs, in GeV
    virtual const std::vector<std::pair<double, double> >& requiredEnergies() const {
      return info().energies();
    }
    /// Declare the list of valid beam energy pairs, in GeV
    virtual Analysis& setRequiredEnergies(const std::vector<std::pair<double, double> >& requiredEnergies) {
      info().setEnergies(requiredEnergies);
      return *this;
    }


    /// Return true if this analysis needs to know the process cross-section.
    bool needsCrossSection() const {
      return info().needsCrossSection();
    }
    /// Declare whether this analysis needs to know the process cross-section from the generator.
    Analysis& setNeedsCrossSection(bool needed=true) {
      info().setNeedsCrossSection(needed);
      return *this;
    }

    //@}


    /// @name Internal metadata modifiying methods
    //@{

    /// Get the actual AnalysisInfo object in which all this metadata is stored (non-const).
    AnalysisInfo& info() {
      assert(_info.get() != 0 && "No AnalysisInfo object :O");
      return *_info;
    }

    /// Set the required beams
    /// @deprecated To be removed in 2.0.0. Use .info file and AnalysisInfo class instead
    virtual Analysis& setBeams(PdgId beam1, PdgId beam2) {
      /// @todo Print out a warning to use setRequiredBeams() instead (and really to use .info files)
      return setRequiredBeams(std::vector<PdgIdPair>(1, make_pair(beam1, beam2)));
    }

    //@}


    /// @name Run conditions
    //@{

    /// Incoming beams for this run
    const ParticlePair& beams() const;

    /// Incoming beam IDs for this run
    const PdgIdPair beamIds() const;

    /// Centre of mass energy for this run
    double sqrtS() const;

    //@}


    /// @name Analysis / beam compatibility testing
    //@{

    /// Check if analysis is compatible with the provided beam particle IDs and energies
    bool isCompatible(const ParticlePair& beams) const;

    /// Check if analysis is compatible with the provided beam particle IDs and energies
    bool isCompatible(PdgId beam1, PdgId beam2, double e1, double e2) const;

    /// Check if analysis is compatible with the provided beam particle IDs and energies
    bool isCompatible(const PdgIdPair& beams, const std::pair<double,double>& energies) const;

    //@}

  public:

    /// Access the controlling AnalysisHandler object.
    AnalysisHandler& handler() const { return *_analysishandler; }

    /// Normalize the given histogram, @a histo. After this call the
    /// histogram will have been transformed to a DataPointSet with the
    /// same name and path. It has the same effect as
    /// @c scale(histo, norm/sumOfWeights).
    /// @param histo The histogram to be normalised.
    /// @param norm The new area of the histogram.
    /// @warning The old histogram will be deleted, and its pointer set to zero.
    void normalize(AIDA::IHistogram1D*& histo, double norm=1.0, bool includeoverflows=true);

    /// Multiplicatively scale the given histogram, @a histo. After this call the
    /// histogram will have been transformed to a DataPointSet with the same name and path.
    /// @param histo The histogram to be scaled.
    /// @param scale The factor used to multiply the histogram bin heights.
    /// @warning The old histogram will be deleted, and its pointer set to zero.
    void scale(AIDA::IHistogram1D*& histo, double scale);

    /// Normalize the given histogram, @a histo. After this call the
    /// histogram will have been transformed to a DataPointSet with the
    /// same name and path. It has the same effect as
    /// @c scale(histo, norm/sumOfWeights).
    /// @param histo The histogram to be normalised.
    /// @param norm The new area of the histogram.
    /// @warning The old histogram will be deleted, and its pointer set to zero.
    void normalize(AIDA::IHistogram2D*& histo, double norm=1.0);

    /// Multiplicatively scale the given histogram, @a histo. After this call the
    /// histogram will have been transformed to a DataPointSet with the same name and path.
    /// @param histo The histogram to be scaled.
    /// @param scale The factor used to multiply the histogram bin heights.
    /// @warning The old histogram will be deleted, and its pointer set to zero.
    void scale(AIDA::IHistogram2D*& histo, double scale);

    /// Set the cross section from the generator
    Analysis& setCrossSection(double xs);


  protected:

    /// Get a Log object based on the name() property of the calling analysis object.
    Log& getLog() const;

    /// Get the process cross-section in pb. Throws if this hasn't been set.
    double crossSection() const;

    /// Get the process cross-section per generated event in pb. Throws if this
    /// hasn't been set.
    double crossSectionPerEvent() const;

    /// Get the number of events seen (via the analysis handler). Use in the
    /// finalize phase only.
    size_t numEvents() const;

    /// Get the sum of event weights seen (via the analysis handler). Use in the
    /// finalize phase only.
    double sumOfWeights() const;
    
    void finishEvent();


  protected:

    /// @name AIDA analysis infrastructure.
    //@{
    /// Access the AIDA analysis factory of the controlling AnalysisHandler object.
    AIDA::IAnalysisFactory& analysisFactory();

    /// Access the AIDA tree of the controlling AnalysisHandler object.
    AIDA::ITree& tree();

    /// Access the AIDA histogram factory of the controlling AnalysisHandler object.
    AIDA::IHistogramFactory& histogramFactory();

    /// Access the AIDA histogram factory of the controlling AnalysisHandler object.
    AIDA::IDataPointSetFactory& datapointsetFactory();

    /// Get the canonical histogram "directory" path for this analysis.
    const std::string histoDir() const;

    /// Get the canonical histogram path for the named histogram in this analysis.
    const std::string histoPath(const std::string& hname) const;

    /// Get the canonical histogram path for the numbered histogram in this analysis.
    const std::string histoPath(size_t datasetId, size_t xAxisId, size_t yAxisId) const;

    /// Get the internal histogram name for given d, x and y (cf. HepData)
    const std::string makeAxisCode(size_t datasetId, size_t xAxisId, size_t yAxisId) const;

    //@}


    /// @name Internal histogram booking (for use by Analysis sub-classes).
    //@{

    /// Get bin edges for a named histo (using ref AIDA caching)
    const BinEdges& binEdges(const std::string& hname) const;

    /// Get bin edges for a numbered histo (using ref AIDA caching)
    const BinEdges& binEdges(size_t datasetId, size_t xAxisId, size_t yAxisId) const;

    /// @brief Get bin edges with logarithmic widths
    ///
    /// @deprecated Prefer logspace. This will disappear in Rivet 2.0.0
    BinEdges logBinEdges(size_t nbins, double lower, double upper);

    /// Book a 1D histogram with @a nbins uniformly distributed across the range @a lower - @a upper .
    /// (NB. this returns a pointer rather than a reference since it will
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly
    /// get the pointer from a reference before they can use it!)
    AIDA::IHistogram1D* bookHistogram1D(const std::string& name,
                                        size_t nbins, double lower, double upper,
                                        const std::string& title="",
                                        const std::string& xtitle="", const std::string& ytitle="");

    /// Book a 1D histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    /// (NB. this returns a pointer rather than a reference since it will
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly
    /// get the pointer from a reference before they can use it!)
    AIDA::IHistogram1D* bookHistogram1D(const std::string& name,
                                        const std::vector<double>& binedges, const std::string& title="",
                                        const std::string& xtitle="", const std::string& ytitle="");

    /// Book a 2D histogram with @a nxbins and @a nybins uniformly
    /// distributed across the ranges @a xlower - @a xupper and @a
    /// ylower - @a yupper respectively along the x- and y-axis.
    /// (NB. this returns a pointer rather than a reference since it
    /// will have to be stored in the analysis class - there's no
    /// point in forcing users to explicitly get the pointer from a
    /// reference before they can use it!)
    AIDA::IHistogram2D*
    bookHistogram2D(const std::string& name,
		    size_t nxbins, double xlower, double xupper,
		    size_t nybins, double ylower, double yupper,
		    const std::string& title="", const std::string& xtitle="",
		    const std::string& ytitle="", const std::string& ztitle="");

    /// Book a 2D histogram with non-uniform bins defined by the
    /// vectorx of bin edges @a xbinedges and @a ybinedges.
    /// (NB. this returns a pointer rather than a reference since it
    /// will have to be stored in the analysis class - there's no
    /// point in forcing users to explicitly get the pointer from a
    /// reference before they can use it!)
    AIDA::IHistogram2D*
    bookHistogram2D(const std::string& name,
		    const std::vector<double>& xbinedges,
		    const std::vector<double>& ybinedges,
		    const std::string& title="", const std::string& xtitle="",
		    const std::string& ytitle="", const std::string& ztitle="");

    /// Book a 1D histogram based on the name in the corresponding AIDA
    /// file. The binnings will be obtained by reading the bundled AIDA data
    /// record file with the same filename as the analysis' name() property.
    AIDA::IHistogram1D* bookHistogram1D(const std::string& name, const std::string& title="",
                                        const std::string& xtitle="", const std::string& ytitle="");

    /// Book a 1D histogram based on the paper, dataset and x/y-axis IDs in the corresponding
    /// HepData record. The binnings will be obtained by reading the bundled AIDA data record file
    /// of the same filename as the analysis' name() property.
    AIDA::IHistogram1D* bookHistogram1D(size_t datasetId, size_t xAxisId, size_t yAxisId,
                                        const std::string& title="",
                                        const std::string& xtitle="", const std::string& ytitle="");

    //@}


    /// @name Internal profile histogram booking (for use by Analysis sub-classes).
    //@{

    /// Book a 1D profile histogram with @a nbins uniformly distributed across the range @a lower - @a upper .
    /// (NB. this returns a pointer rather than a reference since it will
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly
    /// get the pointer from a reference before they can use it!)
    AIDA::IProfile1D* bookProfile1D(const std::string& name,
                                    size_t nbins, double lower, double upper,
                                    const std::string& title="",
                                    const std::string& xtitle="", const std::string& ytitle="");

    /// Book a 1D profile histogram with non-uniform bins defined by the vector of bin edges @a binedges .
    /// (NB. this returns a pointer rather than a reference since it will
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly
    /// get the pointer from a reference before they can use it!)
    AIDA::IProfile1D* bookProfile1D(const std::string& name,
                                    const std::vector<double>& binedges,
                                    const std::string& title="",
                                    const std::string& xtitle="", const std::string& ytitle="");

    /// Book a 1D profile histogram based on the name in the corresponding AIDA
    /// file. The binnings will be obtained by reading the bundled AIDA data
    /// record file with the same filename as the analysis' name() property.
    AIDA::IProfile1D* bookProfile1D(const std::string& name, const std::string& title="",
                                    const std::string& xtitle="", const std::string& ytitle="");

    /// Book a 1D profile histogram based on the paper, dataset and x/y-axis IDs in the corresponding
    /// HepData record. The binnings will be obtained by reading the bundled AIDA data record file
    /// of the same filename as the analysis' name() property.
    AIDA::IProfile1D* bookProfile1D(size_t datasetId, size_t xAxisId, size_t yAxisId,
                                    const std::string& title="",
                                    const std::string& xtitle="", const std::string& ytitle="");
    //@}


    /// @name Internal data point set booking (for use by Analysis sub-classes).
    //@{

    /// Book a 2-dimensional data point set.
    /// (NB. this returns a pointer rather than a reference since it will
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly
    /// get the pointer from a reference before they can use it!)
    AIDA::IDataPointSet* bookDataPointSet(const std::string& name, const std::string& title="",
                                          const std::string& xtitle="", const std::string& ytitle="");


    /// Book a 2-dimensional data point set with equally spaced points in a range.
    /// (NB. this returns a pointer rather than a reference since it will
    /// have to be stored in the analysis class - there's no point in forcing users to explicitly
    /// get the pointer from a reference before they can use it!)
    AIDA::IDataPointSet* bookDataPointSet(const std::string& name,
                                          size_t npts, double lower, double upper,
                                          const std::string& title="",
                                          const std::string& xtitle="", const std::string& ytitle="");

    /// Book a 2-dimensional data point set based on the corresponding AIDA data
    /// file. The binnings (x-errors) will be obtained by reading the bundled
    /// AIDA data record file of the same filename as the analysis' name()
    /// property.
    //AIDA::IDataPointSet* bookDataPointSet(const std::string& name, const std::string& title);

    /// Book a 2-dimensional data point set based on the paper, dataset and x/y-axis IDs in the corresponding
    /// HepData record. The binnings (x-errors) will be obtained by reading the bundled AIDA data record file
    /// of the same filename as the analysis' name() property.
    AIDA::IDataPointSet* bookDataPointSet(size_t datasetId, size_t xAxisId, size_t yAxisId,
                                          const std::string& title="",
                                          const std::string& xtitle="", const std::string& ytitle="");

    //@}


  private:

    /// @name Utility functions
    //@{

    /// Make the histogram directory.
    void _makeHistoDir();

    /// Get the bin edges for this paper from the reference AIDA file, and cache them.
    void _cacheBinEdges() const;

    /// Get the x-axis points for this paper from the reference AIDA file, and cache them.
    void _cacheXAxisData() const;

    //@}


  protected:

    /// Name passed to constructor (used to find .info analysis data file, and as a fallback)
    string _defaultname;

    /// Pointer to analysis metadata object
    shared_ptr<AnalysisInfo> _info;


  private:

    /// @name Cross-section variables
    //@{
    double _crossSection;
    bool _gotCrossSection;
    //@}

    /// The controlling AnalysisHandler object.
    AnalysisHandler* _analysishandler;

    /// Flag to indicate whether the histogram directory is present
    mutable bool _madeHistoDir;

    /// Collection of x-axis point data to speed up many autobookings: the
    /// reference data file should only be read once.
    /// @todo Reduce memory occupancy, or clear after initialisation?
    mutable map<string, vector<DPSXPoint> > _dpsData;

    /// Collection of cached bin edges to speed up many autobookings: the
    /// reference data file should only be read once.
    /// @todo Reduce memory occupancy, or clear after initialisation?
    mutable map<string, BinEdges> _histBinEdges;
    
    mutable vector<AIDA::IHistogram1D *> histos;
    
  private:

    /// The assignment operator is private and must never be called.
    /// In fact, it should not even be implemented.
    Analysis& operator=(const Analysis&);

  };


}


// Include definition of analysis plugin system so that analyses automatically see it when including Analysis.hh
#include "Rivet/AnalysisBuilder.hh"


#endif
