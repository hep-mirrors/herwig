// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

/// @todo Use inline PID functions instead

namespace Rivet {


  /// @brief SLD b-fragmentation measurement
  /// @author Peter Richardson
  class Clustermass : public Analysis {
  public:

    /// Constructor
    Clustermass()
      : Analysis("Clustermass")
    {
    }


    /// @name Analysis methods
    //@{

    /// Book projections and histograms
    void init() {

      const FinalState fs; 
      addProjection(fs, "FS");
      FastJets durhamJets = FastJets(fs, FastJets::DURHAM, 0.7);
      durhamJets.useInvisibles(true);
      addProjection(durhamJets, "DurhamJets");

      _h_R_Durham[0] = bookHisto1D("R2", logspace(50, 0.001, 1.));
      _h_y_Durham[0] = bookHisto1D("y23", logspace(50, 0.001, 1.));
      _h_R_Durham[1] = bookHisto1D("R3", logspace(50, 0.001, 1.));
      _h_y_Durham[1] = bookHisto1D("y34", logspace(50, 0.001, 1.));
      _h_R_Durham[2] = bookHisto1D("R4", logspace(50, 0.001, 1.));
      _h_y_Durham[2] = bookHisto1D("y45", logspace(50, 0.001, 1.));
      _h_R_Durham[3] = bookHisto1D("R5", logspace(50, 0.001, 1.));
      _h_y_Durham[3] = bookHisto1D("y56", logspace(50, 0.001, 1.));
      _h_R_Durham[4] = bookHisto1D("R6", logspace(50, 0.001, 1.));
      _h_y_Durham[4] = bookHisto1D("y67", logspace(50, 0.001, 1.));
      _h_R_Durham[5] = bookHisto1D("R7", logspace(50, 0.001, 1.));



      _histfirst =bookHisto1D("firstclustermass", logspace(100, 0.1, 100.));
      _histlast =bookHisto1D("lastclustermass", logspace(100, 0.1, 100.)); 
    }


    void analyze(const Event& e) {
      // Get event weight for histo filling
      const double weight = e.weight();
      const FastJets& durjet = applyProjection<FastJets>(e, "DurhamJets");
      if (durjet.clusterSeq()) {
        const double y_23 = durjet.clusterSeq()->exclusive_ymerge_max(2);
        const double y_34 = durjet.clusterSeq()->exclusive_ymerge_max(3);
        const double y_45 = durjet.clusterSeq()->exclusive_ymerge_max(4);
        const double y_56 = durjet.clusterSeq()->exclusive_ymerge_max(5);
        const double y_67 = durjet.clusterSeq()->exclusive_ymerge_max(6);
        _h_y_Durham[0]->fill(y_23, weight);
        _h_y_Durham[1]->fill(y_34, weight);
        _h_y_Durham[2]->fill(y_45, weight);
        _h_y_Durham[3]->fill(y_56, weight);
        _h_y_Durham[4]->fill(y_67, weight);

        for (size_t i = 0; i < _h_R_Durham[0]->numBins(); ++i) {
          double ycut = _h_R_Durham[0]->bin(i).xMid();
          double width = _h_R_Durham[0]->bin(i).xWidth();
          if (y_23 < ycut) {
            _h_R_Durham[0]->fillBin(i, weight*width);
          }
        }
        for (size_t i = 0; i < _h_R_Durham[1]->numBins(); ++i) {
          double ycut = _h_R_Durham[1]->bin(i).xMid();
          double width = _h_R_Durham[1]->bin(i).xWidth();
          if (y_34 < ycut && y_23 > ycut) {
            _h_R_Durham[1]->fillBin(i, weight*width);
          }
        }
        for (size_t i = 0; i < _h_R_Durham[2]->numBins(); ++i) {
          double ycut = _h_R_Durham[2]->bin(i).xMid();
          double width = _h_R_Durham[2]->bin(i).xWidth();
          if (y_45 < ycut && y_34 > ycut) {
            _h_R_Durham[2]->fillBin(i, weight*width);
          }
        }
        for (size_t i = 0; i < _h_R_Durham[3]->numBins(); ++i) {
          double ycut = _h_R_Durham[3]->bin(i).xMid();
          double width = _h_R_Durham[3]->bin(i).xWidth();
          if (y_56 < ycut && y_45 > ycut) {
            _h_R_Durham[3]->fillBin(i, weight*width);
          }
        }
        for (size_t i = 0; i < _h_R_Durham[4]->numBins(); ++i) {
          double ycut = _h_R_Durham[4]->bin(i).xMid();
          double width = _h_R_Durham[4]->bin(i).xWidth();
          if (y_67 < ycut&& y_56 > ycut ) {
            _h_R_Durham[4]->fillBin(i, weight*width);
          }
        }
        for (size_t i = 0; i < _h_R_Durham[5]->numBins(); ++i) {
          double ycut = _h_R_Durham[5]->bin(i).xMid();
          double width = _h_R_Durham[5]->bin(i).xWidth();
          if (y_67 > ycut) {
            _h_R_Durham[5]->fillBin(i, weight*width);
          }
        }
      }
    


      foreach (const GenParticle* p, particles(e.genEvent())) {
        const GenVertex* dv = p->end_vertex();
        if (p->pdg_id()==81) {
          if (dv) {
            bool is_last = true;
            for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin() ;
                 pp != dv->particles_out_const_end() ; ++pp) {
              if ((*pp)->pdg_id()==81) {
                is_last = false;
		break;
              }
            }
            if (is_last) {
              _histlast->fill(p->momentum().m(), weight);
            }
          }
	  bool isFirst = true;
          foreach (GenParticle* pp, Rivet::particles(p->production_vertex(), HepMC::parents)) {
            if (abs(pp->pdg_id())==81) {
              
              isFirst = false;
              break;
            }
          }
          if(isFirst)
	     _histfirst->fill(p->momentum().m(), weight);
	  

        }
      }
    }


    // Finalize
    void finalize() {
     const double s = crossSection()/picobarn/sumOfWeights();
      scale(_histlast, s); 
      scale(_histfirst, s);
      for (size_t n = 0; n < 5; ++n) scale(_h_y_Durham[n],s);
      for (size_t n = 0; n < 6; ++n) scale(_h_R_Durham[n], s);
    }


  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the
    /// inclusive single particle distributions' normalisations.

    Histo1DPtr _histlast,_histfirst;
    Histo1DPtr _h_R_Durham[6];
    Histo1DPtr _h_y_Durham[5];

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(Clustermass);

}
