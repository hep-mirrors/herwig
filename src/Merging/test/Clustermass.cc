// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"


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

      _histfirst =bookHisto1D("firstclustermass", logspace(100, 0.1, 100.));
      _histlast =bookHisto1D("lastclustermass", logspace(100, 0.1, 100.)); 
    }


    void analyze(const Event& e) {


      // Get event weight for histo filling
      const double weight = e.weight();

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
    }


  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the
    /// inclusive single particle distributions' normalisations.

    Histo1DPtr _histlast,_histfirst;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(Clustermass);

}
