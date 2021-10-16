// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief MC analysis of quarkonium production
  class MC_Hadron_Onium : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_Hadron_Onium);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      vector<long> states ={441,10441,100441,443,10443,20443,100443,30443,445,100445,20445,447,
			    551,10551,100551,110551,200551,210551,553,10553,20553,100553,30553,110553,
			    120553,130553,200553,210553,220553,300553,555,10555,20555,100555,120555,200555,557,100557,
			    541,10541,543,10543,20543,545,10545,20545,547,100541,100543,30543,
			    4403,5503,5401,5403};
      for( long pid : states ) {
	book(_h_pT[pid],"h_pT_"+to_string(pid),100,0.,100.);
	book(_h_y[pid],"h_y_"+to_string(pid),100,-10.,10. );
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles particles = apply<UnstableParticles>(event,"UFS").particles();
      for(auto p : particles) {
	map<long,Histo1DPtr>::const_iterator loc = _h_pT.find(p.abspid());
	if(loc!=_h_pT.end())
	  loc->second->fill(p.momentum().perp());
	loc = _h_y.find(p.abspid());
	if(loc!=_h_y.end())
	  loc->second->fill(p.momentum().rapidity());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(const auto& kv: _h_pT)
	scale(kv.second, crossSection()/picobarn/sumW());
      for(const auto& kv: _h_y )
	scale(kv.second, crossSection()/picobarn/sumW());
    }

    ///@}


    /// @name Histograms
    ///@{
    map<long,Histo1DPtr> _h_pT,_h_y; 
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MC_Hadron_Onium);

}
