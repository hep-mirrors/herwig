// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief MC analysis for onium spectra in e+e-
  class MC_Simple_Onium : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_Simple_Onium);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      declare(UnstableParticles(),"UFS");
      // Book histograms
      vector<long> states ={441,10441,100441,443,10443,20443,100443,30443,445,100445,20445,447,
			    551,10551,100551,110551,200551,210551,553,10553,20553,100553,30553,110553,
			    120553,130553,200553,210553,220553,300553,555,10555,20555,100555,120555,200555,557,100557,
			    541,10541,543,10543,20543,545,10545,20545,547,100541,100543,30543};
      for( long pid : states ) {
	book(_h[pid],"h_"+to_string(pid),100,0.,1.);
      }
      // book(_h_cc1     , "h_cc1"     , 100, 0.0, 1.0);
      // book(_h_bc0     , "h_bc0"     , 100, 0.0, 1.0);
      // book(_h_bc1     , "h_bc1"     , 100, 0.0, 1.0);
      
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle & p : apply<UnstableParticles>(event,"UFS").particles()) {
	map<long,Histo1DPtr>::const_iterator loc = _h.find(p.abspid());
	if(loc!=_h.end())
	 loc->second->fill(2.*p.momentum().E()/sqrtS());
      // 	else if(p.abspid()==4403) 
      // 	  _h_cc1->fill(2.*p.momentum().E()/sqrtS());
      // 	else if(p.abspid()==5401) 
      // 	  _h_bc0->fill(2.*p.momentum().E()/sqrtS());
      // 	else if(p.abspid()==5403) 
      // 	  _h_bc1->fill(2.*p.momentum().E()/sqrtS());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (const auto & kv : _h) {
	scale(kv.second,0.5/sumW());
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    map<long,Histo1DPtr> _h;
    // Histo1DPtr _h_cc1,_h_bc0,_h_bc1;
  
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MC_Simple_Onium);

}
