// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
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
      book(_h_Bc      , "h_Bc"      , 100, 0.0, 1.0);
      book(_h_BcStar  , "h_BcStar"  , 100, 0.0, 1.0);
      book(_h_Bc0     , "h_Bc0"     , 100, 0.0, 1.0);
      book(_h_Bc2     , "h_Bc2"     , 100, 0.0, 1.0);
      book(_h_Bc3     , "h_Bc3"     , 100, 0.0, 1.0);
      book(_h_Bc1     , "h_Bc1"     , 100, 0.0, 1.0);
      book(_h_Bc1prime, "h_Bc1prime", 100, 0.0, 1.0);
      book(_h_Bc2p    , "h_Bc2p"    , 100, 0.0, 1.0);
      book(_h_Bc2pp   , "h_Bc2pp"   , 100, 0.0, 1.0);
      book(_h_BcStarp , "h_BcStarp" , 100, 0.0, 1.0);
      book(_h_etac    , "h_etac"    , 100, 0.0, 1.0);
      book(_h_etab    , "h_etab"    , 100, 0.0, 1.0);
      book(_h_JPsi    , "h_JPsi"    , 100, 0.0, 1.0);
      book(_h_Upsilon , "h_Upsilon" , 100, 0.0, 1.0);
      book(_h_chic0   , "h_chic0"   , 100, 0.0, 1.0);
      book(_h_chic1   , "h_chic1"   , 100, 0.0, 1.0);
      book(_h_chic2   , "h_chic2"   , 100, 0.0, 1.0);
      book(_h_hc      , "h_hc"      , 100, 0.0, 1.0);
      book(_h_etac2   , "h_etac2"   , 100, 0.0, 1.0);
      book(_h_psi2    , "h_psi2"    , 100, 0.0, 1.0);
      book(_h_psi3    , "h_psi3"    , 100, 0.0, 1.0);
      book(_h_psi3770 , "h_psi3770" , 100, 0.0, 1.0);
      book(_h_cc1     , "h_cc1"     , 100, 0.0, 1.0);
      book(_h_bc0     , "h_bc0"     , 100, 0.0, 1.0);
      book(_h_bc1     , "h_bc1"     , 100, 0.0, 1.0);
      
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for(const Particle & p : apply<UnstableParticles>(event,"UFS").particles()) {
	if(p.abspid()==541) 
	  _h_Bc->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==543) 
	  _h_BcStar->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==10541) 
	  _h_Bc0->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==545) 
	  _h_Bc2->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==547) 
	  _h_Bc3->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==10543) 
	  _h_Bc1->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==20543) 
	  _h_Bc1prime->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==10545) 
	  _h_Bc2p ->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==20545) 
	  _h_Bc2pp->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==30543)
	  _h_BcStarp->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==441) 
	  _h_etac->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==551) 
	  _h_etab->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==443) 
	  _h_JPsi->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==553) 
	  _h_Upsilon->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==10441) 
	  _h_chic0->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==20443) 
	  _h_chic1->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==445) 
	  _h_chic2->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==10443) 
	  _h_hc->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==10445) 
	  _h_etac2->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==20445) 
	  _h_psi2->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==447) 
	  _h_psi3->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==30443) 
	  _h_psi3770->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==4403) 
	  _h_cc1->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==5401) 
	  _h_bc0->fill(2.*p.momentum().E()/sqrtS());
	else if(p.abspid()==5403) 
	  _h_bc1->fill(2.*p.momentum().E()/sqrtS());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_Bc      , 0.5/sumW());
      scale(_h_BcStar  , 0.5/sumW());
      scale(_h_Bc0     , 0.5/sumW());
      scale(_h_Bc2     , 0.5/sumW());
      scale(_h_Bc3     , 0.5/sumW());
      scale(_h_Bc1     , 0.5/sumW());
      scale(_h_Bc1prime, 0.5/sumW());
      scale(_h_BcStarp , 0.5/sumW());
      scale(_h_Bc2p    , 0.5/sumW());
      scale(_h_Bc2pp   , 0.5/sumW());
      
      scale(_h_etac    , 0.5/sumW());
      scale(_h_JPsi    , 0.5/sumW());
      
      scale(_h_etab    , 0.5/sumW());
      scale(_h_Upsilon , 0.5/sumW());
    
      scale(_h_chic0   , 0.5/sumW());
      scale(_h_chic1   , 0.5/sumW());
      scale(_h_chic2   , 0.5/sumW());
      scale(_h_hc      , 0.5/sumW());
      scale(_h_etac2   , 0.5/sumW());
      scale(_h_psi2    , 0.5/sumW());
      scale(_h_psi3770 , 0.5/sumW());
      scale(_h_psi3    , 0.5/sumW());
      scale(_h_cc1     , 0.5/sumW());
      scale(_h_bc0     , 0.5/sumW());
      scale(_h_bc1     , 0.5/sumW());
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_Bc,_h_BcStar,_h_Bc0,_h_Bc2,_h_Bc3,_h_Bc1,_h_Bc1prime,_h_BcStarp,_h_Bc2p,_h_Bc2pp;
    Histo1DPtr _h_etac,_h_JPsi,_h_chic0,_h_chic1,_h_chic2,_h_hc,_h_etac2,_h_psi3770,_h_psi2,_h_psi3;
    Histo1DPtr _h_etab,_h_Upsilon;
    Histo1DPtr _h_cc1,_h_bc0,_h_bc1;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MC_Simple_Onium);

}
