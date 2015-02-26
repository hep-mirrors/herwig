// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Logging.hh"

namespace Rivet {




  /// @brief MC validation analysis for higgs [-> tau tau] + jets events
  class MC_HHJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_HHJETS()
      : MC_JetAnalysis("MC_HHJETS", 4, "Jets")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      IdentifiedFinalState ifs(-10, 10, 0.0*GeV);
      ifs.acceptId(25);
      addProjection(ifs,"IFS");

    
      VetoedFinalState vfs;
      vfs.addVetoPairId(25);
      addProjection(FastJets(vfs, FastJets::ANTIKT, 0.4),
                    "Jets");

      _h_H_mass = bookHisto1D("H_mass", 50, 124.7, 125.3);
      _h_HH_mass = bookHisto1D("HH_mass", 50, 240, 800.0);
      _h_HH_dR = bookHisto1D("HH_dR", 25, 0.5, 7.0);

      _h_H_pT = bookHisto1D("H_pT", 50, 0, 1000.0);
      _h_HH_pT = bookHisto1D("HH_pT", 50, 0, 1000.0);
      _h_H_pT_peak = bookHisto1D("H_pT_peak", 25, 0.0, 25.0);
      _h_H_y = bookHisto1D("H_y", 40, -4.0, 4.0);
      _h_H_phi = bookHisto1D("H_phi", 25, 0.0, TWOPI);
      _h_H_jet1_deta = bookHisto1D("H_jet1_deta", 50, -5.0, 5.0);
      _h_H_jet1_dR = bookHisto1D("H_jet1_dR", 25, 0.5, 7.0);

      MC_JetAnalysis::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {

      const IdentifiedFinalState& ifs
        =applyProjection<IdentifiedFinalState>(e,"IFS");
      ParticleVector allp=ifs.particlesByPt();
      //   cout << "allp.size() " << allp.size() << endl;
      if(allp.size()<1) vetoEvent;

 
      const double weight = e.weight();

      FourMomentum hmom(allp[0].momentum());
      if(allp.size() > 1) { 
	FourMomentum hmom2(allp[1].momentum());
	_h_HH_dR->fill(deltaR(hmom, hmom2), weight);
	_h_HH_pT->fill((hmom+hmom2).pT(), weight);
	_h_HH_mass->fill((hmom+hmom2).mass(), weight);
      }
      _h_H_mass->fill(hmom.mass(),weight);
      _h_H_pT->fill(hmom.pT(),weight);
      _h_H_pT_peak->fill(hmom.pT(),weight);
      _h_H_y->fill(hmom.rapidity(),weight);
      _h_H_phi->fill(hmom.azimuthalAngle(),weight);
   

      // get the jet candidates
      Jets jets;
      foreach (const Jet& jet,
               applyProjection<FastJets>(e, "Jets").jetsByPt(20.0*GeV) ) {
          jets.push_back(jet);
      }
    
      if (jets.size() > 0) {
        _h_H_jet1_deta->fill(hmom.eta()-jets[0].momentum().eta(), weight);
        _h_H_jet1_dR->fill(deltaR(hmom, jets[0].momentum()), weight);
      }

      MC_JetAnalysis::analyze(e);
    }


    /// Finalize
    void finalize() {
      scale(_h_H_mass, crossSection()/sumOfWeights());
      scale(_h_HH_mass, crossSection()/sumOfWeights());
      scale(_h_HH_dR, crossSection()/sumOfWeights());
      scale(_h_H_pT, crossSection()/sumOfWeights());
      scale(_h_HH_pT, crossSection()/sumOfWeights());
      scale(_h_H_pT_peak, crossSection()/sumOfWeights());
      scale(_h_H_y, crossSection()/sumOfWeights());
      scale(_h_H_phi, crossSection()/sumOfWeights());
      scale(_h_H_jet1_deta, crossSection()/sumOfWeights());
      scale(_h_H_jet1_dR, crossSection()/sumOfWeights());

      MC_JetAnalysis::finalize();
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_H_mass;
    Histo1DPtr _h_HH_mass;    
    Histo1DPtr _h_HH_pT;    
    Histo1DPtr _h_HH_dR;
    Histo1DPtr _h_H_pT;
    Histo1DPtr _h_H_pT_peak;
    Histo1DPtr _h_H_y;
    Histo1DPtr _h_H_phi;
    Histo1DPtr _h_H_jet1_deta;
    Histo1DPtr _h_H_jet1_dR;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_HHJETS);

}
