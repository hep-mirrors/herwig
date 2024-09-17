// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief MC analysis of e+e- > e+e- X via gamma gamma
  class MC_EEtoEEX : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MC_EEtoEEX);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(FinalState(),"FS");
      declare(UnstableParticles(), "UFS");
      // book the histograms
      vector<double> thresholds={0.,2.*.1349770,2.*.13957061,2.*.497611,2.*.493677,2.*.93827,
				 2.*1.86483,2.*1.86965,2.*5.27963,2.*5.27932};
      _h_m.resize(thresholds.size(),vector<Histo1DPtr>(4));
      for(unsigned int ix=0;ix<_h_m.size();++ix) {
	book(_h_m[ix][0],"h_m_"+toString(ix)+"_0",100,thresholds[ix],20.);
	book(_h_m[ix][1],"h_m_"+toString(ix)+"_1",100,thresholds[ix],thresholds[ix]+2.);
	if(thresholds[ix]<4.6) book(_h_m[ix][2],"h_m_"+toString(ix)+"_2",100,2.6,4.6);
	if(thresholds[ix]<11.)  book(_h_m[ix][3],"h_m_"+toString(ix)+"_3",100,9.,11.);
      }
      _h_onium.resize(2,vector<Histo1DPtr>(3));
      for(unsigned int ix=0;ix<_h_onium.size();++ix) {
	for(unsigned int iy=0;iy<_h_onium[ix].size();++iy)
	  book(_h_onium[ix][iy],"h_onium_"+toString(ix)+"_"+toString(iy),4,0.5,4.5);
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        } else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    bool findScattered(Particle beam, double& q2, Particle & scat) {
      bool found = false;
      scat = beam;
      while (!scat.children().empty()) {
        found = false;
        for (const Particle & p : scat.children()) {
          if (p.pid()==scat.pid()) {
            scat=p;
            found=true;
            break;
          }
        }
        if (!found) break;
      }
      if (!found) return false;
      q2 = -(beam.momentum() - scat.momentum()).mass2();
      return true;
    }

    void fillHistos(const vector<Histo1DPtr> & histos, double value) {
      for(unsigned int ix=0;ix<histos.size();++ix)
	if(histos[ix]) histos[ix]->fill(value);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find scattered leptons and calc Q2
      const Beam& beams = apply<Beam>(event, "Beams");
      double q12 = -1, q22 = -1;
      pair<Particle,Particle> scattered;
      if (!findScattered(beams.beams().first,  q12,scattered.first )) vetoEvent;
      if (!findScattered(beams.beams().second, q22,scattered.second)) vetoEvent;
      // check the final state
      const FinalState & fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      FourMomentum pTotal;
      for (const Particle& p : fs.particles()) {
	if(p.isSame(scattered.first) || p.isSame(scattered.second)) continue;
        nCount[p.pid()] += 1;
        ++ntotal;
	pTotal+=p.momentum();
      }
      // histograms of the mass of the gamma gamma system
      double W = pTotal.mass();
      // all cases
      for(unsigned int iy=0;iy<_h_m[0].size();++iy) _h_m[0][iy]->fill(W);
      // various 2 particle systems
      if(ntotal==2) {
	if(nCount[111]==2)
	  fillHistos(_h_m[1],W);
	else if(nCount[211]==1 && nCount[-211]==1)
	  fillHistos(_h_m[2],W);
	else if(nCount[130]+nCount[310]+nCount[311]+nCount[-311]==2)
	  fillHistos(_h_m[3],W);
	else if(nCount[321]==1 && nCount[-321]==1)
	  fillHistos(_h_m[4],W);
	else if(nCount[2212]==1 && nCount[-2212]==1)
	  fillHistos(_h_m[5],W);
      }
      // heavy meson 2 body states
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      vector<int> heavy={421,411,511,521};
      bool matched = false;
      for(unsigned int ix=0;ix<heavy.size();++ix) {
	int pid = heavy[ix];
	for (const Particle& p1 : ufs.particles(Cuts::pid==pid)) {
	  // find the children
	  map<long,int> nRes = nCount;
	  int ncount = ntotal;
	  findChildren(p1,nRes,ncount);
	  for(const Particle & p2 : ufs.particles(Cuts::pid==-pid)) {
	    map<long,int> nRes2 = nRes;
	    int ncount2 = ncount;
	    findChildren(p2,nRes2,ncount2);
	    if(ncount2!=0) continue;
	    matched=true;
	    for(auto const & val : nRes2) {
	      if(val.second!=0) {
		matched = false;
		break;
	      }
	    }
	    if(matched) {
	      fillHistos(_h_m[6+ix],W);
	    }
	  }
	  if(matched) break;
	}
	if(matched) break;
      }
      // onium resonances
      vector<int> states = {1,10001,5,10005};
      matched = false;
      for(unsigned int iq=4;iq<6;++iq) {
	for(unsigned int n=0;n<3;++n) {
	  for(unsigned int is=0;is<4;++is) {
	    int id   = iq*110+states[is]+n*100000;
	    for (const Particle& p : ufs.particles(Cuts::pid==id)) {
	      if(p.children().empty()) continue;
	      map<long,int> nRes = nCount;
	      int ncount = ntotal;
	      findChildren(p,nRes,ncount);
	      bool matched = true;
	      for(auto const & val : nRes) {
		if(val.second!=0) {
		  matched = false;
		  break;
		}
	      }
	      if (matched) {
		_h_onium[iq-4][n]->fill(is+1);
		break;
	      }
	    }
	    if(matched) break;
	  }
	  if(matched) break;
	}
	if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // normalize the cross sections
      for(const vector<Histo1DPtr> & hv : _h_m)
	for(const Histo1DPtr & hist : hv)
	  scale(hist, crossSection()/picobarn/sumW());
      for(const vector<Histo1DPtr> & hv : _h_onium)
	for(const Histo1DPtr & hist : hv)
	  scale(hist, crossSection()/picobarn/sumW());
    }

    ///@}


    /// @name Histograms
    ///@{
    vector<vector<Histo1DPtr>> _h_m;
    vector<vector<Histo1DPtr> > _h_onium;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MC_EEtoEEX);

}
