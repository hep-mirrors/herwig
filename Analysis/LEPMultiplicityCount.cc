// -*- C++ -*-
//
// LEPMultiplicityCount.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LEPMultiplicityCount class.
//

#include "LEPMultiplicityCount.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/StandardSelectors.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/EnumParticles.h"
#include "Herwig/Hadronization/Cluster.h"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace Herwig;
using namespace ThePEG;

IBPtr LEPMultiplicityCount::clone() const {
  return new_ptr(*this);
}

IBPtr LEPMultiplicityCount::fullclone() const {
  return new_ptr(*this);
}

LEPMultiplicityCount::LEPMultiplicityCount() : _makeHistograms(false)
{
  // Average particle multiplicities in hadronic Z decay
  // PDG 2006 with 2007 partial update

  // all charged
  _data[0] = MultiplicityInfo(20.76, 0.16, lightMeson); // all charged
  // gamma
  _data[22] = MultiplicityInfo(20.97, 1.17, lightMeson);
  // pi+, pi0, eta
  _data[211] = MultiplicityInfo(17.03, 0.16, lightMeson);
  _data[111] = MultiplicityInfo( 9.76, 0.26, lightMeson);
  _data[221] = MultiplicityInfo( 1.01, 0.08, lightMeson);
  // rho+, rho0, omega, eta'
  _data[213] = MultiplicityInfo( 2.40, 0.49, lightMeson);
  _data[113] = MultiplicityInfo( 1.24, 0.10, lightMeson);
  _data[223] = MultiplicityInfo( 1.02, 0.06, lightMeson);
  _data[331] = MultiplicityInfo( 0.17, 0.05, lightMeson);
  // f_0(980), a_0(980), phi
  _data[10221] = MultiplicityInfo(0.147, 0.011, other);
  _data[9000211] = MultiplicityInfo(0.27, 0.14, other);
  _data[333] = MultiplicityInfo(0.098, 0.006, strangeMeson);
  // f_2(1270), f_1(1285), f_2'(1525), K+, K0
  _data[225] = MultiplicityInfo(0.169, 0.025, other);
  _data[20223] = MultiplicityInfo(0.165, 0.051, other);
  _data[335] = MultiplicityInfo(0.012, 0.006, other);  
  _data[321] = MultiplicityInfo(2.24, 0.04, strangeMeson);
  _data[311] = MultiplicityInfo(2.039, 0.025, lightMeson);
  // K*+(892), K*0(892), K*_2(1430)
  _data[323] = MultiplicityInfo(0.72, 0.05, strangeMeson);
  _data[313] = MultiplicityInfo(0.739, 0.022, strangeMeson);
  _data[315] = MultiplicityInfo(0.073, 0.023, strangeMeson);
  // D+, D0, D_s+
  _data[411] = MultiplicityInfo(0.187, 0.020, other);
  _data[421] = MultiplicityInfo(0.462, 0.026, other);
  _data[431] = MultiplicityInfo(0.131, 0.028, other);
  // D*+(2010), J/Psi(1S), Psi(2S)
  _data[413] = MultiplicityInfo(0.183, 0.008, other);
  _data[443] = MultiplicityInfo(0.0056, 0.0007, other);
  _data[100443] = MultiplicityInfo(0.0023, 0.0007, other);
  // p, Delta++(1232), Lambda, Sigma+, Sigma-, Sigma0
  _data[2212] = MultiplicityInfo(1.046, 0.026, lightBaryon);
  _data[2224] = MultiplicityInfo(0.087, 0.033, lightBaryon);
  _data[3122] = MultiplicityInfo(0.388, 0.009, lightBaryon);
  _data[3222] = MultiplicityInfo(0.107, 0.010, lightBaryon);
  _data[3112] = MultiplicityInfo(0.082, 0.007, lightBaryon);
  _data[3212] = MultiplicityInfo(0.076, 0.010, lightBaryon);
  // Sigma*+, Sigma*-, Xi-, Xi*0, Omega-
  _data[3224] = MultiplicityInfo(0.0239, 0.0021, lightBaryon);  
  _data[3114] = MultiplicityInfo(0.0240, 0.0024, lightBaryon);  
  _data[3312] = MultiplicityInfo(0.0258, 0.0009, lightBaryon);
  _data[3324] = MultiplicityInfo(0.0059, 0.0011, lightBaryon);
  _data[3334] = MultiplicityInfo(0.00164, 0.00028, lightBaryon);
  // Lambda_c+
  _data[4122] = MultiplicityInfo(0.078, 0.024, other);

  // old values from 1.0 paper
//   _data[433] = MultiplicityInfo(0.096, 0.046, other);
  _data[2112] = MultiplicityInfo(0.991, 0.054, lightBaryon);
//   _data[2214] = MultiplicityInfo(0., 0., lightBaryon);
//   _data[2114] = MultiplicityInfo(0., 0., lightBaryon);

// values unknown
  // B mesons
  //  _data[513] = MultiplicityInfo(0.28, 0.04, other);   // flavour averaged value!
  _data[513] = MultiplicityInfo(0., 0., other);
  _data[511] = MultiplicityInfo(0., 0., other); // B0
  _data[521] = MultiplicityInfo(0., 0., other); // B+
  _data[531] = MultiplicityInfo(0., 0., other); // B_s
  _data[541] = MultiplicityInfo(0., 0., other); // B_c
  // B baryons
  _data[5122] = MultiplicityInfo(0., 0., other); // Lambda_b
  _data[5112] = MultiplicityInfo(0., 0., other); // Sig_b-
  _data[5212] = MultiplicityInfo(0., 0., other); // Sig_b0
  _data[5222] = MultiplicityInfo(0., 0., other); // Sig_b+
  _data[5132] = MultiplicityInfo(0., 0., other); // Xi_b-
  _data[5232] = MultiplicityInfo(0., 0., other); // Xi_b0
  _data[5312] = MultiplicityInfo(0., 0., other); // Xi'_b-
  _data[5322] = MultiplicityInfo(0., 0., other); // Xi'_b0
  _data[5332] = MultiplicityInfo(0., 0., other); // Omega_b-
}

namespace {
  bool isLastCluster(tcPPtr p) {
    if ( p->id() != ParticleID::Cluster ) 
      return false;
    for ( size_t i = 0, end = p->children().size();
	  i < end; ++i ) {
      if ( p->children()[i]->id() == ParticleID::Cluster )
	return false;
    }
    return true;
  }

  Energy parentClusterMass(tcPPtr p) {
    if (p->parents().empty()) 
      return -1.0*MeV;

    tcPPtr parent = p->parents()[0];
    if (parent->id() == ParticleID::Cluster) {
      if ( isLastCluster(parent) )
	return parent->mass();
      else
	return p->mass();
    }
    else
      return parentClusterMass(parent);
  }

  bool isPrimaryCluster(tcPPtr p) {
    if ( p->id() != ParticleID::Cluster ) 
      return false;
    if( p->parents().empty())
      return false;
    for ( size_t i = 0, end = p->parents().size();
	  i < end; ++i ) {
      if ( !(p->parents()[i]->dataPtr()->coloured()) )
	return false;
    }
    return true;
  }
}


void LEPMultiplicityCount::analyze(tEventPtr event, long, int, int) {

  set<tcPPtr> particles;
  event->selectFinalState(inserter(particles));

  map <long,long> eventcount;

  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {

    if((*it)->dataPtr()->charged()) 
      ++eventcount[0];
    
    long ID = abs( (*it)->id() );
    ++_finalstatecount[ID];
  }

  // ========

  StepVector steps = event->primaryCollision()->steps();

  particles.clear();
  steps[0]->select(inserter(particles), ThePEG::AllSelector());
  
  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    long ID = (*it)->id();
    ++_collisioncount[ID];
  }

  // =======


  particles.clear();

  if (steps.size() > 2) {
    for ( StepVector::const_iterator it = steps.begin()+2;
	  it != steps.end(); ++it ) {
      (**it).select(inserter(particles), ThePEG::AllSelector());
    }
  }
  
  if( _makeHistograms ) 
    _histograms.insert(make_pair(ParticleID::Cluster, 
				 Histogram(0.0,10.0,200)));
  
  for(set<tcPPtr>::const_iterator it = particles.begin(); 
      it != particles.end(); ++it) {
    long ID = abs( (*it)->id() );
    if(ID==ParticleID::K0) continue;
    if(ID==ParticleID::K_L0||ID==ParticleID::K_S0) ID=ParticleID::K0;
    
    if ( _makeHistograms && isLastCluster(*it) ) {
      _histograms[ParticleID::Cluster] += (*it)->mass()/GeV;
      tcClusterPtr clu = dynamic_ptr_cast<tcClusterPtr>(*it);
      if (clu) {
	_clusters.insert(make_pair(clu->clusterId(), Histogram(0.0,10.0,200)));
	_clusters[clu->clusterId()] += (*it)->mass()/GeV;
      }
    }

    if( _makeHistograms && isPrimaryCluster(*it) ) {
      _primary.insert(make_pair(0, Histogram(0.0,20.0,400)));
      _primary[0] += (*it)->mass()/GeV;
      tcClusterPtr clu = dynamic_ptr_cast<tcClusterPtr>(*it);
      if(clu) {
	_primary.insert(make_pair(clu->clusterId(), Histogram(0.0,20.0,400)));
	_primary[clu->clusterId()] += (*it)->mass()/GeV;
      }
    }
    
    if (_data.find(ID) != _data.end()) {
      ++eventcount[ID];
      
      if (_makeHistograms 
	  && ! (*it)->parents().empty()
	  && (*it)->parents()[0]->id() == ParticleID::Cluster) {
	_histograms.insert(make_pair(ID,Histogram(0.0,10.0,200)));
	_histograms[ID] += parentClusterMass(*it)/GeV;
      }
    }
  }
  
  for(map<long,MultiplicityInfo>::iterator it = _data.begin();
      it != _data.end(); ++it) {
    long currentcount 
      = eventcount.find(it->first) == eventcount.end() ? 0
      : eventcount[it->first];
    it->second.count += currentcount; 
  }
}

void LEPMultiplicityCount::analyze(const tPVector & ) {}

void LEPMultiplicityCount::dofinish() {
  useMe();
  string filename = generator()->filename() + ".mult";
  ofstream outfile(filename.c_str());

  outfile << 
    "\nParticle multiplicities (compared to LEP data):\n"
    "  ID       Name    simMult     obsMult       obsErr     Sigma\n";
  for (map<long,MultiplicityInfo>::const_iterator it = _data.begin();
       it != _data.end();
       ++it) 
    {
      MultiplicityInfo multiplicity = it->second;
      string name = (it->first==0 ? "All chgd" : 
		     generator()->getParticleData(it->first)->PDGName() );

      ios::fmtflags oldFlags = outfile.flags();
      outfile << std::scientific << std::showpoint
	      << std::setprecision(3)
	      << setw(7) << it->first << ' '
	      << setw(9) << name << ' ' 
	      << setw(2) << multiplicity.simMultiplicity() << " | " 
	      << setw(2) << multiplicity.obsMultiplicity << " +/- " 
	      << setw(2) << multiplicity.obsError << ' '
	      << std::showpos << std::setprecision(1)
	      << multiplicity.nSigma() << ' ' 
	      << multiplicity.bargraph()
	      << std::noshowpos;

      outfile << '\n';
      outfile.flags(oldFlags);
    }



  outfile << "\nCount of particles involved in hard process:\n";
  for (map<long,long>::const_iterator it = _collisioncount.begin();
       it != _collisioncount.end(); ++ it) {
    string name = generator()->getParticleData(it->first)->PDGName();
    outfile << name << '\t' << it->second << '\n';
  }



  outfile << "\nFinal state particle count:\n";
  for (map<long,long>::const_iterator it = _finalstatecount.begin();
       it != _finalstatecount.end(); ++ it) {
    string name = generator()->getParticleData(it->first)->PDGName();
    outfile << name << '\t' << it->second << '\n';
  }
  outfile.close();

  if (_makeHistograms) {
    
    Histogram piratio 
      = _histograms[ParticleID::piplus].ratioWith(_histograms[ParticleID::pi0]);
    Histogram Kratio 
      = _histograms[ParticleID::Kplus].ratioWith(_histograms[ParticleID::K0]);

    using namespace HistogramOptions;
    string histofilename = filename + ".top";
    ofstream outfile2(histofilename.c_str());

    for (map<int,Histogram>::const_iterator it = _primary.begin();
	 it != _primary.end(); ++it) {
      ostringstream title1;
      title1 << "Primary Cluster " << it->first;
      string title = title1.str();
      it->second.topdrawOutput(outfile2,Frame|Ylog,"BLACK",title,"",
			       "N (200 bins)","","Cluster mass [GeV]");
    }
    map<long,Histogram>::const_iterator cit = _histograms.find(ParticleID::Cluster);
    string title = generator()->getParticleData(cit->first)->PDGName();
    cit->second.topdrawOutput(outfile2,Frame|Ylog,"BLACK",title,"",
			     "N (200 bins)","","Parent cluster mass [GeV]");

    for (map<int,Histogram>::const_iterator it = _clusters.begin();
	 it != _clusters.end(); ++it) {
      ostringstream title1;
      title1 << "Final Cluster " << it->first;
      string title = title1.str();
      it->second.topdrawOutput(outfile2,Frame|Rawcount|Ylog,"BLACK",title,"",
			       "N (200 bins)","","Cluster mass [GeV]");
    }
    for (map<long,Histogram>::const_iterator it = _histograms.begin();
	 it != _histograms.end(); ++it) {
      string title = generator()->getParticleData(it->first)->PDGName();
      it->second.topdrawOutput(outfile2,Frame|Rawcount|Ylog,"BLACK",title,"",
			       "N (200 bins)","","Parent cluster mass [GeV]");
    }
    piratio.topdrawOutput(outfile2,Frame|Rawcount,"BLACK","pi+ / pi0","",
			  "","","Parent cluster mass [GeV]");
    Kratio.topdrawOutput(outfile2,Frame|Rawcount,"BLACK","K+ / K0","",
			 "","","Parent cluster mass [GeV]");
    outfile2.close();
  }
  AnalysisHandler::dofinish();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<LEPMultiplicityCount,AnalysisHandler>
describeHerwigLEPMultiplicityCount("Herwig::LEPMultiplicityCount", "HwAnalysis.so");

void LEPMultiplicityCount::Init() {

  static ParVector<LEPMultiplicityCount,long> interfaceparticlecodes
    ("ParticleCodes",
     "The PDG code for the particles",
     &LEPMultiplicityCount::_particlecodes,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<LEPMultiplicityCount,double> interfaceMultiplicity
    ("Multiplicity",
     "The multiplicity for the particle",
     &LEPMultiplicityCount::_multiplicity,
     0, 0, 0, 0., 1000., false, false, true);

  static ParVector<LEPMultiplicityCount,double> interfaceError
    ("Error",
     "The error on the multiplicity for the particle",
     &LEPMultiplicityCount::_error,
     0, 0, 0, 0., 1000., false, false, true);

  static ParVector<LEPMultiplicityCount,unsigned int> interfaceSpecies
    ("Species",
     "The type of particle",
     &LEPMultiplicityCount::_species,
     0, 0, other, 0, other, false, false, true);


  static Switch<LEPMultiplicityCount,bool> interfaceHistograms
    ("Histograms",
     "Set to On if detailed histograms are required.",
     &LEPMultiplicityCount::_makeHistograms, false, true, false);
  static SwitchOption interfaceHistogramsYes
    (interfaceHistograms,
     "Yes",
     "Generate histograms of cluster mass dependence.",
     true);
  static SwitchOption interfaceHistogramsNo
    (interfaceHistograms,
     "No",
     "Do not generate histograms.",
     false);

  static ClassDocumentation<LEPMultiplicityCount> documentation
    ("The LEPMultiplicityCount class count the multiplcities of final-state particles"
     " and compares them with LEP data.",
     "The LEP multiplicity analysis uses data from PDG 2006 \\cite{Yao:2006px}.",
     "%\\cite{Yao:2006px}\n"
     "\\bibitem{Yao:2006px}\n"
     "  W.~M.~Yao {\\it et al.}  [Particle Data Group],\n"
     "  %``Review of particle physics,''\n"
     "  J.\\ Phys.\\ G {\\bf 33} (2006) 1.\n"
     "  %%CITATION = JPHGB,G33,1;%%\n"
     );

}

void LEPMultiplicityCount::persistentOutput(PersistentOStream & os) const {
  os << _particlecodes << _multiplicity << _error << _species << _makeHistograms;
}

void LEPMultiplicityCount::persistentInput(PersistentIStream & is, int) {
  is >> _particlecodes >> _multiplicity >> _error >> _species >> _makeHistograms;
}
