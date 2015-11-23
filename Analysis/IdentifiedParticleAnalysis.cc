// -*- C++ -*-
//
// IdentifiedParticleAnalysis.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IdentifiedParticleAnalysis class.
//

#include "IdentifiedParticleAnalysis.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

int IdentifiedParticleAnalysis::getFlavour(const tPVector &pv) {
  tPVector::const_iterator it;
  for(it = pv.begin(); it!=pv.end(); ++it) 
    if (abs((*it)->id()) < 7) break; 
  return abs((*it)->id());
}

void IdentifiedParticleAnalysis::
analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  if ( loop > 0 || state != 0 || !event ) return;
  // get the final-state
  tPVector hadrons=event->getFinalState();
  // get the partons 
  tPVector partons=event->primaryCollision()->steps()[0]->
    getFinalState();
  int flav = getFlavour(partons);
  Energy Emax = 0.5*generator()->maximumCMEnergy();  
  for (tPVector::iterator it = hadrons.begin(); 
       it != hadrons.end(); ++it ) {
    // only looking at charged particles
    if(!(*it)->data().charged()) continue;
    // all particles
    double xp = getX((*it)->momentum(), Emax);
    double pp = (*it)->momentum().vect().mag()/GeV;
    _xpa->addWeighted(xp,event->weight());
    if(abs((*it)->id()) == ParticleID::piplus) {
      _pipma->addWeighted(xp,event->weight());
      _pipm ->addWeighted(pp,event->weight()); 
    } 
    else if(abs((*it)->id()) == ParticleID::Kplus) {
      _kpma->addWeighted(xp,event->weight()); 
      _kpm ->addWeighted(pp,event->weight());
    }
    else if(abs((*it)->id()) == ParticleID::pplus) {
      _ppma->addWeighted(xp,event->weight()); 
      _ppm ->addWeighted(pp,event->weight());
    }
    switch(flav) {
    case 1:
    case 2:
    case 3:
      _xpl ->addWeighted(xp,event->weight()); 
      if(abs((*it)->id()) == ParticleID::piplus)     
	_pipml->addWeighted(xp,event->weight());
      else if(abs((*it)->id()) == ParticleID::Kplus) 
	_kpml ->addWeighted(xp,event->weight());
      else if(abs((*it)->id()) == ParticleID::pplus) 
	_ppml ->addWeighted(xp,event->weight());
      _udsxp->addWeighted(xp,event->weight());
      if (xp > 0) 
	_udsxip->addWeighted( -log(xp),event->weight());
      break;
    case 4: 
      _xpc->addWeighted(xp,event->weight()); 
      if(abs((*it)->id()) == ParticleID::piplus)
	_pipmc->addWeighted(xp,event->weight()); 
      else if(abs((*it)->id()) == ParticleID::Kplus)
	_kpmc ->addWeighted(xp,event->weight()); 
      else if(abs((*it)->id()) == ParticleID::pplus) 
	_ppmc ->addWeighted(xp,event->weight());
      break;
    case 5:
      _xpb ->addWeighted(xp,event->weight()); 
      if(abs((*it)->id()) == ParticleID::piplus) 
	_pipmb ->addWeighted(xp,event->weight());
      else if(abs((*it)->id()) == ParticleID::Kplus) 
	_kpmb ->addWeighted(xp,event->weight());
      else if(abs((*it)->id()) == ParticleID::pplus) 
	_ppmb ->addWeighted(xp,event->weight());
      break;
    default:
      break;
    }
  }
  // finally decaying particles
  set<tcPPtr> allparticles;
  StepVector steps = event->primaryCollision()->steps();
  for ( StepVector::const_iterator it = steps.begin()+2;
	it != steps.end(); ++it ) {
    (**it).select(inserter(allparticles), ThePEG::AllSelector());
  }
  
  for(set<tcPPtr>::const_iterator it = allparticles.begin(); 
      it != allparticles.end(); ++it) {
    // lambda's
    long id = abs( (*it)->id());
    double xe = (*it)->momentum().e()/Emax;
    double xp = (*it)->momentum().vect().mag()/Emax;
    switch(id) {
    case ParticleID::Lambda0:
      _lpm ->addWeighted( xp,event->weight());
      break;
    case ParticleID::Kstarplus:
      _xpKstarplus ->addWeighted( xp,event->weight());
      break;
    case ParticleID::Ximinus:
      _xpXiminus ->addWeighted( xe,event->weight());
      _xiXiminus ->addWeighted( -log(xp),event->weight());
      break;
    case ParticleID::Sigmastarplus:
      _xpSigmaplus ->addWeighted( xe,event->weight());
      _xiSigmaplus ->addWeighted( -log(xp),event->weight());
      break;
    case ParticleID::Sigmastarminus:
      _xpSigmaminus ->addWeighted( xe,event->weight());
      _xiSigmaminus ->addWeighted( -log(xp),event->weight());
      break;
    case ParticleID::Xistar0:
      _xpXi0 ->addWeighted( xe,event->weight());
      _xiXi0 ->addWeighted( -log(xp),event->weight());
      break;
    case 3124: // lambda(1520)
      _xpLambda1520 ->addWeighted( xe,event->weight());
      _xiLambda1520 ->addWeighted(-log(xp),event->weight());
      break;
    case ParticleID::Deltaplus2:
      _xeDelta ->addWeighted(xe,event->weight());
      break;
    case ParticleID::f_0:
      _xpf980 ->addWeighted(xp,event->weight());
      break;
    case ParticleID::f_2:
      _xpf2 ->addWeighted(xp,event->weight());
      break;
    case ParticleID::phi:
      _xpphi ->addWeighted(xp,event->weight());
      break;
    case ParticleID::Kstar0:
      _xpKstar0 ->addWeighted(xp,event->weight());
      break;
    case ParticleID::D0:
      _xeD0 ->addWeighted( xe,event->weight());
      break;
    case ParticleID::Dstarplus:
      _xeDstar ->addWeighted( xe,event->weight());
      break;
    case ParticleID::rho0:
      _xerho0 ->addWeighted( xe,event->weight());
      break;
    case ParticleID::pi0:
      _xepi0 ->addWeighted( xe,event->weight());
      _xipi0 ->addWeighted(-log(xp),event->weight());
    case ParticleID::eta:
      _xeeta ->addWeighted( xe,event->weight());
      _xieta ->addWeighted(-log(xp),event->weight());
    case ParticleID::etaprime:
      _xeetap ->addWeighted( xe,event->weight());
      _xietap ->addWeighted(-log(xp),event->weight());
    case ParticleID::rhoplus:
      _xerhop ->addWeighted( xe,event->weight());
      _xirhop ->addWeighted(-log(xp),event->weight());
    case ParticleID::omega:
      _xeomega ->addWeighted( xe,event->weight());
      _xiomega ->addWeighted(-log(xp),event->weight());
    case ParticleID::a_0plus:
      _xea_0p ->addWeighted( xe,event->weight());
      _xia_0p ->addWeighted(-log(xp),event->weight());
    case ParticleID::K0:
    case ParticleID::K_S0:
    case ParticleID::K_L0:
      _xpK0->addWeighted(xp,event->weight());
    }
  }
}

NoPIOClassDescription<IdentifiedParticleAnalysis> 
IdentifiedParticleAnalysis::initIdentifiedParticleAnalysis;
// Definition of the static class description member.

void IdentifiedParticleAnalysis::Init() {

  static ClassDocumentation<IdentifiedParticleAnalysis> documentation
    ("The IdentifiedParticleAnalysis class compares identified"
     " particle spectra with Z pole data",
     "The LEP IdentifiedParticle analysis uses data from"
     "\\cite{Akers:1994ez,Alexander:1995gq,Alexander:1996qj,Ackerstaff:1998ue,Ackerstaff:1997kj,Abbiendi:2000cv,Ackerstaff:1998ap,Acton:1991aa,Abreu:1998nn,Abreu:1993mn,Barate:1999bg,Barate:1996fi,Abe:1998zs}.",
     "%\\cite{Akers:1994ez}\n"
     "\\bibitem{Akers:1994ez}\n"
     "  R.~Akers {\\it et al.}  [OPAL Collaboration],\n"
     "   ``Measurement of the production rates of charged hadrons in e+ e-\n"
     "  %annihilation at the Z0,''\n"
     "  Z.\\ Phys.\\  C {\\bf 63}, 181 (1994).\n"
     "  %%CITATION = ZEPYA,C63,181;%%\n"
     "%\\cite{Alexander:1995gq}\n"
     "\\bibitem{Alexander:1995gq}\n"
     "  G.~Alexander {\\it et al.}  [OPAL Collaboration],\n"
     "  %``Delta++ production in hadronic Z0 decays,''\n"
     "  Phys.\\ Lett.\\  B {\\bf 358}, 162 (1995).\n"
     "  %%CITATION = PHLTA,B358,162;%%\n"
     "%\\cite{Alexander:1996qj}\n"
     "\\bibitem{Alexander:1996qj}\n"
     "  G.~Alexander {\\it et al.}  [OPAL Collaboration],\n"
     "  %``Strange baryon production in hadronic Z0 decays,''\n"
     "  Z.\\ Phys.\\  C {\\bf 73}, 569 (1997).\n"
     "  %%CITATION = ZEPYA,C73,569;%%\n"
     "%\\cite{Ackerstaff:1998ue}\n"
     "\\bibitem{Ackerstaff:1998ue}\n"
     "  K.~Ackerstaff {\\it et al.}  [OPAL Collaboration],\n"
     "  %``Production of f0(980), f2(1270) and Phi(1020) in hadronic Z0 decay,''\n"
     "  Eur.\\ Phys.\\ J.\\  C {\\bf 4}, 19 (1998)\n"
     "  [arXiv:hep-ex/9802013].\n"
     "  %%CITATION = EPHJA,C4,19;%%\n"
     "%\\cite{Ackerstaff:1997kj}\n"
     "\\bibitem{Ackerstaff:1997kj}\n"
     "  K.~Ackerstaff {\\it et al.}  [OPAL Collaboration],\n"
     "  %``Spin alignment of leading K*(892)0 mesons in hadronic Z0 decays,''\n"
     "  Phys.\\ Lett.\\  B {\\bf 412}, 210 (1997)\n"
     "  [arXiv:hep-ex/9708022].\n"
     "  %%CITATION = PHLTA,B412,210;%%\n"
     "%\\cite{Abbiendi:2000cv}\n"
     "\\bibitem{Abbiendi:2000cv}\n"
     "  G.~Abbiendi {\\it et al.}  [OPAL Collaboration],\n"
     "   ``Multiplicities of pi0, eta, K0 and of charged particles in quark and  gluon\n"
     "  %jets,''\n"
     "  Eur.\\ Phys.\\ J.\\  C {\\bf 17}, 373 (2000)\n"
     "  [arXiv:hep-ex/0007017].\n"
     "  %%CITATION = EPHJA,C17,373;%%\n"
     "%\\cite{Ackerstaff:1998ap}\n"
     "\\bibitem{Ackerstaff:1998ap}\n"
     "  K.~Ackerstaff {\\it et al.}  [OPAL Collaboration],\n"
     "  %``Photon and light meson production in hadronic Z0 decays,''\n"
     "  Eur.\\ Phys.\\ J.\\  C {\\bf 5}, 411 (1998)\n"
     "  [arXiv:hep-ex/9805011].\n"
     "  %%CITATION = EPHJA,C5,411;%%\n"
     "%\\cite{Acton:1991aa}\n"
     "\\bibitem{Acton:1991aa}\n"
     "  P.~D.~Acton {\\it et al.}  [OPAL Collaboration],\n"
     "  %``A Study of charged particle multiplicities in hadronic decays of the Z0,''\n"
     "  Z.\\ Phys.\\  C {\\bf 53}, 539 (1992).\n"
     "  %%CITATION = ZEPYA,C53,539;%%\n"
     "%\\cite{Abreu:1998nn}\n"
     "\\bibitem{Abreu:1998nn}\n"
     "  P.~Abreu {\\it et al.}  [DELPHI Collaboration],\n"
     "   ``Measurement of inclusive rho0, f0(980), f2(1270), K*2(1430)0  and f\'2(1525)\n"
     "  %production in Z0 decays,''\n"
     "  Phys.\\ Lett.\\  B {\\bf 449}, 364 (1999).\n"
     "  %%CITATION = PHLTA,B449,364;%%\n"
     "%\\cite{Abreu:1993mn}\n"
     "\\bibitem{Abreu:1993mn}\n"
     "  P.~Abreu {\\it et al.}  [DELPHI Collaboration],\n"
     "  %``A Measurement of D meson production in Z0 hadronic decays,''\n"
     "  Z.\\ Phys.\\  C {\\bf 59}, 533 (1993)\n"
     "  [Erratum-ibid.\\  C {\\bf 65}, 709 (1995)].\n"
     "  %%CITATION = ZEPYA,C59,533;%%\n"
     "%\\cite{Barate:1999bg}\n"
     "\\bibitem{Barate:1999bg}\n"
     "  R.~Barate {\\it et al.}  [ALEPH Collaboration],\n"
     "  %``Study of charm production in Z decays,''\n"
     "  Eur.\\ Phys.\\ J.\\  C {\\bf 16}, 597 (2000)\n"
     "  [arXiv:hep-ex/9909032].\n"
     "  %%CITATION = EPHJA,C16,597;%%\n"
     "%\\cite{Barate:1996fi}\n"
     "\\bibitem{Barate:1996fi}\n"
     "  R.~Barate {\\it et al.}  [ALEPH Collaboration],\n"
     "  %``Studies of quantum chromodynamics with the ALEPH detector,''\n"
     "  Phys.\\ Rept.\\  {\\bf 294}, 1 (1998).\n"
     "  %%CITATION = PRPLC,294,1;%%\n"
     "%\\cite{Abe:1998zs}\n"
     "\\bibitem{Abe:1998zs}\n"
     "  K.~Abe {\\it et al.}  [SLD Collaboration],\n"
     "   ``Production of pi+, K+, K0, K*0, Phi, p and Lambda0 in hadronic Z0\n"
     "  %decays,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 59}, 052001 (1999)\n"
     "  [arXiv:hep-ex/9805029].\n"
     "  %%CITATION = PHRVA,D59,052001;%%\n"
     );
}

void IdentifiedParticleAnalysis::dofinish() {
  useMe();
  AnalysisHandler::dofinish();
  string fname = generator()->filename() 
    + string("-") + name() + string(".top");
  ofstream output(fname.c_str());
  // chisq
  double chisq,minfrac=0.05;
  unsigned int npoint;
  using namespace HistogramOptions;
  // Histogram for the \f$\xi\f$ distribution for all particles from all quarks
  _xpa->normaliseToData();
  _xpa->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (charged, all quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _xpa->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of charged particles for all events(SLD)",
		      "                                                       ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // light quarks
  _xpl->normaliseToData();
  _xpl->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (charged, light quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _xpl->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of charged particles for all light quark events(SLD)",
		      "                                                                   ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // OPAL
  _udsxp->normaliseToData();
  _udsxp->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL xp (charged, light quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _udsxp->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of charged particles for all light quark events(OPAL)",
		      "                                                                   ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  _udsxip->normaliseToData();
  _udsxip->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL xi (charged, light quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _udsxip->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of charged particles for all light quark events(OPAL)",
		      "                                                                   ",
		      "1/NdN/dX",
		      "       G",
		      "X",
		      "G");
  // charm 
  _xpc->normaliseToData();
  _xpc->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (charged, charm quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _xpc->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of charged particles for all charm events(SLD)",
		      "                                                             ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // bottom
  _xpb->normaliseToData();
  _xpb->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (charged, bottom quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xpb->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of charged particles for all bottom events(SLD)",
		      "                                                              ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // pions all quarks
  _pipma->normaliseToData();
  _pipma->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (pions, all quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _pipma->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of pions for all events(SLD)",
		      "                                           ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // OPAL pions
  _pipm->normaliseToData();
  _pipm->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL momentum (pions, all quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _pipm->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The momentum of pions for all events(OPAL)",
		      "                                          ",
		      "1/NdN/dp/GeV2-13",
		      "            X  X",
		      "p/GeV",
		      " ");
  // light
  _pipml->normaliseToData();
  _pipml->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (pions, light quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _pipml->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of pions for light quark events(SLD)",
		      "                                           ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // charm 
  _pipmc->normaliseToData();
  _pipmc->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (pions, charm quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _pipmc->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of pions for charm events(SLD)",
		      "                                           ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // bottom
  _pipmb->normaliseToData();
  _pipmb->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (pions, bottom quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _pipmb->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of pions for bottom quark events(SLD)",
		      "                                           ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // kaons all quarks
  _kpma->normaliseToData();
  _kpma->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (kaons, all quakrks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _kpma->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of kaons for all events(SLD)",
		      "                                           ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // OPAL kaons
  _kpm->normaliseToData();
  _kpm->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL momentum (kaons, all quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _kpm->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The momentum of kaons for all events(OPAL)",
		      "                                                   ",
		      "1/NdN/dp/GeV2-13",
		      "            X  X",
		      "p/GeV",
		      "     ");
  // light
  _kpml->normaliseToData();
  _kpml->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (kaons, light quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _kpml->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of kaons for light quark events(SLD)",
		      "                                           ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // charm 
  _kpmc->normaliseToData();
  _kpmc->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (kaons, charm quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _kpmc->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of kaons for charm events(SLD)",
		      "                                           ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // bottom
  _kpmb->normaliseToData();
  _kpmb->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (kaons, bottom quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _kpmb->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of kaons for bottom quark events(SLD)",
		      "                                           ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // protons all quarks
  _ppma->normaliseToData();
  _ppma->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (protons, all quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _ppma->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of protons for all events(SLD)",
		      "                                           ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // OPAL protons
  _ppm->normaliseToData();
  _ppm->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL momentum (protons, all quarks ) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _ppm->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The momentum of protons for all events(OPAL)",
		      "                                                   ",
		      "1/NdN/dpGeV2-13",
		      "           X  X",
		      "p/GeV",
		      "     ");
  // light
  _ppml->normaliseToData();
  _ppml->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (protons, light quarks) "
		     << "distribution or " << chisq/npoint << "per degree of freedom\n";
  _ppml->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of protons for light quark events(SLD)",
		      "                                           ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // charm 
  _ppmc->normaliseToData();
  _ppmc->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (protons, charm quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _ppmc->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of protons for charm events(SLD)",
		      "                                           ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // bottom
  _ppmb->normaliseToData();
  _ppmb->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for SLD xp (protons, bottom quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _ppmb->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of protons for bottom quark events(SLD)",
		      "                                           ",
		      "1/NdN/dx0p1",
		      "        X X",
		      "x0p1",
		      " X X");
  // lambda
  _lpm->normaliseToData();
  _lpm->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for ALEPH momentum (lambda, all quarks) "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n"; 
  _lpm->topdrawOutput(output,Frame|Errorbars|Ylog,
		      "RED",
		      "The scaled momentum of L for all events (ALEPH)",
		      "                       F                      ",
		      "1/NdN/dp/GeV2-13",
		      "            X  X",
		      "p/GeV",
		      "     ");
  // K*+
  _xpKstarplus->normaliseToData();
  _xpKstarplus->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for ALEPH momentum K*+ "
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xpKstarplus->topdrawOutput(output,Frame|Errorbars|Ylog,
			      "RED",
			      "The scaled momentum of K2*+3 for all events (ALEPH)",
			      "                        X  X                       ",
			      "1/NdN/dx0p1",
			      "        X X",
			      "x0p1",
			      " X X");
  // xi-
  _xpXiminus->normaliseToData();
  _xpXiminus->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL xi- x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xpXiminus->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum x0E1 of X2-3 for all events (OPAL)",
			    "                     X X    FX X                      ",
			    "1/NdN/dx0E1",
			    "        X X",
			    "x0E1",
			    " X X");
  _xiXiminus->normaliseToData();
  _xiXiminus->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL xi- xi"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xiXiminus->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum X of X2-3 for all events (OPAL)",
			    "                    G    FX X                      ",
			    "1/NdN/dX",
			    "       G",
			    "X",
			    "G");
  // Sigma*+
  _xpSigmaplus->normaliseToData();
  _xpSigmaplus->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL Sigma*+ x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xpSigmaplus->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum x0E1 of S2*+3 for all events (OPAL)",
			    "                     X X    FX  X                      ",
			    "1/NdN/dx0E1",
			    "        X X",
			    "x0E1",
			    " X X");
  _xiSigmaplus->normaliseToData();
  _xiSigmaplus->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL Sigma*+ xi"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xiSigmaplus->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum X of S2*+3 for all events (OPAL)",
			    "                    G    FX  X                      ",
			    "1/NdN/dX",
			    "       G",
			    "X",
			    "G");
  // Sigma*-
  _xpSigmaminus->normaliseToData();
  _xpSigmaminus->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL Sigma*- x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xpSigmaminus->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum x0E1 of S2*-3 for all events (OPAL)",
			    "                     X X    FX  X                      ",
			    "1/NdN/dx0E1",
			    "        X X",
			    "x0E1",
			    " X X");
  _xiSigmaminus->normaliseToData();
  _xiSigmaminus->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL Sigma*- xi"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xiSigmaminus->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum X of S2*-3 for all events (OPAL)",
			    "                    G    FX  X                      ",
			    "1/NdN/dX",
			    "       G",
			    "X",
			    "G");
  // Xi*0
  _xpXi0->normaliseToData();
  _xpXi0->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL Xi*0 x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xpXi0->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum x0E1 of X2*03 for all events (OPAL)",
			    "                     X X    FX  X                      ",
			    "1/NdN/dx0E1",
			    "        X X",
			    "x0E1",
			    " X X");
  _xiXi0->normaliseToData();
  _xiXi0->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL Xi*0 xi"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xiXi0->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum X of X2*03 for all events (OPAL)",
			    "                    G    FX  X                      ",
			    "1/NdN/dX",
			    "       G",
			    "X",
			    "G");
  // lambda(1520)
  _xpLambda1520->normaliseToData();
  _xpLambda1520->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL Lambda(1520) x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xpLambda1520->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum x0E1 of L(1520) for all events (OPAL)",
			    "                     X X    F                            ",
			    "1/NdN/dx0E1",
			    "        X X",
			    "x0E1",
			    " X X");
  _xiLambda1520->normaliseToData();
  _xiLambda1520->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL Lambda(1520) xi"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xiLambda1520->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum X of L(1520) for all events (OPAL)",
			    "                    G    F                            ",
			    "1/NdN/dX",
			    "       G",
			    "X",
			    "G");
  // Delta++
  _xeDelta->normaliseToData();
  _xeDelta->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL Delta++ x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xeDelta->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum x0E1 of D2++3 for all events (OPAL)",
			    "                     X X    FX  X                      ",
			    "1/NdN/dx0E1",
			    "        X X",
			    "x0E1",
			    " X X");
  // f_0
  _xpf980->normaliseToData();
  _xpf980->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL f_0(980) x_p"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xpf980->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum x0p1 of f001(980) for all events (OPAL)",
			    "                     X X     X X                           ",
			    "1/NdN/dx0p1",
			    "        X X",
			    "x0p1",
			    " X X"); 
  // f_2
  _xpf2->normaliseToData();
  _xpf2->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL f_2 x_p"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xpf2->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum x0p1 of f021 for all events (OPAL)",
			    "                     X X     X X                      ",
			    "1/NdN/dx0p1",
			    "        X X",
			    "x0p1",
			    " X X");
  // phi
  _xpphi->normaliseToData();
  _xpphi->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL phi x_p"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xpphi->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum x0p1 of F for all events (OPAL)",
			    "                     X X    G                      ",
			    "1/NdN/dx0p1",
			    "        X X",
			    "x0p1",
			    " X X");
  // K*0
  _xpKstar0->normaliseToData();
  _xpKstar0->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL K*0 x_p"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xpKstar0->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum x0p1 of K2*03 for all events (OPAL)",
			    "                     X X     X  X                      ",
			    "1/NdN/dx0p1",
			    "        X X",
			    "x0p1",
			    " X X");
  // K0
  _xpK0->normaliseToData();
  _xpK0->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL K0 x_p"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xpK0->topdrawOutput(output,Frame|Errorbars|Ylog,
		       "RED",
		       "The scaled momentum x0p1 of K203 for all events (OPAL)",
		       "                     X X     X X                      ",
		       "1/NdN/dx0p1",
		       "        X X",
		       "x0p1",
		       " X X");
  // rho0
  _xerho0->normaliseToData();
  _xerho0->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for DELPHI rho0 x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xerho0->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum x0E1 of R203 for all events (DELPHI)",
			    "                     X X    GX X                        ",
			    "1/NdN/dx0E1",
			    "        X X",
			    "x0E1",
			    " X X");
  // pi0
  _xepi0->normaliseToData();
  _xepi0->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL pi0 x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xepi0->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"The scaled momentum x0E1 of P203 for all events (OPAL)",
			"                     X X    GX X                      ",
			"1/NdN/dx0E1",
			"        X X",
			"x0E1",
			" X X");
  _xipi0->normaliseToData();
  _xipi0->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL pi0 xi"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xipi0->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"The scaled momentum X of P203 for all events (OPAL)",
			"                    G    GX X                      ",
			"1/NdN/dX",
			"       G",
			"X",
			"G");
  // eta
  _xeeta->normaliseToData();
  _xeeta->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL eta x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xeeta->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"The scaled momentum x0E1 of H for all events (OPAL)",
			"                     X X    G                      ",
			"1/NdN/dx0E1",
			"        X X",
			"x0E1",
			" X X");
  _xieta->normaliseToData();
  _xieta->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL eta xi"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xieta->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"The scaled momentum X of H for all events (OPAL)",
			"                    G    G                      ",
			"1/NdN/dX",
			"       G",
			"X",
			"G");
  // eta'
  _xeetap->normaliseToData();
  _xeetap->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL eta' x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xeetap->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"The scaled momentum x0E1 of H' for all events (OPAL)",
			"                     X X    G                       ",
			"1/NdN/dx0E1",
			"        X X",
			"x0E1",
			" X X");
  _xietap->normaliseToData();
  _xietap->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL eta' xi"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xietap->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"The scaled momentum X of H' for all events (OPAL)",
			"                    G    G                       ",
			"1/NdN/dX",
			"       G",
			"X",
			"G");
  // rho+
  _xerhop->normaliseToData();
  _xerhop->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL rho+ x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xerhop->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"The scaled momentum x0E1 of R2+3 for all events (OPAL)",
			"                     X X    GX X                      ",
			"1/NdN/dx0E1",
			"        X X",
			"x0E1",
			" X X");
  _xirhop->normaliseToData();
  _xirhop->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL rho+ xi"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xirhop->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"The scaled momentum X of R2+3 for all events (OPAL)",
			"                    G    GX X                      ",
			"1/NdN/dX",
			"       G",
			"X",
			"G");
  // omega
  _xeomega->normaliseToData();
  _xeomega->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL omega x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xeomega->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"The scaled momentum x0E1 of W for all events (OPAL)",
			"                     X X    G                      ",
			"1/NdN/dx0E1",
			"        X X",
			"x0E1",
			" X X");
  _xiomega->normaliseToData();
  _xiomega->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL omega xi"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xiomega->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"The scaled momentum X of W for all events (OPAL)",
			"                    G    G                      ",
			"1/NdN/dX",
			"       G",
			"X",
			"G");
  // a_0+
  _xea_0p->normaliseToData();
  _xea_0p->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL a_0+ x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xea_0p->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"The scaled momentum x0E1 of A0012+3 for all events (OPAL)",
			"                     X X     X XX X                      ",
			"1/NdN/dx0E1",
			"        X X",
			"x0E1",
			" X X");
  _xia_0p->normaliseToData();
  _xia_0p->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for OPAL a_0+ xi"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xia_0p->topdrawOutput(output,Frame|Errorbars|Ylog,
			"RED",
			"The scaled momentum X of A0012+3 for all events (OPAL)",
			"                    G     X XX X                      ",
			"1/NdN/dX",
			"       G",
			"X",
			"G");
  // D0
  _xeD0->normaliseToData();
  _xeD0->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for DELPHI D0 x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xeD0->topdrawOutput(output,Frame|Errorbars|Ylog,
			    "RED",
			    "The scaled momentum x0E1 of D203 for all events (DELPHI)",
			    "                     X X     X X                        ",
			    "1/NdN/dx0E1",
			    "        X X",
			    "x0E1",
			    " X X");
  // D*+
  _xeDstar->normaliseToData();
  _xeDstar->chiSquared(chisq,npoint,minfrac);
  generator()->log() << "Chi Square = " << chisq << " for " << npoint 
		     << " degrees of freedom for ALEPH Dstar x_E"
		     << "distribution or " << chisq/npoint << " per degree of freedom\n";
  _xeDstar->topdrawOutput(output,Frame|Errorbars|Ylog,
			  "RED",
			  "The scaled momentum x0E1 of D2*+3 for all events (ALEPH)",
			  "                     X X     X  X                        ",
			  "1/NdN/dx0E1",
			  "        X X",
			  "x0E1",
			  " X X");
}

void IdentifiedParticleAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  // SLD data (all charged)
  double xpbins[] = {0.005, 0.008, 0.010, 0.012, 0.014, 
		     0.016, 0.022, 0.027, 0.033, 0.038, 
		     0.044, 0.049, 0.055, 0.060, 0.066, 
		     0.071, 0.077, 0.082, 0.088, 0.099, 
		     0.110, 0.121, 0.143, 0.164, 0.186, 
		     0.208, 0.230, 0.252, 0.274, 0.296, 
		     0.318, 0.351, 0.384, 0.417, 0.450, 
		     0.482, 0.526, 0.570, 0.658, 0.768, 
		     1.000};
  double xpdataa[]={509.2   ,513.9   ,485.5   ,443.4   ,398.7   ,
		    335.8   ,267.2   ,217.8   ,180.7   ,153.2   ,
		    130.6   ,113.4   , 99.09  , 87.57  , 77.56  ,
		    69.78   , 62.84  , 56.72  , 48.90  , 40.78  ,
		    34.39   , 27.35  , 20.53  , 15.64  , 12.25  ,
		    9.67    ,  7.75  , 6.161  , 5.029  , 4.053  ,
		    3.139   , 2.338  , 1.748  , 1.326  , 1.008  ,
		    0.724   , 0.480  , 0.285  , 0.114  , 0.024  };
  double xperrora[]={9.5  ,7.8  ,6.6  ,5.6  ,4.9  ,
		     3.9  ,2.9  ,2.3  ,1.9  ,1.6  ,
		     1.4  ,1.2  ,1.03 ,0.91 ,0.81 ,
		     0.73 ,0.66 ,0.60 ,0.51 ,0.43 ,
		     0.37 ,0.29 ,0.23 ,0.18 ,0.15 ,
		     0.12 ,0.11 ,0.088,0.076,0.065,
		     0.052,0.042,0.034,0.028,0.023,
		     0.018,0.013,0.009,0.005,0.001};
  double xpdatal[] ={507.8 ,505.2 ,465.3 ,421.9 ,371.7 ,
		     315.5 ,250.5 ,200.3 ,167.3 ,140.4 ,
		     121.2 ,105.5 , 91.2 ,81.29 ,72.69 ,
		     65.92 ,58.06 ,53.26 ,45.37 ,38.55 ,
		     32.84 ,26.05 ,19.79 ,15.75 ,12.16 ,
		     10.27 , 8.14 , 6.62 ,5.565 ,4.428 ,
		     3.588 ,2.706 ,2.062 ,1.631 ,1.193 ,
		     0.912 ,0.632 ,0.398 ,0.172 ,0.027};
  double xperrorl[]={ 11.9 ,  9.2 ,  7.5 ,  6.3 ,  5.8 ,
		      4.2 ,  3.2 ,  2.3 ,  2.0 ,  1.6 ,
		      1.3 ,  1.1 ,  1.0 , 0.89 , 0.81 ,
		      0.76 , 0.70 , 0.66 , 0.48 , 0.43 ,
		      0.41 , 0.28 , 0.29 , 0.35 , 0.17 ,
		      0.14 , 0.11 , 0.10 ,0.087 ,0.076 ,
		      0.057 ,0.049 ,0.042 ,0.037 ,0.034 ,
		      0.026 ,0.023 ,0.015 ,0.011 ,0.005 };
  double xpdatac[] ={468.7  ,485.4  ,507.2  ,464.0  ,422.9  ,
		     349.1  ,274.1  ,231.5  ,187.5  ,162.3  ,
		     136.6  ,117.3  , 99.1  ,89.21  ,78.25  ,
		     69.26  ,62.25  ,55.28  ,49.14  ,40.11  ,
		     35.64  ,28.94  ,21.99  ,16.51  ,12.69  ,
		     10.41  , 7.86  , 6.37  ,5.060  ,4.080  ,
		     3.123  ,2.141  ,1.472  ,0.952  ,0.935  ,
		     0.485  ,0.372  ,0.140  ,0.027  ,0.011};
  double xperrorc[]={25.3, 22.6, 20.2, 17.6, 16.9,
		     12.7,  9.5,  7.3,  6.9,  5.2,
		     4.2,  3.6,  3.0, 2.65, 2.42, 
		     2.23, 2.06, 1.93, 1.42, 1.29,
		     1.22, 0.83, 0.86, 1.03, 0.50,
		     0.40, 0.32, 0.28,0.241,0.210,
		     0.156,0.126,0.106,0.090,0.085,
		     0.061,0.055,0.031,0.016,0.010};
  double xpdatab[] ={546.1  ,558.5  ,531.9  ,490.8  ,436.5  ,
		     382.8  ,308.6  ,254.8  ,213.2  ,182.1  ,
		     154.5  ,134.3  ,118.6  ,102.4  ,91.92  ,
		     83.63  ,75.06  ,66.58  ,57.31  ,47.80  ,
		     39.19  ,29.54  ,20.69  ,15.36  ,10.65  ,
		     8.06  , 6.28  , 4.69  ,3.490  ,2.935  ,
		     2.041  ,1.534  ,1.111  ,0.736  ,0.510  ,
		     0.330  ,0.188  ,0.089  ,0.017  ,0.003};
  double xperrorb[]={ 14.2, 10.3,  8.0,  6.3,  6.7,
		      4.4,  3.4,  2.6,  2.1,  1.9, 
		      1.7,  1.6,  1.4,  1.3, 1.22,
		      1.14, 1.10, 1.04, 0.74, 0.71,
		      0.67, 0.48, 0.48, 0.55, 0.26, 
		      0.20, 0.16, 0.13,0.110,0.098,
		      0.068,0.058,0.051,0.042,0.037,
		      0.027,0.022,0.012,0.006,0.002};
  vector<double> bins =vector<double>(xpbins,xpbins+41);
  vector<double> data =vector<double>(xpdataa,xpdataa+40);
  vector<double> error=vector<double>(xperrora,xperrora+40);
  _xpa=new_ptr(Histogram(bins,data,error));
  data =vector<double>(xpdatal,xpdatal+40);
  error=vector<double>(xperrorl,xperrorl+40);
  _xpl=new_ptr(Histogram(bins,data,error));
  data =vector<double>(xpdatac,xpdatac+40);
  error=vector<double>(xperrorc,xperrorc+40);
  _xpc=new_ptr(Histogram(bins,data,error));
  data =vector<double>(xpdatab,xpdatab+40);
  error=vector<double>(xperrorb,xperrorb+40);
  _xpb=new_ptr(Histogram(bins,data,error));
  // SLD pions data
  double pipmbins[] = {0.005, 0.008, 0.010, 0.012, 0.014, 
		       0.016, 0.022, 0.027, 0.033, 0.038, 
		       0.044, 0.049, 0.055, 0.060, 0.066, 
		       0.071, 0.077, 0.082, 0.088, 0.099, 
		       0.110, 0.121, 0.143, 0.164, 0.186, 
		       0.208, 0.230, 0.252, 0.274, 0.296, 
		       0.318, 0.351, 0.384, 0.417, 0.450, 
		       0.482, 0.526, 0.570, 0.658, 0.768, 
		       1.000};
  double pipmdataa[]={ 471.8 , 470.4 , 434.6 , 388.8 , 352.7 ,
		       294.8 , 229.6 , 185.0 , 150.6 , 125.6 ,
		       106.5 , 90.40 , 77.38 , 67.39 , 59.40 ,
		       52.57 , 46.76 , 41.70 , 35.26 , 28.89 ,
		       23.88 , 18.69 , 13.85 , 10.16 , 7.812 ,
		       6.076 , 4.674 , 3.632 , 2.886 , 2.292 ,
		       1.749 , 1.275 , 0.921 , 0.680 , 0.499 ,
		       0.338 , 0.226 , 0.130 ,0.0526 ,0.0113 };
  double pipmerrorastat[]={1.3,    1.1,    1.1,    1.0,    0.9,
			   0.5,    0.5,    0.4,    0.4,    0.4,
			   0.4,   0.35,   0.31,   0.29,   0.27,
			   0.25,   0.24,   0.23,   0.15,   0.13,
			   0.12,   0.08,   0.07,   0.06,  0.050,
			   0.044,  0.039,  0.035,  0.031,  0.028,
			   0.021,  0.018,  0.016,  0.014,  0.013,
			   0.010,  0.009,  0.005, 0.0037, 0.0018};
  double pipmerrorasytm[]={9.2,     6.6,     5.0,     4.0,     3.3,
			   2.2,     1.3,     0.9,     0.7,     0.9,
			   1.1,    1.33,    0.91,    0.70,    0.75,
			   0.60,    0.50,    0.43,    0.36,    0.29, 
			   0.25,    0.19,    0.14,    0.11,   0.069,
			   0.061,   0.053,   0.044,   0.037,   0.031,
			   0.034,   0.028,   0.022,   0.018,   0.014,
			   0.010,   0.007,   0.005,  0.0029,  0.0013};
  double pipmdatal[] ={474.0     ,  467.3    ,  418.2    ,  375.5    ,  327.7    ,
		       275.8    ,  216.0    ,  171.2    ,  140.4    ,  116.4    ,
		       99.9    ,   85.4    ,  72.85    ,  64.51    ,  56.82    ,
		       50.84    ,  45.34    ,  40.71    ,  34.60    ,  28.99    ,
		       24.19    ,  18.97    ,  14.52    ,  11.06    ,   8.67    ,
		       6.79    ,  5.341    ,  4.214    ,  3.452    ,  2.727    ,
		       2.138    ,  1.652    ,  1.164    ,  0.874    ,  0.622    ,
		       0.441    ,  0.300    ,  0.178    ,  0.081    ,  0.016};
  double pipmerrorlstat[]={13.9 , 10.5 ,  8.4 ,  6.9 ,  5.7 ,
			   4.2 ,  3.0 ,  2.2 ,  1.9 ,  1.5 ,
			   1.2 ,  1.0 , 0.89 , 0.79 , 0.72 ,
			   0.66 , 0.61 , 0.56 , 0.40 , 0.35 ,
			   0.31 , 0.22 , 0.17 , 0.14 , 0.12 ,
			   0.10 ,0.085 ,0.073 ,0.064 ,0.056 ,
			   0.042 ,0.036 ,0.031 ,0.027 ,0.024 ,
			   0.019 ,0.017 ,0.010 ,0.007 ,0.003};
  double pipmerrorlsytm[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double pipmdatac[] ={425.5     ,  440.5    ,  453.8    ,  409.2    ,  372.8    ,
		       306.4    ,  234.6    ,  197.4    ,  155.8    ,  132.5    ,
		       109.3    ,   92.9    ,  77.56    ,  68.23    ,  60.06    ,
		       51.76    ,  45.28    ,  40.04    ,  33.50    ,  27.45    , 
		       22.92    ,  18.73    ,  13.72    ,  10.18    ,   7.53    ,
		       5.76    ,  4.381    ,  3.358    ,  2.487    ,  1.947    ,
		       1.436    ,  0.817    ,  0.614    ,  0.386    ,  0.429    ,
		       0.206    ,  0.142    ,  0.066    ,  0.003    ,  0.003};
  double pipmerrorcstat[]={26.6 , 23.2 , 20.0 , 17.2 , 14.6 ,
			   10.9 ,  8.0 ,  6.2 ,  6.0 ,  4.3 ,
			   3.5 ,  2.9 , 2.48 , 2.17 , 1.97  ,
			   1.81 , 1.67 , 1.55 , 1.12 , 0.99 ,
			   0.87 , 0.63 , 0.50 , 0.41 , 0.34 ,
			   0.29 ,0.235 ,0.202 ,0.171 ,0.148 ,
			   0.108 ,0.087 ,0.074 ,0.063 ,0.061 ,
			   0.043 ,0.037 ,0.021 ,0.010 ,0.006};
  double pipmerrorcsytm[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double pipmdatab[] ={478.1   ,488.4  ,463.7  ,432.2  ,382.4  ,
		       333.3  ,261.7  ,214.2  ,175.2  ,145.4  ,
		       121.4  ,103.3  ,89.24  ,75.47  ,65.97  ,
		       59.39  ,52.11  ,45.86  ,38.29  ,30.57  ,
		       24.34  ,18.21  ,12.27  , 8.25  , 5.83  ,
		       4.14  ,2.984  ,2.303  ,1.642  ,1.365  ,
		       0.886  ,0.631  ,0.490  ,0.276  ,0.187  ,
		       0.111  ,0.045  ,0.039  ,0.011  ,0.003};
  double pipmerrorbstat[]={15.8,  11.9,   9.5,   7.7,   6.5,
			   4.6,   3.3,   2.7,   2.3,   1.9, 
			   1.7,   1.5,  1.36,  1.21,  1.12,
			   1.04,  0.97,  0.90,  0.65,  0.58,
			   0.51,  0.36,  0.28,  0.22,  0.18,
			   0.15, 0.128, 0.110, 0.094, 0.085,
			   0.063, 0.052, 0.047, 0.038, 0.033,
			   0.025, 0.019, 0.010, 0.005, 0.002};
  double pipmerrorbsytm[]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double pipmerrora[40],pipmerrorl[40],pipmerrorc[40],pipmerrorb[40];
  for(unsigned int ix=0;ix<40;++ix) {
    pipmerrora[ix]=sqrt(sqr(pipmerrorastat[ix])+sqr(pipmerrorasytm[ix]));
    pipmerrorl[ix]=sqrt(sqr(pipmerrorlstat[ix])+sqr(pipmerrorlsytm[ix]));
    pipmerrorc[ix]=sqrt(sqr(pipmerrorcstat[ix])+sqr(pipmerrorcsytm[ix]));
    pipmerrorb[ix]=sqrt(sqr(pipmerrorbstat[ix])+sqr(pipmerrorbsytm[ix]));
  }
  bins =vector<double>(pipmbins,pipmbins+41);
  data =vector<double>(pipmdataa,pipmdataa+40);
  error=vector<double>(pipmerrora,pipmerrora+40);
  _pipma=new_ptr(Histogram(bins,data,error));
  data =vector<double>(pipmdatal,pipmdatal+40);
  error=vector<double>(pipmerrorl,pipmerrorl+40);
  _pipml=new_ptr(Histogram(bins,data,error));
  data =vector<double>(pipmdatac,pipmdatac+40);
  error=vector<double>(pipmerrorc,pipmerrorc+40);
  _pipmc=new_ptr(Histogram(bins,data,error));
  data =vector<double>(pipmdatab,pipmdatab+40);
  error=vector<double>(pipmerrorb,pipmerrorb+40);
  _pipmb=new_ptr(Histogram(bins,data,error));
  // opal pion data
  double pibinso[] = {0.227, 0.239, 0.251, 0.263, 0.276, 
		      0.290, 0.305, 0.320, 0.336, 0.353, 
		      0.371, 0.390, 0.410, 0.431, 0.453, 
		      0.476, 0.500, 0.525, 0.552, 0.580, 
		      0.610, 0.641, 0.673, 0.708, 0.744, 
		      0.782, 0.822, 0.864, 2.02, 2.12, 
		      2.23, 
		      2.34, 2.59, 2.72, 2.86, 3.01, 
		      3.16, 3.32, 3.49, 3.67, 3.86, 
		      4.06, 4.95, 6.05, 7.39, 9.02, 
		      11.02, 13.46, 16.44, 20.08, 29.95, 
		      45.60};
  double pidatao[]={  9.89   ,  9.98   , 10.37   , 10.38   , 10.42   ,
		      10.40   , 10.53   , 10.66   , 10.53   , 10.56   ,
		      10.46   , 10.23   , 10.29   , 10.02   ,  9.83   ,
		      9.62   ,  9.46   ,  9.23   ,  9.05   ,  8.71   ,
		      8.40   ,  8.12   ,  7.87   ,  7.50   ,  7.26   ,
		      6.95   ,  6.56   ,  0.0	 ,  2.374  ,  2.264  ,
		      2.103  ,  1.944  ,  1.672  ,  1.550  ,  1.420  ,
		      1.328  ,  1.221  ,  1.115  ,  1.035  ,  0.955  ,
		      0.879  ,  0.705  ,  0.478  ,  0.319  ,  0.2052 ,
		      0.1246 ,  0.0717 ,  0.0386 ,  0.0206 ,  0.0058 ,
		      0.0006 };
  double pierrorostat[]={0.07  ,0.07  ,0.07  ,0.07  ,0.07  ,
			 0.07  ,0.06  ,0.06  ,0.06  ,0.06  ,
			 0.06  ,0.06  ,0.06  ,0.05  ,0.05  ,
			 0.05  ,0.05  ,0.05  ,0.05  ,0.04  ,
			 0.04  ,0.04  ,0.04  ,0.04  ,0.03  ,
			 0.03  ,0.04  ,0.0   ,0.031 ,0.013 ,
			 0.018 ,0.010 ,0.009 ,0.013 ,0.008 ,
			 0.008 ,0.008 ,0.007 ,0.007 ,0.006 ,
			 0.006 ,0.003 ,0.002 ,0.001 ,0.0012,
			 0.0008,0.0006,0.0005,0.0004,0.0001,
			 0.0001};
  double pierrorosyst[]={0.29   ,0.29   ,0.31   ,0.31   ,0.31   ,
			 0.31   ,0.31   ,0.31   ,0.31   ,0.31   ,
			 0.31   ,0.30   ,0.30   ,0.29   ,0.29   ,
			 0.28   ,0.28   ,0.27   ,0.26   ,0.25   ,
			 0.25   ,0.24   ,0.23   ,0.22   ,0.21   ,
			 0.20     ,0.20     ,0.0,0.209   ,0.079   ,
			 0.036   ,0.027   ,0.024   ,0.024   ,0.022  ,
			 0.019   ,0.017   ,0.015   ,0.015   ,0.014   ,
			 0.012   ,0.009   ,0.007   ,0.005   ,0.0032  ,
			 0.0024  ,0.0016  ,0.0012  ,0.0013  ,0.0004  ,
			 0.0002 };
  double pierroro[51];
  for(unsigned int ix=0;ix<51;++ix)
    {pierroro[ix]=sqrt(sqr(pierrorostat[ix])+sqr(pierrorosyst[ix]));}
  bins =vector<double>(pibinso,pibinso+52);
  data =vector<double>(pidatao,pidatao+51);
  error=vector<double>(pierroro,pierroro+51);
  _pipm=new_ptr(Histogram(bins,data,error)); 
  // SLD kaons
  double Kbins[] = {0.014, 0.016, 0.022, 0.027, 0.033, 
		    0.038, 0.044, 0.049, 0.055, 0.060, 
		    0.066, 0.071, 0.077, 0.082, 0.088, 
		    0.099, 0.110, 0.121, 0.143, 0.164, 
		    0.186, 0.208, 0.230, 0.252, 0.274, 
		    0.296, 0.318, 0.351, 0.384, 0.417, 
		    0.450, 0.482, 0.526, 0.570, 0.658, 
		    0.768, 1.000};
  double Kdataa[]={28.59 ,  21.57 ,  21.62 ,  19.65 ,  18.02 ,
		   17.27 ,  15.78 , 14.664 , 13.535 , 12.599 ,
		   12.036 , 11.349 , 10.207 ,  9.571 ,  8.671 ,
		   7.784 ,  7.237 ,  5.746 ,  3.959 ,  3.473 , 
		   2.739 ,  2.452 ,  1.903 ,  1.574 ,  1.360 ,
		   1.118 ,  0.890 ,  0.683 ,  0.567 ,  0.433 ,
		   0.351 ,  0.264 ,  0.188 ,  0.122 , 0.0485 ,
		   0.0078};
  double Kerrorastat[]={0.64 ,   0.20 ,   0.19 ,   0.18 ,   0.16 ,
			0.17 ,   0.17 ,  0.194 ,  0.189 ,  0.176 ,
			0.165 ,  0.162 ,  0.164 ,  0.160 ,  0.113 ,
			0.114 ,  0.120 ,  0.089 ,  0.102 ,  0.134 ,
			0.047 ,  0.037 ,  0.030 ,  0.027 ,  0.024 ,
			0.022 ,  0.016 ,  0.014 ,  0.013 ,  0.012 ,
			0.011 ,  0.008 ,  0.008 ,  0.005 , 0.0037 ,
			0.0022};
  double Kerrorasytm[]={9.26,   1.57,   0.80,   0.53,   0.44,
			0.43,   0.47,  0.442,  0.503,  0.558,
			0.635,  0.622,  0.603,  0.566,  0.505,
			0.440,  0.395,  0.369,  0.381,  0.532,
			0.419,  0.163,  0.063,  0.036,  0.026,
			0.020,  0.017,  0.016,  0.015,  0.014,
			0.012,  0.010,  0.008,  0.006, 0.0027,
			0.0011};
  double Kdatal[] ={27.05 ,  20.00 ,  19.74 ,  17.52 ,  16.08 ,
		    15.04 ,  13.54 ,  11.87 ,  11.44 ,  10.64 ,
		    10.24 ,   9.67 ,   8.13 ,   7.98 ,   7.00 ,
		    6.36 ,   5.85 ,   4.89 ,   3.41 ,   2.84 ,
		    2.564 ,  2.401 ,  1.973 ,  1.643 ,  1.481 ,
		    1.211 ,  1.001 ,  0.746 ,  0.666 ,  0.559 ,
		    0.426 ,  0.363 ,  0.261 ,  0.183 ,  0.079 ,
		    0.008};
  double Kerrorlstat[]={1.27 , 0.42 , 0.40 , 0.37 , 0.37 ,
			0.34 , 0.34 , 0.34 , 0.33 , 0.30 ,
			0.29 , 0.29 , 0.27 , 0.28 , 0.19 ,
			0.19 , 0.20 , 0.15 , 0.17 , 0.22 ,
			0.082 ,0.067 ,0.054 ,0.048 ,0.044 ,
			0.039 ,0.029 ,0.025 ,0.023 ,0.022 ,
			0.020 ,0.016 ,0.015 ,0.010 ,0.007 ,
			0.004};
  double Kerrorlsytm[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,0.0,0.0};
  double Kdatac[] ={30.92 ,  22.43 ,  22.04 ,  20.82 ,  16.79 ,
		    16.68 ,  16.46 ,  15.81 ,  12.62 ,  12.24 ,
		    11.42 ,  10.95 ,  10.88 ,   9.62 ,   9.84 ,
		    8.08 ,   8.98 ,   6.59 ,   5.50 ,   5.12 ,
		    3.850 ,  3.087 ,  2.074 ,  1.960 ,  1.681 ,
		    1.368 ,  1.043 ,  0.874 ,  0.600 ,  0.408 ,
		    0.408 ,  0.243 ,  0.173 ,  0.064 ,  0.009 ,
		    0.008};
  double Kerrorcstat[]={3.86 , 1.37 , 1.27 , 1.17 , 1.15 ,
			1.15 , 1.06 , 1.08 , 0.99 , 0.92 ,
			0.87 , 0.85 , 0.84 , 0.81 , 0.59 ,
			0.58 , 0.63 , 0.45 , 0.51 , 0.68 ,
			0.245 ,0.190 ,0.145 ,0.132 ,0.119 ,
			0.104 ,0.076 ,0.068 ,0.058 ,0.050 ,
			0.050 ,0.037 ,0.034 ,0.020 ,0.011 ,
			0.008};
  double Kerrorcsytm[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,0.0,0.0};
  double Kdatab[] ={30.21 ,  23.06 ,  22.89 ,  21.64 ,  21.36 ,
		    21.36 ,  19.90 ,  18.91 ,  18.46 ,  17.43 ,
		    16.92 ,  15.62 ,  15.11 ,  13.18 ,  12.43 ,
		    11.56 ,   9.96 ,   7.17 ,   4.58 ,   4.20 ,
		    2.541 ,  2.009 ,  1.627 ,  1.116 ,  0.830 ,
		    0.640 ,  0.452 ,  0.337 ,  0.245 ,  0.149 ,
		    0.108 ,  0.057 ,  0.061 ,  0.012 ,  0.002 ,
		    0.001};
  double Kerrorbstat[]={1.99 , 0.62 , 0.60 , 0.55 , 0.53 ,
			0.56 , 0.57 , 0.60 , 0.58 , 0.54 ,
			0.53 , 0.52 , 0.52 , 0.50 , 0.36 ,
			0.37 , 0.38 , 0.27 , 0.29 , 0.36 ,
			0.126 ,0.096 ,0.078 ,0.062 ,0.053 ,
			0.045 ,0.032 ,0.028 ,0.024 ,0.020 ,
			0.018 ,0.012 ,0.013 ,0.005 ,0.003 ,
			0.001};
  double Kerrorbsytm[]={0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,
			0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,
			0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,
			0.0 ,0.0 ,0.0 ,0.0 ,0.0 ,0.0};
  double Kerrora[36],Kerrorl[36],Kerrorc[36],Kerrorb[36];
  for(unsigned int ix=0;ix<36;++ix)
    {
      Kerrora[ix]=sqrt(sqr(Kerrorastat[ix])+sqr(Kerrorasytm[ix]));
      Kerrorl[ix]=sqrt(sqr(Kerrorlstat[ix])+sqr(Kerrorlsytm[ix]));
      Kerrorc[ix]=sqrt(sqr(Kerrorcstat[ix])+sqr(Kerrorcsytm[ix]));
      Kerrorb[ix]=sqrt(sqr(Kerrorbstat[ix])+sqr(Kerrorbsytm[ix]));
    }
  bins =vector<double>(Kbins,Kbins+37);
  data =vector<double>(Kdataa,Kdataa+36);
  error=vector<double>(Kerrora,Kerrora+36);
  _kpma=new_ptr(Histogram(bins,data,error));
  data =vector<double>(Kdatal,Kdatal+36);
  error=vector<double>(Kerrorl,Kerrorl+36);
  _kpml=new_ptr(Histogram(bins,data,error));
  data =vector<double>(Kdatac,Kdatac+36);
  error=vector<double>(Kerrorc,Kerrorc+36);
  _kpmc=new_ptr(Histogram(bins,data,error));
  data =vector<double>(Kdatab,Kdatab+36);
  error=vector<double>(Kerrorb,Kerrorb+36);
  _kpmb=new_ptr(Histogram(bins,data,error));
    // OPAL kaons
  double Kbinso[] = {0.271, 0.281, 0.292, 0.304, 0.317, 
		     0.331, 0.346, 0.362, 0.379, 0.397, 
		     0.416, 0.436, 0.457, 0.480, 0.504, 
		     0.528, 0.555, 0.583, 0.612, 0.643, 
		     0.675, 0.709, 0.745, 0.783, 4.05, 
		     4.95, 
		     6.05, 7.38, 9.02, 11.01, 13.45, 
		     16.43, 20.06, 29.93, 45.60};
  double Kdatao[]={0.363 ,0.373 ,0.367 ,0.374 ,0.375 ,
		    0.410 ,0.431 ,0.418 ,0.456 ,0.499 ,
		    0.514 ,0.486 ,0.522 ,0.541 ,0.539 ,
		    0.557 ,0.587 ,0.590 ,0.586 ,0.591 ,
		    0.614 ,0.597 ,0.613 ,0.0   ,0.181 ,
		    0.138 ,0.103 ,0.0767,0.0536,0.0349,
		    0.0220,0.0127,0.0042,0.0004};
  double Kerrorostat[]={0.030 ,0.027 ,0.024 ,0.022 ,0.021 ,
			 0.020 ,0.020 ,0.018 ,0.018 ,0.018 ,
			 0.017 ,0.011 ,0.011 ,0.011 ,0.011 ,
			 0.011 ,0.011 ,0.010 ,0.010 ,0.010 ,
			 0.009 ,0.009 ,0.009 ,0.0   ,0.004 ,
			 0.003 ,0.001 ,0.0010,0.0006,0.0005,
			 0.0003,0.0003,0.0001,0.0001};
  double Kerrorosyst[]={0.028 ,0.028 ,0.025 ,0.023 ,0.019 ,
			 0.019 ,0.019 ,0.020 ,0.016 ,0.017 ,
			 0.016 ,0.011 ,0.011 ,0.012 ,0.018 ,
			 0.020 ,0.014 ,0.048 ,0.016 ,0.044 ,
			 0.024 ,0.025 ,0.031 , 0.0  ,0.015 ,
			 0.009 ,0.006 ,0.0042,0.0029,0.0018,
			 0.0012,0.0007,0.0003,0.0001};
  double Kerroro[34];
  for(unsigned int ix=0;ix<34;++ix)
    {Kerroro[ix]=sqrt(sqr(Kerrorostat[ix])+sqr(Kerrorosyst[ix]));}
  bins =vector<double>(Kbinso,Kbinso+35);
  data =vector<double>(Kdatao,Kdatao+34);
  error=vector<double>(Kerroro,Kerroro+34);
  _kpm=new_ptr(Histogram(bins,data,error));
  // SLD proton data
  double pbins[] = {0.014, 0.016, 0.022, 0.027, 0.033, 
		    0.038, 0.044, 0.049, 0.055, 0.060, 
		    0.066, 0.071, 0.077, 0.082, 0.088, 
		    0.099, 0.110, 0.121, 0.143, 0.164, 
		    0.186, 0.208, 0.230, 0.252, 0.274, 
		    0.296, 0.318, 0.351, 0.384, 0.417, 
		    0.450, 0.482, 0.526, 0.570, 0.658, 
		    0.768, 1.000};
  double pdataa[]={14.51 ,  17.32 ,  13.75 ,  11.12 ,  10.75 ,  9.048 ,  7.669 ,  7.410 ,  6.587 ,  5.788 ,  5.344 ,  4.987 ,  4.278 ,  4.117 ,  3.633 ,  3.036 ,  2.568 ,  2.165 ,  1.931 ,  1.603 ,  0.871 ,  0.912 ,  0.775 ,  0.639 ,  0.511 ,  0.419 ,  0.358 ,  0.254 ,  0.173 ,  0.141 , 0.0950 , 0.0688 , 0.0470 , 0.0241 , 0.0093 , 0.0015};
  double perrorastat[]={ 0.52  ,   0.27  ,   0.29  ,   0.17  ,   0.14  ,  0.123  ,  0.117  ,  0.113  ,  0.109  ,  0.105  ,  0.100  ,  0.104  ,  0.100  ,  0.101  ,  0.072  ,  0.076  ,  0.081  ,  0.069  ,  0.096  ,  0.133  ,  0.045  ,  0.030  ,  0.025  ,  0.022  ,  0.019  ,  0.016  ,  0.011  ,  0.009  ,  0.008  ,  0.007  , 0.0055  , 0.0039  , 0.0032  , 0.0017  , 0.0010  , 0.0003};
  double perrorasytm[]={5.08,  2.58,  2.50,  1.24,  0.47, 0.350, 0.298, 0.294, 0.259, 0.238, 0.228, 0.229, 0.242, 0.253, 0.269, 0.300, 0.357, 0.398, 0.452, 0.594, 0.255, 0.179, 0.062, 0.044, 0.033, 0.024, 0.018, 0.012, 0.008, 0.005,0.0036,0.0027,0.0018,0.0012,0.0006,0.0001};
  double pdatal[] ={13.98 , 17.63 , 13.42 , 10.57 ,  9.98 ,  8.37 ,  7.33 ,  7.79 ,  6.62 ,  5.88 ,  5.39 ,  5.22 ,  4.42 ,  4.44 ,  3.65 ,  3.11 ,  2.73 ,  2.15 ,  1.83 ,  1.84 , 0.905 , 1.065 , 0.822 , 0.762 , 0.628 , 0.486 , 0.446 , 0.306 , 0.230 , 0.197 , 0.145 , 0.108 , 0.070 , 0.036 , 0.013 , 0.003};
  double perrorlstat[]={0.99 , 0.58 , 0.60 , 0.36 , 0.31 , 0.26 , 0.24 , 0.23 , 0.22 , 0.20 , 0.19 , 0.19 , 0.18 , 0.19 , 0.13 , 0.13 , 0.15 , 0.12 , 0.16 , 0.24 ,0.078 ,0.054 ,0.044 ,0.038 ,0.033 ,0.029 ,0.022 ,0.018 ,0.015 ,0.013 ,0.011 ,0.008 ,0.006 ,0.003 ,0.002 ,0.001};
  double perrorlsytm[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double pdatac[] ={13.28, 15.22, 13.32,  9.60, 11.64, 10.07,  8.10,  6.09,  6.54,  6.36,  4.62,  4.43,  4.08,  3.67,  4.07,  2.98,  2.30,  2.39,  1.72,  0.31, 0.561, 0.978, 0.907, 0.652, 0.572, 0.494, 0.454, 0.314, 0.170, 0.103, 0.064, 0.015, 0.044, 0.007, 0.015, 0.001};
  double perrorcstat[]={ 2.94 ,  1.89 ,  1.86 ,  1.16 ,  1.01 ,  0.87 ,  0.76 ,  0.72 ,  0.67 ,  0.63 ,  0.59 ,  0.58 ,  0.55 ,  0.55 ,  0.42 ,  0.41 ,  0.43 ,  0.36 ,  0.50 ,  0.71 , 0.235 , 0.163 , 0.136 , 0.116 , 0.101 , 0.089 , 0.069 , 0.054 , 0.041 , 0.033 , 0.028 , 0.018 , 0.017 , 0.007 , 0.005 , 0.002};
  double perrorcsytm[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double pdatab[] ={13.79 , 17.93 , 16.41 , 12.11 , 10.32 ,  9.52 ,  7.72 ,  6.86 ,  6.19 ,  4.96 ,  4.82 ,  4.57 ,  4.07 ,  3.82 ,  3.29 ,  2.68 ,  2.24 ,  1.84 ,  1.91 ,  1.25 , 0.867 , 0.739 , 0.645 , 0.392 , 0.252 , 0.266 , 0.146 , 0.102 , 0.020 , 0.034 , 0.004 , 0.016 , 0.003 , 0.004 , 0.001 , 0.000};
  double perrorbstat[]={1.49 , 0.78 , 0.93 , 0.52 , 0.42 ,
			0.40 , 0.36 , 0.32 , 0.32 , 0.29 ,
			0.29 , 0.29 , 0.29 , 0.29 , 0.21 ,
			0.22 , 0.23 , 0.19 , 0.27 , 0.36 ,
			0.121 ,0.072 ,0.060 ,0.047 ,0.038 ,
			0.035 ,0.022 ,0.018 ,0.012 ,0.011 ,
			0.007 ,0.006 ,0.003 ,0.002 ,0.001 ,
			0.000};
  double perrorbsytm[]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,0.0,0.0};
  double perrora[36],perrorl[36],perrorc[36],perrorb[36];
  for(unsigned int ix=0;ix<36;++ix)
    {
      perrora[ix]=sqrt(sqr(perrorastat[ix])+sqr(perrorasytm[ix]));
      perrorl[ix]=sqrt(sqr(perrorlstat[ix])+sqr(perrorlsytm[ix]));
      perrorc[ix]=sqrt(sqr(perrorcstat[ix])+sqr(perrorcsytm[ix]));
      perrorb[ix]=sqrt(sqr(perrorbstat[ix])+sqr(perrorbsytm[ix]));
    }
  bins =vector<double>(pbins,pbins+37);
  data =vector<double>(pdataa,pdataa+36);
  error=vector<double>(perrora,perrora+36);
  _ppma=new_ptr(Histogram(bins,data,error));
  data =vector<double>(pdatal,pdatal+36);
  error=vector<double>(perrorl,perrorl+36);
  _ppml=new_ptr(Histogram(bins,data,error));
  data =vector<double>(pdatac,pdatac+36);
  error=vector<double>(perrorc,perrorc+36);
  _ppmc=new_ptr(Histogram(bins,data,error));
  data =vector<double>(pdatab,pdatab+36);
  error=vector<double>(perrorb,perrorb+36);
  _ppmb=new_ptr(Histogram(bins,data,error));
  // OPAL protons
  double pbinso[] = {0.406, 0.421, 0.438, 0.456, 0.475, 
		     0.495, 0.517, 0.541, 0.565, 0.592, 
		     0.620, 0.650, 0.681, 0.714, 0.750, 
		     0.787, 0.826, 0.867, 0.911, 0.957, 
		     1.005, 1.056, 1.109, 1.166, 1.225, 
		     1.287, 1.353, 1.422, 4.05, 4.95, 
		     6.04, 
		     7.38, 9.01, 11.01, 13.44, 16.42, 
		     20.05, 29.90, 45.60};
  double pdatao[]={0.169      ,0.174       ,0.171        ,0.185      ,
		   0.181      ,0.202       ,0.215      ,0.228      ,
		   0.230      ,0.230       ,0.236       ,0.239       ,
		   0.246       ,0.248       ,0.248      ,0.254      ,
		   0.260      ,0.260      ,0.257      ,0.263        ,
		   0.257      ,0.256        ,0.257        ,0.248        ,
		   0.241        ,0.250        ,0.245        ,0.0	   ,
		   0.0606      ,0.0517      ,0.0352      ,0.0214      ,
		   0.0152      ,0.0093      ,0.0045      ,0.0015      ,
		   0.00056    ,0.000015};
  double perrorostat[]={0.010         ,0.009         ,0.009         ,0.008         ,
			0.008         ,0.008         ,0.008         ,0.008         ,
			0.008         ,0.008         ,0.007         ,0.007         ,
			0.007         ,0.007         ,0.004         ,0.004         ,
			0.004         ,0.004         ,0.004         ,0.004         ,
			0.004         ,0.004         ,0.004         ,0.004         ,
			0.004         ,0.004         ,0.004         ,0.0	   ,
			0.0036        ,0.0029        ,0.0011        ,0.0011        ,
			0.0006        ,0.0004        ,0.0003        ,0.0002        ,
			0.00005       ,0.000004};
  double perrorosyst[]={0.013   ,0.015   ,0.012   ,0.015   ,0.006   ,
			0.010   ,0.012   ,0.014   ,0.010   ,0.007   ,
			0.012   ,0.008   ,0.006   ,0.009   ,0.009   ,
			0.008   ,0.007   ,0.007   ,0.006   ,0.007   ,
			0.021   ,0.027   ,0.029   ,0.027   ,0.017   ,
			0.016   ,0.021   ,   0.0  ,0.0129  ,0.0062  ,
			0.0031  ,0.0023  ,0.0021  ,0.0013  ,0.0006  ,
			0.0005  ,0.00014 ,0.000026};
  double perroro[38];
  for(unsigned int ix=0;ix<38;++ix)
    perroro[ix]=sqrt(sqr(perrorostat[ix])+sqr(perrorosyst[ix]));
  bins =vector<double>(pbinso,pbinso+39);
  data =vector<double>(pdatao,pdatao+38);
  error=vector<double>(perroro,perroro+38);
  _ppm=new_ptr(Histogram(bins,data,error)); 
  // OPAL light quarks
  double udsbinso[] = {0.00, 0.01, 0.02, 0.03, 0.04, 
		       0.05, 0.06, 0.07, 0.08, 0.09, 
		       0.10, 0.12, 0.14, 0.16, 0.18, 
		       0.20, 0.25, 0.30, 0.40, 0.50, 
		       0.60, 0.80, 1.00};
  double udsdatao[]={388.  ,390. ,241.  ,176.  ,122.6 ,
		     95.7   ,79.3   ,65.0   ,53.3   ,43.3   ,
		     35.1   ,27.7   ,21.2   ,17.1   ,13.3   ,
		     9.86   ,6.30   ,3.42   ,1.50   ,0.668  ,
		     0.241   ,  0.031};
  double udserrorostat[]={5.  ,5.  ,4.  ,3.,2.7   ,2.2 ,
			  1.9 ,1.6 ,1.6 ,1.5 ,0.9 ,0.7 ,
			  0.7 ,0.6 ,0.6 ,0.26,0.19,0.09,
			  0.05,0.033   ,0.008   ,0.007};
  double udserrorosyst[]={9. ,10.,7. ,5. ,3.9,2.9,
			  2.3,1.7,1.3,1.0,0.7,0.4,
			  0.4,0.3,0.3,0.30   ,0.25   ,
			  0.17   ,0.10   ,0.048  ,0.024  ,0.007};
  double udserroro[22];
  for(unsigned int ix=0;ix<22;++ix)
    {udserroro[ix]=sqrt(sqr(udserrorostat[ix])+sqr(udserrorosyst[ix]));}
  bins =vector<double>(udsbinso,udsbinso+23);
  data =vector<double>(udsdatao,udsdatao+22);
  error=vector<double>(udserroro,udserroro+22);
  _udsxp=new_ptr(Histogram(bins,data,error));
  double udsxibinso[] = {0.0, 0.2, 0.4, 0.6, 0.8, 
			 1.0, 1.2, 1.4, 1.6, 1.8, 
			 2.0, 2.2, 2.4, 2.6, 2.8, 
			 3.0, 3.2, 3.4, 3.6, 3.8, 
			 4.0, 4.2, 4.4, 4.6, 4.8, 
			 5.0, 5.2, 5.4, 5.6, 5.8};
  double udsxidatao[]={0.024    , 0.114    , 0.277    , 0.529    ,  0.86      ,
		       1.31      ,  1.76      ,  2.22      ,  2.70      ,  3.06      ,
		       3.76      ,  4.03      ,  4.48      ,  5.12      ,  5.22      ,
		       5.26      ,  6.24      ,  6.02      ,  5.89      ,  6.04      ,
		       5.85      ,  5.58      ,  5.15      ,  4.21      ,  3.99      ,
		       2.94      ,  2.14      ,  1.93      ,  1.43};
  double udsxierrorostat[]={0.006    ,0.003    ,0.009    ,0.016    ,0.02     ,
			    0.03     ,0.05     ,0.06     ,0.07     ,0.08     ,
			    0.09     ,0.10     ,0.10     ,0.11     ,0.12     ,
			    0.13     ,0.12     ,0.12     ,0.13     ,0.12     ,
			    0.13     ,0.11     ,0.11     ,0.12     ,0.10     ,
			    0.10     ,0.10     ,0.08     ,0.09};
  double udsxierrorosyst[]={0.006,0.011,0.025,0.032,0.05 ,0.06 ,
			    0.07 ,0.06 ,0.06 ,0.09 ,0.11 ,0.13 ,
			    0.18 ,0.16 ,0.17 ,0.19 ,0.21 ,0.20 ,
			    0.26 ,0.20 ,0.20 ,0.14 ,0.09 ,0.24 ,
			    0.14 ,0.15 ,0.12 ,0.13 ,0.23 };
  double udsxierroro[22];
  for(unsigned int ix=0;ix<22;++ix)
    {udsxierroro[ix]=sqrt(sqr(udsxierrorostat[ix])+sqr(udsxierrorosyst[ix]));}
  bins =vector<double>(udsxibinso,udsxibinso+23);
  data =vector<double>(udsxidatao,udsxidatao+22);
  error=vector<double>(udsxierroro,udsxierroro+22);
  _udsxip=new_ptr(Histogram(bins,data,error));
  // lambdas
  double Lbinso[] = {0.012, 0.014, 0.016, 0.018, 0.020, 
		    0.025, 0.030, 0.035, 0.040, 0.050, 
		    0.060, 0.080, 0.100, 0.120, 0.140, 
		    0.160, 0.180, 0.200, 0.250, 0.300, 
		    0.350, 0.400, 0.500, 0.600, 0.700, 
		    0.900};
  double Ldatao[]={2.97    ,3.43    ,3.74    ,3.70    ,3.69    ,3.68    ,
		   3.70    ,3.41    ,3.18    ,2.66    ,2.04    ,1.52    ,
		   1.19    ,0.956   ,0.771   ,0.630   ,0.528   ,0.408   ,
		   0.269   ,0.182   ,0.129   ,0.078   ,0.035   ,0.0118  ,0.0026};
  double Lerrorostat[]={0.35    ,0.30    ,0.29    ,0.21    ,0.18    ,0.16    ,
			0.15    ,0.14    ,0.11    ,0.09    ,0.06    ,0.04    ,
			0.03    ,0.023   ,0.018   ,0.015   ,0.013   ,0.010   ,
			0.008   ,0.007   ,0.006   ,0.005   ,0.003   ,0.0019  ,0.0012};
  double Lerrorosyst[]={0.0,0.0,0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,0.0,0.0,
			0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double Lerroro[25];
  for(unsigned int ix=0;ix<25;++ix) 
    Lerroro[ix]=sqrt(sqr(Lerrorostat[ix])+sqr(Lerrorosyst[ix]));
  bins =vector<double>(Lbinso,Lbinso+26);
  data =vector<double>(Ldatao,Ldatao+25);
  error=vector<double>(Lerroro,Lerroro+25);
  _lpm=new_ptr(Histogram(bins,data,error));
  // K*+/- CERN-PPE-96-186
  double Kstarpmbin[]={0.03    ,0.06    ,0.09    ,0.12    ,0.15    ,
		       0.18    ,0.22    ,0.26    ,0.32    ,0.44    ,1.00};
  double Kstarpmdata[]={5.17   ,3.43    ,2.09    ,2.01    ,1.54   ,
			1.16   ,0.71    ,0.59    ,0.38    ,0.06};
  double Kstarpmstat[]={0.53    ,0.29    ,0.20    ,0.16    ,0.15    ,
			0.12    ,0.09    ,0.06    ,0.04    ,0.01};
  double Kstarpmsyst[]={0.54,0.58,0.22,0.23,0.19,0.22,0.07,0.07,0.03,0.01};
  double Kstarpmerror[10];
  for(unsigned int ix=0;ix<10;++ix) 
    Kstarpmerror[ix] = sqrt(sqr(Kstarpmstat[ix])+sqr(Kstarpmsyst[ix]));
  bins =vector<double>(Kstarpmbin,Kstarpmbin+11);
  data =vector<double>(Kstarpmdata,Kstarpmdata+10);
  error=vector<double>(Kstarpmerror,Kstarpmerror+10);
  _xpKstarplus = new_ptr(Histogram(bins,data,error));
  // xi- ALEXANDER 96   ZP C73,569
  double ximbinsA[]={0.035  ,0.050  ,0.060  ,0.080  ,0.100  ,
		     0.150  ,0.200  ,0.300  ,0.400  ,0.500};
  double ximdataA[]={0.341   ,0.254   ,0.180   ,0.114   ,0.0737  ,
		     0.0438  ,0.0192  ,0.0120  ,0.0084};
  double ximstatA[]={0.018   ,0.014   ,0.007   ,0.005   ,0.0024  ,
		     0.0021  ,0.0015  ,0.0012  ,0.0014};
  double ximsystA[]={0.032   ,0.021   ,0.012   ,0.008   ,0.0048  ,
		     0.0031  ,0.0014  ,0.0012  ,0.0015};
  double ximerrorA[9];
  for(unsigned int ix=0;ix<9;++ix) 
    ximerrorA[ix] = sqrt(sqr(ximstatA[ix])+sqr(ximsystA[ix]));
  bins =vector<double>(ximbinsA,ximbinsA+10);
  data =vector<double>(ximdataA,ximdataA+9);
  error=vector<double>(ximerrorA,ximerrorA+9);
  _xpXiminus = new_ptr(Histogram(bins,data,error));
  double ximbinsB[]={0.6  ,0.8  ,1.0  ,1.2  ,1.4  ,
		     1.6  ,1.8  ,2.0  ,2.2  ,2.4  ,
		     2.6  ,2.8  ,3.0  ,3.2  ,3.4  ,3.6};
  double ximdataB[]={0.0026  ,0.0044  ,0.0044  ,0.0050  ,0.0071  ,
		     0.0075  ,0.0088  ,0.0086  ,0.0094  ,0.0103  ,
		     0.0121  ,0.0102  ,0.0100  ,0.0087  ,0.0088};
  double ximstatB[]={0.0005  ,0.0004  ,0.0004  ,0.0005  ,0.0004  ,
		     0.0004  ,0.0004  ,0.0004  ,0.0004  ,0.0005  ,
		     0.0006  ,0.0006  ,0.0007  ,0.0009  ,0.0011};
  double ximsystB[]={0.0017  ,0.0011  ,0.0005  ,0.0005  ,0.0006  ,
		     0.0006  ,0.0007  ,0.0006  ,0.0007  ,0.0008  ,
		     0.0009  ,0.0008  ,0.0009  ,0.0010  ,0.0013};
  double ximerrorB[15];
  for(unsigned int ix=0;ix<15;++ix) 
    ximerrorB[ix] = sqrt(sqr(ximstatB[ix])+sqr(ximsystB[ix]));
  bins =vector<double>(ximbinsB,ximbinsB+16);
  data =vector<double>(ximdataB,ximdataB+15);
  error=vector<double>(ximerrorB,ximerrorB+15);
  _xiXiminus = new_ptr(Histogram(bins,data,error));
  // sigma*+ALEXANDER 96   ZP C73,569
  double sigmapbinsA[]={0.04  ,0.07  ,0.10  ,0.15  ,0.20  ,0.30  ,0.50};
  double sigmapdataA[]={0.280   ,0.118   ,0.0619  ,0.0403  ,0.0248  ,0.0076};
  double sigmapstatA[]={0.020   ,0.012   ,0.0060  ,0.0048  ,0.0027  ,0.0013};
  double sigmapsystA[]={0.021   ,0.009   ,0.0047  ,0.0031  ,0.0019  ,0.0006};
  double sigmaperrorA[6];
  for(unsigned int ix=0;ix<6;++ix) 
    sigmaperrorA[ix] = sqrt(sqr(sigmapstatA[ix])+sqr(sigmapsystA[ix]));
  bins =vector<double>(sigmapbinsA,sigmapbinsA+7);
  data =vector<double>(sigmapdataA,sigmapdataA+6);
  error=vector<double>(sigmaperrorA,sigmaperrorA+6);
  _xpSigmaplus = new_ptr(Histogram(bins,data,error));
  double sigmapbinsB[]={0.70  ,1.21  ,1.62  ,1.92  ,2.35  ,2.76  ,3.65};
  double sigmapdataB[]={0.0030  ,0.0060  ,0.0068  ,0.0072  ,0.0086  ,0.0095};
  double sigmapstatB[]={0.0005  ,0.0007  ,0.0008  ,0.0007  ,0.0008  ,0.0007};
  double sigmapsystB[]={0.0002  ,0.0005  ,0.0005  ,0.0005  ,0.0007  ,0.0007};
  double sigmaperrorB[6];
  for(unsigned int ix=0;ix<6;++ix) 
    sigmaperrorB[ix] = sqrt(sqr(sigmapstatB[ix])+sqr(sigmapsystB[ix]));
  bins =vector<double>(sigmapbinsB,sigmapbinsB+7);
  data =vector<double>(sigmapdataB,sigmapdataB+6);
  error=vector<double>(sigmaperrorB,sigmaperrorB+6);
  _xiSigmaplus = new_ptr(Histogram(bins,data,error));
  // sigma*+ ALEXANDER 96   ZP C73,569
  double sigmambinsA[]={0.04  ,0.07  ,0.10  ,0.15  ,0.20  ,0.30  ,0.50};
  double sigmamdataA[]={0.291   ,0.116   ,0.0646        ,0.0414  ,0.0235  ,0.0062};
  double sigmamstatA[]={0.021   ,0.013   ,0.0071  ,0.0061  ,0.0040  ,0.0020};
  double sigmamsystA[]={0.022 ,0.009 ,0.000       ,0.0032,0.0018,0.0005};
  double sigmamerrorA[6];
  for(unsigned int ix=0;ix<6;++ix) 
    sigmamerrorA[ix] = sqrt(sqr(sigmamstatA[ix])+sqr(sigmamsystA[ix]));
  bins =vector<double>(sigmambinsA,sigmambinsA+7);
  data =vector<double>(sigmamdataA,sigmamdataA+6);
  error=vector<double>(sigmamerrorA,sigmamerrorA+6);
  _xpSigmaminus = new_ptr(Histogram(bins,data,error));
  double sigmambinsB[]={0.70  ,1.21  ,1.62  ,1.92  ,2.35  ,2.76  ,3.65};
  double sigmamdataB[]={0.0024  ,0.0057  ,0.0070  ,0.0075  ,0.0085  ,0.0098};
  double sigmamstatB[]={0.0008  ,0.0010  ,0.0010  ,0.0008  ,0.0009  ,0.0007};
  double sigmamsystB[]={0.0002,0.0004,0.0005,0.0006,0.0006,0.0007};
  double sigmamerrorB[6];
  for(unsigned int ix=0;ix<6;++ix) 
    sigmamerrorB[ix] = sqrt(sqr(sigmamstatB[ix])+sqr(sigmamsystB[ix]));
  bins =vector<double>(sigmambinsB,sigmambinsB+7);
  data =vector<double>(sigmamdataB,sigmamdataB+6);
  error=vector<double>(sigmamerrorB,sigmamerrorB+6);
  _xiSigmaminus = new_ptr(Histogram(bins,data,error));
  // xi*0 ALEXANDER 96   ZP C73,569
  double xi0binsA[]={0.04  ,0.07  ,0.10  ,0.15  ,0.20  ,0.30};
  double xi0dataA[]={0.064   ,0.0387  ,0.0239  ,0.0144  ,0.0049};
  double xi0statA[]={0.014   ,0.0054  ,0.0028  ,0.0025  ,0.0015};
  double xi0systA[]={0.009 ,0.0042,0.0024,0.0017,0.0006};
  double xi0errorA[5];
  for(unsigned int ix=0;ix<5;++ix) 
    xi0errorA[ix] = sqrt(sqr(xi0statA[ix])+sqr(xi0systA[ix]));
  bins =vector<double>(xi0binsA,xi0binsA+6);
  data =vector<double>(xi0dataA,xi0dataA+5);
  error=vector<double>(xi0errorA,xi0errorA+5);
  _xpXi0 = new_ptr(Histogram(bins,data,error));
  double xi0binsB[]={1.21  ,1.62  ,1.92  ,2.36  ,2.79  ,3.83};
  double xi0dataB[]={0.0012  ,0.0024  ,0.0027  ,0.0027  ,0.0019};
  double xi0statB[]={0.0004  ,0.0004  ,0.0003  ,0.0004  ,0.0004};
  double xi0systB[]={0.0002,0.0003,0.0003,0.0003,0.0003};
  double xi0errorB[5];
  for(unsigned int ix=0;ix<5;++ix) 
    xi0errorB[ix] = sqrt(sqr(xi0statB[ix])+sqr(xi0systB[ix]));
  bins =vector<double>(xi0binsB,xi0binsB+6);
  data =vector<double>(xi0dataB,xi0dataB+5);
  error=vector<double>(xi0errorB,xi0errorB+5);
  _xiXi0 = new_ptr(Histogram(bins,data,error));
  double lambda1520binsA[]={0.035  ,0.039  ,0.045  ,0.052  ,0.060  ,0.100  ,0.300  ,0.500};
  double lambda1520dataA[]={0.447  ,0.139  ,0.150  ,0.163  ,0.000  ,0.033  ,0.008};
  double lambda1520statA[]={0.070    ,0.036   ,0.032   ,0.052   ,0.000   ,0.007   ,0.003};
  double lambda1520systA[]={0.076,0.018,0.018,0.022,0.000,0.006,0.002};
  double lambda1520errorA[7];
  for(unsigned int ix=0;ix<7;++ix) 
    lambda1520errorA[ix] = sqrt(sqr(lambda1520statA[ix])+sqr(lambda1520systA[ix]));
  bins =vector<double>(lambda1520binsA,lambda1520binsA+8);
  data =vector<double>(lambda1520dataA,lambda1520dataA+7);
  error=vector<double>(lambda1520errorA,lambda1520errorA+7);
  _xpLambda1520 = new_ptr(Histogram(bins,data,error));
  double lambda1520binsB[]={0.70  ,1.21  ,2.36  ,3.00  ,3.22  ,3.50  ,3.90,4.54};
  double lambda1520dataB[]={0.0032  ,0.0060  ,0.0000  ,0.0058  ,0.0037  ,0.0020  ,0.0026};
  double lambda1520statB[]={0.0010  ,0.0012  ,0.0000  ,0.0018  ,0.0008  ,0.0005  ,0.0004};
  double lambda1520systB[]={0.0007,0.0010,0.0000,0.0008,0.0005,0.0003,0.0004};
  double lambda1520errorB[7];
  for(unsigned int ix=0;ix<7;++ix) 
    lambda1520errorB[ix] = sqrt(sqr(lambda1520statB[ix])+sqr(lambda1520systB[ix]));
  bins =vector<double>(lambda1520binsB,lambda1520binsB+8);
  data =vector<double>(lambda1520dataB,lambda1520dataB+7);
  error=vector<double>(lambda1520errorB,lambda1520errorB+7);
  _xiLambda1520 = new_ptr(Histogram(bins,data,error));
  // delta++
  double deltabins []={0.05,0.075,0.1,0.15,0.2,0.3,1.0};
  double deltadata []={1.9,2.8,0.38,0.18,0.073 ,0.006};
  double deltaerror[]={0.7,0.8,0.09,0.10,0.043,0.0029};
  bins =vector<double>(deltabins,deltabins+7);
  data =vector<double>(deltadata,deltadata+6);
  error=vector<double>(deltaerror,deltaerror+6);
  _xeDelta = new_ptr(Histogram(bins,data,error));
  // f_0(980)
  double f980bins []={0.00  ,0.06  ,0.12  ,0.14  ,0.16  ,0.20  ,0.25  ,0.35  ,0.50  ,1.00};
  double f980data []={1.04    ,0.57    ,0.30    ,0.20    ,0.21    ,0.13    ,0.085   ,0.046   ,0.0079};
  double f980error[]={0.09    ,0.05    ,0.06    ,0.05    ,0.03    ,0.02    ,0.011   ,0.005   ,0.0009};
  bins =vector<double>(f980bins ,f980bins +10);
  data =vector<double>(f980data ,f980data + 9);
  error=vector<double>(f980error,f980error+ 9);
  _xpf980 = new_ptr(Histogram(bins,data,error));
  // f_2
  double f2bins []={0.00  ,0.06  ,0.12  ,0.14  ,0.16  ,0.20  ,0.25  ,0.35  ,0.50  ,1.00  };
  double f2data []={1.00   ,0.69   ,0.41   ,0.25   ,0.27   ,0.22   ,0.091  ,0.035  ,0.008};
  double f2error[]={0.14   ,0.08   ,0.09   ,0.08   ,0.04   ,0.03   ,0.016  ,0.008  ,0.001};
  bins =vector<double>(f2bins ,f2bins +10);
  data =vector<double>(f2data ,f2data + 9);
  error=vector<double>(f2error,f2error+ 9);
  _xpf2 = new_ptr(Histogram(bins,data,error));
// phi
  double phibins []={0.00  ,0.06  ,0.12  ,0.14  ,0.16  ,0.20  ,0.25  ,0.35  ,0.50  ,1.00};
  double phidata []={0.464  ,0.316  ,0.285  ,0.197  ,0.167  ,0.133  ,0.096  ,0.045  ,0.010};
  double phistat []={0.011  ,0.021  ,0.020  ,0.019  ,0.017  ,0.007  ,0.004  ,0.002  ,0.001};
  double phisyst []={0.005,0.007,0.009,0.006,0.002,0.002,0.001,0.001,0.000};
  double phierror[9];
  for(unsigned int ix=0;ix<9;++ix) 
    phierror[ix] = sqrt(sqr(phistat[ix])+sqr(phisyst[ix]));
  bins =vector<double>(phibins ,phibins +10);
  data =vector<double>(phidata ,phidata + 9);
  error=vector<double>(phierror,phierror+ 9);
  _xpphi = new_ptr(Histogram(bins,data,error));
  // K*0
  double Kstar0bins []={0.0    ,0.01   ,0.03   ,0.1    ,0.125  ,
			0.14   ,0.16   ,0.2    ,0.3    ,0.4    ,0.5    ,0.7    ,1.0 };
  double Kstar0data []={1.22   ,4.96   ,4.14   ,2.35   ,1.99   ,
			1.60   ,1.30   ,0.81   ,0.44   ,0.22   ,0.090  ,0.013  	    };
  double Kstar0stat []={0.15   ,0.17   ,0.20   ,0.16   ,0.15   ,
			0.11   ,0.09   ,0.04   ,0.03   ,0.02   ,0.009  ,0.004  	    };
  double Kstar0syst []={0.04 ,0.15 ,0.19 ,0.13 ,0.09 ,
			0.10 ,0.06 ,0.05 ,0.03 ,0.01 ,0.003,0.003};
  double Kstar0error[12];
  for(unsigned int ix=0;ix<12;++ix) 
    Kstar0error[ix] = sqrt(sqr(Kstar0stat[ix])+sqr(Kstar0syst[ix]));
  bins =vector<double>(Kstar0bins ,Kstar0bins +13);
  data =vector<double>(Kstar0data ,Kstar0data +12);
  error=vector<double>(Kstar0error,Kstar0error+12);
  _xpKstar0 = new_ptr(Histogram(bins,data,error));
  // D0
  double D0bins []={0.15  ,0.25  ,0.35  ,0.45  ,0.55  ,0.65  ,0.75  ,0.85  ,1.00};
  double D0data []={41.1  ,101.4  , 58.4  , 52.9  , 53.2  , 36.0  , 13.5  ,  2.0};
  double D0error[]={11.3  ,13.8  ,10.5  , 8.1  , 8.4  , 8.0  , 3.7  , 2.0};
  bins =vector<double>(D0bins ,D0bins + 9);
  data =vector<double>(D0data ,D0data + 8);
  error=vector<double>(D0error,D0error+ 8);
  _xeD0 = new_ptr(Histogram(bins,data,error));
  // Dstar
  double Dstarbins []={0.10  ,0.15  ,0.20  ,0.25  ,0.30  ,
		       0.35  ,0.40  ,0.45  ,0.50  ,0.55  ,
		       0.60  ,0.65  ,0.70  ,0.75  ,0.80  ,
		       0.85  ,0.90  ,0.95  ,1.00};
  double Dstardata []={7.47  , 9.03  ,10.42  ,10.76  , 9.89  ,
		       8.97  , 8.17  , 6.94  , 6.73  , 5.56  ,
		       4.94  , 3.49  , 3.13  , 2.00  , 1.27  ,
		       0.50  , 0.27  , 0.06};
  double Dstarerror[]={0.63  ,0.49  ,0.44  ,0.43  ,0.38  ,
		       0.35  ,0.32  ,0.28  ,0.27  ,0.24  ,
		       0.22  ,0.18  ,0.17  ,0.14  ,0.11  ,
		       0.07  ,0.05  ,0.03};
  bins =vector<double>(Dstarbins ,Dstarbins +19);
  data =vector<double>(Dstardata ,Dstardata +18);
  error=vector<double>(Dstarerror,Dstarerror+18);
  _xeDstar = new_ptr(Histogram(bins,data,error));
  // rho0
  double rho0bins []={0.05 ,0.1  ,0.2  ,0.3  ,0.4  ,0.6  ,0.8  ,1.0};
  double rho0data []={6.15 ,2.16 ,0.92 ,0.45 ,0.13 ,0.027,0.003};
  double rho0error[]={0.72 ,0.23 ,0.10 ,0.05 ,0.02 ,0.005,0.002};
  bins =vector<double>(rho0bins ,rho0bins +8);
  data =vector<double>(rho0data ,rho0data +7);
  error=vector<double>(rho0error,rho0error+7);
  _xerho0 = new_ptr(Histogram(bins,data,error));
  // pi0
  double pi0binsA[]={0.007  ,0.009  ,0.011  ,0.013  ,0.016  ,
		     0.020  ,0.025  ,0.030  ,0.035  ,0.040  ,
		     0.050  ,0.060  ,0.070  ,0.085  ,0.100  ,
		     0.125  ,0.150  ,0.200  ,0.300  ,0.400  ,0.500};
  double pi0dataA[]={254.     ,266.     ,248.     ,211.     ,178.     ,
		     139.     ,113.     , 94.1    , 77.7    , 62.5    ,
		     45.7    , 34.7    , 26.2    , 19.4    , 13.2    ,
		     9.05   ,  5.36   ,  2.26   ,  0.764  ,  0.455};
  double pi0statA[]={18.     ,12.     , 6.     , 3.     , 2.     ,
		     2.     , 1.     , 0.9    , 0.8    , 0.4    ,
		     0.4    , 0.3    , 0.2    , 0.2    , 0.1    ,
		     0.13   , 0.10   , 0.13   , 0.085  , 0.095};
  double pi0systA[]={48   ,38   ,28   ,18   ,14   ,
		     6    ,5    ,4.0  ,4.3  ,3.9  ,
		     3.0  ,3.0  ,1.8  ,1.4  ,2.9  ,
		     0.76 ,0.69 ,0.38 ,0.309,0.244};
  double pi0errorA[20];
  for(unsigned int ix=0;ix<20;++ix) 
    pi0errorA[ix] = sqrt(sqr(pi0statA[ix])+sqr(pi0systA[ix]));
  bins =vector<double>(pi0binsA,pi0binsA+21);
  data =vector<double>(pi0dataA,pi0dataA+20);
  error=vector<double>(pi0errorA,pi0errorA+20);
  _xepi0 = new_ptr(Histogram(bins,data,error));
  double pi0binsB[]={0.69  ,0.92  ,1.20  ,1.61  ,1.90  ,
		     2.08  ,2.30  ,2.47  ,2.66  ,2.81  ,
		     3.00  ,3.22  ,3.36  ,3.51  ,3.70  ,
		     3.92  ,4.15  ,4.37  ,4.55  ,4.77  ,5.06};
  double pi0dataB[]={0.204  ,0.266  ,0.558  ,0.931  ,1.240  ,
		     1.48   ,1.79   ,2.02   ,2.25   ,2.50   ,
		     2.79   ,2.89   ,3.03   ,3.06   ,3.05   ,
		     3.11   ,2.92   ,2.78   ,2.42   ,1.74};
  double pi0statB[]={0.043  ,0.030  ,0.031  ,0.017  ,0.017  
		     ,0.02   ,0.02   ,0.02   ,0.02   ,0.02   ,
		     0.02   ,0.03   ,0.03   ,0.03   ,0.03   ,
		     0.03   ,0.05   ,0.06   ,0.11   ,0.12};
  double pi0systB[]={0.110 ,0.107 ,0.094 ,0.120 ,0.105 ,
		     0.32  ,0.13  ,0.14  ,0.20  ,0.16  ,
		     0.17  ,0.16  ,0.13  ,0.13  ,0.13  ,
		     0.25  ,0.25  ,0.31  ,0.34  ,0.33};
  double pi0errorB[20];
  for(unsigned int ix=0;ix<20;++ix) 
    pi0errorB[ix] = sqrt(sqr(pi0statB[ix])+sqr(pi0systB[ix]));
  bins =vector<double>(pi0binsB,pi0binsB+21);
  data =vector<double>(pi0dataB,pi0dataB+20);
  error=vector<double>(pi0errorB,pi0errorB+20);
  _xipi0 = new_ptr(Histogram(bins,data,error));
  // eta
  double etabinsA[]={0.025  ,0.035  ,0.050  ,0.075  ,0.100  ,
		     0.125  ,0.150  ,0.200  ,0.300  ,0.400  ,0.500  ,0.600  ,0.800  ,1.000};
  double etadataA[]={10.6   , 7.63  , 5.10  , 3.81  , 2.83  ,
		     2.21   , 1.46  , 0.733 , 0.364 , 0.220 , 0.086   , 0.033   , 0.0013};
  double etastatA[]={1.5    ,0.78   ,0.38   ,0.21    ,0.12  ,
		     0.10    ,0.05    ,0.026   ,0.022   ,0.019   ,0.010   ,0.004   ,0.0004};
  double etasystA[]={2.4   ,1.27  ,0.61  ,0.44  ,0.28  ,
		     0.22  ,0.13  ,0.062 ,0.047 ,0.031 ,0.019 ,0.008 ,0.0011};
  double etaerrorA[13];
  for(unsigned int ix=0;ix<13;++ix) 
    etaerrorA[ix] = sqrt(sqr(etastatA[ix])+sqr(etasystA[ix]));
  bins =vector<double>(etabinsA,etabinsA+14);
  data =vector<double>(etadataA,etadataA+13);
  error=vector<double>(etaerrorA,etaerrorA+13);
  _xeeta = new_ptr(Histogram(bins,data,error));
  double etabinsB[]={0.00  ,0.22  ,0.51  ,0.69  ,0.92  ,
		     1.20  ,1.61  ,1.90  ,2.08  ,2.31  ,2.60  ,3.03  ,3.42  ,3.82};
  double etadataB[]={0.0012  ,0.023   ,0.047   ,0.099   ,0.126   ,
		     0.180   ,0.252   ,0.301   ,0.314   ,0.324   ,0.302   ,0.294   ,0.261};
  double etastatB[]={0.0004  ,0.003   ,0.006   ,0.008   ,0.008   ,
		     0.006   ,0.009   ,0.014   ,0.014   ,0.018   ,0.023   ,0.030   ,0.038};
  double etasystB[]={0.0009 ,0.005  ,0.009  ,0.011  ,0.014  ,
		     0.014  ,0.021  ,0.027  ,0.028  ,0.032  ,0.028  ,0.038  ,0.046};
  double etaerrorB[13];
  for(unsigned int ix=0;ix<13;++ix) 
    etaerrorB[ix] = sqrt(sqr(etastatB[ix])+sqr(etasystB[ix]));
  bins =vector<double>(etabinsB,etabinsB+14);
  data =vector<double>(etadataB,etadataB+13);
  error=vector<double>(etaerrorB,etaerrorB+13);
  _xieta = new_ptr(Histogram(bins,data,error));
  // rho+  
  double rhopbinsA[]={0.016  ,0.025  ,0.035  ,0.050  ,0.075  ,
		      0.100  ,0.125  ,0.150  ,0.200  ,0.300  ,0.400  ,0.600  ,0.800, 1.000};
  double rhopdataA[]={17.3    ,32.3    ,21.3    ,16.7    , 9.89   ,
		      7.11   , 5.90   , 3.60   , 2.02   , 1.03   , 0.430  , 0.075  , 0.013};
  double rhopstatA[]={8.1    ,2.5    ,0.7    ,0.4    ,0.40   ,
		      0.25   ,0.25   ,0.12   ,0.07   ,0.04   ,0.023  ,0.013  ,0.003};
  double rhopsystA[]={12.2 ,9.7  ,4.5  ,1.8  ,1.46 ,1.04 ,0.78 ,0.48 ,0.21 ,0.27 ,0.081,0.032,0.009};
  double rhoperrorA[13];
  for(unsigned int ix=0;ix<13;++ix) 
    rhoperrorA[ix] = sqrt(sqr(rhopstatA[ix])+sqr(rhopsystA[ix]));
  bins =vector<double>(rhopbinsA,rhopbinsA+14);
  data =vector<double>(rhopdataA,rhopdataA+13);
  error=vector<double>(rhoperrorA,rhoperrorA+13);
  _xerhop = new_ptr(Histogram(bins,data,error));
  double rhopbinsB[]={0.0  ,0.5  ,1.0  ,1.5  ,2.0  ,2.5  ,3.0  ,3.5  ,4.0  ,4.5  ,5.0};
  double rhopdataB[]={0.034  ,0.217  ,0.419  ,0.603  ,0.805  ,0.868  ,0.692  ,0.500  ,0.419  ,0.171};
  double rhopstatB[]={0.004  ,0.010  ,0.014  ,0.017  ,0.022  ,0.021  ,0.028  ,0.092  ,0.035  ,0.008};
  double rhopsystB[]={0.019,0.055,0.073,0.078,0.104,0.126,0.165,0.138,0.111,0.081};
  double rhoperrorB[10];
  for(unsigned int ix=0;ix<10;++ix) 
    rhoperrorB[ix] = sqrt(sqr(rhopstatB[ix])+sqr(rhopsystB[ix]));
  bins =vector<double>(rhopbinsB,rhopbinsB+11);
  data =vector<double>(rhopdataB,rhopdataB+10);
  error=vector<double>(rhoperrorB,rhoperrorB+10);
  _xirhop = new_ptr(Histogram(bins,data,error));
  // omega
  double omegabinsA[]={0.025  ,0.035  ,0.050  ,0.075  ,0.100  ,0.125  ,0.150  ,0.200  ,0.300  ,0.400  ,0.600};
  double omegadataA[]={15.2    , 9.88   , 5.82   , 4.12   , 2.74   , 2.23   , 1.45   , 0.789  , 0.335  , 0.130};
  double omegastatA[]={2.4    ,0.84   ,0.35   ,0.25   ,0.16   ,0.14   ,0.09   ,0.049  ,0.037  ,0.027};
  double omegasystA[]={2.1  ,1.48 ,0.75 ,0.54 ,0.32 ,0.24 ,0.17 ,0.099,0.042,0.028};
  double omegaerrorA[10];
  for(unsigned int ix=0;ix<10;++ix) 
    omegaerrorA[ix] = sqrt(sqr(omegastatA[ix])+sqr(omegasystA[ix]));
  bins =vector<double>(omegabinsA,omegabinsA+11);
  data =vector<double>(omegadataA,omegadataA+10);
  error=vector<double>(omegaerrorA,omegaerrorA+10);
  _xeomega = new_ptr(Histogram(bins,data,error));
  double omegabinsB[]={0.51  ,0.92  ,1.21  ,1.61  ,1.90  ,2.09  ,2.32  ,2.62  ,3.06  ,3.49  ,4.01};
  double omegadataB[]={0.064  ,0.116  ,0.193  ,0.250  ,0.301  ,0.299  ,0.344  ,0.330  ,0.344  ,0.293};
  double omegastatB[]={0.013  ,0.013  ,0.012  ,0.016  ,0.018  ,0.018  ,0.021  ,0.020  ,0.029  ,0.046};
  double omegasystB[]={0.014,0.014,0.024,0.029,0.032,0.035,0.045,0.043,0.051,0.040};
  double omegaerrorB[10];
  for(unsigned int ix=0;ix<10;++ix) 
    omegaerrorB[ix] = sqrt(sqr(omegastatB[ix])+sqr(omegasystB[ix]));
  bins =vector<double>(omegabinsB,omegabinsB  +11);
  data =vector<double>(omegadataB,omegadataB  +10);
  error=vector<double>(omegaerrorB,omegaerrorB+10);
  _xiomega = new_ptr(Histogram(bins,data,error));
  // eta'
  double etapbinsA[]={0.050  ,0.070  ,0.100  ,0.125  ,0.150  ,0.200  ,0.300  ,0.400  ,0.600  ,0.800};
  double etapdataA[]={1.01   ,0.462  ,0.460  ,0.293  ,0.354  ,0.137  ,0.088  ,0.034  ,0.013};
  double etapstatA[]={0.38   ,0.180  ,0.144  ,0.099  ,0.068  ,0.028  ,0.020  ,0.010  ,0.006};
  double etapsystA[]={0.14 ,0.073,0.082,0.049,0.054,0.017,0.011,0.006,0.003};
  double etaperrorA[9];
  for(unsigned int ix=0;ix<9;++ix) 
    etaperrorA[ix] = sqrt(sqr(etapstatA[ix])+sqr(etapsystA[ix]));
  bins =vector<double>(etapbinsA,etapbinsA+10);
  data =vector<double>(etapdataA,etapdataA+9);
  error=vector<double>(etaperrorA,etaperrorA+9);
  _xeetap = new_ptr(Histogram(bins,data,error));
  double etapbinsB[]={0.22  ,0.51  ,0.92  ,1.21  ,1.61  ,1.91  ,2.09  ,2.33  ,2.71  ,3.09};
  double etapdataB[]={0.009  ,0.017  ,0.030  ,0.034  ,0.061  ,0.039  ,0.050  ,0.036  ,0.052};
  double etapstatB[]={0.004  ,0.005  ,0.007  ,0.007  ,0.012  ,0.013  ,0.016  ,0.014  ,0.020};
  double etapsystB[]={0.002,0.003,0.004,0.004,0.009,0.007,0.009,0.006,0.007};
  double etaperrorB[9];
  for(unsigned int ix=0;ix<9;++ix) 
    etaperrorB[ix] = sqrt(sqr(etapstatB[ix])+sqr(etapsystB[ix]));
  bins =vector<double>(etapbinsB,etapbinsB+10);
  data =vector<double>(etapdataB,etapdataB+9);
  error=vector<double>(etaperrorB,etaperrorB+9);
  _xietap = new_ptr(Histogram(bins,data,error));
  // a_0+
  double a0binsA[]={0.050  ,0.070  ,0.100  ,0.125  ,0.150  ,0.200  ,0.300  ,0.400  ,0.600  ,0.800  ,1.000};
  double a0dataA[]={1.65    ,1.05    ,0.747   ,0.985   ,0.623   ,0.207   ,0.093   ,0.038   ,0.014   ,0.0040};
  double a0statA[]={1.03    ,0.49    ,0.215   ,0.238   ,0.107   ,0.046   ,0.027   ,0.015   ,0.005   ,0.0018};
  double a0systA[]={0.75  ,0.73  ,0.214 ,0.560 ,0.171 ,0.069 ,0.040 ,0.015 ,0.006 ,0.0024};
  double a0errorA[10];
  for(unsigned int ix=0;ix<10;++ix) 
    a0errorA[ix] = sqrt(sqr(a0statA[ix])+sqr(a0systA[ix]));
  bins =vector<double>(a0binsA,a0binsA+11);
  data =vector<double>(a0dataA,a0dataA+10);
  error=vector<double>(a0errorA,a0errorA+10);
  _xea_0p = new_ptr(Histogram(bins,data,error));
  double a0binsB[]={0.00  ,0.50  ,1.00  ,1.50  ,2.00  ,2.50  ,3.00  ,3.50};
  double a0dataB[]={0.0071  ,0.019   ,0.040   ,0.088   ,0.076   ,0.104   ,0.093};
  double a0statB[]={0.0025  ,0.006   ,0.009   ,0.013   ,0.019   ,0.041   ,0.063};
  double a0systB[]={0.0022,0.007 ,0.012 ,0.023 ,0.030 ,0.041 ,0.050};
  double a0errorB[7];
  for(unsigned int ix=0;ix<7;++ix) 
    a0errorB[ix] = sqrt(sqr(a0statB[ix])+sqr(a0systB[ix]));
  bins =vector<double>(a0binsB,a0binsB+8);
  data =vector<double>(a0dataB,a0dataB+7);
  error=vector<double>(a0errorB,a0errorB+7);
  _xia_0p = new_ptr(Histogram(bins,data,error));
  // K_0
  double K0bins []={0.0114  ,0.020   ,0.030   ,0.040   ,0.050   ,
		    0.060   ,0.070   ,0.080   ,0.090   ,0.100   ,
		    0.125   ,0.150   ,0.200   ,0.250   ,0.300   ,
		    0.350   ,0.400   ,0.450   ,0.500   ,0.600  ,0.800};
  double K0data []={25.731  ,24.617  ,19.349  ,15.500  ,13.170  ,
		    11.144  , 9.360  , 8.470  , 7.010  , 5.734  ,
		    4.488  , 3.100  , 1.945  , 1.266  , 0.860  ,
		    0.579  , 0.394  , 0.253  , 0.163  , 0.051};
  double K0stat []={0.232  ,0.120  ,0.116  ,0.061  ,0.072  ,
		    0.073  ,0.066  ,0.061  ,0.059  ,0.029  ,
		    0.028  ,0.019  ,0.015  ,0.010  ,0.010  ,
		    0.009  ,0.008  ,0.005  ,0.003  ,0.001};
  double K0syst []={1.430,1.300,1.040,0.767,0.690,
		    0.600,0.500,0.468,0.401,0.312,
		    0.247,0.169,0.104,0.071,0.050,
		    0.035,0.026,0.018,0.018,0.010};
  double K0error[20];
  for(unsigned int ix=0;ix<20;++ix) 
    K0error[ix] = sqrt(sqr(K0stat[ix])+sqr(K0syst[ix]));
  bins =vector<double>(K0bins ,K0bins +21);
  data =vector<double>(K0data ,K0data +20);
  error=vector<double>(K0error,K0error+20);
  _xpK0 = new_ptr(Histogram(bins,data,error));
}

      
  

  

   


//   // a_0+
//   double a0binsA[]={};
//   double a0dataA[]={};
//   double a0statA[]={};
//   double a0systA[]={};
//   double a0errorA[7];
//   for(unsigned int ix=0;ix<10;++ix) 
//     a0errorA[ix] = sqrt(sqr(a0statA[ix])+sqr(a0systA[ix]));
//   bins =vector<double>(a0binsA,a0binsA+11);
//   data =vector<double>(a0dataA,a0dataA+10);
//   error=vector<double>(a0errorA,a0errorA+10);
//   _xea_0p = new_ptr(Histogram(bins,data,error));
//   double a0binsB[]={};
//   double a0dataB[]={};
//   double a0statB[]={};
//   double a0systB[]={};
//   double a0errorB[7];
//   for(unsigned int ix=0;ix<7;++ix) 
//     a0errorB[ix] = sqrt(sqr(a0statB[ix])+sqr(a0systB[ix]));
//   bins =vector<double>(a0binsB,a0binsB+8);
//   data =vector<double>(a0dataB,a0dataB+7);
//   error=vector<double>(a0errorB,a0errorB+7);
//   _xia_0p = new_ptr(Histogram(bins,data,error));

