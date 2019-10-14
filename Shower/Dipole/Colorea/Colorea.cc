  // -*- C++ -*-
  //
  // Colorea.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2017 The Herwig Collaboration
  //
  // Herwig is licenced under version 3 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
  //
  // This is the implementation of the non-inlined, non-templated member
  // functions of the Colorea class.
  //

#include "Colorea.h"
#include "eeuugg.h"
#include "eeuuggg.h"
#include "eeuugggg.h"
#include "Herwig/Shower/Dipole/Utility/DipolePartonSplitter.h"
#include "Herwig/Shower/Dipole/Base/DipoleChain.h"
#include "Herwig/Shower/Dipole/Base/Dipole.h"


#include <iterator>

using namespace Herwig;

Colorea::Colorea(){}

list<Dipole>& Colorea::dipoles() { return theChain->dipoles(); }

void Colorea::rearrange( int dipmax , int diplong ){
  
  assert( dipmax >= diplong );
    // if there are only 2 dipoles in the chain
    //there is nothing to do for now.
  if( dipoles().size() < 3 )return;
  
  if( dipoles().size() == 3 ){
      // get the 3 dipoles:
    auto dipi=dipoles().begin();
    auto dipj=dipoles().begin();dipj++;
    auto dipk=dipoles().begin();dipk++;dipk++;
    rearrange3(dipi,dipj,dipk);
    
    if( theChain->circular () ){
        // if the chain is circular,
        // we need to check also the
        // connected ends.
      rearrange3(dipj,dipk,dipi);
      rearrange3(dipk,dipi,dipj);
    }
    return;
  }
  if( dipoles().size() == 4 && dipmax >= 4 ){
      // get the 4 dipoles
    auto dipi=dipoles().begin();
    auto dipj=dipoles().begin();dipj++;
    auto dipk=dipoles().begin();dipk++;dipk++;
    auto dipl=dipoles().begin();dipl++;dipl++;dipl++;
    rearrange4(dipi,dipj,dipk,dipl);
    return;
  }
  if( dipoles().size() == 5 && dipmax >= 5 ){
      // get the 5 dipoles
    auto dipi=dipoles().begin();
    auto dipj=dipoles().begin();dipj++;
    auto dipk=dipoles().begin();dipk++;dipk++;
    auto dipl=dipoles().begin();dipl++;dipl++;dipl++;
    auto dipm=dipoles().begin();dipm++;dipm++;dipm++;dipm++;
    rearrange5(dipi,dipj,dipk,dipl,dipm);
    return;
  }
    // if the chain is longer than dipmax
    // go though the chain with diplong dipoles.
  rearrangeLong( diplong );
  return;
  
}

void Colorea::rearrangeLong( int diplong ){
  if ( diplong == 3 ){
      // get the 3 dipoles:
    auto dipi=dipoles().begin();
    auto dipj=dipoles().begin();dipj++;
    auto dipk=dipoles().begin();dipk++;dipk++;
    while (dipk!=dipoles().end()) {
      rearrange3(dipi,dipj,dipk);
      dipi++;dipj++;dipk++;
    }
    if( theChain->circular() ){
        // if the chain is circular,
        // we need to check also the
        // connected ends.
      dipk=dipoles().begin(); // as dipk was end
      rearrange3( dipi , dipj , dipk );
        // rotate further.
      dipi++;dipk++;
      dipj=dipoles().begin();
      rearrange3(dipi,dipj,dipk);
    }
  }else if (diplong ==4){
      // get the 4 dipoles:
    assert(!theChain->circular ());
    auto dipi=dipoles().begin();
    auto dipj=dipoles().begin();dipj++;
    auto dipk=dipoles().begin();dipk++;dipk++;
    auto dipl=dipoles().begin();dipl++;dipl++;dipl++;
    while (dipl!=dipoles().end()) {
      rearrange4(dipi,dipj,dipk,dipl);
      dipi++;dipj++;dipk++;dipl++;
    }
  }else if (diplong ==5){
      // get the 5 dipoles:
    assert(!theChain->circular ());
    auto dipi=dipoles().begin();
    auto dipj=dipoles().begin();dipj++;
    auto dipk=dipoles().begin();dipk++;dipk++;
    auto dipl=dipoles().begin();dipl++;dipl++;dipl++;
    auto dipm=dipoles().begin();dipm++;dipm++;dipm++;dipm++;
    
    while (dipm!=dipoles().end()) {
      rearrange5(dipi,dipj,dipk,dipl,dipm);
      dipi++;dipj++;dipk++;dipl++;dipm++;
    }
  }else{
    assert(false);
  }
}



void Colorea::rmcol(tColinePtr A,
                    tColinePtr B,
                    list<Dipole>::iterator & dip){
  if(A==B){
    auto AA=A->coloured();
    for (auto col:AA){
      if(col==dip->leftParticle()) A->removeColoured (col);
      if(col==dip->rightParticle())A->removeColoured (col);
    }
    auto BB=A->antiColoured();
    for (auto col:BB){
      if(col==dip->leftParticle()) A->removeAntiColoured (col);
      if(col==dip->rightParticle())A->removeAntiColoured (col);
    }
  }
}



void Colorea::removeColors(list<Dipole>::iterator dipi,
                           list<Dipole>::iterator dipj,
                           list<Dipole>::iterator dipk){
  auto ilc = dipi->leftParticle()->colourInfo()->    colourLine();
  auto irc = dipi->rightParticle()->colourInfo()->   colourLine();
  auto ila = dipi->leftParticle()->colourInfo()->    antiColourLine();
  auto ira = dipi->rightParticle()->colourInfo()->   antiColourLine();
  auto jlc = dipj->leftParticle()->colourInfo()->    colourLine();
  auto jrc = dipj->rightParticle()->colourInfo()->   colourLine();
  auto jla = dipj->leftParticle()->colourInfo()->    antiColourLine();
  auto jra = dipj->rightParticle()->colourInfo()->   antiColourLine();
  auto klc = dipk->leftParticle()->colourInfo()->    colourLine();
  auto krc = dipk->rightParticle()->colourInfo()->   colourLine();
  auto kla = dipk->leftParticle()->colourInfo()->    antiColourLine();
  auto kra = dipk->rightParticle()->colourInfo()->   antiColourLine();
  
  if( ilc && irc ) rmcol(ilc,irc,dipi);
  if( ila && ira ) rmcol(ila,ira,dipi);
  if( jlc && jrc ) rmcol(jlc,jrc,dipj);
  if( jla && jra ) rmcol(jla,jra,dipj);
  if( klc && krc ) rmcol(klc,krc,dipk);
  if( kla && kra ) rmcol(kla,kra,dipk);
  
  return;
}



void Colorea::rearrange3_FF_FF_FF(list<Dipole>::iterator dipi,
                                  list<Dipole>::iterator dipj,
                                  list<Dipole>::iterator dipk){
  produceSwapping( dipi, dipj, dipk,true);
    // compared to IF of FI dipoles we do not need to repair colors as
    // the swapping just modified the gluon momenta.
}





void Colorea::rearrange3_FF_FI_IF(list<Dipole>::iterator dipi,
                                  list<Dipole>::iterator dipj,
                                  list<Dipole>::iterator dipk){
  auto ilc = dipi->leftParticle()->colourInfo()->    colourLine();
  auto irc = dipi->rightParticle()->colourInfo()->   colourLine();
  auto ila = dipi->leftParticle()->colourInfo()->    antiColourLine();
  auto ira = dipi->rightParticle()->colourInfo()->   antiColourLine();
  auto jlc = dipj->leftParticle()->colourInfo()->    colourLine();
  auto jrc = dipj->rightParticle()->colourInfo()->   colourLine();
  auto jla = dipj->leftParticle()->colourInfo()->    antiColourLine();
  auto jra = dipj->rightParticle()->colourInfo()->   antiColourLine();
  auto klc = dipk->leftParticle()->colourInfo()->    colourLine();
  auto krc = dipk->rightParticle()->colourInfo()->   colourLine();
  auto kla = dipk->leftParticle()->colourInfo()->    antiColourLine();
  auto kra = dipk->rightParticle()->colourInfo()->   antiColourLine();
  
  bool didperm = produceSwapping( dipi, dipj, dipk,false);
  
  if (!didperm)return;
  removeColors(dipi,dipj,dipk);
  auto new1Line = new_ptr(ColourLine());
  auto new2Line = new_ptr(ColourLine());
  auto new3Line = new_ptr(ColourLine());
    //Example:
    // before
    //ila (F)--irc (F) =g1(in,jl)= jla (F) -- jra(I)=g2(in,jr)=klc(I) -- krc(F)
    // Should be after:
    //ila (F) -- ira (I) =g2(in,jl)= jlc (I) -- jrc(F) =g1(out,jr)= kla(F) -- krc(F)
  
    //Then:
    // FaFc_FaIa_IcFc    FFFIIF  aa -> FaIa_IcFc_FaFc
    // FcFa_FcIc_IaFa    FFFIIF  cc -> FcIc_IaFa_FcFa
    // FaFc_FaIa_IcIa    FFFIII  aa -> FaIa_IcFc_FaIa
    // FcFa_FcIc_IaIc    FFFIII  cc -> FcIc_IaFa_FcIc
    // IcFc_FaIa_IcFc    IFFIIF  aa -> IcIa_IcFc_FaFc
    // IaFa_FcIc_IaFa    IFFIIF  cc -> IaIc_IaFa_FcFa
    // IcFc_FaIa_IcIa    IFFIII  aa -> IcIa_IcFc_FaIa
    // IaFa_FcIc_IaIc    IFFIII  cc -> IaIc_IaFa_FcIc
    // end particles keep colors
    // central all change color.
    //
  
  if(ilc==irc){
    new1Line->addColoured        (dipi->leftParticle());
    new1Line->addAntiColoured    (dipi->rightParticle());}
  if(ilc==ira){
    new1Line->addColoured        (dipi->leftParticle());
    new1Line->addColoured        (dipi->rightParticle());}
  if(ila==irc){
    new1Line->addAntiColoured    (dipi->leftParticle());
    new1Line->addAntiColoured    (dipi->rightParticle());}
  if(ila==ira){
    new1Line->addAntiColoured    (dipi->leftParticle());
    new1Line->addColoured        (dipi->rightParticle());}
  
  if(jlc==jrc){
    new2Line->addAntiColoured    (dipj->leftParticle());
    new2Line->addAntiColoured    (dipj->rightParticle());}
  if(jlc==jra){
    new2Line->addAntiColoured    (dipj->leftParticle());
    new2Line->addColoured        (dipj->rightParticle());}
  if(jla==jrc){
    new2Line->addColoured        (dipj->leftParticle());
    new2Line->addAntiColoured    (dipj->rightParticle());}
  if(jla==jra){
    new2Line->addColoured        (dipj->leftParticle());
    new2Line->addColoured        (dipj->rightParticle());}
  
  if(klc==krc){
    new3Line->addAntiColoured    (dipk->leftParticle());
    new3Line->addColoured        (dipk->rightParticle());}
  if(klc==kra){
    new3Line->addAntiColoured    (dipk->leftParticle());
    new3Line->addAntiColoured    (dipk->rightParticle());}
  if(kla==krc){
    new3Line->addColoured        (dipk->leftParticle());
    new3Line->addColoured        (dipk->rightParticle());}
  if(kla==kra){
    new3Line->addColoured        (dipk->leftParticle());
    new3Line->addAntiColoured    (dipk->rightParticle());}
  
}



void Colorea::rearrange3(list<Dipole>::iterator dipi,
                         list<Dipole>::iterator dipj,
                         list<Dipole>::iterator dipk){
  
    // We only care about the cenral dipole.
    // If central dipole is FI or IF we need to take care that
    // momentum fractions and pdfs are treated correctly.
    // If central dipole is FF we only swap momenta.
    // Also the color lines of IF and FI swappings
    // needs to be repaired.
  
  bool dipjisFF = dipj->leftFraction() == 1 && dipj->rightFraction() == 1;
  bool dipjisIF = dipj->leftFraction() != 1 && dipj->rightFraction() == 1;
  bool dipjisFI = dipj->leftFraction() == 1 && dipj->rightFraction() != 1;
  bool dipjisII = dipj->leftFraction() != 1 && dipj->rightFraction() != 1;
  
    // if the possibly swaping dipole is II we
    // do nothing as this should not be swapped!(?)
  if(dipjisII) return;
  
  if(dipjisFF){
      // possibly swap gluon momenta
    rearrange3_FF_FF_FF(dipi,dipj,dipk);
  }
    // If central dipole is IF or FI -> keep end colors and swap central colors
  else if( dipjisIF || dipjisFI ){
    rearrange3_FF_FI_IF(dipi,dipj,dipk);
  }
  else {
    assert(false);
  }
  
  
  return;
  
}

  // local helper function to fill madgraph momenta.
void fillMGmom(double * mom,Lorentz5Momentum lmom){
  mom[0]=lmom.t()/GeV;
  mom[1]=lmom.x()/GeV;
  mom[2]=lmom.y()/GeV;
  mom[3]=lmom.z()/GeV;
}


bool Colorea::produceSwapping(list<Dipole>::iterator dipi,
                              list<Dipole>::iterator dipj,
                              list<Dipole>::iterator dipk,bool FF){
  
    // get the process for all tripple dipole swappings.
  static auto pro= eeuugg();
  
    // set mass of q qbar
  pro.setMass( 3 , dipi->leftParticle()->momentum().m()/GeV);
  pro.setMass( 4 , dipk->rightParticle()->momentum().m()/GeV);
  
    // If the central dipole dipj is not a FF-dipole the rearragment is
    // not a simple swapping of gluon momenta. We need to take care that
    // the 'new' dipoles get the correct fraction, pdf and index assignments.
  if(!FF){
      // for the interfacing with madgraph we need new momenta
    double  mom0[4], mom1[4], mom2[4];
    double  mom3[4], mom4[4], mom5[4];
    
      //q
    double prefact=1.;//dipi->leftFraction()!=0.?-1.:1.; ??
    auto PP=dipi->leftParticle()->momentum();
    PP.setX(prefact*PP.x());
    PP.setY(prefact*PP.y());
    PP.setZ(prefact*PP.z());
      // get the momentum sum for the triple dipole system.
      // we need this for the incoming e+e-.
    auto fmomsum=PP;
      // lets rotate the system to be sure that momenta
      // are not misinterpreted as incomming.
    PP=PP.rotateX(0.3);
    fillMGmom(mom2,PP);
    
    
      //qbar
    prefact=1.;//dipk->rightFraction()!=0.?-1.:1.;
    PP=dipk->rightParticle()->momentum();
    PP.setX(prefact*PP.x());
    PP.setY(prefact*PP.y());
    PP.setZ(prefact*PP.z());
    fmomsum+=PP;
    PP=PP.rotateX(0.3);
    fillMGmom(mom3,PP);
    
      //g1
    prefact=1.;//dipj->leftFraction()!=0.?-1.:1.;
    PP=dipj->leftParticle()->momentum();
    PP.setX(prefact*PP.x());
    PP.setY(prefact*PP.y());
    PP.setZ(prefact*PP.z());
    fmomsum+=PP;
    PP=PP.rotateX(0.3);
    fillMGmom(mom4,PP);
    
      //g2
    prefact=1.;//dipj->rightFraction()!=0.?-1.:1.;
    PP=dipj->rightParticle()->momentum();
    PP.setX(prefact*PP.x());
    PP.setY(prefact*PP.y());
    PP.setZ(prefact*PP.z());
    fmomsum+=PP;
    PP=PP.rotateX(0.3);
    fillMGmom(mom5,PP);
    
      //  get the boost to to triple dipole restframe
    auto boostto=fmomsum.findBoostToCM();
      // center of mass energy.
    auto CME=fmomsum.m2();
      // construct incoming momenta.
      // point incoming axis along arbitrary y axis;
    Lorentz5Momentum in1=Lorentz5Momentum(ZERO,sqrt(CME)/2.,ZERO,sqrt(CME)/2.);
    Lorentz5Momentum in2=Lorentz5Momentum(ZERO,-sqrt(CME)/2.,ZERO,sqrt(CME)/2.);
    
    in1=in1.boost(-boostto);
    in1=in1.rotateX(0.3);
    in2=in2.boost(-boostto);
    in2=in2.rotateX(0.3);
    fillMGmom(mom0,in1);
    fillMGmom(mom1,in2);
    
    vector < double * >  momenta{{ mom0, mom1, mom2,
                                   mom3, mom4, mom5}};

      // calculate the permutation from madgraph process.
    auto perm=pro.producePermutation(UseRandom::rnd(),momenta);
      // for tripple dipoles there is only one possible permutation.
    bool didperm=perm[0]!=5;
    
    if(!didperm) return didperm;
    
      // insert the original dipj fractions
    vector<double > tmpfrac{{dipj->leftFraction(),dipj->rightFraction()}};
      // insert the original dipj pdfs
    vector< PDF > tmppdf;
    tmppdf.push_back(dipj->leftPDF());
    tmppdf.push_back(dipj->rightPDF());
      // insert the original dipj particles
    vector<PPtr> tmppart;
    tmppart.push_back(dipj->leftParticle());
    tmppart.push_back(dipj->rightParticle());
      // make sure the dipj are gluons.
    assert(dipj->leftParticle()->id()==21);
    assert(dipj->rightParticle()->id()==21);
    
      // Now set the new or old particles in the dipole
    dipi->rightParticle(tmppart[perm[0]-1-4]);
    dipj->leftParticle (tmppart[perm[0]-1-4]);
    dipj->rightParticle(tmppart[perm[1]-1-4]);
    dipk->leftParticle(tmppart[perm[1]-1-4]);
    
    dipi->rightPDF(tmppdf[perm[0]-1-4]);
    dipj->leftPDF (tmppdf[perm[0]-1-4]);
    dipj->rightPDF(tmppdf[perm[1]-1-4]);
    dipk->leftPDF (tmppdf[perm[1]-1-4]);
    
    dipi->rightFraction(tmpfrac[perm[0]-1-4]);
    dipj->leftFraction (tmpfrac[perm[0]-1-4]);
    dipj->rightFraction(tmpfrac[perm[1]-1-4]);
    dipk->leftFraction (tmpfrac[perm[1]-1-4]);
    
    
    auto firstindex =DipoleIndex(&(dipj->leftParticle()->data()),
                                 &(dipj->rightParticle()->data()),
                                 dipj->leftPDF(),
                                 dipj->rightPDF());
    auto secondindex =DipoleIndex(&(dipj->rightParticle()->data()),
                                  &(dipj->leftParticle()->data()),
                                  dipj->rightPDF(),
                                  dipj->leftPDF());
    dipj->setFirstIndex(firstindex);
    dipj->setSecondIndex(secondindex);
    
    auto firstindex2 =DipoleIndex(&(dipi->leftParticle()->data()),
                                  &(dipi->rightParticle()->data()),
                                  dipi->leftPDF(),
                                  dipi->rightPDF());
    auto secondindex2 =DipoleIndex(&(dipi->rightParticle()->data()),
                                   &(dipi->leftParticle()->data()),
                                   dipi->rightPDF(),
                                   dipi->leftPDF());
    dipi->setFirstIndex(firstindex2);
    dipi->setSecondIndex(secondindex2);
    
    auto firstindex3 =DipoleIndex(&(dipk->leftParticle()->data()),
                                  &(dipk->rightParticle()->data()),
                                  dipk->leftPDF(),
                                  dipk->rightPDF());
    auto secondindex3 =DipoleIndex(&(dipk->rightParticle()->data()),
                                   &(dipk->leftParticle()->data()),
                                   dipk->rightPDF(),
                                   dipk->leftPDF());
    dipk->setFirstIndex(firstindex3);
    dipk->setSecondIndex(secondindex3);
      // check a few things
      //    assert(dipj->leftFraction()!=1.&&(dipj->leftPDF().pdf() )  || dipj->leftFraction()==1. &&!(dipj->leftPDF().pdf()));
      //    assert(dipj->rightFraction()!=1.&&(dipj->rightPDF().pdf() )|| dipj->rightFraction()==1.&&!(dipj->rightPDF().pdf()));
      //    assert(dipi->leftFraction()!=1.&&(dipi->leftPDF().pdf() )  || dipi->leftFraction()==1. &&!(dipi->leftPDF().pdf()));
      //    assert(dipi->rightFraction()!=1.&&(dipi->rightPDF().pdf() )|| dipi->rightFraction()==1.&&!(dipi->rightPDF().pdf()));
      //    assert(dipk->leftFraction()!=1.&&(dipk->leftPDF().pdf() )  || dipk->leftFraction()==1. &&!(dipk->leftPDF().pdf()));
      //    assert(dipk->rightFraction()!=1.&&(dipk->rightPDF().pdf() )|| dipk->rightFraction()==1.&&!(dipk->rightPDF().pdf()));
    
    return didperm;
  }else{
      // if dipj is a FF dipole the color rearrangement is done by
      // swapping the gluon momenta. Let's do this!
    
    double  mom0[4];
    double  mom1[4];
    double  mom2[4];
    double  mom3[4];
    double  mom4[4];
    double  mom5[4];
    vector<Lorentz5Momentum> tmpvec;
    
      //q
    auto PP=dipi->leftParticle()->momentum();
    auto fmomsum=PP;
    PP=PP.rotateX(0.3);
    fillMGmom(mom2,PP);
    
      //qbar
    PP=dipk->rightParticle()->momentum();
    fmomsum+=PP;
    PP=PP.rotateX(0.3);
    fillMGmom(mom3,PP);
    
      //g1
    PP=dipj->leftParticle()->momentum();
    fmomsum+=PP;
    PP=PP.rotateX(0.3);
    fillMGmom(mom4,PP);
    
      //g2
    PP=dipj->rightParticle()->momentum();
    fmomsum+=PP;
    PP=PP.rotateX(0.3);
    fillMGmom(mom5,PP);
    
      //  get the boost to to triple dipole restframe
    auto boostto=fmomsum.findBoostToCM();
    
    auto CME=fmomsum.m2();
      // point incoming axis along arbitrary y axis;
    Lorentz5Momentum in1=Lorentz5Momentum(ZERO,sqrt(CME)/2.,ZERO,sqrt(CME)/2.);
    Lorentz5Momentum in2=Lorentz5Momentum(ZERO,-sqrt(CME)/2.,ZERO,sqrt(CME)/2.);
    
    
    
    assert(dipj->leftParticle()->id()==21);
    assert(dipj->rightParticle()->id()==21);
    
    in1=in1.boost(-boostto);
    in1=in1.rotateX(0.3);
    in2=in2.boost(-boostto);
    in2=in2.rotateX(0.3);
    fillMGmom(mom0,in1);
    fillMGmom(mom1,in2);
    
    vector < double * >  momenta{{ mom0, mom1, mom2,
                                   mom3, mom4, mom5}};
    
      // calculate the permutation from madgraph process.
    auto perm=pro.producePermutation(UseRandom::rnd(),momenta);
    
    bool didperm=perm[0]!=5;
    
      // For final state we can just swap gluon momenta.
    tmpvec.push_back(dipj->leftParticle()->momentum());  // g1
    tmpvec.push_back(dipj->rightParticle()->momentum()); // g2
    
    dipj->leftParticle() ->setMomentum(tmpvec[perm[0]-1-4]);
    dipj->rightParticle()->setMomentum(tmpvec[perm[1]-1-4]);
    
    return didperm;
  }
  
}





void Colorea::rearrange4(list<Dipole>::iterator dipi,
                         list<Dipole>::iterator dipj,
#ifndef NDEBUG
                         list<Dipole>::iterator dipk,
#else
                         list<Dipole>::iterator ,
#endif
                         list<Dipole>::iterator dipl){
  
  assert(dipk->leftParticle()==dipj->rightParticle());
  assert(dipl->leftParticle()==dipk->rightParticle());
  
  double  mom0[4], mom1[4], mom2[4], mom3[4];
  double  mom6[4], mom4[4], mom5[4];
  
    // q       g1       g2      g3      qbar
    // il -- ir=jl -- jr=kl -- kr=ll -- lr
  
    //q
  fillMGmom(mom2,dipi->leftParticle()->momentum());
    //qbar
  fillMGmom(mom3,dipl->rightParticle()->momentum());
    //g1
  fillMGmom(mom4,dipj->leftParticle()->momentum());
    //g2
  fillMGmom(mom5,dipj->rightParticle()->momentum());
    //g3
  fillMGmom(mom6,dipl->leftParticle()->momentum());
  
    // q       g1       g2      g3      qbar
    // il -- ir=jl -- jr=kl -- kr=ll -- lr
  auto fmomsum=dipi->leftParticle()->momentum()+
  dipj->leftParticle()->momentum()+
  dipj->rightParticle()->momentum()+
  dipl->leftParticle()->momentum()+
  dipl->rightParticle()->momentum();
  
  auto boostto=fmomsum.findBoostToCM();
  
  auto CME=fmomsum.m2();
    // point incoming axis along arbitrary y axis;
  Lorentz5Momentum in1=Lorentz5Momentum(ZERO,sqrt(CME)/2.,ZERO,sqrt(CME)/2.);
  Lorentz5Momentum in2=Lorentz5Momentum(ZERO,-sqrt(CME)/2.,ZERO,sqrt(CME)/2.);
  
  in1=in1.boost(-boostto);
  in2=in2.boost(-boostto);
  fillMGmom(mom0,in1);
  fillMGmom(mom1,in2);
  
  vector < double * >  momenta{{ mom0, mom1,
                                 mom2, mom3,
                                 mom4, mom5, mom6}};
  
  static auto pro= eeuuggg();
  
    // set mass of q qbar
  pro.setMass( 3 , dipi->leftParticle()->momentum().m()/GeV);
  pro.setMass( 5 , dipl->rightParticle()->momentum().m()/GeV);
  
  
  auto perm=pro.producePermutation(UseRandom::rnd(),momenta);
  
  vector<Lorentz5Momentum> tmpvec;
  tmpvec.push_back(dipj->leftParticle()->momentum());  // g1
  tmpvec.push_back(dipj->rightParticle()->momentum()); // g2
  tmpvec.push_back(dipl->leftParticle()->momentum()); // g3
  
  dipj->leftParticle() ->setMomentum(tmpvec[perm[0]-1-4]);
  dipj->rightParticle()->setMomentum(tmpvec[perm[1]-1-4]);
  dipl->leftParticle() ->setMomentum(tmpvec[perm[2]-1-4]);
  
}





void Colorea::rearrange5(list<Dipole>::iterator dipi,
                         list<Dipole>::iterator dipj,
#ifndef NDEBUG
                         list<Dipole>::iterator dipk,
#else
                         list<Dipole>::iterator ,
#endif
                         list<Dipole>::iterator dipl,
                         list<Dipole>::iterator dipm){
  
  assert(dipk->leftParticle()==dipj->rightParticle());
  assert(dipl->leftParticle()==dipk->rightParticle());
  
  
  double  mom0[4], mom1[4], mom2[4], mom3[4];
  double  mom4[4], mom5[4], mom6[4], mom7[4];
    // q       g1       g2      g3      g4      qbar
    // il -- ir=jl -- jr=kl -- kr=ll -- lr=ml -- mr
  
    //q
  fillMGmom(mom2,dipi->leftParticle()->momentum());
    //qbar
  fillMGmom(mom3,dipm->rightParticle()->momentum());
    //g1
  fillMGmom(mom4,dipj->leftParticle()->momentum());
    //g2
  fillMGmom(mom5,dipj->rightParticle()->momentum());
    //g3
  fillMGmom(mom6,dipl->leftParticle()->momentum());
    //g4
  fillMGmom(mom7,dipm->leftParticle()->momentum());
  
    // q       g1       g2      g3      g4      qbar
    // il -- ir=jl -- jr=kl -- kr=ll -- lr=ml -- mr
  auto fmomsum=dipi->leftParticle()->momentum()+
  dipj->leftParticle()->momentum()+
  dipj->rightParticle()->momentum()+
  dipl->leftParticle()->momentum()+
  dipl->rightParticle()->momentum()+
  dipm->rightParticle()->momentum();
  
  auto boostto=fmomsum.findBoostToCM();
  
  auto CME=fmomsum.m2();
    // point incoming axis along arbitrary y axis;
  Lorentz5Momentum in1=Lorentz5Momentum(ZERO,sqrt(CME)/2.,ZERO,sqrt(CME)/2.);
  Lorentz5Momentum in2=Lorentz5Momentum(ZERO,-sqrt(CME)/2.,ZERO,sqrt(CME)/2.);
  
  
  in1=in1.boost(-boostto);
  in2=in2.boost(-boostto);
  fillMGmom(mom0,in1);
  fillMGmom(mom1,in2);
  
  vector < double * >  momenta{{ mom0, mom1, mom2, mom3,
                                 mom4, mom5, mom6, mom7}};
  
  static auto pro= eeuugggg();
  
    // set mass of q qbar
  pro.setMass( 3 , dipi->leftParticle()->momentum().m()/GeV);
  pro.setMass( 5 , dipl->rightParticle()->momentum().m()/GeV);
  
  auto perm=pro.producePermutation(UseRandom::rnd(),momenta);
  
  vector<Lorentz5Momentum> tmpvec;
  tmpvec.push_back(dipj->leftParticle()->momentum());  // g1
  tmpvec.push_back(dipj->rightParticle()->momentum()); // g2
  tmpvec.push_back(dipl->leftParticle()->momentum()); // g3
  tmpvec.push_back(dipm->leftParticle()->momentum()); // g4
  
    // set the gluon momenta with permuted momenta.
  dipj->leftParticle() ->setMomentum(tmpvec[perm[0]-1-4]);
  dipj->rightParticle()->setMomentum(tmpvec[perm[1]-1-4]);
  dipl->leftParticle() ->setMomentum(tmpvec[perm[2]-1-4]);
  dipm->leftParticle() ->setMomentum(tmpvec[perm[3]-1-4]);
  
}






