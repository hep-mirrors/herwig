// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEDiffraction class.
//

#include "MEDiffraction.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/StandardXComb.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#include "Herwig/Utilities/Kinematics.h"


MEDiffraction::MEDiffraction()
: HwMEBase(),
  deltaOnly(false),  
  isInRunPhase(false),
  theProtonMass(-MeV) // to be set in doinit
{}


void MEDiffraction::getDiagrams() const {
    //incoming particles
    cPDPair incomingHardons = generator()->eventHandler()->incoming();
      
    tcPDPtr pom = getParticleData(990);
    
    //get incoming particles
    tcPDPtr prt11 = getParticleData(incomingHardons.first->id());
    tcPDPtr prt12 = getParticleData(incomingHardons.second->id());
    
    //get sign of id
    int sign1=0, sign2=0;
    sign1 = (incomingHardons.first->id() > 0) ? 1 : -1;
    sign2 = (incomingHardons.second->id() > 0) ? 1 : -1;
    
    tcPDPtr prt21 = getParticleData(sign1*2214);//Delta+
    tcPDPtr prt22 = getParticleData(sign2*2214);//Delta+    
    
    //for the left side
    tcPDPtr q11 = getParticleData(sign1*2); //u
    tcPDPtr q21 = getParticleData(sign1*1); //d
    //for the right side
    tcPDPtr q12 = getParticleData(sign2*2); //u
    tcPDPtr q22 = getParticleData(sign2*1); //d
    //for the left side
    tcPDPtr dq11 = getParticleData(sign1*2101); //ud_0
    tcPDPtr dq111 = getParticleData(sign1*2103); //ud_1
    tcPDPtr dq21 = getParticleData(sign1*2203); //uu_1
    //for the right side
    tcPDPtr dq12 = getParticleData(sign2*2101); //ud_0
    tcPDPtr dq112 = getParticleData(sign2*2103); //ud_1
    tcPDPtr dq22 = getParticleData(sign2*2203); //uu_1
    
    tcPDPtr gl = getParticleData(21);//gluon
    
    //switch between dissociation decays to different 
    //number of clusters or dissociation into delta only
    //(Maybe can be automated???)
    //(Should be generalized to ppbar, for example!!!)
    switch(dissociationDecay){
      case 0: //one cluster or only delta in the final state 
        if(deltaOnly) //only delta in the final state
        {
            
          switch (diffDirection){
          case 0:
            add(new_ptr((Tree2toNDiagram(3), prt11, pom, prt12, 1, prt21, 2, prt12, -1)));
            break;
          case 1:
            add(new_ptr((Tree2toNDiagram(3), prt11, pom, prt12, 1, prt11, 2, prt22, -1)));
            break;
          case 2:
            add(new_ptr((Tree2toNDiagram(3), prt11, pom, prt12, 1, prt21, 2, prt22, -1)));
            break;
            } 
      
      
        }else
        {
          //switch between direction of dissociated proton for single diffraction or 
          //double diffraction 
          switch (diffDirection){
          case 0: //left
            //u -- ud_0
            add(new_ptr((Tree2toNDiagram(4), prt11, q11, pom, prt12, 3, prt12, 1, dq11, 2, q11, -1)));
            //d -- uu_1
            add(new_ptr((Tree2toNDiagram(4), prt11, q21, pom, prt12, 3, prt12, 1, dq21, 2, q21, -2)));
            break;
          case 1: //right
            //u -- ud_0
            add(new_ptr((Tree2toNDiagram(4), prt11, pom, q12, prt12, 1, prt11, 3, dq12, 2, q12, -1)));
            
            //d -- uu_1
            add(new_ptr((Tree2toNDiagram(4), prt11, pom, q22, prt12, 1, prt11, 3, dq22, 2, q22, -2)));
            break;
          case 2: //double
            //u -- ud_0 left u -- ud_0 right  
            add(new_ptr((Tree2toNDiagram(5), prt11, q11, pom, q12, prt12, 1, dq11, 2, q11, 3, q12, 4, dq12, -1)));
            
            //u -- ud_0 left d -- uu_1 right
            add(new_ptr((Tree2toNDiagram(5), prt11, q11, pom, q22, prt12, 1, dq11, 2, q11, 3, q22, 4, dq22, -2)));
            
            //d -- uu_1 left u -- ud_0 right
            add(new_ptr((Tree2toNDiagram(5), prt11, q21, pom, q12, prt12, 1,dq21, 2, q21, 3, q12, 4, dq12, -3)));
            
            //d -- uu_1 left d -- uu_1 right
            add(new_ptr((Tree2toNDiagram(5), prt11, q21, pom, q22, prt12, 1, dq21, 2, q21, 3, q22, 4, dq22, -4)));
          break;
            }
      
        }
        break;
      case 1: //two clusters (cases with ud_1 not included)
        switch (diffDirection){
          case 0: //left
            //u -- ud_0
            add(new_ptr((Tree2toNDiagram(5), prt11, q11, gl, pom, prt12, 1, dq11, 2, q11, 3, gl, 4, prt12, -1)));
            //d -- uu_1
            add(new_ptr((Tree2toNDiagram(5), prt11, q21, gl, pom, prt12, 1, dq21, 2, q21, 3, gl, 4, prt12, -2)));
            break;
          case 1: //right
            //u -- ud_0
            add(new_ptr((Tree2toNDiagram(5), prt11, pom, gl, q12, prt12, 1, prt11, 2, gl, 3, q12, 4, dq12, -1)));
            //d -- ud_1
            add(new_ptr((Tree2toNDiagram(5), prt11, pom, gl, q22, prt12, 1, prt11, 2, gl, 3, q22, 4, dq22, -2))); 
            break;
          case 2: //double
            //u -- ud_0 left u -- ud_0 right
            add(new_ptr((Tree2toNDiagram(7), prt11, q11, gl, pom, gl, q12, prt12, 1, dq11, 2, q11, 3, gl, 4, 
            gl, 5, q12, 6, dq12, -1)));
            //u -- ud_0 left d -- uu_1 right
            add(new_ptr((Tree2toNDiagram(7), prt11, q11, gl, pom, gl, q22, prt12, 1, dq11, 2, q11, 3, gl, 4, 
            gl, 5, q22, 6, dq22, -2)));
            //d -- uu_1 left u -- ud_0 right
            add(new_ptr((Tree2toNDiagram(7), prt11, q21, gl, pom, gl, q12, prt12, 1, dq21, 2, q21, 3, gl, 4, 
            gl, 5, q12, 6, dq12, -3)));
            //d -- uu_1 left d -- uu_1 right
            add(new_ptr((Tree2toNDiagram(7), prt11, q21, gl, pom, gl, q22, prt12, 1, dq21, 2, q21, 3, gl, 4, 
            gl, 5, q22, 6, dq22, -4)));
            break;
        } 
        break;  
    }    
}

Energy2 MEDiffraction::scale() const {
  return sqr(10*GeV);
}

int MEDiffraction::nDim() const {
  return 0;
}

void MEDiffraction::setKinematics() {
  HwMEBase::setKinematics(); // Always call the base class method first
}

bool MEDiffraction::generateKinematics(const double * ) {
  // generate the masses of the particles
  for (size_t i = 2; i < meMomenta().size(); ++i)
    meMomenta()[i] = Lorentz5Momentum(mePartonData()[i]->generateMass()); 

  /* sample M12, M22 and t, characterizing the diffractive final state */
  const pair<pair<Energy2,Energy2>,Energy2> point = diffractiveMassAndMomentumTransfer();
  const Energy2 M12 (point.first.first);
  const Energy2 M22 (point.first.second);
  const Energy2 t(point.second);
  
 
  /* construct the hadronic momenta in the lab frame */
  const double phi = UseRandom::rnd() * Constants::twopi;
  const Energy cmEnergy = generator()->maximumCMEnergy();
  const Energy2 s = sqr(cmEnergy);
  
  //proton mass
  const Energy2 m2 = sqr( theProtonMass );
  
  const Energy E3 = (s - M22 + M12) / (2.*cmEnergy);
  const Energy E4 = (s + M22 - M12) / (2.*cmEnergy);
  
  //Momentum of outgoing proton and dissociated proton
  const Energy pprime = sqrt(kallen(s, M12, M22)) / (2.*cmEnergy);
  
  //costheta of scattering angle  
  double costheta = s*(s + 2*t - 2*m2 - M12 - M22)
                        / sqrt( kallen(s, M12, M22)*kallen(s, m2, m2) );
  
  assert(abs(costheta)<=1.);
 
  const Energy pzprime = pprime*costheta;
  const Energy pperp = pprime*sqrt(1 - sqr(costheta));

  /* momenta in the lab frame */
  const Lorentz5Momentum p3 = Lorentz5Momentum(pperp*cos(phi), pperp*sin(phi), pzprime, E3);
  const Lorentz5Momentum p4 = Lorentz5Momentum(-pperp*cos(phi), -pperp*sin(phi), -pzprime, E4);
  
  /* decay dissociated proton into quark-diquark */
  //squares of constituent masses of quark and diquark
  const Energy2 mq2(sqr(mq()));
  
  Energy2 Mx2;
  switch(diffDirection){
    case 0:
      Mx2=M12;
      break;
    case 1: 
      Mx2=M22;
      break;
  }
   
    
  /* Select between left/right single diffraction and double diffraction */
  //check if we want only delta for the excited state
  
  
  //pair of momenta for double decay for a two cluster case
  pair<Lorentz5Momentum,Lorentz5Momentum> momPair, momPair1;
  //fraction of momenta
  double frac = UseRandom::rnd();
  
  switch(dissociationDecay){
    case 0:
      if(!deltaOnly)
      {
    
        pair<Lorentz5Momentum,Lorentz5Momentum> decayMomenta;
        pair<Lorentz5Momentum,Lorentz5Momentum> decayMomentaTwo;
        const double phiprime = phi;
        
        //aligned with outgoing dissociated proton
        const double costhetaprime = costheta;
        
        const double sinthetaprime=sqrt(1-sqr(costhetaprime));
        //axis along which diquark from associated proton is aligned
         Axis dir = Axis(sinthetaprime*cos(phiprime), sinthetaprime*sin(phiprime), costhetaprime);
        
        switch (diffDirection){
            case 0://Left single diffraction
              meMomenta()[4].setT(sqrt(mq2+sqr(meMomenta()[4].x())+sqr(meMomenta()[4].y())+sqr(meMomenta()[4].z())));
              ////////////////////////////////////////////////////
              
              do{}
              while(!Kinematics::twoBodyDecay(p3,mqq(),mq(),-dir,decayMomenta.first,decayMomenta.second));
              ///////////
              
              meMomenta()[2].setVect(p4.vect());
              meMomenta()[2].setT(p4.t());
        
              meMomenta()[3].setVect(decayMomenta.first.vect());
              meMomenta()[3].setT(decayMomenta.first.t());
              meMomenta()[4].setVect(decayMomenta.second.vect());
              meMomenta()[4].setT(decayMomenta.second.t());
        
              meMomenta()[2].rescaleEnergy();
              meMomenta()[3].rescaleEnergy();
              meMomenta()[4].rescaleEnergy(); 
              break;
            case 1://Right single diffraction
              meMomenta()[4].setT(sqrt(mq2+sqr(meMomenta()[4].x())+sqr(meMomenta()[4].y())+sqr(meMomenta()[4].z())));
              ////////////////////////////////////////////////////
              
              do{}
              while(!Kinematics::twoBodyDecay(p4,mqq(),mq(),dir,decayMomenta.first,decayMomenta.second));
              
              meMomenta()[2].setVect(p3.vect());
              meMomenta()[2].setT(p3.t());
        
              meMomenta()[3].setVect(decayMomenta.first.vect());
              meMomenta()[3].setT(decayMomenta.first.t());
              meMomenta()[4].setVect(decayMomenta.second.vect());
              meMomenta()[4].setT(decayMomenta.second.t());
        
              meMomenta()[2].rescaleEnergy();
              meMomenta()[3].rescaleEnergy();
              meMomenta()[4].rescaleEnergy(); 
                break;
            case 2://double diffraction
              
              do{}
              while(!Kinematics::twoBodyDecay(p3,mqq(),mq(),-dir,decayMomenta.first,decayMomenta.second));
              
              do{}
              while(!Kinematics::twoBodyDecay(p4,mqq(),mq(),dir,decayMomentaTwo.first,decayMomentaTwo.second));
              
              meMomenta()[2].setVect(decayMomenta.first.vect());
              meMomenta()[2].setT(decayMomenta.first.t());
              meMomenta()[3].setVect(decayMomenta.second.vect());
              meMomenta()[3].setT(decayMomenta.second.t());
              
              meMomenta()[4].setVect(decayMomentaTwo.second.vect());
              meMomenta()[4].setT(decayMomentaTwo.second.t());
              meMomenta()[5].setVect(decayMomentaTwo.first.vect());
              meMomenta()[5].setT(decayMomentaTwo.first.t());
              
              
              meMomenta()[2].rescaleEnergy();
              meMomenta()[3].rescaleEnergy();
              meMomenta()[4].rescaleEnergy(); 
        
              meMomenta()[5].rescaleEnergy();
        
              break;            
      }
      
      }else
      {
            const auto tmp=diffDirection==1?1:0;
            meMomenta()[2+tmp].setVect(p3.vect());
            meMomenta()[2+tmp].setT(p3.t());
            meMomenta()[3-tmp].setVect(p4.vect());
            meMomenta()[3-tmp].setT(p4.t());
    
            meMomenta()[2].rescaleEnergy();
            meMomenta()[3].rescaleEnergy();
 
      }
      break;
    case 1:
      switch(diffDirection){
        case 0: 
          //quark and diquark masses
          meMomenta()[2].setMass(mqq());
          meMomenta()[3].setMass(mq());
          
          //gluon constituent mass
          meMomenta()[4].setMass(getParticleData(21)->constituentMass());
          
          //outgoing proton
          meMomenta()[5].setVect(p4.vect());
          meMomenta()[5].setT(p4.t());
          
          //two body decay of the outgoing dissociation proton
          do{}
            while(!Kinematics::twoBodyDecay(p3,mqq()+mq(),getParticleData(21)->constituentMass(),
                p3.vect().unit(),momPair.first,momPair.second));
            //put gluon back-to-back with quark-diquark   
            //set momenta of quark and diquark
            frac = mqq()/(mqq()+mq());
            meMomenta()[2].setVect(frac*momPair.first.vect());
            meMomenta()[2].setT(sqrt(sqr(frac)*momPair.first.vect().mag2()+sqr(mqq())));
            meMomenta()[3].setVect((1-frac)*momPair.first.vect());
            meMomenta()[3].setT(sqrt(sqr(1-frac)*momPair.first.vect().mag2()+sqr(mq())));
            //set momentum of gluon
            meMomenta()[4].setVect(momPair.second.vect());
            meMomenta()[4].setT(momPair.second.t());
                    
          break;
        case 1: 
          //quark and diquark masses
          meMomenta()[5].setMass(mqq());
          meMomenta()[4].setMass(mq());
          
          //gluon constituent mass
          meMomenta()[3].setMass(getParticleData(21)->constituentMass());
          
          //outgoing proton
          meMomenta()[2].setVect(p3.vect());
          meMomenta()[2].setT(p3.t());
          
          //two body decay of the outgoing dissociation proton
          do{}
            while(!Kinematics::twoBodyDecay(p4,mqq()+mq(),getParticleData(21)->constituentMass(),
                p4.vect().unit(),momPair.first,momPair.second));
            
            //put gluon back-to-back with quark-diquark   
            //set momenta of quark and diquark
            frac = mqq()/(mqq()+mq());
            meMomenta()[5].setVect(frac*momPair.first.vect());
            meMomenta()[5].setT(sqrt(sqr(frac)*momPair.first.vect().mag2()+sqr(mqq())));
            meMomenta()[4].setVect((1-frac)*momPair.first.vect());
            meMomenta()[4].setT(sqrt(sqr(1-frac)*momPair.first.vect().mag2()+sqr(mq())));
            //set momentum of gluon
            meMomenta()[3].setVect(momPair.second.vect());
            meMomenta()[3].setT(momPair.second.t());
          
          
        
          break;
        case 2: 
          //first dissociated proton constituents
          meMomenta()[2].setMass(mqq());
          meMomenta()[3].setMass(mq());
          meMomenta()[4].setMass(getParticleData(21)->constituentMass());
          //second dissociated proton constituents
          meMomenta()[5].setMass(getParticleData(21)->constituentMass());
          meMomenta()[6].setMass(mq());
          meMomenta()[7].setMass(mqq());
          
          
          //two body decay of the outgoing dissociation proton
          do{}
            while(!Kinematics::twoBodyDecay(p3,mqq()+mq(),getParticleData(21)->constituentMass(),
                p3.vect().unit(),momPair.first,momPair.second));
            
            do{}
            while(!Kinematics::twoBodyDecay(p4,mqq()+mq(),getParticleData(21)->constituentMass(),
                p4.vect().unit(),momPair1.first,momPair1.second));    
          
          //put gluon back-to-back with quark-diquark
          frac = mqq()/(mqq()+mq());
          
          //first dissociated proton
          //set momenta of quark and diquark
            
            meMomenta()[2].setVect(frac*momPair.first.vect());
            meMomenta()[2].setT(sqrt(sqr(frac)*momPair.first.vect().mag2()+sqr(mqq())));
            meMomenta()[3].setVect((1-frac)*momPair.first.vect());
            meMomenta()[3].setT(sqrt(sqr(1-frac)*momPair.first.vect().mag2()+sqr(mq())));
            //set momentum of gluon
            meMomenta()[4].setVect(momPair.second.vect());
            meMomenta()[4].setT(momPair.second.t());
            
            //first dissociated proton
          //set momenta of quark and diquark
            
            meMomenta()[7].setVect(frac*momPair1.first.vect());
            meMomenta()[7].setT(sqrt(sqr(frac)*momPair1.first.vect().mag2()+sqr(mqq())));
            meMomenta()[6].setVect((1-frac)*momPair1.first.vect());
            meMomenta()[6].setT(sqrt(sqr(1-frac)*momPair1.first.vect().mag2()+sqr(mq())));
            //set momentum of gluon
            meMomenta()[5].setVect(momPair1.second.vect());
            meMomenta()[5].setT(momPair1.second.t());
          break;
        
      }
      meMomenta()[2].rescaleEnergy();
      meMomenta()[3].rescaleEnergy();
      meMomenta()[4].rescaleEnergy();
      meMomenta()[5].rescaleEnergy();
      if(diffDirection==2){
        meMomenta()[6].rescaleEnergy();
        meMomenta()[7].rescaleEnergy();
      }
      
      break;
  }
  
  jacobian(sqr(cmEnergy)/GeV2);
  return true;
}
//Generate masses of dissociated protons and momentum transfer from probability f(M2,t) 
//(for single diffraction). Sample according to f(M2,t)=f(M2)f(t|M2).
pair<pair<Energy2,Energy2>,Energy2> MEDiffraction::diffractiveMassAndMomentumTransfer() const {
  Energy2 theM12(ZERO),theM22(ZERO), thet(ZERO);
  int count = 0;
  //proton mass squared
  const Energy2 m2 = sqr(theProtonMass);
  //delta mass squared
  const Energy2 md2 = sqr(getParticleData(2214)->mass());
  Energy2 M2;
  bool condition = true;
  do {  
    
    //check if we want only delta 
    if(deltaOnly) {
      switch(diffDirection){
        case 0:
          theM12 = md2;
          theM22 = m2;
          M2 = md2;
          thet = randomt(md2);
          break;
        case 1:
          theM22 = md2;
          theM12 = m2;
          M2 = md2;   
          thet = randomt(md2);
          break;
        case 2:
          theM12 = md2;
          theM22 = md2;
          M2 = md2;
          thet = doublediffrandomt(theM12,theM22);
          break;  
      }

    }
    else {
      switch (diffDirection){
      case 0:
        M2=randomM2();
        thet = randomt(M2);
        theM12=M2;
        
        theM22=m2;
        
        break;
      case 1:
        
        theM12=m2;
        M2=randomM2();
        thet = randomt(M2);
        
        
        theM22=M2; 
        break;
      case 2:
        theM12=randomM2();
        theM22=randomM2();
        M2=(theM12>theM22) ? theM12: theM22;
        
        thet = doublediffrandomt(theM12,theM22);
        
        break;
      }
    }
    count++;
  
  const Energy cmEnergy = generator()->maximumCMEnergy();
  const Energy2 s = sqr(cmEnergy);
    if(generator()->maximumCMEnergy()<sqrt(theM12)+sqrt(theM22)) {
      condition = true;
    }
    else {
  InvEnergy2 slope;
  if(diffDirection==2){
    slope = 2*softPomeronSlope()*log(.1+(sqr(cmEnergy)/softPomeronSlope())/(theM12*theM22));
  }else{
    slope = protonPomeronSlope()
                         + 2*softPomeronSlope()*log(sqr(cmEnergy)/M2);
  }

  
  const double expmax = exp(slope*tmaxfun(s,m2,M2));
  const double expmin = exp(slope*tminfun(s,m2,M2));
  
  
  //without (1-M2/s) constraint
  condition = (UseRandom::rnd()>(protonPomeronSlope()*GeV2)*(expmax-expmin)/(slope*GeV2))
      ||((theM12/GeV2)*(theM22/GeV2)>=(sqr(cmEnergy)/GeV2)/(softPomeronSlope()*GeV2));
    }
  }
  while(condition);
  
  return make_pair (make_pair(theM12,theM22),thet);
}

//Decay of the excited proton to quark-diquark
pair<Lorentz5Momentum,Lorentz5Momentum> MEDiffraction::twoBodyDecayMomenta(Lorentz5Momentum pp) const{
  //Decay of the excited proton
    const Energy2 Mx2(sqr(pp.mass())),mq2(sqr(mq())),mqq2(sqr(mqq()));
    
        const Energy2 psq = ((Mx2-sqr(mq()+mqq()))*(Mx2-sqr(mq()-mqq())))/(4*Mx2);    

    assert(psq/GeV2>0);
    const Energy p(sqrt(psq));
    
    const double phi = UseRandom::rnd() * Constants::twopi;
    const double costheta =1-2*UseRandom::rnd();
    const double sintheta = sqrt(1-sqr(costheta));
    
    Lorentz5Momentum k1=Lorentz5Momentum(p*sintheta*cos(phi), p*sintheta*sin(phi), p*costheta, sqrt(mq2+psq));
    Lorentz5Momentum k2=Lorentz5Momentum(-p*sintheta*cos(phi), -p*sintheta*sin(phi), -p*costheta,sqrt(mqq2+psq));
    
    //find boost to pp center of mass
    const Boost betap3 = (pp).findBoostToCM();
  
    //k1 and k2 calculated at p3 center of mass, so boost back
    k1.boost(-betap3);
    k2.boost(-betap3);
    
    //first is quark, second diquark
    return make_pair(k1,k2);
}



Energy2 MEDiffraction::randomt(Energy2 M2) const {
  assert(protonPomeronSlope()*GeV2 > 0);
  //proton mass
  const Energy2 m2 = sqr( theProtonMass );
  const Energy cmEnergy = generator()->maximumCMEnergy();
  const Energy2 ttmin = tminfun(sqr(cmEnergy),m2,M2);
  const Energy2 ttmax = tmaxfun(sqr(cmEnergy),m2,M2);

  const InvEnergy2 slope = protonPomeronSlope()
                         + 2*softPomeronSlope()*log(sqr(cmEnergy)/M2);
    return log( exp(slope*ttmin) +
              UseRandom::rnd()*(exp(slope*ttmax) - exp(slope*ttmin)) ) / slope;
}

Energy2 MEDiffraction::doublediffrandomt(Energy2 M12, Energy2 M22) const {
  
  const Energy cmEnergy = generator()->maximumCMEnergy();
  
  const double shift = 0.1;
  const InvEnergy2 slope = 2*softPomeronSlope()*log(shift+(sqr(cmEnergy)/softPomeronSlope())/(M12*M22));
  
  const Energy2 ttmin = tminfun(sqr(cmEnergy),M12,M22);
  const Energy2 ttmax = tmaxfun(sqr(cmEnergy),M12,M22);
  double r = UseRandom::rnd();
  Energy2 newVal;
  if(slope*ttmax>slope*ttmin) {
    newVal = ttmax + log( r + (1.-r)*exp(slope*(ttmin-ttmax)) ) / slope;
  }
  else {
    newVal = ttmin + log( 1. - r + r*exp(slope*(ttmax-ttmin))) / slope;
  }
  return newVal;
}

Energy2 MEDiffraction::randomM2() const {
  const double tmp = 1 - softPomeronIntercept();
  
  const Energy cmEnergy = generator()->maximumCMEnergy();
  return sqr(cmEnergy) * pow( pow(M2min()/sqr(cmEnergy),tmp) +
                   UseRandom::rnd() * (pow(M2max()/sqr(cmEnergy),tmp) - pow(M2min()/sqr(cmEnergy),tmp)),
                   1.0/tmp );
}

Energy2 MEDiffraction::tminfun(Energy2 s, Energy2 M12, Energy2 M22) const {
  const Energy2 m2 = sqr( theProtonMass );
  return 0.5/s*(-sqrt(kallen(s, m2, m2)*kallen(s, M12, M22))-sqr(s)+2*s*m2+s*M12+s*M22);
}

Energy2 MEDiffraction::tmaxfun(Energy2 s, Energy2 M12, Energy2 M22) const  {
  const Energy2 m2 = sqr( theProtonMass );
  
  return 0.5/s*(sqrt(kallen(s, m2, m2)*kallen(s, M12, M22))-sqr(s)+2*s*m2+s*M12+s*M22);
}


double MEDiffraction::me2() const{
  return theme2; 
}

CrossSection MEDiffraction::dSigHatDR() const {
  return me2()*jacobian()/sHat()*sqr(hbarc);
}

unsigned int MEDiffraction::orderInAlphaS() const {
  return 0;
}

unsigned int MEDiffraction::orderInAlphaEW() const {
  return 0;
}

Selector<MEBase::DiagramIndex>
MEDiffraction::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
    if(!deltaOnly){
      if(diffDirection<2){
        
        for(unsigned int i = 0; i < diags.size(); i++){
          if(diags[0]->id()==-1) 
            sel.insert(2./3.,i);
          else
            sel.insert(1./3.,i);  
        }
        
      }else{
        for(unsigned int i = 0; i < diags.size(); i++){
          if(diags[0]->id()==-1) 
            sel.insert(4./9.,i);
          else if(diags[0]->id()==-2)
            sel.insert(2./9.,i);  
          else if(diags[0]->id()==-3)
            sel.insert(2./9.,i);
          else
            sel.insert(1./9.,i);    
        }
      }
    }else{
      sel.insert(1.0,0);
    }  
  
  return sel;
}

Selector<const ColourLines *>
MEDiffraction::colourGeometries(tcDiagPtr ) const {
  Selector<const ColourLines *> sel;
  
  int sign1=0, sign2=0;
  sign1 = (generator()->eventHandler()->incoming().first->id() > 0) ? 1 : -1;
  sign2 = (generator()->eventHandler()->incoming().second->id() > 0) ? 1 : -1;
  
  switch(dissociationDecay){
    case 0:
      if(!deltaOnly)
      { 
        if(diffDirection!=2){
          
          if (diffDirection == 0){
            if(sign1>0){
              static ColourLines dqq0=ColourLines("-6 2 7");
              sel.insert(1.0,&dqq0);
            }else{
              static ColourLines dqq0=ColourLines("6 -2 -7");
              sel.insert(1.0,&dqq0);
            }
            
            
          }
          else{
            if(sign2>0){
              static ColourLines dqq1=ColourLines("-6 3 7");
              sel.insert(1.0,&dqq1);
            }else{
              static ColourLines dqq1=ColourLines("6 -3 -7");
              sel.insert(1.0,&dqq1);
            }
          }
          
        }else{
          
          if(sign1>0 && sign2>0){
            static ColourLines ddqq0=ColourLines("-6 2 7, -9 4 8");
            sel.insert(1.0,&ddqq0);
          }else if(sign1<0 && sign2>0){
            static ColourLines ddqq0=ColourLines("6 -2 -7, -9 4 8");
            sel.insert(1.0,&ddqq0);
          }else if(sign1>0&& sign2<0){
            static ColourLines ddqq0=ColourLines("-6 2 7, 9 -4 -8");
            sel.insert(1.0,&ddqq0);
          }else{
            static ColourLines ddqq0=ColourLines("6 -2 -7, 9 -4 -8");
            sel.insert(1.0,&ddqq0);
          }
          
          
        }
    
      }else
      {
      static ColourLines cl("");
      
      sel.insert(1.0, &cl);    
      }
      break;
    case 1:
      switch(diffDirection){
        case 0: 
          static ColourLines clleft("-6 2 3 8, -8 -3 7");
          sel.insert(1.0, &clleft);
          break;
        case 1: 
          static ColourLines clright("-9 4 3 7, -7 -3 8");
          sel.insert(1.0, &clright);
          break;
        case 2: 
          static ColourLines cldouble("-8 2 3 10, -10 -3 9, -13 6 5 11, -11 -5 12");
          sel.insert(1.0, &cldouble);
          break;
      }
      break;
  }
  
  return sel;
}

void MEDiffraction::doinit() {
  HwMEBase::doinit();
  theProtonMass = getParticleData(2212)->mass();
}

void MEDiffraction::doinitrun() {
  HwMEBase::doinitrun();
  isInRunPhase = true;
}

IBPtr MEDiffraction::clone() const {
  return new_ptr(*this);
}

IBPtr MEDiffraction::fullclone() const {
  return new_ptr(*this);
}


ClassDescription<MEDiffraction> MEDiffraction::initMEDiffraction;
// Definition of the static class description member.

void MEDiffraction::persistentOutput(PersistentOStream & os) const {
  os << theme2 << deltaOnly << diffDirection << theprotonPomeronSlope
     << thesoftPomeronIntercept << thesoftPomeronSlope << dissociationDecay
     << ounit(theProtonMass,GeV);
}

void MEDiffraction::persistentInput(PersistentIStream & is, int) {
  is >> theme2 >> deltaOnly >> diffDirection >> theprotonPomeronSlope
     >> thesoftPomeronIntercept >> thesoftPomeronSlope >> dissociationDecay
     >> iunit(theProtonMass,GeV);
}

InvEnergy2 MEDiffraction::protonPomeronSlope() const{
  return theprotonPomeronSlope/GeV2;
}

double MEDiffraction::softPomeronIntercept() const {
  return thesoftPomeronIntercept;
}

InvEnergy2 MEDiffraction::softPomeronSlope() const {
  return thesoftPomeronSlope/GeV2;
}

void MEDiffraction::Init() {

  static ClassDocumentation<MEDiffraction> documentation
    ("There is no documentation for the MEDiffraction class");
  
  static Parameter<MEDiffraction,double> interfaceme2
    ("DiffractionAmplitude",
     "The square of the diffraction amplitude used to determine the "
     "cross section.",
     &MEDiffraction::theme2, 1.0, 0.00001, 100.0,
     false, false, Interface::limited);
   
   static Parameter<MEDiffraction,double> interfaceprotonPomeronSlope
    ("ProtonPomeronSlope",
     "The proton-pomeron slope parameter.",
     &MEDiffraction::theprotonPomeronSlope, 10.1, 0.00001, 100.0,
     false, false, Interface::limited);    
   
   static Parameter<MEDiffraction,double> interfacesoftPomeronIntercept
    ("SoftPomeronIntercept",
     "The soft pomeron intercept.",
     &MEDiffraction::thesoftPomeronIntercept, 1.08, 0.00001, 100.0,
     false, false, Interface::limited);    
   
   static Parameter<MEDiffraction,double> interfacesoftPomeronSlope
    ("SoftPomeronSlope",
     "The soft pomeron slope parameter.",
     &MEDiffraction::thesoftPomeronSlope, 0.25, 0.00001, 100.0,
     false, false, Interface::limited);      
     
  
  static Switch<MEDiffraction, bool> interfaceDeltaOnly
    ("DeltaOnly",
     "proton-proton to proton-delta only",
     &MEDiffraction::deltaOnly, 0, false, false);
  static SwitchOption interfaceDeltaOnly0
    (interfaceDeltaOnly,"No","Final state with Delta only is OFF", 0);
  static SwitchOption interfaceDeltaOnly1
    (interfaceDeltaOnly,"Yes","Final state with Delta only is ON", 1);  
    
    //Select if the left, right or both protons are excited
  static Switch<MEDiffraction, unsigned int> interfaceDiffDirection
    ("DiffDirection",
     "Direction of the excited proton",
     &MEDiffraction::diffDirection, 0, false, false);
  static SwitchOption left
    (interfaceDiffDirection,"Left","Proton moving in the positive z direction", 0);
  static SwitchOption right
    (interfaceDiffDirection,"Right","Proton moving in the negative z direction", 1);
  static SwitchOption both
    (interfaceDiffDirection,"Both","Both protons", 2);
  
  //Select if two or three body decay
  static Switch<MEDiffraction, unsigned int> interfaceDissociationDecay
    ("DissociationDecay",
     "Number of clusters the dissociated proton decays",
     &MEDiffraction::dissociationDecay, 0, false, false);
  static SwitchOption one
    (interfaceDissociationDecay,"One","Dissociated proton decays into one cluster", 0);
  static SwitchOption two
    (interfaceDissociationDecay,"Two","Dissociated proton decays into two clusters", 1);
} 
