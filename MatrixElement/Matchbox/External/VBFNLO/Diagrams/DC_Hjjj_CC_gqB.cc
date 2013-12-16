#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"
#include "DC_Hjjj_CC_gqB.h"

using namespace Herwig;

IBPtr DC_Hjjj_CC_gqB::clone() const {
  return new_ptr(*this);
}
IBPtr DC_Hjjj_CC_gqB::fullclone() const {
  return new_ptr(*this);
}

vector<DiagPtr> DC_Hjjj_CC_gqB::getDiagrams() const {
  // get the particle data objects

  vector<DiagPtr> allPossibleDiagrams;
  for (int genU = 0; genU < 2; genU++) {
    for (int genL = 0; genL < 2; genL++) {
      PDPtr qu, qd, ku, kd;
      if (genU == 0){
        qu = theME->getParticleData(ParticleID::u);
        qd = theME->getParticleData(ParticleID::d);
      }
      else if (genU == 1){
        qu = theME->getParticleData(ParticleID::c);
        qd = theME->getParticleData(ParticleID::s);
      }
      else if (genU == 2){
        qu = theME->getParticleData(ParticleID::t);
        qd = theME->getParticleData(ParticleID::b);
      }
      if (genL == 0){
        ku = theME->getParticleData(ParticleID::u);
        kd = theME->getParticleData(ParticleID::d);
      }
      else if (genL == 1){
        ku = theME->getParticleData(ParticleID::c);
        kd = theME->getParticleData(ParticleID::s);
      }
      else if (genL == 2){
        ku = theME->getParticleData(ParticleID::t);
        kd = theME->getParticleData(ParticleID::b);
      }

      tcPDPtr quB = (*qu).CC();
      tcPDPtr qdB = (*qd).CC();

      tcPDPtr quC = SU2Helper::SU2CC(qu);
      tcPDPtr qdC = SU2Helper::SU2CC(qd);
      tcPDPtr quBC = SU2Helper::SU2CC(quB);
      tcPDPtr qdBC = SU2Helper::SU2CC(qdB);

      tcPDPtr kuB = (*ku).CC();
      tcPDPtr kdB = (*kd).CC();

      tcPDPtr kuC = SU2Helper::SU2CC(ku);
      tcPDPtr kdC = SU2Helper::SU2CC(kd);
      tcPDPtr kuBC = SU2Helper::SU2CC(kuB);
      tcPDPtr kdBC = SU2Helper::SU2CC(kdB);
      PDPtr Wplus = theME->getParticleData(ParticleID::Wplus); 
      PDPtr Wminus = theME->getParticleData(ParticleID::Wminus); 
      PDPtr g = theME->getParticleData(ParticleID::g); 
      PDPtr h0 = theME->getParticleData(ParticleID::h0); 
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qdBC, Wplus, Wplus, kdB, 1, qdC, 4, kdBC, 2, qdB, 3, h0, -1)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qd, Wplus, Wplus, kdB, 2, qdC, 4, kdBC, 1, qdB, 3, h0, -2)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, quBC, Wminus, Wminus, kuB, 1, quC, 4, kuBC, 2, quB, 3, h0, -3)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), g, qu, Wminus, Wminus, kuB, 2, quC, 4, kuBC, 1, quB, 3, h0, -4)));
    }
  }
  return allPossibleDiagrams;
}



Selector<const ColourLines *> DC_Hjjj_CC_gqB::colourGeometries(tcDiagPtr diag) const {
  // colour lines for DC_Hjjj_CC_gqB

  static const ColourLines diag1[1] = { 
    ColourLines("1 6, -1 -2 -8, -5 -7")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("1 2 6, -1 -8, -5 -7")
  }; 
  static const ColourLines diag3[1] = { 
    ColourLines("1 6, -1 -2 -8, -5 -7")
  }; 
  static const ColourLines diag4[1] = { 
    ColourLines("1 2 6, -1 -8, -5 -7")
  }; 

  Selector <const ColourLines *> sel;

  if( diag->id() == -1 )  {
   sel.insert( 1.0,  &(diag1[0]) );
  } 
  else if( diag->id() == -2 )  {
   sel.insert( 1.0,  &(diag2[0]) );
  } 
  else if( diag->id() == -3 )  {
   sel.insert( 1.0,  &(diag3[0]) );
  } 
  else if( diag->id() == -4 )  {
   sel.insert( 1.0,  &(diag4[0]) );
  } 
  return sel;
}
