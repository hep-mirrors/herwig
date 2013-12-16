#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"
#include "DC_Hjj_NC_qq.h"

using namespace Herwig;

IBPtr DC_Hjj_NC_qq::clone() const {
  return new_ptr(*this);
}
IBPtr DC_Hjj_NC_qq::fullclone() const {
  return new_ptr(*this);
}

vector<DiagPtr> DC_Hjj_NC_qq::getDiagrams() const {
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
      PDPtr Z0 = theME->getParticleData(ParticleID::Z0); 
      PDPtr h0 = theME->getParticleData(ParticleID::h0); 
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, kd, 1, qd, 3, kd, 2, h0, -1)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, ku, 1, qd, 3, ku, 2, h0, -2)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, kd, 1, qu, 3, kd, 2, h0, -3)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, ku, 1, qu, 3, ku, 2, h0, -4)));
    }
  }
  return allPossibleDiagrams;
}



Selector<const ColourLines *> DC_Hjj_NC_qq::colourGeometries(tcDiagPtr diag) const {
  // colour lines for DC_Hjj_NC_qq

  static const ColourLines diag1[1] = { 
    ColourLines("1 5, 4 6")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("1 5, 4 6")
  }; 
  static const ColourLines diag3[1] = { 
    ColourLines("1 5, 4 6")
  }; 
  static const ColourLines diag4[1] = { 
    ColourLines("1 5, 4 6")
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
