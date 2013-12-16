#include "Herwig++/MatrixElement/Matchbox/Utility/SU2Helper.h"
#include "DC_Hjjj_NC_qq.h"

using namespace Herwig;

IBPtr DC_Hjjj_NC_qq::clone() const {
  return new_ptr(*this);
}
IBPtr DC_Hjjj_NC_qq::fullclone() const {
  return new_ptr(*this);
}

vector<DiagPtr> DC_Hjjj_NC_qq::getDiagrams() const {
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
      PDPtr g = theME->getParticleData(ParticleID::g); 
      PDPtr h0 = theME->getParticleData(ParticleID::h0); 
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Z0, Z0, kdB, kd, 1, qd, 3, kd, 4, g, 2, h0, -1)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, kd, 3, kd, 1, qd, 5, kd, 5, g, 2, h0, -2)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, qd, Z0, Z0, kd, 2, qd, 4, kd, 1, g, 3, h0, -3)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, kd, 1, qd, 5, qd, 3, kd, 5, g, 2, h0, -4)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, Z0, Z0, kuB, ku, 1, qd, 3, ku, 4, g, 2, h0, -5)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, ku, 3, ku, 1, qd, 5, ku, 5, g, 2, h0, -6)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qd, qd, Z0, Z0, ku, 2, qd, 4, ku, 1, g, 3, h0, -7)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qd, Z0, Z0, ku, 1, qd, 5, qd, 3, ku, 5, g, 2, h0, -8)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Z0, Z0, kdB, kd, 1, qu, 3, kd, 4, g, 2, h0, -9)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, kd, 3, kd, 1, qu, 5, kd, 5, g, 2, h0, -10)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, qu, Z0, Z0, kd, 2, qu, 4, kd, 1, g, 3, h0, -11)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, kd, 1, qu, 5, qu, 3, kd, 5, g, 2, h0, -12)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, Z0, Z0, kuB, ku, 1, qu, 3, ku, 4, g, 2, h0, -13)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, ku, 3, ku, 1, qu, 5, ku, 5, g, 2, h0, -14)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(5), qu, qu, Z0, Z0, ku, 2, qu, 4, ku, 1, g, 3, h0, -15)));
      allPossibleDiagrams.push_back(new_ptr((Tree2toNDiagram(4), qu, Z0, Z0, ku, 1, qu, 5, qu, 3, ku, 5, g, 2, h0, -16)));
    }
  }
  return allPossibleDiagrams;
}



Selector<const ColourLines *> DC_Hjjj_NC_qq::colourGeometries(tcDiagPtr diag) const {
  // colour lines for DC_Hjjj_NC_qq

  static const ColourLines diag1[1] = { 
    ColourLines("1 6, 5 8, 7 -4 -8")
  }; 
  static const ColourLines diag2[1] = { 
    ColourLines("1 6, 4 5 8, 7 -8")
  }; 
  static const ColourLines diag3[1] = { 
    ColourLines("1 8, 5 7, 6 2 -8")
  }; 
  static const ColourLines diag4[1] = { 
    ColourLines("1 5 8, 4 7, 6 -8")
  }; 
  static const ColourLines diag5[1] = { 
    ColourLines("1 6, 5 8, 7 -4 -8")
  }; 
  static const ColourLines diag6[1] = { 
    ColourLines("1 6, 4 5 8, 7 -8")
  }; 
  static const ColourLines diag7[1] = { 
    ColourLines("1 8, 5 7, 6 2 -8")
  }; 
  static const ColourLines diag8[1] = { 
    ColourLines("1 5 8, 4 7, 6 -8")
  }; 
  static const ColourLines diag9[1] = { 
    ColourLines("1 6, 5 8, 7 -4 -8")
  }; 
  static const ColourLines diag10[1] = { 
    ColourLines("1 6, 4 5 8, 7 -8")
  }; 
  static const ColourLines diag11[1] = { 
    ColourLines("1 8, 5 7, 6 2 -8")
  }; 
  static const ColourLines diag12[1] = { 
    ColourLines("1 5 8, 4 7, 6 -8")
  }; 
  static const ColourLines diag13[1] = { 
    ColourLines("1 6, 5 8, 7 -4 -8")
  }; 
  static const ColourLines diag14[1] = { 
    ColourLines("1 6, 4 5 8, 7 -8")
  }; 
  static const ColourLines diag15[1] = { 
    ColourLines("1 8, 5 7, 6 2 -8")
  }; 
  static const ColourLines diag16[1] = { 
    ColourLines("1 5 8, 4 7, 6 -8")
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
  else if( diag->id() == -5 )  {
   sel.insert( 1.0,  &(diag5[0]) );
  } 
  else if( diag->id() == -6 )  {
   sel.insert( 1.0,  &(diag6[0]) );
  } 
  else if( diag->id() == -7 )  {
   sel.insert( 1.0,  &(diag7[0]) );
  } 
  else if( diag->id() == -8 )  {
   sel.insert( 1.0,  &(diag8[0]) );
  } 
  else if( diag->id() == -9 )  {
   sel.insert( 1.0,  &(diag9[0]) );
  } 
  else if( diag->id() == -10 )  {
   sel.insert( 1.0,  &(diag10[0]) );
  } 
  else if( diag->id() == -11 )  {
   sel.insert( 1.0,  &(diag11[0]) );
  } 
  else if( diag->id() == -12 )  {
   sel.insert( 1.0,  &(diag12[0]) );
  } 
  else if( diag->id() == -13 )  {
   sel.insert( 1.0,  &(diag13[0]) );
  } 
  else if( diag->id() == -14 )  {
   sel.insert( 1.0,  &(diag14[0]) );
  } 
  else if( diag->id() == -15 )  {
   sel.insert( 1.0,  &(diag15[0]) );
  } 
  else if( diag->id() == -16 )  {
   sel.insert( 1.0,  &(diag16[0]) );
  } 
  return sel;
}
