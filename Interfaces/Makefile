#LIBOBJECT := libHwInterfaces.so
LIBOBJECTS := libHwAmegic.so libHwKtJet.so

SOURCES := $(wildcard *.cc)

OBJECTS := $(SOURCES:.cc=.o)

#PROGS   := $(SOURCES:.cc=.exe)

HWLIBS = -lHw -lHwHadronization -lHwShower -lHwDecay

EXTRALIBS = -lAmegic -lAmplitude -lAmplTools -lISR -lModel -lPhasespace -lProcess -lString -lVector -lZfunctions

AMEGICPATH = /usera/stephens/projects/AMEGIC++-1.0

KTJETPATH = /var/nwork/pcck/gieseke/mc/ktjet

EXTRAINCLUDES = -I$(AMEGICPATH)/Include -I$(KTJETPATH)

.PHONY: Lib

%.o: %.cc
	$(CXX) -c -ansi $(CXXEXT) $(COMPILER_FLAG) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(DEFS) $(INCLUDES) $(EXTRAINCLUDES) $<

include ../Config/Makefile.common

CXXDEBUGFLAG = -g

lib: $(LIBOBJECTS)

%.exe:

depend:
	$(SRCDIR)/Config/g++dep $(INCLUDES) $(EXTRAINCLUDES) $(SOURCES) $(CHECKSOURCES)

libHwAmegic.so: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(CXXDLFLAGS) AmegicInterface.o -o $@

libHwKtJet.so: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(CXXDEBUGFLAG) $(CXXOPTFLAG) $(CXXDLFLAGS) KtJetInterface.o -o $@

# DO NOT DELETE THIS LINE -- g++dep uses it.
# DO NOT PUT ANYTHING AFTER THIS LINE, IT WILL GO AWAY.


AmegicInterface.o: AmegicInterface.cc \
	 \
	 ../../Herwig++/Interfaces/AmegicInterface.h \
	 ../../Pythia7/Interface/Interfaced.h \
	 ../../Pythia7/Config/Pythia7.h \
	 ../../Pythia7/Config/config.h \
	 ../../Pythia7/Pointer/Ptr.h \
	 ../../Pythia7/Pointer/RCPtr.fh \
	 ../../Pythia7/Pointer/PtrTraits.h \
	 ../../Pythia7/Pointer/RCPtr.h \
	 ../../Pythia7/Pointer/ReferenceCounted.h \
	 ../../Pythia7/Pointer/ReferenceCounted.icc \
	 ../../Pythia7/Pointer/RCPtr.icc \
	 ../../Pythia7/Utilities/Rebinder.fh \
	 ../../Pythia7/Utilities/Interval.fh \
	 ../../Pythia7/Utilities/ClassDescription.fh \
	 ../../Pythia7/Interface/InterfaceBase.fh \
	 ../../Pythia7/Persistency/PersistentOStream.fh \
	 ../../Pythia7/Persistency/PersistentIStream.fh \
	 ../../Pythia7/Config/Units.h \
	 ../../Pythia7/CLHEPWrap/Lorentz5Vector.fh \
	 ../../Pythia7/CLHEPWrap/LorentzVector.fh \
	 ../../Pythia7/CLHEPWrap/ThreeVector.fh \
	 ../../Pythia7/CLHEPWrap/Transverse.fh \
	 ../../Pythia7/CLHEPWrap/SystemOfUnits.h \
	 ../../Pythia7/CLHEPWrap/CLHEP.h \
	 ../../Pythia7/Config/Constants.h \
	 ../../Pythia7/Config/std.h \
	 ../../Pythia7/Config/Pointers.h \
	 ../../Pythia7/Config/Containers.h \
	 ../../Pythia7/Utilities/FixedSizeAllocator.h \
	 ../../Pythia7/Utilities/Debug.h \
	 ../../Pythia7/Utilities/Debug.icc \
	 ../../Pythia7/Utilities/FixedSizeAllocator.icc \
	 ../../Pythia7/Interface/InterfacedBase.h \
	 ../../Pythia7/Misc/Named.h \
	 ../../Pythia7/Misc/Named.icc \
	 ../../Pythia7/Utilities/ClassDescription.h \
	 ../../Pythia7/Utilities/ClassTraits.h \
	 ../../Pythia7/Utilities/DescriptionList.h \
	 ../../Pythia7/Utilities/DescriptionList.icc \
	 ../../Pythia7/Utilities/ClassDescription.icc \
	 ../../Pythia7/Utilities/ClassDescription.tcc \
	 ../../Pythia7/Utilities/HoldFlag.h \
	 ../../Pythia7/Utilities/HoldFlag.icc \
	 ../../Pythia7/Interface/InterfacedBase.xh \
	 ../../Pythia7/Interface/InterfaceBase.h \
	 ../../Pythia7/Interface/InterfaceBase.xh \
	 ../../Pythia7/Misc/Exception.h \
	 ../../Pythia7/Misc/Exception.fh \
	 ../../Pythia7/Utilities/sstream.h \
	 ../../Pythia7/Misc/Exception.icc \
	 ../../Pythia7/Interface/InterfaceBase.icc \
	 ../../Pythia7/Interface/InterfacedBase.icc \
	 ../../Pythia7/Interface/Interfaced.icc \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Amegic.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/MyString.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Flavour.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Message.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Process.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Topology.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Point.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Vertex.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/MyComplex.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Polarisation.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Basic_Sfuncs.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Vector.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Distribution_Functions.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Pfunc.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Single_Amplitude.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Zfunc_Generator.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Prop_Generator.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Color_Generator.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Building_Blocks.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Kabbala.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Composite_Zfuncs.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Elementary_Zfuncs.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/String_Generator.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/String_Tree.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Values.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Mathematica_Interface.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/String_Handler.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Helicity.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Amplitude_Handler.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Amplitude_Generator.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Color.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Amplitude_Base.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Phase_Space_Handler.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Channel.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Channel_Elements.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Channel_Basics.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Selector.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Phase_Space_Integrator.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/ISR_Base.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Poincare.H \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Extra_XS.H \
	 ../../Herwig++/Config/Herwig.h \
	 ../../Pythia7/Repository/Strategy.fh \
	 ../../Herwig++/Interfaces/AmegicInterface.icc \
	 ../../Pythia7/Repository/EventGenerator.h \
	 ../../Pythia7/Repository/EventGenerator.fh \
	 ../../Pythia7/Repository/RandomGenerator.h \
	 ../../Pythia7/CLHEPWrap/RandomEngine.fh \
	 ../../Pythia7/Repository/RandomGenerator.icc \
	 ../../Pythia7/Repository/UseRandom.h \
	 ../../Pythia7/Repository/UseRandom.icc \
	 ../../Pythia7/Repository/CurrentGenerator.fh \
	 ../../Pythia7/Repository/EventGenerator.xh \
	 ../../Pythia7/Repository/EventGenerator.icc \
	 ../../Pythia7/Repository/EventGenerator.tcc \
	 ../../Herwig++/Utilities/HwDebug.h \
	 ../../Herwig++/Utilities/SmplHist.h \
	 ../../Herwig++/Utilities/SmplStat.h \
	 ../../Pythia7/Interface/ClassDocumentation.h \
	 ../../Pythia7/Interface/ClassDocumentation.fh \
	 ../../Pythia7/Interface/ClassDocumentation.icc \
	 ../../Pythia7/Persistency/PersistentOStream.h \
	 ../../Pythia7/Persistency/PersistentOStream.xh \
	 ../../Pythia7/Persistency/PersistentOStream.icc \
	 ../../Pythia7/Persistency/PersistentIStream.h \
	 ../../Pythia7/Persistency/InputDescription.h \
	 ../../Pythia7/Persistency/InputDescription.icc \
	 ../../Pythia7/Persistency/PersistentIStream.icc \
	 ../../Pythia7/Persistency/PersistentIStream.tcc \
	 ../../Pythia7/EventRecord/Particle.h \
	 ../../Pythia7/EventRecord/EventConfig.h \
	 ../../Pythia7/CLHEPWrap/Lorentz5Vector.h \
	 ../../Pythia7/CLHEPWrap/LorentzVector.h \
	 ../../Pythia7/CLHEPWrap/Transverse.h \
	 ../../Pythia7/CLHEPWrap/Transverse.icc \
	 ../../Pythia7/Utilities/Direction.h \
	 ../../Pythia7/Utilities/Direction.xh \
	 ../../Pythia7/Utilities/Direction.icc \
	 ../../Pythia7/Utilities/Math.h \
	 ../../Pythia7/Utilities/Math.icc \
	 ../../Pythia7/CLHEPWrap/Lorentz5Vector.icc \
	 ../../Pythia7/CLHEPWrap/LorentzRotation.h \
	 ../../Pythia7/CLHEPWrap/LorentzRotation.fh \
	 ../../Pythia7/EventRecord/ColourBase.h \
	 ../../Pythia7/EventRecord/EventInfoBase.h \
	 ../../Pythia7/EventRecord/EventInfoBase.icc \
	 ../../Pythia7/EventRecord/ColourLine.h \
	 ../../Pythia7/EventRecord/ColourLine.icc \
	 ../../Pythia7/EventRecord/ColourBase.icc \
	 ../../Pythia7/EventRecord/SpinBase.h \
	 ../../Pythia7/EventRecord/SpinBase.icc \
	 ../../Pythia7/EventRecord/Particle.icc \
	 ../../Pythia7/PDT/ParticleData.h \
	 ../../Pythia7/PDT/PDT.h \
	 ../../Pythia7/PDT/PDT.icc \
	 ../../Pythia7/CLHEPWrap/PhysicalConstants.h \
	 ../../Pythia7/Utilities/Selector.h \
	 ../../Pythia7/Utilities/Selector.icc \
	 ../../Pythia7/Utilities/Selector.tcc \
	 ../../Pythia7/PDT/WidthGenerator.fh \
	 ../../Pythia7/PDT/MassGenerator.fh \
	 ../../Pythia7/PDT/DecayMode.fh \
	 ../../Pythia7/PDT/ParticleData.icc \
	 ../../Pythia7/EventRecord/Particle.tcc \
	 ../../Pythia7/Interface/Parameter.h \
	 ../../Pythia7/Interface/Parameter.xh \
	 ../../Pythia7/Interface/Parameter.fh \
	 ../../Pythia7/Interface/Parameter.icc \
	 ../../Pythia7/Interface/Parameter.tcc \
	 ../../Pythia7/Repository/Repository.h \
	 ../../Pythia7/Repository/BaseRepository.h \
	 ../../Pythia7/Repository/BaseRepository.xh \
	 ../../Pythia7/Repository/BaseRepository.icc \
	 ../../Pythia7/Repository/BaseRepository.tcc \
	 ../../Pythia7/PDT/MatcherBase.h \
	 ../../Pythia7/PDT/MatcherBase.icc \
	 ../../Pythia7/PDT/DecayMode.h \
	 ../../Pythia7/PDT/Decayer.h \
	 ../../Pythia7/PDT/Decayer.fh \
	 ../../Pythia7/Handlers/HandlerBase.h \
	 ../../Pythia7/Handlers/HandlerBase.icc \
	 ../../Pythia7/MatrixElement/Amplitude.h \
	 ../../Pythia7/Config/Complex.h \
	 ../../Pythia7/MatrixElement/Amplitude.icc \
	 ../../Pythia7/PDT/Decayer.icc \
	 ../../Pythia7/PDT/DecayMode.icc \
	 ../../Pythia7/Repository/Repository.icc \
	 ../../Pythia7/Repository/Repository.tcc \
	 /usera/stephens/projects/AMEGIC++-1.0/Include/Run_Parameter.H \
	 ../../Pythia7/Utilities/Rebinder.h \
	 ../../Pythia7/Utilities/Rebinder.icc \
	 ../../Pythia7/Interface/Reference.h \
	 ../../Pythia7/Interface/Reference.xh \
	 ../../Pythia7/Interface/Reference.fh \
	 ../../Pythia7/Interface/Reference.icc \
	 ../../Pythia7/Interface/Reference.tcc \
	 ../../Pythia7/PDT/EnumParticles.h

KtJetInterface.o: KtJetInterface.cc \
	 \
	 ../../Herwig++/Interfaces/KtJetInterface.h \
	 /var/nwork/pcck/gieseke/mc/ktjet/KtJet/KtLorentzVector.h \
	 /var/nwork/pcck/gieseke/mc/ktjet/KtJet/KtUtil.h \
	 /var/nwork/pcck/gieseke/mc/ktjet/KtJet/KtLorentzVector.icc \
	 ../../Pythia7/Config/Pythia7.h \
	 ../../Pythia7/Config/config.h \
	 ../../Pythia7/Pointer/Ptr.h \
	 ../../Pythia7/Pointer/RCPtr.fh \
	 ../../Pythia7/Pointer/PtrTraits.h \
	 ../../Pythia7/Pointer/RCPtr.h \
	 ../../Pythia7/Pointer/ReferenceCounted.h \
	 ../../Pythia7/Pointer/ReferenceCounted.icc \
	 ../../Pythia7/Pointer/RCPtr.icc \
	 ../../Pythia7/Utilities/Rebinder.fh \
	 ../../Pythia7/Utilities/Interval.fh \
	 ../../Pythia7/Utilities/ClassDescription.fh \
	 ../../Pythia7/Interface/InterfaceBase.fh \
	 ../../Pythia7/Persistency/PersistentOStream.fh \
	 ../../Pythia7/Persistency/PersistentIStream.fh \
	 ../../Pythia7/Config/Units.h \
	 ../../Pythia7/CLHEPWrap/Lorentz5Vector.fh \
	 ../../Pythia7/CLHEPWrap/LorentzVector.fh \
	 ../../Pythia7/CLHEPWrap/ThreeVector.fh \
	 ../../Pythia7/CLHEPWrap/Transverse.fh \
	 ../../Pythia7/CLHEPWrap/SystemOfUnits.h \
	 ../../Pythia7/CLHEPWrap/CLHEP.h \
	 ../../Pythia7/Config/Constants.h \
	 ../../Pythia7/Config/std.h \
	 ../../Pythia7/Config/Pointers.h \
	 ../../Pythia7/Config/Containers.h \
	 ../../Pythia7/Utilities/FixedSizeAllocator.h \
	 ../../Pythia7/Utilities/Debug.h \
	 ../../Pythia7/Utilities/Debug.icc \
	 ../../Pythia7/Utilities/FixedSizeAllocator.icc \
	 ../../Herwig++/Config/Herwig.h \
	 ../../Pythia7/Repository/Strategy.fh \
	 ../../Pythia7/EventRecord/Particle.h \
	 ../../Pythia7/EventRecord/EventConfig.h \
	 ../../Pythia7/CLHEPWrap/Lorentz5Vector.h \
	 ../../Pythia7/CLHEPWrap/LorentzVector.h \
	 ../../Pythia7/CLHEPWrap/Transverse.h \
	 ../../Pythia7/CLHEPWrap/Transverse.icc \
	 ../../Pythia7/Utilities/Direction.h \
	 ../../Pythia7/Utilities/Direction.xh \
	 ../../Pythia7/Misc/Exception.h \
	 ../../Pythia7/Misc/Exception.fh \
	 ../../Pythia7/Utilities/sstream.h \
	 ../../Pythia7/Misc/Exception.icc \
	 ../../Pythia7/Utilities/Direction.icc \
	 ../../Pythia7/Utilities/Math.h \
	 ../../Pythia7/Utilities/Math.icc \
	 ../../Pythia7/CLHEPWrap/Lorentz5Vector.icc \
	 ../../Pythia7/CLHEPWrap/LorentzRotation.h \
	 ../../Pythia7/CLHEPWrap/LorentzRotation.fh \
	 ../../Pythia7/Utilities/ClassDescription.h \
	 ../../Pythia7/Misc/Named.h \
	 ../../Pythia7/Misc/Named.icc \
	 ../../Pythia7/Utilities/ClassTraits.h \
	 ../../Pythia7/Utilities/DescriptionList.h \
	 ../../Pythia7/Utilities/DescriptionList.icc \
	 ../../Pythia7/Utilities/ClassDescription.icc \
	 ../../Pythia7/Utilities/ClassDescription.tcc \
	 ../../Pythia7/EventRecord/ColourBase.h \
	 ../../Pythia7/EventRecord/EventInfoBase.h \
	 ../../Pythia7/EventRecord/EventInfoBase.icc \
	 ../../Pythia7/EventRecord/ColourLine.h \
	 ../../Pythia7/EventRecord/ColourLine.icc \
	 ../../Pythia7/EventRecord/ColourBase.icc \
	 ../../Pythia7/EventRecord/SpinBase.h \
	 ../../Pythia7/EventRecord/SpinBase.icc \
	 ../../Pythia7/EventRecord/Particle.icc \
	 ../../Pythia7/PDT/ParticleData.h \
	 ../../Pythia7/PDT/PDT.h \
	 ../../Pythia7/PDT/PDT.icc \
	 ../../Pythia7/CLHEPWrap/PhysicalConstants.h \
	 ../../Pythia7/Interface/Interfaced.h \
	 ../../Pythia7/Interface/InterfacedBase.h \
	 ../../Pythia7/Utilities/HoldFlag.h \
	 ../../Pythia7/Utilities/HoldFlag.icc \
	 ../../Pythia7/Interface/InterfacedBase.xh \
	 ../../Pythia7/Interface/InterfaceBase.h \
	 ../../Pythia7/Interface/InterfaceBase.xh \
	 ../../Pythia7/Interface/InterfaceBase.icc \
	 ../../Pythia7/Interface/InterfacedBase.icc \
	 ../../Pythia7/Interface/Interfaced.icc \
	 ../../Pythia7/Utilities/Selector.h \
	 ../../Pythia7/Utilities/Selector.icc \
	 ../../Pythia7/Utilities/Selector.tcc \
	 ../../Pythia7/PDT/WidthGenerator.fh \
	 ../../Pythia7/PDT/MassGenerator.fh \
	 ../../Pythia7/PDT/DecayMode.fh \
	 ../../Pythia7/PDT/ParticleData.icc \
	 ../../Pythia7/EventRecord/Particle.tcc

# IF YOU PUT ANYTHING HERE IT WILL GO AWAY
