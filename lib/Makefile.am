lib_LTLIBRARIES = libHerwig.la
libHerwig_la_SOURCES =
libHerwig_la_LIBADD = \
$(top_builddir)/UnderlyingEvent/libUA5.la \
$(top_builddir)/Models/StandardModel/libHwStandardModel.la \
$(top_builddir)/Models/RSModel/libHwRSModel.la \
$(top_builddir)/Decay/FormFactors/libHwFormFactor.la \
$(top_builddir)/Decay/Tau/libHwTauDecay.la \
$(top_builddir)/Decay/Baryon/libHwBaryonDecay.la \
$(top_builddir)/Decay/VectorMeson/libHwVMDecay.la \
$(top_builddir)/Decay/Perturbative/libHwPerturbativeDecay.la \
$(top_builddir)/Decay/WeakCurrents/libHwWeakCurrent.la \
$(top_builddir)/Decay/ScalarMeson/libHwSMDecay.la \
$(top_builddir)/Decay/TensorMeson/libHwTMDecay.la \
$(top_builddir)/Decay/Partonic/libHwPartonicDecay.la \
$(top_builddir)/Decay/libHwDecay.la \
$(top_builddir)/PDT/libHwPDT.la \
$(top_builddir)/Shower/libHwShower.la \
$(top_builddir)/Hadronization/libHwHadronization.la \
$(top_builddir)/Helicity/Correlations/libHwCorrelations.la \
$(top_builddir)/Helicity/WaveFunction/libHwWaveFunction.la \
$(top_builddir)/Helicity/Vertex/StandardModel/libHwSMVertex.la \
$(top_builddir)/Helicity/Vertex/Vector/libHwVVertex.la \
$(top_builddir)/Helicity/Vertex/RSModel/libHwRSVertex.la \
$(top_builddir)/Helicity/Vertex/Scalar/libHwSVertex.la \
$(top_builddir)/Helicity/Vertex/Tensor/libHwTVertex.la \
$(top_builddir)/Helicity/Vertex/libHwVertex.la \
$(top_builddir)/Utilities/libHwUtils.la \
$(top_builddir)/MatrixElement/libHwME.la

all-local: done-all-links

done-all-links:
	$(LN_S) ../*/.libs/*.so .
	echo "stamp" > done-all-links

clean-local:
	rm -f *.so done-all-links
