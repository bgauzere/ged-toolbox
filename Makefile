TARGETS_OCTAVE = computeEditDistance.mex costpaths/bagOfSimpleLabeledPaths.mex costwalks/labeledKron.mex hungarianLSAP.mex quadraticTerm.mex utils/DatasetToAdjacency.mex hungarianLSAPE.mex quadraticTermLSAPE.mex LSAPEtoLSAPMapping.mexa64
TARGETS_MATLAB = computeEditDistance.mexa64 costpaths/bagOfSimpleLabeledPaths.mexa64 costwalks/labeledKron.mexa64 hungarianLSAP.mexa64 quadraticTerm.mexa64 utils/DatasetToAdjacency.mexa64 hungarianLSAPE.mexa64 quadraticTermLSAPE.mexa64 computeDelta.mexa64 computeDeltaSub.mexa64 LSAPEtoLSAPMapping.mexa64

octave: MEX_CC = mkoctfile --mex
octave: CXXFLAGS = -I$(IDIR) -Wall -I./
octave: $(TARGETS_OCTAVE)

matlab: EXT=mexa64
matlab: MEX_CC = ~/bin/matlab2013a/2013a/bin/mex
matlab: CXXFLAGS = -I$(IDIR) -I./ -O
matlab: $(TARGETS_MATLAB)

IDIR = ./graph-lib/include
LIBDIR = ./graph-lib/obj
OBJS = $(LIBDIR)/graph.o  $(LIBDIR)/GraphEditDistance.o 




%.mex: %.cpp $(OBJS)
	$(MEX_CC) $(CXXFLAGS) $< $(OBJS)

%.mexa64: %.cpp $(OBJS)
	$(MEX_CC) $(CXXFLAGS) $< $(OBJS)

clean:
	rm -f *mex *mexa64 *.o
