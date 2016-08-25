IDIR = ./include
CXXFLAGS = -I$(IDIR) -Wall  -fpermissive -fPIC -std=c++11
BINDIR = ./bin



ODIR = ./obj
SRCDIR = ./src
_DEPS = graph.h GraphEditDistance.h 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_SRCDEPS = graph.cpp GraphEditDistance.cpp
DEPS_SRC += $(patsubst %,$(SRCDIR)/%,$(_DEPS_SRC))

_OBJ = GraphEditDistance.o graph.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: $(BINDIR)/test_GraphEditDistance $(BINDIR)/contestGraphEditDistance

debug: CXXFLAGS += -DDEBUG -g
debug: all

optim: CXXFLAGS += -O3
optim: all


$(BINDIR)/contestGraphEditDistance: $(OBJ) $(SRCDIR)/contestGraphEditDistance.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

$(BINDIR)/test_GraphEditDistance: $(OBJ) $(SRCDIR)/test_GraphEditDistance.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

$(BINDIR)/%: $(OBJ) $(SRCDIR)/%.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

$(ODIR)/%.o: $(SRCDIR)/%.cpp $(DEPS) $(DEPS_SRC)
	$(CXX) -c -o $@ $< $(CXXFLAGS)


.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~  $(BINDIR)/*
