CXXFLAGS = -Wall
ROOTFLAGS = `$(ROOTSYS)/bin/root-config --cflags --libs`
FASTJETFLAGS = `$(FASTJETSYS)/bin/fastjet-config --cxxflags --libs --plugins`
makeJets.out : src/partonJets.cc
	$(CXX) $(CXXFLAGS) -o bin/makeJets.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

clean :
	rm data/histos.root
