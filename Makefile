CXXFLAGS = -Wall
ROOTFLAGS = `$(ROOTSYS)/bin/root-config --cflags --libs`
FASTJETFLAGS = `$(FASTJETSYS)/bin/fastjet-config --cxxflags --libs --plugins`

makeHistos.out : src/analyzeJets.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/makeHistos.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

writeJetMomenta.out : src/writeJetMomenta.cc src/helpers.h
	$(CXX) $(CXXFLAGS) -o bin/writeJetMomenta.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)

histos : 
	./bin/makeHistos.out events.root histos.root 

data : 
	./bin/writeJetMomenta.out events.root matchedJets.txt
	
.PHONY: data histos 
