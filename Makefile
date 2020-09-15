CXXFLAGS = -Wall
ROOTFLAGS = `$(ROOTSYS)/bin/root-config --cflags --libs`
FASTJETFLAGS = `$(FASTJETSYS)/bin/fastjet-config --cxxflags --libs --plugins`
makeJets.out : src/partonJets.cc
	$(CXX) $(CXXFLAGS) -o bin/makeJets.out $^ $(ROOTFLAGS) $(FASTJETFLAGS)
	
data : 
	./bin/makeJets.out JetNtuple_NoCut.root histosNoCut.root jetInfoNoCut.txt

dataCut : bin/makeJets.out
	./bin/makeJets.out JetNtuple_30GeVCut.root histosCut.root jetInfoCut.txt

dataGenJets : 
	./bin/makeJets.out JetNtuple_Cut_wgenJets.root histosNoCut.root jetInfoNoCut.txt

clean :
	rm data/histos.root

.PHONY: data dataCut clean
