import ROOT 
import copy
from Samples import samp
import numpy as np

class MyAnalysis(object):
   
    def __init__(self, sample):

        """ The Init() function is called when an object MyAnalysis is initialised
        The tree corresponding to the specific sample is picked up 
        and histograms are booked.
        """

        self._tree = ROOT.TTree()        
        if(sample not in samp.keys() and sample != "data"):
            print (RuntimeError("Sample %s not valid. please, choose among these: %s" % (sample, str(samp.keys())) ))
            exit
        self.histograms = {}
        self.sample = sample
        self._file = ROOT.TFile("files/"+sample+".root")
        self._file.cd()
        tree = self._file.Get("events")
        self._tree = tree
        self.nEvents = self._tree.GetEntries()
        print ("Number of entries for " + self.sample + ": " + str(self.nEvents))
        
        ### Book histograms
        self.bookHistos()

    def getTree(self):
        return self._tree

    def getHistos(self):
        return self.histograms

    def bookHistos(self):
        h_nJet = ROOT.TH1F("NJet","#of jets", 3, 3.5, 6.5)
        h_nJet.SetXTitle("%# of jets")
        self.histograms["NJet"] = h_nJet 

        h_nJetFinal = ROOT.TH1F("NJetFinal","#of jets", 7, -0.5, 6.5)
        h_nJetFinal.SetXTitle("%# of jets")
        self.histograms["NJetFinal"] = h_nJetFinal 

        h_MuonIso = ROOT.TH1F("Muon_Iso","Muon Isolation", 25, 0., 3.)
        h_MuonIso.SetXTitle("Muon Isolation")
        self.histograms["Muon_Iso"] = h_MuonIso 

        h_NIsoMu = ROOT.TH1F("NIsoMu","Number of isolated muons", 5, 0.5, 5.5)
        h_NIsoMu.SetXTitle("Number of isolated muons")
        self.histograms["NIsoMu"] = h_NIsoMu 

        h_MuonPt = ROOT.TH1F("Muon_Pt","Muon P_T", 50, 0., 200.)
        h_MuonPt.SetXTitle("Muon P_T")
        self.histograms["Muon_Pt"] = h_MuonPt 

        h_MuonEta = ROOT.TH1F("Muon_eta","Muon eta", 50, 0., 5.)
        h_MuonEta.SetXTitle("Muon eta")
        self.histograms["Muon_eta"] = h_MuonEta 

        h_METpt = ROOT.TH1F("MET_Pt","MET P_T", 25, 0., 300.)
        h_METpt.SetXTitle("MET P_T")
        self.histograms["MET_Pt"] = h_METpt 

        h_JetPt = ROOT.TH1F("Jet_Pt","Jet P_T", 50, 0., 200.)
        h_JetPt.SetXTitle("Jet P_T")
        self.histograms["Jet_Pt"] = h_JetPt 

        h_JetBtag = ROOT.TH1F("Jet_Btag","Jet B tag", 10, 1., 6.)
        h_JetBtag.SetXTitle("Jet B tag")
        self.histograms["Jet_btag"] = h_JetBtag 

        h_NBtag = ROOT.TH1F("NBtag","Jet B tag", 4, 0.5, 4.5)
        h_NBtag.SetXTitle("Number of B tagged jets")
        self.histograms["NBtag"] = h_NBtag 


    def saveHistos(self):
        outfilename = self.sample + "_histos.root"
        outfile = ROOT.TFile(outfilename, "RECREATE")
        outfile.cd()
        for h in self.histograms.values():
            h.Write()
        outfile.Close()

    ### processEvent function implements the actions to perform on each event
    ### This is the place where to implement the analysis strategy: study of most sensitive variables
    ### and signal-like event selection

    def processEvent(self, entry):
        tree = self.getTree()
        tree.GetEntry(entry)
        w = tree.EventWeight

        event_passed = True

        # First, we check if the event passes the cut. Then we fill the histograms.
        # Along the way, I save some variables I need

        # -------- MUON CUTS ----------
        muonPtCut = 25. 
        muonRelIsoCut = 0.03
        nIsoMu = 0
        for m in range(tree.NMuon):
            muon = ROOT.TLorentzVector(tree.Muon_Px[m],tree.Muon_Py[m],tree.Muon_Pz[m],tree.Muon_E[m])
            if(muon.Pt() > muonPtCut and (tree.Muon_Iso[m]/muon.Pt() < muonRelIsoCut)):
                nIsoMu += 1
        if(nIsoMu != 1): event_passed = False


        # -------- NEUTRINO CUTS ----------
        MET = np.sqrt(tree.MET_px**2 + tree.MET_py**2)
        if(MET < 20): event_passed = False

        # -------- JET CUTS ----------
        njets = tree.NJet
        if(njets < 4): event_passed = False
        onejetabove25Pt = False
        onejetBtagged = False
        for j in range(tree.NJet):
            jet = ROOT.TLorentzVector(tree.Jet_Px[j],tree.Jet_Py[j],tree.Jet_Pz[j],tree.Jet_E[j])
            if(jet.Pt() > 25): onejetabove25Pt = True
            if(tree.Jet_btag[j] > 2.0): onejetBtagged = True
        if(not onejetabove25Pt): event_passed = False
        if(not onejetBtagged): event_passed = False

        # -------- TRIGGER ----------
        if(not tree.triggerIsoMu24 == 1): event_passed = False

        
        # fill histograms
        if(not event_passed):
            return
        
        for m in range(tree.NMuon):
            muon = ROOT.TLorentzVector(tree.Muon_Px[m],tree.Muon_Py[m],tree.Muon_Pz[m],tree.Muon_E[m])
            self.histograms["Muon_Iso"].Fill(tree.Muon_Iso[m], w)                
            self.histograms["Muon_eta"].Fill(np.abs(muon.PseudoRapidity()), w)                
            if(muon.Pt() > muonPtCut and (tree.Muon_Iso[m]/muon.Pt()) < muonRelIsoCut):
                self.histograms["Muon_Pt"].Fill(muon.Pt(), w)
        self.histograms["NIsoMu"].Fill(nIsoMu, w)
        self.histograms["MET_Pt"].Fill(MET,w)
        self.histograms["NJet"].Fill(njets,w)


    ### processEvents run the function processEvent on each event stored in the tree    
    def processEvents(self):
        nevts = self.nEvents
        for i in range(nevts):
            self.processEvent(i)
        
        self.saveHistos()




