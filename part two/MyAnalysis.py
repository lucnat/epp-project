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

        h_MuonpEta = ROOT.TH1F("Muonp_eta","Muon+ eta", 500, 0., 2.0)
        h_MuonpEta.SetXTitle("Muon+ eta")
        self.histograms["Muonp_eta"] = h_MuonpEta 

        h_MuonmEta = ROOT.TH1F("Muonm_eta","Muon- eta", 500, 0., 2.0)
        h_MuonmEta.SetXTitle("Muon- eta")
        self.histograms["Muonm_eta"] = h_MuonmEta 

        h_MET_p = ROOT.TH1F("MET_p","MET +", 50, -10., 170.)
        h_MET_p.SetXTitle("MET +")
        self.histograms["MET_p"] = h_MET_p 

        h_MET_m = ROOT.TH1F("MET_m","MET -", 50, -10., 170.)
        h_MET_m.SetXTitle("MET -")
        self.histograms["MET_m"] = h_MET_m 

        mWp = ROOT.TH1F("mWp","Transverse W+ mass", 100, 0., 150.)
        mWp.SetXTitle("W mass (transverse)")
        self.histograms["mWp"] = mWp 

        mWm = ROOT.TH1F("mWm","Transverse W- mass", 100, 0., 150.)
        mWm.SetXTitle("W mass (transverse)")
        self.histograms["mWm"] = mWm 

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
        muonRelIsoCut = 0.1
        etacut = 2.0
        nMu = 0
        for m in range(tree.NMuon):
            muon = ROOT.TLorentzVector(tree.Muon_Px[m],tree.Muon_Py[m],tree.Muon_Pz[m],tree.Muon_E[m])
            if(muon.Pt() > muonPtCut and np.abs(muon.PseudoRapidity()) < etacut and tree.Muon_Iso[m]/muon.Pt() < muonRelIsoCut):
                nMu += 1
        if(not nMu > 0): event_passed = False


        # -------- NEUTRINO CUTS ----------
        MET = np.sqrt(tree.MET_px**2 + tree.MET_py**2)
        if(MET < 20): event_passed = False

        # -------- TRIGGER ----------
        if(not tree.triggerIsoMu24 > 0): event_passed = False

        # onejetBtagged = False
        # for j in range(tree.NJet):
        #     if(tree.Jet_btag[j] > 2.0): onejetBtagged = True
        # if(not onejetBtagged): event_passed = False


        # fill histograms
        if(not event_passed):
            return
        
        for m in range(tree.NMuon):
            muon = ROOT.TLorentzVector(tree.Muon_Px[m],tree.Muon_Py[m],tree.Muon_Pz[m],tree.Muon_E[m])
            if(muon.Pt() > muonPtCut and np.abs(muon.PseudoRapidity()) < etacut and tree.Muon_Iso[m]/muon.Pt() < muonRelIsoCut):
                muon_T = ROOT.TLorentzVector(tree.Muon_Px[m],tree.Muon_Py[m],0.,np.sqrt(tree.Muon_E[m]**2-tree.Muon_Pz[m]**2))
                MET_vec = ROOT.TLorentzVector(tree.MET_px,tree.MET_py,0.,MET)
                mW = (muon_T + MET_vec).M()
                # if(not 86.0 < mW < 86.5):
                if(tree.Muon_Charge[m] > 0):
                    self.histograms["Muonp_eta"].Fill(np.abs(muon.PseudoRapidity()), w)
                    self.histograms["MET_p"].Fill(MET, w)
                    self.histograms["mWp"].Fill(mW,w)
                if(tree.Muon_Charge[m] < 0):
                    self.histograms["Muonm_eta"].Fill(np.abs(muon.PseudoRapidity()), w)
                    self.histograms["MET_m"].Fill(MET, w)
                    self.histograms["mWm"].Fill(mW,w)
 
    ### processEvents run the function processEvent on each event stored in the tree    
    def processEvents(self):
        nevts = self.nEvents
        for i in range(nevts):
            self.processEvent(i)
        
        self.saveHistos()




