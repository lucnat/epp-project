
from MyAnalysis import MyAnalysis
from ROOT import TTree, TFile, double
from Plotter import getBkgHisto, getHisto, plotVar, plotVarNorm, plotShapes, getSigHisto
from pymule import *
from functools import *
import numpy as np
import ctypes

def add(vals):
  # adds values with errors
  res = [0,0]
  for val in vals:
    res[0] += val[0]
    res[1] += val[1]**2
  res[1] = np.sqrt(res[1])
  return np.array(res)

def mult(vals):
  # multiplay values with errors
  res = [1,0]
  for val in vals:
    res[0] *= val[0]
    res[1] += (val[1]/val[0])**2 # add relative errors
  res[1] = np.sqrt(res[1])*res[0]
  return np.array(res)

def inv(val):
  # compute inverse of value with error
  res = 1./val[0]
  res_rel_err = val[1]/val[0]
  res_err = res * res_rel_err
  return np.array([res,res_err])

def compare(a, b):
  # compares two values
  print('relative agreement = ', round(abs((a[0]-b[0])/a[0])*1000000)/1000, '‰')
  print('chi2 = ', chisq([a,b]))

def permille(val):
  return str(round(abs(val[1]/val[0])*1000000.)/1000.) + '‰'

process = False
if process:
  
  ### Instantiation of an object of kind MyAnalysis for each single sample
  TT = MyAnalysis("ttbar")
  TT.processEvents()

  DY = MyAnalysis("dy")
  DY.processEvents()

  QCD = MyAnalysis("qcd")
  QCD.processEvents()

  SingleTop = MyAnalysis("single_top")
  SingleTop.processEvents()

  WJets = MyAnalysis("wjets")
  WJets.processEvents()

  WW = MyAnalysis("ww")
  WW.processEvents()

  ZZ = MyAnalysis("zz")
  ZZ.processEvents()

  WZ = MyAnalysis("wz")
  WZ.processEvents()

  Data = MyAnalysis("data")
  Data.processEvents()


samples = ["qcd", "zz", "wz", "ww",  "single_top", "dy","wjets","ttbar"]

# create pdfs for these
vars = ["NIsoMu", "Muon_Pt", "MET_Pt", "NJet", "Muon_eta", "NJetFinal"]

for v in vars:
    print("Variable: ", v)
    ### plotShapes (variable, samples,logScale )
    plotShapes(v, samples,  False)
    ### plotVar(variable, samples,isData, logScale )
    plotVar(v, samples,  True, False)

def getIntegral(histo):
  minBin = histo.GetXaxis().GetFirst()
  maxBin = histo.GetXaxis().GetLast()
  nDataErr = ctypes.c_double()
  nData = histo.IntegralAndError(minBin, maxBin, nDataErr)
  return nData, nDataErr.value

print('\n\n----------- signal and background ----------')
for var in vars:
  print(var,'signal')
  (sig, sig_err) = getIntegral(getSigHisto(var))
  print(sig, sig_err)
  print(var, 'background')
  (bg, bg_err) = getIntegral(getBkgHisto(var, samples))
  print(bg, bg_err)

print('\n\n----------- perform subtraction ----------')
observations = np.array(getIntegral(getHisto('MET_Pt','data')))
print('observations',observations)
background = np.array(getIntegral(getBkgHisto('MET_Pt', samples)))
print('background',background)
diff = add([observations,-1.*background])
print('diff',diff)

print('\n\n----------- evaluate acceptance and trigger efficiency ----------')
ttbar = np.array(getIntegral(getHisto('MET_Pt','ttbar')))
print('ttbar',ttbar)
acceptance = mult([ttbar,inv([7928.61,np.sqrt(7928.61)])])
print('acceptance = ', acceptance)

print('\n\n----------- cross section ----------')
L = 50
trigger_eff = [0.88, 0.03]
sigma = mult([diff,inv(mult([acceptance,trigger_eff]))]) / L 
print(sigma)
theory = [173.60,11.7]
print('theory', theory)
print('chi2',chisq([sigma,theory]))

