
import sys
from MyAnalysis import MyAnalysis
from ROOT import TTree, TFile, double
from Plotter import getBkgHisto, getHisto, plotVar, plotVarNorm, plotShapes, getSigHisto
from pymule import *
from functools import *
import numpy as np
import ctypes
from matplotlib.pyplot import *

def add(vals):
  # adds values with errors
  res = [0,0]
  for val in vals:
    res[0] += val[0]
    res[1] += val[1]**2
  res[1] = np.sqrt(res[1])
  return np.array(res)

def mult(vals):
  #multiplay values with errosr
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

def percent(val):
  return str(round(abs(val[1]/val[0])*100000.)/1000.) + '‰'


def process_events():
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

# process_events()  # TOGGLE THAT ONE

# samples = ["qcd", "zz", "wz", "ww",  "single_top", "dy","wjets","ttbar"]
samples = ["qcd", "dy","wjets"]

# create pdfs for these
vars = ["Muonp_eta","Muonm_eta","MET_p","MET_m","mWp","mWm"]

for v in vars:
  print("Variable: ", v)  # plotShapes (variable, samples,logScale )
  plotShapes(v, samples,  False)
  plotVar(v, samples,  True, False) # plotVar(variable, samples,isData, logScale )


def getIntegral(histo):
  minBin = histo.GetXaxis().GetFirst()
  maxBin = histo.GetXaxis().GetLast()
  nDataErr = ctypes.c_double()
  nData = histo.IntegralAndError(minBin, maxBin, nDataErr)
  return nData, nDataErr.value

def getIntegralInRange(histo,range):
  minBin = range[0]
  maxBin = range[1]
  nDataErr = ctypes.c_double()
  nData = histo.IntegralAndError(minBin, maxBin, nDataErr)
  return nData, nDataErr.value



print('\n\n----------- signal and background ----------')

# for acceptance (excluding trigger)
mc_wjets_plus=[
  [11092.86897278, 158.73256466],
  [10581.54048157, 154.82274318],
  [11315.79328918, 160.18202564],
  [11311.9407196 , 160.28179247],
  [10581.63168335, 155.09028513 ],
  [10955.44876099, 157.71305959],
  [10504.94567871, 154.6395998],
  [10326.83833313, 153.13294028],
  [10678.53982544, 155.6018195],
  [10062.38841248, 150.93913333]
]
mc_wjets_minus=[
  [8127.68174744, 135.7751357],
  [8082.2029953 , 135.35601153],
  [8314.84244537, 137.25817777],
  [8483.17984009, 138.88843318],
  [7005.9394989 , 126.38788524 ],
  [7175.13549042, 127.74574654],
  [6934.59186554, 125.48050082],
  [6682.91755676, 123.34063736],
  [6339.49119568, 119.96630658],
  [5613.02454376, 112.87557078]
]
mc_wjets_plus_withcuts = [
  [8852.38839722, 142.09539595],
  [8438.86091614, 138.56661597],
  [9065.72822571, 143.73701607],
  [9009.22663879, 143.29596476],
  [8466.79747009, 139.05611548],
  [8786.59803009, 141.55911622],
  [8282.81063843, 137.55675226],
  [8337.2121582 , 137.9173904],
  [8388.09836578, 138.23706529],
  [8129.73973083, 136.05288352]
]
mc_wjets_minus_withcuts = [
  [6613.08628082, 122.81966499],
  [6578.48357391, 122.4850103 ],
  [6803.95197296, 124.53285731],
  [6885.35274506, 125.20749821],
  [5746.2306366 , 114.68535027],
  [5880.58129883, 115.88098179],
  [5768.89819336, 114.68506754],
  [5520.39295959, 112.21725121],
  [5148.12277603, 108.40759907],
  [4521.37570572, 101.43760033]
]



print('compute acceptance and trigger efficiency')
a_p = []
a_m = []
for i in range(0,10):
  acceptance_plus = mult([mc_wjets_plus_withcuts[i],inv(mc_wjets_plus[i])])
  acceptance_minus = mult([mc_wjets_minus_withcuts[i],inv(mc_wjets_minus[i])])
  a_p.append(acceptance_plus)
  a_m.append(acceptance_minus)
print(a_p)
print(a_m)

# for trigger efficiency

# for etamax in np.linspace(0.2,2.0,10):
ratios = []
ratiosp = []
ratiosm = []


for index,maxbin in enumerate(range(50,550,50)):
  range = [maxbin-49,maxbin] # corresponds to eta in [0,0.2]
  print('\n-------------------------------- ')
  print('doing range', np.round(np.array(range)/500*2,1))

  mcsignal_p = np.array(getIntegralInRange(getSigHisto('Muonp_eta'),range))
  background_p = np.array(getIntegralInRange(getBkgHisto('Muonp_eta', samples),range))
  mcsignal_m = np.array(getIntegralInRange(getSigHisto('Muonm_eta'),range))
  background_m = np.array(getIntegralInRange(getBkgHisto('Muonm_eta', samples),range))
  print('mcsignal +',mcsignal_p)
  print('background +',background_p)
  print('mcsignal -',mcsignal_m)
  print('background -', background_m)

  observations_p = np.array(getIntegralInRange(getHisto('Muonp_eta','data'),range))
  observations_m = np.array(getIntegralInRange(getHisto('Muonm_eta','data'),range))
  print('observations +',observations_p)
  print('observations -',observations_m)
  observations_p[1] *= 1.2
  observations_m[1] *= 1.2
  
  signal_p = add([observations_p,-1.*background_p])
  signal_m = add([observations_m,-1.*background_m])
  print('signal +', signal_p)
  print('signal -', signal_m)
  
  ratiop = signal_p[0]/mcsignal_p[0]
  ratiom = signal_m[0]/mcsignal_m[0]
  # print('ratio +', ratiop,'ratio -',ratiom)
  ratiosp.append(ratiop)
  ratiosm.append(ratiom)

  corr_p = mult([signal_p,inv(a_p[index])])
  corr_m = mult([signal_m,inv(a_m[index])])
  print('-----------------')
  print(corr_p)
  print(corr_m)
  print('-----------------')

  sum = add([corr_p,corr_m])
  diff = add([corr_p,-1*corr_m])
  ratio = mult([diff,inv(sum)])
  ratios.append(ratio)
  print('ratio',ratio * 100, '%')

ratios = np.array(ratios)
# plot results
subplots(2,1,gridspec_kw={'height_ratios': [2, 1]})
subplot(2,1,1)
eta = np.array([0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8])+0.1
plt.scatter(eta,ratios[:,0],marker='_',s=np.ones(10)*1400,linewidth=5,color="black")
errorbar(eta,ratios[:,0],ratios[:,1]*1.0,fmt='none',capsize=5,color="black")
grid(color="#ddd")
xticks(eta-0.1)
# xlim([0,2])
# ylim([0.05,0.28])
ylabel('Muon charge asymmetry')
gca().set_axisbelow(True)

subplot(2,1,2)
plt.scatter(eta,ratiosp)
plt.scatter(eta,ratiosm)
legend(['data/MC for $W_+$','data/MC for $W_-$'])
grid(color="#ddd")
xlabel('Muon Pseudorapidity $\eta$')
ylabel('$ratio$')
xticks(eta-0.1)
ylim([0.8,1.2])

show()

