#!/usr/bin/env python

import glob
import ROOT as root
from src.helpers import *

root.gROOT.SetBatch(True)


#####################################
## Configuration (all for 8TeV)

# Scale factor for higgs cross-sections
SIGSF=100.

luminosity = 19.712 # fb^-1 Pixel-Lumi for 22Jan2013 ReRECO 

#originalNEventsDict = {
  # Signal Needs to be checked!!
 # "ggHmumu":100000,
 # "vbfHmumu":100000,
 # "wHmumu":100000,
 # "zHmumu":100000,

 # "DY":30086987,
 # "TT":6921652,

 # "WW":5218045,
 # "WZ":10000283,
 # "ZZ":10320000,
#}

# in pb
xsecDict = {
  "ggHmumu":0.0042944*SIGSF,
  "vbfHmumu":0.00034716*SIGSF,
  "wHmumu":0.000153252*SIGSF,
  "zHmumu":8.6746e-05*SIGSF,

  "DY":3503.71,
  "TT":225.197,

  "WW":54.838,
  "WZ":33.21,
  "ZZ":17.654,
}

colorsDict = {
  "ggHmumu":root.kRed,
  "vbfHmumu":root.kBlue,
  "wHmumu":root.kGreen+1,
  "zHmumu":root.kOrange+7,

  "DY":root.kOrange,
  "TT":root.kGreen-1,

  "WW":root.kCyan,
  "WZ":root.kMagenta,
  "ZZ":root.kPink,


}

#################################################3

# Load files
inputfileNames = glob.glob("*.root")
files = []
for fn in inputfileNames:
  files.append(root.TFile(fn))

# Find histogram names in files
histNames = []
for key in files[0].GetListOfKeys():
  if key.GetClassName() == "TH1F":
    histNames.append(key.GetName())

# Find sample names from file names
inputfileNamesNoSuffix = []
for fn in inputfileNames:
  inputfileNamesNoSuffix.append(fn[:-5])

# split files into MC and data
mcFiles = []
mcSampleNames = []
dataFiles = []
for sampleName,rootFile in zip(inputfileNamesNoSuffix,files):
  if xsecDict.has_key(sampleName):
     mcFiles.append(rootFile)
     mcSampleNames.append(sampleName)
  else:
    dataFiles.append(rootFile)

canvas = root.TCanvas()
canvas.SetLogy(True)
# loop over histograms
for histName in histNames:
  mcHistList = []
  for f,sampleName in zip(mcFiles,mcSampleNames):
    hist = f.Get(histName)
   # hist.Scale(luminosity*1000.*xsecDict[sampleName]/originalNEventsDict[sampleName])
    hist.SetFillColor(colorsDict[sampleName])
    hist.SetFillStyle(1)
    hist.SetLineColor(colorsDict[sampleName])
    hist.SetMarkerStyle(0)
    hist.SetTitle(sampleName)
    mcHistList.append(hist)
  #mcSampleNames = sorted(mcSampleNames,key= lambda name: sorted(mcHistList).index(mcHistList[mcSampleNames.index(name)]))
  mcHistList.sort(key=lambda h: h.Integral())
  mcHistStack = root.THStack("stacked"+histName,"")
  mcHistSum = None
  if len(mcHistList)>0:
    mcHistSum = mcHistList[0].Clone("sumHist_"+histName)
    mcHistSum.Reset()
  for mcH in mcHistList:
    mcHistSum.Add(mcH)
    mcHistStack.Add(mcH)
  dataHist = None
  for f in dataFiles:
    hist = f.Get(histName)
    if dataHist == None:
      dataHist = hist
      dataHist.SetLineColor(root.kBlack)
    else:
      dataHist.Add(hist)

  maximum = mcHistSum.GetMaximum()
  if dataHist != None:
    maximum = max(maximum,dataHist.GetMaximum())
  maximum *= 1.2
  xMin = mcHistSum.GetXaxis().GetBinLowEdge(1)
  xMax = mcHistSum.GetXaxis().GetBinLowEdge(mcHistSum.GetNbinsX()+1)

  # Create 2D hist to have fine control over the axis limits
  axisHist = root.TH2F("axisHist"+histName,"",1,xMin,xMax,1,0.1,maximum)
  axisHist.GetXaxis().SetTitle(mcHistSum.GetXaxis().GetTitle())
  axisHist.GetYaxis().SetTitle("Events/Bin")
  
  # Create a legend
  leg = root.TLegend(.2,.7,.4,.9,"Decay type")
  for hist in mcHistList:  
    leg.AddEntry(hist,hist.GetTitle(),"f")
  
  axisHist.Draw()
  mcHistStack.Draw("hist same")
#  leg.Draw("same")
  if dataHist != None:
    dataHist.Draw("same")
  canvas.RedrawAxis()
  canvas.SaveAs("Hist_"+histName+".png")
    


