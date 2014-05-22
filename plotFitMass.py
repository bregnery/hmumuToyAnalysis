#!/usr/bin/env python

import glob
import ROOT as root
from src.helpers import *

root.gROOT.SetBatch(True)

f = root.TFile("wHmumu.root")

canvas = root.TCanvas()

hist = f.Get("jetMass")
hist2 = f.Get("jjfitM")

#hist.Fit("gaus","","",60,100)
#hist.GetFunction("gaus").SetLineColor(1)
#hist2.Fit("gaus","","",60,100)
#hist2.GetFunction("gaus").SetLineColor(2)


xMin = 0
xMax = 200
yMin = 0
yMax = 6
xTitle = "2 Jet Mass [GeV/c^{2}]"

axisHist = root.TH2F("axisHist","",1,xMin,xMax,1,yMin,yMax)
axisHist.GetXaxis().SetTitle(xTitle)
axisHist.GetYaxis().SetTitle("Events/Bin")
axisHist.Draw()

hist.SetLineColor(1)
hist.SetFillStyle(1)
hist.SetFillColor(1)
hist.SetMarkerStyle(0)
hist2.SetFillStyle(1)
hist2.SetLineColor(2)
hist2.SetFillColor(2)
hist2.SetMarkerStyle(0)

hist.Draw("hist same")
hist2.Draw("hist same")

canvas.SaveAs("Hist_jjfitMComp.png")
