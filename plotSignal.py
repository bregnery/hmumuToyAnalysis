#!/usr/bin/env python

import glob
import ROOT as root
from src.helpers import *

root.gROOT.SetBatch(True)

f = root.TFile("wHmumu.root")

canvas = root.TCanvas()

hist = f.Get("jetMass")
hist2 = f.Get("JetMass1MET")
hist3 = f.Get("JetMass2MET")
hist4 = f.Get("JetMass85MET")

hist.Fit("gaus","","",60,100)
hist.GetFunction("gaus").SetLineColor(1)
hist2.Fit("gaus","","",60,100)
hist2.GetFunction("gaus").SetLineColor(2)
hist3.Fit("gaus","","",60,100)
hist3.GetFunction("gaus").SetLineColor(4)
hist4.Fit("gaus","","",60,100)
hist4.GetFunction("gaus").SetLineColor(6)

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
hist.SetMarkerColor(1)
hist2.SetLineColor(2)
hist2.SetMarkerColor(2)
hist3.SetLineColor(4)
hist3.SetMarkerColor(4)
hist4.SetLineColor(6)
hist4.SetMarkerColor(6)

hist.Draw("same")
hist2.Draw("same")
hist3.Draw("same")
hist4.Draw("same")
#hist.Draw("HIST L same")
#hist2.Draw("HIST L same")
#hist3.Draw("HIST L same")
#hist4.Draw("HIST L same")

canvas.SaveAs("Hist_JetMassImprov.png")
