#!/usr/bin/env python

import glob
import ROOT as root
from src.helpers import *

root.gROOT.SetBatch(True)

f = root.TFile("wHmumu.root")

canvas = root.TCanvas()

hist = f.Get("genJMass")
hist2 = f.Get("jetMass")
hist3 = f.Get("diJetOMass")
hist4 = f.Get("diGenOMass")

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
hist.SetMarkerStyle(0)
hist.SetLineStyle(7)
hist2.SetLineColor(1)
hist2.SetMarkerStyle(0)
hist3.SetLineColor(3)
hist3.SetMarkerStyle(0)
hist4.SetLineColor(3)
hist4.SetMarkerStyle(0)
hist4.SetLineStyle(7)

leg = root.TLegend(.2,.7,.4,.9,"Jet Type")
leg.AddEntry(hist,"Gen Jets","l")
leg.AddEntry(hist2,"Reco Jets","l")
leg.AddEntry(hist3,"Dijet Pt pair","l")
leg.AddEntry(hist4,"Dijet Gen pair","l")

hist.Draw("HIST same")
hist2.Draw("HIST same")
hist3.Draw("HIST same")
hist4.Draw("HIST same")
leg.Draw()

canvas.SaveAs("Hist_jjPairMComp.png")
