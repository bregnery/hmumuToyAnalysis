#!/usr/bin/env python

import glob
import ROOT as root
from src.helpers import *

root.gROOT.SetBatch(True)

f = root.TFile("wHmumu.root")

canvas = root.TCanvas()

hist = f.Get("j1fitPt")

#hist.Fit("gaus","","",-100,100)

xMin = 0
xMax = 200
yMin = 0
yMax = 20
xTitle = "Pt [GeV/c]"

axisHist = root.TH2F("axisHist","",1,xMin,xMax,1,yMin,yMax)
axisHist.GetXaxis().SetTitle(xTitle)
axisHist.GetYaxis().SetTitle("Events/Bin")
axisHist.Draw()


hist.Draw("histo same")

canvas.SaveAs("Hist_j1fitPt.png")