#!/usr/bin/env python

import glob
import ROOT as root
from src.helpers import *

root.gROOT.SetBatch(True)

f = root.TFile("wHmumu.root")

canvas = root.TCanvas()

hist = f.Get("DiffDimuPt")

xMin = -200
xMax = 200
yMin = 0
yMax = 200
xTitle = "Pt Difference [GeV/c^{2}]"

axisHist = root.TH2F("axisHist","",1,xMin,xMax,1,yMin,yMax)
axisHist.GetXaxis().SetTitle(xTitle)
axisHist.GetYaxis().SetTitle("Dimuon Pt [GeV/c^{2}]")
axisHist.Draw()


hist.Draw("COL same")

canvas.SaveAs("Hist_DiffDimu.png")
