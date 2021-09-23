#! /usr/bin/env python
import os 
import shutil
import glob
import math
import array
import sys
import time
import argparse                                                                                                                                                          

import ROOT                                                                                                                                                              
import CMS_lumi, tdrstyle                                                                                                                                                                          
#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetFitFormat('3.2g')
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.04,'X')
ROOT.gStyle.SetLabelSize(0.04,'Y')
ROOT.gStyle.SetTitleSize(0.04,'X')
ROOT.gStyle.SetTitleSize(0.04,'Y')
ROOT.gStyle.SetTitleOffset(1.1,'X')
ROOT.gStyle.SetTitleOffset(1.4,'Y')
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


Vov = 1.50

temperatures = [0, -25, -40]

tRes = {}
err_dcr = {}

# un-irradiated SiPMs - No DCR
fRes0 = ROOT.TFile.Open('../plots/HPK_unirr_52deg_T5C_summaryPlots.root') 
g = fRes0.Get('g_deltaT_energyRatioCorr_bestTh_vs_bar_Vov%.02f_enBin01'%Vov)                                                                 
if (g != None): 
    fitTemp = ROOT.TF1('fitTemp','pol0', 0, 16)
    g.Fit(fitTemp,'QSR')
    tRes[0.0] = [fitTemp.GetParameter(0), fitTemp.GetParError(0)]
    #tRes[0.0] = [g.GetMean(2) , g.GetRMS(2)/math.sqrt(16) ]   
    err_dcr[0.0] = 0.001
    
dcr0 = 0
scaled_dcr = 0

for temp in temperatures:
    
    print '>>>>>> T = ', temp

    if (temp == 0): 
        Vov = 1.75
    else:
        Vov = 1.50

    fDCR = ROOT.TFile.Open('../plots/DCR_vs_Vov_HPK_1E13_52deg_T%dC.root'%temp)

    I_array = 0.5 * ( (fDCR.Get('g_I_vs_Vov_ch3')).Eval(Vov) + (fDCR.Get('g_I_vs_Vov_ch4')).Eval(Vov)) 
    VovEff  = Vov - (I_array*10)/1000 - (I_array/16*68)/1000
    print temp, Vov, VovEff
    dcr = 0.5 * ( (fDCR.Get('g_DCR_vs_VovEff_ch3')).Eval(VovEff) + (fDCR.Get('g_DCR_vs_VovEff_ch4')).Eval(VovEff))
    err_dcr[dcr] = 0.5 * abs(( (fDCR.Get('g_DCR_vs_VovEff_ch3')).Eval(VovEff) - (fDCR.Get('g_DCR_vs_VovEff_ch4')).Eval(VovEff)) )

    print 'dcr = ', dcr

    if (temp == 0 ): 
        dcr0 = dcr
    else:
        scaled_dcr = dcr0 / (2.*(0-temp)/10)
        print temp, VovEff, dcr, scaled_dcr
        #err_dcr[scaled_dcr] = err_dcr[dcr]
        #dcr = scaled_dcr
 

    # get tRes
    fRes = ROOT.TFile.Open('../plots/HPK_1E13_52deg_T%dC_summaryPlots.root'%temp)
    g = fRes.Get('g_deltaT_energyRatioCorr_bestTh_vs_bar_Vov%.02f_enBin01'%Vov)
    if (g == None): continue
    fitTemp = ROOT.TF1('fitTemp','pol0', 0, 16)
    g.Fit(fitTemp,'QSR')
    tRes[dcr] = [fitTemp.GetParameter(0), fitTemp.GetParError(0)]   
    #tRes[dcr] = [g.GetMean(2) , g.GetRMS(2)/math.sqrt(16) ]    


#sorted(tRes.items(), key=lambda x: x[1])
print tRes


g = ROOT.TGraphErrors()

for dcr in sorted(tRes):
    print dcr, err_dcr[dcr],tRes[dcr]
    tres = math.sqrt(tRes[dcr][0]*tRes[dcr][0] - tRes[0.0][0]*tRes[0.0][0])
    tres_err = 0
    if tres > 0:
        tres_err = 1./tres * math.sqrt(pow(tRes[dcr][0]*tRes[dcr][1],2)+pow(tRes[0.0][0]*tRes[0.0][1],2))
    else:
        tres_err = tRes[dcr][1]
    g.SetPoint(g.GetN(), dcr, tres)
#    g.SetPointError(g.GetN()-1, err_dcr[dcr], tRes[dcr][1])
    g.SetPointError(g.GetN()-1, err_dcr[dcr], tres_err)


fitFun = ROOT.TF1('fitFun','[0] * pow(x,[1])',0,100)
fitFun.SetNpx(1000)
fitFun.SetLineWidth(1)
fitFun.SetLineStyle(2)
fitFun.SetLineColor(2)
fitFun.SetParameter(0,30)
fitFun.SetParameter(1,0.5)
g.Fit(fitFun,'RS')


c = ROOT.TCanvas('c','c',600,600)
g.SetMarkerStyle(20)
g.SetMarkerSize(1)
g.GetXaxis().SetTitle('DCR [GHz]')
g.GetYaxis().SetTitle('#sigma_{DCR} [ps]')
g.GetYaxis().SetRangeUser(0,180)
g.GetXaxis().SetRangeUser(0,50)
g.Draw('ap')
#fitFun.Draw('same')
c.Update()
st = g.FindObject("stats")
st.SetX1NDC(0.7) 
st.SetX2NDC(0.93) 
st.SetY1NDC(0.8) #new x end position
st.SetY2NDC(0.93) #new x end position
c.Modified()

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.035)
latex.DrawLatex(0.2,0.88,'V_{ov} = %.02f V'%Vov)
latex.DrawLatex(0.2,0.83,'HPK arrays - 1e13 n_{eq}/cm^{2}')

c.SaveAs('/var/www/html/TOFHIR2A/MTDTB_CERN_Jul21/ModuleCharacterization/c_tResDCR_vs_DCR_Vov%0.02f_HPK_1E13.png'%Vov)
c.SaveAs('/var/www/html/TOFHIR2A/MTDTB_CERN_Jul21/ModuleCharacterization/c_tResDCR_vs_DCR_Vov%0.02f_HPK_1E13.pdf'%Vov)


raw_input('ok?')
