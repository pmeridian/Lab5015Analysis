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
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.05,'X')
ROOT.gStyle.SetLabelSize(0.05,'Y')
ROOT.gStyle.SetTitleSize(0.05,'X')
ROOT.gStyle.SetTitleSize(0.05,'Y')
ROOT.gStyle.SetTitleOffset(1.0,'X')
ROOT.gStyle.SetTitleOffset(1.2,'Y')                                                                                                                
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning 

def Gain(vov, irr):
    gain = (50738.5+95149*vov)  
    if (irr == '2E14'):
        gain = gain*0.7
    return gain

def DCR(I, gain):
    return (I*0.001/16/(1.602E-19*gain)/1000000000)


irr = '1E13'
#irr = '2E14'

Vbd = {'ch3': 37.55, 
       'ch4': 37.55}

inputDir = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2_pedSubtraction_tb/log_stressTest/ASIC0_HPK_non-irr__ASIC2_HPK_1E13__Tsipm0C__52deg/'
outDir = '/var/www/html/TOFHIR2A/MTDTB_CERN_Jul21/ModuleCharacterization/HPK_1E13_52deg_T0C_summaryPlots/DCR/'
outFileName = '../plots/DCR_vs_Vov_HPK_1E13_52deg_T0C.root'                  
outfile = ROOT.TFile(outFileName, 'RECREATE' )  

if (irr == '2E14'):
    Vbd = {'ch3': 38.15,
           'ch4': 38.50} 
    inputDir = '/home/cmsdaq/DAQ/tofhir/sw_daq_tofhir2_pedSubtraction_tb/log_stressTest/ASIC0_HPK_non-irr__ASIC2_HPK_2E14__Tsipm-40C__52deg/'
    outDir = '/var/www/html/TOFHIR2A/MTDTB_CERN_Jul21/ModuleCharacterization/HPK_2E14_52deg_T-40C_summaryPlots/DCR'
    outFileName = '../plots/DCR_vs_Vov_HPK_2E14_52deg_T-40C.root'                  
    outfile = ROOT.TFile(outFileName, 'RECREATE' )  

if (os.path.isdir(outDir) == False):
    os.system('mkdir %s'%outDir)  

fileList = os.listdir(inputDir)

#print fileList

bvs = {'ch3': [],
       'ch4': []}

fname = {}

for f in fileList:
    vbias = f.split('_')[2][2:7]
    if ('ch3' in f): 
        bvs['ch3'].append( f.split('_')[2][2:7])
        fname['ch3', vbias] = f 
    if ('ch4' in f):
        bvs['ch4'].append( f.split('_')[2][2:7])
        fname['ch4', vbias] = f 

print bvs
print fname

g_I_vs_Vov = {}
g_DCR_vs_Vov = {}

g_I_vs_VovEff = {}
g_DCR_vs_VovEff = {}

for ch in ['ch3','ch4']:
    g_I_vs_Vov[ch] = ROOT.TGraphErrors()
    g_DCR_vs_Vov[ch] = ROOT.TGraphErrors()
    g_I_vs_VovEff[ch] = ROOT.TGraphErrors()
    g_DCR_vs_VovEff[ch] = ROOT.TGraphErrors()
    g_I_vs_VovEff[ch].SetMarkerColor(ROOT.kRed) 
    g_DCR_vs_VovEff[ch].SetMarkerColor(ROOT.kRed) 
    if (ch == 'ch3'): 
            g_I_vs_Vov[ch].SetMarkerStyle(20)
            g_DCR_vs_Vov[ch].SetMarkerStyle(20)
            g_I_vs_VovEff[ch].SetMarkerStyle(20)
            g_DCR_vs_VovEff[ch].SetMarkerStyle(20)
    if (ch == 'ch4'): 
            g_I_vs_Vov[ch].SetMarkerStyle(24)
            g_DCR_vs_Vov[ch].SetMarkerStyle(24)
            g_I_vs_VovEff[ch].SetMarkerStyle(24)
            g_DCR_vs_VovEff[ch].SetMarkerStyle(24)
    for bv in bvs[ch]:
        f = ROOT.TFile.Open(inputDir+'/'+fname[ch,bv])
        g_I = f.Get('g_I')
        lastGoodPoint = g_I.GetN()-1
        for j in range(0, g_I.GetN()-1):
            if ( (g_I.GetY()[j]) < 1000.):
                lastGoodPoint = j
                break
        vov = float(bv)-Vbd[ch]

        fitFun = ROOT.TF1('fitFun','pol0',0,10000)
        fitFun.SetLineColor(2)
        fitFun.SetRange(0,g_I.GetX()[lastGoodPoint])
        g_I.Fit('fitFun','QRS')
        ctemp = ROOT.TCanvas('ctemp','ctemp')
        g_I.SetMarkerStyle(20)
        g_I.GetYaxis().SetRangeUser(0, fitFun.GetParameter(0)*1.5)
        g_I.GetYaxis().SetTitle('I_array [#muA]')
        g_I.GetXaxis().SetTitle('time elapsed[s]')
        g_I.Draw('ap')
        ctemp.SaveAs(outDir+'/Iarray_vs_time_%s_Vov%.02f.png'%(ch, vov))

        I_array = fitFun.GetParameter(0)/1000. # uA --> mA
        I_array_err = fitFun.GetParError(0)/1000.
        
        gain = Gain(vov, irr)
        dcr  = DCR(I_array, gain)
        dcr_err = 0.5*(DCR(I_array+I_array_err, gain)-DCR(I_array-I_array_err, gain))
        print 'Vov = ', vov, '   DCR = ', dcr

        g_I_vs_Vov[ch].SetPoint(g_I_vs_Vov[ch].GetN(), vov, I_array)
        g_I_vs_Vov[ch].SetPointError(g_I_vs_Vov[ch].GetN()-1, 0,  I_array_err)

        g_DCR_vs_Vov[ch].SetPoint(g_DCR_vs_Vov[ch].GetN(), vov, dcr)
        g_DCR_vs_Vov[ch].SetPoint(g_DCR_vs_Vov[ch].GetN(), vov, dcr_err)

        dVdrop = (I_array * 10)/1000. + (I_array/16 * 68)/1000.
        vov_eff = vov - dVdrop# deltaVdrop = I_array/R with R = 10 ohm + I_array/16*68 ohm
        gain = Gain(vov_eff, irr)
        dcr  = DCR(I_array, gain)
        dcr_err = 0.5*(DCR(I_array+I_array_err, gain)-DCR(I_array-I_array_err, gain))
        print 'I_array = ', I_array, '   deltaV_Drop = ', dVdrop,  '  Vov_eff = ', vov_eff, '   DCR_eff = ', dcr, '  dcr_err =', dcr_err

        g_I_vs_VovEff[ch].SetPoint(g_I_vs_VovEff[ch].GetN(), vov_eff, I_array)
        g_I_vs_VovEff[ch].SetPointError(g_I_vs_VovEff[ch].GetN()-1, 0, fitFun.GetParError(0)/1000.)

        g_DCR_vs_VovEff[ch].SetPoint(g_DCR_vs_VovEff[ch].GetN(), vov_eff, dcr)
        g_DCR_vs_VovEff[ch].SetPointError(g_DCR_vs_VovEff[ch].GetN()-1, 0, dcr_err)



nn = g_I_vs_Vov['ch3'].GetN()

c0 = ROOT.TCanvas()
c0.SetGridx()
c0.SetGridy()
hdummy = ROOT.TH2F('hdummy','',100,0,3.0,100,0,g_I_vs_Vov['ch3'].GetY()[nn-1]*1.5 )
hdummy.GetXaxis().SetTitle('V_{OV} [V]')
hdummy.GetYaxis().SetTitle('I_{array} (mA)')
hdummy.Draw()
g_I_vs_Vov['ch3'].Draw('psame')
g_I_vs_Vov['ch4'].Draw('psame')


c1 = ROOT.TCanvas()
c1.SetGridx()
c1.SetGridy()
hdummy.Draw()
g_I_vs_VovEff['ch3'].Draw('psame')
g_I_vs_VovEff['ch4'].Draw('psame')


c2 = ROOT.TCanvas()
c2.SetGridx()
c2.SetGridy()
hdummy2 = ROOT.TH2F('hdummy2','',100,0,3.0,100,0,g_DCR_vs_VovEff['ch3'].GetY()[g_DCR_vs_VovEff['ch3'].GetN()-1]*1.5)
hdummy2.GetXaxis().SetTitle('V_{OV} [V]')
hdummy2.GetYaxis().SetTitle('DCR [GHz]')
hdummy2.Draw()
g_DCR_vs_Vov['ch3'].Draw('psame')
g_DCR_vs_Vov['ch4'].Draw('psame')

c3 = ROOT.TCanvas()
c3.SetGridx()
c3.SetGridy()
hdummy2.Draw()
g_DCR_vs_VovEff['ch3'].Draw('psame')
g_DCR_vs_VovEff['ch4'].Draw('psame')



outfile.cd()
g_DCR_vs_Vov['ch3'].Write('g_DCR_vs_Vov_ch3')
g_DCR_vs_Vov['ch4'].Write('g_DCR_vs_Vov_ch4')
g_DCR_vs_VovEff['ch3'].Write('g_DCR_vs_VovEff_ch3')
g_DCR_vs_VovEff['ch4'].Write('g_DCR_vs_VovEff_ch4')
g_I_vs_Vov['ch3'].Write('g_I_vs_Vov_ch3')
g_I_vs_Vov['ch4'].Write('g_I_vs_Vov_ch4')
g_I_vs_VovEff['ch3'].Write('g_I_vs_VovEff_ch3')
g_I_vs_VovEff['ch4'].Write('g_I_vs_VovEff_ch4')
outfile.Close() 

c0.SaveAs(outDir+'/Iarray_vs_Vov.png')
c0.SaveAs(outDir+'/Iarray_vs_Vov.pdf')

c1.SaveAs(outDir+'/Iarray_vs_VovEff.png')
c1.SaveAs(outDir+'/Iarray_vs_VovEff.pdf')

c2.SaveAs(outDir+'/DCR_vs_Vov.png')
c2.SaveAs(outDir+'/DCR_vs_Vov.pdf')

c3.SaveAs(outDir+'/DCR_vs_VovEff.png')
c3.SaveAs(outDir+'/DCR_vs_VovEff.pdf')

raw_input('OK?')

