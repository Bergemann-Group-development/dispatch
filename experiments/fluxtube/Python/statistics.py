# -*- coding: utf-8 -*-
from read_file import read_file
import matplotlib.pyplot as pl
import numpy as np
import os
from array import array
from ROOT import *

gROOT.Reset()

# Set the showing of statistics and fitting results (0 means off, 1111 means all on):
gStyle.SetOptStat(1111)
gStyle.SetOptFit(1111)
#-------------------------------------------------------------------------------
mean = array('f')
mean_e = array('f')
x_e = array('f')
x = array('f')
for i in xrange(10):
    data_format=[1,i]
    
    data_mhd = "pure_mhd_1_patch"
    name_mhd="../../../../plot_data/"+data_mhd+"/{:04d}_{:04d}"
    file_mhd=name_mhd.format(data_format[0],data_format[1])
    f_mhd=read_file(file_mhd,True)
    d_slice_mhd = f_mhd.bz[f_mhd.n[0]/2,:,0]
    
    data_pic = "pure_pic_1_patch"
    name_pic="../../../../plot_data/"+data_pic+"/{:04d}_{:04d}"
    file_pic=name_pic.format(data_format[0],data_format[1])
    f_pic=read_file(file_pic,True)
    d_slice_pic = f_pic.bz[f_pic.n[0]/2,:,0]
    
    perc_diff = abs(d_slice_mhd-d_slice_pic)/d_slice_mhd
    
    #pl.figure(1)
    #pl.plot(perc_diff,'r.')
    #pl.show()
    
    Hist_Diff         = TH1F("Hist_Diff"        , "Hist_Diff"        , 30, 0.0, max(perc_diff)+1)
    for j in xrange(len(perc_diff)):
        Hist_Diff.Fill(perc_diff[j])
    #canvas1 = TCanvas( "canvas1", "Hist_Diff", 20, 20, 1200, 500)
    #Hist_Diff.Draw()
    #canvas1.Update()
    
    mean.append(Hist_Diff.GetMean())
    mean_e.append(Hist_Diff.GetRMS()/sqrt(10))
    x.append(i)
    x_e.append(0)
   # raw_input('Press Enter to exit')

#pl.figure(2)
#pl.errorbar('r.')
#pl.title('Evolution in time of the mean of the difference in density between MHD and PIC')
#pl.xlabel('Time [10e-1 s]')
#pl.ylabel('Mean of density difference [%]')
#pl.show()

gStyle.SetOptStat("")
gStyle.SetOptFit(0)
gStyle.SetStatFormat("8.6g")
gStyle.SetStatX(0.49)
gStyle.SetStatY(0.88)

graph_diff = TGraphErrors(10, x, mean, x_e, mean_e)

canvas1 = TCanvas( "canvas1", "Hist_Diff", 20, 20, 1200, 500)
graph_diff.Draw("AP")

canvas1.Update()


























