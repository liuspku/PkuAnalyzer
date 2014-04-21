#! /usr/bin/env python
import os
import glob
import math
import array
import ROOT
import ntpath
import sys
import subprocess
from subprocess import Popen
from optparse import OptionParser

from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHist,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite

ROOT.gSystem.Load("PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load("PDFs/Util_cxx.so")
#ROOT.gSystem.Load("PDFs/RooRelBWRunningWidth_cxx.so")
ROOT.gSystem.Load("PDFs/HWWLVJRooPdfs_cxx.so")
from ROOT import draw_error_band, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error, RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf, RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf
#

###############################
## Tools ##
###############################

## Set basic TDR style for canvas, pad ..etc ..
def setTDRStyle(): 
    print "setting TDR style"
    tdrStyle =TStyle("tdrStyle","Style for P-TDR");
    #For the canvas:
    tdrStyle.SetCanvasBorderMode(0);
    tdrStyle.SetCanvasColor(kWhite);
    tdrStyle.SetCanvasDefH(600); #Height of canvas
    tdrStyle.SetCanvasDefW(600); #Width of canvas
    tdrStyle.SetCanvasDefX(0); #POsition on screen
    tdrStyle.SetCanvasDefY(0);

    #For the Pad:
    tdrStyle.SetPadBorderMode(0);
    tdrStyle.SetPadColor(kWhite);
    tdrStyle.SetPadGridX(False);
    tdrStyle.SetPadGridY(False);
    tdrStyle.SetGridColor(0);
    tdrStyle.SetGridStyle(3);
    tdrStyle.SetGridWidth(1);

    #For the frame:
    tdrStyle.SetFrameBorderMode(0);
    tdrStyle.SetFrameBorderSize(1);
    tdrStyle.SetFrameFillColor(0);
    tdrStyle.SetFrameFillStyle(0);
    tdrStyle.SetFrameLineColor(1);
    tdrStyle.SetFrameLineStyle(1);
    tdrStyle.SetFrameLineWidth(1);

    #For the histo:
    tdrStyle.SetHistLineColor(1);
    tdrStyle.SetHistLineStyle(0);
    tdrStyle.SetHistLineWidth(1);
    tdrStyle.SetEndErrorSize(2);
    tdrStyle.SetErrorX(0.);
    tdrStyle.SetMarkerStyle(20);

    #For the fit/function:
    tdrStyle.SetOptFit(1);
    tdrStyle.SetFitFormat("5.4g");
    tdrStyle.SetFuncColor(2);
    tdrStyle.SetFuncStyle(1);
    tdrStyle.SetFuncWidth(1);

    #For the date:
    tdrStyle.SetOptDate(0);

    #For the statistics box:
    tdrStyle.SetOptFile(0);
    tdrStyle.SetOptStat(0); #To display the mean and RMS:
    tdrStyle.SetStatColor(kWhite);
    tdrStyle.SetStatFont(42);
    tdrStyle.SetStatFontSize(0.025);
    tdrStyle.SetStatTextColor(1);
    tdrStyle.SetStatFormat("6.4g");
    tdrStyle.SetStatBorderSize(1);
    tdrStyle.SetStatH(0.1);
    tdrStyle.SetStatW(0.15);

    #Margins:
    tdrStyle.SetPadTopMargin(0.05);
    tdrStyle.SetPadBottomMargin(0.13);
    tdrStyle.SetPadLeftMargin(0.18);
    tdrStyle.SetPadRightMargin(0.06);

    #For the Global title:
    tdrStyle.SetOptTitle(0);
    tdrStyle.SetTitleFont(42);
    tdrStyle.SetTitleColor(1);
    tdrStyle.SetTitleTextColor(1);
    tdrStyle.SetTitleFillColor(10);
    tdrStyle.SetTitleFontSize(0.05);

    #For the axis titles:
    tdrStyle.SetTitleColor(1, "XYZ");
    tdrStyle.SetTitleFont(42, "XYZ");
    tdrStyle.SetTitleSize(0.03, "XYZ");
    tdrStyle.SetTitleXOffset(0.9);
    tdrStyle.SetTitleYOffset(1.5);

    #For the axis labels:
    tdrStyle.SetLabelColor(1, "XYZ");
    tdrStyle.SetLabelFont(42, "XYZ");
    tdrStyle.SetLabelOffset(0.007, "XYZ");
    tdrStyle.SetLabelSize(0.03, "XYZ");

    #For the axis:
    tdrStyle.SetAxisColor(1, "XYZ");
    tdrStyle.SetStripDecimals(kTRUE);
    tdrStyle.SetTickLength(0.03, "XYZ");
    tdrStyle.SetNdivisions(510, "XYZ");
    tdrStyle.SetPadTickX(1); #To get tick marks on the opposite side of the frame
    tdrStyle.SetPadTickY(1);

    #Change for log plots:
    tdrStyle.SetOptLogx(0);
    tdrStyle.SetOptLogy(0);
    tdrStyle.SetOptLogz(0);

    #Postscript options:
    tdrStyle.SetPaperSize(20.,20.);
    tdrStyle.cd();

#### Method to make a RooAbsPdf giving label, model name, spectrum, if it is mc or not and a constraint list for the parameters          

#label+self.categoryLabel:  label 
#workspace4fit_: workspace 
#def make_Pdf( label, fit_config, mass_spectrum="_obs0_variable", ConstraintsList=[],ismc = 0):
def make_Pdf( label, workspace, fit_config, mass_spectrum):
    #if TString(mass_spectrum).Contains("_obs0_variable"): rrv_x = workspace.var("rrv_mass_j");
    #if TString(mass_spectrum).Contains("_limit_variable"): rrv_x = workspace.var("rrv_mass_lvj");
    rrv_x = workspace.var("rrv"+mass_spectrum);

    in_model_name=fit_config[1];

    # W mass: 80.385
    if in_model_name == "Voig":
        print "########### Voigtian Pdf for mJ ############"
        rrv_mean_voig=RooRealVar("rrv_mean_voig"+label,"rrv_mean_voig"+label,84,78,88);
        rrv_width_voig=RooRealVar("rrv_width_voig"+label,"rrv_width_voig"+label,7.,1,40);
        rrv_sigma_voig=RooRealVar("rrv_sigma_voig"+label,"rrv_sigma_voig"+label,5,0.01,20);
        model_pdf = RooVoigtian("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,rrv_mean_voig,rrv_width_voig,rrv_sigma_voig);

    # Higgs mass 600-1000
    if in_model_name == "Voig_v1":
        print "########### Voigtian Pdf for Higgs limit_variable ############"
        rrv_mean_voig=RooRealVar("rrv_mean_voig"+label,"rrv_mean_voig"+label,650,550,1200);
        rrv_width_voig=RooRealVar("rrv_width_voig"+label,"rrv_width_voig"+label,100.,10,600);
        rrv_sigma_voig=RooRealVar("rrv_sigma_voig"+label,"rrv_sigma_voig"+label,200,10,400);
        model_pdf = RooVoigtian("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,rrv_mean_voig,rrv_width_voig,rrv_sigma_voig);

    # Bulk mass 600-1000
    if in_model_name == "Voig_v2":
        print "########### Voigtian Pdf for Higgs limit_variable ############"
        rrv_mean_voig=RooRealVar("rrv_mean_voig"+label,"rrv_mean_voig"+label,1000,900,1100);#
        rrv_width_voig=RooRealVar("rrv_width_voig"+label,"rrv_width_voig"+label,2.5,0,10);
        rrv_sigma_voig=RooRealVar("rrv_sigma_voig"+label,"rrv_sigma_voig"+label,40,10,80);

        model_pdf = RooVoigtian("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,rrv_mean_voig,rrv_width_voig,rrv_sigma_voig);

    ## BW for the W mass peak 
    if in_model_name == "BW":            
        print "########### BW Pdf for obs0_variable fit ############"
        rrv_mean_BW=RooRealVar("rrv_mean_BW"+label,"rrv_mean_BW"+label,84,78, 88);
        rrv_width_BW=RooRealVar("rrv_width_BW"+label,"rrv_width_BW"+label,20,1,40);
        model_pdf = RooBreitWigner("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,rrv_mean_BW,rrv_width_BW);

    ## BW relativistic for Higgs mass covoluted with CB 
    if in_model_name == "BWRUNxGausErf":

        print "########### BWRUNxGausErf Pdf for Higgs lvj ############"
        rrv_mean_BWRUN = RooRealVar("rrv_mean_BWRUN"+label,"rrv_mean_BWRUN"+label,1000,950,1050);
        rrv_width_BWRUN= RooRealVar("rrv_width_BWRUN"+label,"rrv_width_BWRUN"+label,200,100,370);

        bwrun = RooBWRunPdf("bwrun"+label+mass_spectrum,"bwrun"+label+mass_spectrum,rrv_x, rrv_mean_BWRUN, rrv_width_BWRUN);

        rrv_mean_cb  = RooRealVar("rrv_mean_cb"+label,"rrv_mean_cb"+label,0);
        rrv_sigma_cb = RooRealVar("rrv_sigma_cb"+label,"rrv_sigma_cb"+label,50,10,300);
        cbshape      = RooGaussian("cbshape"+label,"cbshape"+label, rrv_x,rrv_mean_cb,rrv_sigma_cb);
        fft          = RooFFTConvPdf("fft"+label+mass_spectrum,"fft"+label+mass_spectrum, rrv_x, bwrun, cbshape);

        rrv_offset_erf = RooRealVar("rrv_offset_erf"+label,"rrv_offset_erf"+label,450)#,350,550);
        rrv_width_erf = RooRealVar("rrv_width_erf"+label,"rrv_width_erf"+label,50)#,10,250);
        erf = RooGenericPdf("erf"+label+mass_spectrum,"erf"+label+mass_spectrum, "(1.+TMath::Erf((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset_erf.GetName(), rrv_width_erf.GetName()), RooArgList(rrv_x,rrv_offset_erf,rrv_width_erf) )

        model_pdf = RooProdPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, fft, erf );

    ##  Voig for W mass peak
    if in_model_name == "2Voig":

        print "########### Double Voigtian for obs0_variable fit ############"
        rrv_mean_voig    = RooRealVar("rrv_mean_voig"+label,"rrv_mean_voig"+label,84,78,88);#W mass 80.385
        rrv_shift_2Voig  = RooRealVar("rrv_shift_2Voig"+label,"rrv_shift_2Voig"+label,10.8026)# Z mass: 91.1876; shift=91.1876-80.385=10.8026
        rrv_mean_shifted = RooFormulaVar("rrv_mean_voig2"+label,"@0+@1",RooArgList(rrv_mean_voig,rrv_shift_2Voig));

        rrv_width_voig = RooRealVar("rrv_width_voig"+label,"rrv_width_voig"+label,16.,6,26);
        rrv_sigma_voig = RooRealVar("rrv_sigma_voig"+label,"rrv_sigma_voig"+label,5.,0.,10.);

        rrv_frac = RooRealVar("rrv_frac"+label,"rrv_frac"+label,0.8,0.5,1.);

        model_voig1 = RooVoigtian("model_voig1"+label+mass_spectrum,"model_voig1"+label+mass_spectrum, rrv_x,rrv_mean_voig,rrv_width_voig,rrv_sigma_voig);

        model_voig2 = RooVoigtian("model_voig2"+label+mass_spectrum,"model_voig2"+label+mass_spectrum, rrv_x,rrv_mean_shifted,rrv_width_voig,rrv_sigma_voig);
        model_pdf = RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, RooArgList(model_voig1,model_voig2), RooArgList(rrv_frac));

    ## Gaus for the W peak
    if in_model_name == "Gaus":
        print "########### Gaus for W peak  ############"
        rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label,"rrv_mean_gaus"+label,84,78,88);
        rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label,"rrv_sigma_gaus"+label,7,1,15);
        model_pdf = RooGaussian("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

    ## Gaus for the higgs lineshape
    if in_model_name == "Gaus_v1":

        print "########### Gaus for Higgs limit_variable ############"
        rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label,"rrv_mean_gaus"+label,920,900,1000);
        rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label,"rrv_sigma_gaus"+label,200,100,300);

        model_pdf = RooGaussian("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

    if in_model_name == "BifurGaus_v1":

        print "########### BifurGaus for Higgs limit_variable ############"
        rrv_mean_gaus   = RooRealVar("rrv_mean_gaus"+label,"rrv_mean_gaus"+label,920,900,1000);
        rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label,"rrv_sigma1_gaus"+label,200,100,300);
        rrv_sigma2_gaus = RooRealVar("rrv_sigma2_gaus"+label,"rrv_sigma2_gaus"+label,200,100,300);

        model_pdf = RooBifurGauss("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,rrv_mean_gaus,rrv_sigma1_gaus,rrv_sigma2_gaus);

    ## Crystal Ball for the W mass peak
    if in_model_name == "CB":
        print "########### Cystal Ball for obs0_variable fit ############"
        rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label,"rrv_mean_CB"+label,84,78,88);
        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label,"rrv_sigma_CB"+label,7,4,10);
        rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label,"rrv_alpha_CB"+label,-2,-4,-0.5);
        rrv_n_CB     = RooRealVar("rrv_n_CB"+label,"rrv_n_CB"+label,2,0.,4);
        model_pdf = RooCBShape("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,rrv_mean_CB,rrv_sigma_CB,rrv_alpha_CB,rrv_n_CB);

    ## Sum of two CB 
    if in_model_name == "SCB_v1":
        print "########### Cystal Ball + Crystall Ball ############"
        rrv_mean_SCB   = RooRealVar("rrv_mean_SCB"+label,"rrv_mean_SCB"+label,800,550,1000);
        rrv_sigma_SCB  = RooRealVar("rrv_sigma_SCB"+label,"rrv_sigma_SCB"+label,70,40,300);
        rrv_alpha1_SCB = RooRealVar("rrv_alpha1_SCB"+label,"rrv_alpha1_SCB"+label,-2,-4,-0.5);
        rrv_alpha2_SCB = RooRealVar("rrv_alpha2_SCB"+label,"rrv_alpha2_SCB"+label,2,0.5,4);
        rrv_n1_SCB     = RooRealVar("rrv_n1_SCB"+label,"rrv_n1_SCB"+label,2,0.,4);
        rrv_n2_SCB     = RooRealVar("rrv_n2_SCB"+label,"rrv_n2_SCB"+label,2,0.,4);
        frac           = RooRealVar("rrv_frac_SSCB"+label,"rrv_frac_SSCB"+label,0.5)
        scb1 = RooCBShape("model_pdf_scb1"+label+mass_spectrum,"model_pdf_scb1"+label+mass_spectrum, rrv_x,rrv_mean_SCB,rrv_sigma_SCB,rrv_alpha1_SCB,rrv_n1_SCB);
        scb2 = RooCBShape("model_pdf_scb2"+label+mass_spectrum,"model_pdf_scb2"+label+mass_spectrum, rrv_x,rrv_mean_SCB,rrv_sigma_SCB,rrv_alpha2_SCB,rrv_n2_SCB);
        model_pdf = RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,RooArgList(scb1,scb2),RooArgList(frac))

    ## Double Gaus for Bulk GR signal --> narrow width
    if in_model_name == "2Gaus_sig":

        print "########### Double Gauss for Bulk GR ############"
        label_tstring=TString(label);

        rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label,"rrv_mean1_gaus"+label, 1000, 900, 1100);
        rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label,"rrv_sigma1_gaus"+label,50,20,120);
        gaus1 = RooGaussian("gaus1"+label,"gaus1"+label, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

        rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label,"rrv_deltamean_gaus"+label,0,-50,50);
        rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
        rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label,"rrv_scalesigma_gaus"+label,1,0.,10.);
        rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
        gaus2 = RooGaussian("gaus2"+label,"gaus2"+label, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

        rrv_frac = RooRealVar("rrv_frac"+label,"rrv_frac"+label,0.5,0.,1.);

        model_pdf =RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac),1)

    ## Crystal  ball shape for Bulk GR samples and higgs 
    if in_model_name == "CB_v1":
        print "########### Crystal Ball for Higgs and  Bulk GR  limit_variable ############"
        label_tstring=TString(label);

        rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label,"rrv_mean_CB"+label,920,800,1150);
        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label,"rrv_sigma_CB"+label,200,100,300);
        rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label,"rrv_alpha_CB"+label,4,1,5);
        rrv_n_CB     = RooRealVar("rrv_n_CB"+label,"rrv_n_CB"+label,20.,10,40);

        model_pdf = RooCBShape("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,rrv_mean_CB,rrv_sigma_CB,rrv_alpha_CB,rrv_n_CB);

    ## Crystal  ball shape for Bulk GR samples and higgs 
    if in_model_name == "BWCB":

        print "########### Crystal Ball x Breit Wigner for Bulk Graviton width ############"
        label_tstring=TString(label);

        rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label,"rrv_mean_BW"+label, 1000);
        rrv_width_BW = RooRealVar("rrv_sigma_BW"+label,"rrv_sigma_BW"+label, 50);
        rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label,"rrv_mean_CB"+label,0.,0.,50.);
        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label,"rrv_sigma_CB"+label,50,0,200);

        rrv_mean_BW.setConstant(kTRUE);
        rrv_width_BW.setConstant(kTRUE);

        rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label,"rrv_alpha_CB"+label,2,0,4);
        rrv_n_CB     = RooRealVar("rrv_n_CB"+label,"rrv_n_CB"+label,1.,0.,4.);

        bw      = RooBreitWigner("bw"+label,"bw"+label, rrv_x,rrv_mean_BW,rrv_width_BW);
        cbshape = RooCBShape("cbshape"+label,"cbshape"+label, rrv_x,rrv_mean_CB,rrv_sigma_CB,rrv_alpha_CB,rrv_n_CB);

        model_pdf = RooFFTConvPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x, cbshape, bw);

    if in_model_name == "ArgusBW_v1":

        label_tstring=TString(label);
        rrv_width_BW = RooRealVar("rrv_width_BW"+label,"rrv_width_BW"+label,100,50,600);
        rrv_m0_Argus = RooRealVar("rrv_m0_Argus"+label,"rrv_m0_Argus"+label, 950 );
        rrv_c_Argus  = RooRealVar("rrv_c_Argus"+label,"rrv_c_Argus"+label,-1,-2,-1e-1);
        rrv_frac     = RooRealVar("rrv_frac"+label,"rrv_frac"+label,0.5,0.0,1.);

        bw    = RooBreitWigner("bw"+label,"bw"+label, rrv_x,rrv_m0_Argus,rrv_width_BW);
        argus = RooArgusBG("argus"+label,"argus"+label, rrv_x, rrv_m0_Argus,rrv_c_Argus);
        model_pdf = RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, RooArgList(bw,argus), RooArgList(rrv_frac));

    if in_model_name == "CBBW": # FFT: BreitWigner*CBShape
        print "########### Crystal Ball x Breit Wigner for W mass peak ############"
        rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label,"rrv_mean_CB"+label,84.0,78,88);
        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label,"rrv_sigma_CB"+label,7,4,10);
        rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label,"rrv_alpha_CB"+label,-2,-4,-1);
        rrv_n_CB     = RooRealVar("rrv_n_CB"+label,"rrv_n_CB"+label,0.5,0.,2);
        rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label,"rrv_mean_BW"+label,0);
        rrv_width_BW = RooRealVar("rrv_width_BW"+label,"rrv_width_BW"+label,10,5,20);
        cbshape      = RooCBShape("cbshape"+label,"cbshape"+label, rrv_x,rrv_mean_CB,rrv_sigma_CB,rrv_alpha_CB,rrv_n_CB);
        bw           = RooBreitWigner("bw"+label,"bw"+label, rrv_x,rrv_mean_BW,rrv_width_BW);
        model_pdf    = RooFFTConvPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x, cbshape, bw);

    if in_model_name == "LDGaus": # FFT: Landau*Gaus
        print "########### Landau x Breit Wigner for W mass peak ############"
        rrv_mean_landau  = RooRealVar("rrv_mean_landau"+label,"rrv_mean_landau"+label,84.0,78,88);
        rrv_sigma_landau = RooRealVar("rrv_sigma_landau"+label,"rrv_sigma_landau"+label,7,4,10);
        rrv_mean_gaus    = RooRealVar("rrv_mean_gaus"+label,"rrv_mean_gaus"+label,0);
        rrv_sigma_gaus   = RooRealVar("rrv_sigma_gaus"+label,"rrv_sigma_gaus"+label,16,10,20);
        landau           = RooLandau("landau"+label,"landau"+label, rrv_x,rrv_mean_landau,rrv_sigma_landau);
        gaus             = RooBreitWigner("gaus"+label,"gaus"+label, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);
        model_pdf        = RooFFTConvPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x, landau, gaus);

    '''
    ## Crystal  ball shape for Bulk GR samples and higgs 
    if in_model_name == "DoubleCB_v1":
        label_tstring=TString(label);
        print "########### Double CB for Bulk graviton limit_variable ############"

        rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label,"rrv_mean_CB"+label,700,550,2500);
        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label,"rrv_sigma_CB"+label, 50,20 ,120);
        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label,"rrv_n1_CB"+label, 10.,0.01,35);
        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label,"rrv_alpha2_CB"+label,3.,0.5,6.);
        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label,"rrv_n2_CB"+label,20.,0.01,35);
        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label,"rrv_alpha1_CB"+label,3,0.5,6.);


        rrv_mean_scale_p1 = RooRealVar("CMS_sig_p1_jes","CMS_sig_p1_jes",0);
        rrv_mean_scale_p1.setConstant(kTRUE);
        rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_em","CMS_sig_p1_scale_em",0);
        rrv_mean_scale_p2.setConstant(kTRUE);


        rrv_mean_scale_X1 = RooRealVar("rrv_mean_shift_scale_lep"+label,"rrv_mean_shift_scale_lep"+label,float(self.mean_signal_uncertainty_lep_scale));
        rrv_mean_scale_X1.setConstant(kTRUE);
        rrv_mean_scale_X2 = RooRealVar("rrv_mean_shift_scale_jes"+label,"rrv_mean_shift_scale_jes"+label,float(self.mean_signal_uncertainty_jet_scale));
        rrv_mean_scale_X2.setConstant(kTRUE);

        rrv_total_mean_CB = RooFormulaVar("rrv_total_mean_CB"+label,"@0*(1+@1*@2)*(1+@3*@4)", RooArgList(rrv_mean_CB,rrv_mean_scale_p1,rrv_mean_scale_X1,rrv_mean_scale_p2,rrv_mean_scale_X2));

        rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_em","CMS_sig_p2_scale_em",0);
        rrv_sigma_scale_p1.setConstant(kTRUE);

        rrv_sigma_scale_p2 = RooRealVar("CMS_sig_p2_jer","CMS_sig_p2_jer",0);
        rrv_sigma_scale_p3 = RooRealVar("CMS_sig_p2_jes","CMS_sig_p2_jes",0);
        rrv_sigma_scale_p2.setConstant(kTRUE);
        rrv_sigma_scale_p3.setConstant(kTRUE);

        rrv_mean_sigma_X1 = RooRealVar("rrv_sigma_shift_lep_scale"+label,"rrv_sigma_shift_scale"+label,float(self.sigma_signal_uncertainty_lep_scale));
        rrv_mean_sigma_X2 = RooRealVar("rrv_sigma_shift_jes"+label,"rrv_sigma_shift_scale"+label,float(self.sigma_signal_uncertainty_jet_scale));
        rrv_mean_sigma_X3 = RooRealVar("rrv_sigma_shift_res"+label,"rrv_sigma_shift_res"+label,float(self.sigma_signal_uncertainty_jet_res));
        rrv_mean_sigma_X1.setConstant(kTRUE);
        rrv_mean_sigma_X2.setConstant(kTRUE);
        rrv_mean_sigma_X3.setConstant(kTRUE);

        rrv_total_sigma_CB = RooFormulaVar("rrv_total_sigma_CB"+label,"@0*(1+@1*@2)*(1+@3*@4)*(1+@5*@6)", RooArgList(rrv_sigma_CB,rrv_sigma_scale_p1,rrv_mean_sigma_X1,rrv_sigma_scale_p2,rrv_mean_sigma_X2,rrv_sigma_scale_p3,rrv_mean_sigma_X3));        

        model_pdf = ROOT.RooDoubleCrystalBall("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,rrv_total_mean_CB,rrv_total_sigma_CB,rrv_alpha1_CB,rrv_n1_CB,rrv_alpha2_CB,rrv_n2_CB);

    ## Crystal  ball shape for Bulk GR samples and higgs 
    if in_model_name == "BWDoubleCB":

        label_tstring=TString(label);
        print "########### Double CB x BW for Bulk graviton width ############"

        rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label,"rrv_mean_CB"+label,0.,-100,100);
        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label,"rrv_sigma_CB"+label,90,20,250);
        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label,"rrv_n1_CB"+label, 20.,0.01,105);
        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label,"rrv_alpha2_CB"+label,3.5,0.5,50.5);
        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label,"rrv_n2_CB"+label,20.,0.01,105);
        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label,"rrv_alpha1_CB"+label,3.5,0.5,50.5);
        rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label,"rrv_mean_BW"+label,2100);
        rrv_width_BW = RooRealVar("rrv_width_BW"+label,"rrv_width_BW"+label,630);


        ### fix the Breit-Wigner core to the generated one  
        rrv_mean_BW.setConstant(kTRUE);
        rrv_width_BW.setConstant(kTRUE);                    
        bw           = RooBreitWigner("bw"+label,"bw"+label, rrv_x,rrv_mean_BW,rrv_width_BW);

        ### Double Crystall ball term --> add parameters in order to do systematic on the signal shape inside combiner
        rrv_mean_scale_p1 = RooRealVar("CMS_sig_p1_jes","CMS_sig_p1_jes",0);  ## jes effect on the mean
        rrv_mean_scale_p1.setConstant(kTRUE);

        rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_em","CMS_sig_p1_scale_em",0);
        rrv_mean_scale_p2.setConstant(kTRUE);

        ## set the uncertainty value in other two independent variables 
        rrv_mean_scale_X1 = RooRealVar("rrv_mean_shift_scale_lep"+label,"rrv_mean_shift_scale_lep"+label,float(self.mean_signal_uncertainty_lep_scale));
        rrv_mean_scale_X1.setConstant(kTRUE);

        rrv_mean_scale_X2 = RooRealVar("rrv_mean_shift_scale_jes"+label,"rrv_mean_shift_scale_jes"+label,float(self.mean_signal_uncertainty_jet_scale));
        rrv_mean_scale_X2.setConstant(kTRUE);

        ### total mean
        rrv_total_mean_CB = RooFormulaVar("rrv_total_mean_CB"+label,"@0*(1+@1*@2)*(1+@3*@4)", RooArgList(rrv_mean_CB,rrv_mean_scale_p1,rrv_mean_scale_X1,rrv_mean_scale_p2,rrv_mean_scale_X2));

        ### lepton scale effect on the resolution 
        rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_em","CMS_sig_p2_scale_em",0);
        rrv_sigma_scale_p1.setConstant(kTRUE);

        ### jes and jer effect on the resolution             
        rrv_sigma_scale_p2 = RooRealVar("CMS_sig_p2_jer","CMS_sig_p2_jer",0);
        rrv_sigma_scale_p3 = RooRealVar("CMS_sig_p2_jes","CMS_sig_p2_jes",0);

        rrv_sigma_scale_p2.setConstant(kTRUE);
        rrv_sigma_scale_p3.setConstant(kTRUE);

        rrv_mean_sigma_X1 = RooRealVar("rrv_sigma_shift_lep_scale"+label,"rrv_sigma_shift_scale"+label,float(self.sigma_signal_uncertainty_lep_scale));
        rrv_mean_sigma_X2 = RooRealVar("rrv_sigma_shift_jes"+label,"rrv_sigma_shift_scale"+label,float(self.sigma_signal_uncertainty_jet_scale));
        rrv_mean_sigma_X3 = RooRealVar("rrv_sigma_shift_res"+label,"rrv_sigma_shift_res"+label,float(self.sigma_signal_uncertainty_jet_res));

        rrv_mean_sigma_X1.setConstant(kTRUE);
        rrv_mean_sigma_X2.setConstant(kTRUE);
        rrv_mean_sigma_X3.setConstant(kTRUE);

        ### total resolution 
        rrv_total_sigma_CB = RooFormulaVar("rrv_total_sigma_CB"+label,"@0*(1+@1*@2)*(1+@3*@4)*(1+@5*@6)", RooArgList(rrv_sigma_CB,rrv_sigma_scale_p1,rrv_mean_sigma_X1,rrv_sigma_scale_p2,rrv_mean_sigma_X2,rrv_sigma_scale_p3,rrv_mean_sigma_X3));        

        cbshape = ROOT.RooDoubleCrystalBall("DoubleCB"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,rrv_total_mean_CB,rrv_total_sigma_CB,rrv_alpha1_CB,rrv_n1_CB,rrv_alpha2_CB,rrv_n2_CB)

        ### numerical convolution via FFT
        model_pdf = RooFFTConvPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,bw,cbshape);
        model_pdf.setBufferFraction(1.0)
        '''

    ## ExpN pdf for W+jets bkg fit
    if in_model_name == "ExpN":

        print "########### ExpN funtion for W+jets limit_variable ############"
        rrv_c_ExpN = RooRealVar("rrv_c_ExpN"+label,"rrv_c_ExpN"+label,-3e-3,-1e-1,-1e-5);
        rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label,"rrv_n_ExpN"+label, 1e3, -1e2, 1e4);

        model_pdf = ROOT.RooExpNPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_c_ExpN, rrv_n_ExpN);


    ## levelled exp for W+jets bkg fit
    if in_model_name == "ExpTail":
        print "########### ExpTai = levelled exp funtion for W+jets limit_variable ############"
        rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label,"rrv_s_ExpTail"+label, 110,20,242);
        rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label,"rrv_a_ExpTail"+label, 2.9e-2,-1e-2,7.5e-2);

        model_pdf     = ROOT.RooExpTailPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_s_ExpTail, rrv_a_ExpTail);

    ## sum of two exponential 
    if in_model_name == "2Exp":
        print "########### 2Exp = levelled exp funtion for W+jets limit_variable ############"
        rrv_c0_2Exp   = RooRealVar("rrv_c0_2Exp"+label,"rrv_c0_2Exp"+label, -5e-3, -8e-3,-4e-3);
        rrv_c1_2Exp   = RooRealVar("rrv_c1_2Exp"+label,"rrv_c1_2Exp"+label, -1e-3, -4e-3,-1e-4);
        rrv_frac_2Exp = RooRealVar("rrv_frac_2Exp"+label,"rrv_frac_2Exp"+label, 0., 0., 1e-2);
        model_pdf  = ROOT.Roo2ExpPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_c0_2Exp,rrv_c1_2Exp,rrv_frac_2Exp);

    ## sum of two exponential 
    if in_model_name == "Exp" or in_model_name == "Exp_sr":
        print "########### Exp = levelled exp funtion for W+jets limit_variable ############"
        rrv_c_Exp = RooRealVar("rrv_c_Exp"+label,"rrv_c_Exp"+label,-0.05,-0.1,0.);
        model_pdf = ROOT.RooExponential("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_c_Exp);

    ## Erf times for obs0_variable spectrum
    if in_model_name == "ErfExp" :
        print "########### Erf*Exp for obs0_variable fit  ############"
        rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label,"rrv_c_ErfExp"+label,-0.05,-0.1,-1e-4);
        rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label,"rrv_offset_ErfExp"+label,60.,30.,120);
        rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label,"rrv_width_ErfExp"+label,30.,10, 60.);
        model_pdf         = ROOT.RooErfExpPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

    ## different initial values -> for limit_variable
    if in_model_name == "ErfExp_v1" :
        print "########### Erf*Exp for limit_variable fit  ############"
        rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label,"rrv_c_ErfExp"+label,-0.006,-0.1,0.);
        rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label,"rrv_offset_ErfExp"+label,450.,400.,550.);
        rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label,"rrv_width_ErfExp"+label,70.,10,100.);
        model_pdf         = ROOT.RooErfExpPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

    ## different initial values -> for limit_variable
    if in_model_name == "ErfExp_v2" : 
        print "########### Erf*Exp for limit_variable fit  ############"
        rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label,"rrv_c_ErfExp"+label,-0.005,-0.1,0.);
        rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label,"rrv_offset_ErfExp"+label,450.,400.,500.);
        rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label,"rrv_width_ErfExp"+label, 50.,10,100.);
        rrv_residue_ErfExp = RooRealVar("rrv_residue_ErfExp"+label,"rrv_residue_ErfExp"+label,0.,0.,1.);
        model_pdf = RooGenericPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, "(TMath::Exp(%s*%s) + %s)*(1.+TMath::Erf((%s-%s)/%s))/2. "%(rrv_c_ErfExp.GetName(),rrv_x.GetName(), rrv_residue_ErfExp.GetName(), rrv_x.GetName(),rrv_offset_ErfExp.GetName(), rrv_width_ErfExp.GetName()), RooArgList(rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp,rrv_residue_ErfExp) )

    ## different initial values -> for limit_variable
    if in_model_name == "ErfExp_v3" : #different init-value and range
        print "########### Erf*Exp for limit_variable fit  ############"
        rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label,"rrv_c_ErfExp"+label,-0.005,-0.1,0.);
        rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label,"rrv_offset_ErfExp"+label,450.,400,500.);
        rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label,"rrv_width_ErfExp"+label, 50.,10,100.);
        rrv_residue_ErfExp = RooRealVar("rrv_residue_ErfExp"+label,"rrv_residue_ErfExp"+label,0.,0.,1.);
        rrv_high_ErfExp    = RooRealVar("rrv_high_ErfExp"+label,"rrv_high_ErfExp"+label,1.,0.,400);
        rrv_high_ErfExp.setConstant(kTRUE);
        model_pdf = RooGenericPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, "(TMath::Exp(%s*%s) + %s)* TMath::Power( ((1+TMath::Erf((%s-%s)/%s))/2.), %s )"%(rrv_c_ErfExp.GetName(),rrv_x.GetName(), rrv_residue_ErfExp.GetName(),rrv_x.GetName(),rrv_offset_ErfExp.GetName(), rrv_width_ErfExp.GetName(), rrv_high_ErfExp.GetName()), RooArgList(rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_high_ErfExp,rrv_width_ErfExp,rrv_residue_ErfExp) )

    ## Exp+Gaus or obs0_variable spectrum
    if in_model_name == "ExpGaus":
        print "########### Exp + Gaus for obs0_variable  fit  ############"
        rrv_c_Exp       = RooRealVar("rrv_c_Exp"+label,"rrv_c_Exp"+label,0.05,-0.2,0.2);
        exp             = ROOT.RooExponential("exp"+label,"exp"+label,rrv_x,rrv_c_Exp);

        rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label,"rrv_mean1_gaus"+label,84,78,88);
        rrv_sigma1_gaus = RooRealVar("rrv_smgma1_gaus"+label,"rrv_sigma1_gaus"+label,7,4,10);
        rrv_high        = RooRealVar("rrv_high"+label,"rrv_high"+label,0.5,0.,1.);
        gaus            = RooGaussian("gaus"+label,"gaus"+label, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

        model_pdf       = RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,RooArgList(exp,gaus),RooArgList(rrv_high))

    ## Erf*Exp + Gaus for obs0_variable spectrum 
    if in_model_name == "ErfExpGaus":
        print "########### Erf*Exp + Gaus for obs0_variable  fit  ############"
        rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label,"rrv_c_ErfExp"+label,-0.05,-0.4,0.);
        rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label,"rrv_offset_ErfExp"+label,100.,10.,300.);
        rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label,"rrv_width_ErfExp"+label,30.,10,100.);

        erfExp = ROOT.RooErfExpPdf("erfExp"+label,"erfExp"+label,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label,"rrv_mean_gaus"+label,82,78,87);
        rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label,"rrv_sigma_gaus"+label,7,4,10);
        gaus = RooGaussian("gaus"+label,"gaus"+label, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

        rrv_high  = RooRealVar("rrv_high"+label,"rrv_high"+label,0.7,0.,1.);
        model_pdf = RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))

    ## Erf*Exp + Gaus for obs0_variable spectrum with offset == mean
    if in_model_name == "ErfExpGaus_sp":
        print "########### Erf*Exp + Gaus for obs0_variable  fit  ############"
        rrv_c_ErfExp     = RooRealVar("rrv_c_ErfExp"+label,"rrv_c_ErfExp"+label,-0.05,-0.2,0.);
        rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label,"rrv_width_ErfExp"+label,30.,10,200.);
        erfExp           = ROOT.RooErfExpPdf("erfExp"+label,"erfExp"+label,rrv_x,rrv_c_ErfExp,rrv_mean1_gaus,rrv_width_ErfExp);

        rrv_mean1_gaus   = RooRealVar("rrv_mean1_gaus"+label,"rrv_mean1_gaus"+label,84,78,88);
        rrv_sigma1_gaus  = RooRealVar("rrv_sigma1_gaus"+label,"rrv_sigma1_gaus"+label,7,4,10);
        gaus             = RooGaussian("gaus"+label,"gaus"+label, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

        rrv_high  = RooRealVar("rrv_high"+label,"rrv_high"+label,0.5,0.,1.);
        model_pdf = RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))

    ## Erf*Exp+Gaus or obs0_variable spectrum
    if in_model_name == "ErfExpGaus_v0":
        print "########### Erf*Exp + Gaus for obs0_variable  fit  ############"
        rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label,"rrv_c_ErfExp"+label,-0.05,-0.2,0.);
        rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label,"rrv_offset_ErfExp"+label,100.,10.,140.);
        rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label,"rrv_width_ErfExp"+label,30.,10,100.);
        erfExp = ROOT.RooErfExpPdf("erfExp"+label,"erfExp"+label,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        rrv_mean_gaus     = RooRealVar("rrv_mean_gaus"+label,"rrv_mean_gaus"+label,84,78,88);
        rrv_sigma_gaus    = RooRealVar("rrv_sigma_gaus"+label,"rrv_sigma_gaus"+label,7,4,10);
        gaus              = RooGaussian("gaus"+label,"gaus"+label, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

        rrv_high   = RooRealVar("rrv_high"+label,"rrv_high"+label,0.7,0.,1.);
        model_pdf  = RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))

    ## Erf*Exp+Gaus or obs0_variable spectrum
    if in_model_name == "ErfExpGaus_v1":
        print "########### Erf*Exp + Gaus for limit_variable fit  ############"
        rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label,"rrv_c_ErfExp"+label,-0.007,-0.1,0.);
        rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label,"rrv_offset_ErfExp"+label,800.,10.,1400.);
        rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label,"rrv_width_ErfExp"+label,24.,10,150.);
        erfExp             = ROOT.RooErfExpPdf("erfExp"+label,"erfExp"+label,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        rrv_mean_gaus   = RooRealVar("rrv_mean_gaus"+label,"rrv_mean_gaus"+label,700,500,1200);
        rrv_sigma_gaus  = RooRealVar("rrv_sigma_gaus"+label,"rrv_sigma_gaus"+label,150,10,300);
        gaus            = RooGaussian("gaus"+label,"gaus"+label, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

        rrv_high  = RooRealVar("rrv_high"+label,"rrv_high"+label,0.1,0.,1.);
        model_pdf = RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))

    ## Erf*Exp+Gaus or obs0_variable spectrum
    if in_model_name == "ErfExpGaus_sp_v1":
        print "########### Erf*Exp + Gaus for limit_variable fit  ############"
        rrv_c_ErfExp     = RooRealVar("rrv_c_ErfExp"+label,"rrv_c_ErfExp"+label,-0.007,-0.1,0.);
        rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label,"rrv_width_ErfExp"+label,24.,10,150.);
        rrv_mean_gaus    = RooRealVar("rrv_mean_gaus"+label,"rrv_mean_gaus"+label,900,860,1200);
        erfExp           = ROOT.RooErfExpPdf("erfExp"+label,"erfExp"+label,rrv_x,rrv_c_ErfExp,rrv_mean_gaus,rrv_width_ErfExp);

        rrv_sigma_gaus   = RooRealVar("rrv_sigma_gaus"+label,"rrv_sigma_gaus"+label,150,10,300);
        gaus = RooGaussian("gaus"+label,"gaus"+label, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

        rrv_high  = RooRealVar("rrv_high"+label,"rrv_high"+label,0.1,0.,1.);
        model_pdf = RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))

    ## Erf*Exp+Gaus or obs0_variable spectrum
    if in_model_name == "ErfExpGaus_v2":
        print "########### Erf*Exp + Gaus for obs0_variable fit  ############"
        rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label,"rrv_c_ErfExp"+label,-0.05,-10.,0.);
        rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label,"rrv_offset_ErfExp"+label,100.,10.,140.);
        rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label,"rrv_width_ErfExp"+label,30.,10,100.);
        rrv_mean_gaus     = RooRealVar("rrv_mean_gaus"+label,"rrv_mean_gaus"+label,84,78,88);
        rrv_sigma_gaus    = RooRealVar("rrv_sigma_gaus"+label,"rrv_sigma_gaus"+label,7,4,10);
        rrv_high          = RooRealVar("rrv_high"+label,"rrv_high"+label,200.,0.,1000.);
        model_pdf = ROOT.RooErfExp_Gaus_Pdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp,rrv_mean_gaus,rrv_sigma_gaus,rrv_high );

    ## Erf*Exp + 2Gaus  
    if in_model_name == "ErfExp2Gaus":
        print "########### Erf*Exp + 2Gaus for obs0_variable fit  ############"
        rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label,"rrv_c_ErfExp"+label,-0.05,-0.2,0.);
        rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label,"rrv_offset_ErfExp"+label,100.,10.,140.);
        rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label,"rrv_width_ErfExp"+label,30.,10,100.);
        erfExp = ROOT.RooErfExpPdf("erfExp"+label,"erfExp"+label,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        rrv_mean1_gaus   = RooRealVar("rrv_mean1_gaus"+label,"rrv_mean1_gaus"+label,84,78,88);
        rrv_mean2_gaus   = RooRealVar("rrv_mean2_gaus"+label,"rrv_mean2_gaus"+label,180,170,190);
        rrv_sigma1_gaus  = RooRealVar("rrv_sigma1_gaus"+label,"rrv_sigma1_gaus"+label,7,4,10);
        rrv_sigma2_gaus  = RooRealVar("rrv_sigma2_gaus"+label,"rrv_sigma2_gaus"+label,10,7,15);
        gaus1 = RooGaussian("gaus1"+label,"gaus1"+label, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);
        gaus2 = RooGaussian("gaus2"+label,"gaus2"+label, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

        rrv_high1 = RooRealVar("rrv_high1"+label,"rrv_high1"+label,0.6,0.,1.);
        rrv_high2 = RooRealVar("rrv_high2"+label,"rrv_high2"+label,0.4,0.,1.);
        model_pdf =RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,RooArgList(erfExp,gaus1,gaus2),RooArgList(rrv_high1,rrv_high2))

    ## Gaus + Gaus for obs0_variable spectrum
    if in_model_name == "2Gaus":
        print "########### 2Gaus for obs0_variable fit  ############"
        mean1_tmp      = 8.3141e+01; mean1_tmp_err      = 1.63e-01;
        deltamean_tmp  = 6.9129e+00; deltamean_tmp_err  = 1.24e+00;
        sigma1_tmp     = 7.5145e+00; sigma1_tmp_err     = 1.99e-01;
        scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
        frac_tmp       = 6.7125e-01; frac_tmp_err       = 2.09e-02;

        rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label,"rrv_mean1_gaus"+label,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
        rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label,"rrv_sigma1_gaus"+label,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
        gaus1 = RooGaussian("gaus1"+label,"gaus1"+label, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

        rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label,"rrv_deltamean_gaus"+label,deltamean_tmp,deltamean_tmp-deltamean_tmp_err*4 ,deltamean_tmp+deltamean_tmp_err*4);
        rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
        rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label,"rrv_scalesigma_gaus"+label,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*8, scalesigma_tmp+scalesigma_tmp_err*8);
        rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
        gaus2 = RooGaussian("gaus2"+label,"gaus2"+label, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

        rrv_frac  = RooRealVar("rrv_frac"+label,"rrv_frac"+label,frac_tmp, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);
        model_pdf = RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac),1)

    ## 2Gaus+2Gaus for VV obs0_variable spectrum -> WZ and WW
    if in_model_name == "2_2Gaus":

        print "########### 2Gaus +2Gaus for obs0_variable fit  ############"
        mean1_tmp      = 8.3141e+01; mean1_tmp_err      = 1.63e-01;
        deltamean_tmp  = 6.9129e+00; deltamean_tmp_err  = 1.24e+00;
        sigma1_tmp     = 7.5145e+00; sigma1_tmp_err     = 1.99e-01;
        scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
        frac_tmp       = 6.7125e-01; frac_tmp_err       = 2.09e-02;

        rrv_shift = RooRealVar("rrv_shift"+label,"rrv_shift"+label,10.8026) # Z mass: 91.1876; shift=91.1876-80.385=10.8026

        rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label,"rrv_mean1_gaus"+label,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
        rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label,"rrv_sigma1_gaus"+label,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
        gaus1 = RooGaussian("gaus1"+label,"gaus1"+label, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

        rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label,"rrv_deltamean_gaus"+label,0.,-8,10);
        rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
        rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label,"rrv_scalesigma_gaus"+label,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*4, scalesigma_tmp+scalesigma_tmp_err*4);
        rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
        gaus2 = RooGaussian("gaus2"+label,"gaus2"+label, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

        rrv_frac1 = RooRealVar("rrv_frac1"+label,"rrv_frac1"+label,frac_tmp, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);
        gausguas_1 =RooAddPdf("gausguas_1"+label+mass_spectrum,"gausguas_1"+label+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac1),1)

        rrv_mean3_gaus = RooFormulaVar("rrv_mean3_gaus"+label,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_shift));
        rrv_mean4_gaus = RooFormulaVar("rrv_mean4_gaus"+label,"@0+@1",RooArgList(rrv_mean2_gaus, rrv_shift));
        gaus3 = RooGaussian("gaus3"+label,"gaus3"+label, rrv_x,rrv_mean3_gaus,rrv_sigma1_gaus);
        gaus4 = RooGaussian("gaus4"+label,"gaus4"+label, rrv_x,rrv_mean4_gaus,rrv_sigma2_gaus);
        gausguas_2 = RooAddPdf("gausguas_2"+label+mass_spectrum,"gausguas_2"+label+mass_spectrum,RooArgList(gaus3,gaus4),RooArgList(rrv_frac1),1)

        rrv_frac  = RooRealVar("rrv_frac"+label,"rrv_frac"+label,0.74)#,0.5,1.0)
        model_pdf = RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,RooArgList(gausguas_1,gausguas_2),RooArgList(rrv_frac),1)

    ## Erf*Exp + 2Gaus for obs0_variable spectrum
    if in_model_name == "2Gaus_ErfExp":

        print "########### 2Gaus + Erf*Exp for obs0_variable fit  ############"
        mean1_tmp      = 8.3141e+01; mean1_tmp_err      = 1.63e-01;
        deltamean_tmp  = 6.9129e+00; deltamean_tmp_err  = 1.24e+00;
        sigma1_tmp     = 7.5145e+00; sigma1_tmp_err     = 1.99e-01;
        scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
        frac_tmp       = 6.7125e-01; frac_tmp_err       = 2.09e-02;

        rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label,"rrv_mean1_gaus"+label,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
        rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label,"rrv_sigma1_gaus"+label,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
        gaus1 = RooGaussian("gaus1"+label,"gaus1"+label, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

        rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label,"rrv_deltamean_gaus"+label,deltamean_tmp)#, deltamean_tmp, deltamean_tmp);
        rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
        rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label,"rrv_scalesigma_gaus"+label,scalesigma_tmp)#, scalesigma_tmp, scalesigma_tmp);
        rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
        gaus2 = RooGaussian("gaus2"+label,"gaus2"+label, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

        rrv_frac_2gaus = RooRealVar("rrv_frac_2gaus"+label,"rrv_frac_2gaus"+label,frac_tmp);#, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);

        c0_tmp     = -2.9893e-02 ; c0_tmp_err     = 6.83e-03;
        offset_tmp = 7.9350e+01  ; offset_tmp_err = 9.35e+00;
        width_tmp  = 3.3083e+01  ; width_tmp_err  = 2.97e+00;

        rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label,"rrv_c_ErfExp"+label,c0_tmp, c0_tmp-4e-2, c0_tmp+4e-2 );
        rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label,"rrv_offset_ErfExp"+label, offset_tmp)#, offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
        rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label,"rrv_width_ErfExp"+label, width_tmp, width_tmp-10, width_tmp+10);
        erfexp = ROOT.RooErfExpPdf("erfexp"+label+mass_spectrum,"erfexp"+label+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        rrv_frac = RooRealVar("rrv_frac"+label,"rrv_frac"+label, 0.5,0.,1.);
        model_pdf =RooAddPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,RooArgList(erfexp, gaus1,gaus2),RooArgList(rrv_frac, rrv_frac_2gaus),1)


    ## Erf*Exp+Voig+Gaus for obs0_variable spectrum 
    if in_model_name == "ErfExpVoigGaus":
        print "########### Erf*Exp + Voig + Gaus for obs0_variable fit  ############"
        rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label,"rrv_c_ErfExp"+label,-0.1,-10.,0.);
        rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label,"rrv_offset_ErfExp"+label,100.,10.,140.);
        rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label,"rrv_width_ErfExp"+label,30.,10,100.);
        rrv_mean_voig     = RooRealVar("rrv_mean_voig"+label,"rrv_mean_voig"+label,84,78,88);
        rrv_width_voig    = RooRealVar("rrv_width_voig"+label,"rrv_width_voig"+label,7,1,20);
        rrv_sigma_voig    = RooRealVar("rrv_sigma_voig"+label,"rrv_sigma_voig"+label,5,1,100);
        rrv_high1         = RooRealVar("rrv_high1"+label,"rrv_high1"+label,1,0.,200.);
        rrv_mean_gaus     = RooRealVar("rrv_mean_gaus"+label,"rrv_mean_gaus"+label,174)#,160,187);
        rrv_sigma_gaus    = RooRealVar("rrv_sigma_gaus"+label,"rrv_sigma_gaus"+label,20)#,0.1,100);
        rrv_high2 = RooRealVar("rrv_high2"+label,"rrv_high2"+label,0.)#,0.,0.);
        model_pdf = ROOT.RooErfExp_Voig_Gaus_Pdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp,rrv_mean_voig,rrv_width_voig,rrv_sigma_voig,rrv_high1,rrv_mean_gaus,rrv_sigma_gaus,rrv_high2 );

    ## User1 function 
    if in_model_name == "User1":
        print "########### User 1 Pdf  for limit_variable fit ############"
        rrv_p0 = RooRealVar("rrv_p0_User1"+label,"rrv_p0_User1"+label, 30, 10, 90);
        rrv_p1 = RooRealVar("rrv_p1_User1"+label,"rrv_p1_User1"+label, -4, -9, -2);
        model_pdf=RooUser1Pdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_p0,rrv_p1);

    ## QCD pdf  
    if in_model_name == "QCD":
        print "########### QCD Pdf  for limit_variable fit ############"
        rrv_p0 = RooRealVar("rrv_p0_QCD"+label,"rrv_p0_QCD"+label, 0,-200,200);
        rrv_p1 = RooRealVar("rrv_p1_QCD"+label,"rrv_p1_QCD"+label, 0,-200,200);
        rrv_p2 = RooRealVar("rrv_p2_QCD"+label,"rrv_p2_QCD"+label, 0,-200,200);
        model_pdf = RooQCDPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_p0,rrv_p1,rrv_p2);

    if in_model_name == "QCD_v2":#can replace exp
        print "########### QCD Pdf  for limit_variable fit ############"
        rrv_p0 = RooRealVar("rrv_p0_QCD"+label,"rrv_p0_QCD"+label, -15,-50,0);
        rrv_p1 = RooRealVar("rrv_p1_QCD"+label,"rrv_p1_QCD"+label, 20,0,250);
        rrv_p2 = RooRealVar("rrv_p2_QCD"+label,"rrv_p2_QCD"+label,0,-20,20);
        model_pdf = RooQCDPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_p0,rrv_p1,rrv_p2);

    ## For limit_variable fit -> Pow function can replace exp
    if in_model_name == "Pow" or in_model_name == "Pow_sr" :
        print "########### Pow Pdf  for limit_variable fit ############"
        rrv_c = RooRealVar("rrv_c_Pow"+label,"rrv_c_Pow"+label, -5, -20, 0);
        model_pdf = RooPowPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x, rrv_c );

    ## For limit_variable fit -> Pow function can replace exp
    if in_model_name == "Pow2":
        print "########### Pow2 Pdf  for limit_variable fit ############"
        rrv_c0 = RooRealVar("rrv_c0_Pow2"+label,"rrv_c0_Pow2"+label, 5, 0, 20);
        rrv_c1 = RooRealVar("rrv_c1_Pow2"+label,"rrv_c1_Pow2"+label, 0, -5 , 5);
        model_pdf = RooPow2Pdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x, rrv_c0, rrv_c1 );

    ## For limit_variable fit ->Erf*Pow can replace Erf*Exp
    if in_model_name == "ErfPow_v1":
        print "########### Erf*Pow Pdf  for limit_variable fit ############"
        rrv_c      = RooRealVar("rrv_c_ErfPow"+label,"rrv_c_ErfPow"+label, -5,-10,0);
        rrv_offset = RooRealVar("rrv_offset_ErfPow"+label,"rrv_offset_ErfPow"+label, 450,350,550);
        rrv_width  = RooRealVar("rrv_width_ErfPow"+label,"rrv_width_ErfPow"+label,50,20,90);
        model_pdf  = RooErfPowPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_c,rrv_offset,rrv_width);

    ## For limit_variable fit ->Erf*Pow can replace Erf*Exp -> in the sideband
    if in_model_name == "ErfPow2_v1":
        print "########### Erf*Pow2 Pdf  for limit_variable fit ############"
        rrv_c0 = RooRealVar("rrv_c0_ErfPow2"+label,"rrv_c0_ErfPow2"+label,14,1,30);
        rrv_c1 = RooRealVar("rrv_c1_ErfPow2"+label,"rrv_c1_ErfPow2"+label, 5,-5,10);
        rrv_offset = RooRealVar("rrv_offset_ErfPow2"+label,"rrv_offset_ErfPow2"+label, 450,400,520);
        rrv_width  = RooRealVar("rrv_width_ErfPow2"+label,"rrv_width_ErfPow2"+label,30,10,80);
        model_pdf  = RooErfPow2Pdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);

    ## For limit_variable fit ->Erf*Pow can replace Erf*Exp for sr
    if in_model_name == "ErfPow2_v1_sr":
        print "########### Erf*Pow2 Pdf  for limit_variable fit in the SR  ############"
        rrv_c0 = RooRealVar("rrv_c0_ErfPow2"+label,"rrv_c0_ErfPow2"+label, 4,2, 8);
        rrv_c1 = RooRealVar("rrv_c1_ErfPow2"+label,"rrv_c1_ErfPow2"+label, -0.5,-2,0);
        rrv_offset = RooRealVar("rrv_offset_ErfPow2"+label,"rrv_offset_ErfPow2"+label, 490,440,520);
        rrv_width  = RooRealVar("rrv_width_ErfPow2"+label,"rrv_width_ErfPow2"+label,50,30,80);
        model_pdf = RooErfPow2Pdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);

    ## For limit_variable fit ->Erf*Pow*Exp can replace Erf*Exp 
    if in_model_name == "ErfPowExp_v1":
        print "########### Erf*Pow*Exp Pdf  for limit_variable fit   ############"
        rrv_c0 = RooRealVar("rrv_c0_ErfPowExp"+label,"rrv_c0_ErfPowExp"+label,11,5,20);
        rrv_c1 = RooRealVar("rrv_c1_ErfPowExp"+label,"rrv_c1_ErfPowExp"+label, 0,-2,2);
        rrv_offset = RooRealVar("rrv_offset_ErfPowExp"+label,"rrv_offset_ErfPowExp"+label, 470,420,520);
        rrv_width  = RooRealVar("rrv_width_ErfPowExp"+label,"rrv_width_ErfPowExp"+label,40,30,50);
        model_pdf  = RooErfPowExpPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);

    ## For limit_variable fit ->Erf*Pow*Exp can replace Erf*Exp 
    if in_model_name == "ErfPowExp_v1_sr":
        print "########### Erf*Pow*Exp Pdf for limit_variable fit in SR  ############"
        rrv_c0 = RooRealVar("rrv_c0_ErfPowExp"+label,"rrv_c0_ErfPowExp"+label,6,2,15);
        rrv_c1 = RooRealVar("rrv_c1_ErfPowExp"+label,"rrv_c1_ErfPowExp"+label, -1,-3,2);
        rrv_offset = RooRealVar("rrv_offset_ErfPowExp"+label,"rrv_offset_ErfPowExp"+label, 490,440,520);
        rrv_width  = RooRealVar("rrv_width_ErfPowExp"+label,"rrv_width_ErfPowExp"+label,50,30,70);
        model_pdf=RooErfPowExpPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);

    ## For limit_variable fit ->Erf*Pow*Exp can replace Erf*Exp 
    if in_model_name == "ErfPowExp_v1_0":#difference inital value
        print "########### Erf*Pow*Exp Pdf for limit_variable fit in SR  ############"
        rrv_c0 = RooRealVar("rrv_c0_ErfPowExp"+label,"rrv_c0_ErfPowExp"+label,20,15,40);
        rrv_c1 = RooRealVar("rrv_c1_ErfPowExp"+label,"rrv_c1_ErfPowExp"+label, 1.6,0.5,5);
        rrv_offset = RooRealVar("rrv_offset_ErfPowExp"+label,"rrv_offset_ErfPowExp"+label, 470,420,520);
        rrv_width  = RooRealVar("rrv_width_ErfPowExp"+label,"rrv_width_ErfPowExp"+label,47,30,60);
        model_pdf  = RooErfPowExpPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);

    ## Keys 
    if in_model_name == "Keys":
        print "########### Keys PDF  ############"
        rdataset = workspace.data("rdataset4fit"+label+"_limit_variable");
        rdataset.Print();
        model_pdf = RooKeysPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, rrv_x, rdataset);


    ## Hist 
    if in_model_name == "Hist":
        print "########### Hist PDF  ############"
        rdataset = workspace.data("rdataset4fit"+label+"_limit_variable");
        rdataset.Print();
        rdatahist= rdataset.binnedClone("rdatahist4fit"+label+"_limit_variable","rdatahist4fit"+label+"_limit_variable");
        model_pdf = RooHistPdf("model_pdf"+label+mass_spectrum,"model_pdf"+label+mass_spectrum, RooArgSet(rrv_x), rdatahist);


    ## return the pdf
    getattr(workspace,"import")(model_pdf)
    return workspace.pdf("model_pdf"+label+mass_spectrum)

### in order to make the plot_legend
def legend4Plot(plot, title, left=1, isFill=1, x_offset_low=0., y_offset_low=0., x_offset_high =0., y_offset_high =0., TwoCoulum =1.):
    print "############### draw the plot_legend ########################"
    if left==-1:
        theLeg = TLegend(0.65+x_offset_low, 0.58+y_offset_low, 0.93+x_offset_low, 0.87+y_offset_low, "", "NDC");
        theLeg.SetName("theLegend");
        theLeg.SetLineColor(0);
        theLeg.SetTextFont(42);
        theLeg.SetTextSize(.04);
    else:
        theLeg = TLegend(0.41+x_offset_low, 0.61+y_offset_low, 0.76+x_offset_high, 0.93+y_offset_high, "", "NDC");            
        theLeg.SetName("theLegend");
        if TwoCoulum :
            theLeg.SetNColumns(2);

    theLeg.SetFillColor(0);
    theLeg.SetFillStyle(0);
    theLeg.SetBorderSize(0);
    theLeg.SetLineColor(0);
    theLeg.SetLineWidth(0);
    theLeg.SetLineStyle(0);
    theLeg.SetTextSize(0.040);
    theLeg.SetTextFont(42);

    
    nLeg_total=0;
    nLeg_data=0; nLeg_uncertainty=0;
    leg_list=[];
    objName_before = "";
    for obj in range(int(plot.numItems()) ):
        objName = plot.nameOf(obj);
        if objName == "errorband" : objName = "Uncertainty";
        print objName;
        if not (  plot.getInvisible(objName) or TString(objName).Contains("invisi") or TString(objName).Contains("line") or objName ==objName_before ):
            leg_obj = plot.getObject(obj); 
            leg_name= objName;
            leg_opt = plot.getDrawOptions(objName).Data()

            if leg_opt=="P": leg_opt="PE"

            if   TString(leg_name).Contains("Uncertainty"):
                leg_opt="F"; 
                nLeg_uncertainty=nLeg_total;
            elif TString(leg_name).Contains("sigma"): 
                leg_opt="F"; 
            elif TString(leg_name).Data()=="data" : 
                leg_name="CMS Data "+title; 
                nLeg_data=nLeg_total;
            elif TString(leg_name).Data()=="WJets" :
                leg_name="W+jets"; 
            elif TString(leg_name).Data()=="SingleT" : 
                leg_name="Single t"; 
            elif TString(leg_name).Data()=="TTbar" :
                leg_name="t#bar{t}"; 
            elif TString(leg_name).Data()=="VV" : 
                leg_name="WW/WZ"; 
            leg_list.append( [leg_obj, leg_name, leg_opt]  );
            nLeg_total+=1;
            objName_before=objName;

    #need to draw data and uncertainty at first;
    theLeg.AddEntry(leg_list[nLeg_data][0],leg_list[nLeg_data][1],leg_list[nLeg_data][2]);
    theLeg.AddEntry(leg_list[nLeg_uncertainty][0],leg_list[nLeg_uncertainty][1],leg_list[nLeg_uncertainty][2]);
    for iter in range(nLeg_total):
        if not (iter==nLeg_data or iter==nLeg_uncertainty):
            theLeg.AddEntry(leg_list[iter][0],leg_list[iter][1],leg_list[iter][2]);

    return theLeg;


### fix a given model taking the label, and the region --> for extended pdf --> all the parameter of the pdf + normalization
def fix_Model(self, label, limit_variable_region="_signalregion",mass_spectrum="_limit_variable"):
    print "########### Fixing an Extended Pdf for limit_variable  ############"        
    rdataset = self.workspace4fit_.data("rdataset%s%s_%s%s"%(label,limit_variable_region,self.categoryLabel,mass_spectrum))
    model = self.get_limit_variable_Model(label,limit_variable_region);
    rdataset.Print();
    model.Print();
    parameters = model.getParameters(rdataset);
    par=parameters.createIterator(); par.Reset();
    param=par.Next()
    while (param):
        param.setConstant(kTRUE);
        param=par.Next()

### fix a pdf in a different way --> for RooAbsPdf 
def fix_Pdf(model_pdf,argset_notparameter):
    print "########### Fixing a RooAbsPdf for limit_variable or obs0_variable  ############"        
    parameters = model_pdf.getParameters(argset_notparameter);
    par=parameters.createIterator(); par.Reset();
    param=par.Next()
    while (param):
        param.setConstant(kTRUE);
        param.Print();
        param=par.Next()

### print the parameters of a given pdf --> only non constant ones
def ShowParam_Pdf(model_pdf,argset_notparameter):
    print "########### Show Parameters of a input model  ############"        
    model_pdf.Print()
    parameters = model_pdf.getParameters(argset_notparameter);
    par = parameters.createIterator(); par.Reset();
    param = par.Next()
    while (param):
        if not param.isConstant():
            param.Print();
            if (param.getVal()-param.getMin())< param.getError()*1 or (param.getMax()- param.getVal())< param.getError()*1:
                param.Print();
        param=par.Next()


    ### change a dataset to a histpdf roofit object
def change_dataset_to_histpdf(self,x,dataset):
    print "######## change the dataset into a histpdf  ########"        
    datahist = dataset.binnedClone(dataset.GetName()+"_binnedClone",dataset.GetName()+"_binnedClone")
    histpdf = RooHistPdf(dataset.GetName()+"_histpdf",dataset.GetName()+"_histpdf",RooArgSet(x),datahist)
    dataset.Print();
    histpdf.Print();
    getattr(self.workspace4fit_,"import")(histpdf)

### change from a dataset to a histogramm of Roofit
def change_dataset_to_histogram(self, x,dataset,label=""):
    print "######## change the dataset into a histogramm for obs0_variable distribution ########"        
    datahist=dataset.binnedClone(dataset.GetName()+"_binnedClone",dataset.GetName()+"_binnedClone")
    nbin=int( (x.getMax()-x.getMin())/self.obs0_variable_BinWidth);
    if label=="":
        return datahist.createHistogram("histo_%s"%(dataset.GetName()),x, RooFit.Binning( nbin ,x.getMin(),x.getMax()));
    else:
        return datahist.createHistogram("histo_"+label,x, RooFit.Binning( nbin,x.getMin(),x.getMax()));
