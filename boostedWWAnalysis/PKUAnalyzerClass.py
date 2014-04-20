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

from PKUAnalyzerTools import * 

###############################
## doFit Class Implemetation ##
###############################

class doFit_wj_and_wlvj:

    def __init__(self, in_analyzer_config):
        print "############################################################################"
        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9) ;
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9) ;

        self.analyzer_config=in_analyzer_config;
        print self.analyzer_config;

        ######################
        # plot init
        ######################
        #plot style
        setTDRStyle();
        # color palet for plots
        self.color_palet={ 
                'Uncertainty' : kBlack,
                'Other_Backgrounds' : kBlue
                }
        self.sig_bkg_files=self.analyzer_config["sig_bkg_files"];
        for iter in range(2, 3+self.sig_bkg_files[0]+self.sig_bkg_files[1]):
            self.color_palet[ self.sig_bkg_files[iter][0] ] =self.sig_bkg_files[iter][1] 
        print self.color_palet;
        #plot legend
        self.plot_legend = TLegend(); 

        ######################
        # analysis init
        ######################
        self.workspace4limit_ = RooWorkspace("workspace4limit_","workspace4limit_"); ## for unbin limit calculation.
        self.workspace4fit_ = RooWorkspace("workspace4fit_","workspace4fit_"); ## temporary workspace for analysis

        # category
        self.categoryID    = self.analyzer_config["categoryID"];
        self.categoryLabel = self.analyzer_config["categoryLabel"][0];#lable used for naming;
        self.categoryTitle = self.analyzer_config["categoryLabel"][1];#title used for plot
        
        # limit variable, obs variable
        #narrow the obs0_variable_BinWidth and limit_variable_BinWidth by a factor of 5. Because Higgs-Combination-Tools will generate a binned sample, so need the bin width narrow. So, as a easy selution, we will increase the bin-width by a factor of 5 when ploting m_j m_WW
        self.BinWidth_narrow_factor=1.;## mlvj 100/5=20 mj 5/5=1

        self.obs0_variable_BinWidth=self.analyzer_config["obs0_variable_BinWidth"];
        self.obs0_variable_BinWidth=self.obs0_variable_BinWidth/self.BinWidth_narrow_factor;
        mj_min=self.analyzer_config["obs0_variable_full_range_min"];
        mj_max=self.analyzer_config["obs0_variable_full_range_max"];
        nbins_mj=int( (mj_max - mj_min)/self.obs0_variable_BinWidth );
        ## correct mj max
        mj_max=mj_min+nbins_mj*self.obs0_variable_BinWidth;
        ## define jet mass variable
        rrv_mass_j = RooRealVar("rrv_mass_j","Jet Mass",(mj_min+mj_max)/2.,mj_min,mj_max,"GeV");
        rrv_mass_j.setBins(nbins_mj);
        rrv_mass_j.setRange("lowersideband",self.analyzer_config["obs0_variable_lowersideband_range_min"],self.analyzer_config["obs0_variable_lowersideband_range_max"]);
        rrv_mass_j.setRange("signalregion",self.analyzer_config["obs0_variable_signalregion_range_min"],self.analyzer_config["obs0_variable_signalregion_range_max"]);
        rrv_mass_j.setRange("uppersideband",self.analyzer_config["obs0_variable_uppersideband_range_min"],self.analyzer_config["obs0_variable_uppersideband_range_max"]);
        getattr(self.workspace4fit_,"import")(rrv_mass_j);
        rrv_mass_j.Print();

        self.limit_variable_BinWidth=self.analyzer_config["limit_variable_BinWidth"];
        self.limit_variable_BinWidth=self.limit_variable_BinWidth/self.BinWidth_narrow_factor;
        mlvj_min=self.analyzer_config["limit_variable_full_range_min"];
        mlvj_max=self.analyzer_config["limit_variable_full_range_max"];
        nbins_mlvj=int((mlvj_max-mlvj_min)/self.limit_variable_BinWidth);
        ## correct mlvj max 
        mlvj_max=mlvj_min+nbins_mlvj*self.limit_variable_BinWidth;
        ## define invariant mass WW variable
        rrv_mass_lvj= RooRealVar("rrv_mass_lvj","M_{WW}",(mlvj_min+mlvj_max)/2.,mlvj_min,mlvj_max,"GeV");
        rrv_mass_lvj.setBins(nbins_mlvj);
        rrv_mass_lvj.setRange("signalregion", self.analyzer_config["limit_variable_signalregion_range_min"], self.analyzer_config["limit_variable_signalregion_range_max"] );# for cut&count limit calculation 
        getattr(self.workspace4fit_,"import")(rrv_mass_lvj);
        rrv_mass_lvj.Print();

        ## set the model used for the background parametrization
        self.MODEL_4_mlvj       = self.analyzer_config["fit_model"];
        self.MODEL_4_mlvj_alter = self.analyzer_config["fit_model_alter"];


        ######################
        # input and output init
        ######################
        self.additioninformation=self.analyzer_config["additioninformation"];
        self.signal_scale= self.analyzer_config["signal_scale"];

        self.nsig=self.sig_bkg_files[0];
        self.nbkg=self.sig_bkg_files[1];
        self.sig_list=[];
        self.bkg_list=[];
        #self.allsignals=self.sig_bkg_files[3][0];#first signal as prime signal
        self.allsignals="";#first signal as prime signal
        for iter in range(3, 3+self.nsig):
            self.sig_list.append( self.sig_bkg_files[iter] );
            self.allsignals+=self.sig_bkg_files[iter][0]; 
            #self.sig_list.append( (self.sig_bkg_files[iter][0], self.sig_bkg_files[iter][2]) );
        for iter in range( 3+self.nsig, 3+self.nsig+self.nbkg):
            self.bkg_list.append( self.sig_bkg_files[iter] );
            #self.bkg_list.append( (self.sig_bkg_files[iter][0], self.sig_bkg_files[iter][2]) );
        print self.sig_list;
        print self.bkg_list;

        #result files: The event number, parameters and error write into a txt file.
        #The workspace4limit write into a root file
        self.rlt_DIR="cards_%s_%s/"%(self.additioninformation, self.categoryLabel);
        if not os.path.isdir(self.rlt_DIR):
            os.system("mkdir %s"%(self.rlt_DIR));

        ## extra text file
        self.file_rlt_txt = self.rlt_DIR+"other_wwlvj_%s.txt"%(self.allsignals)
        ## workspace for limit
        self.file_rlt_root = self.rlt_DIR+"wwlvj_%s_workspace.root"%(self.allsignals)
        ## datacard for cut&count and ubninned limit
        self.file_datacard_unbin = self.rlt_DIR+"wwlvj_%s_unbin.txt"%(self.allsignals)
        self.file_datacard_counting = self.rlt_DIR+"wwlvj_%s_counting.txt"%(self.allsignals)

        self.file_out=open(self.file_rlt_txt,"w");
        self.file_out.write("Welcome:\n");
        self.file_out.close()
        self.file_out=open(self.file_rlt_txt,"a+");

        # shape parameter uncertainty for unbin limit setting
        self.FloatingParams=RooArgList("floatpara_list");

    def GetLumi(self):
        if TString(self.categoryLabel).Contains("el"): return 19.7;
        elif TString(self.categoryLabel).Contains("mu"): return 19.7;
        elif TString(self.categoryLabel).Contains("em"): return 19.7;

    def analysis(self):
        print "################### fit_AllSamples_Mj_and_Mlvj #####################"
        self.get_data();
        self.fit_Signals()
        self.fit_Backgrounds()
        self.prepare_limit("SimplyMethod")
        self.read_workspace(1)
 

    def fit_Signals(self):
        print "############### fit_Signals #####################"
        for iter_sig in range(self.nsig):
            print self.sig_list[iter_sig];
            self.fit_SingleChannel(self.sig_list[iter_sig][2], "_%s"%self.sig_list[iter_sig][0], self.sig_list[iter_sig][3]);

    def fit_Backgrounds(self):
        print "############### fit_Backgrounds #####################"
        for iter_bkg in range(self.nbkg):
            print self.bkg_list[iter_bkg]; 
            self.fit_SingleChannel(self.bkg_list[iter_bkg][2], "_%s"%self.bkg_list[iter_bkg][0], self.bkg_list[iter_bkg][3]);

 
    #### Define the steps to fit signal distribution in the mj and mlvj spectra
    def fit_SingleChannel(self, file, label, fit_config,):
        print "############# fit_SingleChannel: %s, %s, %s #################"%(file, label, fit_config)
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(file, label, self.analyzer_config["limit_variable"], self.analyzer_config["obs0_variable"])

        self.fit_limit_variable_SingleChannel(label, "_signalregion", fit_config, 1, 0, 1);
        self.fit_limit_variable_SingleChannel(label, "_lowersideband", fit_config, 0, 0, 1);
        #self.fit_obs_variable_SingleChannel(self.file_SingleT_mc,"_SingleT","ExpGaus");


        #model_narrow="DoubleCB_v1",model_width="BWDoubleCB"  
        #if (TString(self.file_signal).Contains("BulkG_WW_inclusive_M1000_W150") or TString(self.file_signal).Contains("BulkG_WW_inclusive_M1000_W50") or
        #        TString(self.file_signal).Contains("BulkG_WW_inclusive_M1000_W300") or TString(self.file_signal).Contains("BulkG_WW_inclusive_M1500_W225") or
        #        TString(self.file_signal).Contains("BulkG_WW_inclusive_M1500_W450") or TString(self.file_signal).Contains("BulkG_WW_inclusive_M1500_W75") or
        #        TString(self.file_signal).Contains("BulkG_WW_inclusive_M2100_W105") or TString(self.file_signal).Contains("BulkG_WW_inclusive_M2100_W315") or
        #        TString(self.file_signal).Contains("BulkG_WW_inclusive_M2100_W450")):
        #    self.fit_limit_variable_SingleChannel(self.file_signal,"_%s"%(self.allsignals),"_signalregion",model_width, 0, 0, 0, 0);            
        #else:
        #    self.fit_limit_variable_SingleChannel(self.file_signal,"_%s"%(self.allsignals),"_signalregion",model_narrow, 0, 0, 0, 0);

    ### Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
    #def fit_limit_variable_SingleChannel(self,in_file_name, label, in_range, mlvj_model, deco=0, show_constant_parameter=0, logy=0, ismc=0):
    def fit_limit_variable_SingleChannel(self, label, in_range, fit_config, show_constant_parameter=0, logy=0, ismc=0):

        print "############### Fit mlvj single MC sample ",label,"  ",in_range,"  ",fit_config," ##################"
        ## import variable and dataset
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rdataset = self.workspace4fit_.data("rdataset4fit"+label+in_range+"_"+self.categoryLabel+"_mlvj");
        constrainslist =[];

        ## make the extended pdf model
        if fit_config[1]=="Kyes": number_rdata=rdataset.sumEntries();
        else: number_rdata=500;
        model = self.make_Model(label+in_range,fit_config,"_mlvj", number_rdata);

        ## make the fit
        model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE) );
        rfresult = model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult.Print();

        ## set the name of the result of the fit and put it in the workspace   
        rfresult.SetName("rfresult"+label+in_range+"_"+self.categoryLabel+"_mlvj")
        getattr(self.workspace4fit_,"import")(rfresult)

        ## plot the result
        fitting_model=fit_config[1];
        mplot = rrv_mass_lvj.frame(RooFit.Title("M_{lvj"+in_range+"} fitted by "+fitting_model), RooFit.Bins(int(rrv_mass_lvj.getBins()/self.BinWidth_narrow_factor)));
        rdataset.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        ## plot the error band but don't store the canvas (only plotted without -b option
        draw_error_band_extendPdf(rdataset, model, rfresult,mplot,2,"L")
        rdataset.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        model.plotOn( mplot )#, RooFit.VLines()); in order to have the right pull 

        ## get the pull 
        mplot_pull      = self.get_pull(rrv_mass_lvj,mplot);
        parameters_list = model.getParameters(rdataset);
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);

        self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s/limit_variable/"%(self.additioninformation, self.categoryLabel), label, "m_lvj"+in_range+fitting_model, show_constant_parameter, logy);


        ## if the shape parameters has to be decorrelated
        #if deco :
        if fit_config[0] :
            print "################### Decorrelated mlvj single mc shape ################"
            model_pdf = self.workspace4fit_.pdf("model_pdf%s%s_%s_mlvj"%(label,in_range,self.categoryLabel)); ## take the pdf from the workspace
            model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) );
            rfresult_pdf = model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));
            rfresult_pdf.Print();

            ## temp workspace for the pdf diagonalizer
            wsfit_tmp = RooWorkspace("wsfit_tmp"+label+in_range+"_"+self.categoryLabel+"_mlvj");
            Deco      = PdfDiagonalizer("Deco"+label+in_range+"_"+self.categoryLabel+"_"+self.wtagger_label+"_mlvj",wsfit_tmp,rfresult_pdf); ## in order to have a good name 
            print "##################### diagonalize ";
            model_pdf_deco = Deco.diagonalize(model_pdf); ## diagonalize            
            print "##################### workspace for decorrelation ";
            wsfit_tmp.Print("v");
            print "##################### original  parameters ";
            model_pdf.getParameters(rdataset).Print("v");
            print "##################### original  decorrelated parameters ";
            model_pdf_deco.getParameters(rdataset).Print("v");
            print "##################### original  pdf ";
            model_pdf.Print();
            print "##################### decorrelated pdf ";
            model_pdf_deco.Print();

            ## import in the workspace and print the diagonalizerd pdf
            getattr(self.workspace4fit_,"import")(model_pdf_deco);

            ### define a frame for TTbar or other plots
            mplot_deco = rrv_mass_lvj.frame( RooFit.Bins(int(rrv_mass_lvj.getBins()/self.BinWidth_narrow_factor)));

            rdataset.plotOn(mplot_deco, RooFit.Name("Data"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
            model_pdf_deco.plotOn(mplot_deco,RooFit.Name(label),RooFit.LineColor(kBlack));

            mplot_deco.GetYaxis().SetRangeUser(1e-2,mplot_deco.GetMaximum()*1.2);

            rrv_number_dataset=RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset.sumEntries());
            rrv_number_dataset.setError(0.)
            draw_error_band(rdataset, model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,self.color_palet["Uncertainty"],"F"); ## don't store the number in the workspace

            self.plot_legend = self.legend4Plot(mplot_deco,0); ## add the plot_legend                
            mplot_deco.addObject(self.plot_legend);

            self.draw_canvas( mplot_deco, "plots_%s_%s_%s_%s/other/"%(self.additioninformation, self.categoryLabel,self.PS_model, self.wtagger_label), "m_lvj"+label+in_range+in_range+mlvj_model+"_deco",0,logy)

        ### Number of the event in the dataset and lumi scale factor --> set the proper number for bkg extraction or for signal region
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_mlvj").Print();
        self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).Print()
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_mlvj").setVal( self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_mlvj").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_mlvj").setError(self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_mlvj").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).getVal() )

        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_mlvj").Print();

        #### Call the evaluation of the normalization in the signal region for signal, TTbar, VV, SingleT, and WJets after the extrapolation via alpha
        self.get_pdf_signalregion_integral(label) ;


    ##### Function that calculate the normalization inside the mlvj signal region (mass window around the resonance in order to fill datacards)
    def get_pdf_signalregion_integral(self, label, model_name=""):

        print "############### get mlvj normalization inside SR ",label," ",model_name," ##################"
        if model_name == "":
            model = self.workspace4fit_.pdf("model"+label+"_signalregion"+"_"+self.categoryLabel+"_mlvj");
        else:
            model = self.workspace4fit_.pdf(model_name);

        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj");

        fullInt   = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj) );
        signalInt = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj),("signalregion"));

        fullInt_val = fullInt.getVal()
        signalInt_val = signalInt.getVal()/fullInt_val

        ## integal in the signal region
        print "######### integral in SR: ",label+"signalInt=%s"%(signalInt_val)

        print "####### Events Number in MC Dataset:"
        self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_mlvj").Print();
        self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.categoryLabel+"_mlvj").Print();

        print "########## Events Number get from fit:"
        rrv_tmp=self.workspace4fit_.var("rrv_number"+label+"_signalregion"+"_"+self.categoryLabel+"_mlvj");
        print "Events Number in Signal Region from fitting: %s"%(rrv_tmp.getVal()*signalInt_val)

        #### store the info in the output file
        self.file_out.write( "\n%s++++++++++++++++++++++++++++++++++++"%(label) )
        self.file_out.write( "\nEvents Number in All Region from dataset : %s"%(self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.categoryLabel+"_mlvj").getVal()) )
        self.file_out.write( "\nEvents Number in Signal Region from dataset: %s"%(self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_mlvj").getVal()) )
        self.file_out.write( "\nRatio signalregion/all_range from dataset :%s"%(self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_mlvj").getVal()/self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.categoryLabel+"_mlvj").getVal() ) )
        self.file_out.write( "\nEvents Number in All Region from fitting : %s\n"%(rrv_tmp.getVal()) )
        self.file_out.write( "\nEvents Number in Signal Region from fitting: %s\n"%(rrv_tmp.getVal()*signalInt_val) )
        self.file_out.write( "\nRatio signalregion/all_range from fitting :%s"%(signalInt_val ) )

        if not self.workspace4fit_.var("rrv_number_fitting_signalregion"+label+"_"+self.categoryLabel+"_mlvj"):
            rrv_number_fitting_signalregion_mlvj = RooRealVar("rrv_number_fitting_signalregion"+label+"_"+self.categoryLabel+"_mlvj","rrv_number_fitting_signalregion"+label+"_"+
                    self.categoryLabel+"_mlvj", rrv_tmp.getVal()*signalInt_val );
            getattr(self.workspace4fit_,"import")(rrv_number_fitting_signalregion_mlvj);
        else :
            self.workspace4fit_.var("rrv_number_fitting_signalregion"+label+"_"+self.categoryLabel+"_mlvj").setVal(rrv_tmp.getVal()*signalInt_val);

        self.workspace4fit_.var("rrv_number_fitting_signalregion"+label+"_"+self.categoryLabel+"_mlvj").Print();


    #### run selection on data to build the datasets 
    def get_data(self):
        print "############### get_data ########################"
        self.get_mj_and_mlvj_dataset(self.sig_bkg_files[2][2],"_data", self.analyzer_config["limit_variable"], self.analyzer_config["obs0_variable"])
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_dataset_signalregion_data_%s_mlvj"%(self.categoryLabel)).clone("observation_for_counting"))

    ##### Method used to cycle on the events and for the dataset to be fitted
    def get_mj_and_mlvj_dataset(self,in_file_name, label, lvj_mass, jet_mass ):# to get the shape of m_lvj,jet_mass="jet_mass_pr"
        print "################### get_mj_and_mlvj_dataset : ",in_file_name,"  ",label,"  ##################";
        fileIn_name = TString(in_file_name);
        #print fileIn_name;
        fileIn = TFile(fileIn_name.Data());
        treeIn = fileIn.Get("SelectedCandidatesPlain");

        rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j")
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rrv_weight   = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.)

        ##### dataset of m_j -> scaleed and not scaled to lumi 
        rdataset_mj     = RooDataSet("rdataset"+label+"_"+self.categoryLabel+"_mj","rdataset"+label+"_"+self.categoryLabel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_mj = RooDataSet("rdataset4fit"+label+"_"+self.categoryLabel+"_mj","rdataset4fit"+label+"_"+self.categoryLabel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        ##### dataset of m_lvj -> scaled and not scaled to lumi in different region
        rdataset_lowersideband_mlvj = RooDataSet("rdataset"+label+"_lowersideband"+"_"+self.categoryLabel+"_mlvj","rdataset"+label+"_lowersideband"+"_"+self.categoryLabel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset_signalregion_mlvj = RooDataSet("rdataset"+label+"_signalregion"+"_"+self.categoryLabel+"_mlvj","rdataset"+label+"_signalregion"+"_"+self.categoryLabel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset_uppersideband_mlvj = RooDataSet("rdataset"+label+"_uppersideband"+"_"+self.categoryLabel+"_mlvj","rdataset"+label+"_uppersideband"+"_"+self.categoryLabel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset4fit_lowersideband_mlvj = RooDataSet("rdataset4fit"+label+"_lowersideband"+"_"+self.categoryLabel+"_mlvj","rdataset4fit"+label+"_lowersideband"+"_"+self.categoryLabel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset4fit_signalregion_mlvj = RooDataSet("rdataset4fit"+label+"_signalregion"+"_"+self.categoryLabel+"_mlvj","rdataset4fit"+label+"_signalregion"+"_"+self.categoryLabel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset4fit_uppersideband_mlvj = RooDataSet("rdataset4fit"+label+"_uppersideband"+"_"+self.categoryLabel+"_mlvj","rdataset4fit"+label+"_uppersideband"+"_"+self.categoryLabel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### categorize the event in sideband and signal region --> combined dataset 

        data_category = RooCategory("data_category","data_category");
        data_category.defineType("sideband");
        data_category.defineType("signalregion");
        combData = RooDataSet("combData"+label+"_"+self.categoryLabel,"combData"+label+"_"+self.categoryLabel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
        combData4fit = RooDataSet("combData4fit"+label+"_"+self.categoryLabel,"combData4fit"+label+"_"+self.categoryLabel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );

        obs0_variable_lowersideband_range_min= rrv_mass_j.getMin("lowersideband");
        obs0_variable_lowersideband_range_max= rrv_mass_j.getMax("lowersideband");
        obs0_variable_signalregion_range_min= rrv_mass_j.getMin("signalregion");
        obs0_variable_signalregion_range_max= rrv_mass_j.getMax("signalregion");
        obs0_variable_uppersideband_range_min= rrv_mass_j.getMin("uppersideband");
        obs0_variable_uppersideband_range_max= rrv_mass_j.getMax("uppersideband");

        limit_variable_signalregion_range_min = rrv_mass_lvj.getMin("signalregion");
        limit_variable_signalregion_range_max = rrv_mass_lvj.getMax("signalregion");
        print limit_variable_signalregion_range_min, limit_variable_signalregion_range_max

        hnum_4region=TH1D("hnum_4region"+label+"_"+self.categoryLabel,"hnum_4region"+label+"_"+self.categoryLabel,4,-1.5,2.5);# m_j -1: lowersideband; 0:signalregion; 1: uppersideband; 2:total
        hnum_2region=TH1D("hnum_2region"+label+"_"+self.categoryLabel,"hnum_2region"+label+"_"+self.categoryLabel,2,-0.5,1.5);# m_lvj 0: signalregion; 1: total

        tmp_lumi=self.GetLumi()*1000;

        print "###### N entries: ", treeIn.GetEntries()
        for i in range(treeIn.GetEntries()):
            if i % 100000 == 0: print "iEntry: ",i
            treeIn.GetEntry(i);

            if i==0:
                tmp_scale_to_lumi=treeIn.LumiWeight*tmp_lumi;

            tmp_jet_mass=getattr(treeIn, jet_mass);
            tmp_lvj_mass=getattr(treeIn, lvj_mass);

            self.isGoodEvent = 0 ;   
            if treeIn.categories==self.categoryID and tmp_lvj_mass> rrv_mass_lvj.getMin() and tmp_lvj_mass<rrv_mass_lvj.getMax() and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() :
                self.isGoodEvent = 1 ;   

            if self.isGoodEvent == 1:
                ### weigh MC events              
                tmp_event_weight     = treeIn.weight*tmp_lumi;
                tmp_event_weight4fit = treeIn.HLTweight*treeIn.PUweight*treeIn.GenWeight*treeIn.BTagWeight*treeIn.VTagWeight;
                tmp_event_weight4fit = tmp_event_weight4fit*treeIn.LumiWeight*tmp_lumi/tmp_scale_to_lumi;
                #### wtagger_eff_reweight
                if TString(label).Contains("data") :
                    tmp_event_weight=1.;
                    tmp_event_weight4fit=1.;
                #print i, tmp_lvj_mass, tmp_jet_mass, tmp_event_weight

                rrv_mass_lvj.setVal(tmp_lvj_mass);

                if tmp_jet_mass >= obs0_variable_lowersideband_range_min and tmp_jet_mass < obs0_variable_lowersideband_range_max:
                    rdataset_lowersideband_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                    rdataset4fit_lowersideband_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );

                    data_category.setLabel("sideband");
                    combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                    combData4fit.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);

                if tmp_jet_mass >= obs0_variable_signalregion_range_min and tmp_jet_mass < obs0_variable_signalregion_range_max:
                    rdataset_signalregion_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                    rdataset4fit_signalregion_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );

                    data_category.setLabel("signalregion");
                    combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                    combData4fit.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);
                    hnum_2region.Fill(1,tmp_event_weight);

                    if tmp_lvj_mass >=limit_variable_signalregion_range_min and tmp_lvj_mass <limit_variable_signalregion_range_max:
                        hnum_2region.Fill(0,tmp_event_weight);

                if tmp_jet_mass >= obs0_variable_uppersideband_range_min and tmp_jet_mass < obs0_variable_uppersideband_range_max:
                    rdataset_uppersideband_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                    rdataset4fit_uppersideband_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );

                rrv_mass_j.setVal( tmp_jet_mass );
                rdataset_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
                rdataset4fit_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );

                if tmp_jet_mass >=obs0_variable_lowersideband_range_min and tmp_jet_mass <obs0_variable_lowersideband_range_max:
                    hnum_4region.Fill(-1,tmp_event_weight );
                if tmp_jet_mass >=obs0_variable_signalregion_range_min and tmp_jet_mass <obs0_variable_signalregion_range_max :
                    hnum_4region.Fill(0,tmp_event_weight);
                if tmp_jet_mass >=obs0_variable_uppersideband_range_min and tmp_jet_mass <obs0_variable_uppersideband_range_max:
                    hnum_4region.Fill(1,tmp_event_weight);

                hnum_4region.Fill(2,tmp_event_weight);

        ### scaler to lumi for MC in 4fit datasets
        rrv_scale_to_lumi=RooRealVar("rrv_scale_to_lumi"+label+"_"+self.categoryLabel,"rrv_scale_to_lumi"+label+"_"+self.categoryLabel,tmp_scale_to_lumi)
        rrv_scale_to_lumi.Print()
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi)

        ### prepare m_lvj dataset to be compared with the fit results
        rrv_number_dataset_signalregion_mlvj=RooRealVar("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_mlvj","rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_mlvj",hnum_2region.GetBinContent(1));
        rrv_number_dataset_AllRange_mlvj=RooRealVar("rrv_number_dataset_AllRange"+label+"_"+self.categoryLabel+"_mlvj","rrv_number_dataset_AllRange"+label+"_"+self.categoryLabel+"_mlvj",hnum_2region.GetBinContent(2));

        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signalregion_mlvj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_AllRange_mlvj)
        ### import the dataser       
        getattr(self.workspace4fit_,"import")(rdataset_lowersideband_mlvj);
        getattr(self.workspace4fit_,"import")(rdataset_signalregion_mlvj);
        getattr(self.workspace4fit_,"import")(rdataset_uppersideband_mlvj);
        getattr(self.workspace4fit_,"import")(rdataset4fit_lowersideband_mlvj);
        getattr(self.workspace4fit_,"import")(rdataset4fit_signalregion_mlvj);
        getattr(self.workspace4fit_,"import")(rdataset4fit_uppersideband_mlvj);
        getattr(self.workspace4fit_,"import")(combData);
        getattr(self.workspace4fit_,"import")(combData4fit);

        ### write in the output 
        self.file_out.write("\n%s events number in m_lvj from dataset: %s"%(label,rdataset_signalregion_mlvj.sumEntries()))

        ### prepare m_j dataset
        rrv_number_dataset_lowersideband_mj=RooRealVar("rrv_number_dataset_lowersideband"+label+"_"+self.categoryLabel+"_mj","rrv_number_dataset_lowersideband"+label+"_"+self.categoryLabel+"_mj",hnum_4region.GetBinContent(1));
        rrv_number_dataset_signalregion_mj=RooRealVar("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_mj","rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_mj",hnum_4region.GetBinContent(2));
        rrv_number_dataset_uppersideband_mj=RooRealVar("rrv_number_dataset_uppersideband"+label+"_"+self.categoryLabel+"_mj","rrv_number_dataset_uppersideband"+label+"_"+self.categoryLabel+"_mj",hnum_4region.GetBinContent(3));
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_lowersideband_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signalregion_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_uppersideband_mj)

        getattr(self.workspace4fit_,"import")(rdataset_mj)
        getattr(self.workspace4fit_,"import")(rdataset4fit_mj)

        #### print everything
        rdataset_lowersideband_mlvj.Print();
        rdataset_signalregion_mlvj.Print();
        rdataset_uppersideband_mlvj.Print();
        rdataset_mj.Print();
        rdataset4fit_lowersideband_mlvj.Print();
        rdataset4fit_signalregion_mlvj.Print();
        rdataset4fit_uppersideband_mlvj.Print();
        rdataset4fit_mj.Print();
        rrv_number_dataset_signalregion_mlvj.Print()
        rrv_number_dataset_AllRange_mlvj.Print()
        rrv_number_dataset_lowersideband_mj.Print()
        rrv_number_dataset_signalregion_mj.Print()
        rrv_number_dataset_uppersideband_mj.Print()
        print rdataset_signalregion_mlvj.sumEntries()
        print rrv_number_dataset_signalregion_mlvj.getVal()
        print rrv_number_dataset_AllRange_mlvj.getVal()

    ### Define the Extended Pdf for and mJ fit giving: label, fit model name, list constraint and ismc
    #def make_Model(self, label, fit_config, mass_spectrum="_mj", ConstraintsList=[], ismc_wjet=0, area_init_value=500):
    def make_Model(self, label, fit_config, mass_spectrum, area_init_value=500):

        ##### define an extended pdf from a standard Roofit One
        print " "
        print "###############################################"
        print "## Make model : ",label," ",fit_config,"##";
        print "###############################################"
        print " "

        rrv_number = RooRealVar("rrv_number"+label+"_"+self.categoryLabel+mass_spectrum,"rrv_number"+label+"_"+self.categoryLabel+mass_spectrum,area_init_value,0.,1e7);
        ## call the make RooAbsPdf method
        #model_pdf = self.make_Pdf(label, fit_config,mass_spectrum,ConstraintsList,ismc_wjet)
        model_pdf = make_Pdf(label+"_"+self.categoryLabel, self.workspace4fit_, fit_config,mass_spectrum)
        print "######## Model Pdf ########"        
        model_pdf.Print();

        ## create the extended pdf
        model = RooExtendPdf("model"+label+"_"+self.categoryLabel+mass_spectrum,"model"+label+"_"+self.categoryLabel+mass_spectrum, model_pdf, rrv_number );
        print "######## Model Extended Pdf ########"        

        #### put all the parameters ant the shape in the workspace
        getattr(self.workspace4fit_,"import")(rrv_number)
        getattr(self.workspace4fit_,"import")(model) 
        self.workspace4fit_.pdf("model"+label+"_"+self.categoryLabel+mass_spectrum).Print();
        ## return the total extended pdf
        return self.workspace4fit_.pdf("model"+label+"_"+self.categoryLabel+mass_spectrum);


    #### Read the final workspace and produce the latest plots 
    def read_workspace(self, logy=0):
        print "================ read workspace ================";

        ### Taket the workspace for limits  
        file = TFile(self.file_rlt_root) ;
        workspace = file.Get("workspace4limit_") ;
        workspace.Print("v");

        ### iterate on the workspace element parameters
        print "----------- Parameter Workspace -------------";
        parameters_workspace = workspace.allVars();
        par = parameters_workspace.createIterator();
        par.Reset();
        param = par.Next()
        while (param):
            param.Print();
            param=par.Next();
        print "---------------------------------------------";

        workspace.data("data_obs_%s"%(self.categoryLabel)).Print()

        print "----------- Pdf in the Workspace -------------";
        pdfs_workspace = workspace.allPdfs();
        par = pdfs_workspace.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            param.Print();
            param = par.Next()
        print "----------------------------------------------";

        rrv_x = workspace.var("rrv_mass_lvj")
        data_obs = workspace.data("data_obs_%s"%(self.categoryLabel));

        model_pdf_bkgs=[];
        matrix_model_pdf_bkgs=[];
        rate_bkgs=[];
        ral_pdf_bkgs=RooArgList();#ral: RooArgList
        ral_rate_bkgs=RooArgList();
        rate_total_bkgs=0.;
        rate_error2_total_bkgs=0.;#error2: error^2
        for iter in range(self.nbkg):
            pdf= workspace.pdf( "%s_%s"%(self.bkg_list[iter][0],self.categoryLabel) );
            pdf.Print();
            model_pdf_bkgs.append(pdf);
            matrix_model_pdf_bkgs.append(pdf.GetName());
            for jter in range(iter):
                matrix_model_pdf_bkgs[jter]+=",%s"%(pdf.GetName());
            rrv_rate= workspace.var("rate_%s_for_unbin"%(self.bkg_list[iter][0]));
            rrv_rate.Print();
            rate_bkgs.append( rrv_rate);
            ral_pdf_bkgs.add(pdf);
            ral_rate_bkgs.add(rrv_rate);
            rate_total_bkgs+=rrv_rate.getVal();
            rate_error2_total_bkgs+=rrv_rate.getError()*rrv_rate.getError();
        print matrix_model_pdf_bkgs;

        #### Prepare the final plot starting from total background 
        rrv_number_Total_bkgs = RooRealVar("rrv_number_Total_bkgs","rrv_number_Total_bkgs", rate_total_bkgs);
        rrv_number_Total_bkgs.setError( TMath.Sqrt( rate_error2_total_bkgs) );
        rrv_number_Total_bkgs.Print()

        #### Total pdf 
        model_Total_bkgs = RooAddPdf("model_Total_bkgs","model_Total_bkgs", ral_pdf_bkgs, ral_rate_bkgs);
        model_Total_bkgs.Print();

        #### scale factor in order to scale MC to data in the final plot -> in order to avoid the normalization to data which is done by default in rooFit
        scale_number_Total_bkgs = rrv_number_Total_bkgs.getVal()/data_obs.sumEntries()
        print scale_number_Total_bkgs, rrv_number_Total_bkgs.getVal(), data_obs.sumEntries()

        #### create the frame
        mplot = rrv_x.frame(RooFit.Title("check_workspace"), RooFit.Bins(int(rrv_x.getBins()/self.BinWidth_narrow_factor)));
        data_obs.plotOn(mplot , RooFit.Name("data_invisible"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0));

        for iter in range(self.nbkg):
            model_Total_bkgs.plotOn(mplot,RooFit.Normalization(scale_number_Total_bkgs), RooFit.Components(matrix_model_pdf_bkgs[iter]), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet[self.bkg_list[iter][0]]), RooFit.Name(self.bkg_list[iter][0]), RooFit.LineColor(kBlack), RooFit.VLines());
            model_Total_bkgs.plotOn(mplot,RooFit.Normalization(scale_number_Total_bkgs), RooFit.Components(matrix_model_pdf_bkgs[iter]), RooFit.Name(self.bkg_list[iter][0]+"_line"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

        for iter in range(self.nsig):
            model_pdf=workspace.pdf( "%s_%s"%(self.sig_list[iter][0],self.categoryLabel) );
            model_pdf.Print();
            rrv_rate= workspace.var("rate_%s_for_unbin"%(self.sig_list[iter][0]));
            rrv_rate.Print();

            ### signal scale to be visible in the plots
            signal_scale=self.signal_scale;
            scale_number_signal = rrv_rate.getVal()/data_obs.sumEntries()*signal_scale;
            model_pdf.plotOn(mplot,RooFit.Normalization(scale_number_signal),RooFit.Name("%s #times %s"%(self.sig_list[iter][0], signal_scale)),RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet[self.sig_list[0][0]]), RooFit.LineStyle(2), RooFit.VLines());

        #### plot the observed data using poissonian error bar
        self.getData_PoissonInterval(data_obs,mplot);

        model_Total_bkgs.plotOn(mplot,RooFit.Normalization(scale_number_Total_bkgs),RooFit.Invisible());
        mplot_pull=self.get_pull(rrv_x,mplot);

        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        draw_error_band(model_Total_bkgs, rrv_x.GetName(), rrv_number_Total_bkgs,self.FloatingParams,workspace ,mplot,self.color_palet["Uncertainty"],"F");

        mplot.Print();
        self.plot_legend = self.legend4Plot(mplot,0,1,-0.01,-0.05,0.11,0.);
        self.plot_legend.SetTextSize(0.036);
        mplot.addObject(self.plot_legend);

        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);


        parameters_list = RooArgList();
        self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s/limit_variable/"%(self.additioninformation, self.categoryLabel),"check_workspace_for_limit","",0,1);

        if workspace.var("rrv_num_floatparameter_in_last_fitting"):
            self.nPar_float_in_fitTo = int(workspace.var("rrv_num_floatparameter_in_last_fitting").getVal());
        else:
            self.nPar_float_in_fitTo = self.FloatingParams.getSize();
        nBinX = mplot.GetNbinsX();
        ndof  = nBinX-self.nPar_float_in_fitTo;
        print "nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare( self.nPar_float_in_fitTo )*ndof, ndof );

  ### in order to make the plot_legend
    def legend4Plot(self, plot, left=1, isFill=1, x_offset_low=0., y_offset_low=0., x_offset_high =0., y_offset_high =0., TwoCoulum =1.):
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
                    leg_name="CMS Data "+self.categoryTitle; 
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


    def getData_PoissonInterval(self,data_obs,mplot):
        rrv_x = self.workspace4fit_.var("rrv_mass_lvj");
        datahist   = data_obs.binnedClone(data_obs.GetName()+"_binnedClone",data_obs.GetName()+"_binnedClone");
        data_histo = datahist.createHistogram("histo_data",rrv_x) ;
        data_histo.SetName("data");
        data_plot  = RooHist(data_histo);
        data_plot.SetMarkerStyle(20);
        data_plot.SetMarkerSize(1.5);

        alpha = 1 - 0.6827;
        for iPoint  in range(data_plot.GetN()):
            N = data_plot.GetY()[iPoint];
            if N==0 : L = 0;
            else : L = (ROOT.Math.gamma_quantile(alpha/2,N,1.));
            U =  ROOT.Math.gamma_quantile_c(alpha/2,N+1,1);
            data_plot.SetPointEYlow(iPoint, N-L);
            data_plot.SetPointEYhigh(iPoint,U-N);
            data_plot.SetPointEXlow(iPoint,0);        
            data_plot.SetPointEXhigh(iPoint,0);        

        mplot.addPlotable(data_plot,"PE");


    ### in order to get the pull
    def get_pull(self, rrv_x, mplot_orig):

        print "############### draw the pull plot ########################"
        hpull = mplot_orig.pullHist();
        x = ROOT.Double(0.); y = ROOT.Double(0) ;
        for ipoint in range(0,hpull.GetN()):
            hpull.GetPoint(ipoint,x,y);
            if(y == 0): hpull.SetPoint(ipoint,x,10)

        mplot_pull = rrv_x.frame(RooFit.Title("Pull Distribution"), RooFit.Bins(int(rrv_x.getBins()/self.BinWidth_narrow_factor)));
        medianLine = TLine(rrv_x.getMin(),0.,rrv_x.getMax(),0); medianLine.SetLineWidth(2); medianLine.SetLineColor(kRed);
        mplot_pull.addObject(medianLine);
        mplot_pull.addPlotable(hpull,"P");
        mplot_pull.SetTitle("");
        mplot_pull.GetXaxis().SetTitle("");
        mplot_pull.GetYaxis().SetRangeUser(-5,5);
        mplot_pull.GetYaxis().SetTitleSize(0.10);
        mplot_pull.GetYaxis().SetLabelSize(0.10);
        mplot_pull.GetXaxis().SetTitleSize(0.10);
        mplot_pull.GetXaxis().SetLabelSize(0.10);
        mplot_pull.GetYaxis().SetTitleOffset(0.40);
        mplot_pull.GetYaxis().SetTitle("#frac{Data-Fit}{#sigma_{data}}");
        mplot_pull.GetYaxis().CenterTitle();

        return mplot_pull;
       #### in order to make the banner on the plots
    def banner4Plot(self, iswithpull=0):
        print "############### draw the banner ########################"

        if iswithpull:
            if TString(self.categoryLabel).Contains("el"):
                #        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e #nu "%(self.GetLumi())));
                banner = TLatex(0.18,0.96,"CMS                                             L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
            elif TString(self.categoryLabel).Contains("mu"):
                #        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
                banner = TLatex(0.18,0.96,"CMS                                             L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
            elif TString(self.categoryLabel).Contains("em"):
                #        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu,e #nu "%(self.GetLumi())));
                banner = TLatex(0.18,0.96,"CMS                                             L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
            banner.SetNDC(); banner.SetTextSize(0.041);
        else:
            if TString(self.categoryLabel).Contains("el"):
                #        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e #nu "%(self.GetLumi())));
                banner = TLatex(0.18,0.96,"CMS                                         L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
            if TString(self.categoryLabel).Contains("mu"):
                #        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
                banner = TLatex(0.18,0.96,"CMS                                         L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
            if TString(self.categoryLabel).Contains("em"):
                #        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu,e #nu "%(self.GetLumi())));
                banner = TLatex(0.18,0.96,"CM                                          L = 19.7 fb^{-1} at #sqrt{s} = 8 TeV");
            banner.SetNDC(); banner.SetTextSize(0.033);

        return banner;


    ##### Prepare the workspace for the limit and to store info to be printed in the datacard
    #def prepare_limit(self,mode, isTTbarFloating=0, isVVFloating=0, isSingleTFloating=0):
    def prepare_limit(self,analysis_mode):
        print "####################### prepare_limit for %s method ####################"%(analysis_mode);
        self.workspace4fit_.Print("v");

        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_mass_lvj"));

        #### number of events in signal region for every sigs and bkgs for cut-and-couting limit
        for iter in range(self.nsig):
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signalregion_%s_%s_mlvj"%(self.sig_list[iter][0], self.categoryLabel)).clone("rate_%s_for_counting"%(self.sig_list[iter][0]) ))
        for iter in range(self.nbkg):
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signalregion_%s_%s_mlvj"%(self.bkg_list[iter][0], self.categoryLabel)).clone("rate_%s_for_counting"%(self.bkg_list[iter][0]) ))

        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signalregion_%s_%s_mlvj"%(self.allsignals,self.categoryLabel)).clone("rate_%s_for_counting"%(self.allsignals)))
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signalregion_WJets0_%s_mlvj"%(self.categoryLabel)).clone("rate_WJets_for_counting"))
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signalregion_VV_%s_mlvj"%(self.categoryLabel)).clone("rate_VV_for_counting"))
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signalregion_TTbar_%s_mlvj"%(self.categoryLabel)).clone("rate_TTbar_for_counting"))
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signalregion_SingleT_%s_mlvj"%(self.categoryLabel)).clone("rate_SingleT_for_counting"))

        ##### number of signal, Wjets, VV, TTbar and SingleT --> unbin
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_%s_signalregion_%s_mlvj"%(self.allsignals, self.categoryLabel)).clone("rate_%s_for_unbin"%(self.allsignals)));
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_WJets0_signalregion_%s_mlvj"%(self.categoryLabel)).clone("rate_WJets_for_unbin"));
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_VV_signalregion_%s_mlvj"%(self.categoryLabel)).clone("rate_VV_for_unbin"));
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_TTbar_signalregion_%s_mlvj"%(self.categoryLabel)).clone("rate_TTbar_for_unbin"));
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_SingleT_signalregion_%s_mlvj"%(self.categoryLabel)).clone("rate_SingleT_for_unbin"));

        #### Set the error properly -> taking into account lumi, Vtagger and theoretical uncertainty on XS -> for VV, TTbar and SingleT
        #self.workspace4limit_.var("rate_VV_for_unbin").setError(self.workspace4limit_.var("rate_VV_for_unbin").getVal()*TMath.Sqrt( self.lumi_uncertainty*self.lumi_uncertainty + self.rrv_wtagger_eff_reweight_forV.getError()/self.rrv_wtagger_eff_reweight_forV.getVal()*self.rrv_wtagger_eff_reweight_forV.getError()/self.rrv_wtagger_eff_reweight_forV.getVal() +self.XS_VV_uncertainty*self.XS_VV_uncertainty ) );
        #self.workspace4limit_.var("rate_SingleT_for_unbin").setError(self.workspace4limit_.var("rate_SingleT_for_unbin").getVal()*TMath.Sqrt( self.lumi_uncertainty*self.lumi_uncertainty + self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal()*self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal() +self.XS_SingleT_uncertainty*self.XS_SingleT_uncertainty ) );
        #self.workspace4limit_.var("rate_TTbar_for_unbin").setError(self.workspace4limit_.var("rate_TTbar_for_unbin").getVal()*TMath.Sqrt( self.lumi_uncertainty*self.lumi_uncertainty + self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal()*self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal() ))

        #### Get the dataset for data into the signal region
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.data("rdataset_data_signalregion_%s_mlvj"%(self.categoryLabel)).Clone("data_obs_%s"%(self.categoryLabel)))

        for iter_sig in range(self.nsig):
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_%s_signalregion_%s_mlvj"%(self.sig_list[iter_sig][0], self.categoryLabel)).clone("rate_%s_for_unbin"%(self.sig_list[iter_sig][0])));
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_%s_signalregion_%s_mlvj"%(self.sig_list[iter_sig][0], self.categoryLabel)).clone("%s_%s"%(self.sig_list[iter_sig][0], self.categoryLabel)));

        for iter_bkg in range(self.nbkg):
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_%s_signalregion_%s_mlvj"%(self.bkg_list[iter_bkg][0], self.categoryLabel)).clone("rate_%s_for_unbin"%(self.bkg_list[iter_bkg][0])));
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_%s_signalregion_%s_mlvj"%(self.bkg_list[iter_bkg][0], self.categoryLabel)).clone("%s_%s"%(self.bkg_list[iter_bkg][0], self.categoryLabel)));

        #### Take the corrected pdf from the alpha method for the WJets
        #if analysis_mode=="sideband_correction_method1":
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_WJets0_signalregion_%s_after_correct_mlvj"%(self.categoryLabel)).clone("WJets_%s_%s"%(self.categoryLabel, self.wtagger_label)));

        #if isTTbarFloating:
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_TTbar_signalregion_%s_mlvj_Deco_TTbar_signalregion_%s_%s_mlvj"%(self.categoryLabel, self.categoryLabel, self.wtagger_label)).clone("TTbar_%s_%s"%(self.categoryLabel,self.wtagger_label)))
        #else :
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_TTbar_signalregion_%s_mlvj"%(self.categoryLabel, self.categoryLabel, self.wtagger_label)).clone("TTbar_%s_%s"%(self.categoryLabel,self.wtagger_label)))

        #if isSingleTFloating :     
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_SingleT_signalregion_%s_mlvj_Deco_SingleT_signalregion_%s_%s_mlvj"%(self.categoryLabel, self.categoryLabel, self.wtagger_label)).clone("SingleT_%s_%s"%(self.categoryLabel,self.wtagger_label)))
        #else :
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_SingleT_signalregion_%s_mlvj"%(self.categoryLabel)).clone("SingleT_%s_%s"%(self.categoryLabel,self.wtagger_label)))

        #if isVVFloating :    
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_VV_signalregion_%s_mlvj_Deco_VV_signalregion_%s_%s_mlvj"%(self.categoryLabel, self.categoryLabel, self.wtagger_label)).clone("VV_%s_%s"%(self.categoryLabel,self.wtagger_label)))
        #else:
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_VV_signalregion_%s_mlvj"%(self.categoryLabel)).clone("VV_%s_%s"%(self.categoryLabel,self.wtagger_label)))

        #if TString(self.allsignals).Contains("BulkG_WW"):
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_%s_signalregion_%s_mlvj"%(self.allsignals,self.categoryLabel)).clone("BulkWW_%s_%s"%(self.categoryLabel, self.wtagger_label)))
        #else:    
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_%s_signalregion_%s_mlvj"%(self.allsignals,self.categoryLabel)).clone(self.allsignals+"_%s_%s"%(self.categoryLabel, self.wtagger_label)))

        #### Fix all the Pdf parameters 
        #rrv_x = self.workspace4limit_.var("rrv_mass_lvj");

        #self.fix_Pdf(self.workspace4limit_.pdf("TTbar_%s_%s"%(self.categoryLabel,self.wtagger_label)), RooArgSet(rrv_x) );            
        #self.fix_Pdf(self.workspace4limit_.pdf("SingleT_%s_%s"%(self.categoryLabel,self.wtagger_label)), RooArgSet(rrv_x));
        #self.fix_Pdf(self.workspace4limit_.pdf("VV_%s_%s"%(self.categoryLabel,self.wtagger_label)), RooArgSet(rrv_x));
        #self.fix_Pdf(self.workspace4limit_.pdf("WJets_%s_%s"%(self.categoryLabel,self.wtagger_label)), RooArgSet(rrv_x));

        #if TString(self.allsignals).Contains("BulkG_WW"):            
        #    self.fix_Pdf(self.workspace4limit_.pdf("BulkWW_%s_%s"%(self.categoryLabel, self.wtagger_label)), RooArgSet(rrv_x));
        #else:    
        #    self.fix_Pdf(self.workspace4limit_.pdf(self.allsignals+"_%s_%s"%(self.categoryLabel, self.wtagger_label)), RooArgSet(rrv_x));

        #print " ############## Workspace for limit ";
        #parameters_workspace = self.workspace4limit_.allVars();
        #par = parameters_workspace.createIterator();
        #par.Reset();
        #param = par.Next()
        #while (param):
        #    param.Print();
        #    param=par.Next()

        params_list = [];
        #### main modality for the alpha function method
        #if analysis_mode=="sideband_correction_method1":

        #    if self.MODEL_4_mlvj=="ErfExp_v1" or self.MODEL_4_mlvj=="ErfPow_v1" or self.MODEL_4_mlvj=="2Exp" :
        #        ### uncertainty inflation on the Wjets shape from fitting data in lowersideband
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig3"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        ### Add to the parameter list
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig3"%(self.categoryLabel, self.wtagger_label)));

        #        ### Do the same for alpha paramter
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig3"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig4"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig5"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        ### Add to the parameter list
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig3"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig4"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0sim_%s_%s_mlvj_eig5"%(self.categoryLabel, self.wtagger_label)))

        #        ### Do the same for the TTbar
        #        if isTTbarFloating !=0 :
        #            self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #         self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #         self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);

        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)));
        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)));
        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)));


        #    if self.MODEL_4_mlvj=="ErfPow2_v1" or self.MODEL_4_mlvj=="ErfPowExp_v1" :
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig3"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig4"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);

        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig3"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig4"%(self.categoryLabel, self.wtagger_label)));

        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig3"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig4"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig5"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig6"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig7"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);

        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig3"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig4"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig5"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig6"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig7"%(self.categoryLabel, self.wtagger_label)))


        #        if isTTbarFloating !=0 :
        #            self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #         self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #         self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #         self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig3"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);

        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)));
        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)));
        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)));
        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig3"%(self.categoryLabel, self.wtagger_label)));


        #    if self.MODEL_4_mlvj=="Exp" or self.MODEL_4_mlvj=="Pow" :

        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);

        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)));


        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);

        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)))

        #        if isTTbarFloating !=0 :
        #            self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)));

        #    if self.MODEL_4_mlvj=="ExpN" or self.MODEL_4_mlvj=="ExpTail" or self.MODEL_4_mlvj=="Pow2" :

        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);

        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)));

        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig3"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);

        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig2"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_mlvj_eig3"%(self.categoryLabel, self.wtagger_label)))


        #        ### TTbar use exp
        #        if isTTbarFloating !=0:
        #            print "##################### TTbar will float in the limit procedure + final plot ######################";
        #            self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #            params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)));

        #        ### VV use ExpTail:
        #        if isVVFloating !=0:
        #            print "##################### VV will float in the limit procedure + final plot ######################";
        #          self.workspace4limit_.var("Deco_VV_signalregion_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_VV);
        #          self.workspace4limit_.var("Deco_VV_signalregion_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_VV);
        #          params_list.append(self.workspace4limit_.var("Deco_VV_signalregion_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)));
        #          params_list.append(self.workspace4limit_.var("Deco_VV_signalregion_%s_%s_mlvj_eig1"%(self.categoryLabel, self.wtagger_label)));

        #        ### SingleT use Exp:
        #        if isSingleTFloating !=0:
        #            print "##################### SingleT will float in the limit procedure + final plot ######################";
        #          self.workspace4limit_.var("Deco_SingleT_signalregion_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_SingleT);
        #          params_list.append(self.workspace4limit_.var("Deco_SingleT_signalregion_%s_%s_mlvj_eig0"%(self.categoryLabel, self.wtagger_label)));


        ##### add signal shape parameters' uncertainty -> increase the uncertainty on the mean and the sigma since we are using a CB or a Double CB or a BWxDB or BWxCB
        #if self.workspace4limit_.var("rrv_mean_CB_%s_signalregion_%s_%s"%(self.allsignals, self.categoryLabel, self.wtagger_label)):

        #    self.workspace4limit_.var( "rrv_mean_shift_scale_lep_%s_signalregion_%s_%s"%(self.allsignals, self.categoryLabel, self.wtagger_label)).setError(self.mean_signal_uncertainty_lep_scale);
        #   self.workspace4limit_.var( "rrv_mean_shift_scale_jes_%s_signalregion_%s_%s"%(self.allsignals, self.categoryLabel, self.wtagger_label)).setError(self.mean_signal_uncertainty_jet_scale);
        #   self.workspace4limit_.var( "rrv_sigma_shift_lep_scale_%s_signalregion_%s_%s"%(self.allsignals, self.categoryLabel, self.wtagger_label)).setError(self.sigma_signal_uncertainty_lep_scale);
        #   self.workspace4limit_.var( "rrv_sigma_shift_jes_%s_signalregion_%s_%s"%(self.allsignals, self.categoryLabel, self.wtagger_label)).setError(self.sigma_signal_uncertainty_jet_scale);
        #   self.workspace4limit_.var( "rrv_sigma_shift_res_%s_signalregion_%s_%s"%(self.allsignals, self.categoryLabel, self.wtagger_label)).setError(self.sigma_signal_uncertainty_jet_res);

        #   if self.categoryLabel == "mu":
        #       self.workspace4limit_.var("CMS_sig_p1_scale_m").setError(1);
        #    self.workspace4limit_.var("CMS_sig_p2_scale_m").setError(1);
        #    params_list.append(self.workspace4limit_.var("CMS_sig_p1_scale_m"));
        #    params_list.append(self.workspace4limit_.var("CMS_sig_p2_scale_m"));
        #elif self.categoryLabel == "el":
        #    self.workspace4limit_.var("CMS_sig_p1_scale_e").setError(1);
        #    self.workspace4limit_.var("CMS_sig_p2_scale_e").setError(1);
        #    params_list.append(self.workspace4limit_.var("CMS_sig_p1_scale_e"));
        #    params_list.append(self.workspace4limit_.var("CMS_sig_p2_scale_e"));
        #elif self.categoryLabel == "em":
        #    self.workspace4limit_.var("CMS_sig_p1_scale_em").setError(1);
        #    self.workspace4limit_.var("CMS_sig_p2_scale_em").setError(1);
        #    params_list.append(self.workspace4limit_.var("CMS_sig_p1_scale_em"));
        #    params_list.append(self.workspace4limit_.var("CMS_sig_p2_scale_em"));

        #   self.workspace4limit_.var("CMS_sig_p1_jes").setError(1);
        #   self.workspace4limit_.var("CMS_sig_p2_jes").setError(1);
        #   self.workspace4limit_.var("CMS_sig_p2_jer").setError(1);

        #   params_list.append(self.workspace4limit_.var("CMS_sig_p1_jes"));
        #   params_list.append(self.workspace4limit_.var("CMS_sig_p2_jes"));
        #   params_list.append(self.workspace4limit_.var("CMS_sig_p2_jer"));


        #### calculate the shape uncertainty for cut-and-counting
        #self.rrv_counting_uncertainty_from_shape_uncertainty = RooRealVar("rrv_counting_uncertainty_from_shape_uncertainty_%s"%(self.categoryLabel),"rrv_counting_uncertainty_from_shape_uncertainty_%s"%(self.categoryLabel),0);
        #self.rrv_counting_uncertainty_from_shape_uncertainty.setError( Calc_error("WJets_%s_%s"%(self.categoryLabel,self.wtagger_label), "rrv_mass_lvj" ,self.FloatingParams,self.workspace4limit_,"signalregion") );
        #self.rrv_counting_uncertainty_from_shape_uncertainty.Print();

        print " param list ",params_list ;

        ### Print the datacard for unbin and couting analysis
        self.print_limit_datacard("unbin",params_list);
        self.print_limit_datacard("counting");
        ### Save the workspace
        self.save_workspace_to_file();



    #### Method used in order to save the workspace in a output root file
    def save_workspace_to_file(self):
        self.workspace4limit_.writeToFile(self.file_rlt_root);
        self.file_out.close()


    #### Method used to print the general format of the datacard for both counting and unbinned analysis
    def print_limit_datacard(self, datacard_mode, params_list=[]):
        print "############## print_limit_datacard for %s ################"%(datacard_mode)
        if not (datacard_mode == "unbin" or datacard_mode == "counting"):
            print "print_limit_datacard use wrong mode: %s"%(datacard_mode);raw_input("ENTER");

        ### open the datacard    
        datacard_out = open(getattr(self,"file_datacard_%s"%(datacard_mode)),"w");

        ### start to print inside 
        datacard_out.write( "imax %s"%(self.nsig) )
        datacard_out.write( "\njmax %s"%(self.nbkg) )
        datacard_out.write( "\nkmax *" )
        datacard_out.write( "\n--------------- ")

        #if datacard_mode == "unbin":
        #    fnOnly = ntpath.basename(self.file_rlt_root) ## workspace for limit --> output file for the workspace
        #    if TString(self.allsignals).Contains("BulkG_WW"):
        #        datacard_out.write("\nshapes BulkWW  CMS_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.categoryLabel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.categoryLabel, self.wtagger_label));
        #    else:
        #        datacard_out.write("\nshapes %s  CMS_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.allsignals,self.categoryLabel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.categoryLabel, self.wtagger_label));

        #    datacard_out.write("\nshapes WJets  CMS_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.categoryLabel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.categoryLabel, self.wtagger_label));
        #    datacard_out.write("\nshapes TTbar  CMS_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.categoryLabel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.categoryLabel, self.wtagger_label));
        #    datacard_out.write("\nshapes SingleT   CMS_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.categoryLabel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.categoryLabel, self.wtagger_label));
        #    datacard_out.write("\nshapes VV     CMS_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.categoryLabel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.categoryLabel, self.wtagger_label));
        #    datacard_out.write("\nshapes data_obs   CMS_%s1J%s  %s %s:$PROCESS_%s_%s"%(self.categoryLabel,self.wtagger_label,fnOnly,self.workspace4limit_.GetName(), self.categoryLabel, self.wtagger_label));
        #    datacard_out.write( "\n--------------- ")

        datacard_out.write( "\nbin CMS_%s "%(self.categoryLabel));    
        #if datacard_mode == "unbin":
        #    datacard_out.write( "\nobservation %0.2f "%(self.workspace4limit_.data("data_obs_%s_%s"%(self.categoryLabel,self.wtagger_label)).sumEntries()) )
        if datacard_mode == "counting":
            datacard_out.write( "\nobservation %0.2f "%(self.workspace4limit_.var("observation_for_counting").getVal()) )

        datacard_out.write( "\n------------------------------" );

        tmp_bin_string="";
        tmp_processname_string="";
        tmp_processnum_string="";
        tmp_rate_counting="";
        for iter in range( self.nsig ):
            tmp_bin_string        +="CMS_%s "%(self.categoryLabel); 
            tmp_processname_string+="%s "%(self.sig_list[iter][0])
            tmp_processnum_string +="%s "%(1-self.nsig+iter)
            tmp_rate_counting     +="%0.5f "%( self.workspace4limit_.var("rate_%s_for_counting"%(self.sig_list[iter][0])).getVal() );

        for iter in range( self.nbkg ):
            tmp_bin_string        +="CMS_%s "%(self.categoryLabel); 
            tmp_processname_string+="%s "%(self.bkg_list[iter][0])
            tmp_processnum_string +="%s "%(1+iter)
            tmp_rate_counting     +="%0.5f "%( self.workspace4limit_.var("rate_%s_for_counting"%(self.bkg_list[iter][0])).getVal() );




        datacard_out.write( "\nbin "+tmp_bin_string );
        datacard_out.write( "\nprocess "+tmp_processname_string ); ## just one signal sample
        datacard_out.write( "\nprocess "+tmp_processnum_string );

        #### rates for the different process
        #if datacard_mode == "unbin":
        #    if TString(self.allsignals).Contains("BulkG_WW"):                    
        #        datacard_out.write( "\nrate %0.5f %0.3f %0.3f %0.3f %0.3f "%(self.workspace4limit_.var("rate_BulkWW_for_unbin").getVal()*self.xs_rescale, self.workspace4limit_.var("rate_WJets_for_unbin").getVal(), self.workspace4limit_.var("rate_TTbar_for_unbin").getVal(), self.workspace4limit_.var("rate_SingleT_for_unbin").getVal(), self.workspace4limit_.var("rate_VV_for_unbin").getVal() ) )
        #    else:
        #        datacard_out.write( "\nrate %0.5f %0.3f %0.3f %0.3f %0.3f "%(self.workspace4limit_.var("rate_%s_for_unbin"%(self.allsignals)).getVal()*self.xs_rescale, self.workspace4limit_.var("rate_WJets_for_unbin").getVal(), self.workspace4limit_.var("rate_TTbar_for_unbin").getVal(), self.workspace4limit_.var("rate_SingleT_for_unbin").getVal(), self.workspace4limit_.var("rate_VV_for_unbin").getVal() ) )

        if datacard_mode == "counting":
            datacard_out.write( "\nrate "+tmp_rate_counting );

        datacard_out.write( "\n-------------------------------- " )

        #### luminosity nouisance
        #datacard_out.write( "\nlumi_8TeV lnN %0.3f - %0.3f %0.3f %0.3f"%(1.+self.lumi_uncertainty, 1.+self.lumi_uncertainty,1.+self.lumi_uncertainty,1.+self.lumi_uncertainty) )

        #### SingleT XS  nouisance in boosted regime
        #datacard_out.write( "\nCMS_XS_SingleT lnN - - - %0.3f -"%(1+self.XS_SingleT_uncertainty) )

        #### VV XS  nouisance in boosted regime
        #datacard_out.write( "\nCMS_XS_VV lnN - - - - %0.3f"%(1+self.XS_VV_uncertainty) )

        #### WJets Normalization from data fit -> data driven
        #if self.number_WJets_insideband >0:
        #    datacard_out.write( "\nCMS_WJ_norm gmN %0.3f %0.3f - - -"%(self.number_WJets_insideband, getattr(self, "datadriven_alpha_WJets_%s"%(datacard_mode)) ) )
        #else:
        #    datacard_out.write( "\nCMS_WJ_norm_%s_%s lnN - %0.3f - - -"%(self.categoryLabel, self.wtagger_label, 1+ self.workspace4limit_.var("rate_WJets_for_unbin").getError()/self.workspace4limit_.var("rate_WJets_for_unbin").getVal() ) );

        #### Top normalization due to SF in the ttbar CR
        #datacard_out.write( "\nCMS_Top_norm_%s_%s lnN - - %0.3f %0.3f -"%(self.categoryLabel, self.wtagger_label, 1+self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal(), 1+self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal() ) );

        #### V-Tagger SF nouisance
        #if self.wtagger_label == "HP":
        #    datacard_out.write( "\nCMS_eff_vtag_tau21_sf lnN %0.3f/%0.3f - - - %0.3f/%0.3f"%(1+self.rrv_wtagger_eff_reweight_forV.getError(),1-self.rrv_wtagger_eff_reweight_forV.getError(), 1+self.rrv_wtagger_eff_reweight_forV.getError(),1-self.rrv_wtagger_eff_reweight_forV.getError()));
        #elif self.wtagger_label == "LP":
        #    datacard_out.write( "\nCMS_eff_vtag_tau21_sf lnN %0.3f/%0.3f - - - %0.3f/%0.3f"%(1-self.rrv_wtagger_eff_reweight_forV.getError(),1+self.rrv_wtagger_eff_reweight_forV.getError(), 1-self.rrv_wtagger_eff_reweight_forV.getError(),1+self.rrv_wtagger_eff_reweight_forV.getError()));
        #else:
        #    datacard_out.write( "\nCMS_eff_vtag_tau21_sf lnN %0.3f/%0.3f - - - %0.3f/%0.3f"%(1+self.rrv_wtagger_eff_reweight_forV.getError(),1-self.rrv_wtagger_eff_reweight_forV.getError(), 1+self.rrv_wtagger_eff_reweight_forV.getError(),1-self.rrv_wtagger_eff_reweight_forV.getError()));

        #### btag scale factor on the MC background
#       # datacard_out.write( "\nCMS_btagger lnN - - %0.3f %0.3f %0.3f"%(self.categoryLabel, 1+self.btag_scale_uncertainty, 1+self.btag_scale_uncertainty, 1+self.btag_scale_uncertainty ) );

        #### btag scale factor on the MC background
        #datacard_out.write( "\n#CMS_eff_vtag_model lnN %0.3f - - - %0.3f"%(1+self.eff_vtag_model,1+self.eff_vtag_model) );

        #### jet Mass effect only if available -> shapes changing due to the jet mass uncertainty (JEC for CA8/AK7) -> affects also WJets
        #if ( self.workspace4fit_.var("rrv_number_WJets0_massup_in_mj_signalregion_from_fitting_%s"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_WJets0_massdn_in_mj_signalregion_from_fitting_%s"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT_massup_%s_mj"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT_massdown_%s_mj"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar_massup_%s_mj"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar_massdn_%s_mj"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_massup_%s_mj"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_massdn_%s_mj"%(self.categoryLabel))) :

        #    datacard_out.write( "\nJetMass_%s lnN - %0.3f %0.3f %0.3f %0.3f"%(self.categoryLabel, 1+self.WJets_normlization_uncertainty_from_jet_mass, 1+self.TTbar_normlization_uncertainty_from_jet_mass, 1+self.SingleT_normlization_uncertainty_from_jet_mass, 1+self.VV_normlization_uncertainty_from_jet_mass ) )

        #if self.categoryLabel == "mu":
        #    self.channel_short = "m"
        #elif self.categoryLabel =="el":
        #    self.channel_short = "e"
        #elif self.categoryLabel =="em":
        #    self.channel_short = "em"

        #### trigger efficiency
        #datacard_out.write( "\nCMS_trigger_%s lnN %0.3f - %0.3f %0.3f %0.3f"%(self.channel_short, 1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty ) );

        #### Lepton SF
        #datacard_out.write( "\nCMS_eff_%s lnN %0.3f - %0.3f %0.3f %0.3f"%(self.channel_short, 1+self.lep_eff_uncertainty, 1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty ) );

        ############# Evaluated just for signal, in principle also on all the backgrounds with the same topology

        #### Lepton Energy scale
        #datacard_out.write( "\nCMS_scale_%s lnN %0.3f - - - -"%(self.channel_short,1+self.signal_lepton_energy_scale_uncertainty));

        #### Lepton Energy Resolution
        #datacard_out.write( "\nCMS_res_%s lnN %0.3f - - - -"%(self.channel_short,1+self.signal_lepton_energy_res_uncertainty));

        #### CA8 jet energy scale
        #datacard_out.write( "\nCMS_scale_j  lnN %0.3f - - - -"%(1+self.signal_jet_energy_scale_uncertainty));

        #### CA8 jet energy resolution
        #datacard_out.write( "\nCMS_res_j  lnN %0.3f - - - -"%(1+self.signal_jet_energy_res_uncertainty));

        #### btag on the signal
        #datacard_out.write( "\nCMS_btag_eff lnN %0.3f - - - -"%(1+self.signal_btag_uncertainty));


        #### print shapes parameter to be taken int account
        #if datacard_mode == "unbin":
        #    for ipar in params_list:
        #        print "Name %s",ipar.GetName();
        #      if TString(ipar.GetName()).Contains("Deco_TTbar_signalregion"):
        #          datacard_out.write( "\n%s param %0.1f %0.1f "%( ipar.GetName(), ipar.getVal(), ipar.getError() ) )
        #      else:
        #          datacard_out.write( "\n%s param %0.1f %0.1f "%( ipar.GetName(), ipar.getVal(), ipar.getError() ) )
        #if datacard_mode == "counting":
        #    datacard_out.write( "\nShape_%s_%s lnN - - %0.3f - - -"%(self.categoryLabel, self.wtagger_label, 1+self.rrv_counting_uncertainty_from_shape_uncertainty.getError()))

 
    def draw_canvas_with_pull(self, mplot, mplot_pull,parameters_list,in_directory, label, in_model_name="", show_constant_parameter=0, logy=0):# mplot + pull

        print "############### draw the canvas with pull ########################"
        mplot.GetXaxis().SetTitleOffset(1.1);
        mplot.GetYaxis().SetTitleOffset(1.3);
        mplot.GetXaxis().SetTitleSize(0.055);
        mplot.GetYaxis().SetTitleSize(0.055);
        mplot.GetXaxis().SetLabelSize(0.045);
        mplot.GetYaxis().SetLabelSize(0.045);
        mplot_pull.GetXaxis().SetLabelSize(0.14);
        mplot_pull.GetYaxis().SetLabelSize(0.14);
        mplot_pull.GetYaxis().SetTitleSize(0.15);
        mplot_pull.GetYaxis().SetNdivisions(205);


        cMassFit = TCanvas("cMassFit","cMassFit", 600,600);
        # if parameters_list is empty, don't draw pad3
        par_first=parameters_list.createIterator();
        par_first.Reset();
        param_first=par_first.Next()
        doParameterPlot = 0 ;
        if param_first and doParameterPlot != 0:
            pad1=TPad("pad1","pad1",0.,0. ,0.8,0.24);
            pad2=TPad("pad2","pad2",0.,0.24,0.8,1. );
            pad3=TPad("pad3","pad3",0.8,0.,1,1);
            pad1.Draw();
            pad2.Draw();
            pad3.Draw();
        else:
            pad1=TPad("pad1","pad1",0.,0. ,0.99,0.24);
            pad2=TPad("pad2","pad2",0.,0.24,0.99,1. );
            pad1.Draw();
            pad2.Draw();

        pad2.cd();
        mplot.Draw();
        banner = self.banner4Plot(1);
        banner.Draw();

        pad1.cd();
        mplot_pull.Draw();

        if param_first and doParameterPlot != 0:

            pad3.cd();
            latex=TLatex();
            latex.SetTextSize(0.1);
            par=parameters_list.createIterator();
            par.Reset();
            param=par.Next()
            i=0;
            while param:
                if (not param.isConstant() ) or show_constant_parameter:
                    param.Print();
                    icolor=1;#if a paramenter is constant, color is 2
                    if param.isConstant(): icolor=2
                    latex.DrawLatex(0,0.9-i*0.04,"#color[%s]{%s}"%(icolor,param.GetName()) );
                    latex.DrawLatex(0,0.9-i*0.04-0.02," #color[%s]{%4.3e +/- %2.1e}"%(icolor,param.getVal(),param.getError()) );
                    i=i+1;
                param=par.Next();

        ## create the directory where store the plots
        Directory = TString(in_directory+self.allsignals);
        if not Directory.EndsWith("/"):Directory = Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
            os.system("mkdir -p "+Directory.Data());

        rlt_file = TString(Directory.Data()+label);
        print "rlt_file=", rlt_file;
        if not in_model_name=="": in_model_name="_"+in_model_name;
        rlt_file=rlt_file.Append(in_model_name+"_with_pull.png");

        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".png",".pdf");
        cMassFit.SaveAs(rlt_file.Data());

        string_file_name = TString(label+in_model_name);

        if logy:
            mplot.GetYaxis().SetRangeUser(1e-3,mplot.GetMaximum()*200);
            pad2.SetLogy() ;
            pad2.Update();
            cMassFit.Update();
            rlt_file.ReplaceAll(".pdf","_log.pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());

        self.draw_canvas(mplot,in_directory,string_file_name.Data(),0,logy,1);

    #### jusr drawing canvas with no pull
    def draw_canvas(self, in_obj,in_directory, label, is_range=0, logy=0, frompull=0):

        print "############### draw the canvas without pull ########################"
        cMassFit = TCanvas("cMassFit","cMassFit", 600,600);

        if frompull and logy :
            in_obj.GetYaxis().SetRangeUser(1e-2,in_obj.GetMaximum()/200)
        elif not frompull and logy :
            in_obj.GetYaxis().SetRangeUser(0.00001,in_obj.GetMaximum())


        if is_range:
            h2=TH2D("h2","",100,400,1400,4,0.00001,4);
            h2.Draw();
            in_obj.Draw("same")
        else :
            in_obj.Draw()

        in_obj.GetXaxis().SetTitleSize(0.045);
        in_obj.GetXaxis().SetTitleOffset(1.15);
        in_obj.GetXaxis().SetLabelSize(0.04);

        in_obj.GetYaxis().SetTitleSize(0.055);
        in_obj.GetYaxis().SetTitleOffset(1.40);
        in_obj.GetYaxis().SetLabelSize(0.04);

        self.plot_legend.SetTextSize(0.031); 

        banner = self.banner4Plot();
        banner.Draw();

        Directory=TString(in_directory+self.allsignals);
        if not Directory.EndsWith("/"):Directory=Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
            os.system("mkdir -p "+Directory.Data());

        rlt_file=TString(Directory.Data()+label);

        rlt_file = rlt_file.Append(".png");
        cMassFit.SaveAs(rlt_file.Data());

        rlt_file.ReplaceAll(".png",".pdf");
        cMassFit.SaveAs(rlt_file.Data());

        if logy:
            in_obj.GetYaxis().SetRangeUser(1e-3,in_obj.GetMaximum()*200);
            cMassFit.SetLogy() ;
            cMassFit.Update();
            rlt_file.ReplaceAll(".pdf","_log.png");
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());

        #print "Show %s; ENTER to continue!"%(label);raw_input("ENTER");

    #### Method to make a RooAbsPdf giving label, model name, spectrum, if it is mc or not and a constraint list for the parameters          
    '''
    def make_Pdf(self, label, fit_config, mass_spectrum="_mj", ConstraintsList=[],ismc = 0):
        if TString(mass_spectrum).Contains("_mj"): rrv_x = self.workspace4fit_.var("rrv_mass_j");
        if TString(mass_spectrum).Contains("_mlvj"): rrv_x = self.workspace4fit_.var("rrv_mass_lvj");

        in_model_name=fit_config[1];

        # W mass: 80.385
        if in_model_name == "Voig":
            print "########### Voigtian Pdf for mJ ############"
            rrv_mean_voig=RooRealVar("rrv_mean_voig"+label+"_"+self.categoryLabel,"rrv_mean_voig"+label+"_"+self.categoryLabel,84,78,88);
            rrv_width_voig=RooRealVar("rrv_width_voig"+label+"_"+self.categoryLabel,"rrv_width_voig"+label+"_"+self.categoryLabel,7.,1,40);
            rrv_sigma_voig=RooRealVar("rrv_sigma_voig"+label+"_"+self.categoryLabel,"rrv_sigma_voig"+label+"_"+self.categoryLabel,5,0.01,20);
            model_pdf = RooVoigtian("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_mean_voig,rrv_width_voig,rrv_sigma_voig);

        # Higgs mass 600-1000
        if in_model_name == "Voig_v1":
            print "########### Voigtian Pdf for Higgs mlvj ############"
            rrv_mean_voig=RooRealVar("rrv_mean_voig"+label+"_"+self.categoryLabel,"rrv_mean_voig"+label+"_"+self.categoryLabel,650,550,1200);
            rrv_width_voig=RooRealVar("rrv_width_voig"+label+"_"+self.categoryLabel,"rrv_width_voig"+label+"_"+self.categoryLabel,100.,10,600);
            rrv_sigma_voig=RooRealVar("rrv_sigma_voig"+label+"_"+self.categoryLabel,"rrv_sigma_voig"+label+"_"+self.categoryLabel,200,10,400);
            model_pdf = RooVoigtian("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_mean_voig,rrv_width_voig,rrv_sigma_voig);

        # Bulk mass 600-1000
        if in_model_name == "Voig_v2":
            label_tstring=TString(label);
            print "########### Voigtian Pdf for Higgs mlvj ############"
            if label_tstring.Contains("600") and (not label_tstring.Contains("1600") ):
                rrv_mean_voig=RooRealVar("rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,600,500,700);
                rrv_width_voig.setConstant(kTRUE);
                rrv_sigma_voig=RooRealVar("rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,40,10,80);
   
            elif label_tstring.Contains("700") and (not label_tstring.Contains("1700") ):
                rrv_mean_voig=RooRealVar("rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,700,600,800);
                rrv_width_voig=RooRealVar("rrv_width_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2.5,0,10);
                rrv_width_voig.setConstant(kTRUE);
                rrv_sigma_voig=RooRealVar("rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,40,10,80);
   
            elif label_tstring.Contains("800") and (not label_tstring.Contains("1800") ):
                rrv_mean_voig=RooRealVar("rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,800,700,900);
                rrv_width_voig=RooRealVar("rrv_width_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2.5,0,10);
                rrv_width_voig.setConstant(kTRUE);
                rrv_sigma_voig=RooRealVar("rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,40,10,80);
   
            elif label_tstring.Contains("900") and (not label_tstring.Contains("1900") ):
                rrv_mean_voig=RooRealVar("rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,900,800,1000);
                rrv_width_voig=RooRealVar("rrv_width_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2.5,0,10);
                rrv_width_voig.setConstant(kTRUE);
                rrv_sigma_voig=RooRealVar("rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,40,10,90);
   
            elif label_tstring.Contains("1000"):
                rrv_mean_voig=RooRealVar("rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1000,900,1100);#
                rrv_width_voig=RooRealVar("rrv_width_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2.5,0,10);
                rrv_width_voig.setConstant(kTRUE);
                rrv_sigma_voig=RooRealVar("rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,40,10,80);

            if label_tstring.Contains("1100"):
                rrv_mean_voig=RooRealVar("rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1100,1000,1200);
                rrv_width_voig=RooRealVar("rrv_width_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3,0,10);
                rrv_width_voig.setConstant(kTRUE);
                rrv_sigma_voig=RooRealVar("rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,40,10,100);
   
            elif label_tstring.Contains("1200"):
                rrv_mean_voig=RooRealVar("rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1200,1100,1300);
                rrv_width_voig=RooRealVar("rrv_width_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3,0,10);
                rrv_width_voig.setConstant(kTRUE);
                rrv_sigma_voig=RooRealVar("rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_voig"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,40,10,100);

            model_pdf = RooVoigtian("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_mean_voig,rrv_width_voig,rrv_sigma_voig);

        ## BW for the W mass peak 
        if in_model_name == "BW":            
            print "########### BW Pdf for mj fit ############"
            rrv_mean_BW=RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel,"rrv_mean_BW"+label+"_"+self.categoryLabel,84,78, 88);
            rrv_width_BW=RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel,"rrv_width_BW"+label+"_"+self.categoryLabel,20,1,40);
            model_pdf = RooBreitWigner("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_mean_BW,rrv_width_BW);

        ## BW relativistic for Higgs mass covoluted with CB 
        if in_model_name == "BWRUNxGausErf":

            print "########### BWRUNxGausErf Pdf for Higgs lvj ############"
            if label=="_ggH600_signalregion" or label=="_ggH600_lowersideband":
                rrv_mean_BWRUN = RooRealVar("rrv_mean_BWRUN"+label+"_"+self.categoryLabel,"rrv_mean_BWRUN"+label+"_"+self.categoryLabel,600,550,650);
                rrv_width_BWRUN= RooRealVar("rrv_width_BWRUN"+label+"_"+self.categoryLabel,"rrv_width_BWRUN"+label+"_"+self.categoryLabel,60,50,70);

            if label=="_ggH700_signalregion" or label=="_ggH700_lowersideband":
                rrv_mean_BWRUN = RooRealVar("rrv_mean_BWRUN"+label+"_"+self.categoryLabel,"rrv_mean_BWRUN"+label+"_"+self.categoryLabel,700,650,750);
                rrv_width_BWRUN= RooRealVar("rrv_width_BWRUN"+label+"_"+self.categoryLabel,"rrv_width_BWRUN"+label+"_"+self.categoryLabel,100,80,120);

            if label=="_ggH800_signalregion" or label=="_ggH800_lowersideband":
                rrv_mean_BWRUN = RooRealVar("rrv_mean_BWRUN"+label+"_"+self.categoryLabel,"rrv_mean_BWRUN"+label+"_"+self.categoryLabel,800,750,850);
                rrv_width_BWRUN= RooRealVar("rrv_width_BWRUN"+label+"_"+self.categoryLabel,"rrv_width_BWRUN"+label+"_"+self.categoryLabel,150,100,180);

            if label=="_ggH900_signalregion" or label=="_ggH900_lowersideband":
                rrv_mean_BWRUN = RooRealVar("rrv_mean_BWRUN"+label+"_"+self.categoryLabel,"rrv_mean_BWRUN"+label+"_"+self.categoryLabel,900,850,990);
                rrv_width_BWRUN= RooRealVar("rrv_width_BWRUN"+label+"_"+self.categoryLabel,"rrv_width_BWRUN"+label+"_"+self.categoryLabel,200,100,260);

            if label=="_ggH1000_signalregion" or label=="_ggH1000_lowersideband":
                rrv_mean_BWRUN = RooRealVar("rrv_mean_BWRUN"+label+"_"+self.categoryLabel,"rrv_mean_BWRUN"+label+"_"+self.categoryLabel,1000,950,1050);
                rrv_width_BWRUN= RooRealVar("rrv_width_BWRUN"+label+"_"+self.categoryLabel,"rrv_width_BWRUN"+label+"_"+self.categoryLabel,200,100,370);

            bwrun = RooBWRunPdf("bwrun"+label+"_"+self.categoryLabel+mass_spectrum,"bwrun"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x, rrv_mean_BWRUN, rrv_width_BWRUN);

            rrv_mean_cb  = RooRealVar("rrv_mean_cb"+label+"_"+self.categoryLabel,"rrv_mean_cb"+label+"_"+self.categoryLabel,0);
            rrv_sigma_cb = RooRealVar("rrv_sigma_cb"+label+"_"+self.categoryLabel,"rrv_sigma_cb"+label+"_"+self.categoryLabel,50,10,300);
            cbshape      = RooGaussian("cbshape"+label+"_"+self.categoryLabel,"cbshape"+label+"_"+self.categoryLabel, rrv_x,rrv_mean_cb,rrv_sigma_cb);
            fft          = RooFFTConvPdf("fft"+label+"_"+self.categoryLabel+mass_spectrum,"fft"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x, bwrun, cbshape);

            rrv_offset_erf = RooRealVar("rrv_offset_erf"+label+"_"+self.categoryLabel,"rrv_offset_erf"+label+"_"+self.categoryLabel,450)#,350,550);
            rrv_width_erf = RooRealVar("rrv_width_erf"+label+"_"+self.categoryLabel,"rrv_width_erf"+label+"_"+self.categoryLabel,50)#,10,250);
            erf = RooGenericPdf("erf"+label+"_"+self.categoryLabel+mass_spectrum,"erf"+label+"_"+self.categoryLabel+mass_spectrum, "(1.+TMath::Erf((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset_erf.GetName(), rrv_width_erf.GetName()), RooArgList(rrv_x,rrv_offset_erf,rrv_width_erf) )

            model_pdf = RooProdPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, fft, erf );

        ##  Voig for W mass peak
        if in_model_name == "2Voig":

            print "########### Double Voigtian for mj fit ############"
            rrv_mean_voig    = RooRealVar("rrv_mean_voig"+label+"_"+self.categoryLabel,"rrv_mean_voig"+label+"_"+self.categoryLabel,84,78,88);#W mass 80.385
            rrv_shift_2Voig  = RooRealVar("rrv_shift_2Voig"+label+"_"+self.categoryLabel,"rrv_shift_2Voig"+label+"_"+self.categoryLabel,10.8026)# Z mass: 91.1876; shift=91.1876-80.385=10.8026
            rrv_mean_shifted = RooFormulaVar("rrv_mean_voig2"+label+"_"+self.categoryLabel,"@0+@1",RooArgList(rrv_mean_voig,rrv_shift_2Voig));

            rrv_width_voig = RooRealVar("rrv_width_voig"+label+"_"+self.categoryLabel,"rrv_width_voig"+label+"_"+self.categoryLabel,16.,6,26);
            rrv_sigma_voig = RooRealVar("rrv_sigma_voig"+label+"_"+self.categoryLabel,"rrv_sigma_voig"+label+"_"+self.categoryLabel,5.,0.,10.);

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.categoryLabel,"rrv_frac"+label+"_"+self.categoryLabel,0.8,0.5,1.);

            model_voig1 = RooVoigtian("model_voig1"+label+"_"+self.categoryLabel+mass_spectrum,"model_voig1"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_mean_voig,rrv_width_voig,rrv_sigma_voig);

            model_voig2 = RooVoigtian("model_voig2"+label+"_"+self.categoryLabel+mass_spectrum,"model_voig2"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_mean_shifted,rrv_width_voig,rrv_sigma_voig);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, RooArgList(model_voig1,model_voig2), RooArgList(rrv_frac));

        ## Gaus for the W peak
        if in_model_name == "Gaus":
            print "########### Gaus for W peak  ############"
            rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,84,78,88);
            rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.categoryLabel,"rrv_sigma_gaus"+label+"_"+self.categoryLabel,7,1,15);
            model_pdf = RooGaussian("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

        ## Gaus for the higgs lineshape
        if in_model_name == "Gaus_v1":

            print "########### Gaus for Higgs mlvj ############"
            if label=="_ggH600_signalregion" or label=="_ggH600_lowersideband":
                rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,580,550,620);
                rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.categoryLabel,"rrv_sigma_gaus"+label+"_"+self.categoryLabel,65,40,80);

            if label=="_ggH700_signalregion" or label=="_ggH700_lowersideband":
                rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,700,650,750);
                rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.categoryLabel,"rrv_sigma_gaus"+label+"_"+self.categoryLabel,100,40,150);

            if label=="_ggH800_signalregion" or label=="_ggH800_lowersideband":

                rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,800,750,850);
                rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.categoryLabel,"rrv_sigma_gaus"+label+"_"+self.categoryLabel,130,120,140);

            if label=="_ggH900_signalregion" or label=="_ggH900_lowersideband":
                rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,900,850,900);
                rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.categoryLabel,"rrv_sigma_gaus"+label+"_"+self.categoryLabel,160,140,180);

            if label=="_ggH1000_signalregion" or label=="_ggH1000_lowersideband":
                rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,920,900,1000);
                rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.categoryLabel,"rrv_sigma_gaus"+label+"_"+self.categoryLabel,200,100,300);

            model_pdf = RooGaussian("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

        if in_model_name == "BifurGaus_v1":

            print "########### BifurGaus for Higgs mlvj ############"
            if label=="_ggH600_signalregion":
                rrv_mean_gaus   = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,600,550,650);
                rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.categoryLabel,"rrv_sigma1_gaus"+label+"_"+self.categoryLabel,67,40,80);
                rrv_sigma2_gaus = RooRealVar("rrv_sigma2_gaus"+label+"_"+self.categoryLabel,"rrv_sigma2_gaus"+label+"_"+self.categoryLabel,67,40,80);

            if label=="_ggH700_signalregion":
                rrv_mean_gaus   = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,700,650,750);
                rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.categoryLabel,"rrv_sigma1_gaus"+label+"_"+self.categoryLabel,100,40,150);
                rrv_sigma2_gaus = RooRealVar("rrv_sigma2_gaus"+label+"_"+self.categoryLabel,"rrv_sigma2_gaus"+label+"_"+self.categoryLabel,100,40,150);

            if label=="_ggH800_signalregion":
                rrv_mean_gaus   = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,800,750,850);
                rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.categoryLabel,"rrv_sigma1_gaus"+label+"_"+self.categoryLabel,130,120,140);
                rrv_sigma2_gaus = RooRealVar("rrv_sigma2_gaus"+label+"_"+self.categoryLabel,"rrv_sigma2_gaus"+label+"_"+self.categoryLabel,130,120,140);

            if label=="_ggH900_signalregion":
                rrv_mean_gaus   = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,900,850,900);
                rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.categoryLabel,"rrv_sigma1_gaus"+label+"_"+self.categoryLabel,160,140,180);
                rrv_sigma2_gaus = RooRealVar("rrv_sigma2_gaus"+label+"_"+self.categoryLabel,"rrv_sigma2_gaus"+label+"_"+self.categoryLabel,160,140,180);

            if label=="_ggH1000_signalregion":
                rrv_mean_gaus   = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,920,900,1000);
                rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.categoryLabel,"rrv_sigma1_gaus"+label+"_"+self.categoryLabel,200,100,300);
                rrv_sigma2_gaus = RooRealVar("rrv_sigma2_gaus"+label+"_"+self.categoryLabel,"rrv_sigma2_gaus"+label+"_"+self.categoryLabel,200,100,300);

            model_pdf = RooBifurGauss("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_mean_gaus,rrv_sigma1_gaus,rrv_sigma2_gaus);

        ## Crystal Ball for the W mass peak
        if in_model_name == "CB":
            print "########### Cystal Ball for mj fit ############"
            rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,84,78,88);
            rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,7,4,10);
            rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.categoryLabel,"rrv_alpha_CB"+label+"_"+self.categoryLabel,-2,-4,-0.5);
            rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.categoryLabel,"rrv_n_CB"+label+"_"+self.categoryLabel,2,0.,4);
            model_pdf = RooCBShape("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_mean_CB,rrv_sigma_CB,rrv_alpha_CB,rrv_n_CB);

        ## Sum of two CB 
        if in_model_name == "SCB_v1":
            print "########### Cystal Ball + Crystall Ball ############"
            rrv_mean_SCB   = RooRealVar("rrv_mean_SCB"+label+"_"+self.categoryLabel,"rrv_mean_SCB"+label+"_"+self.categoryLabel,800,550,1000);
            rrv_sigma_SCB  = RooRealVar("rrv_sigma_SCB"+label+"_"+self.categoryLabel,"rrv_sigma_SCB"+label+"_"+self.categoryLabel,70,40,300);
            rrv_alpha1_SCB = RooRealVar("rrv_alpha1_SCB"+label+"_"+self.categoryLabel,"rrv_alpha1_SCB"+label+"_"+self.categoryLabel,-2,-4,-0.5);
            rrv_alpha2_SCB = RooRealVar("rrv_alpha2_SCB"+label+"_"+self.categoryLabel,"rrv_alpha2_SCB"+label+"_"+self.categoryLabel,2,0.5,4);
            rrv_n1_SCB     = RooRealVar("rrv_n1_SCB"+label+"_"+self.categoryLabel,"rrv_n1_SCB"+label+"_"+self.categoryLabel,2,0.,4);
            rrv_n2_SCB     = RooRealVar("rrv_n2_SCB"+label+"_"+self.categoryLabel,"rrv_n2_SCB"+label+"_"+self.categoryLabel,2,0.,4);
            frac           = RooRealVar("rrv_frac_SSCB"+label+"_"+self.categoryLabel,"rrv_frac_SSCB"+label+"_"+self.categoryLabel,0.5)
            scb1 = RooCBShape("model_pdf_scb1"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf_scb1"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_mean_SCB,rrv_sigma_SCB,rrv_alpha1_SCB,rrv_n1_SCB);
            scb2 = RooCBShape("model_pdf_scb2"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf_scb2"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_mean_SCB,rrv_sigma_SCB,rrv_alpha2_SCB,rrv_n2_SCB);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(scb1,scb2),RooArgList(frac))

        ## Double Gaus for Bulk GR signal --> narrow width
        if in_model_name == "2Gaus_sig":

            print "########### Double Gauss for Bulk GR ############"
            label_tstring=TString(label);

            if label_tstring.Contains("600") and (not label_tstring.Contains("1600") ):
                rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 600, 500, 700);
            elif label_tstring.Contains("700") and (not label_tstring.Contains("1700") ):
                rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 700, 600, 800);
            elif label_tstring.Contains("800") and (not label_tstring.Contains("1800") ):
                rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 800, 700, 900);
            elif label_tstring.Contains("900") and (not label_tstring.Contains("1900") ):
                rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 900, 800, 1000);
            elif label_tstring.Contains("1000"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 1000, 900, 1100);
            elif label_tstring.Contains("1100"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 1100, 1000, 1200);
            elif label_tstring.Contains("1200"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 1200, 1100, 1300);
            elif label_tstring.Contains("1300"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 1300, 1200, 1400);
            elif label_tstring.Contains("1400"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 1400, 1300, 1500);
            elif label_tstring.Contains("1500"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 1500, 1400, 1600);
            elif label_tstring.Contains("1600"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 1600, 1500, 1700);
            elif label_tstring.Contains("1700"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 1700, 1600, 1800);
            elif label_tstring.Contains("1800"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 1800, 1700, 1900);
            elif label_tstring.Contains("1900"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 1900, 1800, 2000);
            elif label_tstring.Contains("2000"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 2000, 1900, 2100);
            elif label_tstring.Contains("2100"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 2100, 2000, 2200);
            elif label_tstring.Contains("2200"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 2200, 2100, 2300);
            elif label_tstring.Contains("2300"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 2300, 2200, 2400);
            elif label_tstring.Contains("2400"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 2400, 2300, 2500);
            elif label_tstring.Contains("2500"): rrv_mean1_gaus=RooRealVar("rrv_mean_1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel, 2500, 2400, 2600);

            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.categoryLabel,"rrv_sigma1_gaus"+label+"_"+self.categoryLabel,50,20,120);
            gaus1 = RooGaussian("gaus1"+label+"_"+self.categoryLabel,"gaus1"+label+"_"+self.categoryLabel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.categoryLabel,"rrv_deltamean_gaus"+label+"_"+self.categoryLabel,0,-50,50);
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.categoryLabel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.categoryLabel,"rrv_scalesigma_gaus"+label+"_"+self.categoryLabel,1,0.,10.);
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.categoryLabel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.categoryLabel,"gaus2"+label+"_"+self.categoryLabel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.categoryLabel,"rrv_frac"+label+"_"+self.categoryLabel,0.5,0.,1.);

            model_pdf =RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac),1)

        ## Crystal  ball shape for Bulk GR samples and higgs 
        if in_model_name == "CB_v1":
            print "########### Crystal Ball for Higgs and  Bulk GR  mlvj ############"
            label_tstring=TString(label);

            if label_tstring.Contains("H600"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,600,580,620);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,67,40,80);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.categoryLabel,"rrv_alpha_CB"+label+"_"+self.categoryLabel,-1,-4,0);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.categoryLabel,"rrv_n_CB"+label+"_"+self.categoryLabel,20.,10,80 );
            elif label_tstring.Contains("H700"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,700,650,750);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,100,40,150);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.categoryLabel,"rrv_alpha_CB"+label+"_"+self.categoryLabel,-1,-3,-0.1);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.categoryLabel,"rrv_n_CB"+label+"_"+self.categoryLabel,20.,10,40);
            elif label_tstring.Contains("ggH800"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,780,700,850);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,140,120,160);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.categoryLabel,"rrv_alpha_CB"+label+"_"+self.categoryLabel,-1,-4,0);
                rrv_n_CB=RooRealVar("rrv_n_CB"+label+"_"+self.categoryLabel,"rrv_n_CB"+label+"_"+self.categoryLabel,5 , 2, 7);
            elif label_tstring.Contains("vbfH800"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,800,750,850);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,140,120,160);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.categoryLabel,"rrv_alpha_CB"+label+"_"+self.categoryLabel,-1,-4,0);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.categoryLabel,"rrv_n_CB"+label+"_"+self.categoryLabel,5 , 2, 7);
            elif label_tstring.Contains("ggH900"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,880,820,950);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,170,140,200);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.categoryLabel,"rrv_alpha_CB"+label+"_"+self.categoryLabel,1,0,4);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.categoryLabel,"rrv_n_CB"+label+"_"+self.categoryLabel, 2., 0.5,5);
            elif label_tstring.Contains("vbfH900"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,900,880,920);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,170,140,200);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.categoryLabel,"rrv_alpha_CB"+label+"_"+self.categoryLabel,1);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.categoryLabel,"rrv_n_CB"+label+"_"+self.categoryLabel, 2., 0.5,5);
            elif label_tstring.Contains("ggH1000"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,920,800,1150);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,200,100,300);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.categoryLabel,"rrv_alpha_CB"+label+"_"+self.categoryLabel,1,0.1,3);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.categoryLabel,"rrv_n_CB"+label+"_"+self.categoryLabel,2.,0.5,4);
            elif label_tstring.Contains("vbfH1000"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,1000,980,1150);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,200,100,300);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.categoryLabel,"rrv_alpha_CB"+label+"_"+self.categoryLabel,0.72);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.categoryLabel,"rrv_n_CB"+label+"_"+self.categoryLabel,2.,0.5,4);
            else:
                if label_tstring.Contains("M600") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M600_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 600, 550, 650);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 30,10 ,80);

                elif label_tstring.Contains("M700") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M700_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 700, 600, 800);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 30,10 ,80);

                elif label_tstring.Contains("M800") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M800_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 800, 600, 800);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 40,10 ,90);

                elif label_tstring.Contains("M900") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M900_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 900, 600, 800);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 40,10 ,90);

                elif label_tstring.Contains("M1000") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1000_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1000, 900,1100);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("M1100") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1100_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1100,1000,1200);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("M1200") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1200_W") :
                    rrv_mean_CB = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1200,1100,1300);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("M1300") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1300_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1300,1200,1400);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,65,45,120);

                elif label_tstring.Contains("M1400") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1400_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1400,1300,1500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,70,45,130);

                elif label_tstring.Contains("M1500") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1500_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1500,1400,1600);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,75,50,145);

                elif label_tstring.Contains("M1600") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1600_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1600,1500,1700);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,80,55,160);

                elif label_tstring.Contains("M1700") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1700_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1700,1600,1800);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,85,60,175);

                elif label_tstring.Contains("M1800") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1800_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1800,1700,1900);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,90,65,190);

                elif label_tstring.Contains("M1900") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1900_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1900,1800,2000);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,95,70,200);

                elif label_tstring.Contains("M2000") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2000_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2000,1900,2100);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,100,75,210);

                elif label_tstring.Contains("M2100") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2100_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100,2000,2200);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,105,80,225);

                elif label_tstring.Contains("M2200") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2200_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2200,2100,2300);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,110,85,245);

                elif label_tstring.Contains("M2300") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2300_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2300,2200,2400);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,115,90,250);

                elif label_tstring.Contains("M2400") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2400_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2400,2300,2500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,120,90,275);

                elif label_tstring.Contains("M2500") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2500_W") :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2500,2400,2600);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,135,90,285);

                elif label_tstring.Contains("M1000_W150") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1000,500,1500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 150,50 ,500);

                elif label_tstring.Contains("M1000_W300") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1000,500,1500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 300,50 ,800);

                elif label_tstring.Contains("M1000_W50") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1000,500,1500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 50,10 ,200);

                elif label_tstring.Contains("M1500_W75") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1500,1000,2000);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 75,50 ,250);

                elif label_tstring.Contains("M1500_W225") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1500,1000,2000);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 225,150 ,450);

                elif label_tstring.Contains("M1500_W450") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1500,1000,2000);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 450,400 ,700);

                elif label_tstring.Contains("M2100_W105") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100,1500,2500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 105,90 ,300);

                elif label_tstring.Contains("M2100_W315") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100,1500,2500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 315,250 ,600);

                elif label_tstring.Contains("M2100_W630") and label_tstring.Contains("BulkG_WW"):
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100,1500,2500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 630,550 ,900);

                else :
                    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,700,550,2500);
                    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 50,20 ,120);

            rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,4,1,5);
            rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,10,40);

            model_pdf = RooCBShape("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_mean_CB,rrv_sigma_CB,rrv_alpha_CB,rrv_n_CB);

        ## Crystal  ball shape for Bulk GR samples and higgs 
        if in_model_name == "BWCB":

            print "########### Crystal Ball x Breit Wigner for Bulk Graviton width ############"
            label_tstring=TString(label);

            rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 1000);
            rrv_width_BW = RooRealVar("rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 50);
            rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,0.,0.,50.);
            rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,50,0,200);

            #if label_tstring.Contains("M1000_W50") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 1000);
            #    rrv_width_BW = RooRealVar("rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 50);
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,0.,0.,50.);
            #    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,50,0,200);

            #elif label_tstring.Contains("M1000_W150") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 1000);
            #    rrv_width_BW = RooRealVar("rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 150);
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,0.,0.,50.);
            #    rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,50,0,200);

            #elif label_tstring.Contains("M1000_W300") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 1000);
            #        rrv_width_BW = RooRealVar("rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 300);
            #        rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,0.,0.,50.);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,50,0,200);

            #elif label_tstring.Contains("M1500_W75") and label_tstring.Contains("BulkG_WW"): 
            #    rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1500);
            #        rrv_sigma_BW = RooRealVar("rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 175,50 ,300);
            #        rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,0.,0.,50.);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,100,0,300);

            #elif label_tstring.Contains("M1500_W225") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1500,1000,2000);
            #        rrv_sigma_BW = RooRealVar("rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 225,150 ,450);
            #        rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,0.,0.,80.);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,100,0,250);

            #elif label_tstring.Contains("M1500_W450") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1500,1000,2000);
            #        rrv_sigma_BW = RooRealVar("rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 450,150 ,450);
            #        rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,0.,0.,80.);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,100,0,250);

            #elif label_tstring.Contains("M2100_W105") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100,1000,2000);
            #        rrv_sigma_BW = RooRealVar("rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 105,150 ,450);
            #        rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,0.,0.,80.);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,150,0,250);

            #elif label_tstring.Contains("M2100_W315") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100,1000,2000);
            #        rrv_sigma_BW = RooRealVar("rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 315,150 ,450);
            #        rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,0.,0.,80.);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,150,0,250);

            #elif label_tstring.Contains("M2100_W630") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100,1000,2000);
            #        rrv_sigma_BW = RooRealVar("rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 630,150 ,450);
            #        rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,0.,0.,80.);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,150,0,250);

            rrv_mean_BW.setConstant(kTRUE);
            rrv_width_BW.setConstant(kTRUE);

            rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.categoryLabel,"rrv_alpha_CB"+label+"_"+self.categoryLabel,2,0,4);
            rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.categoryLabel,"rrv_n_CB"+label+"_"+self.categoryLabel,1.,0.,4.);

            bw      = RooBreitWigner("bw"+label+"_"+self.categoryLabel,"bw"+label+"_"+self.categoryLabel, rrv_x,rrv_mean_BW,rrv_width_BW);
            cbshape = RooCBShape("cbshape"+label+"_"+self.categoryLabel,"cbshape"+label+"_"+self.categoryLabel, rrv_x,rrv_mean_CB,rrv_sigma_CB,rrv_alpha_CB,rrv_n_CB);

            model_pdf = RooFFTConvPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x, cbshape, bw);

        if in_model_name == "ArgusBW_v1":

            label_tstring=TString(label);
            if label_tstring.Contains("ggH1000"):
                rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel,"rrv_width_BW"+label+"_"+self.categoryLabel,100,50,600);
                rrv_m0_Argus = RooRealVar("rrv_m0_Argus"+label+"_"+self.categoryLabel,"rrv_m0_Argus"+label+"_"+self.categoryLabel, 950 );
                rrv_c_Argus  = RooRealVar("rrv_c_Argus"+label+"_"+self.categoryLabel,"rrv_c_Argus"+label+"_"+self.categoryLabel,-1,-2,-1e-1);
                rrv_frac     = RooRealVar("rrv_frac"+label+"_"+self.categoryLabel,"rrv_frac"+label+"_"+self.categoryLabel,0.5,0.0,1.);
            else:
                rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel,"rrv_width_BW"+label+"_"+self.categoryLabel,200,50,400);
                rrv_m0_Argus = RooRealVar("rrv_m0_Argus"+label+"_"+self.categoryLabel,"rrv_m0_Argus"+label+"_"+self.categoryLabel,1000);
                rrv_c_Argus  = RooRealVar("rrv_c_Argus"+label+"_"+self.categoryLabel,"rrv_c_Argus"+label+"_"+self.categoryLabel,-1,-2,0.1);
                rrv_frac     = RooRealVar("rrv_frac"+label+"_"+self.categoryLabel,"rrv_frac"+label+"_"+self.categoryLabel,0.5,0.0,1.);

            bw    = RooBreitWigner("bw"+label+"_"+self.categoryLabel,"bw"+label+"_"+self.categoryLabel, rrv_x,rrv_m0_Argus,rrv_width_BW);
            argus = RooArgusBG("argus"+label+"_"+self.categoryLabel,"argus"+label+"_"+self.categoryLabel, rrv_x, rrv_m0_Argus,rrv_c_Argus);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, RooArgList(bw,argus), RooArgList(rrv_frac));

        if in_model_name == "CBBW": # FFT: BreitWigner*CBShape
            print "########### Crystal Ball x Breit Wigner for W mass peak ############"
            rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel,"rrv_mean_CB"+label+"_"+self.categoryLabel,84.0,78,88);
            rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel,"rrv_sigma_CB"+label+"_"+self.categoryLabel,7,4,10);
            rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.categoryLabel,"rrv_alpha_CB"+label+"_"+self.categoryLabel,-2,-4,-1);
            rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.categoryLabel,"rrv_n_CB"+label+"_"+self.categoryLabel,0.5,0.,2);
            rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel,"rrv_mean_BW"+label+"_"+self.categoryLabel,0);
            rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel,"rrv_width_BW"+label+"_"+self.categoryLabel,10,5,20);
            cbshape      = RooCBShape("cbshape"+label+"_"+self.categoryLabel,"cbshape"+label+"_"+self.categoryLabel, rrv_x,rrv_mean_CB,rrv_sigma_CB,rrv_alpha_CB,rrv_n_CB);
            bw           = RooBreitWigner("bw"+label+"_"+self.categoryLabel,"bw"+label+"_"+self.categoryLabel, rrv_x,rrv_mean_BW,rrv_width_BW);
            model_pdf    = RooFFTConvPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x, cbshape, bw);

        if in_model_name == "LDGaus": # FFT: Landau*Gaus
            print "########### Landau x Breit Wigner for W mass peak ############"
            rrv_mean_landau  = RooRealVar("rrv_mean_landau"+label+"_"+self.categoryLabel,"rrv_mean_landau"+label+"_"+self.categoryLabel,84.0,78,88);
            rrv_sigma_landau = RooRealVar("rrv_sigma_landau"+label+"_"+self.categoryLabel,"rrv_sigma_landau"+label+"_"+self.categoryLabel,7,4,10);
            rrv_mean_gaus    = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,0);
            rrv_sigma_gaus   = RooRealVar("rrv_sigma_gaus"+label+"_"+self.categoryLabel,"rrv_sigma_gaus"+label+"_"+self.categoryLabel,16,10,20);
            landau           = RooLandau("landau"+label+"_"+self.categoryLabel,"landau"+label+"_"+self.categoryLabel, rrv_x,rrv_mean_landau,rrv_sigma_landau);
            gaus             = RooBreitWigner("gaus"+label+"_"+self.categoryLabel,"gaus"+label+"_"+self.categoryLabel, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);
            model_pdf        = RooFFTConvPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x, landau, gaus);

        ## Crystal  ball shape for Bulk GR samples and higgs 
        if in_model_name == "DoubleCB_v1":
            label_tstring=TString(label);
            print "########### Double CB for Bulk graviton mlvj ############"

            #if label_tstring.Contains("M600") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M600_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 600, 550, 650);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 30,10 ,80);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.,0.5,6.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3,0.5,6.);

            #elif label_tstring.Contains("M700") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M700_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 700, 600, 800);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 30,10 ,80);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.,0.5,6.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3,0.5,6.);

            #elif label_tstring.Contains("M800") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M800_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,820,790,880);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,50,40,70);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 15.,5.,25.);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.64,1.,1.9);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,15.,5.,25.);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,1.,1.9);


            #elif label_tstring.Contains("M900") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M900_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,920,850,950);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,59,45,70);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 25.,2,45);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,25.,0.1,45);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.25,0.5,3.);

            #elif label_tstring.Contains("M1000") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1000_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1020,970,1070);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,55,40,65);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,45);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.4,0.5,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.,0.5,3.5);

            #elif label_tstring.Contains("M1100") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1100_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1120,1080,1150);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,65,55,75);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,25);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,25);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);

            #elif label_tstring.Contains("M1200") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1200_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1220,1200,1250);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,65,55,75);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,30);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,5.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,30);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,5.);

            #elif label_tstring.Contains("M1300") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1300_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1320,1300,1350);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,70,60,75);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.3,0.5,3.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.3,0.5,3.);


            #elif label_tstring.Contains("M1400") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1400_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1420,1400,1440);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,77,65,85);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.5);

            #elif label_tstring.Contains("M1500") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1500_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1515,1500,1530);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,81,71,91);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 15.,0.01,25);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,15.,0.01,25);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.5);

            #elif label_tstring.Contains("M1600") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1600_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1620,1600,1640);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,81,70,90);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);

            #elif label_tstring.Contains("M1700") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1700_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1720,1700,1740);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,90,75,96);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);

            #elif label_tstring.Contains("M1800") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1800_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1820,1800,1840);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,90,75,100);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);

            #elif label_tstring.Contains("M1900") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M1900_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1920,1900,1940);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,95,80,115);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);

            #elif label_tstring.Contains("M2000") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2000_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2020,2000,2040);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,100,80,115);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);

            #elif label_tstring.Contains("M2100") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2100_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2120,2100,2140);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,105,85,115);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);

            #elif label_tstring.Contains("M2200") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2200_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2220,2200,2250);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,115,75,140);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);

            #elif label_tstring.Contains("M2300") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2300_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2320,2300,2340);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,115,95,120);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 15.,0.2,30);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,15.,0.2,20);

            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);

            #elif label_tstring.Contains("M2400") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2400_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2420,2400,2440);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,115,100,125);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);

            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);

            #elif label_tstring.Contains("M2500") and label_tstring.Contains("BulkG_WW") and not label_tstring.Contains("M2500_W") :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2520,2500,2540);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,125,90,145);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.);

            #elif label_tstring.Contains("M1000_W150") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1020,970,1070);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,150,130,175);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,45);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.4,0.5,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.,0.5,3.5);

            #elif label_tstring.Contains("M1000_W300") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1020,970,1070);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 300,50 ,800);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,45);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.4,0.5,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.,0.5,3.5);

            #elif label_tstring.Contains("M1000_W50") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1020,970,1070);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,50,25,1000);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,45);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.4,0.5,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.,0.5,3.5);

            #elif label_tstring.Contains("M1500_W75") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1500,1000,2000);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 75,50 ,250);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.,0.5,6.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3,0.5,6.);

            #elif label_tstring.Contains("M1500_W225") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1500,1000,2000);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 225,150 ,450);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.,0.5,6.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3,0.5,6.);

            #elif label_tstring.Contains("M1500_W450") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1500,1000,2000);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 450,400 ,700);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.,0.5,6.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3,0.5,6.);

            #elif label_tstring.Contains("M2100_W105") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100,1500,2500);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 105,90 ,300);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.,0.5,6.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3,0.5,6.);

            #elif label_tstring.Contains("M2100_W315") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100,1500,2500);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 315,250 ,600);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.,0.5,6.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3,0.5,6.);

            #elif label_tstring.Contains("M2100_W630") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100,2000,2200);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 630,500, 900);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.,0.5,6.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3,0.5,6.);

            #else :
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,700,550,2500);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 50,20 ,120);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.,0.5,6.);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3,0.5,6.);

            rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,700,550,2500);
            rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 50,20 ,120);
            rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,35);
            rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.,0.5,6.);
            rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3,0.5,6.);


            rrv_mean_scale_p1 = RooRealVar("CMS_sig_p1_jes","CMS_sig_p1_jes",0);
            rrv_mean_scale_p1.setConstant(kTRUE);
            if self.categoryLabel == "mu" :             
                rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_m","CMS_sig_p1_scale_m",0);
                rrv_mean_scale_p2.setConstant(kTRUE);
            elif self.categoryLabel == "el" :
                rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_e","CMS_sig_p1_scale_e",0);
                rrv_mean_scale_p2.setConstant(kTRUE);
            elif self.categoryLabel == "em":
                rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_em","CMS_sig_p1_scale_em",0);
                rrv_mean_scale_p2.setConstant(kTRUE);


            rrv_mean_scale_X1 = RooRealVar("rrv_mean_shift_scale_lep"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_shift_scale_lep"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,float(self.mean_signal_uncertainty_lep_scale));
            rrv_mean_scale_X1.setConstant(kTRUE);
            rrv_mean_scale_X2 = RooRealVar("rrv_mean_shift_scale_jes"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_shift_scale_jes"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,float(self.mean_signal_uncertainty_jet_scale));
            rrv_mean_scale_X2.setConstant(kTRUE);

            rrv_total_mean_CB = RooFormulaVar("rrv_total_mean_CB"+label+"_"+self.categoryLabel,"@0*(1+@1*@2)*(1+@3*@4)", RooArgList(rrv_mean_CB,rrv_mean_scale_p1,rrv_mean_scale_X1,rrv_mean_scale_p2,rrv_mean_scale_X2));

            if self.categoryLabel == "mu":
                rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_m","CMS_sig_p2_scale_m",0);
                rrv_sigma_scale_p1.setConstant(kTRUE);
            elif self.categoryLabel == "el":
                rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_e","CMS_sig_p2_scale_e",0);
                rrv_sigma_scale_p1.setConstant(kTRUE);
            elif self.categoryLabel == "em":
                rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_em","CMS_sig_p2_scale_em",0);
                rrv_sigma_scale_p1.setConstant(kTRUE);

            rrv_sigma_scale_p2 = RooRealVar("CMS_sig_p2_jer","CMS_sig_p2_jer",0);
            rrv_sigma_scale_p3 = RooRealVar("CMS_sig_p2_jes","CMS_sig_p2_jes",0);
            rrv_sigma_scale_p2.setConstant(kTRUE);
            rrv_sigma_scale_p3.setConstant(kTRUE);

            rrv_mean_sigma_X1 = RooRealVar("rrv_sigma_shift_lep_scale"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_shift_scale"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,float(self.sigma_signal_uncertainty_lep_scale));
            rrv_mean_sigma_X2 = RooRealVar("rrv_sigma_shift_jes"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_shift_scale"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,float(self.sigma_signal_uncertainty_jet_scale));
            rrv_mean_sigma_X3 = RooRealVar("rrv_sigma_shift_res"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_shift_res"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,float(self.sigma_signal_uncertainty_jet_res));
            rrv_mean_sigma_X1.setConstant(kTRUE);
            rrv_mean_sigma_X2.setConstant(kTRUE);
            rrv_mean_sigma_X3.setConstant(kTRUE);

            rrv_total_sigma_CB = RooFormulaVar("rrv_total_sigma_CB"+label+"_"+self.categoryLabel,"@0*(1+@1*@2)*(1+@3*@4)*(1+@5*@6)", RooArgList(rrv_sigma_CB,rrv_sigma_scale_p1,rrv_mean_sigma_X1,rrv_sigma_scale_p2,rrv_mean_sigma_X2,rrv_sigma_scale_p3,rrv_mean_sigma_X3));        

            model_pdf = ROOT.RooDoubleCrystalBall("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_total_mean_CB,rrv_total_sigma_CB,rrv_alpha1_CB,rrv_n1_CB,rrv_alpha2_CB,rrv_n2_CB);

        ## Crystal  ball shape for Bulk GR samples and higgs 
        if in_model_name == "BWDoubleCB":

            label_tstring=TString(label);
            print "########### Double CB x BW for Bulk graviton width ############"

            #if label_tstring.Contains("M1000_W50") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB   = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,0,-100,100);
            #        rrv_sigma_CB  = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,55,0,200);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,10.,0.01,45);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.4,0.2,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.,0.2,3.5);
            #        rrv_mean_BW   = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1000);
            #        rrv_width_BW  = RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,50);

            #elif label_tstring.Contains("M1000_W150") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,0,-100,100);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,55,0,200);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,10.,0.01,45);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.4,0.5,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.,0.5,3.5);
            #        rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1000);
            #        rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,150);

            #elif label_tstring.Contains("M1000_W300") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,0,-100,100);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,55,0,200);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,10.,0.01,45);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.4,0.5,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,35);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.,0.5,3.5);
            #        rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1000);
            #        rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,300);

            #elif label_tstring.Contains("M1500_W75") and label_tstring.Contains("BulkG_WW"): 
            #    rrv_mean_CB   = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,0,-100,100);
            #        rrv_sigma_CB  = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,75,0,200);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,10.,0.01,45);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.4,0.2,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,45);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.,0.2,3.5);
            #        rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1500);
            #        rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,75);

            #elif label_tstring.Contains("M1500_W225") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB   = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,0,-100,100);
            #        rrv_sigma_CB  = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,75,0,200);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,10.,0.01,45);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.4,0.2,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,45);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.,0.2,3.5);
            #        rrv_mean_BW   = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1500);
            #        rrv_width_BW  = RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,225);

            #elif label_tstring.Contains("M1500_W450") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB   = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,0,-100,100);
            #        rrv_sigma_CB  = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,75,0,200);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,10.,0.01,45);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.4,0.2,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,45);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.,0.2,3.5);
            #        rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1500);
            #        rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,450);

            #elif label_tstring.Contains("M2100_W105") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,0.,-100,100);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,90,20,250);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01, 45);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,45);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.5);
            #        rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100);
            #        rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,105);

            #elif label_tstring.Contains("M2100_W315") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,0.,-100,100);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,90,20,250);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 10.,0.01,45);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,45);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,1.5,0.5,3.5);
            #        rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100);
            #        rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,315);

            #elif label_tstring.Contains("M2100_W630") and label_tstring.Contains("BulkG_WW"):
            #    rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,0.,-100,100);
            #        rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,90,20,250);
            #        rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 20.,0.01,105);
            #        rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.5,0.5,50.5);
            #        rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,105);
            #        rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.5,0.5,50.5);
            #        rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100);
            #        rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,630);
            rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,0.,-100,100);
            rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,90,20,250);
            rrv_n1_CB     = RooRealVar("rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label, 20.,0.01,105);
            rrv_alpha2_CB = RooRealVar("rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.5,0.5,50.5);
            rrv_n2_CB     = RooRealVar("rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_n2_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,20.,0.01,105);
            rrv_alpha1_CB = RooRealVar("rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_alpha1_CB"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,3.5,0.5,50.5);
            rrv_mean_BW  = RooRealVar("rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,2100);
            rrv_width_BW = RooRealVar("rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_width_BW"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,630);


            ### fix the Breit-Wigner core to the generated one  
            rrv_mean_BW.setConstant(kTRUE);
            rrv_width_BW.setConstant(kTRUE);                    
            bw           = RooBreitWigner("bw"+label+"_"+self.categoryLabel,"bw"+label+"_"+self.categoryLabel, rrv_x,rrv_mean_BW,rrv_width_BW);

            ### Double Crystall ball term --> add parameters in order to do systematic on the signal shape inside combiner
            rrv_mean_scale_p1 = RooRealVar("CMS_sig_p1_jes","CMS_sig_p1_jes",0);  ## jes effect on the mean
            rrv_mean_scale_p1.setConstant(kTRUE);

            if self.categoryLabel == "mu" :  ### lep scale effect on the mean       
                rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_m","CMS_sig_p1_scale_m",0);
                rrv_mean_scale_p2.setConstant(kTRUE);

            elif self.categoryLabel == "el" :
                rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_e","CMS_sig_p1_scale_e",0);
                rrv_mean_scale_p2.setConstant(kTRUE);

            elif self.categoryLabel == "em":
                rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_em","CMS_sig_p1_scale_em",0);
                rrv_mean_scale_p2.setConstant(kTRUE);

            ## set the uncertainty value in other two independent variables 
            rrv_mean_scale_X1 = RooRealVar("rrv_mean_shift_scale_lep"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_shift_scale_lep"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,float(self.mean_signal_uncertainty_lep_scale));
            rrv_mean_scale_X1.setConstant(kTRUE);

            rrv_mean_scale_X2 = RooRealVar("rrv_mean_shift_scale_jes"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_mean_shift_scale_jes"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,float(self.mean_signal_uncertainty_jet_scale));
            rrv_mean_scale_X2.setConstant(kTRUE);

            ### total mean
            rrv_total_mean_CB = RooFormulaVar("rrv_total_mean_CB"+label+"_"+self.categoryLabel,"@0*(1+@1*@2)*(1+@3*@4)", RooArgList(rrv_mean_CB,rrv_mean_scale_p1,rrv_mean_scale_X1,rrv_mean_scale_p2,rrv_mean_scale_X2));

            ### lepton scale effect on the resolution 
            if self.categoryLabel == "mu":
                rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_m","CMS_sig_p2_scale_m",0);
                rrv_sigma_scale_p1.setConstant(kTRUE);
            elif self.categoryLabel == "el":
                rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_e","CMS_sig_p2_scale_e",0);
                rrv_sigma_scale_p1.setConstant(kTRUE);
            elif self.categoryLabel == "em":
                rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_em","CMS_sig_p2_scale_em",0);
                rrv_sigma_scale_p1.setConstant(kTRUE);

            ### jes and jer effect on the resolution             
            rrv_sigma_scale_p2 = RooRealVar("CMS_sig_p2_jer","CMS_sig_p2_jer",0);
            rrv_sigma_scale_p3 = RooRealVar("CMS_sig_p2_jes","CMS_sig_p2_jes",0);

            rrv_sigma_scale_p2.setConstant(kTRUE);
            rrv_sigma_scale_p3.setConstant(kTRUE);

            rrv_mean_sigma_X1 = RooRealVar("rrv_sigma_shift_lep_scale"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_shift_scale"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,float(self.sigma_signal_uncertainty_lep_scale));
            rrv_mean_sigma_X2 = RooRealVar("rrv_sigma_shift_jes"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_shift_scale"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,float(self.sigma_signal_uncertainty_jet_scale));
            rrv_mean_sigma_X3 = RooRealVar("rrv_sigma_shift_res"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,"rrv_sigma_shift_res"+label+"_"+self.categoryLabel+"_"+self.wtagger_label,float(self.sigma_signal_uncertainty_jet_res));

            rrv_mean_sigma_X1.setConstant(kTRUE);
            rrv_mean_sigma_X2.setConstant(kTRUE);
            rrv_mean_sigma_X3.setConstant(kTRUE);

            ### total resolution 
            rrv_total_sigma_CB = RooFormulaVar("rrv_total_sigma_CB"+label+"_"+self.categoryLabel,"@0*(1+@1*@2)*(1+@3*@4)*(1+@5*@6)", RooArgList(rrv_sigma_CB,rrv_sigma_scale_p1,rrv_mean_sigma_X1,rrv_sigma_scale_p2,rrv_mean_sigma_X2,rrv_sigma_scale_p3,rrv_mean_sigma_X3));        

            cbshape = ROOT.RooDoubleCrystalBall("DoubleCB"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_total_mean_CB,rrv_total_sigma_CB,rrv_alpha1_CB,rrv_n1_CB,rrv_alpha2_CB,rrv_n2_CB)

            ### numerical convolution via FFT
            model_pdf = RooFFTConvPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,bw,cbshape);
            model_pdf.setBufferFraction(1.0)

        ## ExpN pdf for W+jets bkg fit
        if in_model_name == "ExpN":

            print "########### ExpN funtion for W+jets mlvj ############"
            rrv_c_ExpN = RooRealVar("rrv_c_ExpN"+label+"_"+self.categoryLabel,"rrv_c_ExpN"+label+"_"+self.categoryLabel,-3e-3,-1e-1,-1e-5);
            if(ismc==1):
                rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.categoryLabel,"rrv_n_ExpN"+label+"_"+self.categoryLabel, 1e3, -1e2, 1e4);
            else :
                if self.categoryLabel == "el" :
                    rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.categoryLabel,"rrv_n_ExpN"+label+"_"+self.categoryLabel, 1e3, -1e2, 1e4);
                elif self.wtagger_label == "LP" :
                    rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.categoryLabel,"rrv_n_ExpN"+label+"_"+self.categoryLabel, 1e3, -1e2, 1e4);
                else:
                    rrv_n_ExpN = RooRealVar("rrv_n_ExpN"+label+"_"+self.categoryLabel,"rrv_n_ExpN"+label+"_"+self.categoryLabel, 5e2, 0, 1e3);

            model_pdf = ROOT.RooExpNPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_c_ExpN, rrv_n_ExpN);


        ## levelled exp for W+jets bkg fit
        if in_model_name == "ExpTail":
            print "########### ExpTai = levelled exp funtion for W+jets mlvj ############"
            label_tstring=TString(label);
            if self.wtagger_label == "LP":
                rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.categoryLabel,"rrv_s_ExpTail"+label+"_"+self.categoryLabel, 250,-1.e6,1e6);
                rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.categoryLabel,"rrv_a_ExpTail"+label+"_"+self.categoryLabel, 1e-1,-1.e-2,1e6);
            else:
                if self.categoryLabel == "el" :
                    if ismc == 1 and label_tstring.Contains("lowersideband"):
                        rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.categoryLabel,"rrv_s_ExpTail"+label+"_"+self.categoryLabel, 139,0.,355);
                        rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.categoryLabel,"rrv_a_ExpTail"+label+"_"+self.categoryLabel, 2e-2,-1.e-2,5.5e-2);                     
                    elif ismc == 1 and label_tstring.Contains("signalregion"):
                        rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.categoryLabel,"rrv_s_ExpTail"+label+"_"+self.categoryLabel, 162,18,395);
                        rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.categoryLabel,"rrv_a_ExpTail"+label+"_"+self.categoryLabel, 1.6e-2,-1.e-2,5.5e-2);
                    elif ismc == 0 :  
                        rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.categoryLabel,"rrv_s_ExpTail"+label+"_"+self.categoryLabel, 161,70,240);
                        rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.categoryLabel,"rrv_a_ExpTail"+label+"_"+self.categoryLabel, 8e-3,-1e-2,1.3e-1);
   
                if self.categoryLabel == "mu" or self.categoryLabel == "em":
                    if ismc == 1 and label_tstring.Contains("lowersideband"):
                        rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.categoryLabel,"rrv_s_ExpTail"+label+"_"+self.categoryLabel, 99,10,255);
                        rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.categoryLabel,"rrv_a_ExpTail"+label+"_"+self.categoryLabel, 3e-2,-1e-2,7.5e-2);                        
                    elif ismc == 1 and label_tstring.Contains("signalregion"):
                        rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.categoryLabel,"rrv_s_ExpTail"+label+"_"+self.categoryLabel, 110,20,242);
                        rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.categoryLabel,"rrv_a_ExpTail"+label+"_"+self.categoryLabel, 2.9e-2,-1e-2,7.5e-2);
                    elif ismc == 0 :  
                        rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.categoryLabel,"rrv_s_ExpTail"+label+"_"+self.categoryLabel, 161,40,280);
                        rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.categoryLabel,"rrv_a_ExpTail"+label+"_"+self.categoryLabel, 8e-3,-1e-2,1.3e-1);    

            model_pdf     = ROOT.RooExpTailPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_s_ExpTail, rrv_a_ExpTail);

        ## sum of two exponential 
        if in_model_name == "2Exp":
            print "########### 2Exp = levelled exp funtion for W+jets mlvj ############"
            rrv_c0_2Exp   = RooRealVar("rrv_c0_2Exp"+label+"_"+self.categoryLabel,"rrv_c0_2Exp"+label+"_"+self.categoryLabel, -5e-3, -8e-3,-4e-3);
            rrv_c1_2Exp   = RooRealVar("rrv_c1_2Exp"+label+"_"+self.categoryLabel,"rrv_c1_2Exp"+label+"_"+self.categoryLabel, -1e-3, -4e-3,-1e-4);
            rrv_frac_2Exp = RooRealVar("rrv_frac_2Exp"+label+"_"+self.categoryLabel,"rrv_frac_2Exp"+label+"_"+self.categoryLabel, 0., 0., 1e-2);
            model_pdf  = ROOT.Roo2ExpPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_c0_2Exp,rrv_c1_2Exp,rrv_frac_2Exp);

        ## sum of two exponential 
        if in_model_name == "Exp" or in_model_name == "Exp_sr":
            print "########### Exp = levelled exp funtion for W+jets mlvj ############"
            rrv_c_Exp = RooRealVar("rrv_c_Exp"+label+"_"+self.categoryLabel,"rrv_c_Exp"+label+"_"+self.categoryLabel,-0.05,-0.1,0.);
            model_pdf = ROOT.RooExponential("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_c_Exp);

        ## Erf times for mj spectrum
        if in_model_name == "ErfExp" :
            print "########### Erf*Exp for mj fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.categoryLabel,"rrv_c_ErfExp"+label+"_"+self.categoryLabel,-0.05,-0.1,-1e-4);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfExp"+label+"_"+self.categoryLabel,60.,30.,120);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.categoryLabel,"rrv_width_ErfExp"+label+"_"+self.categoryLabel,30.,10, 60.);
            model_pdf         = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        ## different initial values -> for mlvj
        if in_model_name == "ErfExp_v1" :
            print "########### Erf*Exp for mlvj fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.categoryLabel,"rrv_c_ErfExp"+label+"_"+self.categoryLabel,-0.006,-0.1,0.);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfExp"+label+"_"+self.categoryLabel,450.,400.,550.);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.categoryLabel,"rrv_width_ErfExp"+label+"_"+self.categoryLabel,70.,10,100.);
            model_pdf         = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        ## different initial values -> for mlvj
        if in_model_name == "ErfExp_v2" : 
            print "########### Erf*Exp for mlvj fit  ############"
            rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label+"_"+self.categoryLabel,"rrv_c_ErfExp"+label+"_"+self.categoryLabel,-0.005,-0.1,0.);
            rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfExp"+label+"_"+self.categoryLabel,450.,400.,500.);
            rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label+"_"+self.categoryLabel,"rrv_width_ErfExp"+label+"_"+self.categoryLabel, 50.,10,100.);
            rrv_residue_ErfExp = RooRealVar("rrv_residue_ErfExp"+label+"_"+self.categoryLabel,"rrv_residue_ErfExp"+label+"_"+self.categoryLabel,0.,0.,1.);
            model_pdf = RooGenericPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, "(TMath::Exp(%s*%s) + %s)*(1.+TMath::Erf((%s-%s)/%s))/2. "%(rrv_c_ErfExp.GetName(),rrv_x.GetName(), rrv_residue_ErfExp.GetName(), rrv_x.GetName(),rrv_offset_ErfExp.GetName(), rrv_width_ErfExp.GetName()), RooArgList(rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp,rrv_residue_ErfExp) )

        ## different initial values -> for mlvj
        if in_model_name == "ErfExp_v3" : #different init-value and range
            print "########### Erf*Exp for mlvj fit  ############"
            rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label+"_"+self.categoryLabel,"rrv_c_ErfExp"+label+"_"+self.categoryLabel,-0.005,-0.1,0.);
            rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfExp"+label+"_"+self.categoryLabel,450.,400,500.);
            rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label+"_"+self.categoryLabel,"rrv_width_ErfExp"+label+"_"+self.categoryLabel, 50.,10,100.);
            rrv_residue_ErfExp = RooRealVar("rrv_residue_ErfExp"+label+"_"+self.categoryLabel,"rrv_residue_ErfExp"+label+"_"+self.categoryLabel,0.,0.,1.);
            rrv_high_ErfExp    = RooRealVar("rrv_high_ErfExp"+label+"_"+self.categoryLabel,"rrv_high_ErfExp"+label+"_"+self.categoryLabel,1.,0.,400);
            rrv_high_ErfExp.setConstant(kTRUE);
            model_pdf = RooGenericPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, "(TMath::Exp(%s*%s) + %s)* TMath::Power( ((1+TMath::Erf((%s-%s)/%s))/2.), %s )"%(rrv_c_ErfExp.GetName(),rrv_x.GetName(), rrv_residue_ErfExp.GetName(),rrv_x.GetName(),rrv_offset_ErfExp.GetName(), rrv_width_ErfExp.GetName(), rrv_high_ErfExp.GetName()), RooArgList(rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_high_ErfExp,rrv_width_ErfExp,rrv_residue_ErfExp) )

        ## Exp+Gaus or mj spectrum
        if in_model_name == "ExpGaus":
            print "########### Exp + Gaus for mj  fit  ############"
            rrv_c_Exp       = RooRealVar("rrv_c_Exp"+label+"_"+self.categoryLabel,"rrv_c_Exp"+label+"_"+self.categoryLabel,0.05,-0.2,0.2);
            exp             = ROOT.RooExponential("exp"+label+"_"+self.categoryLabel,"exp"+label+"_"+self.categoryLabel,rrv_x,rrv_c_Exp);

            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel,84,78,88);
            rrv_sigma1_gaus = RooRealVar("rrv_smgma1_gaus"+label+"_"+self.categoryLabel,"rrv_sigma1_gaus"+label+"_"+self.categoryLabel,7,4,10);
            rrv_high        = RooRealVar("rrv_high"+label+"_"+self.categoryLabel,"rrv_high"+label+"_"+self.categoryLabel,0.5,0.,1.);
            gaus            = RooGaussian("gaus"+label+"_"+self.categoryLabel,"gaus"+label+"_"+self.categoryLabel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            model_pdf       = RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(exp,gaus),RooArgList(rrv_high))

        ## Erf*Exp + Gaus for mj spectrum 
        if in_model_name == "ErfExpGaus":
            print "########### Erf*Exp + Gaus for mj  fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.categoryLabel,"rrv_c_ErfExp"+label+"_"+self.categoryLabel,-0.05,-0.4,0.);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfExp"+label+"_"+self.categoryLabel,100.,10.,300.);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.categoryLabel,"rrv_width_ErfExp"+label+"_"+self.categoryLabel,30.,10,100.);

            erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.categoryLabel,"erfExp"+label+"_"+self.categoryLabel,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,82,78,87);
            rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.categoryLabel,"rrv_sigma_gaus"+label+"_"+self.categoryLabel,7,4,10);
            gaus = RooGaussian("gaus"+label+"_"+self.categoryLabel,"gaus"+label+"_"+self.categoryLabel, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

            rrv_high  = RooRealVar("rrv_high"+label+"_"+self.categoryLabel,"rrv_high"+label+"_"+self.categoryLabel,0.7,0.,1.);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))

        ## Erf*Exp + Gaus for mj spectrum with offset == mean
        if in_model_name == "ErfExpGaus_sp":
            print "########### Erf*Exp + Gaus for mj  fit  ############"
            rrv_c_ErfExp     = RooRealVar("rrv_c_ErfExp"+label+"_"+self.categoryLabel,"rrv_c_ErfExp"+label+"_"+self.categoryLabel,-0.05,-0.2,0.);
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.categoryLabel,"rrv_width_ErfExp"+label+"_"+self.categoryLabel,30.,10,200.);
            erfExp           = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.categoryLabel,"erfExp"+label+"_"+self.categoryLabel,rrv_x,rrv_c_ErfExp,rrv_mean1_gaus,rrv_width_ErfExp);

            rrv_mean1_gaus   = RooRealVar("rrv_mean1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel,84,78,88);
            rrv_sigma1_gaus  = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.categoryLabel,"rrv_sigma1_gaus"+label+"_"+self.categoryLabel,7,4,10);
            gaus             = RooGaussian("gaus"+label+"_"+self.categoryLabel,"gaus"+label+"_"+self.categoryLabel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_high  = RooRealVar("rrv_high"+label+"_"+self.categoryLabel,"rrv_high"+label+"_"+self.categoryLabel,0.5,0.,1.);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))

        ## Erf*Exp+Gaus or mj spectrum
        if in_model_name == "ErfExpGaus_v0":
            print "########### Erf*Exp + Gaus for mj  fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.categoryLabel,"rrv_c_ErfExp"+label+"_"+self.categoryLabel,-0.05,-0.2,0.);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfExp"+label+"_"+self.categoryLabel,100.,10.,140.);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.categoryLabel,"rrv_width_ErfExp"+label+"_"+self.categoryLabel,30.,10,100.);
            erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.categoryLabel,"erfExp"+label+"_"+self.categoryLabel,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_mean_gaus     = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,84,78,88);
            rrv_sigma_gaus    = RooRealVar("rrv_sigma_gaus"+label+"_"+self.categoryLabel,"rrv_sigma_gaus"+label+"_"+self.categoryLabel,7,4,10);
            gaus              = RooGaussian("gaus"+label+"_"+self.categoryLabel,"gaus"+label+"_"+self.categoryLabel, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

            rrv_high   = RooRealVar("rrv_high"+label+"_"+self.categoryLabel,"rrv_high"+label+"_"+self.categoryLabel,0.7,0.,1.);
            model_pdf  = RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))

        ## Erf*Exp+Gaus or mj spectrum
        if in_model_name == "ErfExpGaus_v1":
            print "########### Erf*Exp + Gaus for mlvj fit  ############"
            rrv_c_ErfExp       = RooRealVar("rrv_c_ErfExp"+label+"_"+self.categoryLabel,"rrv_c_ErfExp"+label+"_"+self.categoryLabel,-0.007,-0.1,0.);
            rrv_offset_ErfExp  = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfExp"+label+"_"+self.categoryLabel,800.,10.,1400.);
            rrv_width_ErfExp   = RooRealVar("rrv_width_ErfExp"+label+"_"+self.categoryLabel,"rrv_width_ErfExp"+label+"_"+self.categoryLabel,24.,10,150.);
            erfExp             = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.categoryLabel,"erfExp"+label+"_"+self.categoryLabel,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_mean_gaus   = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,700,500,1200);
            rrv_sigma_gaus  = RooRealVar("rrv_sigma_gaus"+label+"_"+self.categoryLabel,"rrv_sigma_gaus"+label+"_"+self.categoryLabel,150,10,300);
            gaus            = RooGaussian("gaus"+label+"_"+self.categoryLabel,"gaus"+label+"_"+self.categoryLabel, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

            rrv_high  = RooRealVar("rrv_high"+label+"_"+self.categoryLabel,"rrv_high"+label+"_"+self.categoryLabel,0.1,0.,1.);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))

        ## Erf*Exp+Gaus or mj spectrum
        if in_model_name == "ErfExpGaus_sp_v1":
            print "########### Erf*Exp + Gaus for mlvj fit  ############"
            rrv_c_ErfExp     = RooRealVar("rrv_c_ErfExp"+label+"_"+self.categoryLabel,"rrv_c_ErfExp"+label+"_"+self.categoryLabel,-0.007,-0.1,0.);
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.categoryLabel,"rrv_width_ErfExp"+label+"_"+self.categoryLabel,24.,10,150.);
            rrv_mean_gaus    = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,900,860,1200);
            erfExp           = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.categoryLabel,"erfExp"+label+"_"+self.categoryLabel,rrv_x,rrv_c_ErfExp,rrv_mean_gaus,rrv_width_ErfExp);

            rrv_sigma_gaus   = RooRealVar("rrv_sigma_gaus"+label+"_"+self.categoryLabel,"rrv_sigma_gaus"+label+"_"+self.categoryLabel,150,10,300);
            gaus = RooGaussian("gaus"+label+"_"+self.categoryLabel,"gaus"+label+"_"+self.categoryLabel, rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

            rrv_high  = RooRealVar("rrv_high"+label+"_"+self.categoryLabel,"rrv_high"+label+"_"+self.categoryLabel,0.1,0.,1.);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))

        ## Erf*Exp+Gaus or mj spectrum
        if in_model_name == "ErfExpGaus_v2":
            print "########### Erf*Exp + Gaus for mj fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.categoryLabel,"rrv_c_ErfExp"+label+"_"+self.categoryLabel,-0.05,-10.,0.);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfExp"+label+"_"+self.categoryLabel,100.,10.,140.);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.categoryLabel,"rrv_width_ErfExp"+label+"_"+self.categoryLabel,30.,10,100.);
            rrv_mean_gaus     = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,84,78,88);
            rrv_sigma_gaus    = RooRealVar("rrv_sigma_gaus"+label+"_"+self.categoryLabel,"rrv_sigma_gaus"+label+"_"+self.categoryLabel,7,4,10);
            rrv_high          = RooRealVar("rrv_high"+label+"_"+self.categoryLabel,"rrv_high"+label+"_"+self.categoryLabel,200.,0.,1000.);
            model_pdf = ROOT.RooErfExp_Gaus_Pdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp,rrv_mean_gaus,rrv_sigma_gaus,rrv_high );

        ## Erf*Exp + 2Gaus  
        if in_model_name == "ErfExp2Gaus":
            print "########### Erf*Exp + 2Gaus for mj fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.categoryLabel,"rrv_c_ErfExp"+label+"_"+self.categoryLabel,-0.05,-0.2,0.);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfExp"+label+"_"+self.categoryLabel,100.,10.,140.);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.categoryLabel,"rrv_width_ErfExp"+label+"_"+self.categoryLabel,30.,10,100.);
            erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.categoryLabel,"erfExp"+label+"_"+self.categoryLabel,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_mean1_gaus   = RooRealVar("rrv_mean1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel,84,78,88);
            rrv_mean2_gaus   = RooRealVar("rrv_mean2_gaus"+label+"_"+self.categoryLabel,"rrv_mean2_gaus"+label+"_"+self.categoryLabel,180,170,190);
            rrv_sigma1_gaus  = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.categoryLabel,"rrv_sigma1_gaus"+label+"_"+self.categoryLabel,7,4,10);
            rrv_sigma2_gaus  = RooRealVar("rrv_sigma2_gaus"+label+"_"+self.categoryLabel,"rrv_sigma2_gaus"+label+"_"+self.categoryLabel,10,7,15);
            gaus1 = RooGaussian("gaus1"+label+"_"+self.categoryLabel,"gaus1"+label+"_"+self.categoryLabel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);
            gaus2 = RooGaussian("gaus2"+label+"_"+self.categoryLabel,"gaus2"+label+"_"+self.categoryLabel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_high1 = RooRealVar("rrv_high1"+label+"_"+self.categoryLabel,"rrv_high1"+label+"_"+self.categoryLabel,0.6,0.,1.);
            rrv_high2 = RooRealVar("rrv_high2"+label+"_"+self.categoryLabel,"rrv_high2"+label+"_"+self.categoryLabel,0.4,0.,1.);
            model_pdf =RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(erfExp,gaus1,gaus2),RooArgList(rrv_high1,rrv_high2))

        ## Gaus + Gaus for mj spectrum
        if in_model_name == "2Gaus":
            print "########### 2Gaus for mj fit  ############"
            mean1_tmp      = 8.3141e+01; mean1_tmp_err      = 1.63e-01;
            deltamean_tmp  = 6.9129e+00; deltamean_tmp_err  = 1.24e+00;
            sigma1_tmp     = 7.5145e+00; sigma1_tmp_err     = 1.99e-01;
            scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
            frac_tmp       = 6.7125e-01; frac_tmp_err       = 2.09e-02;

            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.categoryLabel,"rrv_sigma1_gaus"+label+"_"+self.categoryLabel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.categoryLabel,"gaus1"+label+"_"+self.categoryLabel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.categoryLabel,"rrv_deltamean_gaus"+label+"_"+self.categoryLabel,deltamean_tmp,deltamean_tmp-deltamean_tmp_err*4 ,deltamean_tmp+deltamean_tmp_err*4);
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.categoryLabel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.categoryLabel,"rrv_scalesigma_gaus"+label+"_"+self.categoryLabel,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*8, scalesigma_tmp+scalesigma_tmp_err*8);
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.categoryLabel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.categoryLabel,"gaus2"+label+"_"+self.categoryLabel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac  = RooRealVar("rrv_frac"+label+"_"+self.categoryLabel,"rrv_frac"+label+"_"+self.categoryLabel,frac_tmp, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac),1)

        ## 2Gaus+2Gaus for VV mj spectrum -> WZ and WW
        if in_model_name == "2_2Gaus":

            print "########### 2Gaus +2Gaus for mj fit  ############"
            mean1_tmp      = 8.3141e+01; mean1_tmp_err      = 1.63e-01;
            deltamean_tmp  = 6.9129e+00; deltamean_tmp_err  = 1.24e+00;
            sigma1_tmp     = 7.5145e+00; sigma1_tmp_err     = 1.99e-01;
            scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
            frac_tmp       = 6.7125e-01; frac_tmp_err       = 2.09e-02;

            rrv_shift = RooRealVar("rrv_shift"+label+"_"+self.categoryLabel,"rrv_shift"+label+"_"+self.categoryLabel,10.8026) # Z mass: 91.1876; shift=91.1876-80.385=10.8026

            rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.categoryLabel,"rrv_sigma1_gaus"+label+"_"+self.categoryLabel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.categoryLabel,"gaus1"+label+"_"+self.categoryLabel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.categoryLabel,"rrv_deltamean_gaus"+label+"_"+self.categoryLabel,0.,-8,10);
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.categoryLabel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.categoryLabel,"rrv_scalesigma_gaus"+label+"_"+self.categoryLabel,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*4, scalesigma_tmp+scalesigma_tmp_err*4);
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.categoryLabel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.categoryLabel,"gaus2"+label+"_"+self.categoryLabel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac1 = RooRealVar("rrv_frac1"+label+"_"+self.categoryLabel,"rrv_frac1"+label+"_"+self.categoryLabel,frac_tmp, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);
            gausguas_1 =RooAddPdf("gausguas_1"+label+"_"+self.categoryLabel+mass_spectrum,"gausguas_1"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac1),1)

            rrv_mean3_gaus = RooFormulaVar("rrv_mean3_gaus"+label+"_"+self.categoryLabel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_shift));
            rrv_mean4_gaus = RooFormulaVar("rrv_mean4_gaus"+label+"_"+self.categoryLabel,"@0+@1",RooArgList(rrv_mean2_gaus, rrv_shift));
            gaus3 = RooGaussian("gaus3"+label+"_"+self.categoryLabel,"gaus3"+label+"_"+self.categoryLabel, rrv_x,rrv_mean3_gaus,rrv_sigma1_gaus);
            gaus4 = RooGaussian("gaus4"+label+"_"+self.categoryLabel,"gaus4"+label+"_"+self.categoryLabel, rrv_x,rrv_mean4_gaus,rrv_sigma2_gaus);
            gausguas_2 = RooAddPdf("gausguas_2"+label+"_"+self.categoryLabel+mass_spectrum,"gausguas_2"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(gaus3,gaus4),RooArgList(rrv_frac1),1)

            rrv_frac  = RooRealVar("rrv_frac"+label+"_"+self.categoryLabel,"rrv_frac"+label+"_"+self.categoryLabel,0.74)#,0.5,1.0)
            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(gausguas_1,gausguas_2),RooArgList(rrv_frac),1)

        ## Erf*Exp + 2Gaus for mj spectrum
        if in_model_name == "2Gaus_ErfExp":

            print "########### 2Gaus + Erf*Exp for mj fit  ############"
            mean1_tmp      = 8.3141e+01; mean1_tmp_err      = 1.63e-01;
            deltamean_tmp  = 6.9129e+00; deltamean_tmp_err  = 1.24e+00;
            sigma1_tmp     = 7.5145e+00; sigma1_tmp_err     = 1.99e-01;
            scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
            frac_tmp       = 6.7125e-01; frac_tmp_err       = 2.09e-02;

            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.categoryLabel,"rrv_mean1_gaus"+label+"_"+self.categoryLabel,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.categoryLabel,"rrv_sigma1_gaus"+label+"_"+self.categoryLabel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.categoryLabel,"gaus1"+label+"_"+self.categoryLabel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.categoryLabel,"rrv_deltamean_gaus"+label+"_"+self.categoryLabel,deltamean_tmp)#, deltamean_tmp, deltamean_tmp);
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.categoryLabel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.categoryLabel,"rrv_scalesigma_gaus"+label+"_"+self.categoryLabel,scalesigma_tmp)#, scalesigma_tmp, scalesigma_tmp);
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.categoryLabel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.categoryLabel,"gaus2"+label+"_"+self.categoryLabel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac_2gaus = RooRealVar("rrv_frac_2gaus"+label+"_"+self.categoryLabel,"rrv_frac_2gaus"+label+"_"+self.categoryLabel,frac_tmp);#, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);

            c0_tmp     = -2.9893e-02 ; c0_tmp_err     = 6.83e-03;
            offset_tmp = 7.9350e+01  ; offset_tmp_err = 9.35e+00;
            width_tmp  = 3.3083e+01  ; width_tmp_err  = 2.97e+00;

            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.categoryLabel,"rrv_c_ErfExp"+label+"_"+self.categoryLabel,c0_tmp, c0_tmp-4e-2, c0_tmp+4e-2 );
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfExp"+label+"_"+self.categoryLabel, offset_tmp)#, offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.categoryLabel,"rrv_width_ErfExp"+label+"_"+self.categoryLabel, width_tmp, width_tmp-10, width_tmp+10);
            erfexp = ROOT.RooErfExpPdf("erfexp"+label+"_"+self.categoryLabel+mass_spectrum,"erfexp"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.categoryLabel,"rrv_frac"+label+"_"+self.categoryLabel, 0.5,0.,1.);
            model_pdf =RooAddPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,RooArgList(erfexp, gaus1,gaus2),RooArgList(rrv_frac, rrv_frac_2gaus),1)


        ## Erf*Exp+Voig+Gaus for mj spectrum 
        if in_model_name == "ErfExpVoigGaus":
            print "########### Erf*Exp + Voig + Gaus for mj fit  ############"
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.categoryLabel,"rrv_c_ErfExp"+label+"_"+self.categoryLabel,-0.1,-10.,0.);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfExp"+label+"_"+self.categoryLabel,100.,10.,140.);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.categoryLabel,"rrv_width_ErfExp"+label+"_"+self.categoryLabel,30.,10,100.);
            rrv_mean_voig     = RooRealVar("rrv_mean_voig"+label+"_"+self.categoryLabel,"rrv_mean_voig"+label+"_"+self.categoryLabel,84,78,88);
            rrv_width_voig    = RooRealVar("rrv_width_voig"+label+"_"+self.categoryLabel,"rrv_width_voig"+label+"_"+self.categoryLabel,7,1,20);
            rrv_sigma_voig    = RooRealVar("rrv_sigma_voig"+label+"_"+self.categoryLabel,"rrv_sigma_voig"+label+"_"+self.categoryLabel,5,1,100);
            rrv_high1         = RooRealVar("rrv_high1"+label+"_"+self.categoryLabel,"rrv_high1"+label+"_"+self.categoryLabel,1,0.,200.);
            rrv_mean_gaus     = RooRealVar("rrv_mean_gaus"+label+"_"+self.categoryLabel,"rrv_mean_gaus"+label+"_"+self.categoryLabel,174)#,160,187);
            rrv_sigma_gaus    = RooRealVar("rrv_sigma_gaus"+label+"_"+self.categoryLabel,"rrv_sigma_gaus"+label+"_"+self.categoryLabel,20)#,0.1,100);
            rrv_high2 = RooRealVar("rrv_high2"+label+"_"+self.categoryLabel,"rrv_high2"+label+"_"+self.categoryLabel,0.)#,0.,0.);
            model_pdf = ROOT.RooErfExp_Voig_Gaus_Pdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp,rrv_mean_voig,rrv_width_voig,rrv_sigma_voig,rrv_high1,rrv_mean_gaus,rrv_sigma_gaus,rrv_high2 );

        ## User1 function 
        if in_model_name == "User1":
            print "########### User 1 Pdf  for mlvj fit ############"
            rrv_p0     = RooRealVar("rrv_p0_User1"+label+"_"+self.categoryLabel,"rrv_p0_User1"+label+"_"+self.categoryLabel, 30, 10, 90);
            if self.wtagger_label=="HP":
                rrv_p1 = RooRealVar("rrv_p1_User1"+label+"_"+self.categoryLabel,"rrv_p1_User1"+label+"_"+self.categoryLabel, -4, -9, -2);
            else:
                rrv_p1 = RooRealVar("rrv_p1_User1"+label+"_"+self.categoryLabel,"rrv_p1_User1"+label+"_"+self.categoryLabel, -2, -4, 0.);
            model_pdf=RooUser1Pdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_p0,rrv_p1);

        ## QCD pdf  
        if in_model_name == "QCD":
            print "########### QCD Pdf  for mlvj fit ############"
            rrv_p0 = RooRealVar("rrv_p0_QCD"+label+"_"+self.categoryLabel,"rrv_p0_QCD"+label+"_"+self.categoryLabel, 0,-200,200);
            rrv_p1 = RooRealVar("rrv_p1_QCD"+label+"_"+self.categoryLabel,"rrv_p1_QCD"+label+"_"+self.categoryLabel, 0,-200,200);
            rrv_p2 = RooRealVar("rrv_p2_QCD"+label+"_"+self.categoryLabel,"rrv_p2_QCD"+label+"_"+self.categoryLabel, 0,-200,200);
            model_pdf = RooQCDPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_p0,rrv_p1,rrv_p2);

        if in_model_name == "QCD_v2":#can replace exp
            print "########### QCD Pdf  for mlvj fit ############"
            rrv_p0 = RooRealVar("rrv_p0_QCD"+label+"_"+self.categoryLabel,"rrv_p0_QCD"+label+"_"+self.categoryLabel, -15,-50,0);
            rrv_p1 = RooRealVar("rrv_p1_QCD"+label+"_"+self.categoryLabel,"rrv_p1_QCD"+label+"_"+self.categoryLabel, 20,0,250);
            rrv_p2 = RooRealVar("rrv_p2_QCD"+label+"_"+self.categoryLabel,"rrv_p2_QCD"+label+"_"+self.categoryLabel,0,-20,20);
            model_pdf = RooQCDPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_p0,rrv_p1,rrv_p2);

        ## For mlvj fit -> Pow function can replace exp
        if in_model_name == "Pow" or in_model_name == "Pow_sr" :
            print "########### Pow Pdf  for mlvj fit ############"
            rrv_c = RooRealVar("rrv_c_Pow"+label+"_"+self.categoryLabel,"rrv_c_Pow"+label+"_"+self.categoryLabel, -5, -20, 0);
            model_pdf = RooPowPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x, rrv_c );

        ## For mlvj fit -> Pow function can replace exp
        if in_model_name == "Pow2":
            print "########### Pow2 Pdf  for mlvj fit ############"
            rrv_c0 = RooRealVar("rrv_c0_Pow2"+label+"_"+self.categoryLabel,"rrv_c0_Pow2"+label+"_"+self.categoryLabel, 5, 0, 20);
            rrv_c1 = RooRealVar("rrv_c1_Pow2"+label+"_"+self.categoryLabel,"rrv_c1_Pow2"+label+"_"+self.categoryLabel, 0, -5 , 5);
            model_pdf = RooPow2Pdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x, rrv_c0, rrv_c1 );

        ## For mlvj fit ->Erf*Pow can replace Erf*Exp
        if in_model_name == "ErfPow_v1":
            print "########### Erf*Pow Pdf  for mlvj fit ############"
            rrv_c      = RooRealVar("rrv_c_ErfPow"+label+"_"+self.categoryLabel,"rrv_c_ErfPow"+label+"_"+self.categoryLabel, -5,-10,0);
            rrv_offset = RooRealVar("rrv_offset_ErfPow"+label+"_"+self.categoryLabel,"rrv_offset_ErfPow"+label+"_"+self.categoryLabel, 450,350,550);
            rrv_width  = RooRealVar("rrv_width_ErfPow"+label+"_"+self.categoryLabel,"rrv_width_ErfPow"+label+"_"+self.categoryLabel,50,20,90);
            model_pdf  = RooErfPowPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_c,rrv_offset,rrv_width);

        ## For mlvj fit ->Erf*Pow can replace Erf*Exp -> in the sideband
        if in_model_name == "ErfPow2_v1":
            print "########### Erf*Pow2 Pdf  for mlvj fit ############"
            rrv_c0 = RooRealVar("rrv_c0_ErfPow2"+label+"_"+self.categoryLabel,"rrv_c0_ErfPow2"+label+"_"+self.categoryLabel,14,1,30);
            rrv_c1 = RooRealVar("rrv_c1_ErfPow2"+label+"_"+self.categoryLabel,"rrv_c1_ErfPow2"+label+"_"+self.categoryLabel, 5,-5,10);
            rrv_offset = RooRealVar("rrv_offset_ErfPow2"+label+"_"+self.categoryLabel,"rrv_offset_ErfPow2"+label+"_"+self.categoryLabel, 450,400,520);
            rrv_width  = RooRealVar("rrv_width_ErfPow2"+label+"_"+self.categoryLabel,"rrv_width_ErfPow2"+label+"_"+self.categoryLabel,30,10,80);
            model_pdf  = RooErfPow2Pdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);

        ## For mlvj fit ->Erf*Pow can replace Erf*Exp for sr
        if in_model_name == "ErfPow2_v1_sr":
            print "########### Erf*Pow2 Pdf  for mlvj fit in the SR  ############"
            rrv_c0 = RooRealVar("rrv_c0_ErfPow2"+label+"_"+self.categoryLabel,"rrv_c0_ErfPow2"+label+"_"+self.categoryLabel, 4,2, 8);
            rrv_c1 = RooRealVar("rrv_c1_ErfPow2"+label+"_"+self.categoryLabel,"rrv_c1_ErfPow2"+label+"_"+self.categoryLabel, -0.5,-2,0);
            rrv_offset = RooRealVar("rrv_offset_ErfPow2"+label+"_"+self.categoryLabel,"rrv_offset_ErfPow2"+label+"_"+self.categoryLabel, 490,440,520);
            rrv_width  = RooRealVar("rrv_width_ErfPow2"+label+"_"+self.categoryLabel,"rrv_width_ErfPow2"+label+"_"+self.categoryLabel,50,30,80);
            model_pdf = RooErfPow2Pdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);

        ## For mlvj fit ->Erf*Pow*Exp can replace Erf*Exp 
        if in_model_name == "ErfPowExp_v1":
            print "########### Erf*Pow*Exp Pdf  for mlvj fit   ############"
            rrv_c0 = RooRealVar("rrv_c0_ErfPowExp"+label+"_"+self.categoryLabel,"rrv_c0_ErfPowExp"+label+"_"+self.categoryLabel,11,5,20);
            rrv_c1 = RooRealVar("rrv_c1_ErfPowExp"+label+"_"+self.categoryLabel,"rrv_c1_ErfPowExp"+label+"_"+self.categoryLabel, 0,-2,2);
            rrv_offset = RooRealVar("rrv_offset_ErfPowExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfPowExp"+label+"_"+self.categoryLabel, 470,420,520);
            rrv_width  = RooRealVar("rrv_width_ErfPowExp"+label+"_"+self.categoryLabel,"rrv_width_ErfPowExp"+label+"_"+self.categoryLabel,40,30,50);
            model_pdf  = RooErfPowExpPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);

        ## For mlvj fit ->Erf*Pow*Exp can replace Erf*Exp 
        if in_model_name == "ErfPowExp_v1_sr":
            print "########### Erf*Pow*Exp Pdf for mlvj fit in SR  ############"
            rrv_c0 = RooRealVar("rrv_c0_ErfPowExp"+label+"_"+self.categoryLabel,"rrv_c0_ErfPowExp"+label+"_"+self.categoryLabel,6,2,15);
            rrv_c1 = RooRealVar("rrv_c1_ErfPowExp"+label+"_"+self.categoryLabel,"rrv_c1_ErfPowExp"+label+"_"+self.categoryLabel, -1,-3,2);
            rrv_offset = RooRealVar("rrv_offset_ErfPowExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfPowExp"+label+"_"+self.categoryLabel, 490,440,520);
            rrv_width  = RooRealVar("rrv_width_ErfPowExp"+label+"_"+self.categoryLabel,"rrv_width_ErfPowExp"+label+"_"+self.categoryLabel,50,30,70);
            model_pdf=RooErfPowExpPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);

        ## For mlvj fit ->Erf*Pow*Exp can replace Erf*Exp 
        if in_model_name == "ErfPowExp_v1_0":#difference inital value
            print "########### Erf*Pow*Exp Pdf for mlvj fit in SR  ############"
            rrv_c0 = RooRealVar("rrv_c0_ErfPowExp"+label+"_"+self.categoryLabel,"rrv_c0_ErfPowExp"+label+"_"+self.categoryLabel,20,15,40);
            rrv_c1 = RooRealVar("rrv_c1_ErfPowExp"+label+"_"+self.categoryLabel,"rrv_c1_ErfPowExp"+label+"_"+self.categoryLabel, 1.6,0.5,5);
            rrv_offset = RooRealVar("rrv_offset_ErfPowExp"+label+"_"+self.categoryLabel,"rrv_offset_ErfPowExp"+label+"_"+self.categoryLabel, 470,420,520);
            rrv_width  = RooRealVar("rrv_width_ErfPowExp"+label+"_"+self.categoryLabel,"rrv_width_ErfPowExp"+label+"_"+self.categoryLabel,47,30,60);
            model_pdf  = RooErfPowExpPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);

        ## Keys 
        if in_model_name == "Keys":
            print "########### Keys PDF  ############"
            rdataset = self.workspace4fit_.data("rdataset4fit"+label+"_"+self.categoryLabel+"_mlvj");
            rdataset.Print();
            model_pdf = RooKeysPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, rrv_x, rdataset);


        ## Hist 
        if in_model_name == "Hist":
            print "########### Hist PDF  ############"
            rdataset = self.workspace4fit_.data("rdataset4fit"+label+"_"+self.categoryLabel+"_mlvj");
            rdataset.Print();
            rdatahist= rdataset.binnedClone("rdatahist4fit"+label+"_"+self.categoryLabel+"_mlvj","rdatahist4fit"+label+"_"+self.categoryLabel+"_mlvj");
            model_pdf = RooHistPdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum,"model_pdf"+label+"_"+self.categoryLabel+mass_spectrum, RooArgSet(rrv_x), rdatahist);

        ## return the pdf
        getattr(self.workspace4fit_,"import")(model_pdf)
        return self.workspace4fit_.pdf("model_pdf"+label+"_"+self.categoryLabel+mass_spectrum)
'''
'''
    ########### Gaussian contraint of a parameter of a pdf
    def addConstraint(self, rrv_x, x_mean, x_sigma, ConstraintsList):
        print "########### Add to Contraint List some parameters  ############"
     rrv_x_mean = RooRealVar(rrv_x.GetName()+"_mean",rrv_x.GetName()+"_mean",x_mean );
     rrv_x_sigma = RooRealVar(rrv_x.GetName()+"_sigma",rrv_x.GetName()+"_sigma",x_sigma );
     constrainpdf_x = RooGaussian("constrainpdf_"+rrv_x.GetName(),"constrainpdf_"+rrv_x.GetName(),rrv_x, rrv_x_mean, rrv_x_sigma);
     ## import in the workspace and save the name of constriant pdf
     getattr(self.workspace4fit_,"import")(constrainpdf_x)
     ConstraintsList.append(constrainpdf_x.GetName());

    ### get an mj model from the workspace givin the label
    def get_mj_Model(self,label):
        return self.workspace4fit_.pdf("model"+label+"_"+self.categoryLabel+"_mj")

    ### take the dataset, the model , the parameters in order to fix them as constant --> for extended pdf
    def get_General_mj_Model(self, label ):
        print "########### Fixing a general mj model  ############"
        rdataset_General_mj = self.workspace4fit_.data("rdataset%s_%s_mj"%(label,self.categoryLabel))
        model_General = self.get_mj_Model(label);
        rdataset_General_mj.Print();
        model_General.Print();
        ## get the parameters and cycle on them
        parameters_General = model_General.getParameters(rdataset_General_mj);
        par=parameters_General.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            if (TString(label).Contains("VV") or TString(label).Contains("SingleT") or TString(label).Contains("TTbar")):
                param.Print();
            param.setConstant(kTRUE);
            param=par.Next()
        ## return the pdf after having fixed the paramters
        return self.workspace4fit_.pdf("model%s_%s_mj"%(label,self.categoryLabel))

    ### fix only the ttbar component using the default label --> for extended pdf
    def get_TTbar_mj_Model(self,label="_TTbar"):
        print "########### Fixing only the TTbar mj Shape  ############"
        return self.get_General_mj_Model(label);

    ### fix only the stop component using the default label --> for extended pdf
    def get_SingleT_mj_Model(self,label="_SingleT"):
        print "########### Fixing only the Stop mj Shape  ############"
        return self.get_General_mj_Model(label);

    ### fix only the VV component using the default label --> for extended pdf
    def get_VV_mj_Model(self,label="_VV"):
        print "########### Fixing only the VV mj Shape  ############"
        return self.get_General_mj_Model(label);

    ### fix only the WJets model --> for extended pdf (just fix shape parameters of width, offset of ErfExp and p1 of User1 function
    def get_WJets_mj_Model(self,label):
        print "########### Fixing only the WJets mj Shape --> just the printed parameters  ############"
        rdataset_WJets_mj = self.workspace4fit_.data("rdataset%s_%s_mj"%(label,self.categoryLabel))
        model_WJets = self.get_mj_Model(label);
        rdataset_WJets_mj.Print();
        model_WJets.Print();
        parameters_WJets = model_WJets.getParameters(rdataset_WJets_mj);
        par=parameters_WJets.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            paraName=TString(param.GetName());
            if ( paraName.Contains("rrv_width_ErfExp_WJets") or paraName.Contains("rrv_offset_ErfExp_WJets") or paraName.Contains("rrv_p1_User1_WJets")):
                param.setConstant(kTRUE);
             param.Print();
         else:
             param.setConstant(0);
            param=par.Next()
        return self.workspace4fit_.pdf("model%s_%s_mj"%(label,self.categoryLabel))

    ### fix a given model taking the label, and the region --> for extended pdf --> all the parameter of the pdf + normalization
    def fix_Model(self, label, mlvj_region="_signalregion",mass_spectrum="_mlvj"):
        print "########### Fixing an Extended Pdf for mlvj  ############"        
        rdataset = self.workspace4fit_.data("rdataset%s%s_%s%s"%(label,mlvj_region,self.categoryLabel,mass_spectrum))
        model = self.get_mlvj_Model(label,mlvj_region);
        rdataset.Print();
        model.Print();
        parameters = model.getParameters(rdataset);
        par=parameters.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param=par.Next()

    ### fix a pdf in a different way --> for RooAbsPdf 
    def fix_Pdf(self,model_pdf,argset_notparameter):
        print "########### Fixing a RooAbsPdf for mlvj or mj  ############"        
        parameters = model_pdf.getParameters(argset_notparameter);
        par=parameters.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()

    ### print the parameters of a given pdf --> only non constant ones
    def ShowParam_Pdf(self,model_pdf,argset_notparameter):
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


    #### get a generic mlvj model from the workspace
    def get_mlvj_Model(self,label, mlvj_region):
        return self.workspace4fit_.pdf("model"+label+mlvj_region+"_"+self.categoryLabel+"_mlvj");

    #### get a general mlvj model and fiz the paramters --> for extended pdf
    def get_General_mlvj_Model(self, label, mlvj_region="_signalregion"):
        print "########### Fixing a general mlvj model  ############"
        rdataset_General_mlvj = self.workspace4fit_.data("rdataset%s%s_%s_mlvj"%(label, mlvj_region,self.categoryLabel))
        model_General = self.get_mlvj_Model(label,mlvj_region);
        rdataset_General_mlvj.Print();
        model_General.Print();
        parameters_General = model_General.getParameters(rdataset_General_mlvj);
        par=parameters_General.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()
        return self.get_mlvj_Model(label,mlvj_region);

    ###### get TTbar model mlvj in a region 
    def get_TTbar_mlvj_Model(self, mlvj_region="_signalregion"):
        print "########### Fixing TTbar mlvj model for the region",mlvj_region,"  ############"
        return self.get_General_mlvj_Model("_TTbar",mlvj_region);

    ###### get Single Top model mlvj in a region 
    def get_SingleT_mlvj_Model(self, mlvj_region="_signalregion"):
        print "########### Fixing Stop mlvj model for the region",mlvj_region,"  ############"
        return self.get_General_mlvj_Model("_SingleT",mlvj_region);

    ###### get Signal model mlvj in a region 
    def get_signal_mlvj_Model(self, mlvj_region="_signalregion"):
        print "########### Fixing signal mlvj model for the region",mlvj_region,"  ############"
        return self.get_General_mlvj_Model("_%s"%(self.allsignals),mlvj_region);

    ###### get VV mlvj in a region 
    def get_VV_mlvj_Model(self, mlvj_region="_signalregion"):
        print "########### Fixing VV mlvj for the region",mlvj_region,"  ############"
        return self.get_General_mlvj_Model("_VV",mlvj_region);

    ###### get W+jets mlvj in a region 
    def get_WJets_mlvj_Model(self, mlvj_region="_signalregion"):
        rdataset_WJets_mlvj = self.workspace4fit_.data("rdataset_WJets_%s_mlvj"%(mlvj_region))
        model_WJets = self.get_mlvj_Model("_WJets0",mlvj_region);
        print "######## get Wjet mlvj model for the region --> set constant just the normalization from mj fit",mlvj_region," ########";
        rdataset_WJets_mlvj.Print()
        model_WJets.Print()
        parameters_WJets = model_WJets.getParameters(rdataset_WJets_mlvj);
        par = parameters_WJets.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            paraName=TString(param.GetName());
            param.Print();
            if paraName.Contains("rrv_number_WJets"): ## set the correct normalization for W+jets if we are inside the signal region and fix it as constant
                if self.workspace4fit_.var("rrv_number_WJets_in_mj%s_from_fitting_%s"%(mlvj_region,self.categoryLabel)):
                    self.workspace4fit_.var("rrv_number_WJets_in_mj%s_from_fitting_%s"%(mlvj_region,self.categoryLabel)).Print()
                    param.setVal( self.workspace4fit_.var("rrv_number_WJets_in_mj%s_from_fitting_%s"%(mlvj_region,self.categoryLabel)).getVal() )
                if mlvj_region=="_signalregion": param.setConstant(kTRUE);
            param.Print();
            param=par.Next()
        return self.get_mlvj_Model("_WJets0",mlvj_region);


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
        print "######## change the dataset into a histogramm for mj distribution ########"        
        datahist=dataset.binnedClone(dataset.GetName()+"_binnedClone",dataset.GetName()+"_binnedClone")
        nbin=int( (x.getMax()-x.getMin())/self.obs0_variable_BinWidth);
        if label=="":
            return datahist.createHistogram("histo_%s"%(dataset.GetName()),x, RooFit.Binning( nbin ,x.getMin(),x.getMax()));
        else:
            return datahist.createHistogram("histo_"+label,x, RooFit.Binning( nbin,x.getMin(),x.getMax()));


'''
'''
  ### Method for a single MC fit of the mj spectra giving: file name, label, model name
    def fit_obs_variable_SingleChannel(self,in_file_name, label, in_model_name, additioninformation=""):

        print "############### Fit mj single MC sample",in_file_name," ",label,"  ",in_model_name," ##################"
        ## import variable and dataset
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");
        rdataset_mj = self.workspace4fit_.data("rdataset4fit"+label+"_"+self.categoryLabel+"_mj");
        rdataset_mj.Print();

        ## make the extended model
        model = self.make_Model(label,in_model_name);
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.Extended(kTRUE) );
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult.Print();

        ## Plot the result
        mplot = rrv_mass_j.frame(RooFit.Title(label+" fitted by "+in_model_name), RooFit.Bins(int(rrv_mass_j.getBins()/self.BinWidth_narrow_factor)) );
        rdataset_mj.plotOn( mplot, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        ## draw the error band for an extend pdf
        draw_error_band_extendPdf(rdataset_mj, model, rfresult,mplot,2,"L");
        ## re-draw the dataset
        rdataset_mj.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        ## draw the function
        model.plotOn( mplot );# remove RooFit.VLines() in order to get right pull in the 1st bin

        ## Get the pull
        mplot_pull = self.get_pull(rrv_mass_j, mplot); 
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);

        parameters_list = model.getParameters(rdataset_mj);
        self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s_%s_%s/m_j_fitting%s_wtaggercut%s/"%(self.additioninformation, self.categoryLabel,self.PS_model, self.wtagger_label, additioninformation, self.wtagger_label), label+in_file_name, in_model_name)

        #normalize the number of total events to lumi --> correct the number to scale to the lumi
        self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_mj").setVal( self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_mj").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_mj").setError(self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_mj").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).getVal() )

        if TString(label).Contains("ggH"):
            self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_mj").setVal( self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_mj").getVal() )
            self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_mj").setError(self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_mj").getError() )
        self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_mj").Print();

        ##### apply the correction of the mean and sigma from the ttbar control sample to the SingleT, TTbar and VV 
        par=parameters_list.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            if (TString(label).Contains("VV") or TString(label).Contains("SingleT") or TString(label).Contains("TTbar")):
                #param.Print();
                if TString(param.GetName()).Contains("rrv_mean1_gaus"):
                    param.setRange(param.getMin()+self.mean_shift, param.getMax()+self.mean_shift);
                    param.setVal(param.getVal()+self.mean_shift);
                if TString(param.GetName()).Contains("rrv_deltamean_gaus"):
                    param.setRange(param.getMin()-self.mean_shift, param.getMax()-self.mean_shift);
                    param.setVal(param.getVal()-self.mean_shift);
                if TString(param.GetName()).Contains("rrv_sigma1_gaus"):
                    param.setVal(param.getVal()*self.sigma_scale);
                    param.setRange(param.getMin()*self.sigma_scale, param.getMax()*self.sigma_scale);
                if TString(param.GetName()).Contains("rrv_scalesigma_gaus"):
                    param.setRange(param.getMin()/self.sigma_scale, param.getMax()/self.sigma_scale);
                    param.setVal(param.getVal()/self.sigma_scale);
            param=par.Next()

    ### Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
    def fit_limit_variable_SingleChannel(self,in_file_name, label, in_range, mlvj_model, deco=0, show_constant_parameter=0, logy=0, ismc=0):

        print "############### Fit mlvj single MC sample ",in_file_name," ",label,"  ",mlvj_model,"  ",in_range," ##################"
        ## import variable and dataset
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rdataset = self.workspace4fit_.data("rdataset4fit"+label+in_range+"_"+self.categoryLabel+"_mlvj");
        constrainslist =[];

        ## make the extended pdf model
        model = self.make_Model(label+in_range,mlvj_model,"_mlvj",constrainslist,ismc);

        ## make the fit
        model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE) );
        rfresult = model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult.Print();

        ## set the name of the result of the fit and put it in the workspace   
        rfresult.SetName("rfresult"+label+in_range+"_"+self.categoryLabel+"_mlvj")
        getattr(self.workspace4fit_,"import")(rfresult)

        ## plot the result
        mplot = rrv_mass_lvj.frame(RooFit.Title("M_{lvj"+in_range+"} fitted by "+mlvj_model), RooFit.Bins(int(rrv_mass_lvj.getBins()/self.BinWidth_narrow_factor)));
        rdataset.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        ## plot the error band but don't store the canvas (only plotted without -b option
        draw_error_band_extendPdf(rdataset, model, rfresult,mplot,2,"L")
        rdataset.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        model.plotOn( mplot )#, RooFit.VLines()); in order to have the right pull 

        ## get the pull 
        mplot_pull      = self.get_pull(rrv_mass_lvj,mplot);
        parameters_list = model.getParameters(rdataset);
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);

        self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s_%s_%s/m_lvj_fitting/"%(self.additioninformation, self.categoryLabel,self.PS_model,self.wtagger_label), in_file_name,"m_lvj"+in_range+mlvj_model, show_constant_parameter, logy);


        ## if the shape parameters has to be decorrelated
        if deco :
            print "################### Decorrelated mlvj single mc shape ################"
            model_pdf = self.workspace4fit_.pdf("model_pdf%s%s_%s_mlvj"%(label,in_range,self.categoryLabel)); ## take the pdf from the workspace
            model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) );
            rfresult_pdf = model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));
            rfresult_pdf.Print();

            ## temp workspace for the pdf diagonalizer
            wsfit_tmp = RooWorkspace("wsfit_tmp"+label+in_range+"_"+self.categoryLabel+"_mlvj");
            Deco      = PdfDiagonalizer("Deco"+label+in_range+"_"+self.categoryLabel+"_"+self.wtagger_label+"_mlvj",wsfit_tmp,rfresult_pdf); ## in order to have a good name 
            print "##################### diagonalize ";
            model_pdf_deco = Deco.diagonalize(model_pdf); ## diagonalize            
            print "##################### workspace for decorrelation ";
            wsfit_tmp.Print("v");
            print "##################### original  parameters ";
            model_pdf.getParameters(rdataset).Print("v");
            print "##################### original  decorrelated parameters ";
            model_pdf_deco.getParameters(rdataset).Print("v");
            print "##################### original  pdf ";
            model_pdf.Print();
            print "##################### decorrelated pdf ";
            model_pdf_deco.Print();

            ## import in the workspace and print the diagonalizerd pdf
            getattr(self.workspace4fit_,"import")(model_pdf_deco);

            ### define a frame for TTbar or other plots
            mplot_deco = rrv_mass_lvj.frame( RooFit.Bins(int(rrv_mass_lvj.getBins()/self.BinWidth_narrow_factor)));

            if label=="_TTbar" and in_range=="_signalregion":

                rdataset.plotOn(mplot_deco, RooFit.Name("Powheg Sample"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
                model_pdf_deco.plotOn(mplot_deco,RooFit.Name("TTbar_Powheg"),RooFit.LineColor(kBlack));

                mplot_deco.GetYaxis().SetRangeUser(1e-2,mplot_deco.GetMaximum()*1.2);

                rrv_number_dataset = RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset.sumEntries());
                rrv_number_dataset.setError(0.)
                draw_error_band(rdataset, model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,self.color_palet["Uncertainty"],"F"); ## draw the error band with the area
                self.workspace4fit_.var("rrv_number_TTbar_signalregion_%s_mlvj"%(self.categoryLabel)).Print();
            else:
                rdataset.plotOn(mplot_deco, RooFit.Name("Data"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
                model_pdf_deco.plotOn(mplot_deco,RooFit.Name(label),RooFit.LineColor(kBlack));

                mplot_deco.GetYaxis().SetRangeUser(1e-2,mplot_deco.GetMaximum()*1.2);

                rrv_number_dataset=RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset.sumEntries());
                rrv_number_dataset.setError(0.)
                draw_error_band(rdataset, model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,self.color_palet["Uncertainty"],"F"); ## don't store the number in the workspace

            self.plot_legend = self.legend4Plot(mplot_deco,0); ## add the plot_legend                
            mplot_deco.addObject(self.plot_legend);

            self.draw_canvas( mplot_deco, "plots_%s_%s_%s_%s/other/"%(self.additioninformation, self.categoryLabel,self.PS_model, self.wtagger_label), "m_lvj"+label+in_range+in_range+mlvj_model+"_deco",0,logy)

        ### Number of the event in the dataset and lumi scale factor --> set the proper number for bkg extraction or for signal region
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_mlvj").Print();
        self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).Print()
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_mlvj").setVal( self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_mlvj").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_mlvj").setError(self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_mlvj").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).getVal() )

        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_mlvj").Print();


    #### method to fit the WJets normalization inside the mj signal region -> and write the jets mass sys if available
    def fit_WJetsNorm(self, scaleJetMass = 0): # to get the normalization of WJets in signalregion

        print "############### Fit mj Normalization ##################"
        ## fit the two version of pdf for Wjets shape if available
        self.fit_WJetsNormalization_in_Mj_signalregion("_WJets0");
        self.fit_WJetsNormalization_in_Mj_signalregion("_WJets01");

        ## in case fit also the scaled jet mass distributions in order to have the jet mass scale sys included
        if scaleJetMass :
            self.fit_WJetsNormalization_in_Mj_signalregion("_WJets0_massup","massup");
         self.fit_WJetsNormalization_in_Mj_signalregion("_WJets0_massdn","massdn");
         self.fit_WJetsNormalization_in_Mj_signalregion("_WJets1");

        ## take the normalization numbers
        rrv_WJets0  = self.workspace4fit_.var("rrv_number_WJets0_in_mj_signalregion_from_fitting_%s"%(self.categoryLabel));
        rrv_WJets01 = self.workspace4fit_.var("rrv_number_WJets01_in_mj_signalregion_from_fitting_%s"%(self.categoryLabel));
        rrv_WJets0.Print();
        rrv_WJets01.Print();
        if scaleJetMass :
            rrv_WJets1 = self.workspace4fit_.var("rrv_number_WJets1_in_mj_signalregion_from_fitting_%s"%(self.categoryLabel));
         rrv_WJets1.Print();
         rrv_WJets0massup.Print();
         rrv_WJets0massdn.Print();

        ### total uncertainty combining the result with two different shapes
        total_uncertainty = TMath.Sqrt( TMath.Power(rrv_WJets0.getError(),2) + TMath.Power(rrv_WJets01.getVal()-rrv_WJets0.getVal(),2) );
        rrv_WJets0.setError(total_uncertainty);
        rrv_WJets0.Print();

        ##jet mass uncertainty on WJets normalization and the other bkg component
        if self.workspace4fit_.var("rrv_number_WJets0_massup_in_mj_signalregion_from_fitting_%s"%(self.categoryLabel)) and self.workspace4fit_.var("rrv_number_WJets0_massdn_in_mj_signalregion_from_fitting_%s"%(self.categoryLabel)):            
            rrv_WJets0massup = self.workspace4fit_.var("rrv_number_WJets0_massup_in_mj_signalregion_from_fitting_%s"%(self.categoryLabel));
          rrv_WJets0massdn = self.workspace4fit_.var("rrv_number_WJets0_massdn_in_mj_signalregion_from_fitting_%s"%(self.categoryLabel));
          self.WJets_normlization_uncertainty_from_jet_mass= ( TMath.Abs(rrv_WJets0massup.getVal()-rrv_WJets0.getVal())+TMath.Abs(rrv_WJets0massdn.getVal()-rrv_WJets0.getVal() ) )/2./rrv_WJets0.getVal();

        rrv_SingleT  = self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT__%s_mj"%(self.categoryLabel));

        if self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT_massup_%s_mj"%(self.categoryLabel)) and self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT_massdn_%s_mj"%(self.categoryLabel)) :
            rrv_SingleTmassup = self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT_massup_%s_mj"%(self.categoryLabel));
         rrv_SingleTmassdn = self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT_massdn_%s_mj"%(self.categoryLabel));
         self.SingleT_normlization_uncertainty_from_jet_mass=( TMath.Abs(rrv_SingleTmassup.getVal()-rrv_SingleT.getVal())+TMath.Abs(rrv_SingleTmassdn.getVal()-rrv_SingleT.getVal() ) )/2./rrv_SingleT.getVal();

        rrv_TTbar = self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar__%s_mj"%(self.categoryLabel));
        if self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar_massup_%s_mj"%(self.categoryLabel)) and self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar_massdn_%s_mj"%(self.categoryLabel)):
            rrv_TTbarmassup = self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar_massup_%s_mj"%(self.categoryLabel));
         rrv_TTbarmassdn = self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar_massdn_%s_mj"%(self.categoryLabel));
         self.TTbar_normlization_uncertainty_from_jet_mass=( TMath.Abs(rrv_TTbarmassup.getVal()-rrv_TTbar.getVal())+TMath.Abs(rrv_TTbarmassdn.getVal()-rrv_TTbar.getVal() ) )/2./rrv_TTbar.getVal();

        rrv_VV = self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_%s_mj"%(self.categoryLabel));
        if self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_massup_%s_mj"%(self.categoryLabel)) and self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_massdn_%s_mj"%(self.categoryLabel)):
            rrv_VVmassup = self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_massup_%s_mj"%(self.categoryLabel));
         rrv_VVmassdn = self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_massdn_%s_mj"%(self.categoryLabel));
         self.VV_normlization_uncertainty_from_jet_mass=( TMath.Abs(rrv_VVmassup.getVal()-rrv_VV.getVal())+TMath.Abs(rrv_VVmassdn.getVal()-rrv_VV.getVal() ) )/2./rrv_VV.getVal();

    #### make the mj sideband fit on data ti get the Wjets normaliztion 
    def fit_WJetsNormalization_in_Mj_signalregion(self,label,massscale=""): 

        print "############### Fit mj Normalization: ",label," ",massscale," ##################"
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
        ## get real data in mj distribution --> mass up and down have only an effect on Wjets shape -> effect on the normalization -> evaluated in the MC and fit data
        rdataset_data_mj=self.workspace4fit_.data("rdataset_data_%s_mj"%(self.categoryLabel))

        ### Fix TTbar, VV and SingleT
        model_TTbar = self.get_TTbar_mj_Model("_TTbar"+massscale);
        model_SingleT  = self.get_SingleT_mj_Model("_SingleT"+massscale);
        model_VV    = self.get_VV_mj_Model("_VV"+massscale);
        ## only two parameters are fix, offset and width while the exp is floating , otherwise if shape different User1 or ErfExp everything is flaoting
        model_WJets = self.get_WJets_mj_Model(label);

        ## Total Pdf and fit only in sideband 
        model_data = RooAddPdf("model_data%s_%s_mj"%(massscale,self.categoryLabel),"model_data%s_%s_mj"%(massscale,self.categoryLabel),RooArgList(model_WJets,model_VV,model_TTbar,model_SingleT));
        rfresult = model_data.fitTo( rdataset_data_mj, RooFit.Save(1) , RooFit.Range("lowersideband,uppersideband") ,RooFit.Extended(kTRUE), RooFit.NumCPU(2) );
        rfresult = model_data.fitTo( rdataset_data_mj, RooFit.Save(1) , RooFit.Range("lowersideband,uppersideband") ,RooFit.Extended(kTRUE), RooFit.NumCPU(2), RooFit.Minimizer("Minuit2") );
        rfresult.Print();
        rfresult.covarianceMatrix().Print();
        getattr(self.workspace4fit_,"import")(model_data);

        ## Total numver of event 
        rrv_number_data_mj = RooRealVar("rrv_number_data%s_%s_mj"%(massscale,self.categoryLabel),"rrv_number_data%s_%s_mj"%(massscale,self.categoryLabel),
                self.workspace4fit_.var("rrv_number_TTbar%s_%s_mj"%(massscale,self.categoryLabel)).getVal()+
                self.workspace4fit_.var("rrv_number_SingleT%s_%s_mj"%(massscale,self.categoryLabel)).getVal()+
                self.workspace4fit_.var("rrv_number_VV%s_%s_mj"%(massscale,self.categoryLabel)).getVal()+
                self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.categoryLabel)).getVal());

        rrv_number_data_mj.setError(TMath.Sqrt(self.workspace4fit_.var("rrv_number_TTbar%s_%s_mj"%(massscale,self.categoryLabel)).getError()*
            self.workspace4fit_.var("rrv_number_TTbar%s_%s_mj"%(massscale,self.categoryLabel)).getError()+
            self.workspace4fit_.var("rrv_number_SingleT%s_%s_mj"%(massscale,self.categoryLabel)).getError()*
            self.workspace4fit_.var("rrv_number_SingleT%s_%s_mj"%(massscale,self.categoryLabel)).getError()+
            self.workspace4fit_.var("rrv_number_VV%s_%s_mj"%(massscale,self.categoryLabel)).getError()*
            self.workspace4fit_.var("rrv_number_VV%s_%s_mj"%(massscale,self.categoryLabel)).getError()+
            self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.categoryLabel)).getError()*
            self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.categoryLabel)).getError()));
        getattr(self.workspace4fit_,"import")(rrv_number_data_mj);

        ## if fit on Wjets default with the default shape
        if TString(label).Contains("_WJets0"):

            ## make the final plot
            mplot = rrv_mass_j.frame(RooFit.Title(""), RooFit.Bins(int(rrv_mass_j.getBins()/self.BinWidth_narrow_factor)));
            rdataset_data_mj.plotOn(mplot, RooFit.Name("data_invisible"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

            ## plot solid style 
            model_data.plotOn(mplot,RooFit.Name("VV"), RooFit.Components("model%s_%s_mj,model_SingleT_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj"%(label,self.categoryLabel,self.categoryLabel,self.categoryLabel,self.categoryLabel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack),RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            model_data.plotOn(mplot,RooFit.Name("TTbar"), RooFit.Components("model%s_%s_mj,model_SingleT_%s_mj,model_TTbar_%s_mj"%(label,self.categoryLabel,self.categoryLabel,self.categoryLabel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack),RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            model_data.plotOn(mplot,RooFit.Name("SingleT"), RooFit.Components("model%s_%s_mj,model_SingleT_%s_mj"%(label,self.categoryLabel,self.categoryLabel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["SingleT"]), RooFit.LineColor(kBlack),RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            model_data.plotOn(mplot,RooFit.Name("WJets"), RooFit.Components("model%s_%s_mj"%(label,self.categoryLabel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack),RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            ## plot "dashed" style area
            model_data.plotOn(mplot,RooFit.Name("VV_invisible"), RooFit.Components("model%s_%s_mj,model_SingleT_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj"%(label,self.categoryLabel,self.categoryLabel,self.categoryLabel,self.categoryLabel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3002),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            model_data.plotOn(mplot,RooFit.Name("TTbar_invisible"), RooFit.Components("model%s_%s_mj,model_SingleT_%s_mj,model_TTbar_%s_mj"%(label,self.categoryLabel,self.categoryLabel,self.categoryLabel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3002),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            model_data.plotOn(mplot,RooFit.Name("SingleT_invisible"), RooFit.Components("model%s_%s_mj,model_SingleT_%s_mj"%(label,self.categoryLabel,self.categoryLabel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["SingleT"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3002),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());
            model_data.plotOn(mplot,RooFit.Name("WJets_invisible"), RooFit.Components("model%s_%s_mj"%(label,self.categoryLabel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]),RooFit.FillStyle(3002),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()), RooFit.LineColor(kBlack),RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());


            ### solid line
            model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj"%(label,self.categoryLabel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2) ,RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_SingleT_%s_mj"%(label,self.categoryLabel,self.categoryLabel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2),RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_SingleT_%s_mj,model_TTbar_%s_mj"%(label,self.categoryLabel,self.categoryLabel,self.categoryLabel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2) ,RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_SingleT_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj"%(label,self.categoryLabel,self.categoryLabel,self.categoryLabel,self.categoryLabel)),RooFit.LineColor(kBlack), RooFit.LineWidth(2) ,RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            ### dash line
            model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj"%(label,self.categoryLabel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.LineStyle(kDashed) ,RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_SingleT_%s_mj"%(label,self.categoryLabel,self.categoryLabel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.LineStyle(kDashed) ,RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_SingleT_%s_mj,model_TTbar_%s_mj"%(label,self.categoryLabel,self.categoryLabel,self.categoryLabel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.LineStyle(kDashed) ,RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_SingleT_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj"%(label,self.categoryLabel,self.categoryLabel,self.categoryLabel,self.categoryLabel)),RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.LineStyle(kDashed) ,RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_SingleT_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj"%(label,self.categoryLabel,self.categoryLabel,self.categoryLabel,self.categoryLabel)),RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.LineStyle(kDashed) ,RooFit.NormRange("lowersideband,uppersideband"));

            rdataset_data_mj.plotOn(mplot, RooFit.Name("data"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

            ### draw the error band using the sum of all the entries component MC + fit           
            draw_error_band(rdataset_data_mj, model_data, rrv_number_data_mj,rfresult,mplot,self.color_palet["Uncertainty"],"F");
            rdataset_data_mj.plotOn(mplot, RooFit.Name("data_invisible"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0) );

            ### Get the pull and plot it 
            mplot_pull=self.get_pull(rrv_mass_j,mplot);

            obs0_variable_signalregion_range_min= rrv_mass_j.getMin("signalregion");
            obs0_variable_signalregion_range_max= rrv_mass_j.getMax("signalregion");

            ### signal window zone with vertical lines
            lowerLine = TLine(obs0_variable_signalregion_range_min,0.,obs0_variable_signalregion_range_min,mplot.GetMaximum()*0.9); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kBlack); lowerLine.SetLineStyle(9);
            upperLine = TLine(obs0_variable_signalregion_range_max,0.,obs0_variable_signalregion_range_max,mplot.GetMaximum()*0.9); upperLine.SetLineWidth(2); upperLine.SetLineColor(kBlack); upperLine.SetLineStyle(9);
            mplot.addObject(lowerLine);
            mplot.addObject(upperLine);

            ### plot_legend of the plot
            self.plot_legend = self.legend4Plot(mplot,0,1,-0.10,-0.01,0.10,0.01);
            mplot.addObject(self.plot_legend);
            mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.6);

            parameters_list = model_data.getParameters(rdataset_data_mj);
            self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s_%s_%s/m_j_fitting_wtaggercut%s/"%(self.additioninformation, self.categoryLabel,self.PS_model,self.wtagger_label, self.wtagger_label), "m_j_sideband%s"%(label),"",1)

            ### call the function for getting the normalizatio in signal region for data, TTbar, SingleT, VV and W+jets = label -> store in a output txt file
            self.get_mj_normalization_insignalregion("_data");
            self.get_mj_normalization_insignalregion("_TTbar");
            self.get_mj_normalization_insignalregion("_SingleT");
            self.get_mj_normalization_insignalregion("_VV");
            self.get_mj_normalization_insignalregion(label);

        #### to calculate the WJets's normalization and error in M_J signalregion. The error must contain the shape error: model_WJets have new parameters fitting data
        fullInt   = model_WJets.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j) );
        signalInt = model_WJets.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("signalregion"));
        fullInt_val = fullInt.getVal()
        signalInt_val = signalInt.getVal()/fullInt_val
        ## take the value from the fit (normalization) and multiply it from the ratio of the integrals
        rrv_number_WJets_in_mj_signalregion_from_fitting = RooRealVar("rrv_number%s_in_mj_signalregion_from_fitting_%s"%(label,self.categoryLabel),"rrv_number%s_in_mj_signalregion_from_fitting_%s"%(label,self.categoryLabel),self.workspace4fit_.var("rrv_number%s_%s_mj"%(label,self.categoryLabel)).getVal()*signalInt_val);

        #### Error on the normalization --> from a dedicated function taking into account shape uncertainty
        rrv_number_WJets_in_mj_signalregion_from_fitting.setError( Calc_error_extendPdf(rdataset_data_mj, model_WJets, rfresult,"signalregion") );
        print "########## error on the normaliztion due to shape + norm = %s"%(rrv_number_WJets_in_mj_signalregion_from_fitting.getError());
        getattr(self.workspace4fit_,"import")(rrv_number_WJets_in_mj_signalregion_from_fitting);
        rrv_number_WJets_in_mj_signalregion_from_fitting.Print();


    ##### Counting of the events of each component in the signal region taking the lavel for the model
    def get_mj_normalization_insignalregion(self, label):
        print "################## get mj normalization ",label," ################## ";
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");
        model      = self.workspace4fit_.pdf("model"+label+"_"+self.categoryLabel+"_mj");

        fullInt   = model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j) );
        lowersidebandInt  = model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("lowersideband"));
        signalInt = model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("signalregion"));
        uppersidebandInt  = model.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("uppersideband"));

        fullInt_val   = fullInt.getVal()
        lowersidebandInt_val  = lowersidebandInt.getVal()/fullInt_val
        uppersidebandInt_val  = uppersidebandInt.getVal()/fullInt_val
        signalInt_val = signalInt.getVal()/fullInt_val

        print "########### Events Number in MC Dataset: #############"
        self.workspace4fit_.var("rrv_number_dataset_lowersideband"+label+"_"+self.categoryLabel+"_mj").Print();
        self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_mj").Print();
        self.workspace4fit_.var("rrv_number_dataset_uppersideband"+label+"_"+self.categoryLabel+"_mj").Print();

        print "########### Events Number get from fit: ##############"
        rrv_tmp = self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_mj");
        rrv_tmp.Print();
        print "Events Number in sideband_low :%s"%(rrv_tmp.getVal()*lowersidebandInt_val)
        print "Events Number in Signal Region:%s"%(rrv_tmp.getVal()*signalInt_val)
        print "Events Number in sideband_high:%s"%(rrv_tmp.getVal()*uppersidebandInt_val)
        print "Total Number in sidebands :%s"%(rrv_tmp.getVal()*(lowersidebandInt_val+uppersidebandInt_val) )
        print "Ratio signalregion/sidebands :%s"%(signalInt_val/(lowersidebandInt_val+uppersidebandInt_val) )

        ##### Save numbers in the output text file
        self.file_out.write( "\n%s++++++++++++++++++++++++++++++++++++"%(label) )
        self.file_out.write( "\nEvents Number in sideband_low from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_lowersideband"+label+"_"+self.categoryLabel+"_mj").getVal() ) )
        self.file_out.write( "\nEvents Number in Signal Region from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_mj").getVal() ) )
        self.file_out.write( "\nEvents Number in sideband_high from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_uppersideband"+label+"_"+self.categoryLabel+"_mj").getVal() ) )
        self.file_out.write( "\nTotal Number in sidebands from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_lowersideband"+label+"_"+self.categoryLabel+"_mj").getVal()+ self.workspace4fit_.var("rrv_number_dataset_uppersideband"+label+"_"+self.categoryLabel+"_mj").getVal() ) )
        self.file_out.write( "\nRatio signalregion/sidebands from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_mj").getVal()/(self.workspace4fit_.var("rrv_number_dataset_lowersideband"+label+"_"+self.categoryLabel+"_mj").getVal()+ self.workspace4fit_.var("rrv_number_dataset_uppersideband"+label+"_"+self.categoryLabel+"_mj").getVal()) ) )

        self.file_out.write( "\nEvents Number in sideband_low from fitting:%s"%(rrv_tmp.getVal()*lowersidebandInt_val) )
        self.file_out.write( "\nEvents Number in Signal Region from fitting:%s"%(rrv_tmp.getVal()*signalInt_val) )
        self.file_out.write( "\nEvents Number in sideband_high from fitting:%s"%(rrv_tmp.getVal()*uppersidebandInt_val) )
        self.file_out.write( "\nTotal Number in sidebands from fitting:%s"%(rrv_tmp.getVal()*(lowersidebandInt_val+uppersidebandInt_val) ) )
        self.file_out.write( "\nRatio signalregion/sidebands from fitting:%s"%(signalInt_val/(lowersidebandInt_val+uppersidebandInt_val) ) )

    ##### Method to fit data mlvj shape in the sideband -> first step for the background extraction of the shape
    def fit_mlvj_in_Mj_sideband(self, label, mlvj_region, mlvj_model,logy=0):

        print "############### Fit mlvj in mj sideband: ",label," ",mlvj_region,"  ",mlvj_model," ##################"
        rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j")
        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj")
        rdataset_data_mlvj = self.workspace4fit_.data("rdataset_data%s_%s_mlvj"%(mlvj_region,self.categoryLabel))

        ## get the minor component shapes in the sb low
        model_VV_backgrounds    = self.get_VV_mlvj_Model("_lowersideband");
        number_VV_lowersideband_mlvj    = self.workspace4fit_.var("rrv_number_VV_lowersideband_%s_mlvj"%(self.categoryLabel))
        model_TTbar_backgrounds = self.get_TTbar_mlvj_Model("_lowersideband");
        number_TTbar_lowersideband_mlvj = self.workspace4fit_.var("rrv_number_TTbar_lowersideband_%s_mlvj"%(self.categoryLabel))
        model_SingleT_backgrounds  = self.get_SingleT_mlvj_Model("_lowersideband");
        number_SingleT_lowersideband_mlvj  = self.workspace4fit_.var("rrv_number_SingleT_lowersideband_%s_mlvj"%(self.categoryLabel))

        self.workspace4fit_.var("rrv_number_TTbar_lowersideband_%s_mlvj"%(self.categoryLabel)).Print();
        self.workspace4fit_.var("rrv_number_SingleT_lowersideband_%s_mlvj"%(self.categoryLabel)).Print();
        self.workspace4fit_.var("rrv_number_VV_lowersideband_%s_mlvj"%(self.categoryLabel)).Print();

        ### Make the Pdf for the WJets
        model_pdf_WJets = self.make_Pdf("%s_lowersideband_from_fitting"%(label), mlvj_model,"_mlvj");
        model_pdf_WJets.Print();
        ### inititalize the value to what was fitted with the mc in the sideband
        number_WJets_lowersideband = self.workspace4fit_.var("rrv_number%s_lowersideband_%s_mlvj"%(label,self.categoryLabel)).clone("rrv_number%s_lowersideband_from_fitting_%s_mlvj"%(label,self.categoryLabel));
        model_WJets =RooExtendPdf("model%s_lowersideband_from_fitting_%s_mlvj"%(label,self.categoryLabel),"model%s_lowersideband_from_fitting_%s_mlvj"%(label,self.categoryLabel),model_pdf_WJets,number_WJets_lowersideband);
        model_pdf_WJets.Print();
        number_WJets_lowersideband.Print()

        ## Add the other bkg component fixed to the total model
        model_data = RooAddPdf("model_data%s%s_mlvj"%(label,mlvj_region),"model_data%s%s_mlvj"%(label,mlvj_region),RooArgList(model_WJets,model_VV_backgrounds, model_TTbar_backgrounds, model_SingleT_backgrounds));

        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE));
        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"));
        rfresult.Print();
        rfresult.covarianceMatrix().Print();
        getattr(self.workspace4fit_,"import")(model_data)

        model_WJets.Print();
        model_WJets.getParameters(rdataset_data_mlvj).Print("v");
        self.workspace4fit_.pdf("model_pdf%s_lowersideband_%s_mlvj"%(label,self.categoryLabel)).getParameters(rdataset_data_mlvj).Print("v");

        ### data in the sideband plus error from fit
        rrv_number_data_lowersideband_mlvj = RooRealVar("rrv_number_data_lowersideband_%s_mlvj"%(self.categoryLabel),"rrv_number_data_lowersideband_%s_mlvj"%(self.categoryLabel),
                self.workspace4fit_.var("rrv_number_TTbar_lowersideband_%s_mlvj"%(self.categoryLabel)).getVal()+
                self.workspace4fit_.var("rrv_number_SingleT_lowersideband_%s_mlvj"%(self.categoryLabel)).getVal()+
                self.workspace4fit_.var("rrv_number_VV_lowersideband_%s_mlvj"%(self.categoryLabel)).getVal()+
                self.workspace4fit_.var("rrv_number%s_lowersideband_from_fitting_%s_mlvj"%(label,self.categoryLabel)).getVal() );

        rrv_number_data_lowersideband_mlvj.setError( TMath.Sqrt(self.workspace4fit_.var("rrv_number%s_lowersideband_from_fitting_%s_mlvj"%(label,self.categoryLabel)).getError()*
            self.workspace4fit_.var("rrv_number%s_lowersideband_from_fitting_%s_mlvj"%(label,self.categoryLabel)).getError()+
            self.workspace4fit_.var("rrv_number_TTbar_lowersideband_%s_mlvj"%(self.categoryLabel)).getError()*
            self.workspace4fit_.var("rrv_number_TTbar_lowersideband_%s_mlvj"%(self.categoryLabel)).getError()+
            self.workspace4fit_.var("rrv_number_SingleT_lowersideband_%s_mlvj"%(self.categoryLabel)).getError()*
            self.workspace4fit_.var("rrv_number_SingleT_lowersideband_%s_mlvj"%(self.categoryLabel)).getError()+
            self.workspace4fit_.var("rrv_number_VV_lowersideband_%s_mlvj"%(self.categoryLabel)).getError()*
            self.workspace4fit_.var("rrv_number_VV_lowersideband_%s_mlvj"%(self.categoryLabel)).getError()));

        getattr(self.workspace4fit_,"import")(rrv_number_data_lowersideband_mlvj)

        ### plot for WJets default + default shape
        if TString(label).Contains("_WJets0"):

            mplot = rrv_mass_lvj.frame(RooFit.Title("M_lvj fitted in M_j sideband "), RooFit.Bins(int(rrv_mass_lvj.getBins()/self.BinWidth_narrow_factor)));

            rdataset_data_mlvj.plotOn( mplot , RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0) );

            model_data.plotOn(mplot, RooFit.Components("model%s_lowersideband_from_fitting_%s_mlvj,model_TTbar_lowersideband_%s_mlvj,model_SingleT_lowersideband_%s_mlvj,model_VV_lowersideband_%s_mlvj"%(label,self.categoryLabel,self.categoryLabel,self.categoryLabel,self.categoryLabel)), RooFit.Name("WJets"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_lowersideband_%s_mlvj,model_SingleT_lowersideband_%s_mlvj,model_VV_lowersideband_%s_mlvj"%(self.categoryLabel,self.categoryLabel,self.categoryLabel)),RooFit.Name("VV"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_lowersideband_%s_mlvj,model_SingleT_lowersideband_%s_mlvj"%(self.categoryLabel,self.categoryLabel)), RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components("model_SingleT_lowersideband_%s_mlvj"%(self.categoryLabel)), RooFit.Name("SingleT"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["SingleT"]), RooFit.LineColor(kBlack), RooFit.VLines());

            #solid line
            model_data.plotOn(mplot, RooFit.Components("model%s_lowersideband_from_fitting_%s_mlvj,model_TTbar_lowersideband_%s_mlvj,model_SingleT_lowersideband_%s_mlvj,model_VV_lowersideband_%s_mlvj"%(label,self.categoryLabel,self.categoryLabel,self.categoryLabel,self.categoryLabel)), RooFit.Name("WJets_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_lowersideband_%s_mlvj,model_SingleT_lowersideband_%s_mlvj,model_VV_lowersideband_%s_mlvj"%(self.categoryLabel,self.categoryLabel,self.categoryLabel)),RooFit.Name("VV_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_lowersideband_%s_mlvj,model_SingleT_lowersideband_%s_mlvj"%(self.categoryLabel,self.categoryLabel)), RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components("model_SingleT_lowersideband_%s_mlvj"%(self.categoryLabel)), RooFit.Name("SingleT_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());


            ### draw the error band 
            draw_error_band(rdataset_data_mlvj, model_data,self.workspace4fit_.var("rrv_number_data_lowersideband_%s_mlvj"%(self.categoryLabel)) ,rfresult,mplot,self.color_palet["Uncertainty"],"F");
            model_data.plotOn( mplot , RooFit.VLines(), RooFit.Invisible());
            model_data.plotOn( mplot , RooFit.Invisible());
            self.getData_PoissonInterval(rdataset_data_mlvj,mplot);

            mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);  

            ### Add the plot_legend to the plot 
            self.plot_legend=self.legend4Plot(mplot,0,1,0., 0.06, 0.16, 0.);
            mplot.addObject(self.plot_legend)

            ### calculate the chi2
            self.nPar_float_in_fitTo = rfresult.floatParsFinal().getSize();
            nBinX = mplot.GetNbinsX();
            ndof  = nBinX-self.nPar_float_in_fitTo;
            print mplot.chiSquare();
            print "#################### nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo ,mplot.chiSquare(self.nPar_float_in_fitTo)*ndof, ndof );
            ### write the result in the output
            self.file_out.write("\n fit_mlvj_in_Mj_sideband: nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare( self.nPar_float_in_fitTo )*ndof, ndof ) );

            ### get the pull plot and store the canvas
            mplot_pull = self.get_pull(rrv_mass_lvj,mplot);
            parameters_list = model_data.getParameters(rdataset_data_mlvj);

            self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s_%s_%s/m_lvj_fitting/"%(self.additioninformation, self.categoryLabel,self.PS_model,self.wtagger_label), "m_lvj_lowersideband%s"%(label),"",1,1)

        #### Decorrelate the parameters in order to have a proper shape in the workspace
        wsfit_tmp = RooWorkspace("wsfit_tmp%s_lowersideband_from_fitting_mlvj"%(label));
        Deco      = PdfDiagonalizer("Deco%s_lowersideband_from_fitting_%s_%s_mlvj"%(label,self.categoryLabel,self.wtagger_label),wsfit_tmp,rfresult);
        print"#################### diagonalize data sideband fit "
        model_pdf_WJets_deco = Deco.diagonalize(model_pdf_WJets);
        print"#################### print parameters "
        model_pdf_WJets_deco.Print("v");
        model_pdf_WJets_deco.getParameters(rdataset_data_mlvj).Print("");
        getattr(self.workspace4fit_,"import")(model_pdf_WJets_deco);

        #### Call the alpha evaluation in automatic
        self.get_WJets_mlvj_correction_lowersideband_to_signalregion(label,mlvj_model);

        ### Fix the pdf of signal, TTbar, SingleT and VV in the signal region 
        self.fix_Model("_%s"%(self.allsignals),"_signalregion","_mlvj")
        self.fix_Model("_TTbar","_signalregion","_mlvj")
        self.fix_Model("_SingleT","_signalregion","_mlvj")
        self.fix_Model("_VV","_signalregion","_mlvj")

        ### Call the evaluation of the normalization in the signal region for signal, TTbar, VV, SingleT, and WJets after the extrapolation via alpha
        self.get_pdf_signalregion_integral("_%s"%(self.allsignals));
        self.get_pdf_signalregion_integral("_TTbar");
        self.get_pdf_signalregion_integral("_SingleT");
        self.get_pdf_signalregion_integral("_VV");
        self.get_pdf_signalregion_integral(label,"model_pdf%s_signalregion_%s_after_correct_mlvj"%(label,self.categoryLabel));    


    ##### Function that calculate the normalization inside the mlvj signal region (mass window around the resonance in order to fill datacards)
    def get_pdf_signalregion_integral(self, label, model_name=""):

        print "############### get mlvj normalization inside SR ",label," ",model_name," ##################"
        if model_name == "":
            model = self.workspace4fit_.pdf("model"+label+"_signalregion"+"_"+self.categoryLabel+"_mlvj");
        else:
            model = self.workspace4fit_.pdf(model_name);

        rrv_mass_lvj = self.workspace4fit_.var("rrv_mass_lvj");

        fullInt   = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj) );
        signalInt = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj),("signalregion"));
        highMassInt = model.createIntegral(RooArgSet(rrv_mass_lvj),RooArgSet(rrv_mass_lvj),("high_mass"));

        fullInt_val = fullInt.getVal()
        signalInt_val = signalInt.getVal()/fullInt_val
        highMassInt_val = highMassInt.getVal()/fullInt_val 

        ## integal in the signal region
        print "######### integral in SR: ",label+"signalInt=%s"%(signalInt_val)

        print "####### Events Number in MC Dataset:"
        self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_mlvj").Print();
        self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.categoryLabel+"_mlvj").Print();

        print "########## Events Number get from fit:"
        rrv_tmp=self.workspace4fit_.var("rrv_number"+label+"_signalregion"+"_"+self.categoryLabel+"_mlvj");
        print "Events Number in Signal Region from fitting: %s"%(rrv_tmp.getVal()*signalInt_val)

        #### store the info in the output file
        self.file_out.write( "\n%s++++++++++++++++++++++++++++++++++++"%(label) )
        self.file_out.write( "\nEvents Number in All Region from dataset : %s"%(self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.categoryLabel+"_mlvj").getVal()) )
        self.file_out.write( "\nEvents Number in Signal Region from dataset: %s"%(self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_mlvj").getVal()) )
        self.file_out.write( "\nRatio signalregion/all_range from dataset :%s"%(self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_mlvj").getVal()/self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.categoryLabel+"_mlvj").getVal() ) )
        self.file_out.write( "\nEvents Number in All Region from fitting : %s\n"%(rrv_tmp.getVal()) )
        self.file_out.write( "\nEvents Number in Signal Region from fitting: %s\n"%(rrv_tmp.getVal()*signalInt_val) )
        self.file_out.write( "\nEvents Number in High Mass Region from fitting: %s\n"%(rrv_tmp.getVal()*highMassInt_val) )
        self.file_out.write( "\nRatio signalregion/all_range from fitting :%s"%(signalInt_val ) )

        if not self.workspace4fit_.var("rrv_number_fitting_signalregion"+label+"_"+self.categoryLabel+"_mlvj"):
            rrv_number_fitting_signalregion_mlvj = RooRealVar("rrv_number_fitting_signalregion"+label+"_"+self.categoryLabel+"_mlvj","rrv_number_fitting_signalregion"+label+"_"+
                    self.categoryLabel+"_mlvj", rrv_tmp.getVal()*signalInt_val );
            getattr(self.workspace4fit_,"import")(rrv_number_fitting_signalregion_mlvj);
        else :
            self.workspace4fit_.var("rrv_number_fitting_signalregion"+label+"_"+self.categoryLabel+"_mlvj").setVal(rrv_tmp.getVal()*signalInt_val);

        self.workspace4fit_.var("rrv_number_fitting_signalregion"+label+"_"+self.categoryLabel+"_mlvj").Print();

    ### method to get the alpha function to extrapolate the wjets in the signal region
    def get_WJets_mlvj_correction_lowersideband_to_signalregion(self,label, mlvj_model):

        print" ############# get the extrapolation function alpha from MC : ",label,"   ",mlvj_model," ###############";          
        tmp_Style = self.tdrStyle.Clone("tmp_Style");
        tmp_Style.SetPadRightMargin(0.08);
        tmp_Style.SetPadTickY(0);
        tmp_Style.cd();

        ### take input var and datasets from 4fit collection --> mc not scaled to lumi --> just a shape here 
        rrv_x = self.workspace4fit_.var("rrv_mass_lvj");
        rdataset_WJets_lowersideband_mlvj = self.workspace4fit_.data("rdataset4fit%s_lowersideband_%s_mlvj"%(label,self.categoryLabel))
        rdataset_WJets_signalregion_mlvj = self.workspace4fit_.data("rdataset4fit%s_signalregion_%s_mlvj"%(label,self.categoryLabel))

        ### create a frame for the next plots 
        mplot = rrv_x.frame(RooFit.Title("correlation_pdf"), RooFit.Bins(int(rrv_x.getBins()/self.BinWidth_narrow_factor))) ;
        mplot.GetYaxis().SetTitle("arbitrary units");

        ### model used for Higgs analysis --> parameters in the SR has to be fitted, not yet done in order to take into account correlations between mj and mlvj
        if mlvj_model=="ErfExp_v1":

            rrv_c_sb       = self.workspace4fit_.var("rrv_c_ErfExp%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_offset_sb  = self.workspace4fit_.var("rrv_offset_ErfExp%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_width_sb   = self.workspace4fit_.var("rrv_width_ErfExp%s_lowersideband_%s"%(label,self.categoryLabel));

            rrv_delta_c      = RooRealVar("rrv_delta_c_ErfExp%s_%s"%(label,self.categoryLabel),"rrv_delta_c_ErfExp%s_%s"%(label,self.categoryLabel),0.,
                    -100*rrv_c_sb.getError(),100*rrv_c_sb.getError());
            rrv_delta_offset = RooRealVar("rrv_delta_offset_ErfExp%s_%s"%(label,self.categoryLabel),"rrv_delta_offset_ErfExp%s_%s"%(label,self.categoryLabel),0.,
                    -100*rrv_offset_sb.getError(),100*rrv_offset_sb.getError());
            rrv_delta_width = RooRealVar("rrv_delta_width_ErfExp%s_%s"%(label,self.categoryLabel),"rrv_delta_width_ErfExp%s_%s"%(label,self.categoryLabel),0.,
                    -100*rrv_width_sb.getError(),100*rrv_width_sb.getError());

            rrv_c_sr      = RooFormulaVar("rrv_c_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_c_sb, rrv_delta_c ) );
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_offset_sb, rrv_delta_offset ) );
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_width_sb, rrv_delta_width ) );

            correct_factor_pdf = RooAlpha("correct_factor_pdf","correct_factor_pdf", rrv_x, rrv_c_sr, rrv_offset_sr,rrv_width_sr, rrv_c_sb, rrv_offset_sb, rrv_width_sb, rrv_x.getMin(), rrv_x.getMax());

        if mlvj_model=="ErfPow_v1":

            rrv_c_sb      = self.workspace4fit_.var("rrv_c_ErfPow%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_offset_sb = self.workspace4fit_.var("rrv_offset_ErfPow%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_width_sb  = self.workspace4fit_.var("rrv_width_ErfPow%s_lowersideband_%s"%(label,self.categoryLabel));

            rrv_delta_c      = RooRealVar("rrv_delta_c_ErfPow%s_%s"%(label,self.categoryLabel),"rrv_delta_c_ErfPow%s_%s"%(label,self.categoryLabel),0.,
                    -100*rrv_c_sb.getError(),100*rrv_c_sb.getError());
            rrv_delta_offset = RooRealVar("rrv_delta_offset_ErfPow%s_%s"%(label,self.categoryLabel),"rrv_delta_offset_ErfPow%s_%s"%(label,self.categoryLabel),0.,
                    -100*rrv_offset_sb.getError(),100*rrv_offset_sb.getError());
            rrv_delta_width  = RooRealVar("rrv_delta_width_ErfPow%s_%s"%(label,self.categoryLabel),"rrv_delta_width_ErfPow%s_%s"%(label,self.categoryLabel),0.,
                    -100*rrv_width_sb.getError(),100*rrv_width_sb.getError());

            rrv_c_sr      = RooFormulaVar("rrv_c_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_c_sb, rrv_delta_c ) );
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_offset_sb, rrv_delta_offset ) );
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_width_sb, rrv_delta_width ) );

            correct_factor_pdf = RooAlpha4ErfPowPdf("correct_factor_pdf","correct_factor_pdf", rrv_x, rrv_c_sr, rrv_offset_sr,rrv_width_sr, rrv_c_sb, rrv_offset_sb, rrv_width_sb);

        if mlvj_model=="ErfPow2_v1":

            rrv_c0_sb     = self.workspace4fit_.var("rrv_c0_ErfPow2%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_c1_sb     = self.workspace4fit_.var("rrv_c1_ErfPow2%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_offset_sb = self.workspace4fit_.var("rrv_offset_ErfPow2%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_width_sb  = self.workspace4fit_.var("rrv_width_ErfPow2%s_lowersideband_%s"%(label,self.categoryLabel));

            rrv_delta_c0      = RooRealVar("rrv_delta_c0_ErfPow2%s_%s"%(label,self.categoryLabel),"rrv_delta_c0_ErfPow2%s_%s"%(label,self.categoryLabel),-8, -20 ,0);
            rrv_delta_c1      = RooRealVar("rrv_delta_c1_ErfPow2%s_%s"%(label,self.categoryLabel),"rrv_delta_c1_ErfPow2%s_%s"%(label,self.categoryLabel),0., -5, 5);
            rrv_delta_offset  = RooRealVar("rrv_delta_offset_ErfPow2%s_%s"%(label,self.categoryLabel),"rrv_delta_offset_ErfPow2%s_%s"%(label,self.categoryLabel),30., 1.,80 );
            rrv_delta_width   = RooRealVar("rrv_delta_width_ErfPow2%s_%s"%(label,self.categoryLabel),"rrv_delta_width_ErfPow2%s_%s"%(label,self.categoryLabel),15,1.,100*rrv_width_sb.getError());

            rrv_c0_sr     = RooFormulaVar("rrv_c0_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_c0_sb, rrv_delta_c0 ) );
            rrv_c1_sr     = RooFormulaVar("rrv_c1_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_c1_sb, rrv_delta_c1 ) );
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_offset_sb, rrv_delta_offset ) );
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_width_sb, rrv_delta_width ) );

            correct_factor_pdf = RooAlpha4ErfPow2Pdf("correct_factor_pdf","correct_factor_pdf", rrv_x, rrv_c0_sr, rrv_c1_sr, rrv_offset_sr,rrv_width_sr, rrv_c0_sb, rrv_c1_sb, rrv_offset_sb, rrv_width_sb);

        if mlvj_model=="ErfPowExp_v1": ## take initial value from what was already fitted in the SR

            rrv_c0_sb     = self.workspace4fit_.var("rrv_c0_ErfPowExp%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_c1_sb     = self.workspace4fit_.var("rrv_c1_ErfPowExp%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_offset_sb = self.workspace4fit_.var("rrv_offset_ErfPowExp%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_width_sb  = self.workspace4fit_.var("rrv_width_ErfPowExp%s_lowersideband_%s"%(label,self.categoryLabel));

            rrv_delta_c0  = RooRealVar("rrv_delta_c0_ErfPowExp%s_%s"%(label,self.categoryLabel),"rrv_delta_c0_ErfPowExp%s_%s"%(label,self.categoryLabel),
                    self.workspace4fit_.var("rrv_c0_ErfPowExp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c0_sb.getVal(),
                    self.workspace4fit_.var("rrv_c0_ErfPowExp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c0_sb.getVal()-4*rrv_c0_sb.getError(),
                    self.workspace4fit_.var("rrv_c0_ErfPowExp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c0_sb.getVal()+4*rrv_c0_sb.getError() )

            rrv_delta_c1 = RooRealVar("rrv_delta_c1_ErfPowExp%s_%s"%(label,self.categoryLabel),"rrv_delta_c1_ErfPowExp%s_%s"%(label,self.categoryLabel),
                    self.workspace4fit_.var("rrv_c1_ErfPowExp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c1_sb.getVal(),
                    self.workspace4fit_.var("rrv_c1_ErfPowExp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c1_sb.getVal()-4*rrv_c1_sb.getError(),
                    self.workspace4fit_.var("rrv_c1_ErfPowExp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c1_sb.getVal()+4*rrv_c1_sb.getError() )

            rrv_delta_offset = RooRealVar("rrv_delta_offset_ErfPowExp%s_%s"%(label,self.categoryLabel),"rrv_delta_offset_ErfPowExp%s_%s"%(label,self.categoryLabel),
                    self.workspace4fit_.var("rrv_offset_ErfPowExp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_offset_sb.getVal(),
                    self.workspace4fit_.var("rrv_offset_ErfPowExp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_offset_sb.getVal()-4*rrv_offset_sb.getError(),
                    self.workspace4fit_.var("rrv_offset_ErfPowExp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_offset_sb.getVal()+4*rrv_offset_sb.getError())

            rrv_delta_width = RooRealVar("rrv_delta_width_ErfPowExp%s_%s"%(label,self.categoryLabel),"rrv_delta_width_ErfPowExp%s_%s"%(label,self.categoryLabel),
                    self.workspace4fit_.var("rrv_width_ErfPowExp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_width_sb.getVal(),
                    self.workspace4fit_.var("rrv_width_ErfPowExp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_width_sb.getVal()-4*rrv_width_sb.getError(),
                    self.workspace4fit_.var("rrv_width_ErfPowExp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_width_sb.getVal()+4*rrv_width_sb.getError() )

            rrv_c0_sr     = RooFormulaVar("rrv_c0_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_c0_sb, rrv_delta_c0 ) );
            rrv_c1_sr     = RooFormulaVar("rrv_c1_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_c1_sb, rrv_delta_c1 ) );
            rrv_offset_sr = RooFormulaVar("rrv_offset_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_offset_sb, rrv_delta_offset ) );
            rrv_width_sr  = RooFormulaVar("rrv_width_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_width_sb, rrv_delta_width ) );

            correct_factor_pdf = RooAlpha4ErfPowExpPdf("correct_factor_pdf","correct_factor_pdf", rrv_x, rrv_c0_sr, rrv_c1_sr, rrv_offset_sr,rrv_width_sr, rrv_c0_sb, rrv_c1_sb, rrv_offset_sb, rrv_width_sb);

        if mlvj_model=="Exp":
            rrv_c_sb    = self.workspace4fit_.var("rrv_c_Exp%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_delta_c = RooRealVar("rrv_delta_c_Exp%s_%s"%(label,self.categoryLabel),"rrv_delta_c_Exp%s_%s"%(label,self.categoryLabel),
                    self.workspace4fit_.var("rrv_c_Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c_sb.getVal(),
                    self.workspace4fit_.var("rrv_c_Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c_sb.getVal()-4*rrv_c_sb.getError(),
                    self.workspace4fit_.var("rrv_c_Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c_sb.getVal()+4*rrv_c_sb.getError() )

            correct_factor_pdf = RooExponential("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c);

        if mlvj_model=="2Exp":
            rrv_c0_sb    = self.workspace4fit_.var("rrv_c0_2Exp%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_delta_c0 = RooRealVar("rrv_delta_c0_2Exp%s_%s"%(label,self.categoryLabel),"rrv_delta_c0_2Exp%s_%s"%(label,self.categoryLabel),
                    self.workspace4fit_.var("rrv_c0_2Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c0_sb.getVal(),
                    self.workspace4fit_.var("rrv_c0_2Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c0_sb.getVal()-4*rrv_c0_sb.getError(),
                    self.workspace4fit_.var("rrv_c0_2Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c0_sb.getVal()+4*rrv_c0_sb.getError() )
            rrv_c0_sr = RooFormulaVar("rrv_c0_2Exp_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_c0_sb, rrv_delta_c0 ) );

            rrv_c1_sb = self.workspace4fit_.var("rrv_c1_2Exp%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_delta_c1 = RooRealVar("rrv_delta_c1_2Exp%s_%s"%(label,self.categoryLabel),"rrv_delta_c1_2Exp%s_%s"%(label,self.categoryLabel),
                    self.workspace4fit_.var("rrv_c1_2Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c1_sb.getVal(),
                    self.workspace4fit_.var("rrv_c1_2Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c1_sb.getVal()-4*rrv_c1_sb.getError(),
                    self.workspace4fit_.var("rrv_c1_2Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c1_sb.getVal()+4*rrv_c1_sb.getError() )
            rrv_c1_sr =RooFormulaVar("rrv_c1_2Exp_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_c1_sb, rrv_delta_c1 ) );

            rrv_frac_sb    = self.workspace4fit_.var("rrv_frac_2Exp%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_delta_frac = RooRealVar("rrv_delta_frac_2Exp%s_%s"%(label,self.categoryLabel),"rrv_delta_frac_2Exp%s_%s"%(label,self.categoryLabel),
                    self.workspace4fit_.var("rrv_frac_2Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_frac_sb.getVal(),
                    self.workspace4fit_.var("rrv_frac_2Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_frac_sb.getVal()-4*rrv_frac_sb.getError(),
                    self.workspace4fit_.var("rrv_frac_2Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_frac_sb.getVal()+4*rrv_frac_sb.getError() )
            rrv_frac_sr = RooFormulaVar("rrv_frac_2Exp_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_frac_sb, rrv_delta_frac ) );

            correct_factor_pdf = RooAlpha42ExpPdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_c0_sr,rrv_c1_sr,rrv_frac_sr, rrv_c0_sb,rrv_c1_sb,rrv_frac_sb );

        if mlvj_model=="Pow":

            rrv_c_sb    = self.workspace4fit_.var("rrv_c_Pow%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_delta_c = RooRealVar("rrv_delta_c_Pow%s_%s"%(label,self.categoryLabel),"rrv_delta_c_Pow%s_%s"%(label,self.categoryLabel),0., -100*rrv_c_sb.getError(),100*rrv_c_sb.getError());
            correct_factor_pdf = RooPowPdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c);

        if mlvj_model=="ExpN":
            rrv_c_sb  = self.workspace4fit_.var("rrv_c_ExpN%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_n_sb  = self.workspace4fit_.var("rrv_n_ExpN%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_delta_c = RooRealVar("rrv_delta_c_ExpN%s_%s"%(label,self.categoryLabel),"rrv_delta_c_ExpN%s_%s"%(label,self.categoryLabel),
                    self.workspace4fit_.var("rrv_c_ExpN%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c_sb.getVal(),
                    self.workspace4fit_.var("rrv_c_ExpN%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c_sb.getVal()-4*rrv_c_sb.getError(),
                    self.workspace4fit_.var("rrv_c_ExpN%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c_sb.getVal()+4*rrv_c_sb.getError() )
            rrv_delta_n = RooRealVar("rrv_delta_n_ExpN%s_%s"%(label,self.categoryLabel),"rrv_delta_n_ExpN%s_%s"%(label,self.categoryLabel),
                    self.workspace4fit_.var("rrv_n_ExpN%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_n_sb.getVal(),
                    self.workspace4fit_.var("rrv_n_ExpN%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_n_sb.getVal()-4*rrv_n_sb.getError(),
                    self.workspace4fit_.var("rrv_n_ExpN%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_n_sb.getVal()+4*rrv_n_sb.getError() )

            correct_factor_pdf = RooExpNPdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c, rrv_delta_n);

        if mlvj_model=="ExpTail":
            rrv_s_sb =self.workspace4fit_.var("rrv_s_ExpTail%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_a_sb =self.workspace4fit_.var("rrv_a_ExpTail%s_lowersideband_%s"%(label,self.categoryLabel));

            rrv_delta_s = RooRealVar("rrv_delta_s_ExpTail%s_%s"%(label,self.categoryLabel),"rrv_delta_s_ExpTail%s_%s"%(label,self.categoryLabel),
                    self.workspace4fit_.var("rrv_s_ExpTail%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_s_sb.getVal(),
                    self.workspace4fit_.var("rrv_s_ExpTail%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_s_sb.getVal()-4*rrv_s_sb.getError(),
                    self.workspace4fit_.var("rrv_s_ExpTail%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_s_sb.getVal()+4*rrv_s_sb.getError() )
            rrv_delta_a = RooRealVar("rrv_delta_a_ExpTail%s_%s"%(label,self.categoryLabel),"rrv_delta_a_ExpTail%s_%s"%(label,self.categoryLabel),
                    self.workspace4fit_.var("rrv_a_ExpTail%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_a_sb.getVal(),
                    self.workspace4fit_.var("rrv_a_ExpTail%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_a_sb.getVal()-4*rrv_a_sb.getError(),
                    self.workspace4fit_.var("rrv_a_ExpTail%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_a_sb.getVal()+4*rrv_a_sb.getError() )

            rrv_a_sr = RooFormulaVar("rrv_a_ExpTail_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_a_sb, rrv_delta_a ) );
            rrv_s_sr = RooFormulaVar("rrv_s_ExpTail_sr%s_%s"%(label,self.categoryLabel), "@0+@1",RooArgList(rrv_s_sb, rrv_delta_s ) );

            correct_factor_pdf = RooAlpha4ExpTailPdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_s_sr, rrv_a_sr, rrv_s_sb, rrv_a_sb);

        if mlvj_model=="Pow2":

            rrv_c0_sb    = self.workspace4fit_.var("rrv_c0_Pow2%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_c1_sb    = self.workspace4fit_.var("rrv_c1_Pow2%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_delta_c0 = RooRealVar("rrv_delta_c0_Pow2%s_%s"%(label,self.categoryLabel),"rrv_delta_c0_Pow2%s_%s"%(label,self.categoryLabel),0., -100*rrv_c0_sb.getError(),100*rrv_c0_sb.getError());
            rrv_delta_c1 = RooRealVar("rrv_delta_c1_Pow2%s_%s"%(label,self.categoryLabel),"rrv_delta_c1_Pow2%s_%s"%(label,self.categoryLabel),0., -100*rrv_c1_sb.getError(),100*rrv_c1_sb.getError());
            correct_factor_pdf = RooPow2Pdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c0,rrv_delta_c1);

        ### define the category and do the simultaneous fit taking the combined dataset of events in mlvj sb and sr

        data_category = RooCategory("data_category","data_category");
        data_category.defineType("sideband");
        data_category.defineType("signalregion");
        combData4fit = self.workspace4fit_.data("combData4fit%s_%s"%(label,self.categoryLabel));

        model_pdf_lowersideband_WJets         = self.workspace4fit_.pdf("model_pdf%s_lowersideband_%s_mlvj"%(label,self.categoryLabel));
        model_pdf_signalregion_WJets = RooProdPdf("model_pdf%s_signalregion_%s_mlvj"%(label,self.categoryLabel),"model_pdf%s_signalregion_%s_mlvj"%(label,self.categoryLabel) ,model_pdf_lowersideband_WJets,correct_factor_pdf);

        simPdf = RooSimultaneous("simPdf","simPdf",data_category);
        simPdf.addPdf(model_pdf_lowersideband_WJets,"sideband");
        simPdf.addPdf(model_pdf_signalregion_WJets,"signalregion");
        rfresult=simPdf.fitTo(combData4fit,RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE));
        rfresult=simPdf.fitTo(combData4fit,RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));
        rfresult.Print();
        rfresult.covarianceMatrix().Print();

        ### Decorrelate the parameters in the alpha shape
        wsfit_tmp = RooWorkspace("wsfit_tmp%s_sim_mlvj"%(label));
        print "############### diagonalizer alpha ";
        Deco      = PdfDiagonalizer("Deco%s_sim_%s_%s_mlvj"%(label,self.categoryLabel,self.wtagger_label),wsfit_tmp,rfresult);
        correct_factor_pdf_deco = Deco.diagonalize(correct_factor_pdf);
        correct_factor_pdf_deco.Print();
        correct_factor_pdf_deco.getParameters(rdataset_WJets_signalregion_mlvj).Print("v");
        getattr(self.workspace4fit_,"import")(correct_factor_pdf_deco);

        ## in case of default Wjets with default shape
        if TString(label).Contains("_WJets0"):

            ### only mc plots in the SB region
            mplot_lowersideband = rrv_x.frame(RooFit.Title("WJets sb low"), RooFit.Bins(int(rrv_x.getBins()/self.BinWidth_narrow_factor)));

            rdataset_WJets_lowersideband_mlvj.plotOn(mplot_lowersideband, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
            model_pdf_lowersideband_WJets.plotOn(mplot_lowersideband);
            mplot_pull_sideband = self.get_pull(rrv_x,mplot_lowersideband);
            parameters_list     = model_pdf_lowersideband_WJets.getParameters(rdataset_WJets_lowersideband_mlvj);
            mplot_lowersideband.GetYaxis().SetRangeUser(1e-2,mplot_lowersideband.GetMaximum()*1.2);
            self.draw_canvas_with_pull( mplot_lowersideband, mplot_pull_sideband,parameters_list,"plots_%s_%s_%s_%s/other/"%(self.additioninformation, self.categoryLabel,self.PS_model,self.wtagger_label), "m_lvj%s_lowersideband_sim"%(label),"",1,1)

            ### only mc plots in the SR region
            mplot_signalregion = rrv_x.frame(RooFit.Title("WJets sr"), RooFit.Bins(int(rrv_x.getBins()/self.BinWidth_narrow_factor)));

            rdataset_WJets_signalregion_mlvj.plotOn(mplot_signalregion, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
            model_pdf_signalregion_WJets.plotOn(mplot_signalregion);
            mplot_pull_signalregion = self.get_pull(rrv_x, mplot_signalregion);
            parameters_list = model_pdf_signalregion_WJets.getParameters(rdataset_WJets_signalregion_mlvj);
            mplot_signalregion.GetYaxis().SetRangeUser(1e-2,mplot_signalregion.GetMaximum()*1.2);
            self.draw_canvas_with_pull( mplot_signalregion, mplot_pull_signalregion,parameters_list,"plots_%s_%s_%s_%s/other/"%(self.additioninformation, self.categoryLabel,self.PS_model,self.wtagger_label), "m_lvj%s_signalregion_sim"%(label),"",1,1);

        ### Total plot shape in lowersideband, sr and alpha
        model_pdf_lowersideband_WJets.plotOn(mplot,RooFit.Name("Sideband"),RooFit.LineStyle(10));
        model_pdf_signalregion_WJets.plotOn(mplot, RooFit.LineColor(kRed) ,RooFit.LineStyle(8), RooFit.Name("Signal Region"));
        correct_factor_pdf_deco.plotOn(mplot, RooFit.LineColor(kBlack),RooFit.Name("#alpha") );

        ### plot also what is get from other source if available : alternate PS and shape: 1 PS and 01 is shape or fitting function
        if TString(label).Contains("_WJets0"):
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj"%(self.categoryLabel, self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj"%(self.categoryLabel, self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha: Alternate PS") );

            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_sim_%s_%s_mlvj"%(self.categoryLabel, self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_sim_%s_%s_mlvj"%(self.categoryLabel, self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7),RooFit.Name("#alpha: Alternate Function") );

        paras=RooArgList();
        ### Make a list of paramters as a function of the model after decorrelation 
        if mlvj_model=="ErfExp_v1" or mlvj_model=="ErfPow_v1" or mlvj_model=="2Exp" :
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig0"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig1"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig2"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig3"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig4"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig5"%(label,self.categoryLabel, self.wtagger_label) ));

        if mlvj_model=="ErfPow2_v1" or mlvj_model=="ErfPowExp_v1" :
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig0"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig1"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig2"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig3"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig4"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig5"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig6"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig7"%(label,self.categoryLabel, self.wtagger_label) ));

        if mlvj_model=="Exp" or mlvj_model=="Pow":
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig0"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig1"%(label,self.categoryLabel, self.wtagger_label) ));

        if mlvj_model=="ExpN" or mlvj_model=="ExpTail" or mlvj_model=="Pow2":
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig0"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig1"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig2"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_mlvj_eig3"%(label,self.categoryLabel, self.wtagger_label) ));

        if TString(label).Contains("_WJets0") or TString(label).Contains("_WJets1"): ### draw error band ar 1 and 2 sigma using the decorrelated shape
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj"%(label,self.categoryLabel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,1 ,mplot,kGray+3,"F",3001,"#alpha #pm",20,400);
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj"%(label,self.categoryLabel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,2 ,mplot,kGreen+2,"F",3002,"#alpha #pm",20,400);
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj"%(label,self.categoryLabel, self.wtagger_label),"rrv_mass_lvj", paras, self.workspace4fit_,1 ,mplot,kGray+3,"F",3001,"#alpha_invisible #pm",20,400);

        ### plot on the same canvas
        correct_factor_pdf_deco.plotOn(mplot, RooFit.LineColor(kBlack),RooFit.Name("#alpha_invisible") );

        if TString(label).Contains("_WJets0") : ## add also the plot of alternate ps and function on the canvas
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj"%(self.categoryLabel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj"%(self.categoryLabel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha_invisible: Alternate PS") );
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_sim_%s_%s_mlvj"%(self.categoryLabel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_sim_%s_%s_mlvj"%(self.categoryLabel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7),RooFit.Name("#alpha_invisible: Alternate Function") );

        elif TString(label).Contains("_WJets01") : ## add also the plot of alternate ps and function on the canvas
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj"%(self.categoryLabel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_mlvj"%(self.categoryLabel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha_invisible: Alternate PS") );
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets0_sim_%s_%s_mlvj"%(self.categoryLabel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets0_sim_%s_%s_mlvj"%(self.categoryLabel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7),RooFit.Name("#alpha_invisible: Alternate Function") );

        ### Add the plot_legend
        self.plot_legend=self.legend4Plot(mplot,1,0, -0.01, -0.14, 0.01, -0.06, 0.);
        mplot.addObject(self.plot_legend);

        ## set the Y axis in arbitrary unit 
        if self.allsignals=="ggH600" or self.allsignals=="ggH700": tmp_y_max=0.25
        else: tmp_y_max=0.28
        mplot.GetYaxis().SetRangeUser(0.,tmp_y_max);

        #### Draw another axis with the real value of alpha
        model_pdf_lowersideband_WJets.getVal(RooArgSet(rrv_x)),
        model_pdf_signalregion_WJets.getVal(RooArgSet(rrv_x)),
        correct_factor_pdf_deco.getVal(RooArgSet(rrv_x)),
        tmp_alpha_ratio = ( model_pdf_signalregion_WJets.getVal(RooArgSet(rrv_x))/model_pdf_lowersideband_WJets.getVal(RooArgSet(rrv_x)) );
        tmp_alpha_pdf   = correct_factor_pdf_deco.getVal(RooArgSet(rrv_x)) * mplot.getFitRangeBinW(); ## value of the pdf in each point
        tmp_alpha_scale = tmp_alpha_ratio/tmp_alpha_pdf;

        #add alpha scale axis
        axis_alpha=TGaxis( rrv_x.getMax(), 0, rrv_x.getMax(), tmp_y_max, 0, tmp_y_max*tmp_alpha_scale, 510, "+L");
        axis_alpha.SetTitle("#alpha");
        axis_alpha.SetTitleOffset(0.65);
        axis_alpha.SetTitleSize(0.05);
        axis_alpha.SetLabelSize(0.045);
        axis_alpha.SetTitleFont(42);
        axis_alpha.SetLabelFont(42);
        #axis_alpha.RotateTitle(1);
        mplot.addObject(axis_alpha);

        self.draw_canvas(mplot,"plots_%s_%s_%s_%s/other/"%(self.additioninformation, self.categoryLabel,self.PS_model,self.wtagger_label),"correction_pdf%s_%s_%s_M_lvj_signalregion_to_sideband"%(label,self.PS_model,mlvj_model),0,1);

        correct_factor_pdf_deco.getParameters(rdataset_WJets_lowersideband_mlvj).Print("v");
        model_pdf_WJets_lowersideband_from_fitting_mlvj_Deco = self.workspace4fit_.pdf("model_pdf%s_lowersideband_from_fitting_%s_mlvj_Deco%s_lowersideband_from_fitting_%s_%s_mlvj"%(label,self.categoryLabel,label, self.categoryLabel,self.wtagger_label));
        model_pdf_WJets_lowersideband_from_fitting_mlvj_Deco.Print("v");

        ### Wjets shape in the SR correctedfunction * sb 
        model_pdf_WJets_signalregion_after_correct_mlvj = RooProdPdf("model_pdf%s_signalregion_%s_after_correct_mlvj"%(label,self.categoryLabel),"model_pdf%s_signalregion_%s_after_correct_mlvj"%(label,self.categoryLabel),model_pdf_WJets_lowersideband_from_fitting_mlvj_Deco,self.workspace4fit_.pdf("correct_factor_pdf_Deco%s_sim_%s_%s_mlvj"%(label,self.categoryLabel,self.wtagger_label)) );
        model_pdf_WJets_signalregion_after_correct_mlvj.Print()
        ### fix the parmaters and import in the workspace
        getattr(self.workspace4fit_,"import")(model_pdf_WJets_signalregion_after_correct_mlvj)

        ##### calculate the normalization and alpha for limit datacard
        self.workspace4fit_.var("rrv_number%s_signalregion_%s_mlvj"%(label,self.categoryLabel)).Print();
        self.workspace4fit_.var("rrv_number%s_in_mj_signalregion_from_fitting_%s"%(label,self.categoryLabel)).Print();
        self.workspace4fit_.var("rrv_number%s_signalregion_%s_mlvj"%(label,self.categoryLabel)).setVal(self.workspace4fit_.var("rrv_number%s_in_mj_signalregion_from_fitting_%s"%(label,self.categoryLabel)).getVal());
        self.workspace4fit_.var("rrv_number%s_signalregion_%s_mlvj"%(label,self.categoryLabel)).setError(self.workspace4fit_.var("rrv_number%s_in_mj_signalregion_from_fitting_%s"%(label,self.categoryLabel)).getError());

        self.workspace4fit_.var("rrv_number%s_signalregion_%s_mlvj"%(label,self.categoryLabel)).setConstant(kTRUE);


  ### in order to make the plot_legend
    def legend4Plot(self, plot, left=1, isFill=1, x_offset_low=0., y_offset_low=0., x_offset_high =0., y_offset_high =0., TwoCoulum =1.):
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

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        

        if   self.categoryID == 0: legHeader="(e#nu, 1JLP)";
        elif self.categoryID == 1: legHeader="(e#nu, 1JHP)";
        elif self.categoryID == 2: legHeader="(#mu#nu, 1JLP)";
        elif self.categoryID == 3: legHeader="(#mu#nu, 1JHP)";

        for obj in range(int(plot.numItems()) ):
            objName = plot.nameOf(obj);
          if objName == "errorband" : objName = "Uncertainty";
          print objName;
          if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
                  objName ==objName_before ):
              theObj = plot.getObject(obj);
            objTitle = objName;
            drawoption= plot.getDrawOptions(objName).Data()
            if drawoption=="P":drawoption="PE"
            if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):  objName_before=objName; continue ;
            elif TString(objName).Contains("Graph") :  objName_before=objName; continue ;
        elif TString(objName).Data()=="data" : theLeg.AddEntry(theObj, "CMS Data "+legHeader,"PE");  objName_before=objName;                 
    else: objName_before=objName; continue ;

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        

        for obj in range(int(plot.numItems()) ):
            objName = plot.nameOf(obj);
          if objName == "errorband" : objName = "Uncertainty";
          print objName;
          if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
                  objName ==objName_before ):
              theObj = plot.getObject(obj);
            objTitle = objName;
            drawoption= plot.getDrawOptions(objName).Data()
            if drawoption=="P":drawoption="PE"
            if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):  objName_before=objName; continue ;
            elif TString(objName).Contains("Graph") :  objName_before=objName; continue ;
        elif TString(objName).Data()=="WJets" : theLeg.AddEntry(theObj, "W+jets","F");  objName_before=objName;                 
    else:  objName_before=objName; continue ;

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";        


        for obj in range(int(plot.numItems()) ):
            objName = plot.nameOf(obj);
            if objName == "errorband" : objName = "Uncertainty";
            print objName;
            if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or 
                    objName ==objName_before ):
                theObj = plot.getObject(obj);
                objTitle = objName;
                drawoption= plot.getDrawOptions(objName).Data()
                if drawoption=="P":drawoption="PE"
                if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):
                    theLeg.AddEntry(theObj, objName,"F");
                elif TString(objName).Contains("Graph") :
                    if not (objName_before=="Graph" or objName_before=="Uncertainty"): theLeg.AddEntry(theObj, "Uncertainty","F");
                    else:
                        if TString(objName).Data()=="SingleT" : theLeg.AddEntry(theObj, "Single Top","F");
                        elif TString(objName).Data()=="TTbar" : theLeg.AddEntry(theObj, "t#bar{t}","F");
                    elif TString(objName).Data()=="VV" : theLeg.AddEntry(theObj, "WW/WZ","F");
                elif TString(objName).Data()=="data" :  objName_before=objName; entryCnt = entryCnt+1; continue ;
            elif TString(objName).Data()=="WJets" : objName_before=objName; entryCnt = entryCnt+1; continue;
        elif TString(objName).Contains("vbfH"): theLeg.AddEntry(theObj, (TString(objName).ReplaceAll("vbfH","qqH")).Data() ,"L");
    elif TString(objName).Contains("Uncertainty"): theLeg.AddEntry(theObj, objTitle,drawoption);
elif TString(objName).Contains("Bulk"):
    if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M600") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M600"):
        objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.6 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M700") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M700"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.7 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M800") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M800"):
                           objName_signal_graviton = theObj ; 
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.8 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M900") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M900"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.9 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1000") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1000"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1100") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1100"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.1 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1200") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1200"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.2 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1300") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1300"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.3 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1400") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1400"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.4 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1500") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1500"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.5 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1600") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1600"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.6 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1700") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1700"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.7 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1800") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1800"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.8 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1900") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1900"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.9 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2000") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2000"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=2 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2100") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2100"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.1 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2200") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2200"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.2 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2300") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2300"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.3 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2400") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2400"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.4 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2500") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2500"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.5 TeV #tilde{k}=0.5 (#times100)";
                    else : theLeg.AddEntry(theObj, objTitle,drawoption);
                entryCnt=entryCnt+1;
            objName_before=objName;
        if objName_signal_graviton !="" :
            theLeg.AddEntry(objName_signal_graviton, TString(objNameLeg_signal_graviton).Data() ,"L");
        return theLeg;


    ##### Define the steps to fit WJets MC in the mj and mlvj spectra
    def fit_WJets(self):
        print "######################### fit_WJets ########################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0", "mJJNoKinFit")# to get the shape of m_lvj
        self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets01", "mJJNoKinFit")# to get the shape of m_lvj

        ### Fit in mj depends on the mlvj lower limit -> fitting the turn on at low mass or not
        if self.MODEL_4_mlvj=="ErfPowExp_v1" or self.MODEL_4_mlvj=="ErfPow2_v1" or self.MODEL_4_mlvj=="ErfExp_v1" :
            self.fit_obs_variable_SingleChannel(self.file_WJets0_mc,"_WJets0","ErfExp");
            self.fit_obs_variable_SingleChannel(self.file_WJets0_mc,"_WJets01","User1");
        else:
            self.fit_obs_variable_SingleChannel(self.file_WJets0_mc,"_WJets0","User1");
            self.fit_obs_variable_SingleChannel(self.file_WJets0_mc,"_WJets01","ErfExp");

        #### Fit the mlvj in lowersideband, signal region using two different model as done in the mj
        self.fit_limit_variable_SingleChannel(self.file_WJets0_mc,"_WJets0","_lowersideband",self.MODEL_4_mlvj, 0, 0, 1, 1);
        self.fit_limit_variable_SingleChannel(self.file_WJets0_mc,"_WJets0","_signalregion",self.MODEL_4_mlvj, 0, 0, 1, 1);
        self.fit_limit_variable_SingleChannel(self.file_WJets0_mc,"_WJets01","_lowersideband",self.MODEL_4_mlvj_alter, 0, 0, 1, 1);
        self.fit_limit_variable_SingleChannel(self.file_WJets0_mc,"_WJets01","_signalregion",self.MODEL_4_mlvj_alter, 0, 0, 1, 1);

        print "________________________________________________________________________"


    ##### Define the steps to fit VV MC in the mj and mlvj spectra
    def fit_VV(self):
        print "############################# fit_VV ################################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV", "mJJNoKinFit")

        ### fitting shape as a function of the mlvj region -> signal mass
        if self.MODEL_4_mlvj=="ErfPowExp_v1" or self.MODEL_4_mlvj=="ErfPow2_v1" or self.MODEL_4_mlvj=="ErfExp_v1":
            if self.wtagger_label=="LP":
                self.fit_obs_variable_SingleChannel(self.file_VV_mc,"_VV","ExpGaus");
            else:
                self.fit_obs_variable_SingleChannel(self.file_VV_mc,"_VV","2_2Gaus");
        else:
            if self.wtagger_label=="LP":
                self.fit_obs_variable_SingleChannel(self.file_VV_mc,"_VV","ExpGaus");
            else:
                self.fit_obs_variable_SingleChannel(self.file_VV_mc,"_VV","2_2Gaus");

        if self.MODEL_4_mlvj=="ErfPowExp_v1" or self.MODEL_4_mlvj=="ErfPow2_v1" or self.MODEL_4_mlvj=="ErfExp_v1":
            self.fit_limit_variable_SingleChannel(self.file_VV_mc,"_VV","_lowersideband","ErfExp_v1", 0, 0, 1);
            self.fit_limit_variable_SingleChannel(self.file_VV_mc,"_VV","_signalregion",self.MODEL_4_mlvj, 1, 0, 1);

        else:
            self.fit_limit_variable_SingleChannel(self.file_VV_mc,"_VV","_lowersideband","Exp", 0, 0, 1);
            self.fit_limit_variable_SingleChannel(self.file_VV_mc,"_VV","_signalregion",self.MODEL_4_mlvj, 1, 0, 1);

        print "________________________________________________________________________"

    ##### Define the steps to fit TTbar MC in the mj and mlvj spectra
    def fit_TTbar(self):
        print "################################ fit_TTbar #########################################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar", "mJJNoKinFit")# to get the shape of m_lvj

        if self.MODEL_4_mlvj=="ErfPowExp_v1" or self.MODEL_4_mlvj=="ErfPow2_v1" or self.MODEL_4_mlvj=="ErfExp_v1":
            if self.wtagger_label== "LP": self.fit_obs_variable_SingleChannel(self.file_TTbar_mc,"_TTbar","ExpGaus");
            else:                         self.fit_obs_variable_SingleChannel(self.file_TTbar_mc,"_TTbar","2Gaus_ErfExp");
        else:
            if self.wtagger_label== "LP" : self.fit_obs_variable_SingleChannel(self.file_TTbar_mc,"_TTbar","ExpGaus");
            else:                          self.fit_obs_variable_SingleChannel(self.file_TTbar_mc,"_TTbar","2Gaus_ErfExp");

        if self.MODEL_4_mlvj=="ErfPowExp_v1" or self.MODEL_4_mlvj=="ErfPow2_v1" or self.MODEL_4_mlvj=="ErfExp_v1" :
            self.fit_limit_variable_SingleChannel(self.file_TTbar_mc,"_TTbar","_lowersideband","ErfExp_v1", 0, 0, 1);
            self.fit_limit_variable_SingleChannel(self.file_TTbar_mc,"_TTbar","_signalregion",self.MODEL_4_mlvj,1, 0, 1);

        else:
            self.fit_limit_variable_SingleChannel(self.file_TTbar_mc,"_TTbar","_lowersideband","Exp");
            self.fit_limit_variable_SingleChannel(self.file_TTbar_mc,"_TTbar","_signalregion","Exp",1, 0, 1);

        print "________________________________________________________________________"


    #### Define the steps to fit SingleT MC in the mj and mlvj spectra
    def fit_SingleT(self):
        print "############################## fit_SingleT  #################################"
        self.get_mj_and_mlvj_dataset(self.file_SingleT_mc,"_SingleT", "mJJNoKinFit")
        self.fit_obs_variable_SingleChannel(self.file_SingleT_mc,"_SingleT","ExpGaus");

        if self.MODEL_4_mlvj=="ErfPowExp_v1" or self.MODEL_4_mlvj=="ErfPow2_v1" or self.MODEL_4_mlvj=="ErfExp_v1":
            self.fit_limit_variable_SingleChannel(self.file_SingleT_mc,"_SingleT","_lowersideband","ErfExp_v1", 0, 0, 1);
            self.fit_limit_variable_SingleChannel(self.file_SingleT_mc,"_SingleT","_signalregion","ErfExp_v1", 1, 0, 1);
        else:
            self.fit_limit_variable_SingleChannel(self.file_SingleT_mc,"_SingleT","_lowersideband","Exp", 0, 0, 1);
            self.fit_limit_variable_SingleChannel(self.file_SingleT_mc,"_SingleT","_signalregion","Exp", 1, 0, 1);

        print "________________________________________________________________________"

    ##### Fit of all the MC in both mj and mlvj : Signal, TTbar, SingleT, VV and Wjets
    def fit_AllSamples_Mj_and_Mlvj(self):
        print "################### fit_AllSamples_Mj_and_Mlvj #####################"
        self.fit_Signal()
        self.fit_WJets()
        self.fit_TTbar()
        self.fit_VV()
        self.fit_SingleT()
        print "________________________________________________________________________"


    ##### Analysis with sideband alpha correction 
    def analysis_sideband_correction_method1(self):
        print "##################### Start sideband correction full analysis ##############";
        ### Fit all MC components in both mj and mlvj
        self.fit_AllSamples_Mj_and_Mlvj();
        ### take the real data
        self.get_data()
        ### fit the WJets Normalization into the signal region -> no jet mass fluctuation has been done
        self.fit_WJetsNorm();
        ### fit data in the mlvj low sideband with two different models
        self.fit_mlvj_in_Mj_sideband("_WJets01","_lowersideband",self.MODEL_4_mlvj_alter,1)
        self.fit_mlvj_in_Mj_sideband("_WJets0","_lowersideband",self.MODEL_4_mlvj,1)

        ### Prepare the workspace and datacards     
        self.prepare_limit("sideband_correction_method1",1,0,0)
        ### finale plot and check of the workspace
        self.read_workspace(1)

    ##### Analysis with no shape uncertainty on alpha
    def analysis_sideband_correction_method1_without_shape_and_psmodel_systermatic(self):
        #### fit all the MC samples 
        self.fit_AllSamples_Mj_and_Mlvj()
        #### take the real data
        self.get_data()
        #### fit WJets just with one shape parametrization
        self.fit_WJetsNormalization_in_Mj_signalregion("_WJets0");
        #self.fit_WJetsNormalization_in_Mj_signalregion("_WJets0");
        #### fit sb lo with just one parametrization
        self.fit_mlvj_in_Mj_sideband("_WJets0","_lowersideband", self.MODEL_4_mlvj,1)
        #### prepare limit 
        self.prepare_limit("sideband_correction_method1",1,0,0)
        #### read the workspace
        self.read_workspace(1)
        '''


