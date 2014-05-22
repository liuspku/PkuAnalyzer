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

class doFit:

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
        self.BinWidth_narrow_factor=1.;## limit_variable 100/5=20 obs0_variable 5/5=1

        self.obs0_variable_BinWidth=self.analyzer_config["obs0_variable_BinWidth"];
        self.obs0_variable_BinWidth=self.obs0_variable_BinWidth/self.BinWidth_narrow_factor;
        obs0_variable_min=self.analyzer_config["obs0_variable_full_range_min"];
        obs0_variable_max=self.analyzer_config["obs0_variable_full_range_max"];
        nbins_obs0_variable=int( (obs0_variable_max - obs0_variable_min)/self.obs0_variable_BinWidth );
        ## correct obs0_variable max
        obs0_variable_max=obs0_variable_min+nbins_obs0_variable*self.obs0_variable_BinWidth;
        ## define jet mass variable
        #rrv_obs0_variable = RooRealVar("rrv_obs0_variable","Jet Mass",(obs0_variable_min+obs0_variable_max)/2.,obs0_variable_min,obs0_variable_max,"GeV");
        rrv_obs0_variable = RooRealVar("rrv_obs0_variable", self.analyzer_config["obs0_variable"][1],(obs0_variable_min+obs0_variable_max)/2.,obs0_variable_min,obs0_variable_max,"GeV");
        rrv_obs0_variable.setBins(nbins_obs0_variable);
        rrv_obs0_variable.setRange("lowersideband",self.analyzer_config["obs0_variable_lowersideband_range_min"],self.analyzer_config["obs0_variable_lowersideband_range_max"]);
        rrv_obs0_variable.setRange("signalregion",self.analyzer_config["obs0_variable_signalregion_range_min"],self.analyzer_config["obs0_variable_signalregion_range_max"]);
        rrv_obs0_variable.setRange("uppersideband",self.analyzer_config["obs0_variable_uppersideband_range_min"],self.analyzer_config["obs0_variable_uppersideband_range_max"]);
        getattr(self.workspace4fit_,"import")(rrv_obs0_variable);
        rrv_obs0_variable.Print();

        self.limit_variable_BinWidth=self.analyzer_config["limit_variable_BinWidth"];
        self.limit_variable_BinWidth=self.limit_variable_BinWidth/self.BinWidth_narrow_factor;
        limit_variable_min=self.analyzer_config["limit_variable_full_range_min"];
        limit_variable_max=self.analyzer_config["limit_variable_full_range_max"];
        nbins_limit_variable=int((limit_variable_max-limit_variable_min)/self.limit_variable_BinWidth);
        ## correct limit_variable max 
        limit_variable_max=limit_variable_min+nbins_limit_variable*self.limit_variable_BinWidth;
        ## define invariant mass WW variable
        #rrv_limit_variable= RooRealVar("rrv_limit_variable","M_{WW}",(limit_variable_min+limit_variable_max)/2.,limit_variable_min,limit_variable_max,"GeV");
        rrv_limit_variable= RooRealVar("rrv_limit_variable", self.analyzer_config["limit_variable"][1],(limit_variable_min+limit_variable_max)/2.,limit_variable_min,limit_variable_max,"GeV");
        rrv_limit_variable.setBins(nbins_limit_variable);
        rrv_limit_variable.setRange("signalregion", self.analyzer_config["limit_variable_signalregion_range_min"], self.analyzer_config["limit_variable_signalregion_range_max"] );# for cut&count limit calculation 
        getattr(self.workspace4fit_,"import")(rrv_limit_variable);
        rrv_limit_variable.Print();

        ## set the model used for the background parametrization
        self.MODEL_4_limit_variable       = self.analyzer_config["fit_model"];
        self.MODEL_4_limit_variable_alter = self.analyzer_config["fit_model_alter"];


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
        self.file_rlt_txt = self.rlt_DIR+"other_%s.txt"%(self.allsignals)
        ## workspace for limit
        self.file_rlt_root = self.rlt_DIR+"card4limit_%s_workspace.root"%(self.allsignals)
        ## datacard for cut&count and ubninned limit
        self.file_datacard_unbin = self.rlt_DIR+"card4limit_%s_unbin.txt"%(self.allsignals)
        self.file_datacard_counting = self.rlt_DIR+"card4limit_%s_counting.txt"%(self.allsignals)

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

    def analysis(self, analysis_strategy_name):
        print "################### fit_AllSamples limit_variables #####################"
        self.get_data();
        self.fit_Signals()
        self.fit_Backgrounds()

        if not analysis_strategy_name=="":
            print "Extend Strategy: ", analysis_strategy_name;
            if hasattr(self, analysis_strategy_name ):
                getattr(self, analysis_strategy_name)();
            else:
                raw_input("Couldn't find this strateg! ENTER to continue:");

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

 
    #### Define the steps to fit signal distribution in the obs0_variable and limit_variable spectra
    def fit_SingleChannel(self, file, label, fit_config,):
        print "############# fit_SingleChannel: %s, %s, %s #################"%(file, label, fit_config)
        ### Build the dataset
        self.get_obs0_variable_and_limit_variable_dataset(file, label, self.analyzer_config["limit_variable"][0], self.analyzer_config["obs0_variable"][0])

        self.fit_limit_variable_SingleChannel(label, "_signalregion", fit_config, 1, 0, 1);
        self.fit_limit_variable_SingleChannel(label, "_lowersideband", fit_config, 0, 0, 1);
        self.fit_obs_variable_SingleChannel(label,(0,"ExpGaus"));


        #model_narrow="DoubleCB_v1",model_width="BWDoubleCB"  
        #if (TString(self.file_signal).Contains("BulkG_WW_inclusive_M1000_W150") or TString(self.file_signal).Contains("BulkG_WW_inclusive_M1000_W50") or
        #        TString(self.file_signal).Contains("BulkG_WW_inclusive_M1000_W300") or TString(self.file_signal).Contains("BulkG_WW_inclusive_M1500_W225") or
        #        TString(self.file_signal).Contains("BulkG_WW_inclusive_M1500_W450") or TString(self.file_signal).Contains("BulkG_WW_inclusive_M1500_W75") or
        #        TString(self.file_signal).Contains("BulkG_WW_inclusive_M2100_W105") or TString(self.file_signal).Contains("BulkG_WW_inclusive_M2100_W315") or
        #        TString(self.file_signal).Contains("BulkG_WW_inclusive_M2100_W450")):
        #    self.fit_limit_variable_SingleChannel(self.file_signal,"_%s"%(self.allsignals),"_signalregion",model_width, 0, 0, 0, 0);            
        #else:
        #    self.fit_limit_variable_SingleChannel(self.file_signal,"_%s"%(self.allsignals),"_signalregion",model_narrow, 0, 0, 0, 0);

    ### Define the Extended Pdf for and limit_variable fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
    #def fit_limit_variable_SingleChannel(self,in_file_name, label, in_range, limit_variable_model, deco=0, show_constant_parameter=0, logy=0, ismc=0):
    def fit_limit_variable_SingleChannel(self, label, in_range, fit_config, show_constant_parameter=0, logy=0, ismc=0):

        print "############### Fit limit_variable single MC sample ",label,"  ",in_range,"  ",fit_config," ##################"
        ## import variable and dataset
        rrv_limit_variable = self.workspace4fit_.var("rrv_limit_variable")
        rdataset = self.workspace4fit_.data("rdataset4fit"+label+in_range+"_"+self.categoryLabel+"_limit_variable");
        constrainslist =[];

        ## make the extended pdf model
        if fit_config[1]=="Kyes": number_rdata=rdataset.sumEntries();
        else: number_rdata=500;
        model = self.make_Model(label+in_range,fit_config,"_limit_variable", number_rdata);

        ## make the fit
        model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE) );
        rfresult = model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult.Print();

        ## set the name of the result of the fit and put it in the workspace   
        rfresult.SetName("rfresult"+label+in_range+"_"+self.categoryLabel+"_limit_variable")
        getattr(self.workspace4fit_,"import")(rfresult)

        ## plot the result
        fitting_model=fit_config[1];
        mplot = rrv_limit_variable.frame(RooFit.Title(self.analyzer_config["limit_variable"][1]+in_range+" fitted by "+fitting_model), RooFit.Bins(int(rrv_limit_variable.getBins()/self.BinWidth_narrow_factor)));
        rdataset.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        ## plot the error band but don't store the canvas (only plotted without -b option
        draw_error_band_extendPdf(rdataset, model, rfresult,mplot,2,"L")
        rdataset.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        model.plotOn( mplot )#, RooFit.VLines()); in order to have the right pull 

        ## get the pull 
        mplot_pull      = self.get_pull(rrv_limit_variable,mplot);
        parameters_list = model.getParameters(rdataset);
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);

        self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s/limit_variable/"%(self.additioninformation, self.categoryLabel), label, self.analyzer_config["limit_variable"][0]+in_range+fitting_model, show_constant_parameter, logy);


        ## if the shape parameters has to be decorrelated
        #if deco :
        if fit_config[0] :
            print "################### Decorrelated limit_variable single mc shape ################"
            model_pdf = self.workspace4fit_.pdf("model_pdf%s%s_%s_limit_variable"%(label,in_range,self.categoryLabel)); ## take the pdf from the workspace
            model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) );
            rfresult_pdf = model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));
            rfresult_pdf.Print();

            ## temp workspace for the pdf diagonalizer
            wsfit_tmp = RooWorkspace("wsfit_tmp"+label+in_range+"_"+self.categoryLabel+"_limit_variable");
            Deco      = PdfDiagonalizer("Deco"+label+in_range+"_"+self.categoryLabel+"_"+self.wtagger_label+"_limit_variable",wsfit_tmp,rfresult_pdf); ## in order to have a good name 
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
            mplot_deco = rrv_limit_variable.frame( RooFit.Bins(int(rrv_limit_variable.getBins()/self.BinWidth_narrow_factor)));

            rdataset.plotOn(mplot_deco, RooFit.Name("Data"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
            model_pdf_deco.plotOn(mplot_deco,RooFit.Name(label),RooFit.LineColor(kBlack));

            mplot_deco.GetYaxis().SetRangeUser(1e-2,mplot_deco.GetMaximum()*1.2);

            rrv_number_dataset=RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset.sumEntries());
            rrv_number_dataset.setError(0.)
            draw_error_band(rdataset, model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,self.color_palet["Uncertainty"],"F"); ## don't store the number in the workspace

            self.plot_legend = legend4Plot(mplot_deco, self.categoryTitle,0); ## add the plot_legend                
            mplot_deco.addObject(self.plot_legend);

            self.draw_canvas( mplot_deco, "plots_%s_%s_%s_%s/other/"%(self.additioninformation, self.categoryLabel,self.PS_model, self.wtagger_label), self.analyzer_config["limit_variable"][0]+label+in_range+in_range+limit_variable_model+"_deco",0,logy)

        ### Number of the event in the dataset and lumi scale factor --> set the proper number for bkg extraction or for signal region
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_limit_variable").Print();
        self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).Print()
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_limit_variable").setVal( self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_limit_variable").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_limit_variable").setError(self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_limit_variable").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).getVal() )

        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_limit_variable").Print();

        #### Call the evaluation of the normalization in the signal region for signal, TTbar, VV, SingleT, and WJets after the extrapolation via alpha
        self.get_pdf_signalregion_integral(label) ;


    ### Method for a single MC fit of the obs0_variable spectra giving: file name, label, model name
    def fit_obs_variable_SingleChannel(self, label, fit_config, additioninformation=""):


        in_model_name  =fit_config[1];
        print "############### Fit obs0_variable single MC sample: ",label,"  by ",in_model_name," ##################"
        ## import variable and dataset
        rrv_obs0_variable = self.workspace4fit_.var("rrv_obs0_variable");
        rdataset_obs0_variable = self.workspace4fit_.data("rdataset4fit"+label+"_"+self.categoryLabel+"_obs0_variable");
        rdataset_obs0_variable.Print();

        ## make the extended model
        model = self.make_Model(label, fit_config,"_obs0_variable");
        rfresult = model.fitTo(rdataset_obs0_variable,RooFit.Save(1), RooFit.Extended(kTRUE) );
        rfresult = model.fitTo(rdataset_obs0_variable,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult = model.fitTo(rdataset_obs0_variable,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult.Print();

        ## Plot the result
        mplot = rrv_obs0_variable.frame(RooFit.Title(label+" fitted by "+in_model_name), RooFit.Bins(int(rrv_obs0_variable.getBins()/self.BinWidth_narrow_factor)) );
        rdataset_obs0_variable.plotOn( mplot, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        ## draw the error band for an extend pdf
        draw_error_band_extendPdf(rdataset_obs0_variable, model, rfresult,mplot,2,"L");
        ## re-draw the dataset
        rdataset_obs0_variable.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        ## draw the function
        model.plotOn( mplot );# remove RooFit.VLines() in order to get right pull in the 1st bin

        ## Get the pull
        mplot_pull = self.get_pull(rrv_obs0_variable, mplot); 
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);

        parameters_list = model.getParameters(rdataset_obs0_variable);
        self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s/obs0_variable/"%(self.additioninformation, self.categoryLabel), label, self.analyzer_config["limit_variable"][0]+in_model_name);

        #normalize the number of total events to lumi --> correct the number to scale to the lumi
        self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_obs0_variable").setVal( self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_obs0_variable").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_obs0_variable").setError(self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_obs0_variable").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).getVal() )

        ###### apply the correction of the mean and sigma from the ttbar control sample to the SingleT, TTbar and VV 
        #par=parameters_list.createIterator();
        #par.Reset();
        #param=par.Next()
        #while (param):
        #    if (TString(label).Contains("VV") or TString(label).Contains("SingleT") or TString(label).Contains("TTbar")):
        #        #param.Print();
        #        if TString(param.GetName()).Contains("rrv_mean1_gaus"):
        #            param.setRange(param.getMin()+self.mean_shift, param.getMax()+self.mean_shift);
        #            param.setVal(param.getVal()+self.mean_shift);
        #        if TString(param.GetName()).Contains("rrv_deltamean_gaus"):
        #            param.setRange(param.getMin()-self.mean_shift, param.getMax()-self.mean_shift);
        #            param.setVal(param.getVal()-self.mean_shift);
        #        if TString(param.GetName()).Contains("rrv_sigma1_gaus"):
        #            param.setVal(param.getVal()*self.sigma_scale);
        #            param.setRange(param.getMin()*self.sigma_scale, param.getMax()*self.sigma_scale);
        #        if TString(param.GetName()).Contains("rrv_scalesigma_gaus"):
        #            param.setRange(param.getMin()/self.sigma_scale, param.getMax()/self.sigma_scale);
        #            param.setVal(param.getVal()/self.sigma_scale);
        #    param=par.Next()


    ##### Function that calculate the normalization inside the limit_variable signal region (mass window around the resonance in order to fill datacards)
    def get_pdf_signalregion_integral(self, label, model_name=""):

        print "############### get limit_variable normalization inside SR ",label," ",model_name," ##################"
        if model_name == "":
            model = self.workspace4fit_.pdf("model"+label+"_signalregion"+"_"+self.categoryLabel+"_limit_variable");
        else:
            model = self.workspace4fit_.pdf(model_name);

        rrv_limit_variable = self.workspace4fit_.var("rrv_limit_variable");

        fullInt   = model.createIntegral(RooArgSet(rrv_limit_variable),RooArgSet(rrv_limit_variable) );
        signalInt = model.createIntegral(RooArgSet(rrv_limit_variable),RooArgSet(rrv_limit_variable),("signalregion"));

        fullInt_val = fullInt.getVal()
        signalInt_val = signalInt.getVal()/fullInt_val

        ## integal in the signal region
        print "######### integral in SR: ",label+"signalInt=%s"%(signalInt_val)

        print "####### Events Number in MC Dataset:"
        self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_limit_variable").Print();
        self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.categoryLabel+"_limit_variable").Print();

        print "########## Events Number get from fit:"
        rrv_tmp=self.workspace4fit_.var("rrv_number"+label+"_signalregion"+"_"+self.categoryLabel+"_limit_variable");
        print "Events Number in Signal Region from fitting: %s"%(rrv_tmp.getVal()*signalInt_val)

        #### store the info in the output file
        self.file_out.write( "\n%s++++++++++++++++++++++++++++++++++++"%(label) )
        self.file_out.write( "\nEvents Number in All Region from dataset : %s"%(self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.categoryLabel+"_limit_variable").getVal()) )
        self.file_out.write( "\nEvents Number in Signal Region from dataset: %s"%(self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_limit_variable").getVal()) )
        self.file_out.write( "\nRatio signalregion/all_range from dataset :%s"%(self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_limit_variable").getVal()/self.workspace4fit_.var("rrv_number_dataset_AllRange"+label+"_"+self.categoryLabel+"_limit_variable").getVal() ) )
        self.file_out.write( "\nEvents Number in All Region from fitting : %s\n"%(rrv_tmp.getVal()) )
        self.file_out.write( "\nEvents Number in Signal Region from fitting: %s\n"%(rrv_tmp.getVal()*signalInt_val) )
        self.file_out.write( "\nRatio signalregion/all_range from fitting :%s"%(signalInt_val ) )

        if not self.workspace4fit_.var("rrv_number_fitting_signalregion"+label+"_"+self.categoryLabel+"_limit_variable"):
            rrv_number_fitting_signalregion_limit_variable = RooRealVar("rrv_number_fitting_signalregion"+label+"_"+self.categoryLabel+"_limit_variable","rrv_number_fitting_signalregion"+label+"_"+
                    self.categoryLabel+"_limit_variable", rrv_tmp.getVal()*signalInt_val );
            getattr(self.workspace4fit_,"import")(rrv_number_fitting_signalregion_limit_variable);
        else :
            self.workspace4fit_.var("rrv_number_fitting_signalregion"+label+"_"+self.categoryLabel+"_limit_variable").setVal(rrv_tmp.getVal()*signalInt_val);

        self.workspace4fit_.var("rrv_number_fitting_signalregion"+label+"_"+self.categoryLabel+"_limit_variable").Print();


    #### run selection on data to build the datasets 
    def get_data(self):
        print "############### get_data ########################"
        self.get_obs0_variable_and_limit_variable_dataset(self.sig_bkg_files[2][2],"_data", self.analyzer_config["limit_variable"][0], self.analyzer_config["obs0_variable"][0])
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_dataset_signalregion_data_%s_limit_variable"%(self.categoryLabel)).clone("observation_for_counting"))

    ##### Method used to cycle on the events and for the dataset to be fitted
    def get_obs0_variable_and_limit_variable_dataset(self,in_file_name, label, limit_variable_name, obs0_variable_name ):# to get the shape of m_lvj,obs0_variable_name="jet_mass_pr"
        print "################### get_obs0_variable_and_limit_variable_dataset : ",in_file_name,"  ",label,"  ##################";
        fileIn_name = TString(in_file_name);
        #print fileIn_name;
        fileIn = TFile(fileIn_name.Data());
        treeIn = fileIn.Get("SelectedCandidatesPlain");

        rrv_obs0_variable   = self.workspace4fit_.var("rrv_obs0_variable")
        rrv_limit_variable = self.workspace4fit_.var("rrv_limit_variable")
        rrv_weight   = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.)

        ##### dataset of m_j -> scaleed and not scaled to lumi 
        rdataset_obs0_variable     = RooDataSet("rdataset"+label+"_"+self.categoryLabel+"_obs0_variable","rdataset"+label+"_"+self.categoryLabel+"_obs0_variable",RooArgSet(rrv_obs0_variable,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_obs0_variable = RooDataSet("rdataset4fit"+label+"_"+self.categoryLabel+"_obs0_variable","rdataset4fit"+label+"_"+self.categoryLabel+"_obs0_variable",RooArgSet(rrv_obs0_variable,rrv_weight),RooFit.WeightVar(rrv_weight) );
        ##### dataset of m_lvj -> scaled and not scaled to lumi in different region
        rdataset_lowersideband_limit_variable = RooDataSet("rdataset"+label+"_lowersideband"+"_"+self.categoryLabel+"_limit_variable","rdataset"+label+"_lowersideband"+"_"+self.categoryLabel+"_limit_variable",RooArgSet(rrv_limit_variable,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset_signalregion_limit_variable = RooDataSet("rdataset"+label+"_signalregion"+"_"+self.categoryLabel+"_limit_variable","rdataset"+label+"_signalregion"+"_"+self.categoryLabel+"_limit_variable",RooArgSet(rrv_limit_variable,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset_uppersideband_limit_variable = RooDataSet("rdataset"+label+"_uppersideband"+"_"+self.categoryLabel+"_limit_variable","rdataset"+label+"_uppersideband"+"_"+self.categoryLabel+"_limit_variable",RooArgSet(rrv_limit_variable,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset4fit_lowersideband_limit_variable = RooDataSet("rdataset4fit"+label+"_lowersideband"+"_"+self.categoryLabel+"_limit_variable","rdataset4fit"+label+"_lowersideband"+"_"+self.categoryLabel+"_limit_variable",RooArgSet(rrv_limit_variable,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset4fit_signalregion_limit_variable = RooDataSet("rdataset4fit"+label+"_signalregion"+"_"+self.categoryLabel+"_limit_variable","rdataset4fit"+label+"_signalregion"+"_"+self.categoryLabel+"_limit_variable",RooArgSet(rrv_limit_variable,rrv_weight),RooFit.WeightVar(rrv_weight) );

        rdataset4fit_uppersideband_limit_variable = RooDataSet("rdataset4fit"+label+"_uppersideband"+"_"+self.categoryLabel+"_limit_variable","rdataset4fit"+label+"_uppersideband"+"_"+self.categoryLabel+"_limit_variable",RooArgSet(rrv_limit_variable,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### categorize the event in sideband and signal region --> combined dataset 

        data_category = RooCategory("data_category","data_category");
        data_category.defineType("sideband");
        data_category.defineType("signalregion");
        combData = RooDataSet("combData"+label+"_"+self.categoryLabel,"combData"+label+"_"+self.categoryLabel,RooArgSet(rrv_limit_variable, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
        combData4fit = RooDataSet("combData4fit"+label+"_"+self.categoryLabel,"combData4fit"+label+"_"+self.categoryLabel,RooArgSet(rrv_limit_variable, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );

        obs0_variable_lowersideband_range_min= rrv_obs0_variable.getMin("lowersideband");
        obs0_variable_lowersideband_range_max= rrv_obs0_variable.getMax("lowersideband");
        obs0_variable_signalregion_range_min= rrv_obs0_variable.getMin("signalregion");
        obs0_variable_signalregion_range_max= rrv_obs0_variable.getMax("signalregion");
        obs0_variable_uppersideband_range_min= rrv_obs0_variable.getMin("uppersideband");
        obs0_variable_uppersideband_range_max= rrv_obs0_variable.getMax("uppersideband");

        limit_variable_signalregion_range_min = rrv_limit_variable.getMin("signalregion");
        limit_variable_signalregion_range_max = rrv_limit_variable.getMax("signalregion");
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

            tmp_obs0_variable=getattr(treeIn, obs0_variable_name);
            tmp_limit_variable=getattr(treeIn, limit_variable_name);

            self.isGoodEvent = 0 ;   
            if treeIn.categories==self.categoryID and tmp_limit_variable> rrv_limit_variable.getMin() and tmp_limit_variable<rrv_limit_variable.getMax() and tmp_obs0_variable>rrv_obs0_variable.getMin() and tmp_obs0_variable<rrv_obs0_variable.getMax() :
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
                #print i, tmp_limit_variable, tmp_obs0_variable, tmp_event_weight

                rrv_limit_variable.setVal(tmp_limit_variable);

                if tmp_obs0_variable >= obs0_variable_lowersideband_range_min and tmp_obs0_variable < obs0_variable_lowersideband_range_max:
                    rdataset_lowersideband_limit_variable.add( RooArgSet( rrv_limit_variable ), tmp_event_weight );
                    rdataset4fit_lowersideband_limit_variable.add( RooArgSet( rrv_limit_variable ), tmp_event_weight4fit );

                    data_category.setLabel("sideband");
                    combData.add(RooArgSet(rrv_limit_variable,data_category),tmp_event_weight);
                    combData4fit.add(RooArgSet(rrv_limit_variable,data_category),tmp_event_weight4fit);

                if tmp_obs0_variable >= obs0_variable_signalregion_range_min and tmp_obs0_variable < obs0_variable_signalregion_range_max:
                    rdataset_signalregion_limit_variable.add( RooArgSet( rrv_limit_variable ), tmp_event_weight );
                    rdataset4fit_signalregion_limit_variable.add( RooArgSet( rrv_limit_variable ), tmp_event_weight4fit );

                    data_category.setLabel("signalregion");
                    combData.add(RooArgSet(rrv_limit_variable,data_category),tmp_event_weight);
                    combData4fit.add(RooArgSet(rrv_limit_variable,data_category),tmp_event_weight4fit);
                    hnum_2region.Fill(1,tmp_event_weight);

                    if tmp_limit_variable >=limit_variable_signalregion_range_min and tmp_limit_variable <limit_variable_signalregion_range_max:
                        hnum_2region.Fill(0,tmp_event_weight);

                if tmp_obs0_variable >= obs0_variable_uppersideband_range_min and tmp_obs0_variable < obs0_variable_uppersideband_range_max:
                    rdataset_uppersideband_limit_variable.add( RooArgSet( rrv_limit_variable ), tmp_event_weight );
                    rdataset4fit_uppersideband_limit_variable.add( RooArgSet( rrv_limit_variable ), tmp_event_weight4fit );

                rrv_obs0_variable.setVal( tmp_obs0_variable );
                rdataset_obs0_variable.add( RooArgSet( rrv_obs0_variable ), tmp_event_weight );
                rdataset4fit_obs0_variable.add( RooArgSet( rrv_obs0_variable ), tmp_event_weight4fit );

                if tmp_obs0_variable >=obs0_variable_lowersideband_range_min and tmp_obs0_variable <obs0_variable_lowersideband_range_max:
                    hnum_4region.Fill(-1,tmp_event_weight );
                if tmp_obs0_variable >=obs0_variable_signalregion_range_min and tmp_obs0_variable <obs0_variable_signalregion_range_max :
                    hnum_4region.Fill(0,tmp_event_weight);
                if tmp_obs0_variable >=obs0_variable_uppersideband_range_min and tmp_obs0_variable <obs0_variable_uppersideband_range_max:
                    hnum_4region.Fill(1,tmp_event_weight);

                hnum_4region.Fill(2,tmp_event_weight);

        ### scaler to lumi for MC in 4fit datasets
        rrv_scale_to_lumi=RooRealVar("rrv_scale_to_lumi"+label+"_"+self.categoryLabel,"rrv_scale_to_lumi"+label+"_"+self.categoryLabel,tmp_scale_to_lumi)
        rrv_scale_to_lumi.Print()
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi)

        ### prepare m_lvj dataset to be compared with the fit results
        rrv_number_dataset_signalregion_limit_variable=RooRealVar("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_limit_variable","rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_limit_variable",hnum_2region.GetBinContent(1));
        rrv_number_dataset_AllRange_limit_variable=RooRealVar("rrv_number_dataset_AllRange"+label+"_"+self.categoryLabel+"_limit_variable","rrv_number_dataset_AllRange"+label+"_"+self.categoryLabel+"_limit_variable",hnum_2region.GetBinContent(2));

        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signalregion_limit_variable)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_AllRange_limit_variable)
        ### import the dataser       
        getattr(self.workspace4fit_,"import")(rdataset_lowersideband_limit_variable);
        getattr(self.workspace4fit_,"import")(rdataset_signalregion_limit_variable);
        getattr(self.workspace4fit_,"import")(rdataset_uppersideband_limit_variable);
        getattr(self.workspace4fit_,"import")(rdataset4fit_lowersideband_limit_variable);
        getattr(self.workspace4fit_,"import")(rdataset4fit_signalregion_limit_variable);
        getattr(self.workspace4fit_,"import")(rdataset4fit_uppersideband_limit_variable);
        getattr(self.workspace4fit_,"import")(combData);
        getattr(self.workspace4fit_,"import")(combData4fit);

        ### write in the output 
        self.file_out.write("\n%s events number in %s from dataset: %s"%(label, self.analyzer_config["limit_variable"][0],rdataset_signalregion_limit_variable.sumEntries()))

        ### prepare m_j dataset
        rrv_number_dataset_lowersideband_obs0_variable=RooRealVar("rrv_number_dataset_lowersideband"+label+"_"+self.categoryLabel+"_obs0_variable","rrv_number_dataset_lowersideband"+label+"_"+self.categoryLabel+"_obs0_variable",hnum_4region.GetBinContent(1));
        rrv_number_dataset_signalregion_obs0_variable=RooRealVar("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_obs0_variable","rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_obs0_variable",hnum_4region.GetBinContent(2));
        rrv_number_dataset_uppersideband_obs0_variable=RooRealVar("rrv_number_dataset_uppersideband"+label+"_"+self.categoryLabel+"_obs0_variable","rrv_number_dataset_uppersideband"+label+"_"+self.categoryLabel+"_obs0_variable",hnum_4region.GetBinContent(3));
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_lowersideband_obs0_variable)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signalregion_obs0_variable)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_uppersideband_obs0_variable)

        getattr(self.workspace4fit_,"import")(rdataset_obs0_variable)
        getattr(self.workspace4fit_,"import")(rdataset4fit_obs0_variable)

        #### print everything
        rdataset_lowersideband_limit_variable.Print();
        rdataset_signalregion_limit_variable.Print();
        rdataset_uppersideband_limit_variable.Print();
        rdataset_obs0_variable.Print();
        rdataset4fit_lowersideband_limit_variable.Print();
        rdataset4fit_signalregion_limit_variable.Print();
        rdataset4fit_uppersideband_limit_variable.Print();
        rdataset4fit_obs0_variable.Print();
        rrv_number_dataset_signalregion_limit_variable.Print()
        rrv_number_dataset_AllRange_limit_variable.Print()
        rrv_number_dataset_lowersideband_obs0_variable.Print()
        rrv_number_dataset_signalregion_obs0_variable.Print()
        rrv_number_dataset_uppersideband_obs0_variable.Print()
        print rdataset_signalregion_limit_variable.sumEntries()
        print rrv_number_dataset_signalregion_limit_variable.getVal()
        print rrv_number_dataset_AllRange_limit_variable.getVal()

    ### Define the Extended Pdf for and mJ fit giving: label, fit model name, list constraint and ismc
    #def make_Model(self, label, fit_config, mass_spectrum="_obs0_variable", ConstraintsList=[], ismc_wjet=0, area_init_value=500):
    def make_Model(self, label, fit_config, mass_spectrum, area_init_value=500):

        ##### define an extended pdf from a standard Roofit One
        print " "
        print "###############################################"
        print "## Make model : ",label," ",fit_config,"##";
        print "###############################################"
        print " "

        rrv_number = RooRealVar("rrv_number"+label+"_"+self.categoryLabel+mass_spectrum,"rrv_number"+label+"_"+self.categoryLabel+mass_spectrum,area_init_value,0.,1e7);
        ## call the make RooAbsPdf method
        model_pdf = make_Pdf(label+"_"+self.categoryLabel, self.workspace4fit_, fit_config, mass_spectrum)
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

        rrv_x = workspace.var("rrv_limit_variable")
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
        self.getData_PoissonInterval(rrv_x, data_obs,mplot);

        model_Total_bkgs.plotOn(mplot,RooFit.Normalization(scale_number_Total_bkgs),RooFit.Invisible());
        mplot_pull=self.get_pull(rrv_x,mplot);

        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        draw_error_band(model_Total_bkgs, rrv_x.GetName(), rrv_number_Total_bkgs,self.FloatingParams,workspace ,mplot,self.color_palet["Uncertainty"],"F");

        mplot.Print();
        self.plot_legend = legend4Plot(mplot, self.categoryTitle,0,1,-0.01,-0.05,0.11,0.);
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


    def getData_PoissonInterval(self, rrv_x, data_obs,mplot):
        #rrv_x = self.workspace4fit_.var("rrv_limit_variable");
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

        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_limit_variable"));

        #### number of events in signal region for every sigs and bkgs for cut-and-couting limit
        for iter in range(self.nsig):
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signalregion_%s_%s_limit_variable"%(self.sig_list[iter][0], self.categoryLabel)).clone("rate_%s_for_counting"%(self.sig_list[iter][0]) ))
        for iter in range(self.nbkg):
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signalregion_%s_%s_limit_variable"%(self.bkg_list[iter][0], self.categoryLabel)).clone("rate_%s_for_counting"%(self.bkg_list[iter][0]) ))

        ##### number of signal, Wjets, VV, TTbar and SingleT --> unbin
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_%s_signalregion_%s_limit_variable"%(self.allsignals, self.categoryLabel)).clone("rate_%s_for_unbin"%(self.allsignals)));
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_WJets0_signalregion_%s_limit_variable"%(self.categoryLabel)).clone("rate_WJets_for_unbin"));
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_VV_signalregion_%s_limit_variable"%(self.categoryLabel)).clone("rate_VV_for_unbin"));
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_TTbar_signalregion_%s_limit_variable"%(self.categoryLabel)).clone("rate_TTbar_for_unbin"));
        #getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_SingleT_signalregion_%s_limit_variable"%(self.categoryLabel)).clone("rate_SingleT_for_unbin"));

        #### Set the error properly -> taking into account lumi, Vtagger and theoretical uncertainty on XS -> for VV, TTbar and SingleT
        #self.workspace4limit_.var("rate_VV_for_unbin").setError(self.workspace4limit_.var("rate_VV_for_unbin").getVal()*TMath.Sqrt( self.lumi_uncertainty*self.lumi_uncertainty + self.rrv_wtagger_eff_reweight_forV.getError()/self.rrv_wtagger_eff_reweight_forV.getVal()*self.rrv_wtagger_eff_reweight_forV.getError()/self.rrv_wtagger_eff_reweight_forV.getVal() +self.XS_VV_uncertainty*self.XS_VV_uncertainty ) );
        #self.workspace4limit_.var("rate_SingleT_for_unbin").setError(self.workspace4limit_.var("rate_SingleT_for_unbin").getVal()*TMath.Sqrt( self.lumi_uncertainty*self.lumi_uncertainty + self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal()*self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal() +self.XS_SingleT_uncertainty*self.XS_SingleT_uncertainty ) );
        #self.workspace4limit_.var("rate_TTbar_for_unbin").setError(self.workspace4limit_.var("rate_TTbar_for_unbin").getVal()*TMath.Sqrt( self.lumi_uncertainty*self.lumi_uncertainty + self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal()*self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal() ))

        #### Get the dataset for data into the signal region
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.data("rdataset_data_signalregion_%s_limit_variable"%(self.categoryLabel)).Clone("data_obs_%s"%(self.categoryLabel)))

        for iter_sig in range(self.nsig):
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_%s_signalregion_%s_limit_variable"%(self.sig_list[iter_sig][0], self.categoryLabel)).clone("rate_%s_for_unbin"%(self.sig_list[iter_sig][0])));
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_%s_signalregion_%s_limit_variable"%(self.sig_list[iter_sig][0], self.categoryLabel)).clone("%s_%s"%(self.sig_list[iter_sig][0], self.categoryLabel)));

        for iter_bkg in range(self.nbkg):
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_%s_signalregion_%s_limit_variable"%(self.bkg_list[iter_bkg][0], self.categoryLabel)).clone("rate_%s_for_unbin"%(self.bkg_list[iter_bkg][0])));
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_%s_signalregion_%s_limit_variable"%(self.bkg_list[iter_bkg][0], self.categoryLabel)).clone("%s_%s"%(self.bkg_list[iter_bkg][0], self.categoryLabel)));

        #### Take the corrected pdf from the alpha method for the WJets
        #if analysis_mode=="sideband_correction_method1":
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_WJets0_signalregion_%s_after_correct_limit_variable"%(self.categoryLabel)).clone("WJets_%s_%s"%(self.categoryLabel, self.wtagger_label)));

        #if isTTbarFloating:
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_TTbar_signalregion_%s_limit_variable_Deco_TTbar_signalregion_%s_%s_limit_variable"%(self.categoryLabel, self.categoryLabel, self.wtagger_label)).clone("TTbar_%s_%s"%(self.categoryLabel,self.wtagger_label)))
        #else :
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_TTbar_signalregion_%s_limit_variable"%(self.categoryLabel, self.categoryLabel, self.wtagger_label)).clone("TTbar_%s_%s"%(self.categoryLabel,self.wtagger_label)))

        #if isSingleTFloating :     
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_SingleT_signalregion_%s_limit_variable_Deco_SingleT_signalregion_%s_%s_limit_variable"%(self.categoryLabel, self.categoryLabel, self.wtagger_label)).clone("SingleT_%s_%s"%(self.categoryLabel,self.wtagger_label)))
        #else :
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_SingleT_signalregion_%s_limit_variable"%(self.categoryLabel)).clone("SingleT_%s_%s"%(self.categoryLabel,self.wtagger_label)))

        #if isVVFloating :    
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_VV_signalregion_%s_limit_variable_Deco_VV_signalregion_%s_%s_limit_variable"%(self.categoryLabel, self.categoryLabel, self.wtagger_label)).clone("VV_%s_%s"%(self.categoryLabel,self.wtagger_label)))
        #else:
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_VV_signalregion_%s_limit_variable"%(self.categoryLabel)).clone("VV_%s_%s"%(self.categoryLabel,self.wtagger_label)))

        #if TString(self.allsignals).Contains("BulkG_WW"):
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_%s_signalregion_%s_limit_variable"%(self.allsignals,self.categoryLabel)).clone("BulkWW_%s_%s"%(self.categoryLabel, self.wtagger_label)))
        #else:    
        #    getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_%s_signalregion_%s_limit_variable"%(self.allsignals,self.categoryLabel)).clone(self.allsignals+"_%s_%s"%(self.categoryLabel, self.wtagger_label)))

        #### Fix all the Pdf parameters 
        #rrv_x = self.workspace4limit_.var("rrv_limit_variable");

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

        #    if self.MODEL_4_limit_variable=="ErfExp_v1" or self.MODEL_4_limit_variable=="ErfPow_v1" or self.MODEL_4_limit_variable=="2Exp" :
        #        ### uncertainty inflation on the Wjets shape from fitting data in lowersideband
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig3"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        ### Add to the parameter list
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig3"%(self.categoryLabel, self.wtagger_label)));

        #        ### Do the same for alpha paramter
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig3"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig4"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig5"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        ### Add to the parameter list
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig3"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig4"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0sim_%s_%s_limit_variable_eig5"%(self.categoryLabel, self.wtagger_label)))

        #        ### Do the same for the TTbar
        #        if isTTbarFloating !=0 :
        #            self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #         self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #         self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);

        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)));
        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)));
        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)));


        #    if self.MODEL_4_limit_variable=="ErfPow2_v1" or self.MODEL_4_limit_variable=="ErfPowExp_v1" :
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig3"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig4"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);

        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig3"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig4"%(self.categoryLabel, self.wtagger_label)));

        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig3"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig4"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig5"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig6"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig7"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);

        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig3"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig4"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig5"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig6"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig7"%(self.categoryLabel, self.wtagger_label)))


        #        if isTTbarFloating !=0 :
        #            self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #         self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #         self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #         self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig3"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);

        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)));
        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)));
        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)));
        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig3"%(self.categoryLabel, self.wtagger_label)));


        #    if self.MODEL_4_limit_variable=="Exp" or self.MODEL_4_limit_variable=="Pow" :

        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);

        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)));


        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);

        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)))

        #        if isTTbarFloating !=0 :
        #            self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #         params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)));

        #    if self.MODEL_4_limit_variable=="ExpN" or self.MODEL_4_limit_variable=="ExpTail" or self.MODEL_4_limit_variable=="Pow2" :

        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);
        #        self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_WJets0);

        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)));
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_lowersideband_from_fitting_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)));

        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);
        #        self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig3"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_alpha);

        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig2"%(self.categoryLabel, self.wtagger_label)))
        #        params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_limit_variable_eig3"%(self.categoryLabel, self.wtagger_label)))


        #        ### TTbar use exp
        #        if isTTbarFloating !=0:
        #            print "##################### TTbar will float in the limit procedure + final plot ######################";
        #            self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_TTbar);
        #            params_list.append(self.workspace4limit_.var("Deco_TTbar_signalregion_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)));

        #        ### VV use ExpTail:
        #        if isVVFloating !=0:
        #            print "##################### VV will float in the limit procedure + final plot ######################";
        #          self.workspace4limit_.var("Deco_VV_signalregion_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_VV);
        #          self.workspace4limit_.var("Deco_VV_signalregion_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_VV);
        #          params_list.append(self.workspace4limit_.var("Deco_VV_signalregion_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)));
        #          params_list.append(self.workspace4limit_.var("Deco_VV_signalregion_%s_%s_limit_variable_eig1"%(self.categoryLabel, self.wtagger_label)));

        #        ### SingleT use Exp:
        #        if isSingleTFloating !=0:
        #            print "##################### SingleT will float in the limit procedure + final plot ######################";
        #          self.workspace4limit_.var("Deco_SingleT_signalregion_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)).setError(self.shape_para_error_SingleT);
        #          params_list.append(self.workspace4limit_.var("Deco_SingleT_signalregion_%s_%s_limit_variable_eig0"%(self.categoryLabel, self.wtagger_label)));


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
        #self.rrv_counting_uncertainty_from_shape_uncertainty.setError( Calc_error("WJets_%s_%s"%(self.categoryLabel,self.wtagger_label), "rrv_limit_variable" ,self.FloatingParams,self.workspace4limit_,"signalregion") );
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

        if datacard_mode == "unbin":
            fnOnly = ntpath.basename(self.file_rlt_root) ## workspace for limit --> output file for the workspace
            datacard_out.write( "\nshapes * * %s %s:$PROCESS_%s"%(fnOnly, self.workspace4limit_.GetName(), self.categoryLabel))

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
        if datacard_mode == "unbin":
            datacard_out.write( "\nobservation %0.2f "%(self.workspace4limit_.data("data_obs_%s"%(self.categoryLabel)).sumEntries()) )
        if datacard_mode == "counting":
            datacard_out.write( "\nobservation %0.2f "%(self.workspace4limit_.var("observation_for_counting").getVal()) )


        datacard_out.write( "\n------------------------------" );

        tmp_bin_string="";
        tmp_processname_string="";
        tmp_processnum_string="";
        tmp_rate="";
        for iter in range( self.nsig ):
            tmp_bin_string        +="CMS_%s "%(self.categoryLabel); 
            tmp_processname_string+="%s "%(self.sig_list[iter][0])
            tmp_processnum_string +="%s "%(1-self.nsig+iter)
            if datacard_mode == "unbin":
                tmp_rate     +="%0.5f "%( self.workspace4limit_.var("rate_%s_for_unbin"%(self.sig_list[iter][0])).getVal() );
            if datacard_mode == "counting":
                tmp_rate     +="%0.5f "%( self.workspace4limit_.var("rate_%s_for_counting"%(self.sig_list[iter][0])).getVal() );

        for iter in range( self.nbkg ):
            tmp_bin_string        +="CMS_%s "%(self.categoryLabel); 
            tmp_processname_string+="%s "%(self.bkg_list[iter][0])
            tmp_processnum_string +="%s "%(1+iter)
            if datacard_mode == "unbin":
                tmp_rate     +="%0.5f "%( self.workspace4limit_.var("rate_%s_for_unbin"%(self.bkg_list[iter][0])).getVal() );
            if datacard_mode == "counting":
                tmp_rate     +="%0.5f "%( self.workspace4limit_.var("rate_%s_for_counting"%(self.bkg_list[iter][0])).getVal() );




        datacard_out.write( "\nbin "+tmp_bin_string );
        datacard_out.write( "\nprocess "+tmp_processname_string ); ## just one signal sample
        datacard_out.write( "\nprocess "+tmp_processnum_string );
        datacard_out.write( "\nrate "+tmp_rate );

        #### rates for the different process
        #if datacard_mode == "unbin":
        #    if TString(self.allsignals).Contains("BulkG_WW"):                    
        #        datacard_out.write( "\nrate %0.5f %0.3f %0.3f %0.3f %0.3f "%(self.workspace4limit_.var("rate_BulkWW_for_unbin").getVal()*self.xs_rescale, self.workspace4limit_.var("rate_WJets_for_unbin").getVal(), self.workspace4limit_.var("rate_TTbar_for_unbin").getVal(), self.workspace4limit_.var("rate_SingleT_for_unbin").getVal(), self.workspace4limit_.var("rate_VV_for_unbin").getVal() ) )
        #    else:
        #        datacard_out.write( "\nrate %0.5f %0.3f %0.3f %0.3f %0.3f "%(self.workspace4limit_.var("rate_%s_for_unbin"%(self.allsignals)).getVal()*self.xs_rescale, self.workspace4limit_.var("rate_WJets_for_unbin").getVal(), self.workspace4limit_.var("rate_TTbar_for_unbin").getVal(), self.workspace4limit_.var("rate_SingleT_for_unbin").getVal(), self.workspace4limit_.var("rate_VV_for_unbin").getVal() ) )


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
        #if ( self.workspace4fit_.var("rrv_number_WJets0_massup_in_obs0_variable_signalregion_from_fitting_%s"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_WJets0_massdn_in_obs0_variable_signalregion_from_fitting_%s"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT_massup_%s_obs0_variable"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT_massdown_%s_obs0_variable"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar_massup_%s_obs0_variable"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar_massdn_%s_obs0_variable"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_massup_%s_obs0_variable"%(self.categoryLabel)) and
        #        self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_massdn_%s_obs0_variable"%(self.categoryLabel))) :

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





    def DataDriven4WJetsNorm_fit_jetmass(self):
        print "DataDriven4WJets_fit_jetmass";
        self.fit_WJetsNormalization_in_Mj_signalregion();
         
         
    def SimpleScale(self):
        if self.categoryID <4:
            #el
            SF_WJets=1.75538; SF_WJets_err=0.1502;
            SF_TTbar=1.667;   SF_TTbar_err= 0.222;
        else:
            #mu
            SF_WJets=2.19802; SF_WJets_error=0.1349;
            SF_TTbar=1.306;   SF_TTbar_err=0.165;

        rrv_WJets_counting=self.workspace4fit_.var("rrv_number_fitting_signalregion_WJets_%s_limit_variable"%(self.categoryLabel));
        rrv_WJets_counting.setVal( rrv_WJets_counting.getVal()* SF_WJets );

        rrv_TTbar_counting=self.workspace4fit_.var("rrv_number_fitting_signalregion_TTbar_%s_limit_variable"%(self.categoryLabel));
        rrv_TTbar_counting.Print();
        rrv_TTbar_counting.setVal( rrv_TTbar_counting.getVal()* SF_TTbar  );
        rrv_TTbar_counting.Print();

        rrv_WJets_unbin=self.workspace4fit_.var("rrv_number_WJets_signalregion_%s_limit_variable"%(self.categoryLabel));
        rrv_WJets_unbin.setVal( rrv_WJets_unbin.getVal()* SF_WJets );

        rrv_TTbar_unbin=self.workspace4fit_.var("rrv_number_TTbar_signalregion_%s_limit_variable"%(self.categoryLabel));
        rrv_TTbar_unbin.setVal( rrv_TTbar_unbin.getVal()* SF_TTbar  );


        
    #### make the obs0_variable sideband fit on data to get the Wjets normaliztion 
    #def fit_WJetsNormalization_in_Mj_signalregion(self,label,massscale=""): 
    def fit_WJetsNormalization_in_Mj_signalregion(self): 

        print "############### Fit obs0_variable Normalization.  ##################"
        rrv_obs0_variable = self.workspace4fit_.var("rrv_obs0_variable")
        ## get dataset in obs0_variable distribution 
        rdataset_data_obs0_variable=self.workspace4fit_.data("rdataset_data_%s_obs0_variable"%(self.categoryLabel))


        model_pdf_bkgs=[];
        matrix_model_pdf_bkgs=[];
        rate_bkgs=[];
        ral_pdf_bkgs=RooArgList();#ral: RooArgList
        ral_rate_bkgs=RooArgList();
        rate_total_bkgs=0.;
        rate_error2_total_bkgs=0.;#error2: error^2
        model_WJets=0;  
        for iter in range(self.nbkg):
            pdf= self.workspace4fit_.pdf("model_%s_%s_obs0_variable"%(self.bkg_list[iter][0], self.categoryLabel));
            pdf.Print();
            fix_Pdf(pdf, RooArgSet(rrv_obs0_variable));# fix parameters' value
            model_pdf_bkgs.append(pdf);
            matrix_model_pdf_bkgs.append(pdf.GetName());
            for jter in range(iter):
                matrix_model_pdf_bkgs[jter]+=",%s"%(pdf.GetName());
            rrv_rate= self.workspace4fit_.var("rrv_number_%s_%s_obs0_variable"%(self.bkg_list[iter][0], self.categoryLabel));
            if not self.bkg_list[iter][0]=="WJets": rrv_rate.setConstant(kTRUE);
            else:
                rrv_rate.setRange(rrv_rate.getVal()/2., rrv_rate.getVal()*2); 
                rrv_rate.setConstant(0);
                model_WJets=pdf;
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
        rrv_number_Total_bkgs.Print();

        #### Total pdf 
        model_Total_bkgs = RooAddPdf("model_Total_bkgs","model_Total_bkgs", ral_pdf_bkgs, ral_rate_bkgs);
        model_Total_bkgs.Print();

        ## Total Pdf and fit only in sideband 
        rfresult = model_Total_bkgs.fitTo( rdataset_data_obs0_variable, RooFit.Save(1) , RooFit.Range("lowersideband,uppersideband") ,RooFit.Extended(kTRUE), RooFit.NumCPU(2) );
        rfresult = model_Total_bkgs.fitTo( rdataset_data_obs0_variable, RooFit.Save(1) , RooFit.Range("lowersideband,uppersideband") ,RooFit.Extended(kTRUE), RooFit.NumCPU(2), RooFit.Minimizer("Minuit2") );
        rfresult.Print(); rfresult.covarianceMatrix().Print();
        getattr(self.workspace4fit_,"import")(model_Total_bkgs);

        rate_total_bkgs=0.;
        rate_error2_total_bkgs=0.;#error2: error^2
        for iter in range(self.nbkg):
            rrv_rate= self.workspace4fit_.var("rrv_number_%s_%s_obs0_variable"%(self.bkg_list[iter][0], self.categoryLabel));
            rrv_rate.Print();
            rate_total_bkgs+=rrv_rate.getVal();
            rate_error2_total_bkgs+=rrv_rate.getError()*rrv_rate.getError();
        rrv_number_Total_bkgs = RooRealVar("rrv_number_Total_bkgs","rrv_number_Total_bkgs", rate_total_bkgs);
        rrv_number_Total_bkgs.setError( TMath.Sqrt( rate_error2_total_bkgs) );
        rrv_number_Total_bkgs.Print();

        #raw_input("ENTER");

        #### create the frame
        mplot = rrv_obs0_variable.frame(RooFit.Title("check_workspace"), RooFit.Bins(int(rrv_obs0_variable.getBins()/self.BinWidth_narrow_factor)));
        rdataset_data_obs0_variable.plotOn(mplot , RooFit.Name("data_invisible"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0));

        for iter in range(self.nbkg):
            model_Total_bkgs.plotOn(mplot, RooFit.Components(matrix_model_pdf_bkgs[iter]), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet[self.bkg_list[iter][0]]), RooFit.Name(self.bkg_list[iter][0]), RooFit.LineColor(kBlack), RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());
            model_Total_bkgs.plotOn(mplot, RooFit.Components(matrix_model_pdf_bkgs[iter]), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet[self.bkg_list[iter][0]]), RooFit.Name(self.bkg_list[iter][0]), RooFit.LineColor(kBlack), RooFit.FillStyle(3002), RooFit.Range(rrv_obs0_variable.getMin(),rrv_obs0_variable.getMax()),RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());

            model_Total_bkgs.plotOn(mplot, RooFit.Components(matrix_model_pdf_bkgs[iter]), RooFit.Name(self.bkg_list[iter][0]+"_line"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());
            model_Total_bkgs.plotOn(mplot, RooFit.Components(matrix_model_pdf_bkgs[iter]), RooFit.Name(self.bkg_list[iter][0]+"_line"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_obs0_variable.getMin(),rrv_obs0_variable.getMax()),RooFit.LineStyle(kDashed) ,RooFit.NormRange("lowersideband,uppersideband"), RooFit.VLines());


        #### plot the observed data using poissonian error bar
        #rdataset_data_obs0_variable.plotOn(mplot , RooFit.Name("data"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(1), RooFit.LineColor(0));
        self.getData_PoissonInterval(rrv_obs0_variable, rdataset_data_obs0_variable,mplot);

        model_Total_bkgs.plotOn(mplot,RooFit.Invisible());
        mplot_pull=self.get_pull(rrv_obs0_variable,mplot);

        ### Plot the list of floating parameters and the uncertainty band is draw taking into account this floating list defined in the prepare_limit
        draw_error_band(model_Total_bkgs, rrv_obs0_variable.GetName(), rrv_number_Total_bkgs,self.FloatingParams,self.workspace4fit_ ,mplot,self.color_palet["Uncertainty"],"F");

        mplot.Print();
        self.plot_legend = legend4Plot(mplot, self.categoryTitle,0,1,-0.01,-0.05,0.11,0.);
        self.plot_legend.SetTextSize(0.036);
        mplot.addObject(self.plot_legend);

        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);

        parameters_list = RooArgList();
        self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s/obs0_variable/"%(self.additioninformation, self.categoryLabel),"DataDriven4WJets_fit_jetmass","",0,1);

        #### call the function for getting the normalizatio in signal region for data, TTbar, SingleT, VV and W+jets = label -> store in a output txt file
        #self.get_obs0_variable_normalization_insignalregion("_data");
        #self.get_obs0_variable_normalization_insignalregion("_TTbar");
        #self.get_obs0_variable_normalization_insignalregion("_SingleT");
        #self.get_obs0_variable_normalization_insignalregion("_VV");
        #self.get_obs0_variable_normalization_insignalregion(label);

        ##### to calculate the WJets's normalization and error in M_J signalregion. The error must contain the shape error: model_WJets have new parameters fitting data
        fullInt   = model_WJets.createIntegral(RooArgSet(rrv_obs0_variable),RooArgSet(rrv_obs0_variable) );
        signalInt = model_WJets.createIntegral(RooArgSet(rrv_obs0_variable),RooArgSet(rrv_obs0_variable),("signalregion"));
        fullInt_val = fullInt.getVal()
        signalInt_val = signalInt.getVal()/fullInt_val
        ## take the value from the fit (normalization) and multiply it from the ratio of the integrals
        label="_WJets";
        rrv_number_WJets_in_obs0_variable_signalregion_from_fitting = RooRealVar("rrv_number%s_in_obs0_variable_signalregion_from_fitting_%s"%(label,self.categoryLabel), "rrv_number%s_in_obs0_variable_signalregion_from_fitting_%s"%(label,self.categoryLabel), self.workspace4fit_.var("rrv_number%s_%s_obs0_variable"%(label,self.categoryLabel)).getVal()*signalInt_val);

        ##### Error on the normalization --> from a dedicated function taking into account shape uncertainty
        rrv_number_WJets_in_obs0_variable_signalregion_from_fitting.setError( Calc_error_extendPdf(rdataset_data_obs0_variable, model_WJets, rfresult,"signalregion") );
        print "########## error on the normaliztion due to shape + norm = %s"%(rrv_number_WJets_in_obs0_variable_signalregion_from_fitting.getError());
        getattr(self.workspace4fit_,"import")(rrv_number_WJets_in_obs0_variable_signalregion_from_fitting);
        rrv_number_WJets_in_obs0_variable_signalregion_from_fitting.Print();

        rrv_number_WJets_in_obs0_variable_signalregion_old=self.workspace4fit_.var("rrv_number_WJets_signalregion"+"_"+self.categoryLabel+"_limit_variable");
        rrv_number_WJets_in_obs0_variable_signalregion_old.Print();

        rrv_WJets_counting=self.workspace4fit_.var("rrv_number_fitting_signalregion_WJets_%s_limit_variable"%(self.categoryLabel));
        rrv_WJets_counting.setVal( rrv_WJets_counting.getVal()* rrv_number_WJets_in_obs0_variable_signalregion_from_fitting.getVal()/rrv_number_WJets_in_obs0_variable_signalregion_old.getVal()  );

        rrv_number_WJets_in_obs0_variable_signalregion_old.setVal( rrv_number_WJets_in_obs0_variable_signalregion_from_fitting.getVal() );
        rrv_number_WJets_in_obs0_variable_signalregion_old.setError( rrv_number_WJets_in_obs0_variable_signalregion_from_fitting.getError() );


    #### Method to make a RooAbsPdf giving label, model name, spectrum, if it is mc or not and a constraint list for the parameters          
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

    ### get an obs0_variable model from the workspace givin the label
    def get_obs0_variable_Model(self,label):
        return self.workspace4fit_.pdf("model"+label+"_"+self.categoryLabel+"_obs0_variable")

    ### take the dataset, the model , the parameters in order to fix them as constant --> for extended pdf
    def get_General_obs0_variable_Model(self, label ):
        print "########### Fixing a general obs0_variable model  ############"
        rdataset_General_obs0_variable = self.workspace4fit_.data("rdataset%s_%s_obs0_variable"%(label,self.categoryLabel))
        model_General = self.get_obs0_variable_Model(label);
        rdataset_General_obs0_variable.Print();
        model_General.Print();
        ## get the parameters and cycle on them
        parameters_General = model_General.getParameters(rdataset_General_obs0_variable);
        par=parameters_General.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            if (TString(label).Contains("VV") or TString(label).Contains("SingleT") or TString(label).Contains("TTbar")):
                param.Print();
            param.setConstant(kTRUE);
            param=par.Next()
        ## return the pdf after having fixed the paramters
        return self.workspace4fit_.pdf("model%s_%s_obs0_variable"%(label,self.categoryLabel))

    ### fix only the ttbar component using the default label --> for extended pdf
    def get_TTbar_obs0_variable_Model(self,label="_TTbar"):
        print "########### Fixing only the TTbar obs0_variable Shape  ############"
        return self.get_General_obs0_variable_Model(label);

    ### fix only the stop component using the default label --> for extended pdf
    def get_SingleT_obs0_variable_Model(self,label="_SingleT"):
        print "########### Fixing only the Stop obs0_variable Shape  ############"
        return self.get_General_obs0_variable_Model(label);

    ### fix only the VV component using the default label --> for extended pdf
    def get_VV_obs0_variable_Model(self,label="_VV"):
        print "########### Fixing only the VV obs0_variable Shape  ############"
        return self.get_General_obs0_variable_Model(label);

    ### fix only the WJets model --> for extended pdf (just fix shape parameters of width, offset of ErfExp and p1 of User1 function
    def get_WJets_obs0_variable_Model(self,label):
        print "########### Fixing only the WJets obs0_variable Shape --> just the printed parameters  ############"
        rdataset_WJets_obs0_variable = self.workspace4fit_.data("rdataset%s_%s_obs0_variable"%(label,self.categoryLabel))
        model_WJets = self.get_obs0_variable_Model(label);
        rdataset_WJets_obs0_variable.Print();
        model_WJets.Print();
        parameters_WJets = model_WJets.getParameters(rdataset_WJets_obs0_variable);
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
        return self.workspace4fit_.pdf("model%s_%s_obs0_variable"%(label,self.categoryLabel))


    #### get a generic limit_variable model from the workspace
    def get_limit_variable_Model(self,label, limit_variable_region):
        return self.workspace4fit_.pdf("model"+label+limit_variable_region+"_"+self.categoryLabel+"_limit_variable");

    #### get a general limit_variable model and fiz the paramters --> for extended pdf
    def get_General_limit_variable_Model(self, label, limit_variable_region="_signalregion"):
        print "########### Fixing a general limit_variable model  ############"
        rdataset_General_limit_variable = self.workspace4fit_.data("rdataset%s%s_%s_limit_variable"%(label, limit_variable_region,self.categoryLabel))
        model_General = self.get_limit_variable_Model(label,limit_variable_region);
        rdataset_General_limit_variable.Print();
        model_General.Print();
        parameters_General = model_General.getParameters(rdataset_General_limit_variable);
        par=parameters_General.createIterator(); par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()
        return self.get_limit_variable_Model(label,limit_variable_region);

    ###### get TTbar model limit_variable in a region 
    def get_TTbar_limit_variable_Model(self, limit_variable_region="_signalregion"):
        print "########### Fixing TTbar limit_variable model for the region",limit_variable_region,"  ############"
        return self.get_General_limit_variable_Model("_TTbar",limit_variable_region);

    ###### get Single Top model limit_variable in a region 
    def get_SingleT_limit_variable_Model(self, limit_variable_region="_signalregion"):
        print "########### Fixing Stop limit_variable model for the region",limit_variable_region,"  ############"
        return self.get_General_limit_variable_Model("_SingleT",limit_variable_region);

    ###### get Signal model limit_variable in a region 
    def get_signal_limit_variable_Model(self, limit_variable_region="_signalregion"):
        print "########### Fixing signal limit_variable model for the region",limit_variable_region,"  ############"
        return self.get_General_limit_variable_Model("_%s"%(self.allsignals),limit_variable_region);

    ###### get VV limit_variable in a region 
    def get_VV_limit_variable_Model(self, limit_variable_region="_signalregion"):
        print "########### Fixing VV limit_variable for the region",limit_variable_region,"  ############"
        return self.get_General_limit_variable_Model("_VV",limit_variable_region);

    ###### get W+jets limit_variable in a region 
    def get_WJets_limit_variable_Model(self, limit_variable_region="_signalregion"):
        rdataset_WJets_limit_variable = self.workspace4fit_.data("rdataset_WJets_%s_limit_variable"%(limit_variable_region))
        model_WJets = self.get_limit_variable_Model("_WJets0",limit_variable_region);
        print "######## get Wjet limit_variable model for the region --> set constant just the normalization from obs0_variable fit",limit_variable_region," ########";
        rdataset_WJets_limit_variable.Print()
        model_WJets.Print()
        parameters_WJets = model_WJets.getParameters(rdataset_WJets_limit_variable);
        par = parameters_WJets.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
            paraName=TString(param.GetName());
            param.Print();
            if paraName.Contains("rrv_number_WJets"): ## set the correct normalization for W+jets if we are inside the signal region and fix it as constant
                if self.workspace4fit_.var("rrv_number_WJets_in_obs0_variable%s_from_fitting_%s"%(limit_variable_region,self.categoryLabel)):
                    self.workspace4fit_.var("rrv_number_WJets_in_obs0_variable%s_from_fitting_%s"%(limit_variable_region,self.categoryLabel)).Print()
                    param.setVal( self.workspace4fit_.var("rrv_number_WJets_in_obs0_variable%s_from_fitting_%s"%(limit_variable_region,self.categoryLabel)).getVal() )
                if limit_variable_region=="_signalregion": param.setConstant(kTRUE);
            param.Print();
            param=par.Next()
        return self.get_limit_variable_Model("_WJets0",limit_variable_region);


    ### Define the Extended Pdf for and limit_variable fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
    def fit_limit_variable_SingleChannel(self,in_file_name, label, in_range, limit_variable_model, deco=0, show_constant_parameter=0, logy=0, ismc=0):

        print "############### Fit limit_variable single MC sample ",in_file_name," ",label,"  ",limit_variable_model,"  ",in_range," ##################"
        ## import variable and dataset
        rrv_limit_variable = self.workspace4fit_.var("rrv_limit_variable")
        rdataset = self.workspace4fit_.data("rdataset4fit"+label+in_range+"_"+self.categoryLabel+"_limit_variable");
        constrainslist =[];

        ## make the extended pdf model
        model = self.make_Model(label+in_range,limit_variable_model,"_limit_variable",constrainslist,ismc);

        ## make the fit
        model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE) );
        rfresult = model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult.Print();

        ## set the name of the result of the fit and put it in the workspace   
        rfresult.SetName("rfresult"+label+in_range+"_"+self.categoryLabel+"_limit_variable")
        getattr(self.workspace4fit_,"import")(rfresult)

        ## plot the result
        mplot = rrv_limit_variable.frame(RooFit.Title("M_{lvj"+in_range+"} fitted by "+limit_variable_model), RooFit.Bins(int(rrv_limit_variable.getBins()/self.BinWidth_narrow_factor)));
        rdataset.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        ## plot the error band but don't store the canvas (only plotted without -b option
        draw_error_band_extendPdf(rdataset, model, rfresult,mplot,2,"L")
        rdataset.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        model.plotOn( mplot )#, RooFit.VLines()); in order to have the right pull 

        ## get the pull 
        mplot_pull      = self.get_pull(rrv_limit_variable,mplot);
        parameters_list = model.getParameters(rdataset);
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);

        self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s_%s_%s/m_lvj_fitting/"%(self.additioninformation, self.categoryLabel,self.PS_model,self.wtagger_label), in_file_name,"m_lvj"+in_range+limit_variable_model, show_constant_parameter, logy);


        ## if the shape parameters has to be decorrelated
        if deco :
            print "################### Decorrelated limit_variable single mc shape ################"
            model_pdf = self.workspace4fit_.pdf("model_pdf%s%s_%s_limit_variable"%(label,in_range,self.categoryLabel)); ## take the pdf from the workspace
            model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) );
            rfresult_pdf = model_pdf.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));
            rfresult_pdf.Print();

            ## temp workspace for the pdf diagonalizer
            wsfit_tmp = RooWorkspace("wsfit_tmp"+label+in_range+"_"+self.categoryLabel+"_limit_variable");
            Deco      = PdfDiagonalizer("Deco"+label+in_range+"_"+self.categoryLabel+"_"+self.wtagger_label+"_limit_variable",wsfit_tmp,rfresult_pdf); ## in order to have a good name 
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
            mplot_deco = rrv_limit_variable.frame( RooFit.Bins(int(rrv_limit_variable.getBins()/self.BinWidth_narrow_factor)));

            if label=="_TTbar" and in_range=="_signalregion":

                rdataset.plotOn(mplot_deco, RooFit.Name("Powheg Sample"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
                model_pdf_deco.plotOn(mplot_deco,RooFit.Name("TTbar_Powheg"),RooFit.LineColor(kBlack));

                mplot_deco.GetYaxis().SetRangeUser(1e-2,mplot_deco.GetMaximum()*1.2);

                rrv_number_dataset = RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset.sumEntries());
                rrv_number_dataset.setError(0.)
                draw_error_band(rdataset, model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,self.color_palet["Uncertainty"],"F"); ## draw the error band with the area
                self.workspace4fit_.var("rrv_number_TTbar_signalregion_%s_limit_variable"%(self.categoryLabel)).Print();
            else:
                rdataset.plotOn(mplot_deco, RooFit.Name("Data"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
                model_pdf_deco.plotOn(mplot_deco,RooFit.Name(label),RooFit.LineColor(kBlack));

                mplot_deco.GetYaxis().SetRangeUser(1e-2,mplot_deco.GetMaximum()*1.2);

                rrv_number_dataset=RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset.sumEntries());
                rrv_number_dataset.setError(0.)
                draw_error_band(rdataset, model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,self.color_palet["Uncertainty"],"F"); ## don't store the number in the workspace

            self.plot_legend = self.legend4Plot(mplot_deco,0); ## add the plot_legend                
            mplot_deco.addObject(self.plot_legend);

            self.draw_canvas( mplot_deco, "plots_%s_%s_%s_%s/other/"%(self.additioninformation, self.categoryLabel,self.PS_model, self.wtagger_label), "m_lvj"+label+in_range+in_range+limit_variable_model+"_deco",0,logy)

        ### Number of the event in the dataset and lumi scale factor --> set the proper number for bkg extraction or for signal region
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_limit_variable").Print();
        self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).Print()
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_limit_variable").setVal( self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_limit_variable").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).getVal() )
        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_limit_variable").setError(self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_limit_variable").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.categoryLabel).getVal() )

        self.workspace4fit_.var("rrv_number"+label+in_range+"_"+self.categoryLabel+"_limit_variable").Print();


    #### method to fit the WJets normalization inside the obs0_variable signal region -> and write the jets mass sys if available
    def fit_WJetsNorm(self, scaleJetMass = 0): # to get the normalization of WJets in signalregion

        print "############### Fit obs0_variable Normalization ##################"
        ## fit the two version of pdf for Wjets shape if available
        self.fit_WJetsNormalization_in_Mj_signalregion("_WJets0");
        self.fit_WJetsNormalization_in_Mj_signalregion("_WJets01");

        ## in case fit also the scaled jet mass distributions in order to have the jet mass scale sys included
        if scaleJetMass :
            self.fit_WJetsNormalization_in_Mj_signalregion("_WJets0_massup","massup");
         self.fit_WJetsNormalization_in_Mj_signalregion("_WJets0_massdn","massdn");
         self.fit_WJetsNormalization_in_Mj_signalregion("_WJets1");

        ## take the normalization numbers
        rrv_WJets0  = self.workspace4fit_.var("rrv_number_WJets0_in_obs0_variable_signalregion_from_fitting_%s"%(self.categoryLabel));
        rrv_WJets01 = self.workspace4fit_.var("rrv_number_WJets01_in_obs0_variable_signalregion_from_fitting_%s"%(self.categoryLabel));
        rrv_WJets0.Print();
        rrv_WJets01.Print();
        if scaleJetMass :
            rrv_WJets1 = self.workspace4fit_.var("rrv_number_WJets1_in_obs0_variable_signalregion_from_fitting_%s"%(self.categoryLabel));
         rrv_WJets1.Print();
         rrv_WJets0massup.Print();
         rrv_WJets0massdn.Print();

        ### total uncertainty combining the result with two different shapes
        total_uncertainty = TMath.Sqrt( TMath.Power(rrv_WJets0.getError(),2) + TMath.Power(rrv_WJets01.getVal()-rrv_WJets0.getVal(),2) );
        rrv_WJets0.setError(total_uncertainty);
        rrv_WJets0.Print();

        ##jet mass uncertainty on WJets normalization and the other bkg component
        if self.workspace4fit_.var("rrv_number_WJets0_massup_in_obs0_variable_signalregion_from_fitting_%s"%(self.categoryLabel)) and self.workspace4fit_.var("rrv_number_WJets0_massdn_in_obs0_variable_signalregion_from_fitting_%s"%(self.categoryLabel)):            
            rrv_WJets0massup = self.workspace4fit_.var("rrv_number_WJets0_massup_in_obs0_variable_signalregion_from_fitting_%s"%(self.categoryLabel));
          rrv_WJets0massdn = self.workspace4fit_.var("rrv_number_WJets0_massdn_in_obs0_variable_signalregion_from_fitting_%s"%(self.categoryLabel));
          self.WJets_normlization_uncertainty_from_jet_mass= ( TMath.Abs(rrv_WJets0massup.getVal()-rrv_WJets0.getVal())+TMath.Abs(rrv_WJets0massdn.getVal()-rrv_WJets0.getVal() ) )/2./rrv_WJets0.getVal();

        rrv_SingleT  = self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT__%s_obs0_variable"%(self.categoryLabel));

        if self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT_massup_%s_obs0_variable"%(self.categoryLabel)) and self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT_massdn_%s_obs0_variable"%(self.categoryLabel)) :
            rrv_SingleTmassup = self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT_massup_%s_obs0_variable"%(self.categoryLabel));
         rrv_SingleTmassdn = self.workspace4fit_.var("rrv_number_dataset_signalregion_SingleT_massdn_%s_obs0_variable"%(self.categoryLabel));
         self.SingleT_normlization_uncertainty_from_jet_mass=( TMath.Abs(rrv_SingleTmassup.getVal()-rrv_SingleT.getVal())+TMath.Abs(rrv_SingleTmassdn.getVal()-rrv_SingleT.getVal() ) )/2./rrv_SingleT.getVal();

        rrv_TTbar = self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar__%s_obs0_variable"%(self.categoryLabel));
        if self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar_massup_%s_obs0_variable"%(self.categoryLabel)) and self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar_massdn_%s_obs0_variable"%(self.categoryLabel)):
            rrv_TTbarmassup = self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar_massup_%s_obs0_variable"%(self.categoryLabel));
         rrv_TTbarmassdn = self.workspace4fit_.var("rrv_number_dataset_signalregion_TTbar_massdn_%s_obs0_variable"%(self.categoryLabel));
         self.TTbar_normlization_uncertainty_from_jet_mass=( TMath.Abs(rrv_TTbarmassup.getVal()-rrv_TTbar.getVal())+TMath.Abs(rrv_TTbarmassdn.getVal()-rrv_TTbar.getVal() ) )/2./rrv_TTbar.getVal();

        rrv_VV = self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_%s_obs0_variable"%(self.categoryLabel));
        if self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_massup_%s_obs0_variable"%(self.categoryLabel)) and self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_massdn_%s_obs0_variable"%(self.categoryLabel)):
            rrv_VVmassup = self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_massup_%s_obs0_variable"%(self.categoryLabel));
         rrv_VVmassdn = self.workspace4fit_.var("rrv_number_dataset_signalregion_VV_massdn_%s_obs0_variable"%(self.categoryLabel));
         self.VV_normlization_uncertainty_from_jet_mass=( TMath.Abs(rrv_VVmassup.getVal()-rrv_VV.getVal())+TMath.Abs(rrv_VVmassdn.getVal()-rrv_VV.getVal() ) )/2./rrv_VV.getVal();


    ##### Counting of the events of each component in the signal region taking the lavel for the model
    def get_obs0_variable_normalization_insignalregion(self, label):
        print "################## get obs0_variable normalization ",label," ################## ";
        rrv_obs0_variable = self.workspace4fit_.var("rrv_obs0_variable");
        model      = self.workspace4fit_.pdf("model"+label+"_"+self.categoryLabel+"_obs0_variable");

        fullInt   = model.createIntegral(RooArgSet(rrv_obs0_variable),RooArgSet(rrv_obs0_variable) );
        lowersidebandInt  = model.createIntegral(RooArgSet(rrv_obs0_variable),RooArgSet(rrv_obs0_variable),("lowersideband"));
        signalInt = model.createIntegral(RooArgSet(rrv_obs0_variable),RooArgSet(rrv_obs0_variable),("signalregion"));
        uppersidebandInt  = model.createIntegral(RooArgSet(rrv_obs0_variable),RooArgSet(rrv_obs0_variable),("uppersideband"));

        fullInt_val   = fullInt.getVal()
        lowersidebandInt_val  = lowersidebandInt.getVal()/fullInt_val
        uppersidebandInt_val  = uppersidebandInt.getVal()/fullInt_val
        signalInt_val = signalInt.getVal()/fullInt_val

        print "########### Events Number in MC Dataset: #############"
        self.workspace4fit_.var("rrv_number_dataset_lowersideband"+label+"_"+self.categoryLabel+"_obs0_variable").Print();
        self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_obs0_variable").Print();
        self.workspace4fit_.var("rrv_number_dataset_uppersideband"+label+"_"+self.categoryLabel+"_obs0_variable").Print();

        print "########### Events Number get from fit: ##############"
        rrv_tmp = self.workspace4fit_.var("rrv_number"+label+"_"+self.categoryLabel+"_obs0_variable");
        rrv_tmp.Print();
        print "Events Number in sideband_low :%s"%(rrv_tmp.getVal()*lowersidebandInt_val)
        print "Events Number in Signal Region:%s"%(rrv_tmp.getVal()*signalInt_val)
        print "Events Number in sideband_high:%s"%(rrv_tmp.getVal()*uppersidebandInt_val)
        print "Total Number in sidebands :%s"%(rrv_tmp.getVal()*(lowersidebandInt_val+uppersidebandInt_val) )
        print "Ratio signalregion/sidebands :%s"%(signalInt_val/(lowersidebandInt_val+uppersidebandInt_val) )

        ##### Save numbers in the output text file
        self.file_out.write( "\n%s++++++++++++++++++++++++++++++++++++"%(label) )
        self.file_out.write( "\nEvents Number in sideband_low from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_lowersideband"+label+"_"+self.categoryLabel+"_obs0_variable").getVal() ) )
        self.file_out.write( "\nEvents Number in Signal Region from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_obs0_variable").getVal() ) )
        self.file_out.write( "\nEvents Number in sideband_high from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_uppersideband"+label+"_"+self.categoryLabel+"_obs0_variable").getVal() ) )
        self.file_out.write( "\nTotal Number in sidebands from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_lowersideband"+label+"_"+self.categoryLabel+"_obs0_variable").getVal()+ self.workspace4fit_.var("rrv_number_dataset_uppersideband"+label+"_"+self.categoryLabel+"_obs0_variable").getVal() ) )
        self.file_out.write( "\nRatio signalregion/sidebands from dataset:%s"%(self.workspace4fit_.var("rrv_number_dataset_signalregion"+label+"_"+self.categoryLabel+"_obs0_variable").getVal()/(self.workspace4fit_.var("rrv_number_dataset_lowersideband"+label+"_"+self.categoryLabel+"_obs0_variable").getVal()+ self.workspace4fit_.var("rrv_number_dataset_uppersideband"+label+"_"+self.categoryLabel+"_obs0_variable").getVal()) ) )

        self.file_out.write( "\nEvents Number in sideband_low from fitting:%s"%(rrv_tmp.getVal()*lowersidebandInt_val) )
        self.file_out.write( "\nEvents Number in Signal Region from fitting:%s"%(rrv_tmp.getVal()*signalInt_val) )
        self.file_out.write( "\nEvents Number in sideband_high from fitting:%s"%(rrv_tmp.getVal()*uppersidebandInt_val) )
        self.file_out.write( "\nTotal Number in sidebands from fitting:%s"%(rrv_tmp.getVal()*(lowersidebandInt_val+uppersidebandInt_val) ) )
        self.file_out.write( "\nRatio signalregion/sidebands from fitting:%s"%(signalInt_val/(lowersidebandInt_val+uppersidebandInt_val) ) )

    ##### Method to fit data limit_variable shape in the sideband -> first step for the background extraction of the shape
    def fit_limit_variable_in_Mj_sideband(self, label, limit_variable_region, limit_variable_model,logy=0):

        print "############### Fit limit_variable in obs0_variable sideband: ",label," ",limit_variable_region,"  ",limit_variable_model," ##################"
        rrv_obs0_variable   = self.workspace4fit_.var("rrv_obs0_variable")
        rrv_limit_variable = self.workspace4fit_.var("rrv_limit_variable")
        rdataset_data_limit_variable = self.workspace4fit_.data("rdataset_data%s_%s_limit_variable"%(limit_variable_region,self.categoryLabel))

        ## get the minor component shapes in the sb low
        model_VV_backgrounds    = self.get_VV_limit_variable_Model("_lowersideband");
        number_VV_lowersideband_limit_variable    = self.workspace4fit_.var("rrv_number_VV_lowersideband_%s_limit_variable"%(self.categoryLabel))
        model_TTbar_backgrounds = self.get_TTbar_limit_variable_Model("_lowersideband");
        number_TTbar_lowersideband_limit_variable = self.workspace4fit_.var("rrv_number_TTbar_lowersideband_%s_limit_variable"%(self.categoryLabel))
        model_SingleT_backgrounds  = self.get_SingleT_limit_variable_Model("_lowersideband");
        number_SingleT_lowersideband_limit_variable  = self.workspace4fit_.var("rrv_number_SingleT_lowersideband_%s_limit_variable"%(self.categoryLabel))

        self.workspace4fit_.var("rrv_number_TTbar_lowersideband_%s_limit_variable"%(self.categoryLabel)).Print();
        self.workspace4fit_.var("rrv_number_SingleT_lowersideband_%s_limit_variable"%(self.categoryLabel)).Print();
        self.workspace4fit_.var("rrv_number_VV_lowersideband_%s_limit_variable"%(self.categoryLabel)).Print();

        ### Make the Pdf for the WJets
        model_pdf_WJets = self.make_Pdf("%s_lowersideband_from_fitting"%(label), limit_variable_model,"_limit_variable");
        model_pdf_WJets.Print();
        ### inititalize the value to what was fitted with the mc in the sideband
        number_WJets_lowersideband = self.workspace4fit_.var("rrv_number%s_lowersideband_%s_limit_variable"%(label,self.categoryLabel)).clone("rrv_number%s_lowersideband_from_fitting_%s_limit_variable"%(label,self.categoryLabel));
        model_WJets =RooExtendPdf("model%s_lowersideband_from_fitting_%s_limit_variable"%(label,self.categoryLabel),"model%s_lowersideband_from_fitting_%s_limit_variable"%(label,self.categoryLabel),model_pdf_WJets,number_WJets_lowersideband);
        model_pdf_WJets.Print();
        number_WJets_lowersideband.Print()

        ## Add the other bkg component fixed to the total model
        model_data = RooAddPdf("model_data%s%s_limit_variable"%(label,limit_variable_region),"model_data%s%s_limit_variable"%(label,limit_variable_region),RooArgList(model_WJets,model_VV_backgrounds, model_TTbar_backgrounds, model_SingleT_backgrounds));

        rfresult = model_data.fitTo( rdataset_data_limit_variable, RooFit.Save(1) ,RooFit.Extended(kTRUE));
        rfresult = model_data.fitTo( rdataset_data_limit_variable, RooFit.Save(1) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"));
        rfresult.Print();
        rfresult.covarianceMatrix().Print();
        getattr(self.workspace4fit_,"import")(model_data)

        model_WJets.Print();
        model_WJets.getParameters(rdataset_data_limit_variable).Print("v");
        self.workspace4fit_.pdf("model_pdf%s_lowersideband_%s_limit_variable"%(label,self.categoryLabel)).getParameters(rdataset_data_limit_variable).Print("v");

        ### data in the sideband plus error from fit
        rrv_number_data_lowersideband_limit_variable = RooRealVar("rrv_number_data_lowersideband_%s_limit_variable"%(self.categoryLabel),"rrv_number_data_lowersideband_%s_limit_variable"%(self.categoryLabel),
                self.workspace4fit_.var("rrv_number_TTbar_lowersideband_%s_limit_variable"%(self.categoryLabel)).getVal()+
                self.workspace4fit_.var("rrv_number_SingleT_lowersideband_%s_limit_variable"%(self.categoryLabel)).getVal()+
                self.workspace4fit_.var("rrv_number_VV_lowersideband_%s_limit_variable"%(self.categoryLabel)).getVal()+
                self.workspace4fit_.var("rrv_number%s_lowersideband_from_fitting_%s_limit_variable"%(label,self.categoryLabel)).getVal() );

        rrv_number_data_lowersideband_limit_variable.setError( TMath.Sqrt(self.workspace4fit_.var("rrv_number%s_lowersideband_from_fitting_%s_limit_variable"%(label,self.categoryLabel)).getError()*
            self.workspace4fit_.var("rrv_number%s_lowersideband_from_fitting_%s_limit_variable"%(label,self.categoryLabel)).getError()+
            self.workspace4fit_.var("rrv_number_TTbar_lowersideband_%s_limit_variable"%(self.categoryLabel)).getError()*
            self.workspace4fit_.var("rrv_number_TTbar_lowersideband_%s_limit_variable"%(self.categoryLabel)).getError()+
            self.workspace4fit_.var("rrv_number_SingleT_lowersideband_%s_limit_variable"%(self.categoryLabel)).getError()*
            self.workspace4fit_.var("rrv_number_SingleT_lowersideband_%s_limit_variable"%(self.categoryLabel)).getError()+
            self.workspace4fit_.var("rrv_number_VV_lowersideband_%s_limit_variable"%(self.categoryLabel)).getError()*
            self.workspace4fit_.var("rrv_number_VV_lowersideband_%s_limit_variable"%(self.categoryLabel)).getError()));

        getattr(self.workspace4fit_,"import")(rrv_number_data_lowersideband_limit_variable)

        ### plot for WJets default + default shape
        if TString(label).Contains("_WJets0"):

            mplot = rrv_limit_variable.frame(RooFit.Title("M_lvj fitted in M_j sideband "), RooFit.Bins(int(rrv_limit_variable.getBins()/self.BinWidth_narrow_factor)));

            rdataset_data_limit_variable.plotOn( mplot , RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0), RooFit.MarkerColor(0), RooFit.LineColor(0) );

            model_data.plotOn(mplot, RooFit.Components("model%s_lowersideband_from_fitting_%s_limit_variable,model_TTbar_lowersideband_%s_limit_variable,model_SingleT_lowersideband_%s_limit_variable,model_VV_lowersideband_%s_limit_variable"%(label,self.categoryLabel,self.categoryLabel,self.categoryLabel,self.categoryLabel)), RooFit.Name("WJets"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_lowersideband_%s_limit_variable,model_SingleT_lowersideband_%s_limit_variable,model_VV_lowersideband_%s_limit_variable"%(self.categoryLabel,self.categoryLabel,self.categoryLabel)),RooFit.Name("VV"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_lowersideband_%s_limit_variable,model_SingleT_lowersideband_%s_limit_variable"%(self.categoryLabel,self.categoryLabel)), RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components("model_SingleT_lowersideband_%s_limit_variable"%(self.categoryLabel)), RooFit.Name("SingleT"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["SingleT"]), RooFit.LineColor(kBlack), RooFit.VLines());

            #solid line
            model_data.plotOn(mplot, RooFit.Components("model%s_lowersideband_from_fitting_%s_limit_variable,model_TTbar_lowersideband_%s_limit_variable,model_SingleT_lowersideband_%s_limit_variable,model_VV_lowersideband_%s_limit_variable"%(label,self.categoryLabel,self.categoryLabel,self.categoryLabel,self.categoryLabel)), RooFit.Name("WJets_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_lowersideband_%s_limit_variable,model_SingleT_lowersideband_%s_limit_variable,model_VV_lowersideband_%s_limit_variable"%(self.categoryLabel,self.categoryLabel,self.categoryLabel)),RooFit.Name("VV_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines()) ;

            model_data.plotOn(mplot, RooFit.Components("model_TTbar_lowersideband_%s_limit_variable,model_SingleT_lowersideband_%s_limit_variable"%(self.categoryLabel,self.categoryLabel)), RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components("model_SingleT_lowersideband_%s_limit_variable"%(self.categoryLabel)), RooFit.Name("SingleT_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());


            ### draw the error band 
            draw_error_band(rdataset_data_limit_variable, model_data,self.workspace4fit_.var("rrv_number_data_lowersideband_%s_limit_variable"%(self.categoryLabel)) ,rfresult,mplot,self.color_palet["Uncertainty"],"F");
            model_data.plotOn( mplot , RooFit.VLines(), RooFit.Invisible());
            model_data.plotOn( mplot , RooFit.Invisible());
            self.getData_PoissonInterval(rdataset_data_limit_variable,mplot);

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
            self.file_out.write("\n fit_limit_variable_in_Mj_sideband: nPar=%s, chiSquare=%s/%s"%(self.nPar_float_in_fitTo, mplot.chiSquare( self.nPar_float_in_fitTo )*ndof, ndof ) );

            ### get the pull plot and store the canvas
            mplot_pull = self.get_pull(rrv_limit_variable,mplot);
            parameters_list = model_data.getParameters(rdataset_data_limit_variable);

            self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s_%s_%s/m_lvj_fitting/"%(self.additioninformation, self.categoryLabel,self.PS_model,self.wtagger_label), "m_lvj_lowersideband%s"%(label),"",1,1)

        #### Decorrelate the parameters in order to have a proper shape in the workspace
        wsfit_tmp = RooWorkspace("wsfit_tmp%s_lowersideband_from_fitting_limit_variable"%(label));
        Deco      = PdfDiagonalizer("Deco%s_lowersideband_from_fitting_%s_%s_limit_variable"%(label,self.categoryLabel,self.wtagger_label),wsfit_tmp,rfresult);
        print"#################### diagonalize data sideband fit "
        model_pdf_WJets_deco = Deco.diagonalize(model_pdf_WJets);
        print"#################### print parameters "
        model_pdf_WJets_deco.Print("v");
        model_pdf_WJets_deco.getParameters(rdataset_data_limit_variable).Print("");
        getattr(self.workspace4fit_,"import")(model_pdf_WJets_deco);

        #### Call the alpha evaluation in automatic
        self.get_WJets_limit_variable_correction_lowersideband_to_signalregion(label,limit_variable_model);

        ### Fix the pdf of signal, TTbar, SingleT and VV in the signal region 
        self.fix_Model("_%s"%(self.allsignals),"_signalregion","_limit_variable")
        self.fix_Model("_TTbar","_signalregion","_limit_variable")
        self.fix_Model("_SingleT","_signalregion","_limit_variable")
        self.fix_Model("_VV","_signalregion","_limit_variable")

        ### Call the evaluation of the normalization in the signal region for signal, TTbar, VV, SingleT, and WJets after the extrapolation via alpha
        self.get_pdf_signalregion_integral("_%s"%(self.allsignals));
        self.get_pdf_signalregion_integral("_TTbar");
        self.get_pdf_signalregion_integral("_SingleT");
        self.get_pdf_signalregion_integral("_VV");
        self.get_pdf_signalregion_integral(label,"model_pdf%s_signalregion_%s_after_correct_limit_variable"%(label,self.categoryLabel));    


    ### method to get the alpha function to extrapolate the wjets in the signal region
    def get_WJets_limit_variable_correction_lowersideband_to_signalregion(self,label, limit_variable_model):

        print" ############# get the extrapolation function alpha from MC : ",label,"   ",limit_variable_model," ###############";          
        tmp_Style = self.tdrStyle.Clone("tmp_Style");
        tmp_Style.SetPadRightMargin(0.08);
        tmp_Style.SetPadTickY(0);
        tmp_Style.cd();

        ### take input var and datasets from 4fit collection --> mc not scaled to lumi --> just a shape here 
        rrv_x = self.workspace4fit_.var("rrv_limit_variable");
        rdataset_WJets_lowersideband_limit_variable = self.workspace4fit_.data("rdataset4fit%s_lowersideband_%s_limit_variable"%(label,self.categoryLabel))
        rdataset_WJets_signalregion_limit_variable = self.workspace4fit_.data("rdataset4fit%s_signalregion_%s_limit_variable"%(label,self.categoryLabel))

        ### create a frame for the next plots 
        mplot = rrv_x.frame(RooFit.Title("correlation_pdf"), RooFit.Bins(int(rrv_x.getBins()/self.BinWidth_narrow_factor))) ;
        mplot.GetYaxis().SetTitle("arbitrary units");

        ### model used for Higgs analysis --> parameters in the SR has to be fitted, not yet done in order to take into account correlations between obs0_variable and limit_variable
        if limit_variable_model=="ErfExp_v1":

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

        if limit_variable_model=="ErfPow_v1":

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

        if limit_variable_model=="ErfPow2_v1":

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

        if limit_variable_model=="ErfPowExp_v1": ## take initial value from what was already fitted in the SR

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

        if limit_variable_model=="Exp":
            rrv_c_sb    = self.workspace4fit_.var("rrv_c_Exp%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_delta_c = RooRealVar("rrv_delta_c_Exp%s_%s"%(label,self.categoryLabel),"rrv_delta_c_Exp%s_%s"%(label,self.categoryLabel),
                    self.workspace4fit_.var("rrv_c_Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c_sb.getVal(),
                    self.workspace4fit_.var("rrv_c_Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c_sb.getVal()-4*rrv_c_sb.getError(),
                    self.workspace4fit_.var("rrv_c_Exp%s_signalregion_%s"%(label,self.categoryLabel)).getVal()-rrv_c_sb.getVal()+4*rrv_c_sb.getError() )

            correct_factor_pdf = RooExponential("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c);

        if limit_variable_model=="2Exp":
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

        if limit_variable_model=="Pow":

            rrv_c_sb    = self.workspace4fit_.var("rrv_c_Pow%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_delta_c = RooRealVar("rrv_delta_c_Pow%s_%s"%(label,self.categoryLabel),"rrv_delta_c_Pow%s_%s"%(label,self.categoryLabel),0., -100*rrv_c_sb.getError(),100*rrv_c_sb.getError());
            correct_factor_pdf = RooPowPdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c);

        if limit_variable_model=="ExpN":
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

        if limit_variable_model=="ExpTail":
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

        if limit_variable_model=="Pow2":

            rrv_c0_sb    = self.workspace4fit_.var("rrv_c0_Pow2%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_c1_sb    = self.workspace4fit_.var("rrv_c1_Pow2%s_lowersideband_%s"%(label,self.categoryLabel));
            rrv_delta_c0 = RooRealVar("rrv_delta_c0_Pow2%s_%s"%(label,self.categoryLabel),"rrv_delta_c0_Pow2%s_%s"%(label,self.categoryLabel),0., -100*rrv_c0_sb.getError(),100*rrv_c0_sb.getError());
            rrv_delta_c1 = RooRealVar("rrv_delta_c1_Pow2%s_%s"%(label,self.categoryLabel),"rrv_delta_c1_Pow2%s_%s"%(label,self.categoryLabel),0., -100*rrv_c1_sb.getError(),100*rrv_c1_sb.getError());
            correct_factor_pdf = RooPow2Pdf("correct_factor_pdf","correct_factor_pdf",rrv_x,rrv_delta_c0,rrv_delta_c1);

        ### define the category and do the simultaneous fit taking the combined dataset of events in limit_variable sb and sr

        data_category = RooCategory("data_category","data_category");
        data_category.defineType("sideband");
        data_category.defineType("signalregion");
        combData4fit = self.workspace4fit_.data("combData4fit%s_%s"%(label,self.categoryLabel));

        model_pdf_lowersideband_WJets         = self.workspace4fit_.pdf("model_pdf%s_lowersideband_%s_limit_variable"%(label,self.categoryLabel));
        model_pdf_signalregion_WJets = RooProdPdf("model_pdf%s_signalregion_%s_limit_variable"%(label,self.categoryLabel),"model_pdf%s_signalregion_%s_limit_variable"%(label,self.categoryLabel) ,model_pdf_lowersideband_WJets,correct_factor_pdf);

        simPdf = RooSimultaneous("simPdf","simPdf",data_category);
        simPdf.addPdf(model_pdf_lowersideband_WJets,"sideband");
        simPdf.addPdf(model_pdf_signalregion_WJets,"signalregion");
        rfresult=simPdf.fitTo(combData4fit,RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE));
        rfresult=simPdf.fitTo(combData4fit,RooFit.Save(kTRUE), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));
        rfresult.Print();
        rfresult.covarianceMatrix().Print();

        ### Decorrelate the parameters in the alpha shape
        wsfit_tmp = RooWorkspace("wsfit_tmp%s_sim_limit_variable"%(label));
        print "############### diagonalizer alpha ";
        Deco      = PdfDiagonalizer("Deco%s_sim_%s_%s_limit_variable"%(label,self.categoryLabel,self.wtagger_label),wsfit_tmp,rfresult);
        correct_factor_pdf_deco = Deco.diagonalize(correct_factor_pdf);
        correct_factor_pdf_deco.Print();
        correct_factor_pdf_deco.getParameters(rdataset_WJets_signalregion_limit_variable).Print("v");
        getattr(self.workspace4fit_,"import")(correct_factor_pdf_deco);

        ## in case of default Wjets with default shape
        if TString(label).Contains("_WJets0"):

            ### only mc plots in the SB region
            mplot_lowersideband = rrv_x.frame(RooFit.Title("WJets sb low"), RooFit.Bins(int(rrv_x.getBins()/self.BinWidth_narrow_factor)));

            rdataset_WJets_lowersideband_limit_variable.plotOn(mplot_lowersideband, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
            model_pdf_lowersideband_WJets.plotOn(mplot_lowersideband);
            mplot_pull_sideband = self.get_pull(rrv_x,mplot_lowersideband);
            parameters_list     = model_pdf_lowersideband_WJets.getParameters(rdataset_WJets_lowersideband_limit_variable);
            mplot_lowersideband.GetYaxis().SetRangeUser(1e-2,mplot_lowersideband.GetMaximum()*1.2);
            self.draw_canvas_with_pull( mplot_lowersideband, mplot_pull_sideband,parameters_list,"plots_%s_%s_%s_%s/other/"%(self.additioninformation, self.categoryLabel,self.PS_model,self.wtagger_label), "m_lvj%s_lowersideband_sim"%(label),"",1,1)

            ### only mc plots in the SR region
            mplot_signalregion = rrv_x.frame(RooFit.Title("WJets sr"), RooFit.Bins(int(rrv_x.getBins()/self.BinWidth_narrow_factor)));

            rdataset_WJets_signalregion_limit_variable.plotOn(mplot_signalregion, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
            model_pdf_signalregion_WJets.plotOn(mplot_signalregion);
            mplot_pull_signalregion = self.get_pull(rrv_x, mplot_signalregion);
            parameters_list = model_pdf_signalregion_WJets.getParameters(rdataset_WJets_signalregion_limit_variable);
            mplot_signalregion.GetYaxis().SetRangeUser(1e-2,mplot_signalregion.GetMaximum()*1.2);
            self.draw_canvas_with_pull( mplot_signalregion, mplot_pull_signalregion,parameters_list,"plots_%s_%s_%s_%s/other/"%(self.additioninformation, self.categoryLabel,self.PS_model,self.wtagger_label), "m_lvj%s_signalregion_sim"%(label),"",1,1);

        ### Total plot shape in lowersideband, sr and alpha
        model_pdf_lowersideband_WJets.plotOn(mplot,RooFit.Name("Sideband"),RooFit.LineStyle(10));
        model_pdf_signalregion_WJets.plotOn(mplot, RooFit.LineColor(kRed) ,RooFit.LineStyle(8), RooFit.Name("Signal Region"));
        correct_factor_pdf_deco.plotOn(mplot, RooFit.LineColor(kBlack),RooFit.Name("#alpha") );

        ### plot also what is get from other source if available : alternate PS and shape: 1 PS and 01 is shape or fitting function
        if TString(label).Contains("_WJets0"):
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_limit_variable"%(self.categoryLabel, self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_limit_variable"%(self.categoryLabel, self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha: Alternate PS") );

            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_sim_%s_%s_limit_variable"%(self.categoryLabel, self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_sim_%s_%s_limit_variable"%(self.categoryLabel, self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7),RooFit.Name("#alpha: Alternate Function") );

        paras=RooArgList();
        ### Make a list of paramters as a function of the model after decorrelation 
        if limit_variable_model=="ErfExp_v1" or limit_variable_model=="ErfPow_v1" or limit_variable_model=="2Exp" :
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig0"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig1"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig2"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig3"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig4"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig5"%(label,self.categoryLabel, self.wtagger_label) ));

        if limit_variable_model=="ErfPow2_v1" or limit_variable_model=="ErfPowExp_v1" :
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig0"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig1"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig2"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig3"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig4"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig5"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig6"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig7"%(label,self.categoryLabel, self.wtagger_label) ));

        if limit_variable_model=="Exp" or limit_variable_model=="Pow":
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig0"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig1"%(label,self.categoryLabel, self.wtagger_label) ));

        if limit_variable_model=="ExpN" or limit_variable_model=="ExpTail" or limit_variable_model=="Pow2":
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig0"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig1"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig2"%(label,self.categoryLabel, self.wtagger_label) ));
            paras.add(self.workspace4fit_.var("Deco%s_sim_%s_%s_limit_variable_eig3"%(label,self.categoryLabel, self.wtagger_label) ));

        if TString(label).Contains("_WJets0") or TString(label).Contains("_WJets1"): ### draw error band ar 1 and 2 sigma using the decorrelated shape
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_limit_variable"%(label,self.categoryLabel, self.wtagger_label),"rrv_limit_variable", paras, self.workspace4fit_,1 ,mplot,kGray+3,"F",3001,"#alpha #pm",20,400);
            draw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_limit_variable"%(label,self.categoryLabel, self.wtagger_label),"rrv_limit_variable", paras, self.workspace4fit_,2 ,mplot,kGreen+2,"F",3002,"#alpha #pm",20,400);
            `raw_error_band_shape_Decor("correct_factor_pdf_Deco%s_sim_%s_%s_limit_variable"%(label,self.categoryLabel, self.wtagger_label),"rrv_limit_variable", paras, self.workspace4fit_,1 ,mplot,kGray+3,"F",3001,"#alpha_invisible #pm",20,400);

        ### plot on the same canvas
        correct_factor_pdf_deco.plotOn(mplot, RooFit.LineColor(kBlack),RooFit.Name("#alpha_invisible") );

        if TString(label).Contains("_WJets0") : ## add also the plot of alternate ps and function on the canvas
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_limit_variable"%(self.categoryLabel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_limit_variable"%(self.categoryLabel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha_invisible: Alternate PS") );
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_sim_%s_%s_limit_variable"%(self.categoryLabel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets01_sim_%s_%s_limit_variable"%(self.categoryLabel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7),RooFit.Name("#alpha_invisible: Alternate Function") );

        elif TString(label).Contains("_WJets01") : ## add also the plot of alternate ps and function on the canvas
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_limit_variable"%(self.categoryLabel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets1_sim_%s_%s_limit_variable"%(self.categoryLabel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kMagenta), RooFit.LineStyle(3),RooFit.Name("#alpha_invisible: Alternate PS") );
            if self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets0_sim_%s_%s_limit_variable"%(self.categoryLabel,self.wtagger_label)):
                self.workspace4fit_.pdf("correct_factor_pdf_Deco_WJets0_sim_%s_%s_limit_variable"%(self.categoryLabel,self.wtagger_label)).plotOn(mplot, RooFit.LineColor(kOrange), RooFit.LineStyle(7),RooFit.Name("#alpha_invisible: Alternate Function") );

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

        self.draw_canvas(mplot,"plots_%s_%s_%s_%s/other/"%(self.additioninformation, self.categoryLabel,self.PS_model,self.wtagger_label),"correction_pdf%s_%s_%s_M_lvj_signalregion_to_sideband"%(label,self.PS_model,limit_variable_model),0,1);

        correct_factor_pdf_deco.getParameters(rdataset_WJets_lowersideband_limit_variable).Print("v");
        model_pdf_WJets_lowersideband_from_fitting_limit_variable_Deco = self.workspace4fit_.pdf("model_pdf%s_lowersideband_from_fitting_%s_limit_variable_Deco%s_lowersideband_from_fitting_%s_%s_limit_variable"%(label,self.categoryLabel,label, self.categoryLabel,self.wtagger_label));
        model_pdf_WJets_lowersideband_from_fitting_limit_variable_Deco.Print("v");

        ### Wjets shape in the SR correctedfunction * sb 
        model_pdf_WJets_signalregion_after_correct_limit_variable = RooProdPdf("model_pdf%s_signalregion_%s_after_correct_limit_variable"%(label,self.categoryLabel),"model_pdf%s_signalregion_%s_after_correct_limit_variable"%(label,self.categoryLabel),model_pdf_WJets_lowersideband_from_fitting_limit_variable_Deco,self.workspace4fit_.pdf("correct_factor_pdf_Deco%s_sim_%s_%s_limit_variable"%(label,self.categoryLabel,self.wtagger_label)) );
        model_pdf_WJets_signalregion_after_correct_limit_variable.Print()
        ### fix the parmaters and import in the workspace
        getattr(self.workspace4fit_,"import")(model_pdf_WJets_signalregion_after_correct_limit_variable)

        ##### calculate the normalization and alpha for limit datacard
        self.workspace4fit_.var("rrv_number%s_signalregion_%s_limit_variable"%(label,self.categoryLabel)).Print();
        self.workspace4fit_.var("rrv_number%s_in_obs0_variable_signalregion_from_fitting_%s"%(label,self.categoryLabel)).Print();
        self.workspace4fit_.var("rrv_number%s_signalregion_%s_limit_variable"%(label,self.categoryLabel)).setVal(self.workspace4fit_.var("rrv_number%s_in_obs0_variable_signalregion_from_fitting_%s"%(label,self.categoryLabel)).getVal());
        self.workspace4fit_.var("rrv_number%s_signalregion_%s_limit_variable"%(label,self.categoryLabel)).setError(self.workspace4fit_.var("rrv_number%s_in_obs0_variable_signalregion_from_fitting_%s"%(label,self.categoryLabel)).getError());

        self.workspace4fit_.var("rrv_number%s_signalregion_%s_limit_variable"%(label,self.categoryLabel)).setConstant(kTRUE);


    ##### Analysis with sideband alpha correction 
    def analysis_sideband_correction_method1(self):
        print "##################### Start sideband correction full analysis ##############";
        ### Fit all MC components in both obs0_variable and limit_variable
        self.fit_AllSamples_Mj_and_Mlvj();
        ### take the real data
        self.get_data()
        ### fit the WJets Normalization into the signal region -> no jet mass fluctuation has been done
        self.fit_WJetsNorm();
        ### fit data in the limit_variable low sideband with two different models
        self.fit_limit_variable_in_Mj_sideband("_WJets01","_lowersideband",self.MODEL_4_limit_variable_alter,1)
        self.fit_limit_variable_in_Mj_sideband("_WJets0","_lowersideband",self.MODEL_4_limit_variable,1)

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
        self.fit_limit_variable_in_Mj_sideband("_WJets0","_lowersideband", self.MODEL_4_limit_variable,1)
        #### prepare limit 
        self.prepare_limit("sideband_correction_method1",1,0,0)
        #### read the workspace
        self.read_workspace(1)
        '''


