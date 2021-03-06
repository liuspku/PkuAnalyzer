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

from PKUAnalyzerClass import doFit


############################################
#              Job steering                #
############################################

parser = OptionParser()

parser.add_option('-a', '--additioninformation',action="store",type="string",dest="additioninformation",default="EXO")
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-c', '--categoryID',action="store",type="int",dest="categoryID",default=1)

parser.add_option('-p', '--psmodel',action="store",type="string",dest="psmodel",default="pythia")

parser.add_option('-s','--simple', action='store', dest='simple', default=True, help='pre-limit in simple mode')
parser.add_option('-m','--multi', action='store_true', dest='multi', default=False, help='pre-limit in multi mode')

parser.add_option('--check', action='store_true', dest='check', default=False, help='check the workspace for limit setting')

parser.add_option('--cprime', action="store",type="int",dest="cprime",default=10)
parser.add_option('--BRnew', action="store",type="int",dest="BRnew",default=0)

parser.add_option('--closuretest', action='store',type="int", dest='closuretest', default=0, help='closure test; 0: no test; 1: A1->A2; 2: A->B')
parser.add_option('--fitSignal', action='store',type="int", dest='fitsignal', default=0, help='fit only signal lineshape with a chosen model')

parser.add_option('--inPath', action="store",type="string",dest="inPath",default="./")

parser.add_option('--category', action="store",type="string",dest="category",default="HP")


(options, args) = parser.parse_args()


### run analysis without systematics
def pre_limit_sb_correction_without_systermatic( categoryID, prime_signal_sample, in_limit_variable_signalregion_min, in_limit_variable_signalregion_max, in_obs0_variable_min, in_obs0_variable_max, in_limit_variable_min, in_limit_variable_max, fit_model, fit_model_alter):

    print "#################### pre_limit_sb_correction_without_systermatic: categoryID %s, signal %s, max and min signal region %f-%f, max and min mJ %f-%f, max and min limit_variable %f-%f, fit model %s and alternate %s ######################"%(categoryID,prime_signal_sample,in_limit_variable_signalregion_min,in_limit_variable_signalregion_max,in_obs0_variable_min,in_obs0_variable_max,in_limit_variable_min,in_limit_variable_max,fit_model,fit_model_alter);


    # data, MC of signal and backgrounds
    #file_Directory="AnaSigTree_WH/";
    file_Directory="AnaSigTree_newCategory_Pt300/";
    #sig_bkg_files:
    #0 line is nsig, 1st line is nbkg,
    #2ed line is for data,3 to 3+nsig-1 lines for sig channel;
    #3+nsig to 3+nsig+nbkg-1 for bkg channels.
    #First signal as prime signal. 
    #Signal and Background Line: Name, color number, file, fit_config
    #fit_config = ( decorrelation_or_not, fit shape name, parameter0 initial value, parameter1 initial value, parameter2 initial value, ..).;
    # when shape name= "Keys", it means not fit 
    sig_bkg_files=(
            1, # nsig: number of signal channels
            4, # nbkg: number of background channels
            ("data"   , 1  , file_Directory+"treeEDBR_data_xwh.root"                         , ( 0, "Hist") ), #data file
            ("Wp1000"  , 1  , file_Directory+"treeEDBR_MWp_1000_xwh.root", ( 0, "Hist") ), #sig
            ("WJets"  , 2  , file_Directory+"treeEDBR_WJetsPt180_xwh.root"                   , ( 0, "Hist") ), #bkg
            ("VV"     , 4  , file_Directory+"treeEDBR_VV_xwh.root"                           , ( 0, "Hist") ), #bkg
            ("TTbar"  , 210, file_Directory+"treeEDBR_TTBARpowheg_xwh.root"                  , ( 0, "Hist") ), #bkg
            ("SingleT", 7  , file_Directory+"treeEDBR_SingleTop_xwh.root"                    , ( 0, "Hist") ) #bkg
            );
    category_ID_label={
            0: ("elSRL0", "El SR L0"),
            1: ("elSRL2", "El SR L2"),
            2: ("elSRL0", "El SR L0"),
            3: ("elSBL2", "El SB L2"),
            4: ("muSRL0", "Mu SR L0"),
            5: ("muSRL2", "Mu SR L2"),
            6: ("muSRL0", "Mu SR L0"),
            7: ("muSBL2", "Mu SB L2"),
            };
    
    #analyzer_config
    analyzer_config={ 
            "sig_bkg_files": sig_bkg_files,
            "categoryID": categoryID,# mu/el, HP/LP 
            "categoryLabel": category_ID_label[categoryID],# mu/el, HP/LP 
            "limit_variable": ("mZZ","M_{WW}"),
            "limit_variable_full_range_min": in_limit_variable_min,
            "limit_variable_full_range_max": in_limit_variable_max,
            "limit_variable_BinWidth": 50.,
            "limit_variable_signalregion_range_min": in_limit_variable_signalregion_min,
            "limit_variable_signalregion_range_max": in_limit_variable_signalregion_max,
            "limit_variable_fit_model": fit_model,
            "limit_variable_fit_model": fit_model_alter,
            "obs0_variable": ("mJJNoKinFit","Jet Mass"),
            "obs0_variable_full_range_min": in_obs0_variable_min,
            "obs0_variable_full_range_max": in_obs0_variable_max,
            "obs0_variable_BinWidth": 5.,
            "obs0_variable_signalregion_range_min": 110,
            "obs0_variable_signalregion_range_max": 135,
            "obs0_variable_lowersideband_range_min": in_obs0_variable_min,
            "obs0_variable_lowersideband_range_max": 65,
            "obs0_variable_uppersideband_range_min": 135,
            "obs0_variable_uppersideband_range_max": in_obs0_variable_max,
            "fit_model": "ErfExp_v1", 
            "fit_model_alter": "ErfPow_v1",
            "additioninformation": options.additioninformation,
            "signal_scale":1e-2 
            };
    boostedW_fitter=doFit( analyzer_config);
    #boostedW_fitter.analysis("");
    #boostedW_fitter.analysis("DataDriven4WJetsNorm_fit_jetmass");
    boostedW_fitter.analysis("SimpleScale");

'''
### run full analysis
def pre_limit_sb_correction(method, categoryID, prime_signal_sample="BulkG_c0p2_M1000", in_limit_variable_signalregion_min=500, in_limit_variable_signalregion_max=700,
        in_obs0_variable_min=30, in_obs0_variable_max=140, in_limit_variable_min=400, in_limit_variable_max=1400, fit_model="ErfExp_v1", fit_model_alter="ErfPow_v1"): 

    print "#################### pre_limit_sb_correction: categoryID %s, signal %s, max and min signal region %f-%f, max and min mJ %f-%f, max and min limit_variable %f-%f, fit model %s and alternate %s ######################"%(categoryID,prime_signal_sample,in_limit_variable_signalregion_min,in_limit_variable_signalregion_max,in_obs0_variable_min,in_obs0_variable_max,in_limit_variable_min,in_limit_variable_max,fit_model,fit_model_alter);

    boostedW_fitter=doFit(categoryID, prime_signal_sample, in_limit_variable_signalregion_min,in_limit_variable_signalregion_max,in_obs0_variable_min,in_obs0_variable_max,in_limit_variable_min,in_limit_variable_max,fit_model,fit_model_alter);
    getattr(boostedW_fitter,"analysis_sideband_correction_%s"%(method) )();

### funtion to run just signal lineshape fit
def pre_fitsignal_only(categoryID, prime_signal_sample="BulkG_c0p2_M1000", in_limit_variable_signalregion_min=500, in_limit_variable_signalregion_max=700,
        in_obs0_variable_min=30, in_obs0_variable_max=140, in_limit_variable_min=400, in_limit_variable_max=1400, fit_model_narrow="CB_v1", fit_model_width="BWCB"): 

    print "#################### pre_fitsignal_only: categoryID %s, signal %s, max and min signal region %f-%f, max and min mJ %f-%f, max and min limit_variable %f-%f, fit model narrow %s, fit model width %s  ######################"%(categoryID,prime_signal_sample,in_limit_variable_signalregion_min,in_limit_variable_signalregion_max,in_obs0_variable_min,in_obs0_variable_max,in_limit_variable_min,in_limit_variable_max,fit_model_narrow, fit_model_width);

    boostedW_fitter=doFit(categoryID, prime_signal_sample, in_limit_variable_signalregion_min,in_limit_variable_signalregion_max,in_obs0_variable_min,in_obs0_variable_max,in_limit_variable_min,in_limit_variable_max);
    boostedW_fitter.fit_Signal(fit_model_narrow,fit_model_width);


### function to check the workspace once it has already created
def check_workspace(categoryID, higgs):
    boostedW_fitter = doFit(categoryID,higgs);
    boostedW_fitter.read_workspace()
'''
#### Main Code
if __name__ == '__main__':

    categoryID=options.categoryID;

    if options.simple and ( not options.multi) and ( not options.check) and ( not options.fitsignal):
        print '################# simple mode for %s sample'%(categoryID)
        pre_limit_sb_correction_without_systermatic(categoryID,"Wp_M1000",900,1100,40,150, 200,3000,"ExpN","ExpTail")

'''    if options.check:
        print '################# check workspace for %s sample'%(categoryID);
        check_workspace(categoryID,"BulkG_c0p2_M2000");

    ### real function called by the command line parsing some arguments as argv
    if options.multi  and ( not options.fitsignal):
        print '################# multi mode for %s sample'%(categoryID)
        pre_limit_sb_correction("method1",sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]), sys.argv[9], sys.argv[10] )

    ### only fit signal lineshape
    if options.fitsignal :
        pre_fitsignal_only(sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]),sys.argv[9],sys.argv[10])
'''
