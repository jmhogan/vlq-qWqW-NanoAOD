import os,sys

filelist = {
## put all your hadd files here by name. 
##If there are any _1.root _2.root, etc, put :Nfiles, otherwise :0 if there's only 1

 'DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8_hadd.root':0,
 'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_hadd.root':0,
 'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_1_hadd.root':2,
 #'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_hadd.root':0,
 'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_hadd.root':0,
 'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_hadd.root':0,
 'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_hadd.root':0,
 'TTToHadronic_TuneCP5_13TeV-powheg-pythia8_hadd.root':0,
 'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_1_hadd.root':3,
 #'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_hadd.root':0,
 'TprimeTprime_M-1000_TuneCP5_PSweights_13TeV-madgraph-pythia8_hadd.root':0,
 'TprimeTprime_M-1100_TuneCP5_PSweights_13TeV-madgraph-pythia8_hadd.root':0,
 'TprimeTprime_M-1200_TuneCP5_PSweights_13TeV-madgraph-pythia8_hadd.root':0,
 'TprimeTprime_M-1300_TuneCP5_PSweights_13TeV-madgraph-pythia8_hadd.root':0,
 'TprimeTprime_M-1400_TuneCP5_PSweights_13TeV-madgraph-pythia8_hadd.root':0,
 'TprimeTprime_M-1500_TuneCP5_PSweights_13TeV-madgraph-pythia8_hadd.root':0,
 'TprimeTprime_M-1600_TuneCP5_PSweights_13TeV-madgraph-pythia8_hadd.root':0,
 'TprimeTprime_M-1700_TuneCP5_PSweights_13TeV-madgraph-pythia8_hadd.root':0,
 'TprimeTprime_M-1800_TuneCP5_PSweights_13TeV-madgraph-pythia8_hadd.root':0,
 'TprimeTprime_M-900_TuneCP5_PSweights_13TeV-madgraph-pythia8_hadd.root':0,
 'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8_hadd.root':0,
 'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8_hadd.root':0,
}
 
from ROOT import TFile, TH1
 
for filekey in sorted(filelist.keys()):
    print('-------------------------------------------------------')
    ## Put your file path here and below
    rfile = TFile.Open('root://cmseos.fnal.gov//store/user/jmanagan/NanoAODv6_1lep2018_040920_step1haddsHTSF/'+filekey)
    hist = rfile.Get("nEventsWeighted")

    if filelist[filekey] > 0:
        print 'Opening ',filelist[filekey],'files:'
        for ifile in range(2,filelist[filekey]+1):
            #print 'file #',ifile
            tempfile = TFile.Open('root://cmseos.fnal.gov//store/user/jmanagan/NanoAODv6_1lep2018_040920_step1haddsHTSF/'+filekey.replace('_1_','_'+str(ifile)+'_'))
            temphist = tempfile.Get("nEventsWeighted")
            temphist.SetDirectory(0)
            hist.Add(temphist)
            tempfile.Close()                                  

    adjusted = hist.GetBinContent(3) - hist.GetBinContent(1)
    integral = hist.GetBinContent(3) + hist.GetBinContent(1)
    #newpdf = hist.GetBinContent(2)

    #if 'prime' not in filekey:
    print(str(adjusted)+'. # from integral '+str(integral)+', file '+filekey)
    #else:
    #    print(str(newpdf)+'. # from integral '+str(integral)+', file '+filekey)







