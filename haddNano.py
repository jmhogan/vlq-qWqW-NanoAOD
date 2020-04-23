import os,sys,datetime,time, subprocess, math
from ROOT import *
execfile("/uscms_data/d3/jmanagan/EOSSafeUtils.py")

start_time = time.time()

## make inDir match your step1 directory, make outDir add "hadds" to the end
#inDir='/store/user/ejagt/NanoAODv6_1lep2018_022720_step1/'
inDir='/store/user/jmanagan/NanoAODv6_1lep2018_040920_step1/'
outDir='/store/user/jmanagan/NanoAODv6_1lep2018_040920_step1haddsHTSF/'
scratchDir='/uscmst1b_scratch/lpc1/3DayLifetime/jmanagan/'

if not os.path.exists(scratchDir): os.system('mkdir '+scratchDir)
os.system('eos root://cmseos.fnal.gov/ mkdir -p '+outDir)

dirList = [
# 'NanoAOD_25Oct_SingleEl',
# 'NanoAOD_25Oct_SingleMu',

# 'QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8',
'ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8',
# 'ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8',
# 'ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8',
# 'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8',
# 'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8',


# 'DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',
# 'DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8',
# 'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
# 'TTToHadronic_TuneCP5_13TeV-powheg-pythia8',
# 'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8',
# 'TprimeTprime_M-1000_TuneCP5_PSweights_13TeV-madgraph-pythia8',
# 'TprimeTprime_M-1100_TuneCP5_PSweights_13TeV-madgraph-pythia8',
# 'TprimeTprime_M-1200_TuneCP5_PSweights_13TeV-madgraph-pythia8',
# 'TprimeTprime_M-1300_TuneCP5_PSweights_13TeV-madgraph-pythia8',
# 'TprimeTprime_M-1400_TuneCP5_PSweights_13TeV-madgraph-pythia8',
# 'TprimeTprime_M-1500_TuneCP5_PSweights_13TeV-madgraph-pythia8',
# 'TprimeTprime_M-1600_TuneCP5_PSweights_13TeV-madgraph-pythia8',
# 'TprimeTprime_M-1700_TuneCP5_PSweights_13TeV-madgraph-pythia8',
# 'TprimeTprime_M-1800_TuneCP5_PSweights_13TeV-madgraph-pythia8',
# 'TprimeTprime_M-900_TuneCP5_PSweights_13TeV-madgraph-pythia8',
# 'WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8',
# 'WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8',
]

for sample in dirList:
    outList = ['none']
    #if 'Tprime' in sample: outList = ['BWBW','TZBW','THBW','TZTH','TZTZ','THTH']
    #elif 'Bprime' in sample: outList = ['TWTW','BZTW','BHTW','BZBH','BZBZ','BHBH']
    #elif 'TTTo' in sample: outList = ['Mtt0to700','Mtt700to1000','Mtt1000toInf']

    for outlabel in outList:

        outsample = sample+'_'+outlabel
        if outlabel == 'none': outsample = sample

        rootfiles = EOSlist_root_files(inDir+'/'+outsample)

        print "------------ hadding Sample:",outsample,"---------------"
        print 'N root files in',outsample,'=',len(rootfiles)


        nFilesPerHadd = 900

        onefile = ' root://cmseos.fnal.gov/'+inDir+'/'+outsample+'/'+rootfiles[-1]
        manyfiles = nFilesPerHadd*onefile
        lengthcheck = len('hadd -f -k root://cmseos.fnal.gov/'+outDir+'/'+outsample+'_hadd.root '+manyfiles)
        if lengthcheck > 126000.:
            toolong = lengthcheck - 126000.
            num2remove = math.ceil(toolong/len(onefile))
            nFilesPerHadd = int(nFilesPerHadd - num2remove)
            print 'Length estimate reduced from',lengthcheck,'by',toolong,'via removing',num2remove,'files for nFilesPerHadd of',nFilesPerHadd

        if len(rootfiles) < nFilesPerHadd:
            haddcommand = 'hadd -f -k '+scratchDir+'/'+outsample+'_hadd.root '
            for file in rootfiles:
                haddcommand+=' root://cmseos.fnal.gov/'+inDir+'/'+outsample+'/'+file
            print 'Length of hadd command =',len(haddcommand)
            subprocess.call(haddcommand,shell=True)
            
            xrdcpcommand = 'xrdcp -f '+scratchDir+'/'+outsample+'_hadd.root root://cmseos.fnal.gov/'+outDir+'/'+outsample+'_hadd.root'
            subprocess.call(xrdcpcommand,shell=True)

            #if bool(EOSisfile(outDir+'/'+outsample+'_hadd.root')) != True:
            #    print haddcommand                
        else:
            for i in range(int(math.ceil(len(rootfiles)/float(nFilesPerHadd)))):
                haddcommand = 'hadd -f '+scratchDir+'/'+outsample+'_'+str(i+1)+'_hadd.root '

                begin=i*nFilesPerHadd
                end=begin+nFilesPerHadd
                if end > len(rootfiles): end=len(rootfiles)
                print 'begin:',begin,'end:',end

                for j in range(begin,end):
                    haddcommand+=' root://cmseos.fnal.gov/'+inDir+'/'+outsample+'/'+rootfiles[j]
                print 'Length of hadd command =',len(haddcommand)
                subprocess.call(haddcommand,shell=True)

                xrdcpcommand = 'xrdcp -f '+scratchDir+'/'+outsample+'_'+str(i+1)+'_hadd.root root://cmseos.fnal.gov/'+outDir+'/'+outsample+'_'+str(i+1)+'_hadd.root'
                subprocess.call(xrdcpcommand,shell=True)

                #if bool(EOSisfile(outDir+'/'+outsample+'_'+str(i+1)+'_hadd.root')) != True:
                #    print haddcommand

print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))



