import os,shutil,datetime,time
import getpass, subprocess
from ROOT import *
import time
execfile("/uscms_data/d3/jmanagan/EOSSafeUtils.py")

start_time = time.time()

#IO directories must be full paths

## USER-specific file paths
relbase = '/uscms_data/d3/jmanagan/ErinHitFit/CMSSW_10_2_10/'
tarfile = '/uscms_data/d3/jmanagan/qWqW_nano.tar'
outDir='/store/user/jmanagan/NanoAODv6_1lep2018_051620_step1/'
condorDir='/uscms_data/d3/jmanagan/NanoAODv6_1lep2018_051620_step1/'

runDir=os.getcwd()

listDir = runDir+'/fileListsNano'
if not os.path.exists(listDir): os.system('mkdir '+listDir)

print 'Making tar:'
if os.path.exists(tarfile): print '*********** tar already exists! I ASSUME YOU WANT TO MAKE A NEW ONE! *************'

os.chdir(relbase)
# YOU NEED TO EXCLUDE ANYTHING ELSE THAT MIGHT LIVE IN THE SAME CMSSW RELEASE, MY LIST IS SUPER LONG
print 'tar --exclude="tmp/" --exclude="src/PhysicsTools" --exclude="src/analysisVLQ2019" --exclude="src/bwbw_2018" --exclude="src/vlq-qWqW-NanoAOD/*.log" --exclude="src/bwbw_2018" --exclude="src/.git" -zcf '+tarfile+' ./*'

os.system('tar --exclude="tmp/" --exclude="src/PhysicsTools" --exclude="src/analysisVLQ2019" --exclude="src/bwbw_2018" --exclude="src/vlq-qWqW-NanoAOD/*.log" --exclude="src/bwbw_2018" --exclude="src/.git" -zcf '+tarfile+' ./*')
os.chdir(runDir)

print 'Starting submission'
count=0

dirList = [
    '/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/NANOAODSIM',
    '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/NANOAODSIM',
    '/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/NANOAODSIM',    

    'NanoAOD_25Oct_SingleEl',
    'NanoAOD_25Oct_SingleMu',
    '/TprimeTprime_M-1000_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/TprimeTprime_M-1100_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/TprimeTprime_M-1200_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/TprimeTprime_M-1300_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/TprimeTprime_M-1400_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/TprimeTprime_M-1500_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/TprimeTprime_M-1600_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/TprimeTprime_M-1700_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/TprimeTprime_M-1800_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/TprimeTprime_M-900_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',

    '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v3/NANOAODSIM',
    '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',

    '/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',

    '/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20_ext1-v1/NANOAODSIM',
    '/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',
    '/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv6-Nano25Oct2019_102X_upgrade2018_realistic_v20-v1/NANOAODSIM',



]

for sample in dirList:
    print "------------ Sample:",sample,"---------------"
    outList = ['none']
    #if 'TTTo' in sample: outList = ['Mtt0to700','Mtt700to1000','Mtt1000toInf']
    
    isData = False
    if 'SingleMu' in sample or 'SingleEl' in sample or 'EGamma' in sample: isData = True
    fullName = sample
    if not isData: sample = sample.split('/')[1]
    
    rootfiles = []

    ## MC is being pulled from published datasets via dasgoclient
    if not isData:
        command = '/cvmfs/cms.cern.ch/common/dasgoclient --query="file dataset='+fullName+'" > '+listDir+'/'+sample+'.txt'
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        if 'TTToSemiLeptonic' in sample:
            command = '/cvmfs/cms.cern.ch/common/dasgoclient --query="file dataset='+fullName.replace("v20","v20_ext3")+'" >> '+listDir+'/'+sample+'.txt'
            proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
        elif 'TTToHadronic' in sample:
            command = '/cvmfs/cms.cern.ch/common/dasgoclient --query="file dataset='+fullName.replace("v20","v20_ext2")+'" >> '+listDir+'/'+sample+'.txt'
            proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)

        os.system('sleep 1')
        os.system('touch '+listDir+'/'+sample+'.txt')
        os.system('sleep 30')
        rootlist = open(listDir+'/'+sample+'.txt')

        for line in rootlist:
            rootfiles.append('root://cmsxrootd.fnal.gov/'+line.strip())  ## Can change to whatever is normal for running @ CERN
        rootlist.close()
        print 'Added ',len(rootfiles),'to rootfiles'

    ## Data has been run through a CRAB job for lumimask, so it's pulled from an EOS path
    else:
        inputDir = '/store/user/jmanagan/'
        finalStateYear = 'SingleMuon/NanoTestPost'
        if 'SingleEl' in sample: finalStateYear = 'EGamma/NanoTestPost'
        runlist = EOSlistdir(inputDir+'/'+sample+'/'+finalStateYear+'/')
        print "Running",len(runlist),"crab directories"

        runcounter = 0
        for run in runlist:
            runcounter += 1
            numlist = EOSlistdir(inputDir+'/'+sample+'/'+finalStateYear+'/'+run+'/')
            
            for num in numlist:
                numpath = inputDir+'/'+sample+'/'+finalStateYear+'/'+run+'/'+num
                pathsuffix = numpath.split('/')[-3:]
                pathsuffix = '/'.join(pathsuffix)

                rootlist = EOSlist_root_files(numpath)            
                for rootfile in rootlist:
                    rootfiles.append('root://cmseos.fnal.gov/'+numpath+'/'+rootfile)

    for outlabel in outList:
        tmpcount = 0
    
        outsample = sample+'_'+outlabel
        if outlabel == 'none': outsample = sample
    
        os.system('eos root://cmseos.fnal.gov/ mkdir -p '+outDir+outsample)
        os.system('mkdir -p '+condorDir+outsample)

        for rootfile in rootfiles:
            tmpcount += 1
            #if tmpcount > 1: continue  ### TEST JOB

            ## Need to limit number of events per job for LPC walltime
            n_jobs = 1
            maxEvtsPerJob = 200000
            if 'QCD' in sample: maxEvtsPerJob = 1200000
            if not isData:
                command = '/cvmfs/cms.cern.ch/common/dasgoclient --query="file='+rootfile.replace('root://cmsxrootd.fnal.gov/','')+' | grep file.nevents" '
                proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
                (out, err) = proc.communicate()
                print out
                try: nevents = int(out.split('\n')[0])
                except:
                    try: nevents = int(out.split('\n')[1])
                    except: 
                        print 'ERROR: couldnt isolate the number of events'
                        exit()

                n_jobs = int(nevents) / int(maxEvtsPerJob)
                if int(nevents) % int(maxEvtsPerJob) > 0: n_jobs += 1 ## and extra one to account for the remainder

            print outsample+'_'+str(tmpcount),'=',n_jobs,'jobs'

            ### split based on the number of events
            for i_split in range(n_jobs):
            
                count+=1
                maxEvents = int(maxEvtsPerJob)
                skipEvents = int(maxEvtsPerJob*i_split)
                if i_split == n_jobs-1:
                    maxEvents = nevents - maxEvtsPerJob*(n_jobs-1) ## up to the last event
                
                dict={'RUNDIR':runDir, 'SAMPLE':outsample, 'FILENAME':rootfile, 'OUTPUTDIR':outDir+outsample, 'ID':str(tmpcount)+'_'+str(i_split), 
                      'TARBALL':tarfile, 'SKIP':skipEvents, 'MAX':maxEvents}
                jdfName=condorDir+'/%(SAMPLE)s/%(SAMPLE)s_%(ID)s.job'%dict
                print jdfName
                jdf=open(jdfName,'w')
                jdf.write(
                    """use_x509userproxy = true
universe = vanilla
Executable = %(RUNDIR)s/runCondorAnalysis.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = %(TARBALL)s
Output = %(SAMPLE)s_%(ID)s.out
Error = %(SAMPLE)s_%(ID)s.err
Log = %(SAMPLE)s_%(ID)s.log
Notification = Never
Arguments = "%(FILENAME)s %(SAMPLE)s_%(ID)s.root %(OUTPUTDIR)s %(SKIP)s %(MAX)s"

Queue 1"""%dict)
                jdf.close()
                os.chdir('%s/%s'%(condorDir,outsample))
                os.system('condor_submit %(SAMPLE)s_%(ID)s.job'%dict)
                os.system('sleep 0.5')                                
                os.chdir('%s'%(runDir))
                print count, "jobs submitted!!!"
        
print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))
