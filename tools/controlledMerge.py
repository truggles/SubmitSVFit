# From truggles
# https://github.com/truggles/Z_to_TauTau_13TeV/blob/master/util/svFitMerger.py

import ROOT
import os, glob, subprocess
import multiprocessing


# Check if directory exists, make it if not
def checkDir( dirName ) :
    if not os.path.exists( dirName ) : os.makedirs( dirName )




def mergeSample( sample, channel, ttreePath, originalDir, targetDir ) :
    files = glob.glob(originalDir+'/%s_*_%s.root' % (sample, channel) )
    checkDir( targetDir )
    for file in files :
        print file

    # Add return if no files found
    # this is useful for all the many signal masses
    if len(files) == 0 : 
        print "\n\n"
        print "No files of samples %s found in dir %s" % (sample, originalDir)
        print "\n\n"
        return

    rep = 0
    runningSize = 0
    runningNumFiles = 0
    toMerge = []
    ints = []
    for file_ in files :

        # Merge to ~ 1000 events per file
        f = ROOT.TFile(file_,'r')
        t = f.Get(ttreePath)
        size = t.GetEntries()
        print size,"   ",file_
        runningSize += size
        runningNumFiles += 1
        print "running size: ",runningSize
        toMerge.append( file_ )
        if runningSize > 1000 or runningNumFiles == 500 :
        #if runningSize > 300 or runningNumFiles == 500 :
            runningSize = 0
            runningNumFiles = 0
            mergeList = ["hadd", "-f", targetDir+"/%s_%i_%s.root" % (sample, rep, channel)]
            for f in toMerge :
                mergeList.append( f )
            subprocess.call( mergeList )
            ints = []
            toMerge = []
            rep += 1
    mergeList = ["hadd", "-f", targetDir+"/%s_%i_%s.root" % (sample, rep, channel)]
    for f in toMerge :
        mergeList.append( f )
    if len( mergeList ) > 3 : # greater than 3 means we actually have a file to merge (not empty)
        subprocess.call( mergeList )



if __name__ == '__main__' :
    ''' Start multiprocessing tests '''
    pool = multiprocessing.Pool(processes = 6 )
    multiprocessingOutputs = []
    debug = False
    doAZH = True
    doHTT = False

    ''' SM-HTT Feb 20, 2017 '''
    if doHTT :
    
        dir1 = 'httMay30Rivet'
        originalDir = '/nfs_scratch/truggles/'+dir1
        channel = 'tt'
        ttreePath = 'tt/final/Ntuple'
    
        """ section 1, Need TES, Recoil type 2, no WJets """
        samples = ['DYJets', 'DYJets1', 'DYJets2', 'DYJets3', 'DYJets4', 'DYJetsLow', 'DYJets1Low', 'DYJets2Low',] 
        for mass in [110, 120, 125, 130, 140] :
            samples.append('ggHtoTauTau%i' % mass)
            samples.append('VBFHtoTauTau%i' % mass)
            samples.append('VBFHtoWW2l2nu%i' % mass)
            samples.append('HtoWW2l2nu%i' % mass)
    
        targetDir = '/nfs_scratch/truggles/'+dir1+'_svFit_prep/Recoil2_TES1_WJ0'
        for sample in samples :
            mergeSample( sample, channel, ttreePath, originalDir, targetDir )
    
    
    
    #    """ section 2, Need TES, Recoil type 2, WJets """
    #    samples = ['WJets', 'WJets1', 'WJets2', 'WJets3', 'WJets4']
    #
    #    targetDir = '/nfs_scratch/truggles/'+dir1+'_svFit_prep/Recoil2_TES1_WJ1'
    #    for sample in samples :
    #        mergeSample( sample, channel, ttreePath, originalDir, targetDir )
    
    
    
        """ section 3, Need TES, no recoil, no WJets """
        samples = ['EWKZ2l', 'EWKZ2nu', 'EWKWMinus', 'EWKWPlus', 'T-tchan', 'Tbar-tchan', 'TT', 'Tbar-tW', 'T-tW', 'VV', 'WW1l1nu2q', 'WWW', 'WZ1l1nu2q', 'WZ1l3nu', 'WZ2l2q', 'WZ3l1nu', 'ZZ2l2q', 'ZZ4l'] # Feb17 for Moriond17 
        for mass in [110, 120, 125, 130, 140] :
            samples.append('WMinusHTauTau%i' % mass)
            samples.append('WPlusHTauTau%i' % mass)
            samples.append('ZHTauTau%i' % mass)
            samples.append('ttHTauTau%i' % mass)
    
        targetDir = '/nfs_scratch/truggles/'+dir1+'_svFit_prep/Recoil0_TES1_WJ0'
        for sample in samples :
            mergeSample( sample, channel, ttreePath, originalDir, targetDir )
    
    
    
        
    #    """ section 4, DATA No TES, no recoil, no WJets """
    #    samples = []
    #    #for era in ['B', 'C', 'D', 'E', 'F', 'G', 'H'] :
    #    for era in ['H',] :
    #        samples.append('dataTT-%s' % era)
    #
    #    targetDir = '/nfs_scratch/truggles/'+dir1+'_svFit_prep/Recoil0_TES0_WJ0'
    #    for sample in samples :
    #        mergeSample( sample, channel, ttreePath, originalDir, targetDir )







    # For svFit, skip eeee, mmmm channels
    if doAZH :
        # AZH June 01 Wisconsin -> uwlogin
        azhSamples = ['DYJets', 'DYJets1', 'DYJets2', 'DYJets3', 'DYJets4', 'ggZZ4m', 'ggZZ2e2m', 'ggZZ2e2tau', 'ggZZ4e', 'ggZZ2m2tau', 'ggZZ4tau', 'TT', 'WWW', 'WWZ', 'WZ3l1nu', 'WZZ', 'WZ', 'ZZ4l', 'ZZZ',] # May 31 samples, no ZZ->all, use ZZ4l

        for mass in [110, 120, 125, 130, 140] :
            #azhSamples.append('ggHtoTauTau%i' % mass)
            #azhSamples.append('VBFHtoTauTau%i' % mass)
            #azhSamples.append('WMinusHTauTau%i' % mass)
            #azhSamples.append('WPlusHTauTau%i' % mass)
            azhSamples.append('ZHTauTau%i' % mass)
            #azhSamples.append('ttHTauTau%i' % mass)
        for mass in [125,] :
            azhSamples.append('WMinusHTauTau%i' % mass)
            azhSamples.append('WPlusHTauTau%i' % mass)
            azhSamples.append('ZHWW%i' % mass)

        azhSamples = ['ggZZ4m','ggZZ4e']
        for mass in [220, 240, 260, 280, 300, 320, 340, 350, 400] :
            azhSamples.append('azh%i' % mass)

        #for era in ['B', 'C', 'D', 'E', 'F', 'G', 'H'] :
        #    azhSamples.append('dataEE-%s' % era)
        #    azhSamples.append('dataMM-%s' % era)

        #azhSamples = ['azh340',]

        #originalDir = '/nfs_scratch/truggles/azhMay31skim/'
        originalDir = '/nfs_scratch/truggles/azhJune06skim/'
        targetDir = '/nfs_scratch/truggles/azhMay31AndJune06skim_svFitPrep'
        jobId = ''
        channels = ['eeet','eett','eemt','eeem','emmt','mmtt','mmmt','emmm',] # 8
        #channels = ['emmm',] # 8
        #azhSamples = ['ZZ4l',]
        for channel in channels :
            ttreePath = channel+'/final/Ntuple'
            for sample in azhSamples :
                if debug:
                    mergeSample( sample, channel, ttreePath, originalDir, targetDir )
                else :
                    multiprocessingOutputs.append( pool.apply_async(mergeSample, args=(
                        sample,
                        channel,
                        ttreePath,
                        originalDir,
                        targetDir )))
        if not debug :
            mpResults = [p.get() for p in multiprocessingOutputs]




