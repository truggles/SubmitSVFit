# From truggles
# https://github.com/truggles/Z_to_TauTau_13TeV/blob/master/util/svFitMerger.py

import ROOT
import os, glob, subprocess


# Check if directory exists, make it if not
def checkDir( dirName ) :
    if not os.path.exists( dirName ) : os.makedirs( dirName )




def mergeSample( sample, channel, ttreePath, originalDir, targetDir ) :
    files = glob.glob(originalDir+'/%s_*_%s.root' % (sample, channel) )
    checkDir( targetDir )

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

    ''' SM-HTT Feb 20, 2017 '''

    originalDir = '/nfs_scratch/truggles/httFeb28sigs'
    channel = 'tt'
    ttreePath = 'tt/final/Ntuple'

    """ section 1, Need TES, Recoil type 2, no WJets """
    samples = ['DYJets', 'DYJets1', 'DYJets2', 'DYJets3', 'DYJets4', 'DYJetsLow', 'DYJets1Low', 'DYJets2Low',] 
    samples = []
    for mass in [110, 120, 125, 130, 140] :
        samples.append('ggHtoTauTau%i' % mass)
        samples.append('VBFHtoTauTau%i' % mass)
        samples.append('VBFHtoWW2l2nu%i' % mass)
        samples.append('HtoWW2l2nu%i' % mass)

    targetDir = '/nfs_scratch/truggles/httFebFeb28sigs_svFit_prep/Recoil2_TES1_WJ0'
    for sample in samples :
        mergeSample( sample, channel, ttreePath, originalDir, targetDir )



#    """ section 2, Need TES, Recoil type 2, WJets """
#    samples = ['WJets', 'WJets1', 'WJets2', 'WJets3', 'WJets4']
#
#    targetDir = '/nfs_scratch/truggles/httFebFeb28sigs_svFit_prep/Recoil2_TES1_WJ1'
#    for sample in samples :
#        mergeSample( sample, channel, ttreePath, originalDir, targetDir )



    """ section 3, Need TES, no recoil, no WJets """
    samples = ['EWKWMinus', 'EWKWPlus', 'T-tchan', 'Tbar-tchan', 'TT', 'Tbar-tW', 'T-tW', 'VV', 'WW1l1nu2q', 'WWW', 'WZ1l1nu2q', 'WZ1l3nu', 'WZ2l2q', 'WZ3l1nu', 'ZZ2l2q', 'ZZ4l'] # Feb17 for Moriond17 
    samples = ['EWKWMinus', 'EWKWPlus', 'EWKZ2l', 'EWKZ2nu']
    for mass in [125,] :
        samples.append('WMinusHTauTau%i' % mass)
        samples.append('WPlusHTauTau%i' % mass)
        samples.append('ZHTauTau%i' % mass)
        samples.append('ttHTauTau%i' % mass)

    targetDir = '/nfs_scratch/truggles/httFebFeb28sigs_svFit_prep/Recoil0_TES1_WJ0'
    for sample in samples :
        mergeSample( sample, channel, ttreePath, originalDir, targetDir )



    
#    """ section 4, DATA No TES, no recoil, no WJets """
#    samples = []
#    for era in ['B', 'C', 'D', 'E', 'F', 'G', 'H'] :
#        samples.append('dataTT-%s' % era)
#
#    targetDir = '/nfs_scratch/truggles/httFebFeb28sigs_svFit_prep/Recoil0_TES0_WJ0'
#    for sample in samples :
#        mergeSample( sample, channel, ttreePath, originalDir, targetDir )










