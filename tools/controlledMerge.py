# From truggles
# https://github.com/truggles/Z_to_TauTau_13TeV/blob/master/util/svFitMerger.py

import ROOT
import os, glob, subprocess


# Check if directory exists, make it if not
def checkDir( dirName ) :
    if not os.path.exists( dirName ) : os.makedirs( dirName )




def mergeSample( jobId, sample, channel, ttreePath, originalDir, targetDir ) :
    files = glob.glob(originalDir+'/%s%s_*_%s.root' % (jobId, sample, channel) )
    checkDir( targetDir )

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

    # HTT Aug 23, hdfs -> UW
    #samples = ['DYJetsAMCNLO', 'DYJets', 'DYJets1', 'DYJets2', 'DYJets3', 'DYJets4', 'EWKWPlus', 'EWKWMinus', 'EWKZ2l', 'EWKZ2nu', 'ZZ2l2q', 'WWW', 'WWZ', 'WZZ', 'WZZ', 'ZZZ', 'T-tchan', 'Tbar-tchan', 'TT', 'Tbar-tW', 'T-tW', 'WJets', 'WJets1', 'WJets2', 'WJets3', 'WJets4', 'WW1l1nu2q', 'WZ1l1nu2q', 'WZ1l3nu', 'WZ2l2q', 'ZZ2l2q', 'VV', 'dataTT', 'VBFHtoTauTau120', 'VBFHtoTauTau125', 'VBFHtoTauTau130', 'ggHtoTauTau120', 'ggHtoTauTau125', 'ggHtoTauTau130']
    #samples = ['EWKWPlus', 'EWKWMinus', 'EWKZ2l', 'EWKZ2nu', 'ZZ2l2q', 'WWW', 'WWZ', 'WZZ', 'WZZ', 'ZZZ']
    #originalDir = '/nfs_scratch/truggles/httSept25ewk'
    #targetDir = '/nfs_scratch/truggles/httSept25ewkSkimMerged'
    #jobId = ''
    #channel = 'tt'
    #ttreePath = 'tt/final/Ntuple'
    #for sample in samples :
    #    mergeSample( jobId, sample, channel, ttreePath, originalDir, targetDir )




    ''' SM-HTT Oct 31 '''

    #samples = ['DYJets', 'DYJets1', 'DYJets2', 'DYJets3', 'DYJets4',] # B-F data, add WW, WZ3l1nu, WZ, ZZ4l, ZZ, removed WWZ, WZZ
    #for mass in [120, 125, 130] :
    #    samples.append('ggHtoTauTau%i' % mass)
    #    samples.append('VBFHtoTauTau%i' % mass)
    #    samples.append('WMinusHTauTau%i' % mass)
    #    samples.append('WPlusHTauTau%i' % mass)
    #    samples.append('ZHTauTau%i' % mass)

    originalDir = '/nfs_scratch/truggles/httOct31svFitPrep'
    targetDir = '/nfs_scratch/truggles/httOct31svFitPrepMerged/DYJetsAndHiggs'
    jobId = ''
    channel = 'tt'
    ttreePath = 'tt/final/Ntuple'
    for sample in samples :
        mergeSample( jobId, sample, channel, ttreePath, originalDir, targetDir )


    #samples = ['WJets', 'WJets1', 'WJets2', 'WJets3', 'WJets4'] # B-F data, add WW, WZ3l1nu, WZ, ZZ4l, ZZ, removed WWZ, WZZ
    #targetDir = '/nfs_scratch/truggles/httOct31svFitPrepMerged/WJets'
    #jobId = ''
    #channel = 'tt'
    #ttreePath = 'tt/final/Ntuple'
    #for sample in samples :
    #    mergeSample( jobId, sample, channel, ttreePath, originalDir, targetDir )

    #samples = ['EWKWPlus', 'EWKWMinus', 'EWKZ2l', 'EWKZ2nu', 'WWW', 'ZZZ', 'T-tchan', 'Tbar-tchan', 'TT', 'Tbar-tW', 'T-tW', 'WW1l1nu2q', 'WZ1l1nu2q', 'WZ1l3nu', 'WZ2l2q', 'ZZ2l2q', 'VV', 'dataTT-B', 'dataTT-C', 'dataTT-D', 'dataTT-E', 'dataTT-F', 'WW', 'WZ3l1nu', 'WZ', 'ZZ4l', 'ZZ'] # B-F data, add WW, WZ3l1nu, WZ, ZZ4l, ZZ, removed WWZ, WZZ
    #samples = ['WW', 'WZ3l1nu', 'WZ', 'ZZ4l', 'ZZ'] # B-F data, add WW, WZ3l1nu, WZ, ZZ4l, ZZ, removed WWZ, WZZ
    #targetDir = '/nfs_scratch/truggles/httOct31svFitPrepMerged/AllOthers'
    #jobId = ''
    #channel = 'tt'
    #ttreePath = 'tt/final/Ntuple'
    #for sample in samples :
    #    mergeSample( jobId, sample, channel, ttreePath, originalDir, targetDir )


