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

    # For 2015 ZTT work
    samples = ['dataTT-C', 'dataTT-D']
    originalDir = '/nfs_scratch/truggles/httSept29_2015TauDataFF'
    targetDir = '/nfs_scratch/truggles/httSept29_2015TauDataFFMerged'
    jobId = ''
    channel = 'tt'
    ttreePath = 'tt/final/Ntuple'
    for sample in samples :
        mergeSample( jobId, sample, channel, ttreePath, originalDir, targetDir )



