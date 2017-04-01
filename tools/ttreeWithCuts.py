
# Helpful for making a small ntuple to test svFit changes on
# Can use this to change directory structure as well
# if desired

import ROOT


def ttreeWithCuts( oldFile, oldTreePath, fOutName='ttreeWithCut.root', cut='' ) :
    f = ROOT.TFile( oldFile, 'r' )
    t = f.Get(oldTreePath)
    print "Num Events in initial TTree:",t.GetEntries()
    
    fOut = ROOT.TFile(fOutName,'RECREATE')
    tOut = t.CopyTree( cut )
    print "Num events in new TTree:",tOut.GetEntries()

    fOut.cd()
    
    # Make same directory path
    info = oldTreePath.split('/')
    info.pop() # Get rid of TTree name
    newPath = '/'.join( info )
    fOut.mkdir( newPath )
    fOut.cd( newPath )
    print "New path: %s/%s" % (newPath, tOut.GetName() )

    tOut.Write()
    fOut.Close()



if __name__ == '__main__' :
    oldFile = 'Sync1Feb12/Sync-VBF125_9_tt.root'
    oldTreePath = 'tt/final/Ntuple'
    fOutName = 'vbfCut.root'
    cut = '(t1Pt>200 && t2Pt > 200)'
    ttreeWithCuts( oldFile, oldTreePath, fOutName, cut )
