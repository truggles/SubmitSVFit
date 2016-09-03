# SubmitSVFit
```
cmsrel CMSSW_7_6_3 #for 2016 analysis
cd CMSSW_7_6_3/src/
cmsenv
git cms-init 
git clone --recursive https://github.com/truggles/SubmitSVFit.git
cd SubmitSVFit
source recipe.sh
scram b -j 8
```

There are a number of current versions of the stand alone svFitter included here.
You can find them in ROOT/bin/SVFitStandAlone... with the names of their executables
defined ROOT/bin/BuildFile.xml.

Most current on being worked on: 
ROOT/bin/SVFitStandAloneFSA.cc

To run in interactive mode for example:
```
SVFitStandAloneFSA inputFile=coolInputFile.root newOutputFile=1 newFile=tmpOut.root doES=1
```

 - inputFile = obvious
 - newOutputFile = 0/1
   - 0 = update input file with svFit vars
   - 1 = output new file with original TTree and new svFit vars
 - newFile = name of output file, default is newFile.root if none specified
 - doES = apply energy scale adjustments providing nominal, shift UP and shift DOWN values of svFit
   - 0 = default, no shift
   - 1 = apply shifts

To submit jobs to condor:
```
cd test
python svFitSubmitter.py -dr -sd /hdfs/store/user/truggles/svFitTestSept03 -es=1 --jobName svFitSept03Test
```

 - -dr = dryRun and outputs a command for you to run
 - -sd = select directory, the input directory, this will produce a list of all files in that directory to run over<BR>
       you must have your files in a /hdfs directory.
 - -es = apply energy scale, see above
 - --jobName = applys a job name to append to your file names and output directory structure


To get your files from elsewhere to /hdfs do something like this:
```
gsido mkdir /hdfs/store/user/truggles/mySubmitDir
gsido cp /path/to/submit/files.root /hdfs/store/user/truggles/mySubmitDir
```

It is VERY helpful to make sure that you have ~1000 events per file when running this on Condor.  Anything much larger will take forever,
especially if you run with Energy Shifts. If using FSA ntuples:
 - skim your ntuples without merging any files
 - do a controlled merge that hadds ~1000 events / output file

To do this controlled merge edit the file tools/controlledMerge.py and specify your
 - original directory
 - samples
 - TTreePath
 - output directory
 - and edit the event count / file if you would like to adjust it away from 1,000
Depending on your file naming convention, you may have to edit line 16<BR>

 
Then:
```
python tools/controlledMerge.py
```


