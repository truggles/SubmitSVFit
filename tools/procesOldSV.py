import os, sys, re

#root -l "removeOldSV.C(\"mutau_tree\", \"/scratch/kkaadze/SkimmedNtuples_SMHTT2016/smhmt_20march/DY.root\", \"DY.root\")"
dir = sys.argv[1]
location = os.popen("ls " + dir + "/*root")
files = location.readlines()

for file in files:
    fileName = re.split("/", file)[-1][:-1]
    print fileName
    command = "root -l \"removeOldSV.C(\\\"mutau_tree\\\", \\\"" + file[:-1] + "\\\", \\\"" + fileName + "\\\")\""
    os.system(command)
