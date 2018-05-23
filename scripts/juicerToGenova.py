# Import bits and bobs
import csv
import subprocess
import pandas as pd
import math
import numpy as np
import argparse
import os
import sys

# Arguments list
parser = argparse.ArgumentParser(description='Convert juicer.hic files to GENOVA-input.')
parser.add_argument('-C', dest='chromSizes', required=True,
                   help='Tab-separated file with chromosomes and their lengths.')
parser.add_argument('-JT', dest='juicerTools', required=True,
                   help='Path of juicer-tools.jar')
parser.add_argument('-H', dest='HiC', required=True,
                   help='.HiC file.')
parser.add_argument('-R', dest='res',required=True, type = int,
                   help='Wanted resolution.')
parser.add_argument('-T', dest='tmp', default = "/tmp/",
                   help='tmp-dir.')
parser.add_argument('-O', dest='out',required=True,
                   help='Path and basename of output-file.')
parser.add_argument('-CIS', dest='CIS', default = True, type = bool,
                   help='Only output cis-interactions?')
parser.add_argument('-force', dest='force', default = False, type = bool,
                   help='Overwrite files?')
parser.add_argument('-norm', dest='norm', default = "KR", type = str,required=True,
                   help='KR or raw? (default: KR)')

args = parser.parse_args()
chromSizesFile = args.chromSizes
juicerFile = args.HiC
juicerTools = args.juicerTools
resolution = args.res
tmpDir = args.tmp
baseOut = args.out
SIGout = baseOut + ".sig"
BEDout = baseOut + ".bed"
CISbool = args.CIS
juicerDumpNorm = args.norm
if juicerDumpNorm not in ["KR", "NONE"]:
    sys.exit("!!!\tPlease choose KR or NONE in -norm\t!!!\n")

# Check if file exists: ask for force
if os.path.exists(SIGout) or os.path.exists(BEDout):
    if args.force is False:
        sys.exit("!!!\tOutput files exists.\t!!!\n!!!\tUse force=True to overwrite.\t!!!")
    else:
        void = subprocess.call(['rm',SIGout])
        void = subprocess.call(['rm',BEDout])

# Read chrom.sizes and store as dict: key==chromName, val==chromSize
chromSizesDict = dict()
with open(chromSizesFile,'r') as tsvin:
    tsvin = csv.reader(tsvin, delimiter='\t')
    for row in tsvin:
        chromSizesDict[row[0]] = row[1]

# Make the BED-file
absPandas = pd.DataFrame()
for c,s in chromSizesDict.items():
    sizeOfChrFloored = math.floor(int(s) / resolution)*resolution
    absBed = np.arange(0,sizeOfChrFloored+resolution,resolution)
    absBedStop = absBed + resolution
    absBedStop[-1] = int(s)
    absPandasCHR = pd.DataFrame({"chr" : c, "X" :absBed,"Y" : absBedStop})
    absPandas = absPandas.append(absPandasCHR)

absPandas.sort_values(["chr","X"], inplace=True) # sort on chrom,X
absPandas['idx'] = np.arange(1,len(absPandas)+1) # add index-column

# Write BED-file
absPandas = absPandas[["chr","X","Y","idx"]]
absPandas.to_csv(path_or_buf=BEDout,sep = "\t", header=False , index = False)

# run Juicer-dump and write signal-file
with open(SIGout, "a") as fp:
    wr = csv.writer(fp,dialect = 'excel-tab')
    for c1 in chromSizesDict.keys():
        absPandasC1 = absPandas[(absPandas.chr == c1)]
        C1dic = dict(zip(absPandasC1.X, absPandasC1.idx))
        for c2 in chromSizesDict.keys():
            if CISbool is True:
                if c1 is not c2:
                    continue
            absPandasC2 = absPandas[(absPandas.chr == c2)]
            C2dic = dict(zip(absPandasC2.X, absPandasC2.idx))
            # for all chromosome-chromosome combinations, run juicer-tools dump and store in tmp dir
            #print("Starting", c1,"versus", c2,"\r")
            OUT = "".join([tmpDir,c1,"_", c2,".tmp"])
            cmd = " ".join(['java', '-jar', juicerTools, 'dump','observed',juicerDumpNorm,juicerFile,c1,c2,"BP",str(resolution),OUT])
            #print(cmd)
            try:
                void = subprocess.call(cmd, shell=True)
                df = pd.read_table(OUT,delimiter="\t",header=None, names = ["pos1","pos2","signal"])
                df.dropna(inplace = True, how = 'any')
                # Get idx-numbers
                for e,row in enumerate(df.itertuples(index=True, name='Pandas')):
                    idxX = C2dic[int(row.pos1)]
                    idxY = C1dic[int(row.pos2)]
                    wr.writerow([int(idxX),int(idxY),row.signal])
                void = subprocess.call(['rm',OUT])
            except:
                continue
