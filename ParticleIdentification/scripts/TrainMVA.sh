#!/bin/bash

SIGFILES=''
BKGFILES=''
JOBNAME=''

while getopts "s:b:j:" opt; do
  case $opt in
    s)
      SIGFILES="$SIGFILES $OPTARG"
      ;;
    b)
      BKGFILES="$BKGFILES $OPTARG"
      ;;
    j)
      JOBNAME=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

echo "Signal files are:"
for i in $SIGFILES;do echo $i;done
echo "Background files are:"
for i in $BKGFILES;do echo $i;done
echo "Job name is $JOBNAME"

##Build CINT dictionary for the MVAResult class
echo '
#pragma link C++ class anab::MVAResult+;
#pragma link C++ class vector<anab::MVAResult>+;
#pragma link C++ class map<string,double>+;
' > Linkdef.h

rootcint -f MVAResultDict.cxx -c $LARDATA_DIR/source/AnalysisBase/MVAResult.h Linkdef.h
g++ -shared -o MVAResultDict.so `root-config --ldflags` -fPIC -I$ROOTSYS/include MVAResultDict.cxx
rm -f MVAResultDict.cxx MVAResultDict.h Linkdef.h

ROOTSCRIPT=`mktemp`
echo "
{
gROOT->ProcessLine(\"gSystem->Load(\\\"MVAResultDict.so\\\")\");
gROOT->ProcessLine(\".L TrainMVA.C\");
std::vector<std::string> sigFiles, bkgFiles;" > $ROOTSCRIPT

for i in $SIGFILES
do
#lar -c runPID.fcl -s $i -n 1000
echo "
BuildTree(\"ntuple.root\",\"train_"`basename $i`"\");
sigFiles.push_back(\"train_"`basename $i`"\");" >> $ROOTSCRIPT
done

for i in $BKGFILES
do
#lar -c runPID.fcl -s $i -n 1000
echo "
BuildTree(\"ntuple.root\",\"train_"`basename $i`"\");
bkgFiles.push_back(\"train_"`basename $i`"\");" >> $ROOTSCRIPT
done

echo "
TrainMVA(sigFiles,bkgFiles,\"mvaPlots.root\",\"$JOBNAME\");
}
" >> $ROOTSCRIPT

cat $ROOTSCRIPT

#root -l -b -q $ROOTSCRIPT
