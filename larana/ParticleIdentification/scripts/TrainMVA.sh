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

#Hmm, when using mrb then $LARANA_DIR points to source tree, but
#in installed version then need to add source to path
if [ -d $LARANA_DIR/source ]; then
    cp $LARANA_DIR/source/ParticleIdentification/scripts/TrainMVA.C ./
    cp $LARANA_DIR/source/ParticleIdentification/scripts/runPID_noMVA.fcl ./
else
    cp $LARANA_DIR/ParticleIdentification/scripts/TrainMVA.C ./
    cp $LARANA_DIR/ParticleIdentification/scripts/runPID_noMVA.fcl ./
fi

if [ -d $LARDATA_DIR/source ]; then
    cp $LARDATA_DIR/source/AnalysisBase/MVAPIDResult.h ./
else
    cp $LARDATA_DIR/AnalysisBase/MVAPIDResult.h ./
fi


##Build CINT dictionary for the MVAResult class
echo '
#pragma link C++ class anab::MVAPIDResult+;
#pragma link C++ class vector<anab::MVAPIDResult>+;
#pragma link C++ class map<string,double>+;
' > Linkdef.h

rootcint -f MVAResultDict.cxx -c $LARDATA_DIR/AnalysisBase/MVAPIDResult.h Linkdef.h
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
lar -c ./runPID_noMVA.fcl -s $i -n 100000
TUPLEFILE=ntuple_`basename $i`
mv ntuple.root $TUPLEFILE
echo "
BuildTree(\"${TUPLEFILE}\",\"train_"`basename $i`"\");
sigFiles.push_back(\"train_"`basename $i`"\");" >> $ROOTSCRIPT
done

for i in $BKGFILES
do
mv ntuple.root ntuple_$i
lar -c ./runPID_noMVA.fcl -s $i -n 100000
TUPLEFILE=ntuple_`basename $i`
mv ntuple.root $TUPLEFILE
echo "
BuildTree(\"${TUPLEFILE}\",\"train_"`basename $i`"\");
bkgFiles.push_back(\"train_"`basename $i`"\");" >> $ROOTSCRIPT
done

echo "
TrainMVA(sigFiles,bkgFiles,\"mvaPlots.root\",\"$JOBNAME\");
}
" >> $ROOTSCRIPT

cat $ROOTSCRIPT

root -l -b -q $ROOTSCRIPT
