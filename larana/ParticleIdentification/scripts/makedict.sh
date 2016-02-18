echo '
#pragma link C++ class anab::MVAResult+;
#pragma link C++ class vector<anab::MVAResult>+;
#pragma link C++ class map<string,double>+;
' > Linkdef.h

rootcint -f MVAResultDict.cxx -c $LARDATA_DIR/AnalysisBase/MVAResult.h Linkdef.h
g++ -shared -o MVAResultDict.so `root-config --ldflags` -fPIC -I$ROOTSYS/include MVAResultDict.cxx
rm MVAResultDict.cxx MVAResultDict.h Linkdef.h
