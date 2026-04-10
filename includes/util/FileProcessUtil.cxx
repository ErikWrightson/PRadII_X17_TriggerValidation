#include "FileProcessUtil.h"


namespace FileProcessUtil{
    
    static vector<TString> processFileList(string fileListFileName){
        ifstream file(fileListFileName);
        string line;
        
        vector<TString> list;
        while(getline(file,line)){
            TString l(line);
            list.push_back(l);
        }

        return list;
    }

    TChain* makeChain(vector<TString> names){
        TChain* chain = new TChain("T");
        for(unsigned int i = 0; i < names.size(); i++){
            chain->Add(names.at(i));
        }

        return chain;
    }
}