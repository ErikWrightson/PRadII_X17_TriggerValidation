/**
 * Contains Utility functions that can be used independently.
 * 
 * @author Erik Wrightson <wrightso@jlab.org>
 * @version 06.15.2026
 * @creation 06.15.2026
 */
#include "Utils.h"

/**
 * Prints out the proper usage directions for this program and what each flag means.
 *
 * @param prog - the progam name that is being currently run.
 */
void Utils::printUsage(const char *prog){
    cerr << "Usage: " << prog << " [options]\n"
              << "\t-a Evaluates for all trigger type efficiencies\n"
              << "\t-l Evaluates for LMS Trigger efficiency\n"
              << "\t-s Evaluates for the Total Sum Trigger efficiency\n"
              << "\t-p Evaluates for the Alpha Source Trigger efficiency\n"
              << "\t-m Evaluates for the Master OR (T10) Trigger efficiency\n"
              << "\t-c Evaluates for the VTP Cluster Trigger efficiency\n"
              << "\t-T Evaluates for trigger efficiencies compared to Total Sum events\n"
              << "\t-R Evaluates for trigger efficiencies from events with random trigger\n"
              << "\t-f <filename> of the .root file to evaluate for the various trigger efficiencies\n"
              << "\t-L <fileList.txt> of the .root files to evaluate. Use if wanting to process multiple files. Follow the format in the README\n"
              << "\t-e Evaluates this trigger's performance with respect to itself. Only checks for internal performance on fired events.\n"
              << "\t-N Evaluates this trigger's performance using reconstructed values.\n"
              << "\t-h Show this help\n"
              << "\tNOTE: Either option -f or -L are REQUIRED for running properly.\n";
}

/**
 * Processes the file with the input file name and reads out the root file names contained within and adds them to a vector.
 *
 * @param fileListFileName - the name of the file with the list of ROOT file names.
 */
vector<TString> Utils::processFileList(string fileListFileName){
        ifstream file(fileListFileName);
        string line;
        
        vector<TString> list;
        while(getline(file,line)){
            TString l(line);
            list.push_back(l);
        }

        return list;
}

/**
 * Makes a TChain of the entries in a TString vector assuming they are valid paths to ROOT files.
 *
 * @param names - the vector of ROOT file names to be linked in the chain.
 */
TChain* Utils::makeChain(vector<TString> names, TString treeName){
        TChain* chain = new TChain(treeName);//"recon");//"events");
        for(unsigned int i = 0; i < names.size(); i++){
            chain->Add(names.at(i));
        }

        return chain;
}

void Utils::makeNeighborMap(string dbName, map<string, vector<string>> &m, map<Int_t, vector<Int_t>> &m2){
    ifstream file(dbName);
    string line;

    while(getline(file,line)){
        stringstream ss(line);
        string word;
        
        Int_t wordNum = 0;
        string key;
        vector<string> val;

        Int_t key2;
        vector<Int_t> val2;

        while (ss >> word) {
            if(wordNum == 0){
                key = word;
                string temp = word;
                if(temp.front() == 'G'){
                    temp.erase(0,1);
                    key2 = stoi(temp);
                }
                else{
                    temp.erase(0,1);
                    key2 = 1000+stoi(temp);
                }
            }
            if(wordNum != 0){
                if(wordNum == 1){
                    word.erase(0,2);
                    word.erase(word.length()-2, word.length());
                }
                else if(word.back() == ']'){
                    word.erase(0,1);
                    word.erase(word.length()-2, word.length());
                }
                else{
                    word.erase(0,1);
                    word.erase(word.length()-2, word.length());
                }
                val.push_back(word);

                string temp = word;
                if(temp.front() == 'G'){
                    temp.erase(0,1);
                    val2.push_back(stoi(temp));
                }
                else{
                    temp.erase(0,1);
                    val2.push_back(1000+stoi(temp));
                }
            }
            wordNum++;
        }
        m[key] = val;
        m2[key2] = val2;
    }

    cout<<"Made the quick reference neighbor map.\n";
}

void Utils::makeGainMap(string dbName, map<string, Float_t> &m, map<Int_t, Float_t> &m2){
    ifstream file(dbName);
    string line;

    while(getline(file,line)){
        stringstream ss(line);
        string word;
        
        Int_t wordNum = 0;
        string key;
        Float_t val;

        Int_t key2;

        while (ss >> word) {
            if(wordNum == 0){
                key = word;
                string temp = word;
                if(temp.front() == 'G'){
                    temp.erase(0,1);
                    key2 = stoi(temp);
                }
                else{
                    temp.erase(0,1);
                    key2 = 1000+stoi(temp);
                }
            }
            else{
                val = stof(word);
            }
            wordNum++;
        }
        m[key] = val;
        m2[key2] = val;
    }

    cout<<"Made the quick reference gain map.\n";
}