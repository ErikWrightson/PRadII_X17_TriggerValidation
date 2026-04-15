/**
 * Trigger Validation Processor
 * @author Erik Wrightson
 * @version 04.07.2026
 * @creation 04.05.2026
 */

#include <iostream>
#include <string.h>
#include <getopt.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>


//ROOT Includes that may be handy to have.
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TObject.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TFitResult.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TColor.h>
#include <TString.h>
//Needed for reading in vector types from root files.
#include <TInterpreter.h>

#include "includes/Trigger.h"
#include "includes/LMSTrig.h"


using namespace std;

/**
 * Processes the file with the input file name and reads out the root file names contained within and adds them to a vector.
 *
 * @param fileListFileName - the name of the file with the list of ROOT file names.
 */
vector<TString> processFileList(string fileListFileName){
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
TChain* makeChain(vector<TString> names){
        TChain* chain = new TChain("events");
        for(unsigned int i = 0; i < names.size(); i++){
            chain->Add(names.at(i));
        }

        return chain;
}


/**
 * Prints out the proper usage directions for this program and what each flag means.
 *
 * @param prog - the progam name that is being currently run.
 */
static void printUsage(const char *prog)
{
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
              << "\t-e Evaluates this trigger's performance with respect to itself. Only checks for internal performance on fired events."
              << "\t-h Show this help\n"
              << "\tNOTE: Either option -f or -L are REQUIRED for running properly.\n";
}

/**
 * The main function that launches the trigger analysis.
 *
 * @param argc - the number of input arguments
 * @param argv - an array of the different arguments as an array of char* (strings).
 */
int main (int argc, char **argv){

    bool all = false;
    bool Lms = false;
    bool sum = false;
    bool alpha = false;
    bool mOR = false;
    bool vtp_clust = false;
    bool comp_TotalSum = false;
    bool rand = false;
    bool self = false;

    string fileName;
    string fileListFileName;

	if (argc<2) {
		cout<<"ERR: Incorrect Arguments: " <<endl;
        printUsage(argv[0]);
		
		return -1;
	}

    // ── Parse command-line ───────────────────────────────────────────────
    int opt;
    while ((opt = getopt(argc, argv, "alspmchTf:rL:e")) != -1) {
        switch (opt) {
            case 'a': all=true; Lms =true; sum=true; alpha=true; mOR=true; vtp_clust = true; break;
            case 'l': Lms = true; break;
            case 's': sum = true; break;
            case 'p': alpha = true; break;
            case 'm': mOR = true; break;
            case 'c': vtp_clust = true; break;
            case 'T': comp_TotalSum = true; break;
            case 'R': rand = true; break;
            case 'f': fileName = optarg; break;
            case 'L': fileListFileName = optarg; break;
            case 'e': self = true; break;
            case 'h':
            default: printUsage(argv[0]); return (opt == 'h') ? 0 : 1;
        }
    }

    struct stat buffer;   
    bool existOne = (stat(fileName.c_str(), &buffer) == 0);
    bool existList = (stat(fileListFileName.c_str(), &buffer) == 0);

    if(!existOne && !existList){
        cerr<<"A single valid input file or a filelist txt file was not provided.\n";
        return -3;
    }

    vector<TString> fileNameVec;
    if(existList){
        fileNameVec = processFileList(fileListFileName);
    }
    if(existOne){
        fileNameVec.push_back((TString) fileName);
    }

    if(!(comp_TotalSum || rand || vtp_clust)){
        cerr<<"You must either choose a comparison trigger to use (Random trigger [-R], or Total Sum trigger [-T]), or be evaluating only the Clusters found by the VTP trigger amongst themselves.";
        return -4;
    }

    if(comp_TotalSum && !(Lms || vtp_clust) && !all){
        cerr<<"Getting trigger efficiencies from Total Sum trigger events requires that it be one of the event types this has been implemented for.\n"
            <<"Currently the following options are available:\n\t-L for the LMS option\n\t-c for the VTP clustering trigger.\n";
        return -2;
    }

    TChain* fChain = makeChain(fileNameVec);

    if(Lms){
        LMSTrig trig1 = LMSTrig(fChain);

        trig1.ProcessData(self, rand, comp_TotalSum);

        if(comp_TotalSum){
            Double_t lms_ts_f = trig1.get_tSum_LMS_Found();
            Double_t lms_ts_e = trig1.get_tSum_LMS_Elligible();
            Double_t lms_ts_eff = lms_ts_f/lms_ts_e * 100;

            Double_t lms_ts_both = trig1.get_tSum_LMS_BothFired();

            cout<< "TSum and LMS  fired: " << lms_ts_f <<" TSum that were LMS Elligible: " << lms_ts_e << " LMS Trigger Efficiency from TSum: " << lms_ts_eff << "%\n";
            cout<<"Total amount of time TSum and LMS fired together: " << lms_ts_both << endl;
        }
    }

    
    

    return 0;
}