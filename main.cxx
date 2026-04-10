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
#include "includes/util/FileProcessUtil.h"

using namespace std;
using namespace FileProcessUtil;

static void printUsage(const char *prog)
{
    cerr << "Usage: " << prog << " [options]\n"
              << "\t-a Evaluates for all trigger type efficiencies\n";
              << "\t-l Evaluates for LMS Trigger efficiency\n"
              << "\t-s Evaluates for the Total Sum Trigger efficiency\n"
              << "\t-p Evaluates for the Alpha Source Trigger efficiency\n"
              << "\t-m Evaluates for the Master OR (T10) Trigger efficiency\n"
              << "\t-c Evaluates for the VTP Cluster Trigger efficiency\n"
              << "\t-T Evaluates for trigger efficiencies compared to Total Sum events\n"
              << "\t-R Evaluates for trigger efficiencies from events with random trigger\n"
              << "\t-f <filename> of the .root file to evaluate for the various trigger efficiencies\n"
              << "\t-L <fileList.txt> of the .root files to evaluate. Use if wanting to process multiple files. Follow the format in the README\n"
              << "\t-h Show this help\n"
              << "\tNOTE: Either option -f or -L are REQUIRED for running properly.";
}

int main (int argc, char **argv){

    bool all = false;
    bool Lms = false;
    bool sum = false;
    bool alpha = false;
    bool mOR = false;
    bool vtp_clust = false;
    bool comp_TotalSum = false;
    bool rand = false;

    TString fileName;
    string fileListFileName;

	if (argc<1) {
		cout<<"ERR: Incorrect Arguments: "<<endl;
		
		return -1;
	}

    // ── Parse command-line ───────────────────────────────────────────────
    int opt;
    while ((opt = getopt(argc, argv, "alspmchf:rL:")) != -1) {
        switch (opt) {
            case 'a': all=true; Lms =true; sum=true; alpha=true; mOr=true; vtp_clust = true; break;
            case 'l': Lms = true; break;
            case 's': sum = true; break;
            case 'p': alpha = true; break;
            case 'm': mOr = true; break;
            case 'c': vtp_clust = true; break;
            case 'T': comp_TotalSum = true; break;
            case 'R': rand = true; break;
            case 'f': fileName = optarg; break;
            case 'L': fileListFileName = optarg; break;
            case 'h':
            default: printUsage(argv[0]); return (opt == 'h') ? 0 : 1;
        }
    }

    struct stat buffer;   
    bool existOne = (stat(fileName.c_str(), &buffer) == 0);
    bool existList = (stat(fileListFileName.c_str(), &buffer) == 0);

    if(!existOne && !existList){
        cerr<<"A single valid input file or a filelist txt file was not provided.";
        return -3;
    }

    vector<TString> fileVec;
    if(existList){
        fileNameVec = processFileList();
    }
    if(existOne){
        fileNameVec.push_back(fileName);
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
        
    }

    return 0;
}