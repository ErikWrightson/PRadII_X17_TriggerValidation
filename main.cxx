/**
 * Trigger Validation Processor
 * @author Erik Wrightson
 * @version 06.08.2026
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
#include "includes/ClustTrig.h"
#include "includes/Utils.h"


using namespace std;

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
    bool recon = false;
    bool applyOfflineCuts = true;

    string fileName;
    string fileListFileName;

    TString tName = "events";
    TString outputDirectory = "outfiles/";
    TString rootOutputDirectory = "rootOutfiles/";
    TString fn = "ClusteringTrigger_wVTPInfo_25376_wCuts_Partial";

	if (argc<2) {
		cout<<"ERR: Incorrect Arguments: " <<endl;
        Utils::printUsage(argv[0]);
		
		return -1;
	}
    

    // ── Parse command-line ───────────────────────────────────────────────
    int opt;
    while ((opt = getopt(argc, argv, "alspmchTf:rL:eD:No:C")) != -1) {
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
            case 'D': outputDirectory = optarg; break;
            case 'o': fn = optarg; break;
            case 'N': recon = true; break;
            case 'e': self = true; break;
            case 'C': applyOfflineCuts = false; break;
            case 'h':
            default: Utils::printUsage(argv[0]); return (opt == 'h') ? 0 : 1;
        }
    }

    struct stat buffer;   
    bool existOne = (stat(fileName.c_str(), &buffer) == 0);
    bool existList = (stat(fileListFileName.c_str(), &buffer) == 0);

    if(recon){
        tName = "recon";
    }

    if(!existOne && !existList){
        cerr<<"A single valid input file or a filelist txt file was not provided.\n";
        return -3;
    }

    vector<TString> fileNameVec;
    if(existList){
        fileNameVec = Utils::processFileList(fileListFileName);
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

    TChain* fChain = Utils::makeChain(fileNameVec, tName);

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

    if(vtp_clust){

        if(recon){
            ClustTrig Cl_trig = ClustTrig(fChain);
            if(applyOfflineCuts){
                Cl_trig.ProcessData_OfflineWithThr(self, rand, comp_TotalSum);
            }
            else{
                Cl_trig.ProcessData(self, rand, comp_TotalSum);
            }
            
            if(comp_TotalSum){
                Cl_trig.printTSumPDF(outputDirectory + fn + "_TSum.pdf");
                Cl_trig.SaveAllTSumHistograms(rootOutputDirectory + fn +"_TSum.root");
                Cl_trig.delete_tSum_Histos();
            }

            if(rand){
                Cl_trig.printRandPDF(outputDirectory + fn + "_wMax_Rand.pdf");
                Cl_trig.SaveAllRandHistograms(rootOutputDirectory + fn +"_wMax_Rand.root");
                Cl_trig.delete_rand_Histos();

                if(applyOfflineCuts){
                    Cl_trig.setUseMax(false);
                    Cl_trig.ProcessData_OfflineWithThr(self, rand, comp_TotalSum);
                    Cl_trig.printRandPDF(outputDirectory + fn + "_NoMax_Rand.pdf");
                    Cl_trig.SaveAllRandHistograms(rootOutputDirectory + fn +"_NoMax_Rand.root");
                    Cl_trig.delete_rand_Histos();
                }
            }
        }
        else{
            string neighborDb = "database/module_neighbors.txt";
            map<string, vector<string>> nBMap;
            map<Int_t, vector<Int_t>> nBMap_ModId;
            Utils::makeNeighborMap(neighborDb, nBMap, nBMap_ModId);

            string gDb = "database/calibration_factor_3p5_June7.txt";
            map<string, Float_t> gMap;
            map<Int_t, Float_t> gMap_ModId;
            Utils::makeGainMap(gDb, gMap, gMap_ModId);

            ClustTrig Cl_trig_Raw = ClustTrig(fChain, nBMap, nBMap_ModId, gMap, gMap_ModId);
        }
    }

    if(mOR){

    }

    if(alpha){

    }

    if(sum){
        
    }

    
    

    return 0;
}