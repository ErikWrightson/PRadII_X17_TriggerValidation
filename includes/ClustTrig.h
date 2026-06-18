/**
 * Header file for Cluster Trigger object to evaluate for various VTP-based clustering triggers.
 * @author Erik Wrightson
 * @version 06.08.2026
 * @creation 05.28.2026
 */

#ifndef ClustTrig_H
#define ClustTrig_H

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

//Generally useful includes
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <map>

using namespace std;

class ClustTrig: public Trigger{

    public:

        //Constuctor that ensures the chain is set up via the parent constructor.
        ClustTrig(TChain* c);
        ClustTrig(TChain* c, map<string, vector<string>> m, map<Int_t, vector<Int_t>> m2, map<string, Float_t> g, map<Int_t, Float_t> g2);

        //Processing Function inhereted from parent.
        void ProcessData(bool self, bool rand, bool tSum) override;

        //Process the reconstructed data applying a individual energy threshold on each cluster.
        void ProcessData_OfflineWithThr(bool self, bool rand, bool tSum);

        //Getter methods
        /*Double_t get_Rand_LMS_Elligible();
        Double_t get_Rand_LMS_Found();
        Double_t get_tSum_LMS_Elligible();
        Double_t get_tSum_LMS_Found();
        Double_t get_tSum_LMS_BothFired();*/

        void printTSumPDF(TString pdfName);

    private:

        Double_t en;
        map<string, vector<string>> neighborMap;
        map<Int_t, vector<Int_t>> neighborMap_ModId;

        map<string, Float_t> gainMap;
        map<Int_t, Float_t> gainMap_ModId;

        //Offline 1 Cluster events (max energy cluster).
        TH1F* h_tSum_1Clust_All;
        //VTP 1 Cluster Trigger (max energy cluster).
        TH1F* h_tSum_1Clust_VTP_800;
        TH1F* h_tSum_1Clust_VTP_900;
        TH1F* h_tSum_1Clust_VTP_600;
        TH1F* h_tSum_1Clust_VTP_500;

        //Offline 1 Cluster Trigger for events where only 1 cluster was found.
        TH1F* h_tSum_1ClustOnly_All;
        //VTP 1 Cluster Trigger for events where only 1 cluster was found.
        TH1F* h_tSum_1ClustOnly_VTP_800;
        TH1F* h_tSum_1ClustOnly_VTP_900;
        TH1F* h_tSum_1ClustOnly_VTP_600;
        TH1F* h_tSum_1ClustOnly_VTP_500;

        //Offline 2 Cluster events (max energy combination).
        TH1F* h_tSum_2Clust_All;
        //VTP 2 Cluster events (max energy combination).
        TH1F* h_tSum_2Clust_VTP_800;
        //Offline 2 Cluster events for events where at least 2 clusters were found (Minimum Cluster E).
        TH1F* h_tSum_2Clust_All_MinE;
        //VTP 2 Cluster events for events where at least 2 clusters were found (Minimum Cluster E).
        TH1F* h_tSum_2Clust_VTP_800_MinE;
        //Maximum Time difference of 2 Cluster events for events where at least 2 clusters were found.
        TH1F* h_tSum_2Clust_timeDif;
        //Maximum Time difference of 2 Cluster events for events where at least 2 clusters were found for events missed by the clustering trigger.
        TH1F* h_tSum_2Clust_Missed_timeDif;
        //Maximum Cluster Time for 2 Cluster events for events where at least 2 clusters were found and the clustering trigger did not find it.
        TH1F* h_tSum_2Clust_Missed_MaxTime;
        //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where at least 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
        TH1F* h_tSum_2Clust_Missed_EdeltaT;
        //Number of Blocks for Missed Events with Delta t < 16 ns where at least 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
        TH1I* h_tSum_2Clust_Missed_nBlocks_less;
        //Number of Blocks for Missed Events with Delta t > 16 ns where at least 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
        TH1I* h_tSum_2Clust_Missed_nBlocks_greater;

        //Offline 2 Cluster events for events where only 2 clusters were found.
        TH1F* h_tSum_2ClustOnly_All;
        //VTP 2 Cluster events for events where only 2 clusters were found.
        TH1F* h_tSum_2ClustOnly_VTP_800;
        //Offline 2 Cluster events for events where only 2 clusters were found (Minimum Cluster E).
        TH1F* h_tSum_2ClustOnly_All_MinE;
        //VTP 2 Cluster events for events where only 2 clusters were found (Minimum Cluster E).
        TH1F* h_tSum_2ClustOnly_VTP_800_MinE;
        //Time difference of 2 Cluster events for events where only 2 clusters were found.
        TH1F* h_tSum_2ClustOnly_timeDif;
        //Maximum Time difference of 2 Cluster events for events where only 2 clusters were found for events missed by the clustering trigger.
        TH1F* h_tSum_2ClustOnly_Missed_timeDif;
        //Maximum Cluster Time for 2 Cluster events for events where only 2 clusters were found and the clustering trigger did not find it.
        TH1F* h_tSum_2ClustOnly_Missed_MaxTime;
        //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where only 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
        TH1F* h_tSum_2ClustOnly_Missed_EdeltaT;
        //Number of Blocks for Missed Events with Delta t < 16 ns where only 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
        TH1I* h_tSum_2ClustOnly_Missed_nBlocks_less;
        //Number of Blocks for Missed Events with Delta t > 16 ns where only 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
        TH1I* h_tSum_2ClustOnly_Missed_nBlocks_greater;

        //Offline 3 Cluster events (max energy combination).
        TH1F* h_tSum_3Clust_All;
        //VTP 3 Cluster events (max energy combination).
        TH1F* h_tSum_3Clust_VTP_800;
        //Offline 3 Cluster events for events where at least 3 clusters were found (Minimum Cluster E).
        TH1F* h_tSum_3Clust_All_MinE;
        //VTP 3 Cluster events for events where at least 3 clusters were found (Minimum Cluster E).
        TH1F* h_tSum_3Clust_VTP_800_MinE;
        //Maximum Time difference of 3 Cluster events for events where at least 3 clusters were found.
        TH1F* h_tSum_3Clust_timeDif;
        //Maximum Time difference of 3 Cluster events for events where at least 3 clusters were found for events missed by the clustering trigger.
        TH1F* h_tSum_3Clust_Missed_timeDif;
        //Maximum Cluster Time for 3 Cluster events for events where at least 3 clusters were found and the clustering trigger did not find it.
        TH1F* h_tSum_3Clust_Missed_MaxTime;
        //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where at least 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
        TH1F* h_tSum_3Clust_Missed_EdeltaT;
        //Number of Blocks for Missed Events with Delta t < 16 ns where at least 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
        TH1I* h_tSum_3Clust_Missed_nBlocks_less;
        //Number of Blocks for Missed Events with Delta t > 16 ns where at least 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
        TH1I* h_tSum_3Clust_Missed_nBlocks_greater;

        //Offline 3 Cluster events for events where only 3 clusters were found.
        TH1F* h_tSum_3ClustOnly_All;
        //VTP 3 Cluster events for events where only 3 clusters were found.
        TH1F* h_tSum_3ClustOnly_VTP_800;
        //Offline 3 Cluster events for events where only 3 clusters were found (Minimum Cluster E).
        TH1F* h_tSum_3ClustOnly_All_MinE;
        //VTP 3 Cluster events for events where only 3 clusters were found (Minimum Cluster E).
        TH1F* h_tSum_3ClustOnly_VTP_800_MinE;
        //Time difference of 3 Cluster events for events where only 3 clusters were found.
        TH1F* h_tSum_3ClustOnly_timeDif;
        //Maximum Time difference of 3 Cluster events for events where only 3 clusters were found for events missed by the clustering trigger.
        TH1F* h_tSum_3ClustOnly_Missed_timeDif;
        //Maximum Cluster Time for 3 Cluster events for events where only 3 clusters were found and the clustering trigger did not find it.
        TH1F* h_tSum_3ClustOnly_Missed_MaxTime;
        //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where only 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
        TH1F* h_tSum_3ClustOnly_Missed_EdeltaT;
        //Number of Blocks for Missed Events with Delta t < 16 ns where at only 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
        TH1I* h_tSum_3ClustOnly_Missed_nBlocks_less;
        //Number of Blocks for Missed Events with Delta t > 16 ns where only 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
        TH1I* h_tSum_3ClustOnly_Missed_nBlocks_greater;

        //Offline 2 Cluster events; HyCal XY Positions
        TH2F* h_tSum_2Clust_All_HC_XY;
        //VTP 2 Cluster events; HyCal XY Positions
        TH2F* h_tSum_2Clust_VTP_HC_XY;
        //2 Cluster events missed by the VTP Trigger; HyCal XY, and Theta Positions
        TH2F* h_tSum_2Clust_Miss_HC_XY;
        TH1F* h_tSum_2Clust_Miss_HC_Theta;

        //Number of VTP Clusters - Offline Clusters
        TH1I* h_tSum_VTP_NClustDiff;
        

        //Offline 3 Cluster events; HyCal XY Positions
        TH2F* h_tSum_3Clust_All_HC_XY;
        //VTP 3 Cluster events; HyCal XY Positions
        TH2F* h_tSum_3Clust_VTP_HC_XY;
        //3 Cluster events missed by the VTP Trigger; HyCal XY, and Theta Positions
        TH2F* h_tSum_3Clust_Miss_HC_XY;
        TH1F* h_tSum_3Clust_Miss_HC_Theta;

        void setup_Histos();

        /*Double_t tSum_LMS_Elligible;
        Double_t tSum_LMS_Found;
        Double_t tSum_LMS_BothFire;*/

};

#endif