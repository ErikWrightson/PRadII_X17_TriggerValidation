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

#include "Trigger.h"

using namespace std;

class ClustTrig: public Trigger{

    public:

        struct randInformation{
            vector<array<double, Trigger::nSamples>> doubleInfoVecs;
            vector<array<int, Trigger::nSamples>> integerInfoVecs;
            array<vector<int>,Trigger::nSamples> allIndexVectors;
        };

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
        void delete_tSum_Histos();
        void SaveAllTSumHistograms(TString rootfileName);

        void printRandPDF(TString pdfName);
        void delete_rand_Histos();
        void SaveAllRandHistograms(TString rootfileName);

        void setUseMax(bool set);

    private:

        Double_t en;
        Long64_t entries;
        map<string, vector<string>> neighborMap;
        map<Int_t, vector<Int_t>> neighborMap_ModId;

        map<string, Float_t> gainMap;
        map<Int_t, Float_t> gainMap_ModId;

        vector<Int_t> excludeMod;

        Double_t Xedge[35], Yedge[35];
        bool usemax;
        
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
        //Individual cluster energy for at least 3 or more cluster events
        TH1F* h_tSum_3Clust_Cl_Energy_All;
        TH1F* h_tSum_3Clust_Cl_Energy_Trig;
        //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
        TH1F* h_tSum_3Clust_GEMMatch_All;
        TH1F* h_tSum_3Clust_GEMMatch_Trig;
        //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
        TH1F* h_tSum_3Clust_GEMMatch_All_2Match_Either;
        TH1F* h_tSum_3Clust_GEMMatch_Trig_2Match_Either;
        //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
        TH1F* h_tSum_3Clust_GEMMatch_All_3Match_Either;
        TH1F* h_tSum_3Clust_GEMMatch_Trig_3Match_Either;
        //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
        TH1F* h_tSum_3Clust_GEMMatch_All_3Match_Both;
        TH1F* h_tSum_3Clust_GEMMatch_Trig_3Match_Both;
        //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
        TH1F* h_tSum_3Clust_GEMMatch_All_Cl_E;
        TH1F* h_tSum_3Clust_GEMMatch_Trig_Cl_E;
        //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
        TH1F* h_tSum_3Clust_GEMMatch_All_2Match_Either_Cl_E;
        TH1F* h_tSum_3Clust_GEMMatch_Trig_2Match_Either_Cl_E;
        //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
        TH1F* h_tSum_3Clust_GEMMatch_All_3Match_Either_Cl_E;
        TH1F* h_tSum_3Clust_GEMMatch_Trig_3Match_Either_Cl_E;
        //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
        TH1F* h_tSum_3Clust_GEMMatch_All_3Match_Both_Cl_E;
        TH1F* h_tSum_3Clust_GEMMatch_Trig_3Match_Both_Cl_E;
        //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches and a high/low cluster energy sum.
        TH1F* h_tSum_3Clust_vZ_High;
        TH1F* h_tSum_3Clust_vZ_Low;
        //Maximum cluster energy for only 3 cluster events
        TH1F* h_tSum_3Clust_Cl_Max_All;
        TH1F* h_tSum_3Clust_Cl_Max_Trig;


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
        //Individual cluster energy for only 3 cluster events
        TH1F* h_tSum_3ClustOnly_Cl_Energy_All;
        TH1F* h_tSum_3ClustOnly_Cl_Energy_Trig;
        //Maximum cluster energy for only 3 cluster events
        TH1F* h_tSum_3ClustOnly_Cl_Max_All;
        TH1F* h_tSum_3ClustOnly_Cl_Max_Trig;
        //
        TH1F* h_tSum_3Clust_DCA_Low;
        TH1F* h_tSum_3Clust_delta_X_Low;
        TH1F* h_tSum_3Clust_DCA_High;
        TH1F* h_tSum_3Clust_delta_X_High;

        //Offline 2 Cluster events; HyCal XY Positions
        TH2F* h_tSum_2Clust_All_HC_XY;
        //VTP 2 Cluster events; HyCal XY Positions
        TH2F* h_tSum_2Clust_VTP_HC_XY;
        //2 Cluster events missed by the VTP Trigger; HyCal XY, and Theta Positions
        TH2F* h_tSum_2Clust_Miss_HC_XY;
        TH1F* h_tSum_2Clust_Miss_HC_Theta;
        TH2F* h_tSum_2Clust_Miss_HC_XY_LowE;

        //Number of VTP Clusters
        TH1I* h_tSum_VTP_NClustDiff;
        //Number of VTP Clusters - Offline Clusters for rejected events.
        TH1I* h_tSum_VTP_NClustDiff_Missed;
        //VTP timing of events with non-zero difference in cluster times
        TH1D* h_tSum_VTP_Timing;
        
        //VTP center cluster module ID
        TH1I* h_tSum_VTP_centerID;
        //Offline center cluster module ID
        TH1I* h_tSum_Off_centerID;
        

        //Offline 3 Cluster events; HyCal XY Positions
        TH2F* h_tSum_3Clust_All_HC_XY;
        //VTP 3 Cluster events; HyCal XY Positions
        TH2F* h_tSum_3Clust_VTP_HC_XY;
        //3 Cluster events missed by the VTP Trigger; HyCal XY, and Theta Positions
        TH2F* h_tSum_3Clust_Miss_HC_XY;
        TH1F* h_tSum_3Clust_Miss_HC_Theta;
        TH2F* h_tSum_3Clust_Miss_HC_XY_LowE;



        TH1F* h_rand_3Clust_All_Sum;
        TH1F* h_rand_3Clust_VTP_Sum[nSSPBits];
        TH1F* h_rand_3Clust_Sum_Ratio[nSSPBits];

        TH1F* h_rand_3Clust_All_Min;
        TH1F* h_rand_3Clust_VTP_Min[nSSPBits];
        TH1F* h_rand_3Clust_Min_Ratio[nSSPBits];


        TH1F* h_rand_3Clust_All_Max;
        TH1F* h_rand_3Clust_VTP_Max[nSSPBits];
        TH1F* h_rand_3Clust_Max_Ratio[nSSPBits];

        TH2F* h_rand_3Clust_timeCorrEarly[nSSPBits];
        TH2F* h_rand_3Clust_timeCorrLate[nSSPBits];

        TH2F* h_rand_3Clust_XY_All;
        TH2F* h_rand_3Clust_XY_VTP[nSSPBits];

        TH1F* h_rand_3Clust_deltaT_All;
        TH1F* h_rand_3Clust_deltaT_VTP[nSSPBits];
        TH1F* h_rand_3Clust_deltaT_Ratio[nSSPBits];


        void setup_tSum_Histos();

        void fill_tSum_Histos(Int_t i);
        void fill_tSum_Histos_wThr(Int_t i);

        void setup_rand_Histos();

        void fill_rand_Histos_wThr(Int_t i);

        randInformation ComputeTimeBinnedClusters(bool appThr);

        /*Double_t tSum_LMS_Elligible;
        Double_t tSum_LMS_Found;
        Double_t tSum_LMS_BothFire;*/

};

#endif