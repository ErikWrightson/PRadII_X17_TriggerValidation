/**
 * Form an object with the proper functionality to find and make all the plots for Clustering Trigger Evaluation.
 *
 * @author Erik Wrightson
 * @version 06.08.2026
 * @creation 05.28.2026
 */
#include "Trigger.h"
#include "ClustTrig.h"

/**
 * Creates the Clustering trigger object for analysis on the reconstructed clusters.
 *
 * @param c - the TChain of events to
 */
ClustTrig::ClustTrig(TChain* c):Trigger(c, false, true){
    cout<<"Set All Branch Addresses for Clustering Trigger Object.\n";

    en = 3485.41;

    setup_Histos();
}

/**
 * Creates the LMS trigger object.
 *
 * @param c - the TChain of events to
 */
ClustTrig::ClustTrig(TChain* c, map<string, vector<string>> m, map<Int_t, vector<Int_t>> m2, map<string, Float_t> g, map<Int_t, Float_t> g2):Trigger(c, false){
    cout<<"Set All Branch Addresses for Clustering Trigger Object.\n";

    neighborMap = m;
    neighborMap_ModId = m2;

    gainMap = g;
    gainMap_ModId = g2;
    
    en = 3485.41;

}

/**
 * Processes the events one by one depending on which of the flags are set to be true.
 * For the LMS it looks for efficiencies when other trigger are fired and also evaluates
 * if any channels were dropped.
 *
 * @param self - Whether or not this trigger should have evaluation done within its own events. {True for Yes; False for No}
 * @param rand - Whether or not this trigger should have evalutaion done using a random trigger source. {True for Yes; False for No}
 * @param tSum - Whether or not this trigger should have evaluation done using total sum trigger source. {True for Yes; False for No}
 */
void ClustTrig::ProcessData(bool self, bool rand, bool tSum){
    
    //Get the amount of entries from each file to limit looping through them.
	Long64_t entries = chain->GetEntries();
    
    for(Int_t i = 0; i < entries; i++){
        chain->GetEntry(i);

        trig.Parse(*sspRawPtr);

        // Fill per-bit TRIG_TIME histograms - from Rafo's code
        for (int b = 0; b < nSSPBits; ++b){
            const vector<int>* times = trig.GetTriggersBit(b);
            if (!times) continue;
            for (int t : *times){
                hTrigTime[b]->Fill(t);
            }
        }


        if(self && (trigger_bits & (1<<CL1_FLAG))){

        }
        if(rand && (trigger_bits & (1<<RAND_FC_FLAG))){

        }
        //cout<<trigger_bits<<endl;
        if(tSum && (!trig.GetTriggersBit(TSUM_FLAG-nSSPBits)->empty()) && !(trigger_bits & (1<<LMS_FLAG))){ //(trigType & (1<<SSP_RAWSUM_TFLAG))

            //cout<<trigger_bits<<endl;

            Double_t clustSum = 0;
            Double_t maxTime = -100000;
            Double_t minTime = 10000000;
            Double_t clust_min = 1000000;
            UChar_t clust_min_bl = 244;
            for(Int_t j = 0; j < nClust; j++){
                if((i+1)%10000 == 0 || (entries-i+1)<10000){
                    cout<<"\rClustering Trigger Total Sum Events: " << i+1 << "/" << entries << flush;
                    if(i+1 == entries){
                        cout<<endl;
                    }
                }
                if(cl_E[j] > CL_IND_THR){
                    clustSum += cl_E[j];
                }
                if(cl_time[j] > maxTime){
                    maxTime = cl_time[j];
                }
                if(cl_time[j] < minTime){
                    minTime = cl_time[j];
                }
                if(cl_E[j] < clust_min){
                    clust_min = cl_E[j];
                    clust_min_bl = cl_nblocks[j];
                }

                Double_t theta = TMath::ATan2(TMath::Sqrt(cl_x[j]*cl_x[j]+cl_y[j]*cl_y[j]),cl_z[j])*rad2Deg;

                if(nClust >= 2){
                    if(!trig.GetTriggersBit(CL2_FLAG-nSSPBits)->empty()){
                        h_tSum_2Clust_VTP_HC_XY->Fill(cl_x[j], cl_y[j]);
                    }
                    else{
                        h_tSum_2Clust_Miss_HC_XY->Fill(cl_x[j], cl_y[j]);
                        h_tSum_2Clust_Miss_HC_Theta->Fill(theta);
                    }
                    h_tSum_2Clust_All_HC_XY->Fill(cl_x[j], cl_y[j]);
                }

                if(nClust >= 3){
                    if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                        h_tSum_3Clust_VTP_HC_XY->Fill(cl_x[j], cl_y[j]);
                    }
                    else{
                        h_tSum_3Clust_Miss_HC_XY->Fill(cl_x[j], cl_y[j]);
                        h_tSum_3Clust_Miss_HC_Theta->Fill(theta);
                    }
                    h_tSum_3Clust_All_HC_XY->Fill(cl_x[j], cl_y[j]);
                }
            }

            h_tSum_VTP_NClustDiff->Fill(vtp_cl_n - nClust);

            if(nClust >= 1){

                if(nClust == 1){
                    if(!trig.GetTriggersBit(CL1_FLAG-nSSPBits)->empty()){//trigger_bits & (1<<CL1_FLAG)){
                        h_tSum_1ClustOnly_VTP_800->Fill(clustSum);
                    }
                    if(!trig.GetTriggersBit(CL1_500_FLAG-nSSPBits)->empty()){
                        h_tSum_1ClustOnly_VTP_500->Fill(clustSum);
                    }
                    if(!trig.GetTriggersBit(CL1_600_FLAG-nSSPBits)->empty()){
                        h_tSum_1ClustOnly_VTP_600->Fill(clustSum);
                    }
                    if(!trig.GetTriggersBit(CL1_900_FLAG-nSSPBits)->empty()){
                        h_tSum_1ClustOnly_VTP_900->Fill(clustSum);
                    }
                    h_tSum_1ClustOnly_All->Fill(clustSum);
                }

               if(!trig.GetTriggersBit(CL1_FLAG-nSSPBits)->empty()){
                    h_tSum_1Clust_VTP_800->Fill(clustSum);
               }
               if(!trig.GetTriggersBit(CL1_500_FLAG-nSSPBits)->empty()){
                    h_tSum_1Clust_VTP_500->Fill(clustSum);
               }
               if(!trig.GetTriggersBit(CL1_600_FLAG-nSSPBits)->empty()){
                    h_tSum_1Clust_VTP_600->Fill(clustSum);
               }
               if(!trig.GetTriggersBit(CL1_900_FLAG-nSSPBits)->empty()){
                    h_tSum_1Clust_VTP_900->Fill(clustSum);
               }
               h_tSum_1Clust_All->Fill(clustSum);
               
               if(nClust >= 2){

                    if(nClust == 2){
                        if(!trig.GetTriggersBit(CL2_FLAG-nSSPBits)->empty()){
                            h_tSum_2ClustOnly_VTP_800->Fill(clustSum);
                            h_tSum_2ClustOnly_VTP_800_MinE->Fill(clust_min);
                        }
                        else{
                            h_tSum_2ClustOnly_Missed_MaxTime->Fill(maxTime);
                            h_tSum_2ClustOnly_Missed_timeDif->Fill(maxTime-minTime);

                            if(maxTime-minTime < 16.0){
                                h_tSum_2ClustOnly_Missed_EdeltaT->Fill(clust_min);
                                h_tSum_2ClustOnly_Missed_nBlocks_less->Fill(clust_min_bl);
                            }
                            else{
                                h_tSum_2ClustOnly_Missed_nBlocks_greater->Fill(clust_min_bl);
                            }
                        }
                        h_tSum_2ClustOnly_All->Fill(clustSum);
                        h_tSum_2ClustOnly_All_MinE->Fill(clust_min);
                        h_tSum_2ClustOnly_timeDif->Fill(maxTime - minTime);
                    }

                    if(!trig.GetTriggersBit(CL2_FLAG-nSSPBits)->empty()){
                        h_tSum_2Clust_VTP_800->Fill(clustSum);
                        h_tSum_2Clust_VTP_800_MinE->Fill(clust_min);
                    }
                    else{
                        h_tSum_2Clust_Missed_MaxTime->Fill(maxTime);
                        h_tSum_2Clust_Missed_timeDif->Fill(maxTime-minTime);

                        if(maxTime-minTime < 16.0){
                            h_tSum_2Clust_Missed_EdeltaT->Fill(clust_min);
                            h_tSum_2Clust_Missed_nBlocks_less->Fill(clust_min_bl);
                        }
                        else{
                            h_tSum_2Clust_Missed_nBlocks_greater->Fill(clust_min_bl);
                        }
                    }
                    h_tSum_2Clust_All->Fill(clustSum);
                    h_tSum_2Clust_All_MinE->Fill(clust_min);
                    h_tSum_2Clust_timeDif->Fill(maxTime - minTime);

                    if(nClust >= 3){

                        if(nClust == 3){
                            if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                                h_tSum_3ClustOnly_VTP_800->Fill(clustSum);
                                h_tSum_3ClustOnly_VTP_800_MinE->Fill(clust_min);
                            }
                            else{
                                h_tSum_3ClustOnly_Missed_MaxTime->Fill(maxTime);
                                h_tSum_3ClustOnly_Missed_timeDif->Fill(maxTime-minTime);

                                if(maxTime-minTime < 16.0){
                                    h_tSum_3ClustOnly_Missed_EdeltaT->Fill(clust_min);
                                    h_tSum_3ClustOnly_Missed_nBlocks_less->Fill(clust_min_bl);
                                }
                                else{
                                    h_tSum_3ClustOnly_Missed_nBlocks_greater->Fill(clust_min_bl);
                                }
                            }
                            h_tSum_3ClustOnly_All->Fill(clustSum);
                            h_tSum_3ClustOnly_All_MinE->Fill(clust_min);
                            h_tSum_3ClustOnly_timeDif->Fill(maxTime - minTime);
                        }

                        if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                            h_tSum_3Clust_VTP_800->Fill(clustSum);
                            h_tSum_3Clust_VTP_800_MinE->Fill(clust_min);
                        }
                        else{
                            h_tSum_3Clust_Missed_MaxTime->Fill(maxTime);
                            h_tSum_3Clust_Missed_timeDif->Fill(maxTime-minTime);

                            if(maxTime-minTime < 16.0){
                                h_tSum_3Clust_Missed_EdeltaT->Fill(clust_min);
                                h_tSum_3Clust_Missed_nBlocks_less->Fill(clust_min_bl);
                            }
                            else{
                                h_tSum_3Clust_Missed_nBlocks_greater->Fill(clust_min_bl);
                            }
                        }
                        h_tSum_3Clust_All->Fill(clustSum);
                        h_tSum_3Clust_All_MinE->Fill(clust_min);
                        h_tSum_3Clust_timeDif->Fill(maxTime - minTime);
                    }
               }
            }


        }
    }
}

/**
 * Processes the events one by one depending on which of the flags are set to be true.
 * For the LMS it looks for efficiencies when other trigger are fired and also evaluates
 * if any channels were dropped.
 *
 * @param self - Whether or not this trigger should have evaluation done within its own events. {True for Yes; False for No}
 * @param rand - Whether or not this trigger should have evalutaion done using a random trigger source. {True for Yes; False for No}
 * @param tSum - Whether or not this trigger should have evaluation done using total sum trigger source. {True for Yes; False for No}
 */
void ClustTrig::ProcessData_OfflineWithThr(bool self, bool rand, bool tSum){
    
    //Get the amount of entries from each file to limit looping through them.
	Long64_t entries = chain->GetEntries();
    
    for(Int_t i = 0; i < entries; i++){
        chain->GetEntry(i);

        trig.Parse(*sspRawPtr);

        // Fill per-bit TRIG_TIME histograms - from Rafo's code
        for (int b = 0; b < nSSPBits; ++b){
            const vector<int>* times = trig.GetTriggersBit(b);
            if (!times) continue;
            for (int t : *times){
                hTrigTime[b]->Fill(t);
            }
        }


        if(self && (trigger_bits & (1<<CL1_FLAG))){

        }
        if(rand && (trigger_bits & (1<<RAND_FC_FLAG))){

        }
        //cout<<trigger_bits<<endl;
        if(tSum && (!trig.GetTriggersBit(TSUM_FLAG-nSSPBits)->empty()) && !(trigger_bits & (1<<LMS_FLAG))){//(trigType & (1<<SSP_RAWSUM_TFLAG))

            //cout<<trigger_bits<<endl;

            Double_t clustSum = 0;
            Int_t clustAbove = 0;
            Double_t maxTime = -100000;
            Double_t minTime = 10000000;

            Double_t clust_min = 1000000;
            UChar_t clust_min_bl = 244;
            for(Int_t j = 0; j < nClust; j++){
                if((i+1)%10000 == 0 || (entries-i+1)<10000){
                    cout<<"\rClustering Trigger Total Sum Events: " << i+1 << "/" << entries << flush;
                    if(i+1 == entries){
                        cout<<endl;
                    }
                }
                if(cl_E[j] > CL_IND_THR && cl_nblocks[j] > 1){
                    clustSum += cl_E[j];
                    clustAbove++;
                    if(cl_E[j] < clust_min){
                        clust_min = cl_E[j];
                        clust_min_bl = cl_nblocks[j];
                    }
                    if(cl_time[j] > maxTime){
                        maxTime = cl_time[j];
                    }
                    if(cl_time[j] < minTime){
                        minTime = cl_time[j];
                    }
                }

                Double_t theta = TMath::ATan2(TMath::Sqrt(cl_x[j]*cl_x[j]+cl_y[j]*cl_y[j]),cl_z[j])*rad2Deg;

                if(nClust >= 2 && cl_E[j] > CL_IND_THR){
                    if(!trig.GetTriggersBit(CL2_FLAG-nSSPBits)->empty()){
                        h_tSum_2Clust_VTP_HC_XY->Fill(cl_x[j], cl_y[j]);
                    }
                    else{
                        h_tSum_2Clust_Miss_HC_XY->Fill(cl_x[j], cl_y[j]);
                        h_tSum_2Clust_Miss_HC_Theta->Fill(theta);
                    }
                    h_tSum_2Clust_All_HC_XY->Fill(cl_x[j], cl_y[j]);
                }

                if(nClust >= 3 && cl_E[j] > CL_IND_THR){
                    if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                        h_tSum_3Clust_VTP_HC_XY->Fill(cl_x[j], cl_y[j]);
                    }
                    else{
                        h_tSum_3Clust_Miss_HC_XY->Fill(cl_x[j], cl_y[j]);
                        h_tSum_3Clust_Miss_HC_Theta->Fill(theta);
                    }
                    h_tSum_3Clust_All_HC_XY->Fill(cl_x[j], cl_y[j]);
                }
            }

            h_tSum_VTP_NClustDiff->Fill(vtp_cl_n - clustAbove);

            if(nClust >= 1 && clustAbove >= 1){

                if(nClust == 1 && clustAbove == 1){
                    if(!trig.GetTriggersBit(CL1_FLAG-nSSPBits)->empty()){//trigger_bits & (1<<CL1_FLAG)){
                        h_tSum_1ClustOnly_VTP_800->Fill(clustSum);
                    }
                    if(!trig.GetTriggersBit(CL1_500_FLAG-nSSPBits)->empty()){
                        h_tSum_1ClustOnly_VTP_500->Fill(clustSum);
                    }
                    if(!trig.GetTriggersBit(CL1_600_FLAG-nSSPBits)->empty()){
                        h_tSum_1ClustOnly_VTP_600->Fill(clustSum);
                    }
                    if(!trig.GetTriggersBit(CL1_900_FLAG-nSSPBits)->empty()){
                        h_tSum_1ClustOnly_VTP_900->Fill(clustSum);
                    }
                    h_tSum_1ClustOnly_All->Fill(clustSum);
                }

               if(!trig.GetTriggersBit(CL1_FLAG-nSSPBits)->empty()){
                    h_tSum_1Clust_VTP_800->Fill(clustSum);
               }
               if(!trig.GetTriggersBit(CL1_500_FLAG-nSSPBits)->empty()){
                    h_tSum_1Clust_VTP_500->Fill(clustSum);
               }
               if(!trig.GetTriggersBit(CL1_600_FLAG-nSSPBits)->empty()){
                    h_tSum_1Clust_VTP_600->Fill(clustSum);
               }
               if(!trig.GetTriggersBit(CL1_900_FLAG-nSSPBits)->empty()){
                    h_tSum_1Clust_VTP_900->Fill(clustSum);
               }
               h_tSum_1Clust_All->Fill(clustSum);
               
               if(nClust >= 2 && clustAbove >= 2 && maxTime - minTime < 16){

                    if(nClust == 2 && clustAbove == 2){
                        if(!trig.GetTriggersBit(CL2_FLAG-nSSPBits)->empty()){
                            h_tSum_2ClustOnly_VTP_800->Fill(clustSum);
                            h_tSum_2ClustOnly_VTP_800_MinE->Fill(clust_min);
                        }
                        else{
                            h_tSum_2ClustOnly_Missed_MaxTime->Fill(maxTime);
                            h_tSum_2ClustOnly_Missed_timeDif->Fill(maxTime-minTime);

                            if(maxTime-minTime < 16.0){
                                h_tSum_2ClustOnly_Missed_EdeltaT->Fill(clust_min);
                                h_tSum_2ClustOnly_Missed_nBlocks_less->Fill(clust_min_bl);
                            }
                            else{
                                h_tSum_2ClustOnly_Missed_nBlocks_greater->Fill(clust_min_bl);
                            }
                        }
                        h_tSum_2ClustOnly_All->Fill(clustSum);
                        h_tSum_2ClustOnly_All_MinE->Fill(clust_min);
                        h_tSum_2ClustOnly_timeDif->Fill(maxTime - minTime);
                    }

                    if(!trig.GetTriggersBit(CL2_FLAG-nSSPBits)->empty()){
                        h_tSum_2Clust_VTP_800->Fill(clustSum);
                        h_tSum_2Clust_VTP_800_MinE->Fill(clust_min);
                    }
                    else{
                        h_tSum_2Clust_Missed_MaxTime->Fill(maxTime);
                        h_tSum_2Clust_Missed_timeDif->Fill(maxTime-minTime);

                        if(maxTime-minTime < 16.0){
                            h_tSum_2Clust_Missed_EdeltaT->Fill(clust_min);
                            h_tSum_2Clust_Missed_nBlocks_less->Fill(clust_min_bl);
                        }
                        else{
                            h_tSum_2Clust_Missed_nBlocks_greater->Fill(clust_min_bl);
                        }
                    }
                    h_tSum_2Clust_All->Fill(clustSum);
                    h_tSum_2Clust_All_MinE->Fill(clust_min);
                    h_tSum_2Clust_timeDif->Fill(maxTime - minTime);

                    if(nClust >= 3 && clustAbove >= 3){

                        if(nClust == 3 && clustAbove == 3){
                            if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                                h_tSum_3ClustOnly_VTP_800->Fill(clustSum);
                                h_tSum_3ClustOnly_VTP_800_MinE->Fill(clust_min);
                            }
                            else{
                                h_tSum_3ClustOnly_Missed_MaxTime->Fill(maxTime);
                                h_tSum_3ClustOnly_Missed_timeDif->Fill(maxTime-minTime);

                                if(maxTime-minTime < 16.0){
                                    h_tSum_3ClustOnly_Missed_EdeltaT->Fill(clust_min);
                                    h_tSum_3ClustOnly_Missed_nBlocks_less->Fill(clust_min_bl);
                                }
                                else{
                                    h_tSum_3ClustOnly_Missed_nBlocks_greater->Fill(clust_min_bl);
                                }
                            }
                            h_tSum_3ClustOnly_All->Fill(clustSum);
                            h_tSum_3ClustOnly_All_MinE->Fill(clust_min);
                            h_tSum_3ClustOnly_timeDif->Fill(maxTime - minTime);
                        }

                        if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                            h_tSum_3Clust_VTP_800->Fill(clustSum);
                            h_tSum_3Clust_VTP_800_MinE->Fill(clust_min);
                        }
                        else{
                            h_tSum_3Clust_Missed_MaxTime->Fill(maxTime);
                            h_tSum_3Clust_Missed_timeDif->Fill(maxTime-minTime);

                            if(maxTime-minTime < 16.0){
                                h_tSum_3Clust_Missed_EdeltaT->Fill(clust_min);
                                h_tSum_3Clust_Missed_nBlocks_less->Fill(clust_min_bl);
                            }
                            else{
                                h_tSum_3Clust_Missed_nBlocks_greater->Fill(clust_min_bl);
                            }
                        }
                        h_tSum_3Clust_All->Fill(clustSum);
                        h_tSum_3Clust_All_MinE->Fill(clust_min);
                        h_tSum_3Clust_timeDif->Fill(maxTime - minTime);
                    }
               }
            }


        }
    }
}

/**
 * Sets up the histograms for this object.
 */
void ClustTrig::setup_Histos(){
    //h_tSum_1Clust_Offline = new TH1F("h_tSum1Clust_Offline", "Total Sum Event with at least 1 Cluster Found Offline Energy;E_{Sum};Counts",en+400,0, en+400)
    //h_tSum_1Clust_VTP = new TH1F("h_tSum1Clust_VTP", "Total Sum Event with at least 1 Cluster Found on the VTP Energy;E_{Sum};Counts",en+400,0, en+400)

    //Offline 1 Cluster events (max energy cluster).
    h_tSum_1Clust_All = new TH1F("h_tSum_1Clust_All", "Total Sum Event with at least 1 Cluster Found;E_{Sum};Counts",en+500,0, en+500);
    //VTP 1 Cluster Trigger (max energy cluster).
    h_tSum_1Clust_VTP_800 = new TH1F("h_tSum_1Clust_VTP_800", "Total Sum Event with at least 1 Cluster Found on the VTP with 800 MeV threshold;E_{Sum};Counts",en+500,0, en+500);
    h_tSum_1Clust_VTP_900 = new TH1F("h_tSum_1Clust_VTP_900", "Total Sum Event with at least 1 Cluster Found on the VTP with 900 MeV threshold;E_{Sum};Counts",en+500,0, en+500);
    h_tSum_1Clust_VTP_600 = new TH1F("h_tSum_1Clust_VTP_600", "Total Sum Event with at least 1 Cluster Found on the VTP with 600 MeV threshold;E_{Sum};Counts",en+500,0, en+500);
    h_tSum_1Clust_VTP_500 = new TH1F("h_tSum_1Clust_VTP_500", "Total Sum Event with at least 1 Cluster Found on the VTP with 500 MeV threshold;E_{Sum};Counts",en+500,0, en+500);

    //Offline 1 Cluster Trigger for events where only 1 cluster was found.
    h_tSum_1ClustOnly_All = new TH1F("h_tSum_1ClustOnly_All", "Total Sum Event with ONLY 1 Cluster Found;E_{Sum};Counts",en+500,0, en+500);
    //VTP 1 Cluster Trigger for events where only 1 cluster was found.
    h_tSum_1ClustOnly_VTP_800 = new TH1F("h_tSum_1ClustOnly_VTP_800", "Total Sum Event with ONLY 1 Cluster Found on the VTP with 800 MeV threshold;E_{Sum};Counts",en+500,0, en+500);
    h_tSum_1ClustOnly_VTP_900 = new TH1F("h_tSum_1ClustOnly_VTP_900", "Total Sum Event with at ONLY 1 Cluster Found on the VTP with 900 MeV threshold;E_{Sum};Counts",en+500,0, en+500);
    h_tSum_1ClustOnly_VTP_600 = new TH1F("h_tSum_1ClustOnly_VTP_600", "Total Sum Event with at ONLY 1 Cluster Found on the VTP with 600 MeV threshold;E_{Sum};Counts",en+500,0, en+500);
    h_tSum_1ClustOnly_VTP_500 = new TH1F("h_tSum_1ClustOnly_VTP_500", "Total Sum Event with at ONLY 1 Cluster Found on the VTP with 500 MeV threshold;E_{Sum};Counts",en+500,0, en+500);

    //Offline 2 Cluster events (max energy combination).
    h_tSum_2Clust_All = new TH1F("h_tSum_2Clust_All", "Total Sum Event with at least 2 Clusters Found;E_{Sum};Counts",en+500,0, en+500);
    //VTP 2 Cluster events (max energy combination).
    h_tSum_2Clust_VTP_800 = new TH1F("h_tSum_2Clust_VTP_800", "Total Sum Event with at least 2 Clusters Found on the VTP with 800 MeV sum threshold;E_{Sum};Counts",en+500,0, en+500);
    //Offline 2 Cluster events (minimum cluster energy).
    h_tSum_2Clust_All_MinE = new TH1F("h_tSum_2Clust_All_MinE", "Total Sum Event with at least 2 Clusters Found Minimum Offline Cluster E;E_{Min};Counts",2000,0, 2000);
    //VTP 2 Cluster events (minimum cluster energy).
    h_tSum_2Clust_VTP_800_MinE = new TH1F("h_tSum_2Clust_VTP_800_MinE", "Total Sum Event with at least 2 Clusters Found on the VTP with 800 MeV sum threshold Minimum Offline Cluster E;E_{Min};Counts",2000,0, 2000);
    //Maximum Time difference of 2 Cluster events for events where at least 2 clusters were found.
    h_tSum_2Clust_timeDif = new TH1F("h_tSum_2Clust_timeDif", "Max(cl_time) - Min(cl_time) of seed channels for offline clusters (2 or more);Time [ns]; Counts", 51, -0.5, 50.5);
    //Maximum Time difference of 2 Cluster events for events where at least 2 clusters were found.
    h_tSum_2Clust_Missed_timeDif = new TH1F("h_tSum_2Clust_Missed_timeDif", "Max(cl_time) - Min(cl_time) of seed channels for Missed offline clusters (2 or more);Time [ns]; Counts", 51, -0.5, 50.5);
    //Maximum Cluster Time for 2 Cluster events for events where at least 2 clusters were found and the clustering trigger did not find it.
    h_tSum_2Clust_Missed_MaxTime = new TH1F("h_tSum_2Clust_Missed_MaxTime", "Max(cl_time) of seed channel for MISSED offline clusters (2 or more);Time [ns]; Counts", 401, -0.5, 400.5);
    //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where at least 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    h_tSum_2Clust_Missed_EdeltaT = new TH1F("h_tSum_2Clust_Missed_EdeltaT", "Minimum Cluster Energy of MISSED offline clusters (2 or more);Energy [MeV]; Counts", 500, 0, 500);
    //Number of Blocks for Missed Events with Delta t < 16 ns where at least 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    h_tSum_2Clust_Missed_nBlocks_less = new TH1I("h_tSum_2Clust_Missed_nBlocks_less", "Number of Blocks fired of MISSED offline clusters < 16ns (2 or more);nBlocks; Counts", 11, -0.5, 10.5);
    //Number of Blocks for Missed Events with Delta t > 16 ns where at least 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    h_tSum_2Clust_Missed_nBlocks_greater = new TH1I("h_tSum_2Clust_Missed_nBlocks_greater", "Number of Blocks fired of MISSED offline clusters > 16ns (2 or more);nBlocks; Counts", 11, -0.5, 10.5);

    //Offline 2 Cluster events for events where only 2 clusters were found.
    h_tSum_2ClustOnly_All = new TH1F("h_tSum_2ClustOnly_All", "Total Sum Event with ONLY 2 Clusters Found, Maximum Offline Combination Energy;E_{Sum};Counts",en+500,0, en+500);
    //VTP 2 Cluster events for events where only 2 clusters were found.
    h_tSum_2ClustOnly_VTP_800 = new TH1F("h_tSum_2ClustOnly_VTP_800", "Total Sum Event with ONLY 2 Clusters Found on the VTP with 800 MeV sum threshold;E_{Sum};Counts",en+500,0, en+500);
    //Offline 2 Cluster events (minimum cluster energy).
    h_tSum_2ClustOnly_All_MinE = new TH1F("h_tSum_2ClustOnly_All_MinE", "Total Sum Event with ONLY 2 Clusters Found Minimum Offline Cluster E;E_{Min};Counts",2000,0, 2000);
    //VTP 2 Cluster events (minimum cluster energy).
    h_tSum_2ClustOnly_VTP_800_MinE = new TH1F("h_tSum_2ClustOnly_VTP_800_MinE", "Total Sum Event with ONLY 2 Clusters Found on the VTP with 800 MeV sum threshold Minimum Offline Cluster E;E_{Min};Counts",2000,0, 2000);
    //Maximum Time difference of 2 Cluster events for events where only 2 clusters were found.
    h_tSum_2ClustOnly_timeDif = new TH1F("h_tSum_2ClustOnly_timeDif", "Max(cl_time) - Min(cl_time) of seed channels for offline clusters (2 exactly);Time [ns]; Counts", 51, -0.5, 50.5);
    //Maximum Time difference of 2 Cluster events for events where only 2 clusters were found.
    h_tSum_2ClustOnly_Missed_timeDif = new TH1F("h_tSum_2ClustOnly_Missed_timeDif", "Max(cl_time) - Min(cl_time) of seed channels for Missed offline clusters (2 exactly);Time [ns]; Counts", 51, -0.5, 50.5);
    //Maximum Cluster Time for 2 Cluster events for events where only 2 clusters were found and the clustering trigger did not find it.
    h_tSum_2ClustOnly_Missed_MaxTime = new TH1F("h_tSum_2ClustOnly_Missed_MaxTime", "Max(cl_time) of seed channel for MISSED offline clusters (2 exactly);Time [ns]; Counts", 401, -0.5, 400.5);
    //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where only 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    h_tSum_2ClustOnly_Missed_EdeltaT = new TH1F("h_tSum_2ClustOnly_Missed_EdeltaT", "Minimum Cluster Energy of MISSED offline clusters (2 exactly);Energy [MeV]; Counts", 500, 0, 500);
    //Number of Blocks for Missed Events with Delta t < 16 ns where only 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    h_tSum_2ClustOnly_Missed_nBlocks_less = new TH1I("h_tSum_2ClustOnly_Missed_nBlocks_less", "Number of Blocks fired of MISSED offline clusters < 16ns (2 exactly);nBlocks; Counts", 11, -0.5, 10.5);
    //Number of Blocks for Missed Events with Delta t > 16 ns where only 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    h_tSum_2ClustOnly_Missed_nBlocks_greater = new TH1I("h_tSum_2ClustOnly_Missed_nBlocks_greater", "Number of Blocks fired of MISSED offline clusters > 16ns (2 exactly);nBlocks; Counts", 11, -0.5, 10.5);

    //Offline 3 Cluster events (max energy combination).
    h_tSum_3Clust_All = new TH1F("h_tSum_3Clust_All", "Total Sum Event with at least 3 Clusters Found;E_{Sum};Counts",en+500,0, en+500);
    //VTP 3 Cluster events (max energy combination).
    h_tSum_3Clust_VTP_800 = new TH1F("h_tSum_3Clust_VTP_800", "Total Sum Event with at least 3 Clusters Found on the VTP with 800 MeV sum threshold;E_{Sum};Counts",en+500,0, en+500);
    //Offline 3 Cluster events (minimum cluster energy).
    h_tSum_3Clust_All_MinE = new TH1F("h_tSum_3Clust_All_MinE", "Total Sum Event with at least 3 Clusters Found Minimum Offline Cluster E;E_{Min};Counts",2000,0, 2000);
    //VTP 3 Cluster events (minimum cluster energy).
    h_tSum_3Clust_VTP_800_MinE = new TH1F("h_tSum_3Clust_VTP_800_MinE", "Total Sum Event with at least 3 Clusters Found on the VTP with 800 MeV sum threshold Minimum Offline Cluster E;E_{Min};Counts",2000,0, 2000);
    //Maximum Time difference of 3 Cluster events for events where at least 3 clusters were found.
    h_tSum_3Clust_timeDif = new TH1F("h_tSum_3Clust_timeDif", "Max(cl_time) - Min(cl_time) of seed channels for offline clusters (3 or more);Time [ns]; Counts", 51, -0.5, 50.5);
    //Maximum Time difference of 3 Cluster events for events where at least 3 clusters were found.
    h_tSum_3Clust_Missed_timeDif = new TH1F("h_tSum_3Clust_Missed_timeDif", "Max(cl_time) - Min(cl_time) of seed channels for Missed offline clusters (3 or more);Time [ns]; Counts", 51, -0.5, 50.5);
    //Maximum Cluster Time for 3 Cluster events for events where at least 3 clusters were found and the clustering trigger did not find it.
    h_tSum_3Clust_Missed_MaxTime = new TH1F("h_tSum_3Clust_Missed_MaxTime", "Max(cl_time) of seed channel for MISSED offline clusters (3 or more);Time [ns]; Counts", 401, -0.5, 400.5);
    //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where at least 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    h_tSum_3Clust_Missed_EdeltaT = new TH1F("h_tSum_3Clust_Missed_EdeltaT", "Minimum Cluster Energy of MISSED offline clusters (3 or more);Energy [MeV]; Counts", 500, 0, 500);
    //Number of Blocks for Missed Events with Delta t < 16 ns where at least 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    h_tSum_3Clust_Missed_nBlocks_less = new TH1I("h_tSum_3Clust_Missed_nBlocks_less", "Number of Blocks fired of MISSED offline clusters < 16ns (3 or more);nBlocks; Counts", 11, -0.5, 10.5);
    //Number of Blocks for Missed Events with Delta t > 16 ns where at least 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    h_tSum_3Clust_Missed_nBlocks_greater = new TH1I("h_tSum_3Clust_Missed_nBlocks_greater", "Number of Blocks fired of MISSED offline clusters > 16ns (3 or more);nBlocks; Counts", 11, -0.5, 10.5);

    //Offline 3 Cluster events for events where only 3 clusters were found.
    h_tSum_3ClustOnly_All = new TH1F("h_tSum_3ClustOnly_All", "Total Sum Event with ONLY 3 Clusters Found, Maximum Offline Combination Energy;E_{Sum};Counts",en+500,0, en+500);
    //VTP 3 Cluster events for events where only 3 clusters were found.
    h_tSum_3ClustOnly_VTP_800 = new TH1F("h_tSum_3ClustOnly_VTP_800", "Total Sum Event with ONLY 3 Clusters Found on the VTP with 800 MeV sum threshold;E_{Sum};Counts",en+500,0, en+500);
    //Offline 3 Cluster events (minimum cluster energy).
    h_tSum_3ClustOnly_All_MinE = new TH1F("h_tSum_3ClustOnly_All_MinE", "Total Sum Event with ONLY 3 Clusters Found Minimum Offline Cluster E;E_{Min};Counts",2000,0, 2000);
    //VTP 3 Cluster events (minimum cluster energy).
    h_tSum_3ClustOnly_VTP_800_MinE = new TH1F("h_tSum_3ClustOnly_VTP_800_MinE", "Total Sum Event with ONLY 3 Clusters Found on the VTP with 800 MeV sum threshold Minimum Offline Cluster E;E_{Min};Counts",2000,0, 2000);
    //Maximum Time difference of 2 Cluster events for events where only 3 clusters were found.
    h_tSum_3ClustOnly_timeDif = new TH1F("h_tSum_3ClustOnly_timeDif", "Max(cl_time) - Min(cl_time) of seed channels for offline clusters (3 exactly);Time [ns]; Counts", 51, -0.5, 50.5);
    //Maximum Time difference of 2 Cluster events for events where only 3 clusters were found.
    h_tSum_3ClustOnly_Missed_timeDif = new TH1F("h_tSum_3ClustOnly_Missed_timeDif", "Max(cl_time) - Min(cl_time) of seed channels for Missed offline clusters (3 exactly);Time [ns]; Counts", 51, -0.5, 50.5);
    //Maximum Cluster Time for 3 Cluster events for events where only 3 clusters were found and the clustering trigger did not find it.
    h_tSum_3ClustOnly_Missed_MaxTime = new TH1F("h_tSum_3ClustOnly_Missed_MaxTime", "Max(cl_time) of seed channel for MISSED offline clusters (3 exactly);Time [ns]; Counts", 401, -0.5, 400.5);
    //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where only 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    h_tSum_3ClustOnly_Missed_EdeltaT = new TH1F("h_tSum_3ClustOnly_Missed_EdeltaT", "Minimum Cluster Energy of MISSED offline clusters (3 exactly);Energy [MeV]; Counts", 500, 0, 500);
    //Number of Blocks for Missed Events with Delta t < 16 ns where only 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    h_tSum_3ClustOnly_Missed_nBlocks_less = new TH1I("h_tSum_3ClustOnly_Missed_nBlocks_less", "Number of Blocks fired of MISSED offline clusters < 16ns (3 exactly);nBlocks; Counts", 11, -0.5, 10.5);
    //Number of Blocks for Missed Events with Delta t > 16 ns where only 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    h_tSum_3ClustOnly_Missed_nBlocks_greater = new TH1I("h_tSum_3ClustOnly_Missed_nBlocks_greater", "Number of Blocks fired of MISSED offline clusters > 16ns (3 exactly);nBlocks; Counts", 11, -0.5, 10.5);

    //Offline 2 Cluster events; HyCal XY Positions
    h_tSum_2Clust_All_HC_XY = new TH2F("h_tSum_2Clust_All_HC_XY", "Total Sum Events with 2 Offline Clusters Found;x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);
    //VTP 2 Cluster events; HyCal XY Positions
    h_tSum_2Clust_VTP_HC_XY = new TH2F("h_tSum_2Clust_VTP_HC_XY", "Total Sum Events with 2 VTP Clusters Found;x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);
    //2 Cluster events missed by the VTP Trigger; HyCal XY, and Theta Positions
    h_tSum_2Clust_Miss_HC_XY = new TH2F("h_tSum_2Clust_Miss_HC_XY", "Total Sum Events with 2 Clusters Found Missed by the VTP;x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);
    h_tSum_2Clust_Miss_HC_Theta = new TH1F("h_tSum_2Clust_Miss_HC_Theta", "Total Sum Events with 2 Clusters Found Missed by the VTP;#theta [#circ];Counts",50, 0, 5);

    //Offline 2 Cluster events; HyCal XY Positions
    h_tSum_3Clust_All_HC_XY = new TH2F("h_tSum_3Clust_All_HC_XY", "Total Sum Events with 3 Offline Clusters Found;x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);
    //VTP 3 Cluster events; HyCal XY Positions
    h_tSum_3Clust_VTP_HC_XY = new TH2F("h_tSum_3Clust_VTP_HC_XY", "Total Sum Events with 3 VTP Clusters Found;x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);
    //3 Cluster events missed by the VTP Trigger; HyCal XY, and Theta Positions
    h_tSum_3Clust_Miss_HC_XY = new TH2F("h_tSum_3Clust_Miss_HC_XY", "Total Sum Events with 3 Clusters Found Missed by the VTP;x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);
    h_tSum_3Clust_Miss_HC_Theta = new TH1F("h_tSum_3Clust_Miss_HC_Theta", "Total Sum Events with 3 Clusters Found Missed by the VTP;#theta [#circ];Counts",50, 0, 5);

    h_tSum_VTP_NClustDiff = new TH1I("h_tSum_VTP_NClustDiff", "(No. VTP Clusters) - (No. Offline Clusters);Difference;Count",11, -5.5, 5.5);
}

/**
 * Prints out the trigger performance plots from the 
 */
void ClustTrig::printTSumPDF(TString pdfName){

    h_tSum_1Clust_All->SetLineColor(kBlack);
    h_tSum_1Clust_VTP_800->SetLineColor(kBlue);
    h_tSum_1Clust_VTP_500->SetLineColor(kGreen);
    h_tSum_1Clust_VTP_600->SetLineColor(kRed);
    h_tSum_1Clust_VTP_900->SetLineColor(kCyan);

    TH1F* h_tSum_1ClustRatio_900 = (TH1F*) h_tSum_1Clust_VTP_900->Clone();
    h_tSum_1ClustRatio_900->Divide(h_tSum_1Clust_All);
    h_tSum_1ClustRatio_900->SetLineColor(kCyan);
    h_tSum_1ClustRatio_900->SetTitle("Trigger Efficiency;E_{Sum};Efficiency");

    TH1F* h_tSum_1ClustRatio_800 = (TH1F*) h_tSum_1Clust_VTP_800->Clone();
    h_tSum_1ClustRatio_800->Divide(h_tSum_1Clust_All);
    h_tSum_1ClustRatio_800->SetLineColor(kBlue);
    h_tSum_1ClustRatio_800->SetTitle("Trigger Efficiency;E_{Sum};Efficiency");

    TH1F* h_tSum_1ClustRatio_600 = (TH1F*) h_tSum_1Clust_VTP_600->Clone();
    h_tSum_1ClustRatio_600->Divide(h_tSum_1Clust_All);
    h_tSum_1ClustRatio_600->SetLineColor(kRed);
    h_tSum_1ClustRatio_600->SetTitle("Trigger Efficiency;E_{Sum};Efficiency");

    TH1F* h_tSum_1ClustRatio_500 = (TH1F*) h_tSum_1Clust_VTP_500->Clone();
    h_tSum_1ClustRatio_500->Divide(h_tSum_1Clust_All);
    h_tSum_1ClustRatio_500->SetLineColor(kGreen);
    h_tSum_1ClustRatio_500->SetTitle("Trigger Efficiency;E_{Sum};Efficiency");


    h_tSum_1ClustOnly_All->SetLineColor(kBlack);
    h_tSum_1ClustOnly_VTP_800->SetLineColor(kBlue);
    h_tSum_1ClustOnly_VTP_500->SetLineColor(kGreen);
    h_tSum_1ClustOnly_VTP_600->SetLineColor(kRed);
    h_tSum_1ClustOnly_VTP_900->SetLineColor(kCyan);

    TH1F* h_tSum_1ClustOnlyRatio_900 = (TH1F*) h_tSum_1ClustOnly_VTP_900->Clone();
    h_tSum_1ClustOnlyRatio_900->Divide(h_tSum_1ClustOnly_All);
    h_tSum_1ClustOnlyRatio_900->SetLineColor(kCyan);
    h_tSum_1ClustOnlyRatio_900->SetTitle("Exactly 1 Cluster Trigger Efficiency;E_{Sum};Efficiency");

    TH1F* h_tSum_1ClustOnlyRatio_800 = (TH1F*) h_tSum_1ClustOnly_VTP_800->Clone();
    h_tSum_1ClustOnlyRatio_800->Divide(h_tSum_1ClustOnly_All);
    h_tSum_1ClustOnlyRatio_800->SetLineColor(kBlue);
    h_tSum_1ClustOnlyRatio_800->SetTitle("Exactly 1 Cluster Trigger Efficiency;E_{Sum};Efficiency");

    TH1F* h_tSum_1ClustOnlyRatio_600 = (TH1F*) h_tSum_1ClustOnly_VTP_600->Clone();
    h_tSum_1ClustOnlyRatio_600->Divide(h_tSum_1ClustOnly_All);
    h_tSum_1ClustOnlyRatio_600->SetLineColor(kRed);
    h_tSum_1ClustOnlyRatio_600->SetTitle("Exactly 1 Cluster Trigger Efficiency;E_{Sum};Efficiency");

    TH1F* h_tSum_1ClustOnlyRatio_500 = (TH1F*) h_tSum_1ClustOnly_VTP_500->Clone();
    h_tSum_1ClustOnlyRatio_500->Divide(h_tSum_1ClustOnly_All);
    h_tSum_1ClustOnlyRatio_500->SetLineColor(kGreen);
    h_tSum_1ClustOnlyRatio_500->SetTitle("Exactly 1 Cluster Trigger Efficiency;E_{Sum};Efficiency");

    h_tSum_2Clust_All->SetLineColor(kBlack);
    h_tSum_2Clust_VTP_800->SetLineColor(kBlue);

    TH1F* h_tSum_2ClustRatio = (TH1F*) h_tSum_2Clust_VTP_800->Clone();
    h_tSum_2ClustRatio->Divide(h_tSum_2Clust_All);
    h_tSum_2ClustRatio->SetLineColor(kBlack);
    h_tSum_2ClustRatio->SetTitle("Trigger Efficiency;E_{Sum};Efficiency");

    h_tSum_2Clust_All_MinE->SetLineColor(kBlack);
    h_tSum_2Clust_VTP_800_MinE->SetLineColor(kBlue);

    TH1F* h_tSum_2ClustRatio_MinE = (TH1F*) h_tSum_2Clust_VTP_800_MinE->Clone();
    h_tSum_2ClustRatio_MinE->Divide(h_tSum_2Clust_All_MinE);
    h_tSum_2ClustRatio_MinE->SetLineColor(kBlack);
    h_tSum_2ClustRatio_MinE->SetTitle("Trigger Efficiency;E_{Min};Efficiency");

    h_tSum_2ClustOnly_All->SetLineColor(kBlack);
    h_tSum_2ClustOnly_VTP_800->SetLineColor(kBlue);

    TH1F* h_tSum_2ClustOnlyRatio = (TH1F*) h_tSum_2ClustOnly_VTP_800->Clone();
    h_tSum_2ClustOnlyRatio->Divide(h_tSum_2ClustOnly_All);
    h_tSum_2ClustOnlyRatio->SetLineColor(kBlack);
    h_tSum_2ClustOnlyRatio->SetTitle("Exactly 2 Clusters Trigger Efficiency;E_{Sum};Efficiency");

    h_tSum_2ClustOnly_All_MinE->SetLineColor(kBlack);
    h_tSum_2ClustOnly_VTP_800_MinE->SetLineColor(kBlue);
    
    TH1F* h_tSum_2ClustOnlyRatio_MinE = (TH1F*) h_tSum_2ClustOnly_VTP_800_MinE->Clone();
    h_tSum_2ClustOnlyRatio_MinE->Divide(h_tSum_2ClustOnly_All_MinE);
    h_tSum_2ClustOnlyRatio_MinE->SetLineColor(kBlack);
    h_tSum_2ClustOnlyRatio_MinE->SetTitle("Trigger Efficiency;E_{Min};Efficiency");

    h_tSum_3Clust_All->SetLineColor(kBlack);
    h_tSum_3Clust_VTP_800->SetLineColor(kBlue);

    TH1F* h_tSum_3ClustRatio = (TH1F*) h_tSum_3Clust_VTP_800->Clone();
    h_tSum_3ClustRatio->Divide(h_tSum_3Clust_All);
    h_tSum_3ClustRatio->SetLineColor(kBlack);
    h_tSum_3ClustRatio->SetTitle("Trigger Efficiency;E_{Sum};Efficiency");

    h_tSum_3Clust_All_MinE->SetLineColor(kBlack);
    h_tSum_3Clust_VTP_800_MinE->SetLineColor(kBlue);

    TH1F* h_tSum_3ClustRatio_MinE = (TH1F*) h_tSum_3Clust_VTP_800_MinE->Clone();
    h_tSum_3ClustRatio_MinE->Divide(h_tSum_3Clust_All_MinE);
    h_tSum_3ClustRatio_MinE->SetLineColor(kBlack);
    h_tSum_3ClustRatio_MinE->SetTitle("Trigger Efficiency;E_{Min};Efficiency");

    h_tSum_3ClustOnly_All->SetLineColor(kBlack);
    h_tSum_3ClustOnly_VTP_800->SetLineColor(kBlue);

    TH1F* h_tSum_3ClustOnlyRatio = (TH1F*) h_tSum_3ClustOnly_VTP_800->Clone();
    h_tSum_3ClustOnlyRatio->Divide(h_tSum_3ClustOnly_All);
    h_tSum_3ClustOnlyRatio->SetLineColor(kBlack);
    h_tSum_3ClustOnlyRatio->SetTitle("Exactly 3 Clusters Trigger Efficiency;E_{Sum};Efficiency");

    h_tSum_3ClustOnly_All_MinE->SetLineColor(kBlack);
    h_tSum_3ClustOnly_VTP_800_MinE->SetLineColor(kBlue);
    
    TH1F* h_tSum_3ClustOnlyRatio_MinE = (TH1F*) h_tSum_3ClustOnly_VTP_800_MinE->Clone();
    h_tSum_3ClustOnlyRatio_MinE->Divide(h_tSum_3ClustOnly_All_MinE);
    h_tSum_3ClustOnlyRatio_MinE->SetLineColor(kBlack);
    h_tSum_3ClustOnlyRatio_MinE->SetTitle("Trigger Efficiency;E_{Min};Efficiency");

    TCanvas *c = new TCanvas("c", "TSum_ClustTrig_Canvas",1000,1000);
    auto legend = new TLegend(0.1,0.8,0.4,0.9);
    auto legend2 = new TLegend(0.1,0.8,0.4,0.9);

    //Page 1
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_1Clust_All->Draw("HIST 0");
    h_tSum_1Clust_VTP_800->Draw("HIST SAME");
    h_tSum_1Clust_VTP_500->Draw("HIST SAME");
    h_tSum_1Clust_VTP_600->Draw("HIST SAME");
    h_tSum_1Clust_VTP_900->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_1Clust_All,"All","l");
	legend->AddEntry(h_tSum_1Clust_VTP_900,"900 MeV","l");
	legend->AddEntry(h_tSum_1Clust_VTP_800,"800 MeV","l");
	legend->AddEntry(h_tSum_1Clust_VTP_600,"600 MeV","l");
	legend->AddEntry(h_tSum_1Clust_VTP_500,"500 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_1ClustRatio_900->Draw("HIST");
    h_tSum_1ClustRatio_800->Draw("HIST SAME");
    h_tSum_1ClustRatio_600->Draw("HIST SAME");
    h_tSum_1ClustRatio_500->Draw("HIST SAME");
    
    legend2->SetHeader("Legend","C"); // option "C" allows to center the header
	legend2->AddEntry(h_tSum_1ClustRatio_900,"900 MeV","l");
	legend2->AddEntry(h_tSum_1ClustRatio_800,"800 MeV","l");
	legend2->AddEntry(h_tSum_1ClustRatio_600,"600 MeV","l");
	legend2->AddEntry(h_tSum_1ClustRatio_500,"500 MeV","l");
    legend2->Draw();
    c->Print(pdfName + "(");
    c->Clear();
    legend->Clear();
    legend2->Clear();

    //Page 2
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_1ClustOnly_All->Draw("HIST");
    h_tSum_1ClustOnly_VTP_800->Draw("HIST SAME");
    h_tSum_1ClustOnly_VTP_500->Draw("HIST SAME");
    h_tSum_1ClustOnly_VTP_600->Draw("HIST SAME");
    h_tSum_1ClustOnly_VTP_900->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_1ClustOnly_All,"All","l");
	legend->AddEntry(h_tSum_1ClustOnly_VTP_900,"900 MeV","l");
	legend->AddEntry(h_tSum_1ClustOnly_VTP_800,"800 MeV","l");
	legend->AddEntry(h_tSum_1ClustOnly_VTP_600,"600 MeV","l");
	legend->AddEntry(h_tSum_1ClustOnly_VTP_500,"500 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_1ClustOnlyRatio_900->Draw("HIST");
    h_tSum_1ClustOnlyRatio_800->Draw("HIST SAME");
    h_tSum_1ClustOnlyRatio_600->Draw("HIST SAME");
    h_tSum_1ClustOnlyRatio_500->Draw("HIST SAME");
    
    legend2->SetHeader("Legend","C"); // option "C" allows to center the header
	legend2->AddEntry(h_tSum_1ClustOnlyRatio_900,"900 MeV","l");
	legend2->AddEntry(h_tSum_1ClustOnlyRatio_800,"800 MeV","l");
	legend2->AddEntry(h_tSum_1ClustOnlyRatio_600,"600 MeV","l");
	legend2->AddEntry(h_tSum_1ClustOnlyRatio_500,"500 MeV","l");
    legend2->Draw();
    c->Print(pdfName);
    c->Clear();
    legend->Clear();
    legend2->Clear();

    //Page 3
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2Clust_All->Draw("HIST");
    h_tSum_2Clust_VTP_800->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2Clust_All,"All","l");
	legend->AddEntry(h_tSum_2Clust_VTP_800,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustRatio->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 4
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2Clust_All->Draw("HIST");
    h_tSum_2Clust_VTP_800->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2Clust_All,"All","l");
	legend->AddEntry(h_tSum_2Clust_VTP_800,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustRatio->SetAxisRange(0.9,1.0, "Y");
    h_tSum_2ClustRatio->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 5
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2Clust_All_MinE->Draw("HIST");
    h_tSum_2Clust_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2Clust_All_MinE,"All","l");
	legend->AddEntry(h_tSum_2Clust_VTP_800_MinE,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 6
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2Clust_All_MinE->Draw("HIST");
    h_tSum_2Clust_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2Clust_All_MinE,"All","l");
	legend->AddEntry(h_tSum_2Clust_VTP_800_MinE,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustRatio_MinE->SetAxisRange(0.9, 1.0, "Y");
    h_tSum_2ClustRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 7
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_2Clust_timeDif->Draw("HIST");
    c->cd(2);
    h_tSum_2Clust_Missed_timeDif->Draw("HIST");
    c->cd(3);
    h_tSum_2Clust_Missed_MaxTime->Draw("HIST");
    c->cd(4);
    h_tSum_2Clust_Missed_EdeltaT->Draw("HIST");
    c->Print(pdfName);
    c->Clear();

    //Page 8
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_2Clust_Missed_nBlocks_less->Draw("HIST");
    c->cd(2);
    h_tSum_2Clust_Missed_nBlocks_greater->Draw("HIST");
    c->Print(pdfName);
    c->Clear();


    //Page 9
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2ClustOnly_All->Draw("HIST");
    h_tSum_2ClustOnly_VTP_800->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2ClustOnly_All,"All","l");
	legend->AddEntry(h_tSum_2ClustOnly_VTP_800,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustOnlyRatio->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 10
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2ClustOnly_All->Draw("HIST");
    h_tSum_2ClustOnly_VTP_800->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2ClustOnly_All,"All","l");
	legend->AddEntry(h_tSum_2ClustOnly_VTP_800,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustOnlyRatio->SetAxisRange(0.9, 1.0, "Y");
    h_tSum_2ClustOnlyRatio->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 11
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2ClustOnly_All_MinE->Draw("HIST");
    h_tSum_2ClustOnly_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2ClustOnly_All_MinE,"All","l");
	legend->AddEntry(h_tSum_2ClustOnly_VTP_800_MinE,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustOnlyRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 12
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2ClustOnly_All_MinE->Draw("HIST");
    h_tSum_2ClustOnly_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2ClustOnly_All_MinE,"All","l");
	legend->AddEntry(h_tSum_2ClustOnly_VTP_800_MinE,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustOnlyRatio_MinE->SetAxisRange(0.9, 1.0, "Y");
    h_tSum_2ClustOnlyRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 13
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_2ClustOnly_timeDif->Draw("HIST");
    c->cd(2);
    h_tSum_2ClustOnly_Missed_timeDif->Draw("HIST");
    c->cd(3);
    h_tSum_2ClustOnly_Missed_MaxTime->Draw("HIST");
    c->cd(4);
    h_tSum_2ClustOnly_Missed_EdeltaT->Draw("HIST");
    c->Print(pdfName);
    c->Clear();

    //Page 14
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_2ClustOnly_Missed_nBlocks_less->Draw("HIST");
    c->cd(2);
    h_tSum_2ClustOnly_Missed_nBlocks_greater->Draw("HIST");
    c->Print(pdfName);
    c->Clear();

    //Page 15
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    gPad->SetLogz(1);
    h_tSum_2Clust_All_HC_XY->Draw("COLZ");
    c->cd(2);
    gPad->SetLogz(1);
    h_tSum_2Clust_VTP_HC_XY->Draw("COLZ");
    c->cd(3);
    gPad->SetLogz(1);
    h_tSum_2Clust_Miss_HC_XY->Draw("COLZ");
    c->cd(4);
    gPad->SetLogy(1);
    h_tSum_2Clust_Miss_HC_Theta->Draw("HIST");
    c->Print(pdfName);
    c->Clear();
    gPad->SetLogz(0);
    gPad->SetLogy(0);

    //Page 16
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3Clust_All->Draw("HIST");
    h_tSum_3Clust_VTP_800->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3Clust_All,"All","l");
	legend->AddEntry(h_tSum_3Clust_VTP_800,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_3ClustRatio->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 17
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3Clust_All_MinE->Draw("HIST");
    h_tSum_3Clust_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3Clust_All_MinE,"All","l");
	legend->AddEntry(h_tSum_3Clust_VTP_800_MinE,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_3ClustRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 18
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3Clust_All_MinE->Draw("HIST");
    h_tSum_3Clust_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3Clust_All_MinE,"All","l");
	legend->AddEntry(h_tSum_3Clust_VTP_800_MinE,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_3ClustRatio_MinE->SetAxisRange(0.9,1.0,"Y");
    h_tSum_3ClustRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 19
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_3Clust_timeDif->Draw("HIST");
    c->cd(2);
    h_tSum_3Clust_Missed_timeDif->Draw("HIST");
    c->cd(3);
    h_tSum_3Clust_Missed_MaxTime->Draw("HIST");
    c->cd(4);
    h_tSum_3Clust_Missed_EdeltaT->Draw("HIST");
    c->Print(pdfName);
    c->Clear();

    //Page 20
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_3Clust_Missed_nBlocks_less->Draw("HIST");
    c->cd(2);
    h_tSum_3Clust_Missed_nBlocks_greater->Draw("HIST");
    c->Print(pdfName);
    c->Clear();

    //Page 20
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    gPad->SetLogz(1);
    h_tSum_3Clust_All_HC_XY->Draw("COLZ");
    c->cd(2);
    gPad->SetLogz(1);
    h_tSum_3Clust_VTP_HC_XY->Draw("COLZ");
    c->cd(3);
    gPad->SetLogz(1);
    h_tSum_3Clust_Miss_HC_XY->Draw("COLZ");
    c->cd(4);
    gPad->SetLogy(1);
    h_tSum_3Clust_Miss_HC_Theta->Draw("HIST");
    c->Print(pdfName);
    c->Clear();
    gPad->SetLogy(0);
    gPad->SetLogz(0);

    //Page 21
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3ClustOnly_All->Draw("HIST");
    h_tSum_3ClustOnly_VTP_800->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3ClustOnly_All,"All","l");
	legend->AddEntry(h_tSum_3ClustOnly_VTP_800,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_3ClustOnlyRatio->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 22
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3ClustOnly_All_MinE->Draw("HIST");
    h_tSum_3ClustOnly_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3ClustOnly_All_MinE,"All","l");
	legend->AddEntry(h_tSum_3ClustOnly_VTP_800_MinE,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_3ClustOnlyRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 23
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3ClustOnly_All_MinE->Draw("HIST");
    h_tSum_3ClustOnly_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3ClustOnly_All_MinE,"All","l");
	legend->AddEntry(h_tSum_3ClustOnly_VTP_800_MinE,"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_3ClustOnlyRatio_MinE->SetAxisRange(0.9, 1.0, "Y");
    h_tSum_3ClustOnlyRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 24
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_3ClustOnly_timeDif->Draw("HIST");
    c->cd(2);
    h_tSum_3ClustOnly_Missed_timeDif->Draw("HIST");
    c->cd(3);
    h_tSum_3ClustOnly_Missed_MaxTime->Draw("HIST");
    c->cd(4);
    h_tSum_3ClustOnly_Missed_EdeltaT->Draw("HIST");
    c->Print(pdfName);
    c->Clear();

    //Page 25
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_3ClustOnly_Missed_nBlocks_less->Draw("HIST");
    c->cd(2);
    h_tSum_3ClustOnly_Missed_nBlocks_greater->Draw("HIST");
    c->Print(pdfName);
    c->Clear();

    //Page 26
    //-----------------------------------------------------------------------------------
    c->Divide(2,4);
    for(Int_t b = 0; b < nSSPBits; b++){
        c->cd(b+1);
        hTrigTime[b]->Draw("HIST");
        hTrigTime[b]->SetStats(0);
    }
    c->Print(pdfName);
    c->Clear();

    c->Print(pdfName);

    c->cd(1);
    h_tSum_VTP_NClustDiff->Draw("HIST");
    c->Print(pdfName + ")");

}