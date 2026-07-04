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
ClustTrig::ClustTrig(TChain* c):Trigger(c, true, true){
    cout<<"Set All Branch Addresses for Clustering Trigger Object.\n";

    en = 2239.51;//3485.41;

    vector<Int_t> excludeMod = {1526, 1527, 1528, 1529, 1560, 1563, 1594, 1597, 1628, 1629, 1630, 1631};

    //Set the bounds of the Histograms we are going to visualize.
   Xedge[0] = -17.*20.77; // mm
   Yedge[0] = -17.*20.75; // mm
   for(int i=1; i<35; i++){
      Xedge[i] = Xedge[i-1] + 20.77;
      Yedge[i] = Yedge[i-1] + 20.75;
   }
   usemax = true;

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
    
    en = 2239.51;//3485.41;

    usemax = true;

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

    if(tSum){
        setup_tSum_Histos();
    }
    if(rand){
        setup_rand_Histos();
    }
    
    //Get the amount of entries from each file to limit looping through them.
	entries = chain->GetEntries();
    
    for(Int_t i = 0; i < entries; i++){
        chain->GetEntry(i);

        if((i+1)%10000 == 0 || (entries-i+1)<10000){
            cout<<"\rClustering Trigger Events: " << i+1 << "/" << entries << flush;
            if(i+1 == entries){
                cout<<endl;
            }
        }

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
        if(tSum && (!trig.GetTriggersBit(TSUM_FLAG-nSSPBits)->empty()) && !(trigger_bits & (1<<LMS_FLAG))){ //(trigType & (1<<SSP_RAWSUM_TFLAG))
            fill_tSum_Histos(i);
        }
    }
}

/**
 * Fill the Total Sum event histograms for each event.
 *
 * @param i - the current event number
 */
void ClustTrig::fill_tSum_Histos(Int_t i){
    Double_t clustSum = 0;
    Double_t maxTime = -100000;
    Double_t minTime = 10000000;
    Double_t clust_min = 1000000;
    UChar_t clust_min_bl = 244;
    Int_t min_ind = -1;
    vector<Int_t> GEM_matches;
    vector<Int_t> GEM_matches_either;
    Int_t GEM_matches_count = 0;
    Int_t GEM_matches_count_both = 0;
    Utils::Point p[MAX_CLUSTERS];
    Double_t clust_max = -1;
    for(Int_t j = 0; j < nClust; j++){
        
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
            min_ind = j;
        }
        if(cl_E[j] > clust_max){
            clust_max = cl_E[j];
        }
        h_tSum_Off_centerID->Fill(cl_center[j]);

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

            if(((match_flag[j] & (1<<0)) || (match_flag[j] & (1<<1))) || ((match_flag[j] & (1<<2)) || (match_flag[j] & (1<<3)))){
                GEM_matches_count++;
                GEM_matches_either.push_back(j);
                if(((match_flag[j] & (1<<0)) || (match_flag[j] & (1<<1))) && ((match_flag[j] & (1<<2)) || (match_flag[j] & (1<<3)))){
                    GEM_matches_count_both++;
                    GEM_matches.push_back(j);

                    vector<Float_t> x_vec;
                    vector<Float_t> y_vec;
                    vector<Float_t> z_vec;
                    /*x_vec.push_back(cl_x[j]);
                    y_vec.push_back(cl_y[j]);
                    z_vec.push_back(cl_z[j]);*/ //Currently the HyCal position still needs correction so they are not being used in the line of best fit yet.
                    for(Int_t q = 0; q < MAX_GEMS; q++){
                        if(mgz[j][q]>0){
                            x_vec.push_back(mgx[j][q]);
                            y_vec.push_back(mgy[j][q]);
                            z_vec.push_back(mgz[j][q]);
                        }
                    }

                    Utils::LineOfBestFit line = Utils::FitLine(x_vec, y_vec, z_vec);
                    p[j] = Utils::ClosestApproachToZAxis(line);
                }
            }
        }
    }

    h_tSum_VTP_NClustDiff->Fill(vtp_cl_n - nClust);
    bool check = false;

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
                h_tSum_VTP_NClustDiff_Missed->Fill(vtp_cl_n - nClust);
                h_tSum_2Clust_Miss_HC_XY_LowE->Fill(cl_x[min_ind],cl_y[min_ind]);

                check = true;

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

                for(Int_t b = 0; b < nClust; b++){
                    if(cl_E[b]>CL_IND_THR){
                        h_tSum_3Clust_Cl_Energy_All->Fill(cl_E[b]);
                        if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                            h_tSum_3Clust_Cl_Energy_Trig->Fill(cl_E[b]);
                        }
                        if(nClust == 3){
                            h_tSum_3ClustOnly_Cl_Energy_All->Fill(cl_E[b]);
                            if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                                h_tSum_3ClustOnly_Cl_Energy_Trig->Fill(cl_E[b]);
                            }
                        }
                    }
                }

                for(UInt_t c = 0; c < GEM_matches.size(); c++){
                    if(GEM_matches_count_both >= 2){
                        h_tSum_3Clust_GEMMatch_All_Cl_E->Fill(cl_E[GEM_matches.at(c)]);
                        if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                            h_tSum_3Clust_GEMMatch_Trig_Cl_E->Fill(cl_E[GEM_matches.at(c)]);
                        }

                        if(GEM_matches_count_both >= 3){
                            h_tSum_3Clust_GEMMatch_All_3Match_Both_Cl_E->Fill(cl_E[GEM_matches.at(c)]);
                            if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                                h_tSum_3Clust_GEMMatch_Trig_3Match_Both_Cl_E->Fill(cl_E[GEM_matches.at(c)]);
                            }
                        }
                    }

                    Float_t gemx_1 = TMath::Max(mgx[GEM_matches.at(c)][0], mgx[GEM_matches.at(c)][2]);
                    Float_t gemx_2 = TMath::Max(mgx[GEM_matches.at(c)][1], mgx[GEM_matches.at(c)][3]);

                    if(clustSum > 2000.0 && clustSum < 2400.0){
                        h_tSum_3Clust_vZ_High->Fill(p[GEM_matches.at(c)].z);
                        h_tSum_3Clust_DCA_High->Fill(TMath::Sqrt(p[GEM_matches.at(c)].x*p[GEM_matches.at(c)].x + p[GEM_matches.at(c)].y*p[GEM_matches.at(c)].y));
                        h_tSum_3Clust_delta_X_High->Fill(gemx_1-gemx_2);
                    }
                    else if(clustSum < 2000.0){
                        h_tSum_3Clust_vZ_Low->Fill(p[GEM_matches.at(c)].z);
                        h_tSum_3Clust_DCA_Low->Fill(TMath::Sqrt(p[GEM_matches.at(c)].x*p[GEM_matches.at(c)].x + p[GEM_matches.at(c)].y*p[GEM_matches.at(c)].y));
                        h_tSum_3Clust_delta_X_Low->Fill(gemx_1-gemx_2);
                    }    
                }

                for(UInt_t d = 0; d < GEM_matches_either.size(); d++){
                    if(GEM_matches_count >= 2){
                        h_tSum_3Clust_GEMMatch_All_2Match_Either_Cl_E->Fill(cl_E[GEM_matches_either.at(d)]);
                        if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                            h_tSum_3Clust_GEMMatch_Trig_2Match_Either_Cl_E->Fill(cl_E[GEM_matches_either.at(d)]);
                        }

                        if(GEM_matches_count >= 3){
                            h_tSum_3Clust_GEMMatch_All_3Match_Either_Cl_E->Fill(cl_E[GEM_matches_either.at(d)]);
                            if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                                h_tSum_3Clust_GEMMatch_Trig_3Match_Either_Cl_E->Fill(cl_E[GEM_matches_either.at(d)]);
                            }
                        }
                    }
                }
                        

                if(GEM_matches_count_both >= 2){
                    h_tSum_3Clust_GEMMatch_All->Fill(clustSum);
                    if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                        h_tSum_3Clust_GEMMatch_Trig->Fill(clustSum);
                    }

                    if(GEM_matches_count_both >= 3){
                        h_tSum_3Clust_GEMMatch_All_3Match_Both->Fill(clustSum);
                        if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                            h_tSum_3Clust_GEMMatch_Trig_3Match_Both->Fill(clustSum);
                        }
                    }
                }

                if(GEM_matches_count >= 2){
                    h_tSum_3Clust_GEMMatch_All_2Match_Either->Fill(clustSum);
                    if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                        h_tSum_3Clust_GEMMatch_Trig_2Match_Either->Fill(clustSum);
                    }

                    if(GEM_matches_count >= 3){
                        h_tSum_3Clust_GEMMatch_All_3Match_Either->Fill(clustSum);
                        if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                            h_tSum_3Clust_GEMMatch_Trig_3Match_Either->Fill(clustSum);
                        }
                    }
                }

                if(nClust == 3){
                    if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                        h_tSum_3ClustOnly_VTP_800->Fill(clustSum);
                        h_tSum_3ClustOnly_VTP_800_MinE->Fill(clust_min);
                        h_tSum_3ClustOnly_Cl_Max_Trig->Fill(clust_max);
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
                    h_tSum_3ClustOnly_Cl_Max_All->Fill(clust_max);
                }

                if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                    h_tSum_3Clust_VTP_800->Fill(clustSum);
                    h_tSum_3Clust_VTP_800_MinE->Fill(clust_min);
                    h_tSum_3Clust_Cl_Max_Trig->Fill(clust_max);
                }
                else{
                    h_tSum_3Clust_Missed_MaxTime->Fill(maxTime);
                    h_tSum_3Clust_Missed_timeDif->Fill(maxTime-minTime);
                    h_tSum_3Clust_Miss_HC_XY_LowE->Fill(cl_x[min_ind],cl_y[min_ind]);

                    if(!check){
                        h_tSum_VTP_NClustDiff_Missed->Fill(vtp_cl_n - nClust);
                    }

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
                h_tSum_3Clust_Cl_Max_All->Fill(clust_max);
            }
        }
    }

    for(Int_t a = 0; a < vtp_cl_n; a++){
        h_tSum_VTP_Timing->Fill(vtp_cl_time[a]);
        h_tSum_VTP_centerID->Fill(vtp_cl_center[a]);
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

    if(tSum){
        setup_tSum_Histos();
    }
    if(rand){
        setup_rand_Histos();
    }
    
    //Get the amount of entries from each file to limit looping through them.
	Long64_t entries = chain->GetEntries();
    
    for(Int_t i = 0; i < entries; i++){
        chain->GetEntry(i);

        if((i+1)%10000 == 0 || (entries-i+1)<10000){
            cout<<"\rClustering Trigger Events: " << i+1 << "/" << entries << flush;
            if(i+1 == entries){
                cout<<endl;
            }
        }

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
            fill_rand_Histos_wThr(i);

            for(Int_t m = 0; m < nSSPBits - 1; m++){
                h_rand_3Clust_Sum_Ratio[m] = (TH1F*) h_rand_3Clust_VTP_Sum[m]->Clone();
                h_rand_3Clust_Sum_Ratio[m]->Divide(h_rand_3Clust_All_Sum);
                h_rand_3Clust_Sum_Ratio[m]->SetTitle("Trigger Efficiency");
                TString sumRatioName  = TString::Format("h_rand_3Clust_Sum_Ratio_bit%d", m);
                h_rand_3Clust_Sum_Ratio[m]->SetName(sumRatioName);

                h_rand_3Clust_Max_Ratio[m] = (TH1F*) h_rand_3Clust_VTP_Max[m]->Clone();
                h_rand_3Clust_Max_Ratio[m]->Divide(h_rand_3Clust_All_Max);
                h_rand_3Clust_Max_Ratio[m]->SetTitle("Trigger Efficiency");
                TString maxRatioName  = TString::Format("h_rand_3Clust_Max_Ratio_bit%d", m);
                h_rand_3Clust_Max_Ratio[m]->SetName(maxRatioName);

                h_rand_3Clust_Min_Ratio[m] = (TH1F*) h_rand_3Clust_VTP_Min[m]->Clone();
                h_rand_3Clust_Min_Ratio[m]->Divide(h_rand_3Clust_All_Min);
                h_rand_3Clust_Min_Ratio[m]->SetTitle("Trigger Efficiency");
                TString minRatioName  = TString::Format("h_rand_3Clust_Min_Ratio_bit%d", m);
                h_rand_3Clust_Min_Ratio[m]->SetName(minRatioName);

                h_rand_3Clust_deltaT_Ratio[m] = (TH1F*) h_rand_3Clust_deltaT_VTP[m]->Clone();
                h_rand_3Clust_deltaT_Ratio[m]->Divide(h_rand_3Clust_deltaT_All);
                h_rand_3Clust_deltaT_Ratio[m]->SetTitle("Trigger Efficiency");
                TString dtRatioName  = TString::Format("h_rand_3Clust_deltaT_Ratio_bit%d",m);
                h_rand_3Clust_deltaT_Ratio[m]->SetName(dtRatioName);
            }

        }
        if(tSum && (!trig.GetTriggersBit(TSUM_FLAG-nSSPBits)->empty()) && !(trigger_bits & (1<<LMS_FLAG))){//(trigType & (1<<SSP_RAWSUM_TFLAG))
            fill_tSum_Histos_wThr(i);
        }
    }
}

/**
 * Fill the various histograms making cuts on the basic clustering thresholds.
 *
 * @param i - the current event number
 */
void ClustTrig::fill_tSum_Histos_wThr(Int_t i){
    Double_t clustSum = 0;
    Int_t clustAbove = 0;
    Double_t maxTime = -100000;
    Double_t minTime = 10000000;

    Double_t clust_min = 1000000;
    Int_t min_ind = -1;
    UChar_t clust_min_bl = 244;

    vector<Int_t> GEM_matches;
    vector<Int_t> GEM_matches_either;
    Int_t GEM_matches_count = 0;
    Int_t GEM_matches_count_both = 0;
    Utils::Point p[MAX_CLUSTERS];
    Double_t clust_max = -1;
    for(Int_t j = 0; j < nClust; j++){
        if(cl_E[j] > CL_IND_THR && cl_nblocks[j] > 1){
            clustSum += cl_E[j];
            clustAbove++;
            if(cl_E[j] > clust_max){
                clust_max = cl_E[j];
            }
            if(cl_E[j] < clust_min){
                clust_min = cl_E[j];
                clust_min_bl = cl_nblocks[j];
                min_ind = j;
            }
            if(cl_time[j] > maxTime){
                maxTime = cl_time[j];
            }
            if(cl_time[j] < minTime){
                minTime = cl_time[j];
            }

            h_tSum_Off_centerID->Fill(cl_center[j]);
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

            if(((match_flag[j] & (1<<0)) || (match_flag[j] & (1<<1))) || ((match_flag[j] & (1<<2)) || (match_flag[j] & (1<<3)))){
                GEM_matches_count++;
                GEM_matches_either.push_back(j);
                if(((match_flag[j] & (1<<0)) || (match_flag[j] & (1<<1))) && ((match_flag[j] & (1<<2)) || (match_flag[j] & (1<<3)))){
                    GEM_matches_count_both++;
                    GEM_matches.push_back(j);

                    vector<Float_t> x_vec;
                    vector<Float_t> y_vec;
                    vector<Float_t> z_vec;
                    /*x_vec.push_back(cl_x[j]);
                    y_vec.push_back(cl_y[j]);
                    z_vec.push_back(cl_z[j]);*/ //Currently the HyCal position still needs correction so they are not being used in the line of best fit yet.
                    for(Int_t q = 0; q < MAX_GEMS; q++){
                        if(mgz[j][q]>0){
                            x_vec.push_back(mgx[j][q]);
                            y_vec.push_back(mgy[j][q]);
                            z_vec.push_back(mgz[j][q]);
                        }
                    }

                    Utils::LineOfBestFit line = Utils::FitLine(x_vec, y_vec, z_vec);
                    p[j] = Utils::ClosestApproachToZAxis(line);
                }
            }
        }
    }

    h_tSum_VTP_NClustDiff->Fill(vtp_cl_n - clustAbove);
    bool check = false;

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
                h_tSum_VTP_NClustDiff_Missed->Fill(vtp_cl_n - clustAbove);
                h_tSum_2Clust_Miss_HC_XY_LowE->Fill(cl_x[min_ind],cl_y[min_ind]);

                check = true;

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

                for(Int_t b = 0; b < nClust; b++){
                    if(cl_E[b]>CL_IND_THR){
                        h_tSum_3Clust_Cl_Energy_All->Fill(cl_E[b]);
                        if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                            h_tSum_3Clust_Cl_Energy_Trig->Fill(cl_E[b]);
                        }
                        if(clustAbove == 3){
                            h_tSum_3ClustOnly_Cl_Energy_All->Fill(cl_E[b]);
                            if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                                h_tSum_3ClustOnly_Cl_Energy_Trig->Fill(cl_E[b]);
                            }
                        }
                    }
                }

                for(UInt_t c = 0; c < GEM_matches.size(); c++){
                    if(GEM_matches_count_both >= 2){
                        h_tSum_3Clust_GEMMatch_All_Cl_E->Fill(cl_E[GEM_matches.at(c)]);
                        if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                            h_tSum_3Clust_GEMMatch_Trig_Cl_E->Fill(cl_E[GEM_matches.at(c)]);
                        }

                        if(GEM_matches_count_both >= 3){
                            h_tSum_3Clust_GEMMatch_All_3Match_Both_Cl_E->Fill(cl_E[GEM_matches.at(c)]);
                            if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                                h_tSum_3Clust_GEMMatch_Trig_3Match_Both_Cl_E->Fill(cl_E[GEM_matches.at(c)]);
                            }
                        }
                    }

                    Float_t gemx_1 = TMath::Max(mgx[GEM_matches.at(c)][0], mgx[GEM_matches.at(c)][2]);
                    Float_t gemx_2 = TMath::Max(mgx[GEM_matches.at(c)][1], mgx[GEM_matches.at(c)][3]);

                    if(clustSum > 2000.0 && clustSum < 2400.0){
                        h_tSum_3Clust_vZ_High->Fill(p[GEM_matches.at(c)].z);
                        h_tSum_3Clust_DCA_High->Fill(TMath::Sqrt(p[GEM_matches.at(c)].x*p[GEM_matches.at(c)].x + p[GEM_matches.at(c)].y*p[GEM_matches.at(c)].y));
                        h_tSum_3Clust_delta_X_High->Fill(gemx_1-gemx_2);
                    }
                    else if(clustSum < 2000.0){
                        h_tSum_3Clust_vZ_Low->Fill(p[GEM_matches.at(c)].z);
                        h_tSum_3Clust_DCA_Low->Fill(TMath::Sqrt(p[GEM_matches.at(c)].x*p[GEM_matches.at(c)].x + p[GEM_matches.at(c)].y*p[GEM_matches.at(c)].y));
                        h_tSum_3Clust_delta_X_Low->Fill(gemx_1-gemx_2);
                    }    
                }

                for(UInt_t d = 0; d < GEM_matches_either.size(); d++){
                    if(GEM_matches_count >= 2){
                        h_tSum_3Clust_GEMMatch_All_2Match_Either_Cl_E->Fill(cl_E[GEM_matches_either.at(d)]);
                        if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                            h_tSum_3Clust_GEMMatch_Trig_2Match_Either_Cl_E->Fill(cl_E[GEM_matches_either.at(d)]);
                        }

                        if(GEM_matches_count >= 3){
                            h_tSum_3Clust_GEMMatch_All_3Match_Either_Cl_E->Fill(cl_E[GEM_matches_either.at(d)]);
                            if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                                h_tSum_3Clust_GEMMatch_Trig_3Match_Either_Cl_E->Fill(cl_E[GEM_matches_either.at(d)]);
                            }
                        }
                    }
                }

                if(GEM_matches_count_both >= 2){
                    h_tSum_3Clust_GEMMatch_All->Fill(clustSum);
                    if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                        h_tSum_3Clust_GEMMatch_Trig->Fill(clustSum);
                    }

                    if(GEM_matches_count_both >= 3){
                        h_tSum_3Clust_GEMMatch_All_3Match_Both->Fill(clustSum);
                        if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                            h_tSum_3Clust_GEMMatch_Trig_3Match_Both->Fill(clustSum);
                        }
                    }
                }

                if(GEM_matches_count >= 2){
                    h_tSum_3Clust_GEMMatch_All_2Match_Either->Fill(clustSum);
                    if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                        h_tSum_3Clust_GEMMatch_Trig_2Match_Either->Fill(clustSum);
                    }

                    if(GEM_matches_count >= 3){
                        h_tSum_3Clust_GEMMatch_All_3Match_Either->Fill(clustSum);
                        if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                            h_tSum_3Clust_GEMMatch_Trig_3Match_Either->Fill(clustSum);
                        }
                    }
                }

                if(nClust == 3 && clustAbove == 3){
                    if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                        h_tSum_3ClustOnly_VTP_800->Fill(clustSum);
                        h_tSum_3ClustOnly_VTP_800_MinE->Fill(clust_min);
                        h_tSum_3ClustOnly_Cl_Max_Trig->Fill(clust_max);
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
                    h_tSum_3ClustOnly_Cl_Max_All->Fill(clust_max);
                }

                if(!trig.GetTriggersBit(CL3_FLAG-nSSPBits)->empty()){
                    h_tSum_3Clust_VTP_800->Fill(clustSum);
                    h_tSum_3Clust_VTP_800_MinE->Fill(clust_min);
                    h_tSum_3Clust_Cl_Max_Trig->Fill(clust_max);
                }
                else{
                    h_tSum_3Clust_Missed_MaxTime->Fill(maxTime);
                    h_tSum_3Clust_Missed_timeDif->Fill(maxTime-minTime);
                    h_tSum_3Clust_Miss_HC_XY_LowE->Fill(cl_x[min_ind],cl_y[min_ind]);
                    if(!check){
                        h_tSum_VTP_NClustDiff_Missed->Fill(vtp_cl_n - clustAbove);
                    }

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
                h_tSum_3Clust_Cl_Max_All->Fill(clust_max);
            }
        }
    }

    for(Int_t a = 0; a < vtp_cl_n; a++){
        h_tSum_VTP_Timing->Fill(vtp_cl_time[a]);
        h_tSum_VTP_centerID->Fill(vtp_cl_center[a]);  
    }
}

/**
 * Sets up the histograms for this object when comparing to the total sum trigger.
 */
void ClustTrig::setup_tSum_Histos(){
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
    h_tSum_2ClustOnly_All = new TH1F("h_tSum_2ClustOnly_All", "Total Sum Event with ONLY 2 Clusters Found;E_{Sum};Counts",en+500,0, en+500);
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
    h_tSum_3Clust_Cl_Max_All = new TH1F("h_tSum_3Clust_Cl_Max_All", "Total Sum Event with at least 3 Clusters, by Max;E_{Max};Counts",en+500,0, en+500);
    h_tSum_3Clust_Cl_Max_Trig = new TH1F("h_tSum_3Clust_Cl_Max_Trig", "Total Sum Event with at least 3 Clusters, by Max;E_{Max};Counts",en+500,0, en+500);

    //Offline 3 Cluster events for events where only 3 clusters were found.
    h_tSum_3ClustOnly_All = new TH1F("h_tSum_3ClustOnly_All", "Total Sum Event with ONLY 3 Clusters Found;E_{Sum};Counts",en+500,0, en+500);
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
    //Maximum cluster energy for only 3 cluster events
    h_tSum_3ClustOnly_Cl_Max_All = new TH1F("h_tSum_3ClustOnly_Cl_Max_All", "Total Sum Event with ONLY 3 Clusters, by Max;E_{Max};Counts",en+500,0, en+500);
    h_tSum_3ClustOnly_Cl_Max_Trig = new TH1F("h_tSum_3ClustOnly_Cl_Max_Trig", "Total Sum Event with ONLY 3 Clusters, by Max;E_{Max};Counts",en+500,0, en+500);

    //Offline 2 Cluster events; HyCal XY Positions
    h_tSum_2Clust_All_HC_XY = new TH2F("h_tSum_2Clust_All_HC_XY", "Total Sum Events with 2 Offline Clusters Found;x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);
    //VTP 2 Cluster events; HyCal XY Positions
    h_tSum_2Clust_VTP_HC_XY = new TH2F("h_tSum_2Clust_VTP_HC_XY", "Total Sum Events with 2 VTP Clusters Found;x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);
    //2 Cluster events missed by the VTP Trigger; HyCal XY, and Theta Positions
    h_tSum_2Clust_Miss_HC_XY = new TH2F("h_tSum_2Clust_Miss_HC_XY", "Total Sum Events with 2 Clusters Found Missed by the VTP;x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);
    h_tSum_2Clust_Miss_HC_Theta = new TH1F("h_tSum_2Clust_Miss_HC_Theta", "Total Sum Events with 2 Clusters Found Missed by the VTP;#theta [#circ];Counts",50, 0, 5);
    h_tSum_2Clust_Miss_HC_XY_LowE = new TH2F("h_tSum_2Clust_Miss_HC_XY_LowE", "Total Sum Events with 2 Clusters Found Missed by the VTP (Lowest Energy Cluster Only);x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);

    //Offline 2 Cluster events; HyCal XY Positions
    h_tSum_3Clust_All_HC_XY = new TH2F("h_tSum_3Clust_All_HC_XY", "Total Sum Events with 3 Offline Clusters Found;x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);
    //VTP 3 Cluster events; HyCal XY Positions
    h_tSum_3Clust_VTP_HC_XY = new TH2F("h_tSum_3Clust_VTP_HC_XY", "Total Sum Events with 3 VTP Clusters Found;x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);
    //3 Cluster events missed by the VTP Trigger; HyCal XY, and Theta Positions
    h_tSum_3Clust_Miss_HC_XY = new TH2F("h_tSum_3Clust_Miss_HC_XY", "Total Sum Events with 3 Clusters Found Missed by the VTP;x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);
    h_tSum_3Clust_Miss_HC_Theta = new TH1F("h_tSum_3Clust_Miss_HC_Theta", "Total Sum Events with 3 Clusters Found Missed by the VTP;#theta [#circ];Counts",50, 0, 5);
    h_tSum_3Clust_Miss_HC_XY_LowE = new TH2F("h_tSum_3Clust_Miss_HC_XY_LowE", "Total Sum Events with 3 Clusters Found Missed by the VTP (Lowest Energy Cluster Only);x [mm];y [mm]",1000,-500, 500, 1000, -500, 500);

    //Number of VTP Clusters - Offline Clusters
    h_tSum_VTP_NClustDiff = new TH1I("h_tSum_VTP_NClustDiff", "(No. VTP Clusters) - (No. Offline Clusters);Difference;Count",11, -5.5, 5.5);
    //Number of VTP Clusters - Offline Clusters
    h_tSum_VTP_NClustDiff_Missed = new TH1I("h_tSum_VTP_NClustDiff_Missed", "(No. VTP Clusters) - (No. Offline Clusters) Missed Events;Difference;Count",11, -5.5, 5.5);

    h_tSum_VTP_Timing = new TH1D("h_tSum_VTP_Timing", "Timing of All VTP Clusters; Samples [4ns per];Count", kTrigTimeBins, kTrigTimeMin, kTrigTimeMax);

    //VTP center cluster module ID
    h_tSum_VTP_centerID = new TH1I("h_tSum_VTP_centerID","Cluster Center Module ID; ID; Count",1157, 1000.5, 2256.5);
    //Offline center cluster module ID
    h_tSum_Off_centerID = new TH1I("h_tSum_Off_centerID","Cluster Center Module ID; ID; Count",1157, 1000.5, 2256.5);

    //Individual cluster energy for at least 3 or more cluster events
    h_tSum_3Clust_Cl_Energy_All = new TH1F("h_tSum_3Clust_Cl_Energy_All", "Total Sum Event with at least 3 Clusters Found, Every Cluster Energy;E_{Cl};Counts",en+500,0, en+500);
    h_tSum_3Clust_Cl_Energy_Trig = new TH1F("h_tSum_3Clust_Cl_Energy_Trig", "Total Sum Event with at least 3 Clusters Found, Every Cluster Energy;E_{Cl};Counts",en+500,0, en+500);

    //Individual cluster energy for only 3 cluster events
    h_tSum_3ClustOnly_Cl_Energy_All = new TH1F("h_tSum_3ClustOnly_Cl_Energy_All", "Total Sum Event with ONLY 3 Clusters Found, Every Cluster Energy;E_{Cl};Counts",en+500,0, en+500);
    h_tSum_3ClustOnly_Cl_Energy_Trig = new TH1F("h_tSum_3ClustOnly_Cl_Energy_Trig", "Total Sum Event with ONLY 3 Clusters Found, Every Cluster Energy;E_{Cl};Counts",en+500,0, en+500);

    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    h_tSum_3Clust_GEMMatch_All = new TH1F("h_tSum_3Clust_GEMMatch_All", "Total Sum Event with at least 3 Clusters Found (2 GEM Matches on Both GEM), Cluster Sum Energy;E_{Sum};Counts",en+500,0, en+500);
    h_tSum_3Clust_GEMMatch_Trig = new TH1F("h_tSum_3Clust_GEMMatch_Trig", "Total Sum Event with at least 3 Clusters Found (2 GEM Matches on Both GEM), Cluster Sum Energy;E_{Sum};Counts",en+500,0, en+500);

    //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
    h_tSum_3Clust_GEMMatch_All_2Match_Either = new TH1F("h_tSum_3Clust_GEMMatch_All_2Match_Either", "Total Sum Event with at least 3 Clusters Found (2 GEM Matches on Either GEM), Cluster Sum Energy;E_{Sum};Counts",en+500,0, en+500);
    h_tSum_3Clust_GEMMatch_Trig_2Match_Either = new TH1F("h_tSum_3Clust_GEMMatch_Trig_2Match_Either", "Total Sum Event with at least 3 Clusters Found (2 GEM Matches on Either GEM), Cluster Sum Energy;E_{Sum};Counts",en+500,0, en+500);
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    h_tSum_3Clust_GEMMatch_All_3Match_Either = new TH1F("h_tSum_3Clust_GEMMatch_All_3Match_Either", "Total Sum Event with at least 3 Clusters Found (3 GEM Matches on Either GEM), Cluster Sum Energy;E_{Sum};Counts",en+500,0, en+500);
    h_tSum_3Clust_GEMMatch_Trig_3Match_Either = new TH1F("h_tSum_3Clust_GEMMatch_Trig_3Match_Either", "Total Sum Event with at least 3 Clusters Found (3 GEM Matches on Either GEM), Cluster Sum Energy;E_{Sum};Counts",en+500,0, en+500);
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    h_tSum_3Clust_GEMMatch_All_3Match_Both = new TH1F("h_tSum_3Clust_GEMMatch_All_3Match_Both", "Total Sum Event with at least 3 Clusters Found (3 GEM Matches on Both GEM), Cluster Sum Energy;E_{Sum};Counts",en+500,0, en+500);
    h_tSum_3Clust_GEMMatch_Trig_3Match_Both = new TH1F("h_tSum_3Clust_GEMMatch_Trig_3Match_Both", "Total Sum Event with at least 3 Clusters Found (3 GEM Matches on Both GEM), Cluster Sum Energy;E_{Sum};Counts",en+500,0, en+500);

    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    h_tSum_3Clust_GEMMatch_All_Cl_E = new TH1F("h_tSum_3Clust_GEMMatch_All_Cl_E", "Total Sum Event with at least 3 Clusters Found (2 GEM Matches on Both GEM), Single Cluster Energy;E_{Cl};Counts",en+500,0, en+500);
    h_tSum_3Clust_GEMMatch_Trig_Cl_E = new TH1F("h_tSum_3Clust_GEMMatch_Trig_Cl_E", "Total Sum Event with at least 3 Clusters Found (2 GEM Matches on Both GEM), Single Cluster Energy;E_{Cl};Counts",en+500,0, en+500);

    //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
    h_tSum_3Clust_GEMMatch_All_2Match_Either_Cl_E = new TH1F("h_tSum_3Clust_GEMMatch_All_2Match_Either_Cl_E", "Total Sum Event with at least 3 Clusters Found (2 GEM Matches on Either GEM), Single Cluster Energy;E_{Cl};Counts",en+500,0, en+500);
    h_tSum_3Clust_GEMMatch_Trig_2Match_Either_Cl_E = new TH1F("h_tSum_3Clust_GEMMatch_Trig_2Match_Either_Cl_E", "Total Sum Event with at least 3 Clusters Found (2 GEM Matches on Either GEM), Single Cluster Energy;E_{Cl};Counts",en+500,0, en+500);
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    h_tSum_3Clust_GEMMatch_All_3Match_Either_Cl_E = new TH1F("h_tSum_3Clust_GEMMatch_All_3Match_Either_Cl_E", "Total Sum Event with at least 3 Clusters Found (3 GEM Matches on Either GEM), Single Cluster Energy;E_{Cl};Counts",en+500,0, en+500);
    h_tSum_3Clust_GEMMatch_Trig_3Match_Either_Cl_E = new TH1F("h_tSum_3Clust_GEMMatch_Trig_3Match_Either_Cl_E", "Total Sum Event with at least 3 Clusters Found (3 GEM Matches on Either GEM), Single Cluster Energy;E_{Cl};Counts",en+500,0, en+500);
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    h_tSum_3Clust_GEMMatch_All_3Match_Both_Cl_E = new TH1F("h_tSum_3Clust_GEMMatch_All_3Match_Both_Cl_E", "Total Sum Event with at least 3 Clusters Found (3 GEM Matches on Both GEM), Single Cluster Energy;E_{Cl};Counts",en+500,0, en+500);
    h_tSum_3Clust_GEMMatch_Trig_3Match_Both_Cl_E = new TH1F("h_tSum_3Clust_GEMMatch_Trig_3Match_Both_Cl_E", "Total Sum Event with at least 3 Clusters Found (3 GEM Matches on Both GEM), Single Cluster Energy;E_{Cl};Counts",en+500,0, en+500);

    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches and a high/low cluster energy sum.
    h_tSum_3Clust_vZ_High = new TH1F("h_tSum_3Clust_vZ_High","Z Vertex (PCA);z (mm);Counts",1000,-6000,7000);
    h_tSum_3Clust_vZ_Low = new TH1F("h_tSum_3Clust_vZ_Low","Z Vertex (PCA);z (mm);Counts",1000,-6000,7000);

    h_tSum_3Clust_DCA_Low = new TH1F("h_tSum_3Clust_DCA_Low","DCA E_{Sum} < 2000 MeV;r (mm);Counts",120,0,60);
    h_tSum_3Clust_delta_X_Low = new TH1F("h_tSum_3Clust_delta_X_Low","GEM #Deltax  E_{Sum} < 2000 MeV; #Deltax (mm);Counts",121,-15.5,15.5);
    h_tSum_3Clust_DCA_High = new TH1F("h_tSum_3Clust_DCA_High","DCA 2000 MeV < E_{Sum} < 2400 MeV;r (mm);Counts",120,0,60);
    h_tSum_3Clust_delta_X_High = new TH1F("h_tSum_3Clust_delta_X_High","GEM #Deltax 2000 MeV < E_{Sum} < 2400 MeV;#Deltax (mm);Counts",121,-15.5,15.5);
}

/**
 * Prints out the trigger performance plots from the Total Sum trigger.\
 *
 * @param pdfName - the name of the PDF to print to
 */
void ClustTrig::printTSumPDF(TString pdfName){

    TString Bit0_Thr = "750 MeV";
    TString Bit1_Thr = "750 MeV";
    TString Bit2_Thr = "750 MeV";
    TString Bit3_Thr = "850 MeV";
    TString Bit4_Thr = "1000 MeV";
    TString Bit5_Thr = "600 MeV";
    TString Bit6_Thr = "900 MeV";

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
	legend->AddEntry(h_tSum_1Clust_VTP_900,Bit6_Thr, "l");
	legend->AddEntry(h_tSum_1Clust_VTP_800,Bit1_Thr,"l");//"800 MeV","l");
	legend->AddEntry(h_tSum_1Clust_VTP_600,Bit5_Thr,"l");
	legend->AddEntry(h_tSum_1Clust_VTP_500,Bit4_Thr,"l");//"500 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_1ClustRatio_900->Draw("HIST");
    h_tSum_1ClustRatio_800->Draw("HIST SAME");
    h_tSum_1ClustRatio_600->Draw("HIST SAME");
    h_tSum_1ClustRatio_500->Draw("HIST SAME");
    
    legend2->SetHeader("Legend","C"); // option "C" allows to center the header
	legend2->AddEntry(h_tSum_1ClustRatio_900,Bit6_Thr,"l");
	legend2->AddEntry(h_tSum_1ClustRatio_800,Bit1_Thr,"l");//"800 MeV","l");
	legend2->AddEntry(h_tSum_1ClustRatio_600,Bit5_Thr,"l");
	legend2->AddEntry(h_tSum_1ClustRatio_500,Bit4_Thr,"l");//"500 MeV","l");
    legend2->Draw();
    c->Print(pdfName + "(");
    c->Clear();
    legend->Clear();
    legend2->Clear();

    //Page 2
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
	legend->AddEntry(h_tSum_1Clust_VTP_900,Bit6_Thr, "l");
	legend->AddEntry(h_tSum_1Clust_VTP_800,Bit1_Thr,"l");//"800 MeV","l");
	legend->AddEntry(h_tSum_1Clust_VTP_600,Bit5_Thr,"l");
	legend->AddEntry(h_tSum_1Clust_VTP_500,Bit4_Thr,"l");//"500 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_1ClustRatio_900->SetAxisRange(0.9,1.0,"Y");
    h_tSum_1ClustRatio_900->Draw("HIST");
    h_tSum_1ClustRatio_800->Draw("HIST SAME");
    h_tSum_1ClustRatio_600->Draw("HIST SAME");
    h_tSum_1ClustRatio_500->Draw("HIST SAME");
    
    legend2->SetHeader("Legend","C"); // option "C" allows to center the header
	legend2->AddEntry(h_tSum_1ClustRatio_900,Bit6_Thr,"l");
	legend2->AddEntry(h_tSum_1ClustRatio_800,Bit1_Thr,"l");//"800 MeV","l");
	legend2->AddEntry(h_tSum_1ClustRatio_600,Bit5_Thr,"l");
	legend2->AddEntry(h_tSum_1ClustRatio_500,Bit4_Thr,"l");//"500 MeV","l");
    legend2->Draw();
    c->Print(pdfName + "(");
    c->Clear();
    legend->Clear();
    legend2->Clear();

    //Page 3
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
	legend->AddEntry(h_tSum_1ClustOnly_VTP_900,Bit6_Thr,"l");
	legend->AddEntry(h_tSum_1ClustOnly_VTP_800,Bit1_Thr,"l");//"800 MeV","l");
	legend->AddEntry(h_tSum_1ClustOnly_VTP_600,Bit5_Thr,"l");
	legend->AddEntry(h_tSum_1ClustOnly_VTP_500,Bit4_Thr,"l");//"500 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_1ClustOnlyRatio_900->Draw("HIST");
    h_tSum_1ClustOnlyRatio_800->Draw("HIST SAME");
    h_tSum_1ClustOnlyRatio_600->Draw("HIST SAME");
    h_tSum_1ClustOnlyRatio_500->Draw("HIST SAME");
    
    legend2->SetHeader("Legend","C"); // option "C" allows to center the header
	legend2->AddEntry(h_tSum_1ClustOnlyRatio_900,Bit6_Thr,"l");
	legend2->AddEntry(h_tSum_1ClustOnlyRatio_800,Bit1_Thr,"l");//"800 MeV","l");
	legend2->AddEntry(h_tSum_1ClustOnlyRatio_600,Bit5_Thr,"l");
	legend2->AddEntry(h_tSum_1ClustOnlyRatio_500,Bit4_Thr,"l");//"500 MeV","l");
    legend2->Draw();
    c->Print(pdfName);
    c->Clear();
    legend->Clear();
    legend2->Clear();

    //Page 4
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
	legend->AddEntry(h_tSum_1ClustOnly_VTP_900,Bit6_Thr,"l");
	legend->AddEntry(h_tSum_1ClustOnly_VTP_800,Bit1_Thr,"l");//"800 MeV","l");
	legend->AddEntry(h_tSum_1ClustOnly_VTP_600,Bit5_Thr,"l");
	legend->AddEntry(h_tSum_1ClustOnly_VTP_500,Bit4_Thr,"l");//"500 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_1ClustOnlyRatio_900->SetAxisRange(0.9,1.0,"Y");
    h_tSum_1ClustOnlyRatio_900->Draw("HIST");
    h_tSum_1ClustOnlyRatio_800->Draw("HIST SAME");
    h_tSum_1ClustOnlyRatio_600->Draw("HIST SAME");
    h_tSum_1ClustOnlyRatio_500->Draw("HIST SAME");
    
    legend2->SetHeader("Legend","C"); // option "C" allows to center the header
	legend2->AddEntry(h_tSum_1ClustOnlyRatio_900,Bit6_Thr,"l");
	legend2->AddEntry(h_tSum_1ClustOnlyRatio_800,Bit1_Thr,"l");//"800 MeV","l");
	legend2->AddEntry(h_tSum_1ClustOnlyRatio_600,Bit5_Thr,"l");
	legend2->AddEntry(h_tSum_1ClustOnlyRatio_500,Bit4_Thr,"l");//"500 MeV","l");
    legend2->Draw();
    c->Print(pdfName);
    c->Clear();
    legend->Clear();
    legend2->Clear();

    //Page 5
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2Clust_All->Draw("HIST");
    h_tSum_2Clust_VTP_800->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2Clust_All,"All","l");
	legend->AddEntry(h_tSum_2Clust_VTP_800,Bit2_Thr,"l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustRatio->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 6
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2Clust_All->Draw("HIST");
    h_tSum_2Clust_VTP_800->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2Clust_All,"All","l");
	legend->AddEntry(h_tSum_2Clust_VTP_800,Bit2_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustRatio->SetAxisRange(0.9,1.0, "Y");
    h_tSum_2ClustRatio->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 7
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2Clust_All_MinE->Draw("HIST");
    h_tSum_2Clust_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2Clust_All_MinE,"All","l");
	legend->AddEntry(h_tSum_2Clust_VTP_800_MinE,Bit2_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 8
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2Clust_All_MinE->Draw("HIST");
    h_tSum_2Clust_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2Clust_All_MinE,"All","l");
	legend->AddEntry(h_tSum_2Clust_VTP_800_MinE,Bit2_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustRatio_MinE->SetAxisRange(0.9, 1.0, "Y");
    h_tSum_2ClustRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 9
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

    //Page 10
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_2Clust_Missed_nBlocks_less->Draw("HIST");
    c->cd(2);
    h_tSum_2Clust_Missed_nBlocks_greater->Draw("HIST");
    c->Print(pdfName);
    c->Clear();


    //Page 11
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2ClustOnly_All->Draw("HIST");
    h_tSum_2ClustOnly_VTP_800->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2ClustOnly_All,"All","l");
	legend->AddEntry(h_tSum_2ClustOnly_VTP_800,Bit2_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustOnlyRatio->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 12
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2ClustOnly_All->Draw("HIST");
    h_tSum_2ClustOnly_VTP_800->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2ClustOnly_All,"All","l");
	legend->AddEntry(h_tSum_2ClustOnly_VTP_800,Bit2_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustOnlyRatio->SetAxisRange(0.9, 1.0, "Y");
    h_tSum_2ClustOnlyRatio->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 13
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2ClustOnly_All_MinE->Draw("HIST");
    h_tSum_2ClustOnly_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2ClustOnly_All_MinE,"All","l");
	legend->AddEntry(h_tSum_2ClustOnly_VTP_800_MinE,Bit2_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustOnlyRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 14
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_2ClustOnly_All_MinE->Draw("HIST");
    h_tSum_2ClustOnly_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_2ClustOnly_All_MinE,"All","l");
	legend->AddEntry(h_tSum_2ClustOnly_VTP_800_MinE,Bit2_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_2ClustOnlyRatio_MinE->SetAxisRange(0.9, 1.0, "Y");
    h_tSum_2ClustOnlyRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 15
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

    //Page 16
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_2ClustOnly_Missed_nBlocks_less->Draw("HIST");
    c->cd(2);
    h_tSum_2ClustOnly_Missed_nBlocks_greater->Draw("HIST");
    c->Print(pdfName);
    c->Clear();

    //Page 17
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

    //Page 18
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    gPad->SetLogz(1);
    h_tSum_2Clust_Miss_HC_XY_LowE->Draw("COLZ");
    gPad->SetLogz(0);
    c->cd(2);
    TH1F* x_temp = (TH1F*) h_tSum_2Clust_VTP_HC_XY->ProjectionX();
    x_temp->Divide((TH1F*) h_tSum_2Clust_All_HC_XY->ProjectionX());
    x_temp->SetTitle("2 Cluster X Position Ratio of VTP/All");
    x_temp->Draw("HIST");
    c->cd(3);
    TH1F* y_temp = (TH1F*) h_tSum_2Clust_VTP_HC_XY->ProjectionY();
    y_temp->Divide((TH1F*) h_tSum_2Clust_All_HC_XY->ProjectionY());
    y_temp->SetTitle("2 Cluster Y Position Ratio of VTP/All");
    y_temp->Draw("HIST");
    c->Print(pdfName);
    c->Clear();
    gPad->SetLogy(0);
    gPad->SetLogz(0);

    //Page 19
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3Clust_All->Draw("HIST");
    h_tSum_3Clust_VTP_800->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3Clust_All,"All","l");
	legend->AddEntry(h_tSum_3Clust_VTP_800,Bit3_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_3ClustRatio->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 20
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3Clust_All_MinE->Draw("HIST");
    h_tSum_3Clust_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3Clust_All_MinE,"All","l");
	legend->AddEntry(h_tSum_3Clust_VTP_800_MinE,Bit3_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_3ClustRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 21
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3Clust_All_MinE->Draw("HIST");
    h_tSum_3Clust_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3Clust_All_MinE,"All","l");
	legend->AddEntry(h_tSum_3Clust_VTP_800_MinE,Bit3_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_3ClustRatio_MinE->SetAxisRange(0.9,1.0,"Y");
    h_tSum_3ClustRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 22
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3Clust_Cl_Max_All->Draw("HIST");
    h_tSum_3Clust_Cl_Max_Trig->SetLineColor(kBlue);
    h_tSum_3Clust_Cl_Max_Trig->Draw("HIST SAME");

    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3Clust_Cl_Max_All,"All","l");
	legend->AddEntry(h_tSum_3Clust_Cl_Max_Trig,Bit3_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    TH1F* rat = (TH1F*) h_tSum_3Clust_Cl_Max_Trig->Clone();
    rat->Divide(h_tSum_3Clust_Cl_Energy_All);
    rat->SetTitle("Trigger Efficiency as a fucntion of Maximum Cluster Energy");
    rat->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();
    
    //Page 23
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

    //Page 24
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_3Clust_Missed_nBlocks_less->Draw("HIST");
    c->cd(2);
    h_tSum_3Clust_Missed_nBlocks_greater->Draw("HIST");
    c->Print(pdfName);
    c->Clear();

    //Page 25
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

    //Page 26
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    gPad->SetLogz(1);
    h_tSum_3Clust_Miss_HC_XY_LowE->Draw("COLZ");
    gPad->SetLogz(0);
    c->cd(2);
    TH1F* x_temp2 = (TH1F*) h_tSum_3Clust_VTP_HC_XY->ProjectionX();
    x_temp2->Divide((TH1F*) h_tSum_3Clust_All_HC_XY->ProjectionX());
    x_temp2->SetTitle("3 Cluster X Position Ratio of VTP/All");
    x_temp2->Draw("HIST");
    c->cd(3);
    TH1F* y_temp2 = (TH1F*) h_tSum_3Clust_VTP_HC_XY->ProjectionY();
    y_temp2->Divide((TH1F*) h_tSum_3Clust_All_HC_XY->ProjectionY());
    y_temp2->SetTitle("3 Cluster Y Position Ratio of VTP/All");
    y_temp2->Draw("HIST");
    c->Print(pdfName);
    c->Clear();
    gPad->SetLogy(0);
    gPad->SetLogz(0);

    //Page 27
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3Clust_Cl_Energy_All->Draw("HIST");
    h_tSum_3Clust_Cl_Energy_Trig->SetLineColor(kBlue);
    h_tSum_3Clust_Cl_Energy_Trig->Draw("HIST SAME");

    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3Clust_Cl_Energy_All,"All","l");
	legend->AddEntry(h_tSum_3Clust_Cl_Energy_Trig,Bit3_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    gPad->SetLogy(1);
    h_tSum_3Clust_Cl_Energy_All->Draw("HIST");
    h_tSum_3Clust_Cl_Energy_Trig->SetLineColor(kBlue);
    h_tSum_3Clust_Cl_Energy_Trig->Draw("HIST SAME");

    legend2->SetHeader("Legend","C"); // option "C" allows to center the header
	legend2->AddEntry(h_tSum_3Clust_Cl_Energy_All,"All","l");
	legend2->AddEntry(h_tSum_3Clust_Cl_Energy_Trig,Bit3_Thr,"l");//"800 MeV","l");
	legend2->Draw();

    c->Print(pdfName);
    gPad->SetLogy(0);
    c->Clear();
    legend->Clear();
    legend->Clear();
    legend2->Clear();

    //Page 28
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3Clust_GEMMatch_All_2Match_Either->Draw("HIST");
    h_tSum_3Clust_GEMMatch_Trig_2Match_Either->SetLineColor(kBlue);
    h_tSum_3Clust_GEMMatch_Trig_2Match_Either->Draw("HIST SAME");
    legend2->SetHeader("Legend","C"); // option "C" allows to center the header
	legend2->AddEntry(h_tSum_3Clust_GEMMatch_All_2Match_Either,"All","l");
	legend2->AddEntry(h_tSum_3Clust_GEMMatch_Trig_2Match_Either,Bit3_Thr,"l");//"800 MeV","l");
    legend2->Draw();

    c->cd(2);
    h_tSum_3Clust_GEMMatch_All_2Match_Either_Cl_E->Draw("HIST");
    h_tSum_3Clust_GEMMatch_Trig_2Match_Either_Cl_E->SetLineColor(kBlue);
    h_tSum_3Clust_GEMMatch_Trig_2Match_Either_Cl_E->Draw("HIST SAME");
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3Clust_GEMMatch_All_2Match_Either_Cl_E,"All","l");
	legend->AddEntry(h_tSum_3Clust_GEMMatch_Trig_2Match_Either_Cl_E,Bit3_Thr,"l");//"800 MeV","l");
    legend->Draw();
    c->Print(pdfName);
    legend->Clear();
    legend2->Clear();
    c->Clear();

    //Page 29
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3Clust_GEMMatch_All->Draw("HIST");
    h_tSum_3Clust_GEMMatch_Trig->SetLineColor(kBlue);
    h_tSum_3Clust_GEMMatch_Trig->Draw("HIST SAME");
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3Clust_GEMMatch_All,"All","l");
	legend->AddEntry(h_tSum_3Clust_GEMMatch_Trig,Bit3_Thr,"l");//"800 MeV","l");
    legend->Draw();

    c->cd(2);
    h_tSum_3Clust_GEMMatch_All_Cl_E->Draw("HIST");
    h_tSum_3Clust_GEMMatch_Trig_Cl_E->SetLineColor(kBlue);
    h_tSum_3Clust_GEMMatch_Trig_Cl_E->Draw("HIST SAME");
    legend2->SetHeader("Legend","C"); // option "C" allows to center the header
	legend2->AddEntry(h_tSum_3Clust_GEMMatch_All_Cl_E,"All","l");
	legend2->AddEntry(h_tSum_3Clust_GEMMatch_Trig_Cl_E,Bit3_Thr,"l");//"800 MeV","l");
    legend2->Draw();

    c->Print(pdfName);
    c->Clear();
    legend->Clear();
    legend2->Clear();

    //Page 30
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3Clust_GEMMatch_All_3Match_Either->Draw("HIST");
    h_tSum_3Clust_GEMMatch_Trig_3Match_Either->SetLineColor(kBlue);
    h_tSum_3Clust_GEMMatch_Trig_3Match_Either->Draw("HIST SAME");
    legend2->SetHeader("Legend","C"); // option "C" allows to center the header
	legend2->AddEntry(h_tSum_3Clust_GEMMatch_All_3Match_Either,"All","l");
	legend2->AddEntry(h_tSum_3Clust_GEMMatch_Trig_3Match_Either,Bit3_Thr,"l");//"800 MeV","l");
    legend2->Draw();

    c->cd(2);
    h_tSum_3Clust_GEMMatch_All_3Match_Either_Cl_E->Draw("HIST");
    h_tSum_3Clust_GEMMatch_Trig_3Match_Either_Cl_E->SetLineColor(kBlue);
    h_tSum_3Clust_GEMMatch_Trig_3Match_Either_Cl_E->Draw("HIST SAME");
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3Clust_GEMMatch_All_3Match_Either_Cl_E,"All","l");
	legend->AddEntry(h_tSum_3Clust_GEMMatch_Trig_3Match_Either_Cl_E,Bit3_Thr,"l");//"800 MeV","l");
    legend->Draw();
    c->Print(pdfName);
    legend->Clear();
    legend2->Clear();
    c->Clear();

    //Page 31
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3Clust_GEMMatch_All_3Match_Both->Draw("HIST");
    h_tSum_3Clust_GEMMatch_Trig_3Match_Both->SetLineColor(kBlue);
    h_tSum_3Clust_GEMMatch_Trig_3Match_Both->Draw("HIST SAME");
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3Clust_GEMMatch_All_3Match_Both,"All","l");
	legend->AddEntry(h_tSum_3Clust_GEMMatch_Trig_3Match_Both,Bit3_Thr,"l");//"800 MeV","l");
    legend->Draw();

    c->cd(2);
    h_tSum_3Clust_GEMMatch_All_3Match_Both_Cl_E->Draw("HIST");
    h_tSum_3Clust_GEMMatch_Trig_3Match_Both_Cl_E->SetLineColor(kBlue);
    h_tSum_3Clust_GEMMatch_Trig_3Match_Both_Cl_E->Draw("HIST SAME");
    legend2->SetHeader("Legend","C"); // option "C" allows to center the header
	legend2->AddEntry(h_tSum_3Clust_GEMMatch_All_3Match_Both_Cl_E,"All","l");
	legend2->AddEntry(h_tSum_3Clust_GEMMatch_Trig_3Match_Both_Cl_E,Bit3_Thr,"l");//"800 MeV","l");
    legend2->Draw();

    c->Print(pdfName);
    c->Clear();
    legend->Clear();
    legend2->Clear();

    //Page 32
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3Clust_vZ_High->Draw("HIST");
    h_tSum_3Clust_vZ_Low->SetLineColor(kBlue);
    h_tSum_3Clust_vZ_Low->Draw("HIST SAME");

    legend2->SetHeader("Legend","C"); // option "C" allows to center the header
	legend2->AddEntry(h_tSum_3Clust_vZ_High,"2000<E_{Sum}<2400","l");
	legend2->AddEntry(h_tSum_3Clust_vZ_Low,"E_{Sum}<2000","l");//"800 MeV","l");
	legend2->Draw();

    c->Print(pdfName);
    c->Clear();
    legend2->Clear();

    //Page 33
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_3Clust_DCA_High->Draw("HIST");

    c->cd(2);
    h_tSum_3Clust_DCA_Low->Draw("HIST");

    c->cd(3);
    gPad->SetLogy(1);
    h_tSum_3Clust_delta_X_High->Draw("HIST");

    c->cd(4);
    gPad->SetLogy(1);
    h_tSum_3Clust_delta_X_Low->Draw("HIST");

    c->Print(pdfName);
    c->Clear();

    //Page 34
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3ClustOnly_All->Draw("HIST");
    h_tSum_3ClustOnly_VTP_800->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3ClustOnly_All,"All","l");
	legend->AddEntry(h_tSum_3ClustOnly_VTP_800,Bit3_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_3ClustOnlyRatio->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 35
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3ClustOnly_All_MinE->Draw("HIST");
    h_tSum_3ClustOnly_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3ClustOnly_All_MinE,"All","l");
	legend->AddEntry(h_tSum_3ClustOnly_VTP_800_MinE,Bit3_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_3ClustOnlyRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 36
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3ClustOnly_All_MinE->Draw("HIST");
    h_tSum_3ClustOnly_VTP_800_MinE->Draw("HIST SAME");
    
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3ClustOnly_All_MinE,"All","l");
	legend->AddEntry(h_tSum_3ClustOnly_VTP_800_MinE,Bit3_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    h_tSum_3ClustOnlyRatio_MinE->SetAxisRange(0.9, 1.0, "Y");
    h_tSum_3ClustOnlyRatio_MinE->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 37
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3ClustOnly_Cl_Max_All->Draw("HIST");
    h_tSum_3ClustOnly_Cl_Max_Trig->SetLineColor(kBlue);
    h_tSum_3ClustOnly_Cl_Max_Trig->Draw("HIST SAME");

    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3ClustOnly_Cl_Max_All,"All","l");
	legend->AddEntry(h_tSum_3ClustOnly_Cl_Max_Trig,Bit3_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    TH1F* rat2 = (TH1F*) h_tSum_3ClustOnly_Cl_Max_Trig->Clone();
    rat2->Divide(h_tSum_3ClustOnly_Cl_Energy_All);
    rat2->SetTitle("Trigger Efficiency as a fucntion of Maximum Cluster Energy");
    rat2->Draw("HIST");

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 38
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_3ClustOnly_Cl_Energy_All->Draw("HIST");
    h_tSum_3ClustOnly_Cl_Energy_Trig->SetLineColor(kBlue);
    h_tSum_3ClustOnly_Cl_Energy_Trig->Draw("HIST SAME");

    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_tSum_3ClustOnly_Cl_Energy_All,"All","l");
	legend->AddEntry(h_tSum_3ClustOnly_Cl_Energy_Trig,Bit3_Thr,"l");//"800 MeV","l");
	legend->Draw();

    c->cd(2);
    gPad->SetLogy(1);
    h_tSum_3ClustOnly_Cl_Energy_All->Draw("HIST");
    h_tSum_3ClustOnly_Cl_Energy_Trig->SetLineColor(kBlue);
    h_tSum_3ClustOnly_Cl_Energy_Trig->Draw("HIST SAME");

    legend2->SetHeader("Legend","C"); // option "C" allows to center the header
	legend2->AddEntry(h_tSum_3ClustOnly_Cl_Energy_All,"All","l");
	legend2->AddEntry(h_tSum_3ClustOnly_Cl_Energy_Trig,Bit3_Thr,"l");//"800 MeV","l");
	legend2->Draw();

    c->Print(pdfName);
    gPad->SetLogy(0);
    c->Clear();
    legend->Clear();
    legend2->Clear();

    //Page 39
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

    //Page 40
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_3ClustOnly_Missed_nBlocks_less->Draw("HIST");
    c->cd(2);
    h_tSum_3ClustOnly_Missed_nBlocks_greater->Draw("HIST");
    c->Print(pdfName);
    c->Clear();

    //Page 41
    //-----------------------------------------------------------------------------------
    c->Divide(2,4);
    for(Int_t b = 0; b < nSSPBits; b++){
        c->cd(b+1);
        hTrigTime[b]->Draw("HIST");
        hTrigTime[b]->SetStats(0);
    }
    c->Print(pdfName);
    c->Clear();

    //Page 42
    //-----------------------------------------------------------------------------------
    c->Print(pdfName);

    //Page 43
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_tSum_VTP_NClustDiff->Draw("HIST");

    c->cd(2);
    h_tSum_VTP_NClustDiff_Missed->Draw("HIST");

    c->cd(3);
    h_tSum_VTP_Timing->Draw("HIST");

    c->cd(4);
    gPad->SetLogy(1);
    h_tSum_VTP_centerID->Draw("HIST");
    h_tSum_Off_centerID->SetLineColor(kBlue);
    h_tSum_Off_centerID->Draw("HIST SAME");
    c->Print(pdfName + ")");
    gPad->SetLogy(0);
    c->Clear();

    delete c;

}

void ClustTrig::SaveAllTSumHistograms(TString rootfileName){

    TObjArray* arr = new TObjArray(0,0);

    //Offline 1 Cluster events (max energy cluster).
    (*arr).Add(h_tSum_1Clust_All);
    //VTP 1 Cluster Trigger (max energy cluster).
    (*arr).Add(h_tSum_1Clust_VTP_800);
    (*arr).Add(h_tSum_1Clust_VTP_900);
    (*arr).Add(h_tSum_1Clust_VTP_600);
    (*arr).Add(h_tSum_1Clust_VTP_500);

        //Offline 1 Cluster Trigger for events where only 1 cluster was found.
    (*arr).Add(h_tSum_1ClustOnly_All);
        //VTP 1 Cluster Trigger for events where only 1 cluster was found.
    (*arr).Add(h_tSum_1ClustOnly_VTP_800);
    (*arr).Add(h_tSum_1ClustOnly_VTP_900);
    (*arr).Add(h_tSum_1ClustOnly_VTP_600);
    (*arr).Add(h_tSum_1ClustOnly_VTP_500);

    //Offline 2 Cluster events (max energy combination).
    (*arr).Add(h_tSum_2Clust_All);
    //VTP 2 Cluster events (max energy combination).
    (*arr).Add(h_tSum_2Clust_VTP_800);
    //Offline 2 Cluster events for events where at least 2 clusters were found (Minimum Cluster E).
    (*arr).Add(h_tSum_2Clust_All_MinE);
    //VTP 2 Cluster events for events where at least 2 clusters were found (Minimum Cluster E).
    (*arr).Add(h_tSum_2Clust_VTP_800_MinE);
    //Maximum Time difference of 2 Cluster events for events where at least 2 clusters were found.
    (*arr).Add(h_tSum_2Clust_timeDif);
    //Maximum Time difference of 2 Cluster events for events where at least 2 clusters were found for events missed by the clustering trigger.
    (*arr).Add(h_tSum_2Clust_Missed_timeDif);
    //Maximum Cluster Time for 2 Cluster events for events where at least 2 clusters were found and the clustering trigger did not find it.
    (*arr).Add(h_tSum_2Clust_Missed_MaxTime);
    //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where at least 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    (*arr).Add(h_tSum_2Clust_Missed_EdeltaT);
    //Number of Blocks for Missed Events with Delta t < 16 ns where at least 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    (*arr).Add(h_tSum_2Clust_Missed_nBlocks_less);
    //Number of Blocks for Missed Events with Delta t > 16 ns where at least 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    (*arr).Add(h_tSum_2Clust_Missed_nBlocks_greater);

    //Offline 2 Cluster events for events where only 2 clusters were found.
    (*arr).Add(h_tSum_2ClustOnly_All);
    //VTP 2 Cluster events for events where only 2 clusters were found.
    (*arr).Add(h_tSum_2ClustOnly_VTP_800);
    //Offline 2 Cluster events for events where only 2 clusters were found (Minimum Cluster E).
    (*arr).Add(h_tSum_2ClustOnly_All_MinE);
        //VTP 2 Cluster events for events where only 2 clusters were found (Minimum Cluster E).
    (*arr).Add(h_tSum_2ClustOnly_VTP_800_MinE);
        //Time difference of 2 Cluster events for events where only 2 clusters were found.
    (*arr).Add(h_tSum_2ClustOnly_timeDif);
        //Maximum Time difference of 2 Cluster events for events where only 2 clusters were found for events missed by the clustering trigger.
    (*arr).Add(h_tSum_2ClustOnly_Missed_timeDif);
        //Maximum Cluster Time for 2 Cluster events for events where only 2 clusters were found and the clustering trigger did not find it.
    (*arr).Add(h_tSum_2ClustOnly_Missed_MaxTime);
        //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where only 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    (*arr).Add(h_tSum_2ClustOnly_Missed_EdeltaT);
        //Number of Blocks for Missed Events with Delta t < 16 ns where only 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    (*arr).Add(h_tSum_2ClustOnly_Missed_nBlocks_less);
        //Number of Blocks for Missed Events with Delta t > 16 ns where only 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    (*arr).Add(h_tSum_2ClustOnly_Missed_nBlocks_greater);

    //Offline 3 Cluster events (max energy combination).
    (*arr).Add(h_tSum_3Clust_All);
    //VTP 3 Cluster events (max energy combination).
    (*arr).Add(h_tSum_3Clust_VTP_800);
    //Offline 3 Cluster events for events where at least 3 clusters were found (Minimum Cluster E).
    (*arr).Add(h_tSum_3Clust_All_MinE);
    //VTP 3 Cluster events for events where at least 3 clusters were found (Minimum Cluster E).
    (*arr).Add(h_tSum_3Clust_VTP_800_MinE);
    //Maximum Time difference of 3 Cluster events for events where at least 3 clusters were found.
    (*arr).Add(h_tSum_3Clust_timeDif);
    //Maximum Time difference of 3 Cluster events for events where at least 3 clusters were found for events missed by the clustering trigger.
    (*arr).Add(h_tSum_3Clust_Missed_timeDif);
    //Maximum Cluster Time for 3 Cluster events for events where at least 3 clusters were found and the clustering trigger did not find it.
    (*arr).Add(h_tSum_3Clust_Missed_MaxTime);
    //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where at least 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    (*arr).Add(h_tSum_3Clust_Missed_EdeltaT);
    //Number of Blocks for Missed Events with Delta t < 16 ns where at least 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    (*arr).Add(h_tSum_3Clust_Missed_nBlocks_less);
    //Number of Blocks for Missed Events with Delta t > 16 ns where at least 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    (*arr).Add(h_tSum_3Clust_Missed_nBlocks_greater);
    //Individual cluster energy for at least 3 or more cluster events
    (*arr).Add(h_tSum_3Clust_Cl_Energy_All);
    (*arr).Add(h_tSum_3Clust_Cl_Energy_Trig);
    //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
    (*arr).Add(h_tSum_3Clust_GEMMatch_All);
    (*arr).Add(h_tSum_3Clust_GEMMatch_Trig);
    //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
    (*arr).Add(h_tSum_3Clust_GEMMatch_All_2Match_Either);
    (*arr).Add(h_tSum_3Clust_GEMMatch_Trig_2Match_Either);
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    (*arr).Add(h_tSum_3Clust_GEMMatch_All_3Match_Either);
    (*arr).Add(h_tSum_3Clust_GEMMatch_Trig_3Match_Either);
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    (*arr).Add(h_tSum_3Clust_GEMMatch_All_3Match_Both);
    (*arr).Add(h_tSum_3Clust_GEMMatch_Trig_3Match_Both);
    //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
    (*arr).Add(h_tSum_3Clust_GEMMatch_All_Cl_E);
    (*arr).Add(h_tSum_3Clust_GEMMatch_Trig_Cl_E);
    //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
    (*arr).Add(h_tSum_3Clust_GEMMatch_All_2Match_Either_Cl_E);
    (*arr).Add(h_tSum_3Clust_GEMMatch_Trig_2Match_Either_Cl_E);
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    (*arr).Add(h_tSum_3Clust_GEMMatch_All_3Match_Either_Cl_E);
    (*arr).Add(h_tSum_3Clust_GEMMatch_Trig_3Match_Either_Cl_E);
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    (*arr).Add(h_tSum_3Clust_GEMMatch_All_3Match_Both_Cl_E);
    (*arr).Add(h_tSum_3Clust_GEMMatch_Trig_3Match_Both_Cl_E);
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches and a high/low cluster energy sum.
    (*arr).Add(h_tSum_3Clust_vZ_High);
    (*arr).Add(h_tSum_3Clust_vZ_Low);
    //Maximum cluster energy for only 3 cluster events
    (*arr).Add(h_tSum_3Clust_Cl_Max_All);
    (*arr).Add(h_tSum_3Clust_Cl_Max_Trig);


    //Offline 3 Cluster events for events where only 3 clusters were found.
    (*arr).Add(h_tSum_3ClustOnly_All);
    //VTP 3 Cluster events for events where only 3 clusters were found.
    (*arr).Add(h_tSum_3ClustOnly_VTP_800);
    //Offline 3 Cluster events for events where only 3 clusters were found (Minimum Cluster E).
    (*arr).Add(h_tSum_3ClustOnly_All_MinE);
    //VTP 3 Cluster events for events where only 3 clusters were found (Minimum Cluster E).
    (*arr).Add(h_tSum_3ClustOnly_VTP_800_MinE);
    //Time difference of 3 Cluster events for events where only 3 clusters were found.
    (*arr).Add(h_tSum_3ClustOnly_timeDif);
    //Maximum Time difference of 3 Cluster events for events where only 3 clusters were found for events missed by the clustering trigger.
    (*arr).Add(h_tSum_3ClustOnly_Missed_timeDif);
    //Maximum Cluster Time for 3 Cluster events for events where only 3 clusters were found and the clustering trigger did not find it.
    (*arr).Add(h_tSum_3ClustOnly_Missed_MaxTime);
    //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where only 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    (*arr).Add(h_tSum_3ClustOnly_Missed_EdeltaT);
    //Number of Blocks for Missed Events with Delta t < 16 ns where at only 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    (*arr).Add(h_tSum_3ClustOnly_Missed_nBlocks_less);
    //Number of Blocks for Missed Events with Delta t > 16 ns where only 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    (*arr).Add(h_tSum_3ClustOnly_Missed_nBlocks_greater);
    //Individual cluster energy for only 3 cluster events
    (*arr).Add(h_tSum_3ClustOnly_Cl_Energy_All);
    (*arr).Add(h_tSum_3ClustOnly_Cl_Energy_Trig);
    //Maximum cluster energy for only 3 cluster events
    (*arr).Add(h_tSum_3ClustOnly_Cl_Max_All);
    (*arr).Add(h_tSum_3ClustOnly_Cl_Max_Trig);

    //Offline 2 Cluster events; HyCal XY Positions
    (*arr).Add(h_tSum_2Clust_All_HC_XY);
    //VTP 2 Cluster events; HyCal XY Positions
    (*arr).Add(h_tSum_2Clust_VTP_HC_XY);
    //2 Cluster events missed by the VTP Trigger; HyCal XY, and Theta Positions
    (*arr).Add(h_tSum_2Clust_Miss_HC_XY);
    (*arr).Add(h_tSum_2Clust_Miss_HC_Theta);
    (*arr).Add(h_tSum_2Clust_Miss_HC_XY_LowE);

    //Number of VTP Clusters
    (*arr).Add(h_tSum_VTP_NClustDiff);
    //Number of VTP Clusters - Offline Clusters for rejected events.
    (*arr).Add(h_tSum_VTP_NClustDiff_Missed);
    //VTP timing of events with non-zero difference in cluster times
    (*arr).Add(h_tSum_VTP_Timing);
        
    //VTP center cluster module ID
    (*arr).Add(h_tSum_VTP_centerID);
    //Offline center cluster module ID
    (*arr).Add(h_tSum_Off_centerID);
        

    //Offline 3 Cluster events; HyCal XY Positions
    (*arr).Add(h_tSum_3Clust_All_HC_XY);
    //VTP 3 Cluster events; HyCal XY Positions
    (*arr).Add(h_tSum_3Clust_VTP_HC_XY);
    //3 Cluster events missed by the VTP Trigger; HyCal XY, and Theta Positions
    (*arr).Add(h_tSum_3Clust_Miss_HC_XY);
    (*arr).Add(h_tSum_3Clust_Miss_HC_Theta);
    (*arr).Add(h_tSum_3Clust_Miss_HC_XY_LowE);

    //3 Cluster events
    (*arr).Add(h_tSum_3Clust_DCA_Low);
    (*arr).Add(h_tSum_3Clust_delta_X_Low);
    (*arr).Add(h_tSum_3Clust_DCA_High);
    (*arr).Add(h_tSum_3Clust_delta_X_High);

    TFile file1(rootfileName,"RECREATE");
    (*arr).Write();
	file1.Close();
}

/**
 * Deletes the total sum histograms.
 */
void ClustTrig::delete_tSum_Histos(){
    //Offline 1 Cluster events (max energy cluster).
    delete h_tSum_1Clust_All;
    //VTP 1 Cluster Trigger (max energy cluster).
    delete h_tSum_1Clust_VTP_800;
    delete h_tSum_1Clust_VTP_900;
    delete h_tSum_1Clust_VTP_600;
    delete h_tSum_1Clust_VTP_500;

    //Offline 1 Cluster Trigger for events where only 1 cluster was found.
    delete h_tSum_1ClustOnly_All;
    //VTP 1 Cluster Trigger for events where only 1 cluster was found.
    delete h_tSum_1ClustOnly_VTP_800;
    delete h_tSum_1ClustOnly_VTP_900;
    delete h_tSum_1ClustOnly_VTP_600;
    delete h_tSum_1ClustOnly_VTP_500;

    //Offline 2 Cluster events (max energy combination).
    delete h_tSum_2Clust_All;
    //VTP 2 Cluster events (max energy combination).
    delete h_tSum_2Clust_VTP_800;
    //Offline 2 Cluster events for events where at least 2 clusters were found (Minimum Cluster E).
    delete h_tSum_2Clust_All_MinE;
    //VTP 2 Cluster events for events where at least 2 clusters were found (Minimum Cluster E).
    delete h_tSum_2Clust_VTP_800_MinE;
    //Maximum Time difference of 2 Cluster events for events where at least 2 clusters were found.
    delete h_tSum_2Clust_timeDif;
    //Maximum Time difference of 2 Cluster events for events where at least 2 clusters were found for events missed by the clustering trigger.
    delete h_tSum_2Clust_Missed_timeDif;
    //Maximum Cluster Time for 2 Cluster events for events where at least 2 clusters were found and the clustering trigger did not find it.
    delete h_tSum_2Clust_Missed_MaxTime;
    //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where at least 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    delete h_tSum_2Clust_Missed_EdeltaT;
    //Number of Blocks for Missed Events with Delta t < 16 ns where at least 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    delete h_tSum_2Clust_Missed_nBlocks_less;
    //Number of Blocks for Missed Events with Delta t > 16 ns where at least 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    delete h_tSum_2Clust_Missed_nBlocks_greater;

    //Offline 2 Cluster events for events where only 2 clusters were found.
    delete h_tSum_2ClustOnly_All;
    //VTP 2 Cluster events for events where only 2 clusters were found.
    delete h_tSum_2ClustOnly_VTP_800;
    //Offline 2 Cluster events for events where only 2 clusters were found (Minimum Cluster E).
    delete h_tSum_2ClustOnly_All_MinE;
    //VTP 2 Cluster events for events where only 2 clusters were found (Minimum Cluster E).
    delete h_tSum_2ClustOnly_VTP_800_MinE;
    //Time difference of 2 Cluster events for events where only 2 clusters were found.
    delete h_tSum_2ClustOnly_timeDif;
    //Maximum Time difference of 2 Cluster events for events where only 2 clusters were found for events missed by the clustering trigger.
    delete h_tSum_2ClustOnly_Missed_timeDif;
    //Maximum Cluster Time for 2 Cluster events for events where only 2 clusters were found and the clustering trigger did not find it.
    delete h_tSum_2ClustOnly_Missed_MaxTime;
    //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where only 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    delete h_tSum_2ClustOnly_Missed_EdeltaT;
    //Number of Blocks for Missed Events with Delta t < 16 ns where only 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    delete h_tSum_2ClustOnly_Missed_nBlocks_less;
    //Number of Blocks for Missed Events with Delta t > 16 ns where only 2 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    delete h_tSum_2ClustOnly_Missed_nBlocks_greater;

    //Offline 3 Cluster events (max energy combination).
    delete h_tSum_3Clust_All;
    //VTP 3 Cluster events (max energy combination).
    delete h_tSum_3Clust_VTP_800;
    //Offline 3 Cluster events for events where at least 3 clusters were found (Minimum Cluster E).
    delete h_tSum_3Clust_All_MinE;
    //VTP 3 Cluster events for events where at least 3 clusters were found (Minimum Cluster E).
    delete h_tSum_3Clust_VTP_800_MinE;
    //Maximum Time difference of 3 Cluster events for events where at least 3 clusters were found.
    delete h_tSum_3Clust_timeDif;
    //Maximum Time difference of 3 Cluster events for events where at least 3 clusters were found for events missed by the clustering trigger.
    delete h_tSum_3Clust_Missed_timeDif;
    //Maximum Cluster Time for 3 Cluster events for events where at least 3 clusters were found and the clustering trigger did not find it.
    delete h_tSum_3Clust_Missed_MaxTime;
    //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where at least 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    delete h_tSum_3Clust_Missed_EdeltaT;
    //Number of Blocks for Missed Events with Delta t < 16 ns where at least 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    delete h_tSum_3Clust_Missed_nBlocks_less;
    //Number of Blocks for Missed Events with Delta t > 16 ns where at least 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    delete h_tSum_3Clust_Missed_nBlocks_greater;
    //Individual cluster energy for at least 3 or more cluster events
    delete h_tSum_3Clust_Cl_Energy_All;
    delete h_tSum_3Clust_Cl_Energy_Trig;
    //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
    delete h_tSum_3Clust_GEMMatch_All;
    delete h_tSum_3Clust_GEMMatch_Trig;
    //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
    delete h_tSum_3Clust_GEMMatch_All_2Match_Either;
    delete h_tSum_3Clust_GEMMatch_Trig_2Match_Either;
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    delete h_tSum_3Clust_GEMMatch_All_3Match_Either;
    delete h_tSum_3Clust_GEMMatch_Trig_3Match_Either;
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    delete h_tSum_3Clust_GEMMatch_All_3Match_Both;
    delete h_tSum_3Clust_GEMMatch_Trig_3Match_Both;
    //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
    delete h_tSum_3Clust_GEMMatch_All_Cl_E;
    delete h_tSum_3Clust_GEMMatch_Trig_Cl_E;
    //Cluster sum energy for at least 3 or more cluster events with at least 2 GEM matches.
    delete h_tSum_3Clust_GEMMatch_All_2Match_Either_Cl_E;
    delete h_tSum_3Clust_GEMMatch_Trig_2Match_Either_Cl_E;
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    delete h_tSum_3Clust_GEMMatch_All_3Match_Either_Cl_E;
    delete h_tSum_3Clust_GEMMatch_Trig_3Match_Either_Cl_E;
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches.
    delete h_tSum_3Clust_GEMMatch_All_3Match_Both_Cl_E;
    delete h_tSum_3Clust_GEMMatch_Trig_3Match_Both_Cl_E;
    //Cluster sum energy for at least 3 or more cluster events with at least 3 GEM matches and a high/low cluster energy sum.
    delete h_tSum_3Clust_vZ_High;
    delete h_tSum_3Clust_vZ_Low;
    //Maximum cluster energy for only 3 cluster events
    delete h_tSum_3Clust_Cl_Max_All;
    delete h_tSum_3Clust_Cl_Max_Trig;


    //Offline 3 Cluster events for events where only 3 clusters were found.
    delete h_tSum_3ClustOnly_All;
    //VTP 3 Cluster events for events where only 3 clusters were found.
    delete h_tSum_3ClustOnly_VTP_800;
    //Offline 3 Cluster events for events where only 3 clusters were found (Minimum Cluster E).
    delete h_tSum_3ClustOnly_All_MinE;
    //VTP 3 Cluster events for events where only 3 clusters were found (Minimum Cluster E).
    delete h_tSum_3ClustOnly_VTP_800_MinE;
    //Time difference of 3 Cluster events for events where only 3 clusters were found.
    delete h_tSum_3ClustOnly_timeDif;
    //Maximum Time difference of 3 Cluster events for events where only 3 clusters were found for events missed by the clustering trigger.
    delete h_tSum_3ClustOnly_Missed_timeDif;
    //Maximum Cluster Time for 3 Cluster events for events where only 3 clusters were found and the clustering trigger did not find it.
    delete h_tSum_3ClustOnly_Missed_MaxTime;
    //Minimum Cluster Energy for Missed Events with Delta t < 16 ns where only 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    delete h_tSum_3ClustOnly_Missed_EdeltaT;
    //Number of Blocks for Missed Events with Delta t < 16 ns where at only 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    delete h_tSum_3ClustOnly_Missed_nBlocks_less;
    //Number of Blocks for Missed Events with Delta t > 16 ns where only 3 clusters were found, the event was missed by the clustering trigger, and the delta T was less than 16ns.
    delete h_tSum_3ClustOnly_Missed_nBlocks_greater;
    //Individual cluster energy for only 3 cluster events
    delete h_tSum_3ClustOnly_Cl_Energy_All;
    delete h_tSum_3ClustOnly_Cl_Energy_Trig;
    //Maximum cluster energy for only 3 cluster events
    delete h_tSum_3ClustOnly_Cl_Max_All;
    delete h_tSum_3ClustOnly_Cl_Max_Trig;
    //
    delete h_tSum_3Clust_DCA_Low;
    delete h_tSum_3Clust_delta_X_Low;
    delete h_tSum_3Clust_DCA_High;
    delete h_tSum_3Clust_delta_X_High;

    //Offline 2 Cluster events; HyCal XY Positions
    delete h_tSum_2Clust_All_HC_XY;
    //VTP 2 Cluster events; HyCal XY Positions
    delete h_tSum_2Clust_VTP_HC_XY;
    //2 Cluster events missed by the VTP Trigger; HyCal XY, and Theta Positions
    delete h_tSum_2Clust_Miss_HC_XY;
    delete h_tSum_2Clust_Miss_HC_Theta;
    delete h_tSum_2Clust_Miss_HC_XY_LowE;

    //Number of VTP Clusters
    delete h_tSum_VTP_NClustDiff;
    //Number of VTP Clusters - Offline Clusters for rejected events.
    delete h_tSum_VTP_NClustDiff_Missed;
    //VTP timing of events with non-zero difference in cluster times
    delete h_tSum_VTP_Timing;
        
    //VTP center cluster module ID
    delete h_tSum_VTP_centerID;
    //Offline center cluster module ID
    delete h_tSum_Off_centerID;
        

    //Offline 3 Cluster events; HyCal XY Positions
    delete h_tSum_3Clust_All_HC_XY;
    //VTP 3 Cluster events; HyCal XY Positions
    delete h_tSum_3Clust_VTP_HC_XY;
    //3 Cluster events missed by the VTP Trigger; HyCal XY, and Theta Positions
    delete h_tSum_3Clust_Miss_HC_XY;
    delete h_tSum_3Clust_Miss_HC_Theta;
    delete h_tSum_3Clust_Miss_HC_XY_LowE;
}

/**
 * Fills all the random trigger event histograms.
 *
 * @param i - the current event number
 */
void ClustTrig::fill_rand_Histos_wThr(Int_t i){

    randInformation passed_3Clust = ComputeTimeBinnedClusters(true);
    array<double,Trigger::nSamples> cl_sums = passed_3Clust.doubleInfoVecs.at(0);
    array<double,Trigger::nSamples> cl_min = passed_3Clust.doubleInfoVecs.at(1);
    array<double,Trigger::nSamples> cl_max = passed_3Clust.doubleInfoVecs.at(2);
    array<double,Trigger::nSamples> cl_earlyTime = passed_3Clust.doubleInfoVecs.at(3);
    array<double,Trigger::nSamples> cl_lateTime = passed_3Clust.doubleInfoVecs.at(4);

    array<int,Trigger::nSamples> cl_count = passed_3Clust.integerInfoVecs.at(0);
    array<int,Trigger::nSamples> cl_minInd = passed_3Clust.integerInfoVecs.at(1);
    array<int,Trigger::nSamples> cl_maxInd = passed_3Clust.integerInfoVecs.at(2);
    array<int,Trigger::nSamples> cl_earlyInd = passed_3Clust.integerInfoVecs.at(3);
    array<int,Trigger::nSamples> cl_lateInd = passed_3Clust.integerInfoVecs.at(4);

    array<vector<int>, Trigger::nSamples> cl_ind =  passed_3Clust.allIndexVectors;

    Int_t ind_MaxEsum = -1;

    for(Int_t j = 10; j < Trigger::nSamples-10; j++){
        if(cl_count[j] < 3){continue;}

        if(ind_MaxEsum > -1 && cl_sums[j] > cl_sums[ind_MaxEsum]){
            ind_MaxEsum = j;
        }
    }

    if(ind_MaxEsum > -1){
        h_rand_3Clust_All_Sum->Fill(cl_sums[ind_MaxEsum]);
        h_rand_3Clust_All_Min->Fill(cl_min[ind_MaxEsum]);
        h_rand_3Clust_All_Max->Fill(cl_max[ind_MaxEsum]);
        Float_t dt = cl_time[cl_lateInd[ind_MaxEsum]] - cl_time[cl_earlyInd[ind_MaxEsum]];
        h_rand_3Clust_deltaT_All->Fill(dt);

        for(UInt_t k = 0; k < cl_ind[ind_MaxEsum].size(); k++){
            int index = cl_ind[ind_MaxEsum].at(k);
            h_rand_3Clust_XY_All->Fill(cl_x[index], cl_y[index]);
        }

        for(Int_t b = 0; b < nSSPBits - 1; b++){
            if(!trig.GetTriggersBit(b)->empty()){
                h_rand_3Clust_VTP_Sum[b]->Fill(cl_sums[ind_MaxEsum]);
                h_rand_3Clust_VTP_Min[b]->Fill(cl_min[ind_MaxEsum]);
                h_rand_3Clust_VTP_Max[b]->Fill(cl_max[ind_MaxEsum]);

                h_rand_3Clust_deltaT_VTP[b]->Fill(dt);

                const vector<int>* times = trig.GetTriggersBit(CL3_FLAG-nSSPBits);
                for (int t : *times){
                    h_rand_3Clust_timeCorrEarly[b]->Fill(cl_earlyTime[ind_MaxEsum],t);
                    h_rand_3Clust_timeCorrLate[b]->Fill(cl_lateTime[ind_MaxEsum],t);
                }

                for(UInt_t l = 0; l < cl_ind[ind_MaxEsum].size(); l++){
                    int index2 = cl_ind[ind_MaxEsum].at(l);
                    h_rand_3Clust_XY_VTP[b]->Fill(cl_x[index2], cl_y[index2]);
                }
            }
        }
    }

    for(Int_t a = 0; a < vtp_cl_n; a++){
        h_tSum_VTP_Timing->Fill(vtp_cl_time[a]);
        h_tSum_VTP_centerID->Fill(vtp_cl_center[a]);  
    }
}

void ClustTrig::setup_rand_Histos(){

    for (int b = 0; b < nSSPBits; ++b) {
        TString sumName  = TString::Format("h_rand_3Clust_VTP_Sum_bit%d", b);
        TString sumTitle = TString::Format( "Random Trigger Event with at least 3 Clusters Found by E_{Sum} - bit %d;E_{Sum};Counts", b);
        h_rand_3Clust_VTP_Sum[b] = new TH1F(sumName, sumTitle, en+500,0, en+500);

        TString minName  = TString::Format("h_rand_3Clust_VTP_Min_bit%d", b);
        TString minTitle = TString::Format( "Random Trigger Event with at least 3 Clusters Found by E_{Min} - bit %d;E_{Min};Counts", b);
        h_rand_3Clust_VTP_Min[b] = new TH1F(minName, minTitle, en+500,0, en+500);

        TString maxName  = TString::Format("h_rand_3Clust_VTP_Max_bit%d", b);
        TString maxTitle = TString::Format( "Random Trigger Event with at least 3 Clusters Found by E_{Max} - bit %d;E_{Min};Counts", b);
        h_rand_3Clust_VTP_Max[b] = new TH1F(maxName, maxTitle, en+500,0, en+500);

        TString timeEarlyName  = TString::Format("h_rand_3Clust_timeCorrEarly_bit%d", b);
        TString timeEarlyTitle = TString::Format( "Time of Trigger Firing v. Time of Sum Threshold Crossing by earliest cluster - bit %d;Clust Time [4ns], Trigger Time [4ns]", b);
        h_rand_3Clust_timeCorrEarly[b] = new TH2F(timeEarlyName, timeEarlyTitle, kTrigTimeBins, kTrigTimeMin, kTrigTimeMax, kTrigTimeBins, kTrigTimeMin, kTrigTimeMax);

        TString timeLateName  = TString::Format("h_rand_3Clust_timeCorrLate_bit%d", b);
        TString timeLateTitle = TString::Format( "Time of Trigger Firing v. Time of Sum Threshold Crossing by latest cluster - bit %d;Clust Time [4ns], Trigger Time [4ns]", b);
        h_rand_3Clust_timeCorrLate[b] = new TH2F(timeLateName, timeLateTitle, kTrigTimeBins, kTrigTimeMin, kTrigTimeMax, kTrigTimeBins, kTrigTimeMin, kTrigTimeMax);

        TString XY_Name  = TString::Format("h_rand_3Clust_XY_VTP_bit%d", b);
        TString XY_Title = TString::Format( "Offline Clusters Found when VTP did fired - bit %d;x [mm];y [mm]", b);
        h_rand_3Clust_XY_VTP[b] = new TH2F(XY_Name, XY_Title, 34, Xedge, 34, Yedge);

        TString dtName = TString::Format("h_rand_3Clust_deltaT_VTP_bit%d", b);
        TString dtTitle = TString::Format("#Deltat (Max_{cl_time} - Min_{cl_time}) - bit %d", b);
        h_rand_3Clust_deltaT_VTP[nSSPBits] = new TH1F(dtName, dtTitle, 32, 0, 32);
        
    }

    h_rand_3Clust_All_Sum = new TH1F("h_rand_3Clust_All_Sum", "Random Trigger Event with at least 3 Clusters Found by E_{Sum};E_{Sum};Counts",en+500,0, en+500);
    h_rand_3Clust_All_Min = new TH1F("h_rand_3Clust_All_Min", "Random Trigger Event with at least 3 Clusters Found by E_{Min};E_{Min};Counts",en+500,0, en+500);
    h_rand_3Clust_All_Max = new TH1F("h_rand_3Clust_All_Max", "Random Trigger Event with at least 3 Clusters Found by E_{Max};E_{Max};Counts",en+500,0, en+500);
    h_rand_3Clust_XY_All = new TH2F("h_rand_3Clust_XY_All", "Offline Clusters Found when VTP should have fired;x [mm];y [mm]",34, Xedge, 34, Yedge);
    h_rand_3Clust_deltaT_All = new TH1F("h_rand_3Clust_deltaT_All", "#Deltat (Max_{cl_time} - Min_{cl_time}) All", 32, 0, 32);
}

/**
 * Forms a vector of each time sample window and the elligible energy sums within each time window.
 *
 * @param cl_mult - cluster multiplicity to search for good events in
 * @param appThr - bool to decide if the blocks and individual energy threshold should be applied.
 */
ClustTrig::randInformation ClustTrig::ComputeTimeBinnedClusters(bool appThr){
    /*array<double, nSamples> sums{};
    if (n_ch <= 0) return sums;

    for (int c = 0; c < n_ch; ++c) {
        for (int p = 0; p < nMaxPeaksPerCh; ++p) {
            const double e = double(peak_integral[c][p]);
            if (e == 0.0) continue;          // skip empty slots
            const double t = peak_time[c][p];

            // i contributes when i*step <= t < i*step + window
            //   ⇒  (t - window)/step  <  i  <=  t/step
            int iMin = static_cast<int>(std::floor((t - kSampleWindow) / kSampleStep)) + 1;
            int iMax = static_cast<int>(std::floor(t / kSampleStep));
            if (iMin < 0)         iMin = 0;
            if (iMax >= nSamples) iMax = nSamples - 1;

            for (int i = iMin; i <= iMax; ++i)
                sums[i] += e;
        }
    }
    return sums;*/

    array<double, Trigger::nSamples> cl_sums{};
    array<double, Trigger::nSamples> cl_min{};
    array<double, Trigger::nSamples> cl_max{};
    array<double, Trigger::nSamples> cl_earlyTime{};
    array<double, Trigger::nSamples> cl_lateTime{};

    array<int, Trigger::nSamples> cl_count{};
    array<int, Trigger::nSamples> cl_minInd{};
    array<int, Trigger::nSamples> cl_maxInd{};
    array<int, Trigger::nSamples> cl_earlyInd{};
    array<int, Trigger::nSamples> cl_lateInd{};

    array<vector<int>, Trigger::nSamples> cl_ind{};

    bool first = true;
    for(int c = 0; c < nClust; c++){
        double e = cl_E[c];
        if(appThr && (cl_E[c] < CL_IND_THR || cl_nblocks[c] < 2 || (usemax && cl_E[c] > CL_IND_MAX_THR) || binary_search(excludeMod.begin(), excludeMod.end(), cl_center[c]))) continue; 
        double t = cl_time[c];

        int iMin = static_cast<int>(std::floor((t - kSampleWindow) / kSampleStep)) + 1;
        int iMax = static_cast<int>(std::floor(t / kSampleStep));
        if (iMin < 0){
            iMin = 0;
        }

        for (int i = iMin; i <= iMax; ++i){
            cl_sums[i] += e;
            cl_ind[i].push_back(c);
            if(first){
                first = false;
                cl_min[i] = e;
                cl_minInd[i] = c;

                cl_max[i] = e;
                cl_maxInd[i] = c;
                
                cl_earlyTime[i] = t;
                cl_earlyInd[i] = c;
                
                cl_lateTime[i] = t;
                cl_lateInd[i] = c;

            }
            else{
                if(cl_min[i] > e){
                    cl_min[i] = e;
                    cl_minInd[i] = c;
                }
                if(cl_max[i] < e){
                    cl_max[i] = e;
                    cl_maxInd[i] = c;
                }
                if(cl_earlyTime[i] > t){
                    cl_earlyTime[i] = t;
                    cl_earlyInd[i] = c;
                }
                if(cl_lateTime[i] < t){
                    cl_lateTime[i] = t;
                    cl_lateInd[i] = c;
                }
            }
            cl_count[i] += 1;
        }

    }

    vector<array<double, Trigger::nSamples>> retDoubleVec;
    retDoubleVec.push_back(cl_sums);
    retDoubleVec.push_back(cl_min);
    retDoubleVec.push_back(cl_max);
    retDoubleVec.push_back(cl_earlyTime);
    retDoubleVec.push_back(cl_lateTime);

    vector<array<int, Trigger::nSamples>> retIntegerVec;
    retIntegerVec.push_back(cl_count);
    retIntegerVec.push_back(cl_minInd);
    retIntegerVec.push_back(cl_maxInd);
    retIntegerVec.push_back(cl_earlyInd);
    retIntegerVec.push_back(cl_lateInd);

    randInformation infoBySampleWindow;
    infoBySampleWindow.doubleInfoVecs = retDoubleVec;
    infoBySampleWindow.integerInfoVecs = retIntegerVec;
    infoBySampleWindow.allIndexVectors = cl_ind;

    return infoBySampleWindow;
}

/**
 * Saves all the relevant random trigger histograms.
 *
 * @param rootfileName - the name of the ROOT file to save these histograms too
 */
void ClustTrig::SaveAllRandHistograms(TString rootfileName){
    TObjArray* arr = new TObjArray(0,0);

    (*arr).Add(h_rand_3Clust_All_Sum);
    (*arr).Add(h_rand_3Clust_All_Max);
    (*arr).Add(h_rand_3Clust_All_Min);
    (*arr).Add(h_rand_3Clust_XY_All);
    (*arr).Add(h_rand_3Clust_deltaT_All);

    for(Int_t b = 0; b < nSSPBits; b++){
        (*arr).Add(h_rand_3Clust_VTP_Sum[b]);
        (*arr).Add(h_rand_3Clust_VTP_Max[b]);
        (*arr).Add(h_rand_3Clust_VTP_Min[b]);
        (*arr).Add(h_rand_3Clust_timeCorrEarly[b]);
        (*arr).Add(h_rand_3Clust_timeCorrLate[b]);
        (*arr).Add(h_rand_3Clust_XY_VTP[b]);
        (*arr).Add(h_rand_3Clust_deltaT_VTP[b]);

        (*arr).Add(h_rand_3Clust_Sum_Ratio[b]);
        (*arr).Add(h_rand_3Clust_Max_Ratio[b]);
        (*arr).Add(h_rand_3Clust_Min_Ratio[b]);
        (*arr).Add(h_rand_3Clust_deltaT_Ratio[b]);
    }

    TFile file1(rootfileName,"RECREATE");
    (*arr).Write();
	file1.Close();
}

/**
 * Deletes the random trigger histograms.
 */
void ClustTrig::delete_rand_Histos(){
    delete h_rand_3Clust_All_Sum;
    delete h_rand_3Clust_All_Max;
    delete h_rand_3Clust_All_Min;
    delete h_rand_3Clust_XY_All;
    delete h_rand_3Clust_deltaT_All;

    for(Int_t b = 0 ; b < nSSPBits; b++){
        delete h_rand_3Clust_VTP_Sum[b];
        delete h_rand_3Clust_VTP_Max[b];
        delete h_rand_3Clust_VTP_Min[b];
        delete h_rand_3Clust_deltaT_VTP[b];
        delete h_rand_3Clust_timeCorrEarly[b];
        delete h_rand_3Clust_timeCorrLate[b];
        delete h_rand_3Clust_XY_VTP[b];

        delete h_rand_3Clust_Sum_Ratio[b];
        delete h_rand_3Clust_Max_Ratio[b];
        delete h_rand_3Clust_Min_Ratio[b];
        delete h_rand_3Clust_deltaT_Ratio[b];
    }
}

/**
 * Prints out the trigger performance plots from the Random trigger.
 *
 * @param pdfName - the name of the PDF to print to
 */
void ClustTrig::printRandPDF(TString pdfName){

    TString Bit_Thr[nSSPBits-1] = {"650 MeV","750 MeV", "950 MeV", "1300 MeV", "1100 MeV", "1400 MeV", "1500 MeV"};
    Color_t Bit_Color[nSSPBits-1] = {kBlue, kRed, kGreen,kCyan, kMagenta, kSpring, kPink};

    TCanvas *c = new TCanvas("c", "Rand_ClustTrig_Canvas",1000,1000);
    auto legend = new TLegend(0.1,0.8,0.4,0.9);

    //Page 1
    //-------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_rand_3Clust_All_Sum->SetLineColor(kBlack);
    h_rand_3Clust_All_Sum->Draw("HIST");

    legend->SetHeader("Legend","C"); // option "C" allows to center the header
    legend->AddEntry(h_rand_3Clust_All_Sum,"All","l");
    for(Int_t b = 0; b < nSSPBits - 1; b++){
        h_rand_3Clust_VTP_Sum[b]->SetLineColor(Bit_Color[b]);
        h_rand_3Clust_VTP_Sum[b]->Draw("SAME HIST");

        legend->AddEntry(h_rand_3Clust_VTP_Sum[b],Bit_Thr[b], "l");
    }
	legend->Draw();

    c->cd(2);
    for(Int_t m = 0; m < nSSPBits - 1; m++){
        h_rand_3Clust_Sum_Ratio[m]->SetLineColor(Bit_Color[m]);
        if(m == 0){
            h_rand_3Clust_Sum_Ratio[m]->Draw("HIST");
        }
        else{
            h_rand_3Clust_Sum_Ratio[m]->Draw("HIST SAME");
        }
    }
    

    c->Print(pdfName + "(");
    c->Clear();
    legend->Clear();

    //Page 2
    //-------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_rand_3Clust_All_Sum->Draw("HIST");

    legend->SetHeader("Legend","C"); // option "C" allows to center the header
    legend->AddEntry(h_rand_3Clust_All_Sum,"All","l");
    for(Int_t l = 0; l < nSSPBits - 1; l++){
        h_rand_3Clust_VTP_Sum[l]->Draw("SAME HIST");
        legend->AddEntry(h_rand_3Clust_VTP_Sum[l],Bit_Thr[l], "l");
    }
	legend->Draw();

    c->cd(2);
    for(Int_t q = 0; q < nSSPBits - 1; q++){
        h_rand_3Clust_Sum_Ratio[q]->SetAxisRange(0.9,1.0,"Y");
        if(q == 0){
            h_rand_3Clust_Sum_Ratio[q]->Draw("HIST");
        }
        else{
            h_rand_3Clust_Sum_Ratio[q]->Draw("HIST SAME");
        }
    }

    c->Print(pdfName + "(");
    c->Clear();
    legend->Clear();

    //Page 3
    //-------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_rand_3Clust_All_Min->SetLineColor(kBlack);
    h_rand_3Clust_All_Min->Draw("HIST");

    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_rand_3Clust_All_Min,"All","l");
    for(Int_t n = 0; n < nSSPBits - 1; n++){
        h_rand_3Clust_VTP_Min[n]->SetLineColor(Bit_Color[n]);
        h_rand_3Clust_VTP_Min[n]->Draw("SAME HIST");
        legend->AddEntry(h_rand_3Clust_VTP_Min[n],Bit_Thr[n], "l");
    }
	legend->Draw();

    c->cd(2);
    for(Int_t p = 0; p < nSSPBits - 1; p++){
        h_rand_3Clust_Min_Ratio[p]->SetLineColor(Bit_Color[p]);

        if(p == 0){
            h_rand_3Clust_Min_Ratio[p]->Draw("HIST");
        }
        else{
            h_rand_3Clust_Min_Ratio[p]->Draw("HIST SAME");
        }
    }

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 4
    //-------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_rand_3Clust_All_Min->Draw("HIST");

    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_rand_3Clust_All_Min,"All","l");
    for(Int_t w = 0; w < nSSPBits - 1; w++){
        h_rand_3Clust_VTP_Min[w]->Draw("SAME HIST");
        legend->AddEntry(h_rand_3Clust_VTP_Min[w],Bit_Thr[w], "l");
    }
	legend->Draw();

    c->cd(2);
    h_rand_3Clust_Min_Ratio[0]->SetAxisRange(0.9,1.0,"Y");

    for(Int_t r = 0; r < nSSPBits - 1; r++){
        if(r == 0){
            h_rand_3Clust_Min_Ratio[r]->Draw("HIST");
        }
        else{
            h_rand_3Clust_Min_Ratio[r]->Draw("HIST SAME");
        }
    }

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 5
    //-------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_rand_3Clust_All_Max->SetLineColor(kBlack);
    h_rand_3Clust_All_Max->Draw("HIST");
    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_rand_3Clust_All_Max,"All","l");
    for(Int_t t = 0; t < nSSPBits - 1; t++){
        h_rand_3Clust_VTP_Max[t]->SetLineColor(Bit_Color[t]);
        h_rand_3Clust_VTP_Max[t]->Draw("SAME HIST");
        legend->AddEntry(h_rand_3Clust_VTP_Max[t],Bit_Thr[t], "l");
    }
	legend->Draw();

    c->cd(2);
    for(Int_t y = 0; y < nSSPBits; y++){
        h_rand_3Clust_Max_Ratio[y]->SetLineColor(Bit_Color[y]);
        if(y == 0){
            h_rand_3Clust_Max_Ratio[y]->Draw("HIST");       
        }
        else{
            h_rand_3Clust_Max_Ratio[y]->Draw("HIST SAME");
        }
    }

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 6
    //-------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_rand_3Clust_All_Max->Draw("HIST");

    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_rand_3Clust_All_Max,"All","l");
    for(Int_t u = 0; u < nSSPBits - 1; u++){
        h_rand_3Clust_VTP_Max[u]->Draw("SAME HIST");
        legend->AddEntry(h_rand_3Clust_VTP_Max[u],Bit_Thr[u], "l");
    }
	legend->Draw();

    c->cd(2);
    h_rand_3Clust_Max_Ratio[0]->SetAxisRange(0.9,1.0,"Y");

    for(Int_t s = 0; s < nSSPBits - 1; s++){
        if(s == 0){
            h_rand_3Clust_Max_Ratio[s]->Draw("HIST");
        }
        else{
            h_rand_3Clust_Max_Ratio[s]->Draw("HIST SAME");
        }
    }

    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 7
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_rand_3Clust_deltaT_All->SetLineColor(kBlack);
    h_rand_3Clust_deltaT_All->Draw("HIST");

    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_rand_3Clust_deltaT_All,"All","l");
    for(Int_t g = 0; g < nSSPBits - 1; g++){
        h_rand_3Clust_deltaT_VTP[g]->Draw("SAME HIST");
        legend->AddEntry(h_rand_3Clust_deltaT_VTP[g],Bit_Thr[g], "l");
    }
	legend->Draw();

    c->cd(2);
    h_rand_3Clust_deltaT_Ratio[0]->SetAxisRange(0.9,1.0,"Y");
    for(Int_t h = 0; h < nSSPBits - 1; h++){
        if(h == 0){
            h_rand_3Clust_deltaT_Ratio[h]->Draw("HIST");
        }
        else{
            h_rand_3Clust_deltaT_Ratio[h]->Draw("HIST SAME");
        }
    }
    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 7
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_rand_3Clust_deltaT_All->SetLineColor(kBlack);
    h_rand_3Clust_deltaT_All->Draw("HIST");

    legend->SetHeader("Legend","C"); // option "C" allows to center the header
	legend->AddEntry(h_rand_3Clust_deltaT_All,"All","l");
    for(Int_t x = 0; x < nSSPBits - 1; x++){
        h_rand_3Clust_deltaT_VTP[x]->Draw("SAME HIST");
        legend->AddEntry(h_rand_3Clust_deltaT_VTP[x],Bit_Thr[x], "l");
    }
	legend->Draw();

    c->cd(2);
    for(Int_t h = 0; h < nSSPBits - 1; h++){
        h_rand_3Clust_deltaT_Ratio[h]->SetLineColor(Bit_Color[h]);
        if(h == 0){
            h_rand_3Clust_deltaT_Ratio[h]->Draw("HIST");
        }
        else{
            h_rand_3Clust_deltaT_Ratio[h]->Draw("HIST SAME");
        }
    }
    c->Print(pdfName);
    c->Clear();
    legend->Clear();

    //Page 8
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    h_rand_3Clust_timeCorrEarly[3]->Draw("COLZ");

    c->cd(2);
    h_rand_3Clust_timeCorrLate[3]->Draw("COLZ");

    c->Print(pdfName);
    c->Clear();

    //Page 8
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    gPad->SetLogz(1);
    h_rand_3Clust_XY_All->Draw("COLZ");

    c->cd(2);
    gPad->SetLogz(1);
    h_rand_3Clust_XY_VTP[0]->Draw("COLZ");

    c->cd(3);
    gPad->SetLogz(1);
    h_rand_3Clust_XY_VTP[1]->Draw("COLZ");

    c->cd(4);
    gPad->SetLogz(1);
    h_rand_3Clust_XY_VTP[2]->Draw("COLZ");

    c->Print(pdfName);
    c->Clear();

    //Page 9
    //-----------------------------------------------------------------------------------
    c->Divide(2,2);
    c->cd(1);
    gPad->SetLogz(1);
    h_rand_3Clust_XY_VTP[3]->Draw("COLZ");

    c->cd(2);
    gPad->SetLogz(1);
    h_rand_3Clust_XY_VTP[4]->Draw("COLZ");

    c->cd(3);
    gPad->SetLogz(1);
    h_rand_3Clust_XY_VTP[5]->Draw("COLZ");

    c->cd(4);
    gPad->SetLogz(1);
    h_rand_3Clust_XY_VTP[6]->Draw("COLZ");

    c->Print(pdfName);
    c->Clear();

    //Page 9
    //-----------------------------------------------------------------------------------
    c->Divide(2,4);
    for(Int_t b = 0; b < nSSPBits; b++){
        c->cd(b+1);
        hTrigTime[b]->Draw("HIST");
        hTrigTime[b]->SetStats(0);
    }
    c->Print(pdfName + ")");
    c->Clear();
}

/**
 * Sets the value that determines whether the maximum energy threshold is used or not.
 * 
 * @param set - true if this object should use the upper cluster limit or not for evaluation.
 */
void ClustTrig::setUseMax(bool set){
    usemax = set;
}