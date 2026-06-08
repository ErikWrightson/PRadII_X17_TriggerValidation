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
 * Creates the LMS trigger object.
 *
 * @param c - the TChain of events to
 */
ClustTrig::ClustTrig(TChain* c):Trigger(c, false, true){
    cout<<"Set All Branch Addresses for Clustering Trigger Object.\n";

    en = 3485.41;

    setup_Histos();
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

        if(self && (trigger_bits & (1<<CL1_FLAG))){

        }
        if(rand && (trigger_bits & (1<<RAND_FC_FLAG))){

        }
        if(tSum && (trigger_bits & (1<<TSUM_FLAG)) && !(trigger_bits & (1<<LMS_FLAG))){ //(trigType & (1<<SSP_RAWSUM_TFLAG))

            Double_t clustSum = 0;
            for(Int_t j = 0; j < nClust; j++){
                if(cl_E[j] > CL_IND_THR){
                    clustSum += cl_E[j];
                }
            }
            
            if(nClust >= 1){

                if(nClust == 1){
                    if(trigger_bits & (1<<CL1_FLAG)){
                        h_tSum_1ClustOnly_VTP_800->Fill(clustSum);
                    }
                    if(trigger_bits & (1<<CL1_500_FLAG)){
                        h_tSum_1ClustOnly_VTP_500->Fill(clustSum);
                    }
                    if(trigger_bits & (1<<CL1_600_FLAG)){
                        h_tSum_1ClustOnly_VTP_600->Fill(clustSum);
                    }
                    if(trigger_bits & (1<<CL1_900_FLAG)){
                        h_tSum_1ClustOnly_VTP_900->Fill(clustSum);
                    }
                    h_tSum_1ClustOnly_All->Fill(clustSum);
                }

               if(trigger_bits & (1<<CL1_FLAG)){
                    h_tSum_1Clust_VTP_800->Fill(clustSum);
               }
               if(trigger_bits & (1<<CL1_500_FLAG)){
                    h_tSum_1Clust_VTP_500->Fill(clustSum);
               }
               if(trigger_bits & (1<<CL1_600_FLAG)){
                    h_tSum_1Clust_VTP_600->Fill(clustSum);
               }
               if(trigger_bits & (1<<CL1_900_FLAG)){
                    h_tSum_1Clust_VTP_900->Fill(clustSum);
               }
               h_tSum_1Clust_All->Fill(clustSum);
               
               if(nClust >= 2){

                    if(nClust == 2){
                        if(trigger_bits & (1<<CL2_FLAG)){
                            h_tSum_2ClustOnly_VTP_800->Fill(clustSum);
                        }
                        h_tSum_2ClustOnly_All->Fill(clustSum);
                    }

                    if(trigger_bits & (1<<CL2_FLAG)){
                        h_tSum_2Clust_VTP_800->Fill(clustSum);
                    }
                    h_tSum_2Clust_All->Fill(clustSum);

                    if(nClust >= 3){

                        if(nClust == 3){
                            if(trigger_bits & (1<<CL3_FLAG)){
                                h_tSum_3ClustOnly_VTP_800->Fill(clustSum);
                            }
                            h_tSum_3ClustOnly_All->Fill(clustSum);
                        }

                        if(trigger_bits & (1<<CL3_FLAG)){
                            h_tSum_3Clust_VTP_800->Fill(clustSum);
                        }
                        h_tSum_3Clust_All->Fill(clustSum);
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
    h_tSum_1Clust_All = new TH1F("h_tSum_1Clust_All", "Total Sum Event with at least 1 Cluster Found, Maximum Offline Energy;E_{Sum};Counts",en+1000,0, en+1000);
    //VTP 1 Cluster Trigger (max energy cluster).
    h_tSum_1Clust_VTP_800 = new TH1F("h_tSum_1Clust_VTP_800", "Total Sum Event with at least 1 Cluster Found on the VTP with 800 MeV threshold;E_{Sum};Counts",en+1000,0, en+1000);
    h_tSum_1Clust_VTP_900 = new TH1F("h_tSum_1Clust_VTP_900", "Total Sum Event with at least 1 Cluster Found on the VTP with 900 MeV threshold;E_{Sum};Counts",en+1000,0, en+1000);
    h_tSum_1Clust_VTP_600 = new TH1F("h_tSum_1Clust_VTP_600", "Total Sum Event with at least 1 Cluster Found on the VTP with 600 MeV threshold;E_{Sum};Counts",en+1000,0, en+1000);
    h_tSum_1Clust_VTP_500 = new TH1F("h_tSum_1Clust_VTP_500", "Total Sum Event with at least 1 Cluster Found on the VTP with 500 MeV threshold;E_{Sum};Counts",en+1000,0, en+1000);

    //Offline 1 Cluster Trigger for events where only 1 cluster was found.
    h_tSum_1ClustOnly_All = new TH1F("h_tSum_1ClustOnly_All", "Total Sum Event with ONLY 1 Cluster Found, Maximum Offline Energy;E_{Sum};Counts",en+1000,0, en+1000);
    //VTP 1 Cluster Trigger for events where only 1 cluster was found.
    h_tSum_1ClustOnly_VTP_800 = new TH1F("h_tSum_1ClustOnly_VTP_800", "Total Sum Event with ONLY 1 Cluster Found on the VTP with 800 MeV threshold;E_{Sum};Counts",en+1000,0, en+1000);
    h_tSum_1ClustOnly_VTP_900 = new TH1F("h_tSum_1ClustOnly_VTP_900", "Total Sum Event with at ONLY 1 Cluster Found on the VTP with 900 MeV threshold;E_{Sum};Counts",en+1000,0, en+1000);
    h_tSum_1ClustOnly_VTP_600 = new TH1F("h_tSum_1ClustOnly_VTP_600", "Total Sum Event with at ONLY 1 Cluster Found on the VTP with 600 MeV threshold;E_{Sum};Counts",en+1000,0, en+1000);
    h_tSum_1ClustOnly_VTP_500 = new TH1F("h_tSum_1ClustOnly_VTP_500", "Total Sum Event with at ONLY 1 Cluster Found on the VTP with 500 MeV threshold;E_{Sum};Counts",en+1000,0, en+1000);

    //Offline 2 Cluster events (max energy combination).
    h_tSum_2Clust_All = new TH1F("h_tSum_2Clust_All", "Total Sum Event with at least 2 Clusters Found, Maximum Offline Combination Energy;E_{Sum};Counts",en+1000,0, en+1000);
    //VTP 2 Cluster events (max energy combination).
    h_tSum_2Clust_VTP_800 = new TH1F("h_tSum_2Clust_VTP_800", "Total Sum Event with at least 2 Clusters Found on the VTP with 800 MeV sum threshold;E_{Sum};Counts",en+1000,0, en+1000);

    //Offline 2 Cluster events for events where only 2 clusters were found.
    h_tSum_2ClustOnly_All = new TH1F("h_tSum_2ClustOnly_All", "Total Sum Event with ONLY 2 Clusters Found, Maximum Offline Combination Energy;E_{Sum};Counts",en+1000,0, en+1000);
    //VTP 2 Cluster events for events where only 2 clusters were found.
    h_tSum_2ClustOnly_VTP_800 = new TH1F("h_tSum_2ClustOnly_VTP_800", "Total Sum Event with ONLY 2 Clusters Found on the VTP with 800 MeV sum threshold;E_{Sum};Counts",en+1000,0, en+1000);

    //Offline 3 Cluster events (max energy combination).
    h_tSum_3Clust_All = new TH1F("h_tSum_3Clust_All", "Total Sum Event with at least 3 Clusters Found, Maximum Offline Combination Energy;E_{Sum};Counts",en+1000,0, en+1000);
    //VTP 3 Cluster events (max energy combination).
    h_tSum_3Clust_VTP_800 = new TH1F("h_tSum_3Clust_VTP_800", "Total Sum Event with at least 3 Clusters Found on the VTP with 800 MeV sum threshold;E_{Sum};Counts",en+1000,0, en+1000);

    //Offline 3 Cluster events for events where only 3 clusters were found.
    h_tSum_3ClustOnly_All = new TH1F("h_tSum_3ClustOnly_All", "Total Sum Event with ONLY 3 Clusters Found, Maximum Offline Combination Energy;E_{Sum};Counts",en+1000,0, en+1000);
    //VTP 3 Cluster events for events where only 3 clusters were found.
    h_tSum_3ClustOnly_VTP_800 = new TH1F("h_tSum_2ClustOnly_VTP_800", "Total Sum Event with ONLY 3 Clusters Found on the VTP with 800 MeV sum threshold;E_{Sum};Counts",en+1000,0, en+1000);
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

    h_tSum_2ClustOnly_All->SetLineColor(kBlack);
    h_tSum_2ClustOnly_VTP_800->SetLineColor(kBlue);

    TH1F* h_tSum_2ClustOnlyRatio = (TH1F*) h_tSum_2ClustOnly_VTP_800->Clone();
    h_tSum_2ClustOnlyRatio->Divide(h_tSum_2ClustOnly_All);
    h_tSum_2ClustOnlyRatio->SetLineColor(kBlack);
    h_tSum_2ClustOnlyRatio->SetTitle("Exactly 2 Clusters Trigger Efficiency;E_{Sum};Efficiency");

    h_tSum_3Clust_All->SetLineColor(kBlack);
    h_tSum_3Clust_VTP_800->SetLineColor(kBlue);

    TH1F* h_tSum_3ClustRatio = (TH1F*) h_tSum_3Clust_VTP_800->Clone();
    h_tSum_3ClustRatio->Divide(h_tSum_3Clust_All);
    h_tSum_3ClustRatio->SetLineColor(kBlack);
    h_tSum_3ClustRatio->SetTitle("Trigger Efficiency;E_{Sum};Efficiency");

    h_tSum_3ClustOnly_All->SetLineColor(kBlack);
    h_tSum_3ClustOnly_VTP_800->SetLineColor(kBlue);

    TH1F* h_tSum_3ClustOnlyRatio = (TH1F*) h_tSum_3ClustOnly_VTP_800->Clone();
    h_tSum_3ClustOnlyRatio->Divide(h_tSum_3ClustOnly_All);
    h_tSum_3ClustOnlyRatio->SetLineColor(kBlack);
    h_tSum_3ClustOnlyRatio->SetTitle("Exactly 3 Clusters Trigger Efficiency;E_{Sum};Efficiency");

    TCanvas *c = new TCanvas("c", "TSum_ClustTrig_Canvas",1000,1000);
    auto legend = new TLegend(0.1,0.8,0.4,0.9);
    auto legend2 = new TLegend(0.1,0.8,0.4,0.9);

    //Page 1
    //-----------------------------------------------------------------------------------
    c->Divide(1,2);
    c->cd(1);
    h_tSum_1Clust_All->Draw("HIST");
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

    //Page 5
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

    //Page 6
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

}