/**
 * Form an object with the proper functionality to find and make all the plots for LMS Trigger Evaluation.
 *
 * @author Erik Wrightson
 * @version 04.08.2026
 * @creation 04.08.2026
 */
#include "Trigger.h"
#include "LMSTrig.h"

/**
 * Creates the LMS trigger object.
 *
 * @param c - the TChain of events to
 */
LMSTrig::LMSTrig(TChain* c):Trigger(c){
    cout<<"Set All Branch Addresses for LMS Trigger Object.\n";
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
void LMSTrig::ProcessData(bool self, bool rand, bool tSum) override{
    //Get the amount of entries from each file to limit looping through them.
	Long64_t entries = chain->GetEntries();
    
    for(Int_t i = 0; i < entries; i++){
        chain->GetEntry(i);

        if(self){

        }
        if(rand){

        }
        if(tSum){

        }
    }
}