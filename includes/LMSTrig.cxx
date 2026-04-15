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
LMSTrig::LMSTrig(TChain* c):Trigger(c, false){
    cout<<"Set All Branch Addresses for LMS Trigger Object.\n";
    rand_LMS_Elligible = 0;
    rand_LMS_Found = 0;

    tSum_LMS_Elligible = 0;
    tSum_LMS_Found = 0;
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
void LMSTrig::ProcessData(bool self, bool rand, bool tSum){
    //Get the amount of entries from each file to limit looping through them.
	Long64_t entries = chain->GetEntries();

    rand_LMS_Elligible = 0;
    rand_LMS_Found = 0;

    tSum_LMS_Elligible = 0;
    tSum_LMS_Found = 0;

    tSum_LMS_BothFire = 0;
    
    for(Int_t i = 0; i < entries; i++){
        chain->GetEntry(i);

        if(self && (trigType & (1<<LMS_TFLAG)) && (trigger_bits & LMS_FLAG)){

        }
        if(rand && (trigger_bits & (1<<RAND_FC_FLAG))){

        }
        if(tSum && (trigger_bits & (1<<TSUM_FLAG))){ //(trigType & (1<<SSP_RAWSUM_TFLAG))
            Int_t crysFire = 0;
            
            if(trigger_bits & (1<<LMS_FLAG)){tSum_LMS_BothFire++;}

            if(lms_numChan==3){
                
                for(int x = 0; x < numChan; x++){
                    //Checks that all 3 LMS Reference PMT Channels fired and that at least 3 crystals fired.
                    //Technically all should fire, but limiting to three allows for the least amount of channels to be checked.
                    if(crysFire >= 3){
                        tSum_LMS_Elligible++;
                        if((trigger_bits & (1<<LMS_FLAG))){tSum_LMS_Found++;} //(trigType & EXT_FLAG) &&
                        break;
                    }

                    if(modId[x] > 1000 && modId[x]<3000){crysFire++;}
                }
            }

        }
    }
}

//Getter methods
Double_t LMSTrig::get_Rand_LMS_Elligible(){return rand_LMS_Elligible;}
Double_t LMSTrig::get_Rand_LMS_Found(){return rand_LMS_Found;}
Double_t LMSTrig::get_tSum_LMS_Elligible(){return tSum_LMS_Elligible;}
Double_t LMSTrig::get_tSum_LMS_Found(){return tSum_LMS_Found;}
Double_t LMSTrig::get_tSum_LMS_BothFired(){return tSum_LMS_BothFire;}