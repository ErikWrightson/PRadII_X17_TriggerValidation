/**
 * Trigger parent class that instantiates and sets up the trees that all the trigger types have in common.
 *
 * @author Erik Wrightson
 * @version 04.10.2026
 * @creation 04.10.2026
 */
#include "Trigger.h"

/**
 * Creates a trigger object and then sets the known branches that are in the PRadII/X17 level 1 replay trees.
 *
 * @param c - the TChain to set up the branching and trees for.
 */
Trigger::Trigger(TChain* c){
    chain = (TChain*) c;

    chain->SetMakeClass(1);

    //General Event Data
    chain->SetBranchAddress("event_num",       &eventNum);
    chain->SetBranchAddress("trigger_type",    &trigType);
    chain->SetBranchAddress("trigger",         &trigger);
    chain->SetBranchAddress("timestamp",       &time);

    //HyCal Information
    chain->SetBranchAddress("hycal.nch",       &numChan);
    chain->SetBranchAddress("hycal.crate",     crate);
    chain->SetBranchAddress("hycal.slot",      slot);
    chain->SetBranchAddress("hycal.channel",   channel);
    chain->SetBranchAddress("hycal.module_id", modId);
    chain->SetBranchAddress("hycal.nsamples",  nSamps);
    chain->SetBranchAddress("hycal.samples",   samps);
    chain->SetBranchAddress("hycal.ped_mean",  pedMean);
    chain->SetBranchAddress("hycal.ped_rms",   pedRMS);
    chain->SetBranchAddress("hycal.integral",  integral);
    
    //Waveform Information
    if (chain->GetBranch("hycal.npeaks")) {
        chain->SetBranchAddress("hycal.npeaks",        nPeaks);
        chain->SetBranchAddress("hycal.peak_height",   peakHeight);
        chain->SetBranchAddress("hycal.peak_time",     peakTime);
        chain->SetBranchAddress("hycal.peak_integral", peakIntegral);
    }
    
    //GEM Information
    if(chain->GetBranch("gem.nch")){
        chain->Branch("gem.nch",          gem_numChan);
        chain->Branch("gem.mpd_crate",    mpdCrate);
        chain->Branch("gem.mpd_fiber",    mpdFiber);
        chain->Branch("gem.apv",          apv);
        chain->Branch("gem.strip",        strip);
        chain->Branch("gem.ssp_samples",  sspSamp);
        chain->Branch("n_ssp_triggers",   numSSPTrigs);
        chain->Branch("ssp_trigger_tags", sspTrigTags);
    }
    
}

/**
 * Processes the events one by one depending on which of the flags are set to be true.
 *
 * @param self - Whether or not this trigger should have evaluation done within its own events. {True for Yes; False for No}
 * @param rand - Whether or not this trigger should have evalutaion done using a random trigger source. {True for Yes; False for No}
 * @param tSum - Whether or not this trigger should have evaluation done using total sum trigger source. {True for Yes; False for No}
 */
virtual void Trigger::ProcessData(bool self, bool rand, bool tSum){
    cout<<"This should be overriden by all child classes.";
}