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
 * @param gem - boolean flag for if any GEM information is needed for this trigger.
 */
Trigger::Trigger(TChain* c, bool gem){
    chain = (TChain*) c;

    chain->SetMakeClass(1);

    //General Event Data
    chain->SetBranchAddress("event_num",       &eventNum);
    chain->SetBranchAddress("trigger_type",    &trigType);
    chain->SetBranchAddress("trigger_bits",    &trigger_bits);
    chain->SetBranchAddress("timestamp",       &time);

    //HyCal Information
    chain->SetBranchAddress("hycal.nch",       &numChan);
    /*chain->SetBranchAddress("hycal.crate",     crate);
    chain->SetBranchAddress("hycal.slot",      slot);
    chain->SetBranchAddress("hycal.channel",   channel);*/
    chain->SetBranchAddress("hycal.module_id", modId);
    chain->SetBranchAddress("hycal.nsamples",  nSamps);
    chain->SetBranchAddress("hycal.samples",   samps);
    chain->SetBranchAddress("hycal.ped_mean",  pedMean);
    chain->SetBranchAddress("hycal.ped_rms",   pedRMS);
    chain->SetBranchAddress("hycal.integral",  integral);
    
    //HyCal Waveform Information
    if (chain->GetBranch("hycal.npeaks")) {
        chain->SetBranchAddress("hycal.npeaks",        nPeaks);
        chain->SetBranchAddress("hycal.peak_height",   peakHeight);
        chain->SetBranchAddress("hycal.peak_time",     peakTime);
        chain->SetBranchAddress("hycal.peak_integral", peakIntegral);
    }
    
    //GEM Information
    if(gem && chain->GetBranch("gem.nch")){
        chain->Branch("gem.nch",          gem_numChan);
        chain->Branch("gem.mpd_crate",    mpdCrate);
        chain->Branch("gem.mpd_fiber",    mpdFiber);
        chain->Branch("gem.apv",          apv);
        chain->Branch("gem.strip",        strip);
        chain->Branch("gem.ssp_samples",  sspSamp);
        chain->Branch("n_ssp_triggers",   numSSPTrigs);
        chain->Branch("ssp_trigger_tags", sspTrigTags);
    }


    //LMS Information
    chain->SetBranchAddress("lms.nch",       &lms_numChan);
    /*chain->SetBranchAddress("lms.crate",     lms_crate);
    chain->SetBranchAddress("lms.slot",      lms_slot);
    chain->SetBranchAddress("lms.channel",   lms_channel);*/
    chain->SetBranchAddress("lms.id", lms_modId);
    chain->SetBranchAddress("lms.nsamples",  lms_nSamps);
    chain->SetBranchAddress("lms.samples",   lms_samps);
    chain->SetBranchAddress("lms.ped_mean",  lms_pedMean);
    chain->SetBranchAddress("lms.ped_rms",   lms_pedRMS);
    chain->SetBranchAddress("lms.integral",  lms_integral);
    
    //LMS Waveform Information
    if (chain->GetBranch("lms.npeaks")) {
        chain->SetBranchAddress("lms.npeaks",        lms_nPeaks);
        chain->SetBranchAddress("lms.peak_height",   lms_peakHeight);
        chain->SetBranchAddress("lms.peak_time",     lms_peakTime);
        chain->SetBranchAddress("lms.peak_integral", lms_peakIntegral);
    }

    //Veto Information
    chain->SetBranchAddress("veto.nch",       &veto_numChan);
    /*chain->SetBranchAddress("veto.crate",     veto_crate);
    chain->SetBranchAddress("veto.slot",      veto_slot);
    chain->SetBranchAddress("veto.channel",   veto_channel);*/
    chain->SetBranchAddress("veto.id", veto_modId);
    chain->SetBranchAddress("veto.nsamples",  veto_nSamps);
    chain->SetBranchAddress("veto.samples",   veto_samps);
    chain->SetBranchAddress("veto.ped_mean",  veto_pedMean);
    chain->SetBranchAddress("veto.ped_rms",   veto_pedRMS);
    chain->SetBranchAddress("veto.integral",  veto_integral);
    
    //Veto Waveform Information
    if (chain->GetBranch("veto.npeaks")) {
        chain->SetBranchAddress("veto.npeaks",        veto_nPeaks);
        chain->SetBranchAddress("veto.peak_height",   veto_peakHeight);
        chain->SetBranchAddress("veto.peak_time",     veto_peakTime);
        chain->SetBranchAddress("veto.peak_integral", veto_peakIntegral);
    }
    
}

/**
 * Processes the events one by one depending on which of the flags are set to be true.
 *
 * @param self - Whether or not this trigger should have evaluation done within its own events. {True for Yes; False for No}
 * @param rand - Whether or not this trigger should have evalutaion done using a random trigger source. {True for Yes; False for No}
 * @param tSum - Whether or not this trigger should have evaluation done using total sum trigger source. {True for Yes; False for No}
 */
void Trigger::ProcessData(bool self, bool rand, bool tSum){
    cout<<"This should be overriden by all child classes.";
}