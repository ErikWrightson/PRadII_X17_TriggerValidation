/**
 * Header file for Trigger parent class.
 * This allows one definition to be made for all of the attributes that each trigger type share, like loading the TChain and Branches.
 *
 * @author Erik Wrightson
 * @version 04.10.2026
 * @creation 04.10.2026
 */

#ifndef Trigger_H
#define Trigger_H

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

class Trigger{

    public:

        //Statically set the flags for each trigger type.
        static constexpr Int_t EXT_FLAG = 0x0FF0000;
        static constexpr Int_t RAND_FLAG = 0x01000000;
        static constexpr Int_t MAX_SAMPLES  = 200;   // samples per channel per event
        static constexpr Int_t MAX_CHANNELS = 64;   // channels per slot (16 for FADC250, 64 for ADC1881M)
        static constexpr Int_t MAX_ROCS     = 10;    // number of ROC crates
        static constexpr Int_t MAX_PEAKS    = 8;    // max peaks per channel waveform
        static constexpr Int_t APV_STRIP_SIZE     = 128;   // channels per APV25 chip
        static constexpr Int_t SSP_TIME_SAMPLES   = 6;     // fixed by SSP firmware
        static constexpr Int_t MAX_APVS_PER_MPD   = 16;    // APV slots per MPD

        //
        static constexpr Int_t MAX_HC_CHANS = 1728;
        static constexpr Int_t MAX_GEM_CHANS = 1000;

        TChain* chain;

        //Locations to store the event data
        UInt_t eventNum;
        UInt8_t trigType;
        UInt_t trigger;
        Long64_t time;

        //Locations to store the HyCal data.
        Int_t numChan;
        UInt8_t crate[MAX_HC_CHANS];
        UInt8_t slot[MAX_HC_CHANS];
        UInt8_t channel[MAX_HC_CHANS];
        UShort_t modId[MAX_HC_CHANS];
        UInt8_t nSamps[MAX_HC_CHANS];
        UShort_t samps[MAX_HC_CHANS][MAX_SAMPLES];
        Float_t pedMean[MAX_HC_CHANS];
        Float_t pedRMS[MAX_HC_CHANS];
        Float_t integral[MAX_HC_CHANS];
        
        //Only sometimes used for storing the HyCal waveform information.
        UInt8_t nPeaks[MAX_HC_CHANS];
        Float_t peakHeight[MAX_HC_CHANS][MAX_PEAKS];
        Float_t peakTime[MAX_HC_CHANS][MAX_PEAKS];
        Float_t peakIntegral[MAX_HC_CHANS][MAX_PEAKS];

        //GEM Variables
        UInt_t gem_numChan;
        UInt8_t mpdCrate[gem_numChan];
        UInt8_t mpdFiber[gem_numChan];
        UInt8_t apv[gem_numChan];
        UInt8_t strip[gem_numChan];
        Short_t sspSamp[gem_numChan][SSP_TIME_SAMPLES];
        Int_t numSSPTrigs;
        UInt_t sspTrigTags[numSSPTrigs][SSP_TIME_SAMPLES];

        //Parent Constructor
        Trigger(TChain* c);

        //Virutal processor function to be overriden by children.
        virtual void ProcessData(bool self, bool rand, bool tSum);


}

#endif