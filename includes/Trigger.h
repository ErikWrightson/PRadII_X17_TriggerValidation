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

        //Trigger Type Flags
        static constexpr UInt_t SSP_RAWSUM_TFLAG = 8; // SSP Raw Sum trigger type flag
        static constexpr UInt_t CLUSTER_TFLAG = 9;    //Clustering trigger type flag
        static constexpr UInt_t PULSER_TFLAG = 15;    //Pulser trigger type flag
        static constexpr UInt_t LMS_TFLAG = 24;       //LMS trigger type flag
        static constexpr UInt_t ALPHA_TFLAG = 25;     //Alpha trigger type flag
        static constexpr UInt_t MOR_TFLAG = 27;       //Master OR trigger type flag
        
        //Trigger Bits Flags
        static constexpr UInt_t TSUM_FLAG = 8; //Total Sum Bit flag
        //VTP Clustering Trigger Flags
        static constexpr UInt_t CL1_FLAG = 9;  // Single Cluster trigger flag
        static constexpr UInt_t CL2_FLAG = 10;  // Two Cluster trigger flag
        static constexpr UInt_t CL3_FLAG = 11;  // Three Cluster trigger flag
        static constexpr UInt_t PULSER_FLAG = 15; // Pulser flag
        static constexpr UInt_t LMS_FLAG = 24;    // LMS LED pulse trigger flag
        static constexpr UInt_t ALPHA_FLAG = 25;  // LMS Alpha pulse flag
        static constexpr UInt_t RAND_FC_FLAG = 26;    // Random Trigger from the Faraday Cup flag
        static constexpr UInt_t MOR_FLAG = 27;
        

        //Limits for arrays.
        static constexpr Int_t MAX_SAMPLES  = 200;       // samples per channel per event
        static constexpr Int_t MAX_CHANNELS = 64;        // channels per slot (16 for FADC250, 64 for ADC1881M)
        static constexpr Int_t MAX_ROCS     = 10;        // number of ROC crates
        static constexpr Int_t MAX_PEAKS    = 8;         // max peaks per channel waveform
        static constexpr Int_t APV_STRIP_SIZE     = 128; // channels per APV25 chip
        static constexpr Int_t SSP_TIME_SAMPLES   = 6;   // fixed by SSP firmware
        static constexpr Int_t MAX_APVS_PER_MPD   = 16;  // APV slots per MPD

        //
        static constexpr Int_t MAX_HC_CHANS = 1728;      //Maximum number of HyCal Channels that could fire.
        static constexpr Int_t MAX_GEM_CHANS = 1000;     //Maximum GEM channels that could fire.
        static constexpr Int_t MAX_NUM_SSP_TRIGS = 1000; //Maximum number of SSP triggers that could fire.

        TChain* chain;

        //Locations to store the event data
        UInt_t eventNum;
        UChar_t trigType;
        UInt_t trigger;
        Long64_t time;

        //Locations to store the HyCal data.
        Int_t numChan;
        UChar_t crate[MAX_HC_CHANS];
        UChar_t slot[MAX_HC_CHANS];
        UChar_t channel[MAX_HC_CHANS];
        UShort_t modId[MAX_HC_CHANS];
        UChar_t nSamps[MAX_HC_CHANS];
        UShort_t samps[MAX_HC_CHANS][MAX_SAMPLES];
        Float_t pedMean[MAX_HC_CHANS];
        Float_t pedRMS[MAX_HC_CHANS];
        Float_t integral[MAX_HC_CHANS];
        
        //Only sometimes used for storing the HyCal waveform information.
        UChar_t nPeaks[MAX_HC_CHANS];
        Float_t peakHeight[MAX_HC_CHANS][MAX_PEAKS];
        Float_t peakTime[MAX_HC_CHANS][MAX_PEAKS];
        Float_t peakIntegral[MAX_HC_CHANS][MAX_PEAKS];

        //GEM Variables
        UChar_t gem_numChan;
        UChar_t mpdCrate[MAX_GEM_CHANS];
        UChar_t mpdFiber[MAX_GEM_CHANS];
        UChar_t apv[MAX_GEM_CHANS];
        UChar_t strip[MAX_GEM_CHANS];
        Short_t sspSamp[MAX_GEM_CHANS][SSP_TIME_SAMPLES];
        Int_t numSSPTrigs;
        UInt_t sspTrigTags[MAX_NUM_SSP_TRIGS][SSP_TIME_SAMPLES];

        //Parent Constructor
        Trigger(TChain* c, bool gem);

        //Virutal processor function to be overriden by children.
        virtual void ProcessData(bool self, bool rand, bool tSum);


};

#endif