/**
 * Header file for LMS Trigger object to evaluate for LMS triggers.
 * @author Erik Wrightson
 * @version 04.08.2026
 * @creation 04.08.2026
 */

#ifndef LMS_Trig_H
#define LMSTrig_H

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

class LMSTrig: public Trigger{

    public:

        //Statically set the flags for this trigger type.
        static const Int_t EXT_FLAG = 0x0FF0000;
        static const Int_t RAND_FLAG = 0x01000000;

        //Constuctor that ensures the chain is set up via the parent constructor.
        LMSTrig(TChain* c);

        //Processing Function inhereted from parent.
        void ProcessData(bool self, bool rand, bool tSum) override;

    private:
        Double_t rand_LMS_Elligible;
        Double_t rand_LMS_Found;

        Double_t tSum_LMS_Elligible;
        Double_t tSum_LMS_Found;


};

#endif