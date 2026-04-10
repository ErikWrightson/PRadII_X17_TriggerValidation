/**
 * Header file for looping through the root data and generating histograms needed.
 * @author Erik Wrightson
 * @version 04.07.2026
 * @creation 04.07.2026
 */

#ifndef FileProcessUtil_H
#define FileProcessUtil_H

#include <iostream>
#include <string.h>
#include <fstream>
#include <TString.h>
#include <TChain.h>

using namespace std;

namespace FileProcessUtil{

    /**
     * Processes a file with the given name and proper formatting and returns a vector with the 
     * names of each file that is listed within.
     *
     * @param fileListFileName - the name of the file that contains the list of files to process.
     */
    static vector<TString> processFileList(string fileListFileName);

    /**
     * Makes TChain of the given vector assuming the file names are correct.
     *
     * @param names - the vector of file names to add to the TChain.
     */
     TChain* makeChain(vector<TString> names);
}
#endif