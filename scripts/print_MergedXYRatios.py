import ROOT
import sys

from ROOT import kBlack
from ROOT import kBlue
from ROOT import kMagenta
from ROOT import kPink
from ROOT import kRed
from ROOT import kGreen
from ROOT import kCyan
from ROOT import kSpring
from ROOT import gPad

def extract_histograms(directory, path=""):
    """Recursively traverses a ROOT directory to find all histograms."""
    histograms = {}
    
    # Loop over all keys inside the current directory
    for key in directory.GetListOfKeys():
        name = key.GetName()
        classname = key.GetClassName()
        full_path = f"{path}/{name}" if path else name
        
        # Retrieve the object from the file
        obj = key.ReadObj()
        
        # Check if the object is a subdirectory
        if isinstance(obj, ROOT.TDirectoryFile):
            # Recurse into the subdirectory
            histograms.update(extract_histograms(obj, full_path))
            
        # Check if the object is a histogram (inherits from TH1)
        elif obj.InheritsFrom("TH1"):
            # Detach histogram from the directory to keep it in memory after file closes
            obj.SetDirectory(0)  
            histograms[full_path] = obj
            print(f"Found histogram [{classname}]: {full_path}")
            
    return histograms

# Main execution
if __name__ == "__main__":

    #if(sys.argc == 2):
        #print("Please input only the root file to print.")
    file_path = sys.argv[1]

    pdfName = file_path
    pdfName.replace(".root", "")
    pdfName = pdfName[:pdfName.rfind("/")]+"_XY_Ratios.pdf"
    
    # Open the ROOT file in read mode
    root_file = ROOT.TFile.Open(file_path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"Error: Could not open file {file_path}")
        exit(1)
        
    print(f"Scanning {file_path} for histograms...")
    all_histograms = extract_histograms(root_file)

    # Close the file safely
    root_file.Close()
    
    print(f"\nSuccessfully read {len(all_histograms)} histograms.")

    xyHistos = [hist for name, hist in all_histograms.items() if "XY" in name]#[y for y in all_histograms if "XY" in y.GetName()]
    xyAllHist = [u for u in xyHistos if "All" in u.GetName()]
    xyVTPHist = [i for i in xyHistos if "VTP" in i.GetName()]
    xyRatioHist = []
    h_XY_Ratio_ProjX = []
    h_XY_Ratio_ProjY = []

    h_XY_All_noMid = xyAllHist[0].Clone("h_XY_All_noMid")
    for j in range(16,20):
            for k in range(16,20):
                h_XY_All_noMid.SetBinContent(j, k, 0.0)
                h_XY_All_noMid.SetBinError(j, k, 0.0)
    h_XY_All_ProjX = h_XY_All_noMid.ProjectionX("h_XY_All_ProjX")
    h_XY_All_ProjY = h_XY_All_noMid.ProjectionY("h_XY_All_ProjY")

    bitThr = ["650", "750", "950", "1300", "1100", "1400", "1500"]

    xyTemp2 = None
    for i in range(0,len(xyVTPHist)):
        xyTemp = xyVTPHist[i].Clone()
        
        xyTemp2 = xyVTPHist[i].Clone()
        for j in range(16,20):
            for k in range(16,20):
                xyTemp2.SetBinContent(j, k, 0.0)
                xyTemp2.SetBinError(j, k, 0.0)
        
        xyTemp_ProjX = xyTemp2.ProjectionX("h_XY_Ratio_ProjX_bit" + str(i))
        xyTemp_ProjX.Divide(xyTemp_ProjX, h_XY_All_ProjX, 1.0, 1.0, "B")
        xyTemp_ProjX.SetTitle("X Projection Weighted Efficiency Bit - " + str(i) + " " + bitThr[i] + "MeV")
        h_XY_Ratio_ProjX.append(xyTemp_ProjX)

        xyTemp_ProjY = xyTemp2.ProjectionY("h_XY_Ratio_ProjY_bit" + str(i))
        xyTemp_ProjY.Divide(xyTemp_ProjY, h_XY_All_ProjY, 1.0, 1.0, "B")
        xyTemp_ProjY.SetTitle("y Projection Weighted Efficiency Bit - " + str(i) + " " + bitThr[i] + "MeV")
        h_XY_Ratio_ProjY.append(xyTemp_ProjY)

        xyTemp.Divide(xyVTPHist[i], xyAllHist[0], 1.0, 1.0, "B")
        xyTemp.SetTitle("XY Trigger Efficiency Bit " + str(i) + " " + bitThr[i] + "MeV")
        xyRatioHist.append(xyTemp)
    
    c = ROOT.TCanvas("c1", "Merged Canvas", 1000, 1000)
    legend = ROOT.TLegend(0.1,0.8,0.4,0.9)

    
    '''alt = 0
    for i in range(0,len(xyVTPHist)):
        h_XY_Ratio_ProjX.append(xyRatioHist[i].ProjectionX("h_XY_Ratio" + str(i) + "_ProjX"))
        h_XY_Ratio_ProjX[i].Scale(1.0/34)
        h_XY_Ratio_ProjX[i].SetTitle("X Projection Efficiency Bit - " + str(i) + " " + bitThr[i] + "MeV")
        
        midbinX1 = h_XY_Ratio_ProjX[i].GetBinContent(17) * 34.0/32.0
        h_XY_Ratio_ProjX[i].SetBinContent(17, midbinX1)
        midbinX1_err = h_XY_Ratio_ProjX[i].GetBinError(17) * 34.0/32.0
        h_XY_Ratio_ProjX[i].SetBinError(17, midbinX1_err)
        
        midbinX2 = h_XY_Ratio_ProjX[i].GetBinContent(18) * 34.0/32.0
        h_XY_Ratio_ProjX[i].SetBinContent(18, midbinX2)
        midbinX2_err = h_XY_Ratio_ProjX[i].GetBinError(18) * 34.0/32.0
        h_XY_Ratio_ProjX[i].SetBinError(18, midbinX2_err)

        h_XY_Ratio_ProjY.append(xyRatioHist[i].ProjectionY("h_XY_Ratio" + str(i) + "_ProjY"))
        h_XY_Ratio_ProjY[i].Scale(1.0/34)
        h_XY_Ratio_ProjY[i].SetTitle("Y Projection Efficiency Bit - " + str(i) + " " + bitThr[i] + "MeV")

        midbinY1 = h_XY_Ratio_ProjY[i].GetBinContent(17) * 34.0/32.0
        h_XY_Ratio_ProjY[i].SetBinContent(17, midbinY1)
        midbinY1_err = h_XY_Ratio_ProjY[i].GetBinError(17) * 34.0/32.0
        h_XY_Ratio_ProjY[i].SetBinError(17, midbinY1_err)
        
        midbinY2 = h_XY_Ratio_ProjY[i].GetBinContent(18) * 34.0/32.0
        h_XY_Ratio_ProjY[i].SetBinContent(18, midbinY2)
        midbinY2_err = h_XY_Ratio_ProjY[i].GetBinError(18) * 34.0/32.0
        h_XY_Ratio_ProjY[i].SetBinError(18, midbinY2_err)'''
    
    Bit_Color = [kBlue, kRed, kGreen,kCyan, kMagenta, kSpring, 28]

    #-----------------------------------------------------------------
    c.Divide(2,2)
    c.cd(1)
    gPad.SetLogz(1)
    xyAllHist[0].Draw("COLZ")
    c.cd(2)
    gPad.SetLogz(1)
    xyVTPHist[0].Draw("COLZ")
    c.cd(3)
    gPad.SetLogz(1)
    xyVTPHist[1].Draw("COLZ")
    c.cd(4)
    gPad.SetLogz(1)
    xyVTPHist[2].Draw("COLZ")
    c.Print(pdfName + "(")
    c.Clear()

    #-----------------------------------------------------------------
    c.Divide(2,2)
    c.cd(1)
    gPad.SetLogz(1)
    xyVTPHist[3].Draw("COLZ")
    c.cd(2)
    gPad.SetLogz(1)
    xyVTPHist[4].Draw("COLZ")
    c.cd(3)
    gPad.SetLogz(1)
    xyVTPHist[5].Draw("COLZ")
    c.cd(4)
    gPad.SetLogz(1)
    xyVTPHist[6].Draw("COLZ")
    c.Print(pdfName)
    c.Clear()

    #-----------------------------------------------------------------
    for i in range(0, len(xyRatioHist)):
        c.Divide(2,2)
        c.cd(1)
        xyRatioHist[i].SetStats(0)
        xyRatioHist[i].Draw("COLZ")
        c.cd(3)
        h_XY_Ratio_ProjX[i].SetMarkerStyle(20)
        h_XY_Ratio_ProjX[i].Draw("P E")
        c.cd(4)
        h_XY_Ratio_ProjY[i].SetMarkerStyle(20)
        h_XY_Ratio_ProjY[i].Draw("P E")
        c.Print(pdfName)
        c.Clear()

    c.Divide(2,2)
    c.cd(1)
    gPad.SetLogz(1)
    xyTemp2.Draw("COLZ")
    c.cd(2)
    gPad.SetLogz(1)
    h_XY_All_noMid.Draw("COLZ")
    c.Print(pdfName + ")")