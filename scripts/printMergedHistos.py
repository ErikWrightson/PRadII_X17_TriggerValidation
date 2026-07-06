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
    pdfName = pdfName[:pdfName.rfind("/")]+"AfterMerge.pdf"
    
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

    sumHistos = [hist for name, hist in all_histograms.items() if "Sum" in name]
    sumAllHist = [w for w in sumHistos if "All" in w.GetName()]
    sumVTPHist = [q for q in sumHistos if "VTP" in q.GetName()]
    sumRatioHist = []
    
    minHistos = [hist for name, hist in all_histograms.items() if "Min" in name]#[e for e in all_histograms if "Min" in e.GetName()]
    minAllHist = [r for r in minHistos if "All" in r.GetName()]
    minVTPHist = [t for t in minHistos if "VTP" in t.GetName()]
    minRatioHist = []

    maxHistos = [hist for name, hist in all_histograms.items() if "Max" in name]#[e for e in all_histograms if "Max" in e.GetName()]
    maxAllHist = [r for r in maxHistos if "All" in r.GetName()]
    maxVTPHist = [t for t in maxHistos if "VTP" in t.GetName()]
    maxRatioHist = []

    xyHistos = [hist for name, hist in all_histograms.items() if "XY" in name]#[y for y in all_histograms if "XY" in y.GetName()]
    xyAllHist = [u for u in xyHistos if "All" in u.GetName()]
    xyVTPHist = [i for i in xyHistos if "VTP" in i.GetName()]
    xyRatioHist = []

    dtHistos = [hist for name, hist in all_histograms.items() if "deltaT" in name]#[o for o in all_histograms if "deltaT" in o.GetName()]
    dtAllHist = [p for p in dtHistos if "All" in p.GetName()]
    dtVTPHist = [a for a in dtHistos if "VTP" in a.GetName()]
    dtRatioHist = []

    timeCorrHist = [hist for name, hist in all_histograms.items() if "Corr" in name]#[o for o in all_histograms if "Corr" in o.GetName()]
    
    for i in range(0,len(dtVTPHist)):
        sumTemp = sumVTPHist[i].Clone()
        sumTemp.Divide(sumVTPHist[i], sumAllHist[0], 1.0, 1.0, "B")
        sumTemp.SetTitle("Trigger Efficiency")
        sumRatioHist.append(sumTemp)

        minTemp = minVTPHist[i].Clone()
        minTemp.Divide(minVTPHist[i], minAllHist[0], 1.0, 1.0, "B")
        minTemp.SetTitle("Trigger Efficiency")
        minRatioHist.append(minTemp)

        maxTemp = maxVTPHist[i].Clone()
        maxTemp.Divide(maxVTPHist[i], maxAllHist[0], 1.0, 1.0, "B")
        maxTemp.SetTitle("Trigger Efficiency")
        maxRatioHist.append(maxTemp)

        dtTemp = dtVTPHist[i].Clone()
        dtTemp.Divide(dtVTPHist[i], dtAllHist[0], 1.0, 1.0, "B")
        dtTemp.SetTitle("Trigger Efficiency")
        dtRatioHist.append(dtTemp)
    
    c = ROOT.TCanvas("c1", "Merged Canvas", 1000, 1000)
    legend = ROOT.TLegend(0.1,0.8,0.4,0.9)
    bitThr = ["650", "750", "950", "1300", "1100", "1400", "1500"]
    Bit_Color = [kBlue, kRed, kGreen,kCyan, kMagenta, kSpring, 28]

    #-----------------------------------------------------------
    c.Divide(1,2)
    c.cd(1)
    sumAllHist[0].SetLineColor(kBlack)
    sumAllHist[0].Draw("HIST")
    legend.AddEntry(sumAllHist[0], "All", "l")
    for ind, s in enumerate(sumVTPHist):
        s.SetLineColor(Bit_Color[ind])
        s.Draw("HIST SAME")
        legend.AddEntry(s, bitThr[ind], "l")
    legend.Draw()  

    c.cd(2)
    sumRatioHist[0].SetLineColor(Bit_Color[0])
    sumRatioHist[0].Draw("HIST")
    sumRatioHist[0].Draw("E1 SAME")
    for d in range(1,len(sumRatioHist)):
        sumRatioHist[d].SetLineColor(Bit_Color[d])
        sumRatioHist[d].Draw("HIST SAME")

    c.Print(pdfName + "(")
    c.Clear()
    legend.Clear()

    #-----------------------------------------------------------
    c.Divide(1,2)
    c.cd(1)
    sumAllHist[0].SetLineColor(kBlack)
    sumAllHist[0].Draw("HIST")
    legend.AddEntry(sumAllHist[0], "All", "l")
    for ind, s in enumerate(sumVTPHist):
        s.Draw("HIST SAME")
        legend.AddEntry(s, bitThr[ind], "l")
    legend.Draw()  

    c.cd(2)
    sumRatioHist[0].SetAxisRange(0.9,1.03,"Y")
    sumRatioHist[0].Draw("HIST")
    sumRatioHist[0].Draw("E1 SAME")
    for d in range(1,len(sumRatioHist)):
        sumRatioHist[d].Draw("HIST SAME")

    c.Print(pdfName)
    c.Clear()
    legend.Clear()

    #-----------------------------------------------------------
    c.Divide(1,2)
    c.cd(1)
    minAllHist[0].SetAxisRange(0,1300,"X")
    minAllHist[0].SetLineColor(kBlack)
    minAllHist[0].Draw("HIST")
    legend.AddEntry(minAllHist[0], "All", "l")
    for ind, s in enumerate(minVTPHist):
        s.SetLineColor(Bit_Color[ind])
        s.Draw("HIST SAME")
        legend.AddEntry(s, bitThr[ind], "l")
    legend.Draw()  

    c.cd(2)
    minRatioHist[0].SetAxisRange(0,1300,"X")
    minRatioHist[0].SetLineColor(Bit_Color[0])
    minRatioHist[0].Draw("HIST")
    minRatioHist[0].Draw("E1 SAME")
    for d in range(1,len(minRatioHist)):
        minRatioHist[d].SetLineColor(Bit_Color[d])
        minRatioHist[d].Draw("HIST SAME")

    c.Print(pdfName + "(")
    c.Clear()
    legend.Clear()

    #-----------------------------------------------------------
    c.Divide(1,2)
    c.cd(1)
    minAllHist[0].Draw("HIST")
    legend.AddEntry(minAllHist[0], "All", "l")
    for ind, s in enumerate(minVTPHist):
        s.Draw("HIST SAME")
        legend.AddEntry(s, bitThr[ind], "l")
    legend.Draw()  

    c.cd(2)
    minRatioHist[0].SetAxisRange(0.9,1.03,"Y")
    minRatioHist[0].Draw("HIST")
    minRatioHist[0].Draw("E1 SAME")
    for d in range(1,len(minRatioHist)):
        minRatioHist[d].Draw("HIST SAME")

    c.Print(pdfName)
    c.Clear()
    legend.Clear()

     #-----------------------------------------------------------
    c.Divide(1,2)
    c.cd(1)
    maxAllHist[0].SetLineColor(kBlack)
    maxAllHist[0].Draw("HIST")
    legend.AddEntry(maxAllHist[0], "All", "l")
    for ind, s in enumerate(maxVTPHist):
        s.SetLineColor(Bit_Color[ind])
        s.Draw("HIST SAME")
        legend.AddEntry(s, bitThr[ind], "l")
    legend.Draw()  

    c.cd(2)
    maxRatioHist[0].SetLineColor(Bit_Color[0])
    maxRatioHist[0].Draw("HIST")
    maxRatioHist[0].Draw("E1 SAME")
    for d in range(1,len(maxRatioHist)):
        maxRatioHist[d].SetLineColor(Bit_Color[d])
        maxRatioHist[d].Draw("HIST SAME")

    c.Print(pdfName + "(")
    c.Clear()
    legend.Clear()

    #-----------------------------------------------------------
    c.Divide(1,2)
    c.cd(1)
    maxAllHist[0].Draw("HIST")
    legend.AddEntry(maxAllHist[0], "All", "l")
    for ind, s in enumerate(maxVTPHist):
        s.Draw("HIST SAME")
        legend.AddEntry(s, bitThr[ind], "l")
    legend.Draw()  

    c.cd(2)
    maxRatioHist[0].SetAxisRange(0.9,1.03,"Y")
    maxRatioHist[0].Draw("HIST")
    maxRatioHist[0].Draw("E1 SAME")
    for d in range(1,len(maxRatioHist)):
        maxRatioHist[d].Draw("HIST SAME")

    c.Print(pdfName)
    c.Clear()
    legend.Clear()

     #-----------------------------------------------------------
    c.Divide(1,2)
    c.cd(1)
    dtAllHist[0].SetLineColor(kBlack)
    dtAllHist[0].Draw("HIST")
    legend.AddEntry(dtAllHist[0], "All", "l")
    for ind, s in enumerate(dtVTPHist):
        s.SetLineColor(Bit_Color[ind])
        s.Draw("HIST SAME")
        legend.AddEntry(s, bitThr[ind], "l")
    legend.Draw()  

    c.cd(2)
    dtRatioHist[0].SetLineColor(Bit_Color[0])
    dtRatioHist[0].Draw("HIST")
    dtRatioHist[0].Draw("E1 SAME")
    for d in range(1,len(dtRatioHist)):
        dtRatioHist[d].SetLineColor(Bit_Color[d])
        dtRatioHist[d].Draw("HIST SAME")

    c.Print(pdfName + "(")
    c.Clear()
    legend.Clear()

    #-----------------------------------------------------------
    c.Divide(1,2)
    c.cd(1)
    dtAllHist[0].Draw("HIST")
    legend.AddEntry(dtAllHist[0], "All", "l")
    for ind, s in enumerate(dtVTPHist):
        s.Draw("HIST SAME")
        legend.AddEntry(s, bitThr[ind], "l")
    legend.Draw()  

    c.cd(2)
    dtRatioHist[0].SetAxisRange(0.9,1.03,"Y")
    dtRatioHist[0].Draw("HIST")
    dtRatioHist[0].Draw("E1 SAME")
    for d in range(1,len(dtRatioHist)):
        dtRatioHist[d].Draw("HIST SAME")

    c.Print(pdfName)
    c.Clear()
    legend.Clear()

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
    c.Print(pdfName)
    c.Clear()

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

    c.Divide(1,2)
    c.cd(1)
    timeCorrHist[0].Draw("COLZ")
    c.cd(2)
    timeCorrHist[1].Draw("COLZ")
    c.Print(pdfName + ")")
    c.Clear()


