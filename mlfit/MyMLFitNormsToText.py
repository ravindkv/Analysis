import re
from sys import argv, stdout, stderr, exit
from optparse import OptionParser
import math

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
import ROOT
ROOT.gROOT.SetBatch(True)

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("-u", "--uncertainties", default=False, action="store_true", help="Report the uncertainties from the fit(s) too")

(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    exit(1)

errors = False
if options.uncertainties: 
    errors = True

file = ROOT.TFile.Open(args[0]);
prefit = file.Get("norm_prefit")
fit_s = file.Get("norm_fit_s")
fit_b = file.Get("norm_fit_b")
if prefit == None: stderr.write("Missing fit_s in %s. Did you run FitDiagnostics in a recent-enough version of combine and with --saveNorm?\n" % file);
if fit_s  == None: raise RuntimeError, "Missing fit_s in %s. Did you run FitDiagnostics with --saveNorm?" % file;
if fit_b  == None: raise RuntimeError, "Missing fit_b in %s. Did you run FitDiagnostics with --saveNorm?" % file;

iter = fit_s.createIterator()
#Headline = "%-30s %-30s     pre-fit   signal+background Fit  bkg-only Fit"%("Channel","Process") if (prefit and errors) else "%-30s %-30s  signal+background Fit  bkg-only Fit"%("Channel","Process")
if prefit and errors :
 headrow  = ["Channel","Process","Pre-fit","S+B Fit","B-Only Fit"]
 headline = ("{:5} {:5} {:^25} {:^25} {:^25}").format(*headrow)
elif prefit: 
 headrow  = ["Channel","Process","Pre-fit","S+B Fit","B-Only Fit"]
 headline = ("{:40} {:25} {:>20} {:>20} {:>20}").format(*headrow)
else : 
 headrow = ["Channel","Process","S+B Fit","B-Only Fit"]
 headline = ("{:40} {:25} {:>20} {:>20}").format(*headrow)

line = "".join(["-" for i in range(len(headline))])
print headline
print line

row_ttbar   = ["SM $\\ttbar$ + jets", "muL", "eL", "muM", "eM", "muT", "eT"]
row_stop    = ["Single \\PQt", "mu L", "e L", "muM", "eM", "muT", "eT"]
row_qcd     = ["QCD multijet", "mu L", "e L", "muM", "eM", "muT", "eT"]
row_wjet    = ["\\PW + jets", "mu L", "e L", "muM", "eM", "muT", "eT"]
row_zjet    = ["$\\PZ/\\PGg$ + jets", "mu L", "e L", "muM", "eM", "muT", "eT"]
row_vv      = ["VV", "mu L", "e L", "muM", "eM", "muT", "eT"]
def modifyAllRow(ch, proc, val, err):
    modifyRow(ch, proc, val, err, "ttbar", row_ttbar)
    modifyRow(ch, proc, val, err, "stop", row_stop)
    modifyRow(ch, proc, val, err, "qcd", row_qcd)
    modifyRow(ch, proc, val, err, "wjet", row_wjet)
    modifyRow(ch, proc, val, err, "zjet", row_zjet)
    modifyRow(ch, proc, val, err, "vv", row_vv)

def modifyRow (ch, proc, val, err, procName, rowName):
    val_ = str(math.ceil(val)).split(".0")[0]
    err_ = str(math.ceil(err)).split(".0")[0]
    if(proc==procName and ch == "ch1"): rowName[1] = val_+" $\\pm$ "+ err_
    if(proc==procName and ch == "ch2"): rowName[3] = val_+" $\\pm$ "+ err_
    if(proc==procName and ch == "ch3"): rowName[5] = val_+" $\\pm$ "+ err_
    if(proc==procName and ch == "ch4"): rowName[2] = val_+" $\\pm$ "+ err_
    if(proc==procName and ch == "ch5"): rowName[4] = val_+" $\\pm$ "+ err_
    if(proc==procName and ch == "ch6"): rowName[6] = val_+" $\\pm$ "+ err_
    return rowName

while True:
    norm_s = iter.Next()
    if norm_s == None: break;
    norm_b = fit_b.find(norm_s.GetName())
    norm_p = prefit.find(norm_s.GetName()) if prefit else None
    # we have to replace any non-standard characters with "_" otherwise the matching will screw up 
    proc_chan_name = (norm_s.GetName()).replace(".","_").replace(":","_").replace(",","_")
    m = re.match(r"(\w+)/(\w+)", proc_chan_name);
    if m == None: m = re.match(r"n_exp_(?:final_)?(?:bin)+(\.\w+)_proc_(\.\w+)", proc_chan_name);
    if m == None: raise RuntimeError, "Non-conforming object name %s" % norm_s.GetName()
    if norm_b == None: raise RuntimeError, "Missing normalization %s for background fit" % norm_s.GetName()
    if prefit and norm_p and errors:
        row = ["%-5s"%m.group(1), "%-5s"%m.group(2), "%10.3f +/- %-10.3f"%(norm_p.getVal(), norm_p.getError()), "%10.3f +/- %-10.3f"%(norm_s.getVal(), norm_s.getError()),"%10.3f +/- %-10.3f"%(norm_b.getVal(), norm_b.getError())]
        #modifyAllRow(m.group(1), m.group(2), norm_p.getVal(), norm_p.getError())
        #modifyAllRow(m.group(1), m.group(2), norm_s.getVal(), norm_s.getError())
        modifyAllRow(m.group(1), m.group(2), norm_b.getVal(), norm_b.getError())
	print("{:<5} {:5} {:10} {:10} {:10}").format(*row)
        #print "%-30s %-30s % 7.3f +/- % 7.3f % 7.3f +/- % 7.3f  % 7.3f +/- % 7.3f" % 
    else:
        if norm_p and prefit:
            row = ["%-40s"%m.group(1), "%-25s"%m.group(2), "%10.3f"%(norm_p.getVal()), "%10.3f"%(norm_s.getVal()),"%10.3f"%(norm_b.getVal())]
	    print("{:<40} {:25} {:>20} {:>20} {:>20}").format(*row)
            #print "%-30s %-30s %7.3f %7.3f %7.3f" % (m.group(1), m.group(2), norm_p.getVal(),  norm_s.getVal(),  norm_b.getVal())
        else:
            row = ["%-40s"%m.group(1), "%-25s"%m.group(2), "%10.3f"%(norm_s.getVal()),"%10.3f"%(norm_b.getVal())]
	    print("{:<40} {:25} {:>20} {:>20}").format(*row)
            #print "%-30s %-30s %7.3f %7.3f" % (m.group(1), m.group(2), norm_s.getVal(), norm_b.getVal())

row_AllBkg = ["All background", "mu L", "e L", "muM", "eM", "muT", "eT"]
def addAllBkg(i):
    y_ttbar = int(row_ttbar[i].split("$\\pm$")[0].strip())
    e_ttbar = int(row_ttbar[i].split("$\\pm$")[1].strip())
    y_stop  = int(row_stop[i].split("$\\pm$")[0].strip())
    e_stop  = int(row_stop[i].split("$\\pm$")[1].strip())
    y_qcd   = int(row_qcd[i].split("$\\pm$")[0].strip())
    e_qcd   = int(row_qcd[i].split("$\\pm$")[1].strip())
    y_wjet  = int(row_wjet[i].split("$\\pm$")[0].strip())
    e_wjet  = int(row_wjet[i].split("$\\pm$")[1].strip())
    y_zjet  = int(row_zjet[i].split("$\\pm$")[0].strip())
    e_zjet  = int(row_zjet[i].split("$\\pm$")[1].strip())
    y_vv    = int(row_vv[i].split("$\\pm$")[0].strip())
    e_vv    = int(row_vv[i].split("$\\pm$")[1].strip())
    yield_bkg = y_ttbar+ y_stop + y_qcd + y_wjet + y_zjet + y_vv
    error_bkg = math.sqrt(e_ttbar**2 + e_stop**2 + e_qcd**2 + e_wjet**2 + e_zjet**2 + e_vv**2)
    return str(math.ceil(yield_bkg)).split(".0")[0] +" $\\pm$ "+ str(math.ceil(error_bkg)).split(".0")[0]
for i in range(6):
    row_AllBkg[i+1] = addAllBkg(i+1)

def printLatexRow(rowName):
    row_ = ""
    max_ = len(rowName) -1
    for i in range(max_):
        row_ = row_ + rowName[i] + " & "
    print row_ + rowName[len(rowName)-1] + " \\\\ "

print "\\hline"
printLatexRow( row_ttbar)
printLatexRow( row_stop )
printLatexRow( row_qcd  )
printLatexRow( row_wjet )
printLatexRow( row_zjet )
printLatexRow( row_vv   )
print "\\hline"
printLatexRow( row_AllBkg)
print "\\hline"
