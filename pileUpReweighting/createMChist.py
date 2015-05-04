import ROOT
import sys

# from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Pileup_MC_Gen_Scenarios
Summer2012_S7 = [
	2.344E-05,
	2.344E-05,
	2.344E-05,
	2.344E-05,
	4.687E-04,
	4.687E-04,
	7.032E-04,
	9.414E-04,
	1.234E-03,
	1.603E-03,
	2.464E-03,
	3.250E-03,
	5.021E-03,
	6.644E-03,
	8.502E-03,
	1.121E-02,
	1.518E-02,
	2.033E-02,
	2.608E-02,
	3.171E-02,
	3.667E-02,
	4.060E-02,
	4.338E-02,
	4.520E-02,
	4.641E-02,
	4.735E-02,
	4.816E-02,
	4.881E-02,
	4.917E-02,
	4.909E-02,
	4.842E-02,
	4.707E-02,
	4.501E-02,
	4.228E-02,
	3.896E-02,
	3.521E-02,
	3.118E-02,
	2.702E-02,
	2.287E-02,
	1.885E-02,
	1.508E-02,
	1.166E-02,
	8.673E-03,
	6.190E-03,
	4.222E-03,
	2.746E-03,
	1.698E-03,
	9.971E-04,
	5.549E-04,
	2.924E-04,
	1.457E-04,
	6.864E-05,
	3.054E-05,
	1.282E-05,
	5.081E-06,
	1.898E-06,
	6.688E-07,
	2.221E-07,
	6.947E-08,
	2.047E-08 ]

Summer2012_S10 = [
	2.560E-06,
	5.239E-06,
	1.420E-05,
	5.005E-05,
	1.001E-04,
	2.705E-04,
	1.999E-03,
	6.097E-03,
	1.046E-02,
	1.383E-02,
	1.685E-02,
	2.055E-02,
	2.572E-02,
	3.262E-02,
	4.121E-02,
	4.977E-02,
	5.539E-02,
	5.725E-02,
	5.607E-02,
	5.312E-02,
	5.008E-02,
	4.763E-02,
	4.558E-02,
	4.363E-02,
	4.159E-02,
	3.933E-02,
	3.681E-02,
	3.406E-02,
	3.116E-02,
	2.818E-02,
	2.519E-02,
	2.226E-02,
	1.946E-02,
	1.682E-02,
	1.437E-02,
	1.215E-02,
	1.016E-02,
	8.400E-03,
	6.873E-03,
	5.564E-03,
	4.457E-03,
	3.533E-03,
	2.772E-03,
	2.154E-03,
	1.656E-03,
	1.261E-03,
	9.513E-04,
	7.107E-04,
	5.259E-04,
	3.856E-04,
	2.801E-04,
	2.017E-04,
	1.439E-04,
	1.017E-04,
	7.126E-05,
	4.948E-05,
	3.405E-05,
	2.322E-05,
	1.570E-05,
	5.005E-06]

def generateHistoFromList( list, name, title="pileup" ):
	histo = ROOT.TH1F( name, title, len(list), 0, len(list) )
	for bin, content in enumerate( list ):
		histo.SetBinContent( bin+1, content ) # bin0 is underflow!
	return histo

def writeObjectsToFile( list, filename ):
	# list contains objects which will be written to a TFile
	outfile = ROOT.TFile( filename, "recreate" )
	outfile.cd()
	for entry in list:
		entry.Write()
	outfile.Close()


if __name__ == "__main__":

	import sys
	if len(sys.argv) < 2:
		sys.exit('Usage: %s outputFilename.root' % sys.argv[0])
	filename = sys.argv[1]

	s10Hist = generateHistoFromList( Summer2012_S10, "pileupScenarioS10" )
	s7Hist = generateHistoFromList( Summer2012_S7, "pileupScenarioS7" )
	writeObjectsToFile( [ s10Hist, s7Hist ], filename )


