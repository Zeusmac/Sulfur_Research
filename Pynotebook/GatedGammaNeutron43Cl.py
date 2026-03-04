import ROOT as rt
import ctypes


#44SGated
class DecayFit:

	def __init__(self,X0,bounds,hist,fitgraph):
		self.X0 = X0
		self.bounds = bounds
		self.hist = hist
		self.figraph = fitgraph

	def decay_parent(x, par):
		"""
		a = amplitude
		b = decay constant (lambda)
		bkgrd = constant background
		"""
		a = par[0]
		b = par[1]
		bkgrd = par[2]
		t = x[0]
		return a * rt.TMath.Exp(-b * t) + bkgrd

	def fit(func, self):
    
		MaxX = hist.GetXaxis().GetXmax();
		#MinX = hist.GetXaxis().GetXmin();
		# Axis range
		MinX = 1
		#MaxX = 
		print("Min X: ", MinX)
		print("Max X: ", MaxX)

		# Canvas
		#fitgraph = rt.TCanvas("fits", "fits", 800, 600)
		fitgraph.Clear()

		# Define TF1
		decayp = rt.TF1("decayp", func, MinX, MaxX, len(self.X0))

		decayp.SetParName(0, "source #")
		decayp.SetParName(1, "Source lambda")
		decayp.SetParName(2, "shift")


		for i in range(0,len(X0)):# intial conditions
			decayp.SetParameter(i, self.X0[i])

		for i in range(0, len(bounds)): #Fit bounds
			decayp.SetParLimits(i, self.bounds[i][0], self.bounds[i][1]) 

		# Print initial parameters
		for i in range(decayp.GetNpar()):
			print(f"Initial p{i} = {decayp.GetParameter(i)}")

		# Perform fit
		hist.Fit("decayp", "RS", "", MinX, MaxX)
		#hist.Print("V")

		# Print parameter limits
		for i in range(decayp.GetNpar()):
			par_min = ctypes.c_double()
			par_max = ctypes.c_double()
			decayp.GetParLimits(i, par_min, par_max)
			print(f"Parameter {i} ({decayp.GetParName(i)}) limits: [{par_min}, {par_max}]")

		# Half-life calculation
		lambda_fit = decayp.GetParameter(1)
		lambda_err = decayp.GetParError(1)

		decayhl = rt.TMath.Log(2) / lambda_fit
		decay_error = (rt.TMath.Log(2) / (lambda_fit**2)) * lambda_err

		print(f"Decay: {decayhl:.6f} ms ± {decay_error:.6f} ms\n")

		# Chi2 / NDF
		chi2 = decayp.GetChisquare()
		ndf = decayp.GetNDF()
		print(f"Chi^2/NDF: {chi2/ndf:.6f}\n")

		# Constant background function
		bckgrnd = rt.TF1("bckgrnd", "[0]", MinX, MaxX)
		bckgrnd.SetParameter(0, decayp.GetParameter(2))
		bckgrnd.SetLineColor(rt.kBlack)

		# Legend
		legend = rt.TLegend(0.7, 0.4, 0.9, 0.7)
		legend.AddEntry(hist, "Data", "E")
		legend.AddEntry(decayp, "Fit", "l")
		legend.AddEntry(bckgrnd, "backgrnd", "l")

		# Styling
		hist.SetTitle("Decay of 44S Gated on 329keV(Neutron);Time (ms);Counts")
		hist.GetXaxis().CenterTitle(True)
		hist.GetYaxis().CenterTitle(True)
		hist.SetLineColor(rt.kBlue)
		decayp.SetLineColor(rt.kRed)

		# Draw
		hist.Draw("E")
		decayp.Draw("same")
		bckgrnd.Draw("same")
		fitgraph.Update()
		fitgraph.Draw()
		hist.Print("V")
		# legend.Draw()  # Uncomment if desired

		fitgraph.SaveAs("GammaN43Cl.png")

###end of functions###

file = rt.TFile("e21062_44S.root")
dtimeN = file.Get("dtimeN")

bounds = [[0,100], # source
          [0,1],   #lamb
          [0,50]] #background shift

X0 = [50,
      .005,
      20]

gamma_canvas = rt.TCanvas("gamma_canvas_1", "gamma_329KeV")
gamma_329 = dtimeN.ProjectionX("Decay_Curve_329KeV",327,331)
#gamma_329.Draw()
f = DecayFit(X0,bounds,gamma_329,gamma_canvas)
f.Decayp(decay_parent)
#gamma_canvas.Update()
gamma_canvas.Draw()

