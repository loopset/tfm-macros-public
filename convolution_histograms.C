//our header include
#include "S_header.C"

//ROOT includes
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TColor.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THStack.h"
#include "TSpline.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "RtypesCore.h"
#include "TVirtualPad.h"
#include "Math/Vector3D.h"
#include "TGraphAsymmErrors.h"
//C++ includes
#include <Rtypes.h>
#include <TRandom.h>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>//for std::nan and so
#include <iterator>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <algorithm>//reverse and iterator operations
#include <cstdlib>


void convolution_histograms(std::string energy="13MeV",
					  std::string ann="-15.84fm",
					  std::string particle="n")
{
	ROOT::EnableImplicitMT();
	//customizing output
	std::cout<<std::fixed<<std::setprecision(2);
	//get energy from string
	auto Tn { energy_from_string(energy)};

	//files with angles configuration
	std::string angles_file {TString::Format("%s/%s/%s/%s/%s/angles_configuration.dat", directory.c_str(), subdir.c_str(),
											 energy.c_str(), ann.c_str(),
											 particle.c_str()).Data()};
	auto [index, theta1_file, phi1_file, theta2_file, phi2_file] = read_angles_file(angles_file.c_str());
	//and nominal ones
	//nominal angles file
	std::string nominal_angles_file {TString::Format("%s/%s/%s/%s/%s/angles_configuration_nominal.dat", directory.c_str(), subdir.c_str(),
													 energy.c_str(), ann.c_str(),
													 particle.c_str()).Data()};
	auto [index_nom, theta1_nom_file, phi1_nom_file, theta2_nom_file, phi2_nom_file] = read_angles_file(nominal_angles_file.c_str());

	//file with convolutions TTree
	auto* infile = new TFile(TString::Format("%s/%s/%s/%s/%s/convolution_tree.root", directory.c_str(), subdir.c_str(),
							  energy.c_str(), ann.c_str(),
							  particle.c_str()).Data(), "read");
	auto* convolution_tree = (TTree*)infile->Get("convolution_tree");
	Double_t T1_tree, T2_tree, t1_tree, t2_tree, p1_tree, p2_tree, sigma_tree, S0_tree, S1_tree, S2_tree;
	Double_t t1_nom_tree, t2_nom_tree, p1_nom_tree, p2_nom_tree;
	convolution_tree->SetBranchAddress("T1_tree", &T1_tree);
	convolution_tree->SetBranchAddress("T2p_tree", &T2_tree);
	convolution_tree->SetBranchAddress("theta1_tree", &t1_tree);
	convolution_tree->SetBranchAddress("theta2_tree", &t2_tree);
	convolution_tree->SetBranchAddress("phi1_tree", &p1_tree);
	convolution_tree->SetBranchAddress("phi2_tree", &p2_tree);
	convolution_tree->SetBranchAddress("theta1_nom_tree", &t1_nom_tree);
	convolution_tree->SetBranchAddress("theta2_nom_tree", &t2_nom_tree);
	convolution_tree->SetBranchAddress("phi1_nom_tree", &p1_nom_tree);
	convolution_tree->SetBranchAddress("phi2_nom_tree", &p2_nom_tree);
	convolution_tree->SetBranchAddress("sigma_tree", &sigma_tree);
	convolution_tree->SetBranchAddress("S0_tree", &S0_tree);
	convolution_tree->SetBranchAddress("S1_tree", &S1_tree);
	convolution_tree->SetBranchAddress("S2_tree", &S2_tree);
	
	//NOMINAL theta angles to make convolution
	std::vector<Double_t> theta1_nominal {20.};
	std::vector<Double_t> theta2_nominal;
	if (particle == "n") theta2_nominal = theta1_nominal;
	else if (particle == "p") theta2_nominal = { 12.};
  	Double_t phi1_nominal { 180.};
	Double_t phi2_nominal { (particle == "n") ? 180. : 90.};

	//final histograms
	std::vector<TGraph*> theoretical_graphs;
	std::vector<TH2D*> final_maps0, final_maps1, final_maps2;
	//parameters of histograms
	Int_t binning_T { 120};
	Int_t binning_xs{ 100};
	Double_t xs_upper_limit { 25.};
	Double_t S_upper_limit { 17.};
	std::cout<<BOLDMAGENTA<<"Bin width: "<<S_upper_limit / binning_T<<RESET<<'\n';

	TStopwatch t; t.Start();
	//and fill new histograms
	Int_t counter { 0};
	for(auto loop1 : theta1_nominal)
	{
		Int_t theta2_break { 0};
		for(auto loop2 : theta2_nominal)
		{
			//force t1 = t2 only for neutrons
			if(particle == "n")
			{
				loop2 = loop1;
				if(theta2_break >= 1) break;
			}
			//nominal values
			std::cout<<BOLDRED<<"t1 nominal: "<<loop1<<" t2 nominal: "<<loop2<<RESET<<'\n';
			//run over values with desired nominal values
			final_maps0.push_back(new TH2D("", "xs on S0_tree", binning_T, 0., S_upper_limit, binning_xs, 0., xs_upper_limit));
			final_maps1.push_back(new TH2D("", "xs on S1_tree", binning_T, 0., S_upper_limit, binning_xs, 0., xs_upper_limit));
			final_maps2.push_back(new TH2D("", "xs on S2_tree", binning_T, 0., S_upper_limit, binning_xs, 0., xs_upper_limit));

			//iterate over entries
			auto entries { convolution_tree->GetEntries()};
			for(Int_t j = 0; j < (Int_t)entries; j++)
			{
				convolution_tree->GetEntry(j);
				if(j % (Int_t)100000 == 0)std::cout<<BOLDGREEN<< static_cast<Double_t>(j) / entries * 100<<" % read"<<RESET<<'\n';

				//check if nominal is ok
				if(!( (equal_doubles(t1_nom_tree, loop1)) &&
					  (equal_doubles(t2_nom_tree, loop2)) )) continue;
				//std::cout<<"Here"<<'\n';
				final_maps0[counter]->Fill(S0_tree, sigma_tree);
				final_maps1[counter]->Fill(S1_tree, sigma_tree);
				final_maps2[counter]->Fill(S2_tree, sigma_tree);
			
			}
		
			//and get theoretical xs
			std::string graph_file { TString::Format("%s/%s/%s/%s/%s/theoretical_xs/xs_%.2f_%.2f_%.2f_%.2f.dat", directory.c_str(), subdir.c_str(),
													 energy.c_str(), ann.c_str(),
													 particle.c_str(),
													 loop1, phi1_nominal, loop2, phi2_nominal).Data()};
			//get theoretical xs file if it does not exist
			if (!file_exists(graph_file))
			{
				gSystem->Exec(TString::Format("%s/xs_on_demand.sh -a %.2f -b %.2f -c %.2f -d %.2f -e %s -f %s -g %s",
											  directory.c_str(),
											  loop1, phi1_nominal, loop2, phi2_nominal,
											  particle.c_str(),
											  ann.c_str(), energy.c_str()));
			}
			theoretical_graphs.push_back(new TGraph(TString(graph_file)));
			theta2_break++;
			counter++;
		}
	}

	//get experimental points, from binning of histograms
	//quantile information
	//number of quantiles
	Int_t nq { 5};
	//1-sigma bands
	Double_t xq[5] { 0.0, 0.32, 0.50, 0.68, 1.00};
	Int_t lower { 1}; Int_t mean { 2}; Int_t upper { 3};
	std::vector<std::vector<TH2D*>> maps { final_maps0, final_maps1, final_maps2};
	
	auto* outfile = new TFile(TString::Format("%s/%s/%s/%s/%s/graphs_and_histos.root", directory.c_str(), subdir.c_str(),
							  energy.c_str(), ann.c_str(),
							  particle.c_str()).Data(), "recreate");
	outfile->cd();

	Int_t g_index { 0};
	for(auto& m : maps)
	{
		for(Int_t i = 0; i < counter; i++)
		{
			auto* push_graph = new TGraphAsymmErrors();
			Int_t aux { 0};
			Int_t point{ 0};
			for(Int_t bin = 1; bin < m[i]->GetNbinsX(); bin++)//exclude overflow bins
			{
				//project each bin!
				auto* projection = (TH1D*)m[i]->ProjectionY(TString::Format("%d_%d", i, bin), bin, bin);
				if(projection->GetEntries() <= 0.)
				{
					aux++;
					continue;
				}
				Double_t yq[nq];
				projection->GetQuantiles(nq, yq, xq);
				Double_t X { m[i]->GetXaxis()->GetBinCenter(bin)};
				//Double_t Y { yq[mean]};
				Double_t Y { projection->GetMean()};
				Double_t EY { projection->GetStdDev()};
				Double_t EYLow { EY};//{ std::abs(yq[lower] - Y)};
				Double_t EYUp  { EY};//{ std::abs(yq[upper] - Y)};
				push_graph->SetPoint(point, X, Y);
				push_graph->SetPointEYlow(point, EYLow);
				push_graph->SetPointEYhigh(point, EYUp);

				point++;
				delete projection;
			}
			std::cout<<BOLDRED<<"Skipped "<<aux<<" projections"<<RESET<<'\n';
			push_graph->Write(TString::Format("map%d_%d", g_index, i));
			delete push_graph;
		}
		g_index++;
	}

	for(Int_t i = 0; i < counter; i++)
	{
		final_maps0[i]->Write(TString::Format("histo0_%d", i));
		final_maps1[i]->Write(TString::Format("histo1_%d", i));
		final_maps2[i]->Write(TString::Format("histo2_%d", i));
		theoretical_graphs[i]->Write(TString::Format("theoretical_%d", i));
	}
	
	outfile->Close();
	infile->Close();

	t.Stop(); t.Print();
}
