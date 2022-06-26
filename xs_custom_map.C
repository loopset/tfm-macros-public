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
#include "TSpline.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "RtypesCore.h"
#include "TVirtualPad.h"
#include "Math/Vector3D.h"
//C++ includes
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>//for std::nan and so
#include <iterator>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <algorithm>//reverse and iterator operations

void xs_custom_map(std::string energy="13MeV",
				   std::string ann="-18.55fm_new",
				   std::string particle="n",
				   //std::string drawing="normal",
				   Bool_t rewrite_S=false,
				   Double_t S_step=0.05)
{
	ROOT::EnableImplicitMT();
	//Get energy from string
	const Double_t Tn { energy_from_string(energy)};
	std::cout<<BOLDGREEN<<"Calculation for "<<energy<<" and "<<ann<<RESET<<'\n';

	//Angles configuration file
	std::string angles_file {TString::Format("%s/%s/%s/%s/%s/angles_configuration.dat", directory.c_str(), subdir.c_str(),
							  energy.c_str(), ann.c_str(),
							  particle.c_str()).Data()};
	//nominal angles file
	std::string nominal_angles_file {TString::Format("%s/%s/%s/%s/%s/angles_configuration_nominal.dat", directory.c_str(), subdir.c_str(),
							  energy.c_str(), ann.c_str(),
							  particle.c_str()).Data()};
	//std::vectors with angle configurations
	auto [index, theta1, phi1, theta2, phi2] = read_angles_file(angles_file.c_str());
	//nominal
	auto [index_nom, theta1_nom, phi1_nom, theta2_nom, phi2_nom] = read_angles_file(nominal_angles_file.c_str());

	//TTree for storing output
	auto* fout = new TFile(TString::Format("%s/%s/%s/%s/%s/xs_map.root", directory.c_str(), subdir.c_str(),
							  energy.c_str(), ann.c_str(),
							  particle.c_str()).Data(), "recreate");
	//variables for tree
	Double_t T1, T2p, T2n, sigma;
	Double_t S_max_kin, S_max_theo;
	Double_t S_T1max, S_T1min;
	Double_t angle_t1, angle_t2, angle_p1, angle_p2;
	Double_t nom_t1, nom_t2, nom_p1, nom_p2;
	Int_t region;//tags kinematical region according to S-T1 plot
	auto* t1 = new TTree("Map_xs", "Map_xs");
	t1->Branch("angle_t1", &angle_t1, "angle_t1/D");
	t1->Branch("angle_t2", &angle_t2, "angle_t2/D");
	t1->Branch("angle_p1", &angle_p1, "angle_p1/D");
	t1->Branch("angle_p2", &angle_p2, "angle_p2/D");
	t1->Branch("nom_t1", &nom_t1, "nom_t1/D");
	t1->Branch("nom_t2", &nom_t2, "nom_t2/D");
	t1->Branch("nom_p1", &nom_p1, "nom_p1/D");
	t1->Branch("nom_p2", &nom_p2, "nom_p2/D");
	t1->Branch("T1", &T1, "T1/D");
	t1->Branch("T2p", &T2p, "T2p/D");
	t1->Branch("T2n", &T2n, "T2n/D");
	t1->Branch("sigma", &sigma, "sigma/D");
	t1->Branch("S_max_kin", &S_max_kin, "S_max_kin/D");
	t1->Branch("S_max_theo", &S_max_theo, "S_max_theo/D");
	t1->Branch("S_T1max", &S_T1max, "S_T1max/D");
	t1->Branch("S_T1min", &S_T1min, "S_T1min/D");
	t1->Branch("region", &region, "region/I");

	//TH2 for output
	auto* h = new TH3D("h", "3d map of xs", 100, 0., Tn, 100, 0., Tn, 100, 0., 25.);

	//Count time of execution
	TStopwatch t;
	t.Start(); 
	for (Int_t i=0; i<=index.back(); ++i)
	{
	
		std::cout<<GREEN<<"xs for: "<<i<<" t1: "<< theta1[i] <<" p1: "<<phi1[i] << " t2: "<<theta2[i]
				 <<" p2: "<<phi2[i]<<RESET<<std::endl;

		//Array with kinematic configuration, in RADIANS
		Double_t K[6] {Tn, deg_to_rad(theta1[i]), deg_to_rad(theta2[i]),
					   deg_to_rad(phi1[i]), deg_to_rad(phi2[i]), static_cast<Double_t>((particle == "n") ? true : false)};

		//Get S curve! Integrate if it does not exist (ignore vT2)
		auto [vS, vT1, _] = get_S_file(energy, theta1[i], theta2[i], phi1[i], phi2[i], particle, "primary", rewrite_S);
		//and now TSpline with integrated kinematics
		auto* S_spline = new TSpline3("S_spline", &(vS[0]), &(vT1[0]), vS.size(),"b2,e2", 0, 0);
		auto* S_func = new TF1("S_func", [&](Double_t* x, Double_t* p){return S_spline->Eval(x[0]);},
							   0., S_spline->GetXmax(), 1);
		
		//now theoretical cross section
		//program theo_xs_getter.C must be executed before!
		std::string graph_file {TString::Format("%s/%s/%s/%s/%s/theoretical_xs/xs_%.2f_%.2f_%.2f_%.2f.dat", directory.c_str(), subdir.c_str(),
																	 energy.c_str(), ann.c_str(),
																	 particle.c_str(),
																	 theta1[i], phi1[i], theta2[i], phi2[i]).Data()};
		//abort if it does not exist
		if(!file_exists(graph_file))
		{
			std::cout<<BOLDRED<<"Error opening file: "<<graph_file<<RESET<<std::endl;
			std::abort();
		}
		auto* xs_graph = new TGraphErrors(TString(graph_file));
		auto* xs_spline = new TSpline3("xs_spline", xs_graph->GetX(), xs_graph->GetY(),
									   xs_graph->GetN(), "b2,e2", 0, 0);
		auto* xs_func = new TF1("xs_func", [&](Double_t* x, Double_t* p){ return xs_spline->Eval(x[0]);}, 0., xs_spline->GetXmax(), 1);
		//max value of S fiven by fortran file
		S_max_theo = xs_spline->GetXmax();
		//max value of S given by our integrated file
		S_max_kin = S_spline->GetXmax();
		
		//kinematical points of interest, when branch changes
		std::vector<Double_t>::iterator index_max { std::max_element( vT1.begin(), vT1.end())};
		std::vector<Double_t>::iterator index_min { std::min_element( vT1.begin(), vT1.end())};
	    S_T1max = vS[std::distance(vT1.begin(), index_max)];
		S_T1min = vS[std::distance(vT1.begin(), index_min)];
		//std::cout<<"- -> +: "<<S_T1max<<" + -> -: "<<S_T1min<<" S_max: "<<other_S_max<<'\n';

		//Build vector for points
		//Decide upper limit comparing S_max values. What matters is S_max_kin
		Double_t upper_limit { (S_max_theo >= S_max_kin) ? S_max_kin : S_max_theo};
		std::vector<Double_t> points{ arange(0., upper_limit, S_step)};//consider other step values to reduce size of .root file
		for(auto& s : points)
		{
			sigma = xs_func->Eval(s);
			T1    = S_func->Eval(s);
			if(std::isnan(T1)) std::abort();//if T1 is nan somewhere, we should take a look at
			//recognize branch from S_T1max,min values
			if( s <= S_T1max )//from 1 to 2
			{
				T2n = T2n_calculator(&T1, K);
				h ->Fill(T1, T2n, sigma);
				T2p = 0.;
				region = 1;
			}
			else if( (s > S_T1max) && (s < S_T1min))//from 1 to zero or 4
			{
				T2p = T2p_calculator(&T1, K);
				h ->Fill(T1, T2p, sigma);
				T2n = 0.;
				region = 2;
			}
			else //from 4 or zero to 3
			{
				T2n = T2n_calculator(&T1, K);
				h ->Fill(T1, T2n, sigma);
				T2p = 0.;
				region = 3;				
			}

			//tag angles
			angle_t1 = theta1[i]; angle_t2 = theta2[i];
			angle_p1 = phi1[i]; angle_p2 = phi2[i];

			//nominal angles
			nom_t1 = theta1_nom[i]; nom_t2 = theta2_nom[i];
			nom_p1 = phi1_nom[i];   nom_p2 = phi2_nom[i];
		
			t1->Fill();
		}
		delete xs_func;
		delete S_func;
		delete xs_spline;
		delete xs_graph;
		delete S_spline;
	}
	t.Stop();
	t.Print();

	fout->cd();
	h->Write();
	t1->Write();
	fout->Close();
	
	/*
	auto* c1 {new TCanvas("c1", "3D map of xs", 700,500)};
	c1->cd();
	h->SetMarkerStyle(1);
	h->SetMarkerColor((particle == "n") ? 8 : 6);
	h->Draw("");
	h->GetXaxis()->SetTitle("T_{n1} [MeV]");
	if (particle == "n") h->GetYaxis()->SetTitle("T_{n2} [MeV]");
	else if (particle == "p") h->GetYaxis()->SetTitle("T_{p} [MeV]");
	h->GetZaxis()->SetTitle("xs [b/MeV/sr^{2}]");
	h->GetZaxis()->SetTitleOffset(1.6);
	*/
	

}
