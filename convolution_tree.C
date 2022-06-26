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

void convolution_tree(std::string energy="13MeV",
					  std::string ann="-15.84fm",
					  std::string particle="n",
					  Int_t iterations=30,
					  Double_t step=0.1)
{
	//enable multithreading
	ROOT::EnableImplicitMT();
	//get energy from string
	auto Tn { energy_from_string(energy)};

	//set precision for output
	std::cout<<std::fixed<<std::setprecision(2);
	//file with cross section map
	auto* infile = new TFile(TString::Format("%s/%s/%s/%s/%s/xs_map.root", directory.c_str(), subdir.c_str(),
							  energy.c_str(), ann.c_str(),
							  particle.c_str()).Data(), "read");
	auto* t1 = (TTree*)infile->Get("Map_xs");

	//file with angles configuration
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

	//variables to read TTree
	Double_t T1, T2p, T2n, sigma;
	Double_t theta1, theta2, phi1, phi2;
	Double_t theta1_nom, theta2_nom, phi1_nom, phi2_nom;
	//Double_t S_T1max, S_T1min;
	Int_t region;
	t1->SetBranchAddress("T1", &T1); t1->SetBranchAddress("T2p", &T2p);
	t1->SetBranchAddress("T2n", &T2n); t1->SetBranchAddress("sigma", &sigma);
	t1->SetBranchAddress("angle_t1", &theta1); t1->SetBranchAddress("angle_t2", &theta2);
	t1->SetBranchAddress("angle_p1", &phi1); t1->SetBranchAddress("angle_p2", &phi2);
	t1->SetBranchAddress("nom_t1", &theta1_nom); t1->SetBranchAddress("nom_t2", &theta2_nom);
	t1->SetBranchAddress("nom_p1", &phi1_nom); t1->SetBranchAddress("nom_p2", &phi2_nom);
	//t1->SetBranchAddress("S_T1max", &S_T1max); t1->SetBranchAddress("S_T1min", &S_T1min);
	t1->SetBranchAddress("region", &region);

	//NOMINAL theta angles to make convolution
	std::vector<Double_t> theta1_nominal {20.};
	std::vector<Double_t> theta2_nominal;
	if (particle == "n") theta2_nominal = theta1_nominal;
	else if (particle == "p") theta2_nominal = { 12.};
  	std::vector<Double_t> phi1_nominal { 180.};
	std::vector<Double_t> phi2_nominal { (particle == "n") ? 180. : 90.};

	//configuring convolution
	Double_t fwhm_factor { 1.};//actually, all quated fwhm are sigmas!!
   
	//final histogram binning
	Int_t binning_T { 100};
	Int_t binning_xs { 100};
	Double_t xs_upper_limit { 25.};
	auto* final_convolution = new TH3D("final_convolution", "3d convolution", binning_T, 0., Tn, binning_T, 0., Tn, binning_xs, 0., xs_upper_limit);

	//and new TFile with TTree!
	auto* outfile = new TFile(TString::Format("%s/%s/%s/%s/%s/convolution_tree.root", directory.c_str(), subdir.c_str(),
							  energy.c_str(), ann.c_str(),
							  particle.c_str()).Data(), "recreate");
	outfile->cd();
	auto* convolution_tree = new TTree("convolution_tree", "convolution_tree");
	Double_t T1_tree, T2p_tree, T2n_tree, t1_tree, t2_tree, p1_tree, p2_tree, sigma_tree, S1_tree, S2_tree, S0_tree;
	Double_t t1_nom_tree, t2_nom_tree, p1_nom_tree, p2_nom_tree;
	Int_t region_tree;
	Double_t distance_cluster;
	convolution_tree->Branch("T1_tree", &T1_tree, "T1_tree/D");
	convolution_tree->Branch("T2p_tree", &T2p_tree, "T2p_tree/D");
	convolution_tree->Branch("T2n_tree", &T2n_tree, "T2n_tree/D");
	convolution_tree->Branch("theta1_tree", &t1_tree, "theta1_tree/D");
	convolution_tree->Branch("theta2_tree", &t2_tree, "theta2_tree/D");
	convolution_tree->Branch("phi1_tree", &p1_tree, "phi1_tree/D");
	convolution_tree->Branch("phi2_tree", &p2_tree, "phi2_tree/D");
	convolution_tree->Branch("theta1_nom_tree", &t1_nom_tree, "theta1_nom_tree/D");
	convolution_tree->Branch("theta2_nom_tree", &t2_nom_tree, "theta2_nom_tree/D");
	convolution_tree->Branch("phi1_nom_tree", &p1_nom_tree, "phi1_nom_tree/D");
	convolution_tree->Branch("phi2_nom_tree", &p2_nom_tree, "phi2_nom_tree/D");
	convolution_tree->Branch("sigma_tree", &sigma_tree, "sigma_tree/D");
	convolution_tree->Branch("S0_tree", &S0_tree, "S0_tree/D");
	convolution_tree->Branch("d_cluster", &distance_cluster, "d_cluster/D");
	convolution_tree->Branch("S1_tree", &S1_tree, "S1_tree/D");
	convolution_tree->Branch("S2_tree", &S2_tree, "S2_tree/D");
	convolution_tree->Branch("region_tree", &region_tree, "region_tree/I");

	Int_t counter { 0};
	TStopwatch t; t.Start();
	//and fill new TTrees
	for(auto loop1 : theta1_nominal)
	{
		Int_t theta2_break { 0};
		for(auto loop2 : theta2_nominal)
		{
			if(particle == "n")
			{
				loop2 = loop1;
				if(theta2_break >= 1) break;
			}
			std::cout<<BOLDRED<<"t1 nominal: "<<loop1<<" t2_nominal: "<<loop2<<RESET<<'\n';
			//Check if angles in index file fall inside uncertainty bands around nominal values
			for(Int_t i = 0; i <= (Int_t)index.back(); i++)
			{
				if(!(equal_doubles(loop1, theta1_nom_file[i]) && (equal_doubles(loop2, theta2_nom_file[i])))) continue;
				std::cout<<BOLDGREEN<<"Convolution number "<<counter<< RESET <<'\n';
				std::cout<<"t1: "<<theta1_file[i]<<" p1: "<<phi1_file[i]<<" t2: "<<theta2_file[i]<<" p2: "<<phi2_file[i]<<'\n';

				//functions to convert to S1 and S2
				auto [vS, vT1, vT2] = get_S_file(energy, theta1_file[i], theta2_file[i],
											phi1_file[i], phi2_file[i],
											particle, "primary", false);
				//get maximum and minimum values in S, which determine +/- branch
				//for T1 vs S
				std::vector<Double_t>::iterator index_max { std::max_element( vT1.begin(), vT1.end())};
				std::vector<Double_t>::iterator index_min { std::min_element( vT1.begin(), vT1.end())};
				Double_t S_T1max { vS[std::distance(vT1.begin(), index_max)]};
				Double_t S_T1min { vS[std::distance(vT1.begin(), index_min)]};
				Double_t T1max { TMath::MaxElement(vT1.size(), &(vT1[0]))};
				Double_t T1min { TMath::MinElement(vT1.size(), &(vT1[0]))};
				
				//vector of kinematic configuration (function takes angles as degrees)
				Double_t K[6] {Tn, deg_to_rad(theta1_file[i]), deg_to_rad(theta2_file[i]),
							   deg_to_rad(phi1_file[i]), deg_to_rad(phi2_file[i]),
							   static_cast<Double_t>((particle == "n") ? true : false)};
				//for T2 vs S
				std::vector<Double_t>::iterator index_max2 { std::max_element( vT2.begin(), vT2.end())};
				std::vector<Double_t>::iterator index_min2 { std::min_element( vT2.begin(), vT2.end())};
				Double_t T2_S_T1max { vT2[std::distance(vT1.begin(), index_max)]};
				Double_t T2_S_T1min { vT2[std::distance(vT1.begin(), index_min)]};
				Double_t S_T2max { vS[std::distance(vT2.begin(), index_max2)]};
				Double_t S_T2min { vS[std::distance(vT2.begin(), index_min2)]};
				Double_t T2max { TMath::MaxElement(vT2.size(), &(vT2[0]))};
				Double_t T2min { TMath::MinElement(vT2.size(), &(vT2[0]))};
				//std::cout<<"T2_S_T1max: "<<T2_S_T1max<<'\n';
				auto spline_T1 = new TSpline3("spline_T1", &(vS[0]), &(vT1[0]), vS.size(),"b2,e2", 0, 0);
				auto function_T1  = new TF1("function_T1", [=](Double_t* x, Double_t* p){return spline_T1->Eval(x[0]);}, 0., spline_T1->GetXmax(), 1);
				auto spline_T2 = new TSpline3("spline_T2", &(vS[0]), &(vT2[0]), vS.size(),"b2,e2", 0, 0);
				auto function_T2  = new TF1("function_T2", [=](Double_t* x, Double_t* p){return spline_T2->Eval(x[0]);}, 0., spline_T2->GetXmax(), 1);
				/*
				//3rd particle kinematics to skip ambivalence on T2 vs S
				auto* function_T3 = new TF1("function_T3", T3_calculator, 0., Tn, 6);
				function_T3->SetParameters(K);
				//T3 depends on T1 only, we must obtain first the boundary at T2max -> T1_S_T2max
				Double_t T1_T2max { vT1[std::distance(vT2.begin(), index_max2)]};
				Double_t limit { function_T3->Eval(T1_T2max)};
				std::cout<<"T1_T2max: "<<T1_T2max<<" limit: "<<limit<<'\n';
				std::cout<<"Spline T1 min: "<<function_T1->GetMinimumX()<<" Spline T2 max: "<<function_T2->GetMaximumX()<<'\n';
				std::cout<<"S_T1min: "<<S_T1min<<'\n';
				*/
				function_T2->SetNpx(150);
				//loop over entries
				for(Int_t j = 0; j < (Int_t)t1->GetEntries(); j++)
				{
					t1->GetEntry(j);
					if(!( (equal_doubles(theta1, theta1_file[i])) && (equal_doubles(theta2, theta2_file[i]))
						  && (equal_doubles(phi1, phi1_file[i])) && (equal_doubles(phi2, phi2_file[i])) )) continue;
					//reassure we are in nominal values
					if(!( (equal_doubles(theta1_nom, theta1_nom_file[i])) && (equal_doubles(theta2_nom, theta2_nom_file[i]))
						  && (equal_doubles(phi1_nom, phi1_nom_file[i])) && (equal_doubles(phi2_nom, phi2_nom_file[i])) )) { std::cout<<BOLDRED<<"Issue tagging values on tree"<<RESET<<'\n'; std::abort();}

					//and perform # iterations convolutions
					for(Int_t it = 0; it < iterations; it++)
					{
						//make convolution here
						Double_t fwhm_x { (it == 0) ? 0. : monster_resolution1->Eval(T1) * T1 / fwhm_factor};
						T1_tree = T1 + gRandom->Gaus(0., fwhm_x);

						//distinguish between regions
						if(region == 1)
						{
							//Double_t fwhm_y { (it == 0) ? 0. : monster_resolution2->Eval(T2n) * T2n / fwhm_factor};
							Double_t fwhm_y;
							if(it == 0) fwhm_y = 0.;
							else if (particle == "n") fwhm_y = monster_resolution2->Eval(T2n) * T2n / fwhm_factor;
							else if (particle == "p") fwhm_y = actar_sigma->Eval(T2n);
							else { std::cout<<BOLDRED<<"Error calculating sigma y"<<RESET<<'\n'; std::abort();}
							
							T2n_tree = T2n + gRandom->Gaus(0., fwhm_y);
							T2p_tree = 0.;

							//check if T1 is smaller than T1max
							if((T1_tree <= T1max) && (T1_tree >= vT1.front())) S1_tree = function_T1->GetX(T1_tree, 0., S_T1max);
							else S1_tree = -20.;

							//for T2n
							if((T2n_tree <= T2_S_T1max) && (T2n_tree >= T2min)) S2_tree = function_T2->GetX(T2n_tree, 0., S_T1max);
							else S2_tree = -40.;

							region_tree = region;
						}
						else if (region == 2)
						{
							//Double_t fwhm_y { (it == 0) ? 0. : monster_resolution2->Eval(T2p) * T2p / fwhm_factor};
							Double_t fwhm_y;
							if(it == 0) fwhm_y = 0.;
							else if (particle == "n") fwhm_y = monster_resolution2->Eval(T2p) * T2p / fwhm_factor;
							else if (particle == "p") fwhm_y = actar_sigma->Eval(T2p);
							else { std::cout<<BOLDRED<<"Error calculating sigma y"<<RESET<<'\n'; std::abort();}
							
							T2p_tree = T2p + gRandom->Gaus(0., fwhm_y);
							//std::cout<<"fwhm_y: "<<fwhm_y<<'\n';
							T2n_tree = 0.;

							//check if T1 is smaller than T1max
							if((T1_tree <= T1max) && (T1_tree >= T1min)) S1_tree = function_T1->GetX(T1_tree, S_T1max, S_T1min);
							else S1_tree = -20.;

							//for T2p we have an ambiguity
							Double_t root1, root2;
							if((T2p_tree > T2_S_T1max) && (T2p_tree < T2max))
							{
								root1 = function_T2->GetX(T2p_tree, S_T1max, S_T2max);
								S2_tree = root1;
								//let's skip the ambiguity by now""
								// //but if we are in the narrow band around the maximum...
								// if((T2p_tree <= T2max) && (T2p_tree >= T2_S_T1min))
								// {
								// 	root2 = function_T2->GetX(T2p_tree, S_T2max, S_T1min);
								// 	S2_tree = (gRandom->Uniform(0., 1.) >= 0.5) ? root1 : root2;
								// }
								// else S2_tree = root1;
							}
							else S2_tree = -40.;

							region_tree = region;
						}
						else if (region == 3)
						{
							//Double_t fwhm_y { (it == 0) ? 0. : monster_resolution2->Eval(T2n) * T2n / fwhm_factor};
							Double_t fwhm_y;
							if(it == 0) fwhm_y = 0.;
							else if (particle == "n") fwhm_y = monster_resolution2->Eval(T2n) * T2n / fwhm_factor;
							else if (particle == "p") fwhm_y = actar_sigma->Eval(T2n);
							else { std::cout<<BOLDRED<<"Error calculating sigma y"<<RESET<<'\n'; std::abort();}
							
							T2n_tree = T2n + gRandom->Gaus(0., fwhm_y);
							T2p_tree = 0.;

							//check if T1 is smaller than T1max
							if((T1_tree <= vT1.back()) && (T1_tree > T1min)) S1_tree = function_T1->GetX(T1_tree, S_T1min, spline_T1->GetXmax());
							else S1_tree = -20.;

							//for T2n
							if((T2n_tree > vT2.back()) && (T2n_tree < T2_S_T1min)) S2_tree = function_T2->GetX(T2n_tree, S_T1min, spline_T1->GetXmax());
							else S2_tree = -40.;

							region_tree = region;
						}
						else
						{
							std::cout<<"Issue when selecting branch on convolution"<<'\n';
							std::abort();
						}

						//and here S0 projection
						std::vector<Double_t> S_step { arange(vS.front(), vS.back(), step, true)};
						std::vector<Double_t> aux_step;
						for(auto& s : S_step)
						{
							Double_t T1_step { function_T1->Eval(s)};
							Double_t T2_step { function_T2->Eval(s)};	
							Double_t distance { S_distance(T1_step, T2_step, T1_tree, T2p_tree + T2n_tree)};
							//if(!(distance <= step / 2)) continue;
							aux_step.push_back(distance);
						}
						Int_t index_min_step { static_cast<Int_t>(TMath::LocMin(aux_step.size(), &(aux_step[0])))};
						S0_tree = S_step[index_min_step];
						distance_cluster = aux_step[index_min_step];
						
						sigma_tree = sigma;

						//fill histogram
						final_convolution->Fill(T1_tree, T2p_tree + T2n_tree, sigma_tree);
						//mark angles
						t1_tree = theta1_file[i]; t2_tree = theta2_file[i];
						p1_tree = phi1_file[i]; p2_tree = phi2_file[i];
						//again nominal angles
						t1_nom_tree = theta1_nom_file[i]; t2_nom_tree = theta2_nom_file[i];
						p1_nom_tree = phi1_nom_file[i];   p2_nom_tree = phi2_nom_file[i];
						
						convolution_tree->Fill();
					}				
				}
				//if(counter == 10) break;
				//delete function_T3;
				delete function_T2;
				delete function_T1;
				delete spline_T2;
				delete spline_T1;
				counter++;
			}
			//variable that makes t1=t2 for neutrons
			theta2_break++;
		}
	}
	convolution_tree->Write();
	final_convolution->Write();
	outfile->Close();
	t.Stop(); t.Print();
	infile->Close();
	

}
