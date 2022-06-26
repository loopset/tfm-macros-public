#include "S_header.C"

#include <RtypesCore.h>
#include <ios>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TF1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TSpline.h"
#include "TStyle.h"
#include <iomanip>
/* Function that runs bash script to get
   theoretical cross section 
*/


void theo_xs_getter(std::string energy="13MeV",
					std::string ann="-15.84fm",
					std::string particle="n",
					Bool_t force_rewrite=false)
{
	//steps for running uncertainty
	Double_t step_theta_neutron { 1.};
	Double_t step_theta_proton  { 0.5};
	Double_t u_theta1 { deltan1};
	Double_t u_theta2 { (particle=="n") ? deltan2 : deltap};
	
	std::vector<Double_t> theta1_nominal {20.};
	std::vector<Double_t> theta2_nominal;
	if (particle == "n") theta2_nominal = theta1_nominal;
	else if (particle == "p") theta2_nominal = { 12.};//arange(thetap_min, thetap_max, 1.);
	
	std::vector<Double_t> phi1_nominal { 180.};
	std::vector<Double_t> phi2_nominal { (particle == "n") ? 180. : 90.};


	//Configurations file
	std::ofstream fout;
	std::ofstream fout_nominal;
	fout.open(TString::Format("%s/%s/%s/%s/%s/angles_configuration.dat", directory.c_str(), subdir.c_str(),
							  energy.c_str(), ann.c_str(),
							  particle.c_str()).Data());
	fout_nominal.open(TString::Format("%s/%s/%s/%s/%s/angles_configuration_nominal.dat", directory.c_str(), subdir.c_str(),
							  energy.c_str(), ann.c_str(),
							  particle.c_str()).Data());

	//count time
	TStopwatch t;
	t.Start();
	std::cout<<std::fixed<<std::setprecision(2);//fix printing precision
	Int_t index { 0};
	if(particle == "n")
	{
		std::cout<<"Entering calculation for n interaction"<<'\n';
		std::cout<<"with uncertainties u(1): "<<deltan1<<" and u(2): "<<deltan2<<'\n';
		for(auto& nom : theta1_nominal)
		{
			//Theta_1 and Theta_2 are the same in this approach!!!!!! (for n_n interaction only)
			std::vector<Double_t> theta1 {arange(nom - u_theta1, nom + u_theta1, step_theta_neutron, true)};
			std::vector<Double_t> theta2 {arange(nom - u_theta2, nom + u_theta2, step_theta_neutron, true)};

			for(auto& t1 : theta1)
			{
				for(auto& t2 : theta2)
				{
					for(auto& phi : phi1_nominal)
					{
						std::vector<Double_t> phi1 { phi + phi_uncertainty(t1, nom, R, d1), phi - phi_uncertainty(t1, nom, R, d1)};
						std::vector<Double_t> phi2 { phi + phi_uncertainty(t2, nom, R, d2), phi - phi_uncertainty(t2, nom, R, d2)};
						for(auto& p1 : phi1)
						{
							for(auto& p2 : phi2)
							{
								std::string xs_file {TString::Format("%s/%s/%s/%s/%s/theoretical_xs/xs_%.2f_%.2f_%.2f_%.2f.dat", directory.c_str(), subdir.c_str(),
																	 energy.c_str(), ann.c_str(),
																	 particle.c_str(),
																	 t1, p1, t2, p2).Data()};
								if(!file_exists(xs_file) | force_rewrite)
								{
									cout<<index<<"    "<<t1<<"    "<<p1<<"    "<<t2<<"    "<<p2<<'\n';
									gSystem->Exec(TString::Format("%s/xs_on_demand.sh -a %.2f -b %.2f -c %.2f -d %.2f -e %s -f %s -g %s",
																  directory.c_str(),
																  t1, p1, t2, p2,
																  particle.c_str(),
																  ann.c_str(), energy.c_str()));
								}
								fout << std::fixed << std::setprecision(2) << std::endl;//ensuring that doubles are written with 2 decimals: very important!
								fout_nominal << std::fixed << std::setprecision(2) << std::endl;
								fout<<index<<"    "<<t1<<"    "<<p1<<"    "<<t2<<"    "<<p2<<'\n';
								fout_nominal<<index<<"    "<<nom<<"    "<<phi<<"    "<<nom<<"    "<<phi<<'\n';
								++index;
							}
						}
					}
				}
			}
		}
	}
	else if (particle == "p")
	{
		std::cout<<"Entering calculation for p interaction"<<'\n';
		std::cout<<"with uncertainties u(1): "<<deltan1<<" and u(2): "<<deltap<<'\n';
		for(auto& nom1 : theta1_nominal)
		{
			//Theta_1 is measured at one of the neutron detectors
			std::vector<Double_t> theta1 { arange(nom1 - u_theta1, nom1 + u_theta1, step_theta_neutron, true)};
			for(auto& nom2 : theta2_nominal)
			{
				//Theta_2 is the proton measured at ACTAR TPC <-- Now measured in all combinations of theta1 vs theta2!!
				std::vector<Double_t> theta2 { arange(nom2 - u_theta2, nom2 + u_theta2, step_theta_proton, true)};
				for(auto& t1 : theta1)
				{
					for(auto& t2 : theta2)
					{
						for(auto& nomp1 : phi1_nominal)
						{
							std::vector<Double_t> phi1 { nomp1 + phi_uncertainty(t1, nom1, R, d1), nomp1 - phi_uncertainty(t1, nom1, R, d1)};
							for(auto& nomp2 : phi2_nominal)
							{
								std::vector<Double_t> phi2 { nomp2 + deltap, nomp2 - deltap};
								for(auto& p1 : phi1)
								{
									for(auto& p2 : phi2)
									{
										std::string xs_file {TString::Format("%s/%s/%s/%s/%s/theoretical_xs/xs_%.2f_%.2f_%.2f_%.2f.dat", directory.c_str(), subdir.c_str(),
																	 energy.c_str(), ann.c_str(),
																	 particle.c_str(),
																			 t1, p1, t2, p2).Data()};
										if(!file_exists(xs_file) | force_rewrite)
										{
											cout<<index<<"    "<<t1<<"    "<<p1<<"    "<<t2<<"    "<<p2<<'\n';
											gSystem->Exec(TString::Format("%s/xs_on_demand.sh -a %.2f -b %.2f -c %.2f -d %.2f -e %s -f %s -g %s",
																  directory.c_str(),
																  t1, p1, t2, p2,
																  particle.c_str(),
																  ann.c_str(), energy.c_str()));
										}
										//very important to fix precision while writing to files
										fout << std::fixed << std::setprecision(2) << std::endl;
										fout_nominal << std::fixed << std::setprecision(2) << std::endl;
										fout<<index<<"    "<<t1<<"    "<<p1<<"    "<<t2<<"    "<<p2<<'\n';
										fout_nominal<<index<<"    "<<nom1<<"    "<<nomp1<<"    "<<nom2<<"    "<<nomp2<<'\n';
										++index;
									}
								}
							}
						}
					}
				}
			}
		}
	}

	t.Stop();
	t.Print();
	fout.close(); fout_nominal.close();
}
