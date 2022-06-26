#ifndef S_HEADER
#define S_HEADER
//ROOT includes
#include "TSystem.h"
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
#include "TH2D.h"
#include "TH3D.h"
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
#include <sys/stat.h>
#include <string>
#include <vector>
#include <algorithm>//reverse and iterator operations

//colored output :)
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

//struct types to pass to functions because ROOT does not like
//C++ tuples, pairs and so.
struct two_ints
{
	Int_t a;
	Int_t b;
};
	
struct two_doubles
{
	Double_t a;
	Double_t b;
};

struct two_vectors
{
	std::vector<Double_t> v1;
	std::vector<Double_t> v2;
};

struct three_vectors
{
	std::vector<Double_t> v1;
	std::vector<Double_t> v2;
	std::vector<Double_t> v3;
};

struct four_vectors{
	std::vector<Double_t> v1;
	std::vector<Double_t> v2;
	std::vector<Double_t> v3;
	std::vector<Double_t> v4;
};

struct five_vectors
{
	std::vector<Double_t> index;
	std::vector<Double_t> theta1;
	std::vector<Double_t> phi1;
	std::vector<Double_t> theta2;
	std::vector<Double_t> phi2;
};
//LOCATE THE CURRENT PATH: S_header.C must be with the major programs on the same directory (and level)
std::string directory { gSystem->pwd()};
std::string subdir    { "scattering_lengths"};
std::string subdirS   { "arc_lengths"};

//CONVERT std::string with energy to Double_t
Double_t energy_from_string(std::string energy)
{
	//energy must be of the form "13MeV"
	auto index { energy.find("M")};
	auto out   { energy.erase(index)};
	Double_t val { std::atof(out.c_str())};
	return val;
}
//TEMPLATES FOR NUMPY-LIKE OPERATIONS ON VECTORS
template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

template <typename T>
std::vector<T> arange(T start, T stop, T step = 1., Bool_t last=false)
{
	std::vector<T> out;
	for (T value = start; (last == false) ?  value < stop : value <= stop; value += step)
	{//does NOT include last value of range!
		out.push_back(value);
	}
	return out;
}

//UTILITY PHYSICAL FUNCTIONS: angle/rad conversion, comparing doubles...

Double_t deg_to_rad (Double_t x)
{
	return x*TMath::Pi()/180.;
}

Double_t rad_to_deg (Double_t x)
{
	return x*180./TMath::Pi();
}

Bool_t equal_doubles(Double_t first, Double_t second, Double_t tol = 1.E-6)
{
    return std::abs(first - second) < tol;
}

//CHECK IF FILE EXISTS
//from: https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exists-using-standard-c-c11-14-17-c
Bool_t file_exists (const std::string& name) {
	//struct stat   buffer;   
	//return (stat (name.c_str(), &buffer) == 0);
	//use ROOT capabilities to check if file exists
	//AccessPathName returns FALSE if file EXISTS
	return !(gSystem->AccessPathName(name.c_str()));
}


//KINEMATIC FORMULAS
//but first constants
//deuterium
const Double_t Md {2.01410177784 * 931.494 - 0.510998950}; //MeV
//proton = m2,m3
const Double_t Mh { 1.007276466 * 931.494}; //MeV
//neutron = mn, m1, m2
const Double_t Mn {1.0086649159 * 931.494};  //MeV
//Q value
const Double_t Q {Mn +Md -2.*Mn -Mh}; //MeV

//UNCERTAINTIES IN ANGLES
//Distances to neutron detectors
const Double_t R {10}; //radius of MONSTER cells, cm
const Double_t d1 {150}; //distance to first detector, cm
const Double_t d2 {250}; //distance to second cell, cm
//uncertainties in angles for NEUTRON
const Double_t deltan1 {3.81}; //degree
const Double_t deltan2 {2.29}; //degree
//uncertainties in angles for PROTON
const Double_t thetap_min { rad_to_deg(TMath::ATan(2. / 12.5))};//nearest angle
const Double_t thetap_max { thetap_min + 5.};
const Double_t deltap { 0.5}; //degree

//uncertainty in phi angles
Double_t phi_uncertainty(Double_t theta, Double_t theta_nom, Double_t R, Double_t d)
{
	//Implicity convert to radians!
	//Outputs degrees!
	
	Double_t val {TMath::ATan(TMath::Sqrt(TMath::Power(R/d, 2) -
										  TMath::Power(TMath::Tan(deg_to_rad(theta-theta_nom)),2)))};
	return rad_to_deg(val);
}


Double_t T2p_calculator(Double_t* x, Double_t* p)
{
	Double_t T1 = x[0];//T1 is the independent variable
	Double_t Tn = p[0], theta1 = p[1], theta2 = p[2], phi1 = p[3], phi2 = p[4];
	Double_t particle = p[5];
	Double_t phi { phi1 - phi2};//sign criterium for Delta phi
	//assigning masses
	Double_t m1 { Mn}; Double_t m2 { (particle == true) ? Mn : Mh}; Double_t m3 { (particle == true) ? Mh : Mn};
	Double_t A {m2 + m3};
	Double_t cos_theta12 {TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Cos(phi)+TMath::Cos(theta1)*TMath::Cos(theta2)};
	Double_t B {2.*cos_theta12*TMath::Sqrt(m1*m2*T1) - 2.*TMath::Cos(theta2)*TMath::Sqrt(Mn*m2*Tn)};
	Double_t C {(m1 + m3)*T1 - m3*Q +Tn*(Mn - m3) -2.*TMath::Cos(theta1)*TMath::Sqrt(Mn*m1*Tn*T1)};
	Double_t D {B*B-4.*A*C};
	if(D>=0.)
	{
		Double_t Tpos {(-B+TMath::Sqrt(D))/(2.*A)};
		if(Tpos>=0.) return Tpos*Tpos;
		else return 0.;
	}
	else return 0.;
}

Double_t T2n_calculator(Double_t* x, Double_t* p)
{
	Double_t T1 = x[0];//T1 is the independent variable
	Double_t Tn = p[0], theta1 = p[1], theta2 = p[2], phi1 = p[3], phi2 = p[4];
	Double_t particle = p[5];
	Double_t phi { phi1 - phi2};//sign criterium for Delta phi
	//assigning masses
	Double_t m1 { Mn}; Double_t m2 { (particle == true) ? Mn : Mh}; Double_t m3 { (particle == true) ? Mh : Mn};
	Double_t A {m2 + m3};
	Double_t cos_theta12 {TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Cos(phi)+TMath::Cos(theta1)*TMath::Cos(theta2)};
	Double_t B {2.*cos_theta12*TMath::Sqrt(m1*m2*T1) - 2.*TMath::Cos(theta2)*TMath::Sqrt(Mn*m2*Tn)};
	Double_t C {(m1 + m3)*T1 - m3*Q +Tn*(Mn - m3) -2.*TMath::Cos(theta1)*TMath::Sqrt(Mn*m1*Tn*T1)};
	Double_t D {B*B-4.*A*C};
	if(D>=0.)
	{
		Double_t Tneg {(-B-TMath::Sqrt(D))/(2.*A)};
		if(Tneg>=0.) return Tneg*Tneg;
		else return 0.;
	}
	else return 0.;
}

Double_t T3_calculator(Double_t* x, Double_t* p)
{
	//Function to pass as TF1 to compute S integral over proton energy vs Tn1
	Double_t T1 = x[0];//T1 is the independent variable
	Double_t Tn = p[0], theta1 = p[1], theta2 = p[2], phi1 = p[3], phi2 = p[4];
	Double_t particle = p[5];
	Double_t phi { phi1 - phi2};//sign criterium for Delta phi
	//assigning masses
	Double_t m1 { Mn}; Double_t m2 { (particle == true) ? Mn : Mh}; Double_t m3 { (particle == true) ? Mh : Mn};
	Double_t A {m2 + m3};
	Double_t cos_theta12 {TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Cos(phi)+TMath::Cos(theta1)*TMath::Cos(theta2)};
	Double_t B {2.*cos_theta12*TMath::Sqrt(m1*m2*T1) - 2.*TMath::Cos(theta2)*TMath::Sqrt(Mn*m2*Tn)};
	Double_t C {(m1 + m3)*T1 - m3*Q +Tn*(Mn - m3) -2.*TMath::Cos(theta1)*TMath::Sqrt(Mn*m1*Tn*T1)};
	Double_t D {B*B-4.*A*C};
	if(D>=0.)
	{
		Double_t Tpos {(-B+TMath::Sqrt(D))/(2.*A)};//only includes positive branch by now
		if(Tpos>=0.)
		{
			Double_t T2 { Tpos*Tpos};
			return Tn + Q - T1 - T2;
		} 
		else return 0.;
	}
	else return 0.;
}

#define NAN_VALUE "999"//type of NaN value for kinematical functions
two_doubles T2_null_points(Double_t Tn, Double_t theta1, std::string particle)//just depends on these two variables
{
	//assigning masses
	Double_t m1 { Mn}; Double_t m2 { (particle == "n") ? Mn : Mh}; Double_t m3 { (particle == "n") ? Mh : Mn};
	Double_t D { m1 + m3};
	Double_t E { 2*TMath::Cos(theta1)*TMath::Sqrt(Mn*m1*Tn)};
	//implicit change of sign in E because eq is: Dx**2 -Ex + F!
	Double_t F { (Mn-m3)*Tn -m3*Q};
	Double_t Delta {E*E-4.*D*F};
	if(Delta > 0.)
	{
		Double_t pos { (E+TMath::Sqrt(Delta))/(2.*D)};
		Double_t neg { (E-TMath::Sqrt(Delta))/(2.*D)};
		//return NaN values when conditions are not fullfilled!
		two_doubles out {(pos <= 0.) ? std::nan(NAN_VALUE) : pos*pos, (neg <= 0.) ? std::nan(NAN_VALUE) : neg*neg};
		return out;
	}
	two_doubles out {std::nan(NAN_VALUE), std::nan(NAN_VALUE)};
	return out;
}

two_doubles T2p_boundaries(Double_t Tn,
						   Double_t theta1, Double_t theta2,
						   Double_t phi1, Double_t phi2,
						   std::string particle)
{
	Double_t phi { phi1 - phi2};//criterium of signs for Delta phi
	//assigning masses
	Double_t m1 { Mn}; Double_t m2 { (particle == "n") ? Mn : Mh}; Double_t m3 { (particle == "n") ? Mh : Mn};
	Double_t cos_theta12 {TMath::Sin(theta1)*TMath::Sin(theta2)*TMath::Cos(phi)+TMath::Cos(theta1)*TMath::Cos(theta2)};
	Double_t B1 { 2.*cos_theta12*TMath::Sqrt(m1*m2)};
	Double_t B2 { 2.*TMath::Cos(theta2)*TMath::Sqrt(Mn*m2*Tn)};
	Double_t A { m2 + m3};
	Double_t C1 { m1 + m3};
	Double_t C2 { 2.*TMath::Cos(theta1)*TMath::Sqrt(Mn*m1*Tn)};
	Double_t C3 { (Mn - m3)*Tn - m3*Q};
	Double_t I { B1*B1 - 4.* A*C1};
	Double_t II { 4*A*C2 - 2.*B1*B2};
	Double_t III { B2*B2-4.*A*C3};
	Double_t Delta { II*II - 4.*I*III};
	if(Delta > 0.)
	{
		Double_t pos { (-II+TMath::Sqrt(Delta))/(2.*I)};
		Double_t neg { (-II-TMath::Sqrt(Delta))/(2.*I)};
		//also return NaN values when this goes wrong
		two_doubles out {(pos <= 0.) ? std::nan(NAN_VALUE) : pos*pos, (neg <= 0.) ? std::nan(NAN_VALUE) : neg*neg};
		return out;
	}
	two_doubles out {std::nan(NAN_VALUE), std::nan(NAN_VALUE)};
	return out;
}

Double_t S_explicit_integrand(Double_t* x, Double_t* p)//ready to be passed as TF1
{
	Double_t T1 = x[0];//T1 is the independent variable
	Double_t Tn = p[0], theta1 = p[1], theta2 = p[2], phi1 = p[3], phi2 = p[4],
		branch = p[5], particle = p[6], line = p[7];
	//Assigning masses
	Double_t m1 { Mn}; Double_t m2 { (particle == true) ? Mn : Mh}; Double_t m3 { (particle == true) ? Mh : Mn};
	//Momenta
	Double_t K[6] {Tn, theta1, theta2, phi1, phi2, particle};//vector to call T2n,p functions
	Double_t p1 { TMath::Sqrt(2.*m1*T1)};
	Double_t p2;
	//Whether is positive or negative branch
	if (branch) p2 = TMath::Sqrt(2.*m2*T2p_calculator(x, K));//if true = 1
	else p2 = TMath::Sqrt(2.*m2*T2n_calculator(x, K));//if false = 0
	//if T2n,p is null, return 0.!
	if (p2 <= 0.) return 0.;

	Double_t pn { TMath::Sqrt(2.*Mn*Tn)};
	
	//T3Vectors
	//following results from https://sci-hub.hkvisa.net/10.1103/PhysRevC.79.014609
	ROOT::Math::XYZVector versor_pn (0., 0., 1.);
	ROOT::Math::XYZVector versor_p1 (TMath::Sin(theta1)*TMath::Cos(phi1),
									 TMath::Sin(theta1)*TMath::Sin(phi1),
									 TMath::Cos(theta1));
	ROOT::Math::XYZVector versor_p2 (TMath::Sin(theta2)*TMath::Cos(phi2),
									 TMath::Sin(theta2)*TMath::Sin(phi2),
									 TMath::Cos(theta2));
	auto vector_pn { pn * versor_pn};
	auto vector_p1 { p1 * versor_p1};
	auto vector_p2 { p2 * versor_p2};

	Double_t numerator   { - m2 * (p1 * (m1 + m3) - m1 * (vector_pn - vector_p2).Dot(versor_p1))};
	Double_t denominator { m1 * (p2 * (m2 + m3) - m2 * (vector_pn - vector_p1).Dot(versor_p2))};

	//integrand
	Double_t integrand;
	//if is primary (T2) kinematic line
	if(line) integrand = TMath::Sqrt(1. + TMath::Power( (m1 * p2 * numerator) / (m2 * p1 *denominator), 2));//n_n
	//else integrand = TMath::Sqrt(1. + TMath::Power(1 + numerator / denominator, 2));//p_n
	else
	{
		Double_t p3dp3dp1 { p1 + (vector_p2 - vector_pn).Dot(versor_p1) + numerator / denominator * (p2 + (vector_p1 - vector_pn).Dot(versor_p2))};
		integrand = TMath::Sqrt(1 + TMath::Power((m1 * 1.) / (m3 * p1) * p3dp3dp1, 2));
	}

	return integrand;
}

//function which converts T1 to T2, keeping S
std::vector<Double_t> convert_T1_T2(std::vector<Double_t> vS, std::vector<Double_t> vT1, Double_t* K)
{
	//get maximum and minimum values in S, which determine +/- branch
	std::vector<Double_t>::iterator index_max { std::max_element( vT1.begin(), vT1.end())};
	std::vector<Double_t>::iterator index_min { std::min_element( vT1.begin(), vT1.end())};
	Double_t S_T1max { vS[std::distance(vT1.begin(), index_max)]};
	Double_t S_T1min { vS[std::distance(vT1.begin(), index_min)]};
	std::vector<Double_t> vT2;
	//std::cout<<"Converting T1 to T2"<<'\n';
	for(Int_t j = 0; j < vT1.size(); j++)
	{
		if(vS[j] <= S_T1max) vT2.push_back(T2n_calculator(&vT1[j], K));
		else if ((vS[j] > S_T1max) && (vS[j] < S_T1min)) vT2.push_back(T2p_calculator(&vT1[j], K));
		else vT2.push_back(T2n_calculator(&vT1[j], K));
	}
	if(vT2.size() != vT1.size())
	{
		std::cout<<"vT1 and vT2 do not share dimensions";
		std::abort();
	}
	return vT2;
}

three_vectors get_S_curve(std::string energy,
						Double_t theta1, Double_t theta2,
						Double_t phi1, Double_t phi2,
						std::string particle, std::string line)
{
	//get energy from string
	auto Tn { energy_from_string(energy)};
	
	//Kinematical points of interest along S curve
	auto [s1, s3] = T2_null_points(Tn, deg_to_rad(theta1), particle);
	auto [s4, s2] = T2p_boundaries(Tn, deg_to_rad(theta1), deg_to_rad(theta2),
								   deg_to_rad(phi1), deg_to_rad(phi2), particle);
	//Check if kinematics exist for this configuration
	if(std::isnan(s1) && std::isnan(s2))
	{
		std::cout<<"Kinematics does not exist for this configuration"<<std::endl;
		std::vector<Double_t> a, b, c;
		three_vectors out { a, b, c};
		return out;
	}
	//std::cout<<"s1: "<<s1<<" MeV and s3: "<<s3<<" MeV\n";
	//std::cout<<"s2: "<<s2<<" MeV and s4: "<<s4<<'\n';
	
	//Else define function
	auto* fS = new TF1("fS", S_explicit_integrand, 0., s2, 8);
	fS->SetParameters(Tn, deg_to_rad(theta1), deg_to_rad(theta2),
					  deg_to_rad(phi1), deg_to_rad(phi2),
					  0, //branch: 0 por negative and 1 for positive
					  (particle == "n") ? true : false,//true for neutron, false for proton
					  (line == "primary") ? true : false );//primary = true(1), secondary = false (0)
	//vector of kinematic configuration (function takes angles as degrees)
	Double_t K[6] {Tn,deg_to_rad(theta1), deg_to_rad(theta2),
				   deg_to_rad(phi1), deg_to_rad(phi2),
				   static_cast<Double_t>((particle == "n") ? true : false)};
	
	//Define epsilon values to avoid integration errors
	Double_t roundoff { 1.E-4}; Double_t zero { 1.E-3};
	//Step of integration
	Double_t step { 1.E-2};

	//1 to 2
	std::vector<Double_t> vT12;
	std::vector<Double_t> vS12;
	if (std::isnan(s1))
	{
		vT12 =  arange( (std::isnan(s4)) ? zero : s4 + roundoff, s2, step);
		for(auto& T1 : vT12)
		{
			vS12.push_back(fS->Integral((std::isnan(s4)) ? zero : s4 + roundoff, T1));
		}
	}
	else
	{
		vT12 = arange(s1, s2, step);
		for(auto& T1 : vT12)
		{
			vS12.push_back(fS->Integral(s1, T1));
		}
	}

	//2 to zero or 4
	fS->SetParameter(5, 1);//switch to positive branch
	std::vector<Double_t> vT24(arange( (std::isnan(s4)) ? zero : s4+roundoff, s2, step));
	std::reverse( vT24.begin(),  vT24.end());
	std::vector<Double_t> vS24;
	for (auto& T1 : vT24)
	{
		vS24.push_back(vS12.back() + fS->Integral(T1, s2-roundoff)); 
	}

	//zero or 4 to 3
	fS->SetParameter(5, 0);
	std::vector<Double_t> vT43;
	std::vector<Double_t> vS43;
	if(std::isnan(s3))
	{
		;//just do nothing in this case
	}
	else
	{
		vT43 = arange((std::isnan(s4)) ? zero + roundoff : s4 + 2*roundoff, s3, step);
		for (auto& T1 : vT43)
		{
			vS43.push_back(vS24.back() + fS->Integral((std::isnan(s4)) ? zero : s4 + roundoff, T1));
		}
	}

	//and now merge all vectors
	std::vector<Double_t> vT1234;
	std::vector<Double_t> vS1234;
	vT1234.reserve(vT12.size() + vT24.size() + vT43.size()); vT1234.insert(vT1234.end(), vT12.begin(), vT12.end());//continues bellow
	vT1234.insert(vT1234.end(), vT24.begin(), vT24.end()); vT1234.insert(vT1234.end(), vT43.begin(), vT43.end());
	vS1234.reserve(vS12.size() + vS24.size() + vT43.size()); vS1234.insert(vS1234.end(), vS12.begin(), vS12.end());
	vS1234.insert(vS1234.end(), vS24.begin(), vS24.end()); vS1234.insert(vS1234.end(), vS43.begin(), vS43.end());

	//and finally compute vT2 from the vS1234 and vT1234
	auto vT_second = convert_T1_T2(vS1234, vT1234, K);
	
	//finally, write to file
	std::ofstream fout(TString::Format("%s/%s/%s/%s/%s/S_t1_%.2f_p1_%.2f_t2_%.2f_p2_%.2f.dat",
									   directory.c_str(), subdirS.c_str(),
									   energy.c_str(), particle.c_str(),
									   line.c_str(),
									   theta1, phi1, theta2, phi2).Data());
	//now X is S (MeV) and Y is T1 (MeV)
	fout << std::fixed << std::setprecision(7) << std::endl;//make sure we write with enough precision to avoid duplicates
	for(Int_t i = 0; i < vS1234.size(); i++)
	{
		fout << vS1234[i] << " " << vT1234[i] << " " << vT_second[i] <<'\n';
	}
	fout.close();
	//delete
	delete fS;	
	//and return the vectors just in case
	three_vectors out { vS1234, vT1234, vT_second};
	return out;
	
}

//AND FINALLY UTILITIES TO READ FILES!
//Angles index file from theo_xs_getter.C
five_vectors read_angles_file(const char* file)
{
	std::vector<Double_t> index, theta1, phi1, theta2, phi2;
	std::fstream infile;
	infile.open(file);
	Double_t ai,at1,ap1,at2,ap2; //auxiliars for storing values on while loop
	while (infile >> ai >> at1 >> ap1 >> at2 >> ap2)
	{
		index.push_back(ai);
		theta1.push_back(at1); phi1.push_back(ap1);
		theta2.push_back(at2); phi2.push_back(ap2);
	}
	infile.close();
	//class variable construction
	five_vectors out{index, theta1, phi1, theta2, phi2};
	return out;
}

//read file from get_S_curve()
three_vectors get_S_file(std::string energy,
					   Double_t theta1, Double_t theta2,
					   Double_t phi1, Double_t phi2,
					   std::string particle,
					   std::string line,
					   Bool_t force_rewrite=false)
{
	//get energy from string
	auto Tn { energy_from_string(energy)};
	
	three_vectors out;
	std::string name { TString::Format("%s/%s/%s/%s/%s/S_t1_%.2f_p1_%.2f_t2_%.2f_p2_%.2f.dat",
									   directory.c_str(), subdirS.c_str(),
									   energy.c_str(), particle.c_str(),
									   line.c_str(),
									   theta1, phi1, theta2, phi2).Data()};

	//check if file exists or if we impose a rewrite
	if(file_exists(name) && force_rewrite == false)
	{
		std::ifstream infile (name);
		std::vector<Double_t> S, T1, T2;
		Double_t auxS, auxT1, auxT2;
		while(infile >> auxS >> auxT1 >> auxT2)
		{
			S.push_back(auxS); T1.push_back(auxT1); T2.push_back(auxT2);
		}
		infile.close();
		//check for duplicates in S
		const auto duplicate = std::adjacent_find(S.begin(), S.end());
		if (duplicate != S.end())
		{
			std::cout << BOLDRED <<"Duplicate " << *duplicate <<" at t1: "
					  <<theta1<<" p1:"<<phi1<<" t2: "<<theta2<<" and p2: "<<phi2<< RESET <<"\n";
			std::abort();
		}
		out = { S, T1, T2};
		return out;
	}
	else
	{
		out = get_S_curve(energy, theta1, theta2, phi1, phi2, particle, line);
		return out;
	}
}

//OTHER UTILITIES
TH1D* histogram_from_file(TGraphErrors* g)
{
	Int_t npts { g->GetN()};
	std::cout<<"N theoretical values: "<< npts <<std::endl;
	//Now create arrays with values
	Double_t* xep { new Double_t[npts]};
	Double_t* yxs { new Double_t[npts]};
	xep = g->GetX();
	yxs = g->GetY();
	Int_t nbins {npts};
	//maximum and minimum values
	Double_t xmin {xep[0]};
	Double_t xmax {xep[npts-1]};
	//Histogram
	auto* hr = new TH1D("", "Hist for xs", nbins, xmin, xmax);
	for(Int_t i = 0; i<npts; i++)
	{
		hr->SetBinContent(i+1, yxs[i]);
	}
	return hr;
}

//find bins for energy within uncertainty while doing the projection
two_ints get_first_last_bins(TH2D* h, Double_t T, Double_t DeltaT)
{
	two_ints out {h->GetYaxis()->FindFixBin(T-DeltaT), h->GetYaxis()->FindFixBin(T+DeltaT)};
	return out;
}//with FindFixBin we do not include under/overflows

//the same but for TH3D histograms
two_ints get_first_last_bins_3d(TH3D* h, Double_t T, Double_t DeltaT)
{
	two_ints out {h->GetYaxis()->FindFixBin(T-DeltaT), h->GetYaxis()->FindFixBin(T+DeltaT)};
	return out;
}

//RESOLUTION FUNCTIONS
//by interpolation of graph (old version)
std::vector<Double_t> resolution_x { 0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20.};
std::vector<Double_t> resolution_y { 0.025, 0.0385, 0.048, 0.056, 0.064, 0.0720, 0.0780, 0.0825, 0.086, 0.092, 0.097 };//given in Delta E /E for L=2m
auto* resolution_spline = new TSpline3("resolution_spline", &(resolution_x[0]), &(resolution_y[0]), resolution_x.size(), "b2,e2", 0., 0.);
Double_t resolution_c_functor(Double_t* x, Double_t* p){ return resolution_spline->Eval(x[0]);}
auto* resolution_func   = new TF1("resolution_func", resolution_c_functor, 0., 20., 1);
//Real resolution from MONSTER applying formula and data from paper
Double_t tof_resolution_1(Double_t* x, Double_t* p)
{
	Double_t T = x[0];
	Double_t L = 1.5; //m
	Double_t Deltat { 1.5e-9}; //s, from MONSTER paper
	Double_t DeltaL { 0.025}; //m, from MONSTER paper
	Double_t v { TMath::Sqrt(2.*T / Mn) * 3.e8}; //they are always neutrons, multiplied by c
	Double_t time { L / v};
	Double_t res { 2.*TMath::Sqrt(TMath::Power(DeltaL / L, 2) + TMath::Power(Deltat / time, 2))};
	return res;
}

Double_t tof_resolution_2(Double_t* x, Double_t* p)
{
	Double_t T = x[0];
	Double_t L = 2.5; //m
	Double_t Deltat { 1.5e-9}; //s, from MONSTER paper
	Double_t DeltaL { 0.025}; //m, from MONSTER paper
	Double_t v { TMath::Sqrt(2.*T / Mn) * 3.e8}; //they are always neutrons, multiplied by c
	Double_t time { L / v};
	Double_t res { 2.*TMath::Sqrt(TMath::Power(DeltaL / L, 2) + TMath::Power(Deltat / time, 2))};
	return res;
}
auto* monster_resolution1 = new TF1("monster_resolution1", tof_resolution_1, 0., 20., 1);
auto* monster_resolution2 = new TF1("monster_resolution2", tof_resolution_2, 0., 20., 1);

//resolution on ACTAR TPC Si detectors
//1. energy straggling in gas
const Double_t straggling_gas { 0.0141}; //MeV
//2. silicon resolution
auto* silicon_sigma = new TF1("silicon_sigma", "0.0213 / 2.35 * TMath::Sqrt(x)", 0., 20.);
Double_t quadratic_sum_actar(Double_t* x, Double_t* p){ return straggling_gas + silicon_sigma->Eval(x[0]);}
auto* actar_sigma = new TF1("actar_sigma", quadratic_sum_actar, 0., 20., 1);

//Computes convolution displacing bin centers by gaussian of the desired sigma.
//Repeats the process iterations times. New values for histogram (bins, init, end) may be given, but better keep same values as projection
TH1D* histogram_convolution(TH1D* projection, //Double_t sigma, //sigma is provided by us now
							Double_t init, Double_t end,
							Int_t bins = 100,
							Int_t iterations = 100)
{
	auto* rand = new TRandom3();
	auto* proj_sum = new TH1D("", "Convolution", bins, init, end);
	auto* proj_new = new TH1D("", "Iteration histogram", bins, init, end);
	proj_sum->Sumw2();
	for(Int_t it = 0; it < iterations; it++)
	{
		//Iterate over projection bins
		for(Int_t i = 0; i < projection->GetXaxis()->GetNbins(); i++)
		{
			Double_t content { projection->GetBinContent(i)};
			if(!(content>0.)) continue; //continue if bin is empty
			Double_t error { projection->GetBinError(i)};
			Double_t center{ projection->GetBinCenter(i)};
			Double_t sigma { resolution_func->Eval(center) * center};
			Double_t new_center { projection->GetBinCenter(i) +
				rand->Gaus(0., sigma)};
			Int_t new_bin { proj_new->FindBin(new_center)};
			proj_new->AddBinContent(new_bin, content);
			//and now add error in an accumulative way (squared)
			if(std::abs(center - new_center) <= sigma) proj_new->SetBinError(new_bin, error);
			else proj_new->SetBinError(new_bin, error);//by default, all displaced bins get the same error... maybe a good idea
			//changes error as function of distance to original bin
		}
		//Append to sum histogram
		proj_sum->Add(proj_new);
		//Prepare to fill again
		proj_new->Reset();
	}
	//proj_sum->Scale(1. / proj_sum->Integral(), "");
	//Manually threatment of errors
	//Erros should be divided by \sqrt{iterations} as they are a standard deviation
	Double_t integral { proj_sum->Integral()};
	for(Int_t j = 0; j < proj_sum->GetNbinsX(); j++) {proj_sum->SetBinContent(j, proj_sum->GetBinContent(j) / integral); proj_sum->SetBinError(j, proj_sum->GetBinError(j) / std::sqrt(iterations));}
	proj_sum->ResetStats();
	//deletes
	delete rand;
	delete proj_new;
	return proj_sum;
}

//CONVOLUTIONS OVER HISTOGRAMS (deprecated in our code)
TH2D* histogram_convolution_2d(TH2D* projection, //Double_t sigma, //sigma is provided by us now
							Double_t init, Double_t end,
							Int_t bins = 100,
							Int_t iterations = 100)
{
	auto* rand = new TRandom3();
	auto* proj_sum = new TH2D("", "Convolution", bins, init, end, bins, init, 25.);
	auto* proj_new = new TH2D("", "Iteration histogram", bins, init, end, bins, init, 25.);
	auto* proj_aux = new TH2D("", "Summing histogram", bins, init, end, bins, init, 25.);
	proj_sum->Sumw2();
	proj_new->Sumw2();
	for(Int_t it = 0; it < iterations; it++)
	{
		//Iterate over projection bins
		for(Int_t j = 0; j < projection->GetNbinsY(); j++)
		{
			for(Int_t i = 0; i < projection->GetNbinsX(); i++)
			{
				//std::cout<<"i,j: "<<i<<" "<<j<<'\n';
				Double_t content { projection->GetBinContent(i, j)};
				if(content <= 0.) continue;
				//std::cout<<content<<'\n';
				Double_t x_center { projection->GetXaxis()->GetBinCenter(i)};
				Double_t y_center { projection->GetYaxis()->GetBinCenter(j)};
				Double_t sigma { resolution_func->Eval(x_center) * x_center};
				Double_t new_x_center {x_center + rand->Gaus(0., sigma)};
				Int_t new_x_bin { proj_aux->GetXaxis()->FindBin(new_x_center)};
				//std::cout<<"new_x_bin: "<<new_x_bin<<'\n';
				proj_aux->SetBinContent(new_x_bin, j, content);
				proj_new->Add(proj_aux);
				proj_aux->Reset();
				//std::cout<<"Nentries: "<<proj_new->GetEntries()<<'\n';
			}
		}
		proj_sum->Add(proj_new);
		proj_new->Reset();
	}
	//and finally divide by iterations to return get mean in BinContent (I Think we can do this with root built-in kIsAverage static function)
	for(Int_t j = 0; j < proj_sum->GetNbinsY(); j++)
	{
		for (Int_t i = 0; i < proj_sum->GetNbinsX(); i++)
		{
			Double_t content { proj_sum->GetBinContent(i, j)};
			Double_t error   { proj_sum->GetBinError(i, j)};
			proj_sum->SetBinContent(i, j, content / iterations);
			proj_sum->SetBinError(i, j, content / TMath::Sqrt(iterations));
		}
	}
	//proj_sum->ResetStats();
	delete rand;
	delete proj_aux;
	delete proj_new;
	return proj_sum;
}

//TH3D Convolution
TH3D* histogram_convolution_3d(TH3D* hist, //Double_t sigma, //sigma is provided by us now
							Double_t init, Double_t end,
							Int_t bins = 100,
							Int_t iterations = 100)
{
	auto* rand = new TRandom3();
	auto* out  = (TH3D*)hist->Clone();
	out->Reset();
	auto* hist_iteration = (TH3D*) out->Clone();
	//auto* hist_aux       = (TH3D*) out->Clone();
	//auto* hist_T1        = (TH3D*) out->Clone();
	//auto* hist_T2        = (TH3D*) out->Clone();
	out->Sumw2();
	hist_iteration->Sumw2();
	for(Int_t it = 0; it < iterations; it++)
	{
		//std::cout<<"Iteration: "<<it<<'\n';
		for(Int_t k = 0; k < hist->GetNbinsZ(); k++)
		{
			for(Int_t j = 0; j < hist->GetNbinsY(); j++)
			{
				for(Int_t i = 0; i < hist->GetNbinsX(); i++)
				{
					//std::cout<<"i,j: "<<i<<" "<<j<<'\n';
					Double_t content { hist->GetBinContent(i, j, k)};
					if(content <= 0.) continue;
					//std::cout<<content<<'\n';
					Double_t x_center { hist->GetXaxis()->GetBinCenter(i)};
					Double_t y_center { hist->GetYaxis()->GetBinCenter(j)};
					Double_t z_center { hist->GetZaxis()->GetBinCenter(k)};
					Double_t x_sigma { (it == 0) ? 0. : resolution_func->Eval(x_center) * x_center};
					Double_t y_sigma { (it == 0) ? 0. : resolution_func->Eval(y_center) * y_center};
					Double_t x_new_center {x_center + rand->Gaus(0., x_sigma)};
					Double_t y_new_center {y_center + rand->Gaus(0., y_sigma)};
					//Int_t x_new_bin { hist_aux->GetXaxis()->FindBin(x_new_center)};
					//Int_t y_new_bin { hist_aux->GetYaxis()->FindBin(y_new_center)};
					//std::cout<<"new_x_bin: "<<new_x_bin<<'\n';
					//hist_aux->SetBinContent(x_new_bin, y_new_bin, k, content);
					//hist_iteration->Add(hist_aux);
					//hist_aux->Reset();
					//std::cout<<"Nentries: "<<proj_new->GetEntries()<<'\n';
					hist_iteration->Fill(x_new_center, y_new_center, z_center);
				}
			}
		}
		out->Add(hist_iteration);
		hist_iteration->Reset();
	}
	//and finally divide by iterations to return get mean in BinContent (I Think we can do this with root built-inkIsAverage static function)
	/*
	for(Int_t k = 0; k < out->GetNbinsZ(); k++)
	{
		for(Int_t j = 0; j < out->GetNbinsY(); j++)
		{
			for (Int_t i = 0; i < out->GetNbinsX(); i++)
			{
				Double_t content { out->GetBinContent(i, j, k)};
				Double_t error   { out->GetBinError(i, j, k)};
				out->SetBinContent(i, j, content / iterations);
				out->SetBinError(i, j, content / TMath::Sqrt(iterations));
			}
		}
	}
	*/
	//proj_sum->ResetStats();
	delete rand;
	//delete hist_aux;
	delete hist_iteration;
	return out;
}

//function to compute S distance when clustering (deprecated)
Double_t S_distance(Double_t x0, Double_t y0, Double_t x1, Double_t y1)
{
	Double_t val { TMath::Sqrt(TMath::Power(x0 - x1, 2) + TMath::Power(y0- y1, 2))};
	return val;
}
//FUNCTIONS TO OBTAIN RESULTS
//get all data from TGraphAsymmErrors
four_vectors get_graph_points(TGraphAsymmErrors* g)
{
	std::vector<Double_t> X, Y, EYLow, EYHigh;
	Int_t npoints { g->GetN()};
	Double_t* x = g->GetX();
	Double_t* y = g->GetY();
	Double_t* eylow = g->GetEYlow();
	Double_t* eyhigh = g->GetEYhigh();
	X = std::vector<Double_t>(x, x+npoints);
	Y = std::vector<Double_t>(y, y+npoints);
	EYLow = std::vector<Double_t>(eylow, eylow+npoints);
	EYHigh = std::vector<Double_t>(eyhigh, eyhigh+npoints);

	four_vectors out { X, Y, EYLow, EYHigh};
	return out;
}
//include statistical error in yield
void include_statistical_error(TGraphAsymmErrors* g, Bool_t rewrite_error=false)
{
	for(Int_t i = 0; i < g->GetN(); i++)
	{
		Double_t Y { g->GetPointY(i)};
		Double_t EYLow { g->GetErrorYlow(i)};
		Double_t EYHigh{ g->GetErrorYhigh(i)};
		Double_t sqrt_error { TMath::Sqrt(Y)};
		if(!rewrite_error)
		{
			g->SetPointEYlow(i, EYLow + sqrt_error);
			g->SetPointEYhigh(i, EYHigh + sqrt_error);
		}
		else//if true replace errors statistical ones
		{
			g->SetPointEYlow(i, sqrt_error);
			g->SetPointEYhigh(i, sqrt_error);
		}
	}
}

//get quadrature ratio
two_vectors ratio_quadrature(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2)
{
	std::vector<Double_t> ratio, sigma, axis;
	//g2 is here the lower ann value; check depending on anns vector sorting!
	auto [x1, y1, eylow1, eyhigh1] = get_graph_points(g1);
	auto [x2, y2, eylow2, eyhigh2] = get_graph_points(g2);
	Int_t npoints1 { g1->GetN()}; Int_t npoints2 { g2->GetN()};
	for(Int_t i = 0; i < npoints1; i++)
	{
		Double_t a { x1[i]};
		//check if x values agree
		if(!(std::find_if(x2.begin(), x2.end(), [a](Double_t b){return equal_doubles(a, b);}) != x2.end())) continue;
		Double_t numerator { std::abs(y2[i] - y1[i])};
		Double_t denominator {TMath::Sqrt(TMath::Power(eylow1[i], 2 ) + TMath::Power(eylow2[i], 2 ))};
		Double_t division { numerator / denominator};
		if(!std::isfinite(division)) continue;//std::isfinite checks if no NaN or Infs are in value
		//otherwise push to ratio and axis vectors
		axis.push_back(x1[i]); ratio.push_back(division);
	}
	two_vectors out { axis, ratio};
	return out;		
}

//calculate ratio without sqrt in denominator, returning itself and the good X axis
two_vectors ratio_calculator(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2)
{
	std::vector<Double_t> ratio, axis;
	//g2 is here the lower ann value; check depending on anns vector sorting!
	auto [x1, y1, eylow1, eyhigh1] = get_graph_points(g2);
	auto [x2, y2, eylow2, eyhigh2] = get_graph_points(g1);
	Int_t npoints1 { g1->GetN()}; Int_t npoints2 { g2->GetN()};
	for(Int_t i = 0; i < npoints1; i++)
	{
		Double_t a { x1[i]};
		//check if x values agree
		if(!(std::find_if(x2.begin(), x2.end(), [a](Double_t b){return equal_doubles(a, b);}) != x2.end())) continue;
		Double_t numerator { y2[i] - y1[i]};
		Double_t denominator { eylow2[i] + eyhigh1[i]};
		Double_t division { numerator / denominator};
		if(!std::isfinite(division)) continue;//std::isfinite checks if no NaN or Infs are in value
		//otherwise push to ratio and axis vectors
		axis.push_back(x1[i]); ratio.push_back(division);
	}
	two_vectors out { axis, ratio};
	return out;	
}

//perform Welch t-test if it were possible
TH1D* welch_test(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2)
{
	std::vector<Double_t> ratio, axis;
	//now up and low errors are the sameeee
	auto [x1, y1, eylow1, eyhigh1] = get_graph_points(g1);
	auto [x2, y2, eylow2, eyhigh2] = get_graph_points(g2);
	Int_t npoints1 { g1->GetN()}; Int_t npoints2 { g2->GetN()};
	auto* h = new TH1D("h", "Welch test", 120, 0., 17.);
	for(Int_t i = 0; i < npoints1; i++)
	{
		Double_t a { x1[i]};
		//check if x values agree
		if(!(std::find_if(x2.begin(), x2.end(), [a](Double_t b){return equal_doubles(a, b);}) != x2.end())) continue;

		Int_t iterations { static_cast<Int_t>(1.E4)};
		auto* h1 = new TH1D("h1", "g1 points", 10000, 0., 5000.);
		auto* h2 = new TH1D("h2", "g2 points", 10000, 0., 5000.);
		Double_t point1 { y1[i]};
		Double_t sigma1 { eylow1[i]};
		Double_t point2 { y2[i]};
		Double_t sigma2 { eylow2[i]};

		for(Int_t j = 0; j < iterations; j++)
		{
			h1->Fill(point1 + gRandom->Gaus(0., sigma1));
			h2->Fill(point2 + gRandom->Gaus(0., sigma2));
		}
		Double_t X1 { h1->GetMean()}; Double_t X2 { h2->GetMean()};
		Double_t S1 { h1->GetStdDev() / TMath::Sqrt(h1->GetEntries())};
		Double_t S2 { h2->GetStdDev() / TMath::Sqrt(h2->GetEntries())};
		delete h1;
		delete h2;
		Double_t DeltaX { std::fabs(X1-X2)};
		Double_t sigmaDelta { TMath::Sqrt(S1*S1 + S2*S2)};
		//std::cout<<"sigmaDelta: "<<sigmaDelta<<'\n';
		Double_t t { DeltaX / sigmaDelta};
		if(!(std::isfinite(t))) continue;
		std::cout<<"t: "<<t<<'\n';
		Double_t nu { TMath::Power(sigmaDelta, 4) / ((TMath::Power(S1, 4) / (iterations-1)) + (TMath::Power(S2, 4) / (iterations - 1)))};
		if(!(std::isfinite(nu))) continue;
		Double_t t0 { TMath::StudentQuantile(0.0025, (Int_t)nu, kFALSE)};;
		//std::cout<<"Delta: "<<DeltaX<<" sigmaDelta: "<<sigmaDelta<<" t: "<<t<<'\n';
		if(t>=t0)h->Fill(x1[i]);

	}
	return h;
}

//confidence interval plot
TGraphErrors* difference_plotter(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2)
{
	std::vector<Double_t> mean, sigma, axis;
	//g2 is here the lower ann value; check depending on anns vector sorting!
	auto [x1, y1, eylow1, eyhigh1] = get_graph_points(g1);
	auto [x2, y2, eylow2, eyhigh2] = get_graph_points(g2);
	Int_t npoints1 { g1->GetN()}; Int_t npoints2 { g2->GetN()};
	for(Int_t i = 0; i < npoints1; i++)
	{
		Double_t a { x1[i]};
		//check if x values agree
		if(!(std::find_if(x2.begin(), x2.end(), [a](Double_t b){return equal_doubles(a, b);}) != x2.end())) continue;
		Double_t diff { std::abs(y2[i] - y1[i])};
		Double_t sigma_diff { TMath::Sqrt(TMath::Power(eylow1[i], 2) + TMath::Power(eylow2[i], 2))};
		if(!std::isfinite(sigma_diff)) continue;//std::isfinite checks if no NaN or Infs are in value
		//otherwise push to ratio and axis vectors
		axis.push_back(x1[i]);
		mean.push_back(diff);
		sigma.push_back(1.0 * sigma_diff);
	}
	auto* graph = new TGraphErrors(axis.size(), &(axis[0]), &(mean[0]), NULL, &(sigma[0]));
	return graph;	
}

#endif
