#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "json.hpp"

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraphPolar.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TLatex.h"
#include "TFile.h"
#include "TStyle.h"
#include "TImage.h"

using namespace std;
using json = nlohmann::json;

//Graph object used in fitting and plotting
TGraphPolar* phi_vs_r_data;

//Number of fit parameters: a, b, c
const int NPARS = 3;

//Colors to be used in plotting
const int colors[10] = {4, 3, 2, 1, 6, 5, 7, 46, 41, 28}; //blue, green, red, black, magenta, yellow, cyan, gold, copper, bronze

//Fit function, r versus phi
double FitFunc(const double *x, const double *pars)
{
  double phi = x[0];
  double cn = pars[0];
  double an = pars[1];
  double bn = pars[2];
  double nsymm = pars[3];
  double phi0 = pars[4];
  return cn + an*cos(nsymm*(phi-phi0)) + bn*cos(2*nsymm*(phi-phi0));
}

//Chi-square function to be minimized
double FuncToMin(const double *pars)
{
  double chi2 = 0;
  double r_uncertainty = 0.01;
  for(int ip=0;ip<phi_vs_r_data->GetN();ip++)
  {
    double phi = phi_vs_r_data->GetX()[ip];
    double rdata = phi_vs_r_data->GetY()[ip];
    double rcalc = FitFunc(&phi,pars);
    double chi = (rdata-rcalc)/r_uncertainty;
    chi2 += chi*chi;
  }
  return chi2;
}

//Main routine
int main(int argc, char *argv[])
{
  //Command line arguments
  if(argc<3) { cerr << "Usage: " << argv[0] << " <symmetry order, 3 or 4, anything else means circle> <input/output file basename> <phi0 in deg>" << endl; return 1; }
  int nsymm = atoi(argv[1]);
  if(nsymm!=3 && nsymm!=4) { cerr << "Allowed symmetry orders: 3 or 4. Setting it to -1 (meaning: circular plate)" << endl; nsymm = -1; }
  string outfilebase = argv[2];
  string infilename = outfilebase + ".json";
  double phi0 = 0;
  if(argc>3) phi0 = atoi(argv[3])*M_PI/180.;
  
  
  // read a JSON file
  ifstream infile(infilename);
  json jsondata;
  infile >> jsondata;
  //cout << "JSON size: " << jsondata.size() << endl;
    
  int Ndata = jsondata["datasetColl"].size();
  cout << "number of datasets: " << Ndata << endl;
  vector<vector<pair<double,double>>> phi_r_pair_vecs;
  for(int idataset=0;idataset<Ndata;idataset++)
  {
    //Store data in vecotr of pairs
    int Npoints = jsondata["datasetColl"][idataset]["data"].size();
    cout << "  dataset size: " << Npoints << endl;
    vector<pair<double,double>> phi_r_pair_vec;
    vector<double> rvec;
    vector<double> phivec;
    for(int ipoint=0;ipoint<Npoints;ipoint++)
    {
      phi_r_pair_vec.emplace_back(
        jsondata["datasetColl"][idataset]["data"][ipoint]["value"][1],  //phi
        jsondata["datasetColl"][idataset]["data"][ipoint]["value"][0]); //r
    }
    //Sort and save
    sort(phi_r_pair_vec.begin(), phi_r_pair_vec.end());
    phi_r_pair_vecs.push_back(phi_r_pair_vec);
    
    /*
    //Write to text flle
    string outfilename = outfilebase + "_curve" + to_string(idataset+1) + ".txt";
    ofstream outfile(outfilename);
    for(auto& it : phi_r_pair_vec)
    {
      outfile << it.first << "\t" << it.second << endl;
    }
    outfile.close();
    */
  }
  
  //Graph objects for drawing
  vector<TGraphPolar*> phi_vs_r_data_plot;
  vector<TGraphPolar*> phi_vs_r_fit;
  const int Ncurves = phi_r_pair_vecs.size();
  //TGraphPolar* phi_vs_r_fit[Ncurves];
  
  //Parameters and errors to be saved
  vector<double> anval(Ncurves,0);
  vector<double> bnval(Ncurves,0);
  vector<double> cnval(Ncurves,0);
  vector<double> anerr(Ncurves,0);
  vector<double> bnerr(Ncurves,0);
  vector<double> cnerr(Ncurves,0);
  
  int icurve = -1;
  for(auto& phi_r_pair_vec : phi_r_pair_vecs)
  {
    icurve++;
    //Create Graph object from data in input file
    vector<double> phivector, rvector;
    for(auto& phi_r_pair : phi_r_pair_vec)
    {
      phivector.push_back(phi_r_pair.first*M_PI/180.); //convert to radians
      rvector.push_back(abs(phi_r_pair.second)); //take absolute value as sometimes negatives are read in the json file
    }
    phi_vs_r_data = new TGraphPolar(phivector.size(),&phivector[0],&rvector[0]);
    
    //Save data graph to array
    phi_vs_r_data_plot.push_back((TGraphPolar*)phi_vs_r_data->Clone());
    
    //Minuit2 minimizer initialization
    ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad ); //Possible methods: kMigrad, kSimplex, kCombined, kScan, kFumili
    min.SetMaxFunctionCalls(1000000);
    min.SetMaxIterations(100000);
    min.SetTolerance(0.001);
    ROOT::Math::Functor ftr(&FuncToMin,NPARS+2);  //Plus two to acommodate not fitted (constant) parameter nsymm (3 for triangle, 4 for square, -1 for circle) and phi0 (shift of 0 degrees in the plots)
    min.SetFunction(ftr);
    
    // Set the variables to be minimized
    min.SetLowerLimitedVariable(0, "cn", 1, 0.01, 0.0);
    min.SetLimitedVariable(1, "an", 0.0, 0.01, -0.9, 0.9);
    min.SetLimitedVariable(2, "bn", 0.0, 0.01, -0.9, 0.9);
    min.SetFixedVariable(3,"nsymm",nsymm);
    min.SetFixedVariable(4,"phi0",phi0);
    //min.FixVariable(1);
    if(nsymm==-1) { min.FixVariable(1); min.FixVariable(2); } //Circular pattern for circular data
    if(icurve==0) min.FixVariable(2); //No higher order corrections for the first curve

    //Fixations for some of the square plots from the 2022 Shu paper
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig4" && icurve==3) min.FixVariable(2);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig5" && icurve==3) min.FixVariable(2);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig6" && icurve==4) min.FixVariable(2);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig7" && icurve==5) min.FixVariable(2);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig8" && icurve==1) min.FixVariable(2);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig8" && icurve==3) min.FixVariable(2);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig8" && icurve==5) min.FixVariable(2);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig9" && icurve==5) min.FixVariable(2);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig10" && icurve==5) min.FixVariable(2);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig11" && icurve==6) min.FixVariable(2);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig12" && icurve==5) min.FixVariable(2);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig12" && icurve==6) min.FixVariable(2);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig12" && icurve==7) min.FixVariable(2);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig13" && icurve==6) min.FixVariable(1);
    if(outfilebase=="shu_plots_and_data/shu2022_wpd_fig13" && icurve==6) min.FixVariable(2);
    
    //Fixations for some of the triangle plots from the 2022 Shu paper
    if(outfilebase=="triangles_shu/shuTR01_wpd" && icurve==1) min.FixVariable(2);
    if(outfilebase=="triangles_shu/shuTR01_wpd" && icurve==2) min.FixVariable(2);
    if(outfilebase=="triangles_shu/shuTR04_wpd" && icurve==3) min.FixVariable(2);
    
    //Perform minimalization and error calculation
    min.Minimize(); 
    min.PrintResults();
    min.Hesse();
    
    //Save parameters for later
    cnval.at(icurve) = min.X()[0];
    anval.at(icurve) = min.X()[1];
    bnval.at(icurve) = min.X()[2];
    cnerr.at(icurve) = min.Errors()[0];
    anerr.at(icurve) = min.Errors()[1];
    bnerr.at(icurve) = min.Errors()[2];
    
    //Create Graph object from fitted curve
    int Nphi = 500;
    vector<double> phicalcvector(Nphi,0), rcalcvector(Nphi,0);
    for(int ip = 0; ip<Nphi; ip++)
    {
      double phi = (ip+0.5)*2*M_PI/Nphi;
      phicalcvector.at(ip) = phi;
      rcalcvector.at(ip) = FitFunc(&phi,min.X());
    } 
    phi_vs_r_fit.push_back(new TGraphPolar(phicalcvector.size(),&phicalcvector[0],&rcalcvector[0]));
  }
  //Creating a canvas to draw on
  TCanvas *c = new TCanvas("c","c",1024,1024);

  for(int icurve=0; icurve<Ncurves; icurve++)
  { 
    //Draw data points
    double rmax = 1.5;
    if(nsymm==3) rmax=0.7; //triangle plots are smaller
    phi_vs_r_data_plot.at(icurve)->SetTitle(infilename.c_str());
    phi_vs_r_data_plot.at(icurve)->Draw("P");
    c->Update();
    phi_vs_r_data_plot.at(icurve)->SetMinPolar(0);
    phi_vs_r_data_plot.at(icurve)->SetMaxPolar(2*M_PI);
    phi_vs_r_data_plot.at(icurve)->SetMinRadial(0);
    phi_vs_r_data_plot.at(icurve)->SetMaxRadial(rmax);
    phi_vs_r_data_plot.at(icurve)->SetMarkerStyle(20);
    phi_vs_r_data_plot.at(icurve)->SetMarkerSize(2.);
    phi_vs_r_data_plot.at(icurve)->SetMarkerColor(colors[icurve]);
    phi_vs_r_data_plot.at(icurve)->Draw("P");
    
    //Draw fitted curve
    phi_vs_r_fit.at(icurve)->SetLineColor(colors[icurve]);
    phi_vs_r_fit.at(icurve)->SetLineWidth(4);
    phi_vs_r_fit.at(icurve)->Draw("L");
  }
  //Save image
  c->Print(Form("%s_all_fit.png",outfilebase.c_str()));
  
  //Create overlaid version -- NOT WORKING!
  //TImage *img = TImage::Open(Form("%s.png",outfilebase.c_str()));
  //img->Draw("x");
  //phi_vs_r_data_plot.at(icurve)->Draw("P");
  //c->Print(Form("%s_overlay_try.png",outfilebase.c_str()));
  
  //Print parameters
  ofstream outfile(Form("%s_fits.out",outfilebase.c_str()),ofstream::app);
  outfile << outfilebase << " an";  for(int icurve=0; icurve<Ncurves; icurve++) outfile << "\t" << anval.at(icurve) << "\t" << anerr.at(icurve); outfile << endl;
  outfile << outfilebase << " bn";  for(int icurve=0; icurve<Ncurves; icurve++) outfile << "\t" << bnval.at(icurve) << "\t" << bnerr.at(icurve); outfile << endl;
  outfile << outfilebase << " cn";  for(int icurve=0; icurve<Ncurves; icurve++) outfile << "\t" << cnval.at(icurve) << "\t" << cnerr.at(icurve); outfile << endl;
  
  //Clean up and return
  delete phi_vs_r_data;
  for(int icurve=0; icurve<Ncurves; icurve++)
    delete phi_vs_r_fit.at(icurve);
  delete c;
  return 0;
} 
