#include <TSpectrum.h>
#include <TFile.h>
#include <TTree.h>
#include <TBrowser.h>
#include <TRandom.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>

#define PI 3.14159

void Pi0TreeTask(){

//Opening the file
TFile *file = TFile::Open("/home/rbaker/2014work/Pi0TreeTask/Pi0_mc_out_250.root");
if(file == 0){
	cout << "Error: File not opened" << endl;
	return;
	}
//Setting up the tree
TTree *tree = (TTree *) file->Get("h1");

//Defining all variables and arrays (arrays not needed)
float X_vtx, Y_vtx, Z_vtx, Px_bm, Py_bm, Pz_bm, Pt_bm, En_bm, Px_l0114, Py_l0114, Pz_l0114, Pt_l0114, En_l0114, Px_l0201, Py_l0201, Pz_l0201, Pt_l0201, En_l0201, Px_l0301, Py_l0301, Pz_l0301, Pt_l0301, En_l0301, phiRAD0, phiRAD1, phi0, phi1, thetaRAD0, thetaRAD1, theta0, theta1, ph0_ph1_mass, proton_rest_mass;

/*float X_vtx[10000], Y_vtx[10000], Z_vtx[10000], Px_bm[10000], Py_bm[10000], Pz_bm[10000], Pt_bm[10000], En_bm[10000], Px_l0114[10000], Py_l0114[10000], Pz_l0114[10000], Pt_l0114[10000], En_l0114[10000], Px_l0201[10000], Py_l0201[10000], Pz_l0201[10000], Pt_l0201[10000], En_l0201[10000], Px_l0301[10000], Py_l0301[10000], Pz_l0301[10000], Pt_l0301[10000], En_l0301[10000];
*/
//Setting up the branches
tree->SetBranchAddress("X_vtx", &X_vtx);
tree->SetBranchAddress("Y_vtx", &Y_vtx);
tree->SetBranchAddress("Z_vtx", &Z_vtx);
tree->SetBranchAddress("Px_bm", &Px_bm);
tree->SetBranchAddress("Py_bm", &Py_bm);
tree->SetBranchAddress("Pz_bm", &Pz_bm);
tree->SetBranchAddress("Pt_bm", &Pt_bm);
tree->SetBranchAddress("En_bm", &En_bm);
tree->SetBranchAddress("Px_l0114", &Px_l0114);
tree->SetBranchAddress("Py_l0114", &Py_l0114);
tree->SetBranchAddress("Pz_l0114", &Pz_l0114);
tree->SetBranchAddress("Pt_l0114", &Pt_l0114);
tree->SetBranchAddress("En_l0114", &En_l0114);
tree->SetBranchAddress("Px_l0201", &Px_l0201);
tree->SetBranchAddress("Py_l0201", &Py_l0201);
tree->SetBranchAddress("Pz_l0201", &Pz_l0201);
tree->SetBranchAddress("Pt_l0201", &Pt_l0201);
tree->SetBranchAddress("En_l0201", &En_l0201);
tree->SetBranchAddress("Px_l0301", &Px_l0301);
tree->SetBranchAddress("Py_l0301", &Py_l0301);
tree->SetBranchAddress("Pz_l0301", &Pz_l0301);
tree->SetBranchAddress("Pt_l0301", &Pt_l0301);
tree->SetBranchAddress("En_l0301", &En_l0301);

//Setting up the vectors
TLorentzVector v_beam, v_proton, v_proton_recoil, v_photon0, v_photon1, v_photon_add, v_proton_rest;

//Old code used to analyze the exact output of each variable
/*

tree->GetEntry(0);
cout << "X_vtx: " << X_vtx[0] << endl;
cout << "Y_vtx: " << Y_vtx[0] << endl;
cout << "Z_vtx: " << Z_vtx[0] << endl;
cout << "Px_bm: " << Px_bm[0] << endl;
cout << "Py_bm: " << Py_bm[0] << endl;
cout << "Pz_bm: " << Pz_bm[0] << endl;
cout << "Pt_bm: " << Pt_bm[0] << endl;
cout << "En_bm: " << En_bm[0] << endl;
cout << "Px_proton: " << Px_l0114[0] << endl;
cout << "Py_proton: " << Py_l0114[0] << endl;
cout << "Pz_proton: " << Pz_l0114[0] << endl;
cout << "Pt_proton: " << Pt_l0114[0] << endl;
cout << "En_proton: " << En_l0114[0] << endl;
cout << "Px_photon0: " << Px_l0201[0] << endl;
cout << "Py_photon0: " << Py_l0201[0] << endl;
cout << "Pz_photon0: " << Pz_l0201[0] << endl;
cout << "Pt_photon0: " << Pt_l0201[0] << endl;
cout << "En_photon0: " << En_l0201[0] << endl;
cout << "Px_photon1: " << Px_l0301[0] << endl;
cout << "Py_photon1: " << Py_l0301[0] << endl;
cout << "Pz_photon1: " << Pz_l0301[0] << endl;
cout << "Pt_photon1: " << Pt_l0301[0] << endl;
cout << "En_photon1: " << En_l0301[0] << endl;
*/

//Setting up the graphs
TH1D *graph0 = new TH1D("En_bm", "Energy of the incoming photon beam (GeV)", 200, 0.24, 0.26);
TH1D *graph1 = new TH1D("ph0_ph1_mass", "Photon Invariant Mass Reconstruction", 500, 0, 220);
TH1D *graph2 = new TH1D("proton_mass", "Proton Invariant Mass Reconstruciton", 500, 900, 1000);
TH1D *graph3 = new TH1D("theta", "Theta for the two Photons", 180, 0, 180);
TH1D *graph4 = new TH1D("phi", "Phi for the two Photons", 180, 0, 180);

//Retrieve the number of entries in the tree
int N_entries = (int) tree->GetEntries();
//Loop over each entry
for(int i = 0; i < N_entries; i++){
	//Get the individual entry
	tree->GetEntry(i);
	//Fill graph0 with the energy info
	graph0->Fill(En_bm);
	//Set the values for the beam, proton, recoil proton, and 2 photon vectors
	double P;
	P = Pt_bm*1000.0;
	v_beam.SetXYZM(Px_bm*P, Py_bm*P, Pz_bm*P, 0.0);
	v_proton.SetXYZM(0., 0., 0., 938.27);
	P = Pt_l0114*1000.0;
	v_proton_recoil.SetXYZM(Px_l0114*P, Py_l0114*P, Pz_l0114*P, 938.27);
	P = Pt_l0201*1000.0;
	v_photon0.SetXYZM(Px_l0201*P, Py_l0201*P, Pz_l0201*P, 0.0);
	P = Pt_l0301*1000.0;
	v_photon1.SetXYZM(Px_l0301*P, Py_l0301*P, Pz_l0301*P, 0.0);

	//Add the 2 photons together to find the invariant mass, then fill the graph with the result
	v_photon_add = v_photon0 + v_photon1;
	ph0_ph1_mass = v_photon_add.M();
	if(ph0_ph1_mass < 0){
		ph0_ph1_mass = -1*ph0_ph1_mass;
		}
	graph1->Fill(ph0_ph1_mass);

        //Calculating phi and theta for photon0, then filling the graph
	thetaRAD0 = v_photon0.Theta();
	theta0 = thetaRAD0*180/PI;
	graph3->Fill(theta0);
	phiRAD0 = v_photon0.Phi();
	phi0 = phiRAD0*180/PI;
	graph4->Fill(phi0);

	//Calculating phi and theta for photon1, then filling the graph
	thetaRAD1 = v_photon1.Theta();
        theta1 = thetaRAD1*180/PI;
        graph3->Fill(theta1);
        phiRAD1 = v_photon1.Phi();
        phi1 = phiRAD1*180/PI;
        graph4->Fill(phi1);

	//Calculating proton rest mass
	v_proton_rest = v_beam + v_proton - v_photon_add;
	proton_rest_mass = v_proton_rest.M();
	if(proton_rest_mass < 0){
		proton_rest_mass = -1*proton_rest_mass;
		}
	graph2->Fill(proton_rest_mass);
	}

//Drawing the graphs
c0 = new TCanvas("c0", "c0", 700, 1300);
graph0->Draw(); 
c1 = new TCanvas("c1", "c1", 700, 1300);
graph1->Draw();
c2 = new TCanvas("c2", "c2", 700, 1300);
graph2->Draw();
c3 = new TCanvas("c3", "c3", 700, 1300);
graph3->Draw();
c4 = new TCanvas("c4", "c4", 700, 1300);
graph4->Draw();

//Using TSpectrum to output the exact peak position of the mass graphs
TSpectrum *s1 = new TSpectrum();
double nfound1 = s1->Search(graph2, 1, "goff", 0.7);
float *xpeaks1 = s1->GetPositionX();
for(int p=0; p<nfound1; p++){
	float xp1 = xpeaks1[p];
       	int bin1 = graph2->GetXaxis()->FindBin(xp1);
        float yp1 = graph2->GetBinContent(bin1);
	cout << "Proton Invariant Mass Peak: " << xp1 << " MeV" << endl;
        }

TSpectrum *s2 = new TSpectrum();
double nfound2 = s2->Search(graph1, 1, "goff", 0.7);
float *xpeaks2 = s2->GetPositionX();
for(int p=0; p<nfound2; p++){
	float xp2 = xpeaks2[p];
	int bin2 = graph1->GetXaxis()->FindBin(xp2);
	float yp2 = graph1->GetBinContent(bin2);
	cout << "Pion Mass Reconstruction Peak: " << xp2 << " MeV" << endl;
	}                                         

}
