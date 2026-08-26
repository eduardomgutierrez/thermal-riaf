#include "adafParameters.h"
#include "adafFunctions.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <fmath/physics.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>
#include <fparameters/SpaceIterator.h>
#include "State.h"

using namespace std;

///////////////////////

void adafParameters() 
{	
	// Check if we should read precomputed ADAF from sol_spin file FIRST
	readPrecomputedADAF = GlobalConfig.get<int>("readPrecomputedADAF", 0);
	R_high = GlobalConfig.get<double>("R_high", 1.0);
	precompPhiIsOmega = 0;
	
	ifstream adafFile,adafParams;
	double accRateNorm = 0.0;
	
	if (readPrecomputedADAF) {
		// Read black hole mass and spin from parameters.json
		blackHoleMass = GlobalConfig.get<double>("blackHoleMass", 10.0);  // Solar masses
		blackHoleSpin = GlobalConfig.get<double>("spin", 0.0);
		cout << "Using precomputed ADAF mode with M = " << blackHoleMass << " Msun, a = " << blackHoleSpin << endl;
		
		// Read sol_spin file
		string solSpinFile = GlobalConfig.get<string>("solSpinFile", "sol_spin_p50.dats");
		ifstream solFile(solSpinFile);
		if (!solFile.is_open()) {
			cerr << "Error: Cannot open sol_spin file: " << solSpinFile << endl;
			readPrecomputedADAF = 0;
		} else {
			// Skip header line
			string header;
			getline(solFile, header);
			
			// Columns: r, density (g/cm^3), vel/c, b^r (G), b^phi (G), B (G), H/r, beta, Temp (K), Omega, v_phi/c
			precomp_r.clear();
			precomp_density.clear();
			precomp_v.clear();
			precomp_vphi.clear();
			precomp_B.clear();
			precomp_HR.clear();
			precomp_beta.clear();
			precomp_Temp.clear();
			
			double r_val, dens_val, v_val, br_val, bphi_val, B_val, HR_val, beta_val, Temp_val, omega_val, vphi_val;
			while (solFile >> r_val >> dens_val >> v_val >> br_val >> bphi_val >> B_val >> HR_val >> beta_val >> Temp_val >> omega_val >> vphi_val) {
				precomp_r.push_back(r_val);
				precomp_density.push_back(dens_val);
				precomp_v.push_back(v_val);
				precomp_vphi.push_back(vphi_val);
				precomp_B.push_back(B_val);
				precomp_HR.push_back(HR_val);
				precomp_beta.push_back(beta_val);
				precomp_Temp.push_back(Temp_val);
			}
			solFile.close();
			cout << "Read " << precomp_r.size() << " points from precomputed ADAF file: " << solSpinFile << endl;
			cout << "Using precomputed azimuthal profile from v_phi/c column" << endl;
			cout << "Using R_high = " << R_high << " for Te = Ti/R prescription" << endl;
		}
		
		// Set unused parameters to safe defaults (not used in precomputed mode)
		s = 0.0;
		magFieldPar = 0.5;  // Will be overridden by precompBeta locally
		alpha = 0.1;
		j = 1.0;
		delta = 0.3;
		rTr = 3.0;
		powerIndex = 0.1;
		
		// Create minimal fallback arrays for logr/logTi/logTe/logv
		// These are only used if precomputed interpolation fails
		nRaux = 2;
		logr.resize(nRaux, 0.0);
		logTi.resize(nRaux, 0.0);
		logTe.resize(nRaux, 0.0);
		logv.resize(nRaux, 0.0);
		// Use precomputed range for fallback
		double r_inner_rg = std::max(precomp_r.back(), 2.0);
		double r_outer_rg = precomp_r.front();
		logr[0] = log(r_inner_rg / 2.0);  // Convert r_g to R_S (r_g = R_S/2)
		logr[1] = log(r_outer_rg / 2.0);
		logTi[0] = log(precomp_Temp.back());
		logTi[1] = log(precomp_Temp.front());
		logTe[0] = log(precomp_Temp.back());
		logTe[1] = log(precomp_Temp.front());
		logv[0] = log(std::abs(precomp_v.back()));
		logv[1] = log(std::abs(precomp_v.front()));
	} else {
		// Original mode: read from adafParameters.txt
		// For non-precomputed case, spin is 0 (Schwarzschild)
		blackHoleSpin = 0.0;
		adafFile.open("adafFile.txt"); adafParams.open("adafParameters.txt");
		adafFile >> nRaux;
		adafParams >> blackHoleMass >> accRateNorm >> s >> magFieldPar >> alpha >> j >> delta >> rTr >> powerIndex;
		adafParams.close();
		
		logr.resize(nRaux,0.0);
		logTi.resize(nRaux,0.0);
		logTe.resize(nRaux,0.0);
		logv.resize(nRaux,0.0);
		
		for (size_t i=0;i<nRaux;i++) {
			adafFile >> logr[i] >> logTi[i] >> logTe[i] >> logv[i];
		}
		adafFile.close();
	}
    
	blackHoleMass *= solarMass;
	schwRadius = 2.0*gravitationalConstant*blackHoleMass / cLight2;
	gravRadius = schwRadius / 2.0;  // r_g = GM/c^2 = R_S/2
	
	// Compute horizon radius: r_H = r_g * (1 + sqrt(1 - a^2)) for Kerr black hole
	// For a=0 (Schwarzschild): r_H = 2*r_g = R_S
	double a = blackHoleSpin;
	if (std::abs(a) >= 1.0) a = 0.998 * (a > 0 ? 1.0 : -1.0);  // Limit to near-extremal
	horizonRadius = gravRadius * (1.0 + sqrt(1.0 - a*a));
	cout << "Gravitational radius r_g = " << gravRadius << " cm" << endl;
	cout << "Horizon radius r_H = " << horizonRadius << " cm = " << horizonRadius/gravRadius << " r_g" << endl;
	
	if (readPrecomputedADAF && !precomp_r.empty()) {
		// Compute accRateOut from precomputed data at the outer boundary
		// Mdot = 4*pi*r*H*rho*|v|  at outer radius
		double r_out_rg = precomp_r.front();  // Outer radius in r_g
		double r_out_cm = r_out_rg * gravRadius;  // Convert to cm
		double rho_out = precomp_density.front();  // g/cm^3
		double v_out = std::abs(precomp_v.front()) * cLight;  // cm/s
		double HR_out = precomp_HR.front();
		double H_out = r_out_cm * HR_out;  // cm
		accRateOut = 4.0 * pi * r_out_cm * H_out * rho_out * v_out;  // g/s
		
		// Compute accRateIn from precomputed data at the inner boundary
		double r_in_rg = precomp_r.back();  // Inner radius in r_g
		double r_in_cm = r_in_rg * gravRadius;  // Convert to cm
		double rho_in = precomp_density.back();  // g/cm^3
		double v_in = std::abs(precomp_v.back()) * cLight;  // cm/s
		double HR_in = precomp_HR.back();
		double H_in = r_in_cm * HR_in;  // cm
		double accRateIn = 4.0 * pi * r_in_cm * H_in * rho_in * v_in;  // g/s
		
		double eddAccRate = 1.39e18 * blackHoleMass/solarMass;
		cout << "Computed accretion rates from precomputed data:" << endl;
		cout << "  At r_out = " << r_out_rg << " r_g: Mdot_out = " << accRateOut << " g/s" 
		     << " (" << accRateOut/eddAccRate << " Mdot_Edd)" << endl;
		cout << "  At r_in  = " << r_in_rg << " r_g: Mdot_in  = " << accRateIn << " g/s"
		     << " (" << accRateIn/eddAccRate << " Mdot_Edd)" << endl;
		cout << "Outer accretion power = " << accRateOut*cLight2 << " erg s^{-1}" << endl;
		cout << "Inner accretion power = " << accRateIn*cLight2 << " erg s^{-1}" << endl;
	} else {
		double eddAccRate = 1.39e18 * blackHoleMass/solarMass;
		accRateOut = accRateNorm*eddAccRate;
		cout << "Total outer accretion power = " << accRateOut*cLight2 << " erg s^{-1}" << endl;
		cout << "Total inner accretion power = " << accRateOut*pow(1.0/exp(logr.back()),s)*cLight2 
			 << " erg s^{-1}" << endl;
	}
    
	rTr = rTr * schwRadius;
	rOutCD = GlobalConfig.get<double>("rOutCD") * schwRadius;
	eMeanMolecularWeight = GlobalConfig.get<double>("mu_e");
	iMeanMolecularWeight = GlobalConfig.get<double>("mu_i");

	nR = GlobalConfig.get<int>("model.particle.default.dim.radius.samples");
	nE = GlobalConfig.get<int>("model.particle.photon.dim.energy.samples");
	nRcd = GlobalConfig.get<int>("model.particle.default.dim.radius_cd.samples");
    
	// Calculate radial step - use precomputed range if available
	if (readPrecomputedADAF && !precomp_r.empty()) {
		// precomp_r is ordered from large to small, so front() is outer, back() is inner
		// Enforce inner radius >= horizon radius (r_H/r_g depends on spin)
		double r_H_rg = horizonRadius / gravRadius;  // Horizon radius in units of r_g
		double r_inner_rg = std::max(precomp_r.back(), r_H_rg*1.01);
		paso_r = pow(precomp_r.front()/r_inner_rg, 1.0/(nR));
		cout << "Using sol_spin radial range: r_in = " << r_inner_rg << " r_g, r_out = " 
			 << precomp_r.front() << " r_g" << endl;
		cout << "Horizon at r_H = " << r_H_rg << " r_g" << endl;
	} else {
		paso_r = pow(exp(logr.back())/exp(logr.front()),1.0/(nR));
	}
	paso_rCD = pow(rOutCD/rTr,1.0/(nRcd));

	logMinEnergy = GlobalConfig.get<double>("model.particle.photon.dim.energy.min");
	logMaxEnergy = GlobalConfig.get<double>("model.particle.photon.dim.energy.max");

	inclination = GlobalConfig.get<double>("inclination");
    
    calculateComptonScatt = GlobalConfig.get<int>("calculateComptonScatt");
	height_method = GlobalConfig.get<int>("height_method");
    
    calculateThermal = GlobalConfig.get<int>("calculateThermal");
	numProcesses = GlobalConfig.get<int>("thermal.numProcesses");
    if (calculateThermal) {
        calculateComptonRedMatrix = GlobalConfig.get<int>("thermal.compton.calculateRedMatrix");
        comptonMethod = GlobalConfig.get<int>("thermal.compton.method");
        if (1) {

            nGammaCompton = GlobalConfig.get<size_t>("thermal.compton.redMatrixParams.nGammaCompton");
            nTempCompton = GlobalConfig.get<size_t>("thermal.compton.redMatrixParams.nTempCompton");
            nNuPrimCompton = GlobalConfig.get<size_t>("thermal.compton.redMatrixParams.nNuPrimCompton");
            nNuCompton = GlobalConfig.get<size_t>("thermal.compton.redMatrixParams.nNuCompton");
            gammaMinCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.gammaMinCompton");
            gammaMaxCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.gammaMaxCompton");
            tempMinCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.tempMinCompton");
            tempMaxCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.tempMaxCompton");
            nuPrimMinCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.nuPrimMinCompton");
            nuPrimMaxCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.nuPrimMaxCompton");
            nuMinCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.nuMinCompton");
            nuMaxCompton = GlobalConfig.get<double>("thermal.compton.redMatrixParams.nuMaxCompton");
            
            ofstream fileSizes;
            fileSizes.open("comptonRedMatrix/sizesVecCompton.dat",ios::out);
            fileSizes << nGammaCompton << "\t" << gammaMinCompton << "\t" << gammaMaxCompton << endl
                      << nTempCompton << "\t" << tempMinCompton << "\t" << tempMaxCompton << endl
                      << nNuPrimCompton << "\t" << nuPrimMinCompton << "\t" << nuPrimMaxCompton << endl
                      << nNuCompton << "\t" << nuMinCompton << "\t" << nuMaxCompton << endl;
            fileSizes.close();
        } else {
            ifstream fileSizes;
            fileSizes.open("comptonRedMatrix/sizesVecCompton.dat",ios::in);
            fileSizes >> nGammaCompton >> gammaMinCompton >> gammaMaxCompton
                      >> nTempCompton >> tempMinCompton >> tempMaxCompton
                      >> nNuPrimCompton >> nuPrimMinCompton >> nuPrimMaxCompton
                      >> nNuCompton >> nuMinCompton >> nuMaxCompton;
            fileSizes.close();
        }
        calculatePhotonDensityGap = GlobalConfig.get<int>("thermal.calculatePhotonDensityGap");
    }
	calculateJetEmission = GlobalConfig.get<int>("nonThermal.calculateJetEmission");
	
    calculateNonThermal = GlobalConfig.get<int>("calculateNonThermal");
    if (calculateNonThermal) {
		accMethod = GlobalConfig.get<double>("nonThermal.acc_method");
		calculateNTprotons = GlobalConfig.get<int>("nonThermal.protons");
		calculateNTelectrons = GlobalConfig.get<int>("nonThermal.electrons");
        calculateLosses = GlobalConfig.get<int>("nonThermal.calculateLosses");
		calculateFlare = GlobalConfig.get<int>("nonThermal.calculateFlare");
        calculateNTdistributions = GlobalConfig.get<int>("nonThermal.calculateDistributions");
        calculateNonThermalLum = GlobalConfig.get<int>("nonThermal.calculateLuminosities");
		calculateNonThermalHE = GlobalConfig.get<int>("nonThermal.calculateHighEnergyProcesses");
        calculateNeutronInj = GlobalConfig.get<int>("nonThermal.neutrons.calculateInjection");
		calculateNeutronDis = GlobalConfig.get<int>("nonThermal.neutrons.calculatePropagation");
		calculateJetDecay = GlobalConfig.get<int>("nonThermal.neutrons.calculateJetDecay");
		calculateSecondaries = GlobalConfig.get<int>("nonThermal.calculateSecondaries");
		calculateNeutrinos = GlobalConfig.get<int>("nonThermal.calculateNeutrinos");
	}
}