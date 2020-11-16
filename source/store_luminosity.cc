#include "galprop_classes.h"
#include "galprop_internal.h"

#include <ErrorLogger.h>

#include <cassert>
#include <string>
#include <cstring>
#include <valarray>

using namespace std;

#include <CCfits/CCfits>
#include <CCfits/FITS.h>

void Galprop::store_luminosity() {

  INFO("Entry");

  const string fName = "!" + configure.fOutputDirectory + configure.fOutputPrefix + "luminosity_" + galdef.galdef_ID + ".gz";

  CCfits::FITS luminosityFile(fName, CCfits::Write);

  vector<string> columnNamesCR(6, "Energy"), columnFormatCR(6, "D"), columnUnitsCR(6, "erg s^-1");

  columnNamesCR[1] = "Protons";
  columnNamesCR[2] = "Helium Nuclei";
  columnNamesCR[3] = "Other Nuclei";
  columnNamesCR[4] = "Primary Leptons";
  columnNamesCR[5] = "Secondary Leptons";
  
  columnUnitsCR[0] = "MeV";

  CCfits::Table* crTable = luminosityFile.addTable("Differential Cosmic-ray Luminosity", galaxy.n_pgrid, columnNamesCR, columnFormatCR, columnUnitsCR);

  CCfits::Column& crEnergyCol = crTable->column("Energy");
  CCfits::Column& crProtonCol = crTable->column("Protons");
  CCfits::Column& crHeliumCol = crTable->column("Helium Nuclei");
  CCfits::Column& crNucleiCol = crTable->column("Other Nuclei");
  CCfits::Column& crPrimaryLeptonCol = crTable->column("Primary Leptons");
  CCfits::Column& crSecondaryLeptonCol = crTable->column("Secondary Leptons");

  crEnergyCol.write(galaxy.p, 1);
  crProtonCol.write(galaxy.fProtonLuminosity, 1);
  crHeliumCol.write(galaxy.fHeliumLuminosity, 1);
  crNucleiCol.write(galaxy.fNucleiLuminosity, 1);
  crPrimaryLeptonCol.write(galaxy.fPrimaryElectronLuminosity, 1);
  crSecondaryLeptonCol.write(galaxy.fSecondaryElectronLuminosity, 1);

  crTable->addKey("NCRBINS", galaxy.n_pgrid, "Number of cosmic ray energy/momentum bins");
  const double crKeyV = galaxy.p_factor;
  crTable->addKey("CRBINSIZE", crKeyV, "Cosmic ray energy/momentum bin size");

  vector<string> columnNamesPhotons(7, "Energy"), columnFormatPhotons(7, "D"), columnUnitsPhotons(7, "erg s^-1");

  columnNamesPhotons[1] = "Frequency";
  columnNamesPhotons[2] = "Synchrotron";
  columnNamesPhotons[3] = "Interstellar Radiation Field";
  columnNamesPhotons[4] = "Neutral Pion Decay";
  columnNamesPhotons[5] = "Bremsstrahlung";
  columnNamesPhotons[6] = "Inverse Compton";

  columnUnitsPhotons[0] = "eV";
  columnUnitsPhotons[1] = "Hz";

  // Here we combine the information calculated for different frequency/energy
  // ranges where different processes are dominant into a combined SED 
  // spanning synchrotron -> gamma ray energies

  valarray<double> syncEnergy(0., galaxy.nu_synch.size()), gammaEnergy(0., galaxy.n_E_gammagrid), targetEnergy(0., galaxy.nu_ISRF.size());
    
  targetEnergy = h_planck*erg_to_eV*galaxy.nu_ISRF;
  syncEnergy = galaxy.nu_synch*h_planck*erg_to_eV; // Hz -> eV
  gammaEnergy = galaxy.E_gamma*1e6; // MeV -> eV

  vector< pair<double, int> > rangeSED;

  for (int i = 0; i < syncEnergy.size(); ++i)
    rangeSED.push_back(pair<double, int>(syncEnergy[i], 1));
  
  for (int i = 0; i < targetEnergy.size(); ++i)
    rangeSED.push_back(pair<double, int>(targetEnergy[i], 2));

  for (int i = 0; i < gammaEnergy.size(); ++i)
    rangeSED.push_back(pair<double, int>(gammaEnergy[i], 3));

  stable_sort(rangeSED.begin(), rangeSED.end());
  
  //for (int i = 0; i < rangeSED.size(); ++i)
  //cout << i << " " << rangeSED[i].first << " " << rangeSED[i].second << endl;

  vector<double> energy, freq, sync, isrf, pi0, brem, ic;

  for (int i = 0; i < rangeSED.size(); ++i) {

    int syncIndex = -1, isrfIndex = -1, gammaIndex = -1;

    if (1 == rangeSED[i].second) {

      const double* syncVal = std::lower_bound(&syncEnergy[0], &syncEnergy[syncEnergy.size()-1], rangeSED[i].first);

      syncIndex = (syncVal - &syncEnergy[0] >= 0 && syncVal - &syncEnergy[0] < syncEnergy.size() ? syncVal - &syncEnergy[0] : -1);

    }

    if (2 == rangeSED[i].second) {

      const double* isrfVal = std::lower_bound(&targetEnergy[0], &targetEnergy[targetEnergy.size()-1], rangeSED[i].first);

      isrfIndex = (isrfVal - &targetEnergy[0] >= 0 && isrfVal - &targetEnergy[0] < targetEnergy.size() ? isrfVal - &targetEnergy[0] : -1);

    }

    if (3 == rangeSED[i].second) {

      const double* gammaVal = std::lower_bound(&gammaEnergy[0], &gammaEnergy[gammaEnergy.size()-1], rangeSED[i].first);

      gammaIndex = (gammaVal - &gammaEnergy[0] >= 0 && gammaVal - &gammaEnergy[0] < gammaEnergy.size() ? gammaVal - &gammaEnergy[0] : -1);

    }

    energy.push_back(rangeSED[i].first);
    freq.push_back(energy.back()/(utl::kPlanck_SI/utl::e_SI));
    sync.push_back(0.);
    isrf.push_back(0.);
    pi0.push_back(0.);
    brem.push_back(0.);
    ic.push_back(0.);
 
    if (1 == rangeSED[i].second)
      sync.back() = galaxy.fSynchrotronLuminosity[syncIndex];

    if (2 == rangeSED[i].second)
      isrf.back() = galaxy.fISRFLuminosity[isrfIndex];

    if (3 == rangeSED[i].second) {

      pi0.back() = galaxy.fPi0DecayLuminosity[gammaIndex];
      brem.back() = galaxy.fBremsstrahlungLuminosity[gammaIndex];
      ic.back() = galaxy.fICLuminosity[gammaIndex];
      
    }

    // Below deals with case when we have entries for two processes for the same energy

    if (i < rangeSED.size()-1) 
      if (rangeSED[i].first == rangeSED[i+1].first) {

	if (1 == rangeSED[i+1].second)
	  sync.back() = galaxy.fSynchrotronLuminosity[syncIndex];
	
	if (2 == rangeSED[i+1].second)
	  isrf.back() = galaxy.fISRFLuminosity[isrfIndex];
	
	if (3 == rangeSED[i+1].second) {
	  
	  pi0.back() = galaxy.fPi0DecayLuminosity[gammaIndex];
	  brem.back() = galaxy.fBremsstrahlungLuminosity[gammaIndex];
	  ic.back() = galaxy.fICLuminosity[gammaIndex];
	  
	}

      }

    //cout << i << " " << rangeSED[i].first << " " << rangeSED[i].second << " " << syncIndex << " " << isrfIndex << " " << gammaIndex << endl;

    //cout << i << " " << rangeSED[i].first << " " << energy.back() << " " << freq.back() << " " << sync.back() << " " << isrf.back() << " " << pi0.back() << " " << brem.back() << " " << ic.back() << endl;

  }

  CCfits::Table* photonTable = luminosityFile.addTable("Differential Photon Luminosity", energy.size(), columnNamesPhotons, columnFormatPhotons, columnUnitsPhotons);

  CCfits::Column& energyCol = photonTable->column("Energy");
  CCfits::Column& freqCol = photonTable->column("Frequency");
  CCfits::Column& syncCol = photonTable->column("Synchrotron");
  CCfits::Column& isrfCol = photonTable->column("Interstellar Radiation Field");
  CCfits::Column& pi0Col = photonTable->column("Neutral Pion Decay");
  CCfits::Column& bremCol = photonTable->column("Bremsstrahlung");
  CCfits::Column& icCol = photonTable->column("Inverse Compton");

  energyCol.write(energy, 1);
  freqCol.write(freq, 1);
  syncCol.write(sync, 1);
  isrfCol.write(isrf, 1);
  pi0Col.write(pi0, 1);
  bremCol.write(brem, 1);
  icCol.write(ic, 1);
  
  photonTable->addKey("NPHOTONBINS", energy.size(), "Number of photon energy bins");

  //cout << *crTable << " " << crTable->numCols() << endl;
  //cout << *photonTable << " " << photonTable->numCols() << endl;

  //exit(0);

  INFO("Exit");

}
