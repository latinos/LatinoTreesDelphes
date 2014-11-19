/****************************************************************************************************************************

        -this program create a sample on minimum bias events for SLHC-14TeV
        -reference for the Tune:pp option used can be found here: https://cds.cern.ch/record/1697700/files/GEN-14-001-pas.pdf 

        -source minbias_setup_slc6, then:
        -compile with ---> gcc -I$HEPMC/include/ -L$HEPMC/lib -I$PYTHIA8DATA/../include -L$PYTHIA8DATA/../lib -I$LHAPDF/include -L$LHAPDF/lib 
                               -lLHAPDF -lHepMC -lpythia8tohepmc -lpythia8 -o genMinBias_14TeV genMinBias_14TeV.cpp
****************************************************************************************************************************/

#include "Pythia8/Pythia.h"
#include "Pythia8/Pythia8ToHepMC.h"
#include "Pythia8/LHAPDFInterface.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/IO_AsciiParticles.h"

using namespace Pythia8;  
using namespace std;

//********************************************************************************************

int main(int argc, char** argv) 
{
    if(argc < 3)
    {
        cout << "missing argument!" << endl
             << "-----------------" << endl
             << "genMinBias_14TeV events_number output_file" << endl
             << "-----------------" << endl << endl;
        return 0;
    }
    int nEvents = atoi(argv[1]);
    string outFile = argv[2];
    
    //-----Pythia setup-----
    Pythia pythia;     
    //---general setup
    pythia.readString("Random:seed = 0");
    pythia.readString("Random:setSeed = on"); 
    pythia.readString("Beams:idA = 2212");
    pythia.readString("Beams:idB = 2212");
    pythia.readString("Beams:eA = 7000.");
    pythia.readString("Beams:eB = 7000.");
    pythia.readString("Beams:eCM = 14000.");
    pythia.readString("HadronLevel:Hadronize = on");
    //---CMS tuned min bias
    pythia.readString("ParticleDecays:limitTau0 = on");
    pythia.readString("ParticleDecays:tauMax = 10");
    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("SoftQCD:singleDiffractive = on");
    pythia.readString("SoftQCD:doubleDiffractive = on");
    pythia.readString("Tune:ee = 3");
    pythia.readString("Tune:pp = 2");
    //pythia.readString("Tune:pp = 15");
    //---Pythia initialization
    pythia.init();

    HepMC::IO_GenEvent hepmc_file_out(outFile, std::ios::out);
    HepMC::Pythia8ToHepMC ToHepMC;
	
    for (int iEvent = 0; iEvent < nEvents; ++iEvent)
    {
	if(iEvent%(nEvents/10) == 0) 
	{
	    cout << "Events:  " << iEvent << endl;
	}
	if(!pythia.next())
	{
	    cout << "CRASH! ---> skip" << endl;
	    continue;
	}
	// construct new HepMC event setting units.	
	HepMC::GenEvent* hepmcevt = new HepMC::GenEvent(HepMC::Units::GEV, HepMC::Units::MM);
	// fill the event including PDF infos
	ToHepMC.fill_next_event( pythia, hepmcevt );
	hepmc_file_out << hepmcevt;
        delete hepmcevt;
    }
    // Done.                           
    return 0;
}
