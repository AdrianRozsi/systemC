
//g++ -DSC_INCLUDE_FX -w -g -I. -I/sysc/systemc-2.3.2/include -L. -L/sysc/systemc-2.3.2/lib-linux64 -o main *.cpp -lsystemc -lm


#include <systemc>
#include "vp.hpp"
#include "hard.hpp"
#include "soft.hpp"
#include <iostream>
#include <fstream>
#include "AudioFile.h"
#include "interconnect.hpp"
#include <tlm_utils/simple_initiator_socket.h>
#include <tlm_utils/simple_target_socket.h>
extern int thr;
#include <chrono>
using namespace sc_core;
using namespace tlm;
using namespace std::chrono;  // Biztosítjuk, hogy elérhetők legyenek az időzítési funkciók
using namespace std;

sc_core::sc_time offset, delay;
int counter;
const sc_dt::uint64 BASS = 0;
const sc_dt::uint64 MID = 1;
const sc_dt::uint64 TREBLE = 2;

int resize(int samples);

int sc_main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: ./program someFile.wav preset (bass, mid, treble)" << std::endl;
        return 0;
    }

    std::string wavFile = argv[1];
    std::string preset = argv[2];

    
    Soft soft("soft_instance", wavFile, preset);
    // Declare and instantiate Interconnect if necessary
Interconnect interconnect("interconnect");



// Bind the appropriate ports if Soft has TLM ports

    soft.processAudio();
    
    
    std::cout << "Data written to file successfully." << std::endl;

    // Start the SystemC simulation
    Vp vp("VP", argv[1], argv[2]);
    sc_start(10, sc_core::SC_MS);  // Adjust simulation time as needed

    return 0;
}
