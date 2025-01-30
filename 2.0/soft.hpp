#ifndef SOFT_HPP_
#define SOFT_HPP_
#include "AudioFile.h"
#include <iostream>
#include <fstream>
#include <systemc>
#include <tlm_utils/simple_initiator_socket.h>
#include <stdio.h>
#include <stdlib.h>
#include "typedefs.hpp"
#include "utils.hpp"
#include <complex>

using Complex = std::complex<double>;




class Soft : public sc_core::sc_module
{
public:
    void setParameters(const std::string& wavFile, const std::string& preset) {
        this->wavFile = wavFile;
        this->preset = preset;
    }



    std::string wavFile;   // WAV fájl neve
    std::string preset;    // Kiválasztott preset (BASS, MID, TREBLE)

   Soft(sc_core::sc_module_name name, const std::string& wavFile, const std::string& preset);
    ~Soft();  // Destruktor deklaráció
  void processAudio();

  tlm_utils::simple_initiator_socket<Soft> interconnect_socket;
  
  void IFFT_Function(); // IFFT számítást végző függvény
  vector<Complex> IFFT(const vector<Complex>& in); // IFFT algoritmus
  
  void FFT_Function();
  vector<Complex> FFT(const vector<Complex>& in);
   

  vector<double> getFreq(int n, int sampleRate);

  vector<sc_dt::uint64> gauss_freq;
  vector<vector<int>> gauss_amp;
  vector<double> generate_gaussian(vector<double> freq, uint8_t preset);

  vector<Complex> normalizePeak(vector<Complex> fftArr, uint8_t preset);
  vector<Complex> normalizePower(vector<Complex> modifiedFftArr, vector<Complex> fftArr);
  char* coeffs;

protected:
   void gen();
   //void gen(); //main software function
   pl_t pl; //payload
   sc_core::sc_time offset; //time
   //std::vector<num_t> coeffs(COEFFS_SIZE, 0);
   std::ifstream coeffsFile; 
   std::ifstream samplesFile; //input files
   std::ifstream audioFile; //input files - dodato
   std::ofstream output; 
     num_t result[10*SAMPLES_SIZE];
  num_t samplesH[10*SAMPLES_SIZE];
  num_t coeffsH[COEFFS_SIZE];
  num_t fft_res[COEFFS_SIZE];
  num_t inputSamples[BUFFER_LEN];
  num_t gain = 1.0;
  bool start, done;
   //FILE *fd_in;
   //unsigned short input[SAMPLES_SIZE];
   string line;
   int count = 0;
   //unsigned short temp;
  
  num_t read_bram(sc_dt::uint64 addr, unsigned char type); 
  void write_bram(sc_dt::uint64 addr, num_t val, unsigned char type); 

  num_t read_hard(sc_dt::uint64 addr);
  void write_hard(sc_dt::uint64 addr, int val);

};

#endif // SOFT_HPP_
