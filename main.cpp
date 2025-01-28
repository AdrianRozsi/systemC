
//g++ -DSC_INCLUDE_FX -w -g -I. -I/sysc/systemc-2.3.2/include -L. -L/sysc/systemc-2.3.2/lib-linux64 -o main *.cpp -lsystemc -lm


#include <systemc>
#include "vp.hpp"
#include "hard.hpp"
#include "soft.hpp"
#include <iostream>
#include <fstream>
#include "AudioFile.h"
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

int sc_main(int argc, char* argv[]) {
    auto startMain = high_resolution_clock::now();
    uint8_t preset;
    string arg1 = argv[1]; // get input file
    string arg2;
    if(argc!=3){
  	std::cout<<"Start with: "<<std::endl;
  	std::cout << "./filename 'wavFileNam' 'presets' " << std::endl;
    std::cout << "presets can be: bass, mid, trebble" << std::endl;
    std::cout << "example: ./program someFile.wav bass" << std::endl;
    return 0;
  }

    if(argc > 2) arg2 = argv[2];

    

    if (arg2 == "bass")
        preset = BASS;
    else if (arg2 == "mid")
        preset = MID;
    else if (arg2 == "treble")
        preset = TREBLE;
    else
    {
        cout << "Argument for presets incorrect" << endl;
        cout << "call with -h argument: ./filename -h" << endl;
        return 1;
    }
    AudioFile<double> audioFile;              // audiFile (vector)
    vector<complex<double>> audioFileComplex; // audiofile converted to complex
                                              // so it can be used in in fft function
    vector<complex<double>> transform;
    size_t oldSize;
    size_t newSize;
    float sampleRate;

    bool loadedOK = audioFile.load(arg1); // read file

if (!loadedOK)
    {
        cout << "problem with reading the wav file";
        cout << "call with -h argument: ./filename -h" << endl;

        return 1;
    }
    oldSize = audioFile.getNumSamplesPerChannel(); // save original size
    sampleRate = audioFile.getSampleRate();

    // resize to the closest 2 exponential, fill with 0's
    int resize(int samples);
    newSize = resize(audioFile.getNumSamplesPerChannel());
    audioFile.samples[0].resize(newSize);
    for (int i = oldSize; i < newSize; i++)
        audioFile.samples[0][i] = 0;

   // convert double to complex double
    // fill imaginary numbers with 0
    for (int i = 0; i < audioFile.getNumSamplesPerChannel(); i++)
        audioFileComplex.push_back(audioFile.samples[0][i]);


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    char* coeffs;
    char* samples;
    Soft soft("soft_instance", coeffs, samples);

    auto start = high_resolution_clock::now();
    transform = soft.FFT(audioFileComplex); // do the fft
    auto stop = high_resolution_clock::now();

    auto fftDuration = duration_cast<milliseconds>(stop - start);
 
    // To get the value of duration use the count()
    // member function on the duration object
    cout <<"fft duration in milliseconds: "<< fftDuration.count() << endl;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


vector<double> freq(audioFile.getNumSamplesPerChannel());
    vector<double> gauss;
    vector<complex<double>> filteredTransform(transform);
    vector<complex<double>> ifftOutput;
    vector<double> filteredData(oldSize);

    // get  array with the list of frequencies
    freq = soft.getFreq(audioFile.getNumSamplesPerChannel(), audioFile.getSampleRate());

    // generate the gaussian function
    gauss = soft.generate_gaussian(freq, preset);

    // apply the gaussian function to the signal
    for (int i = 0; i < freq.size(); i++)
        filteredTransform[i] += filteredTransform[i] * gauss[i];

    // normalize the signal so we dont have to much gain
    filteredTransform = soft.normalizePower(filteredTransform, transform);

    AudioFile<double> outputAudio;

    start = high_resolution_clock::now();
    // do the ifft
    ifftOutput = soft.IFFT(filteredTransform);
    stop = high_resolution_clock::now();


    auto ifftDuration = duration_cast<milliseconds>(stop - start);
    cout <<"ifft duration in milliseconds: "<< ifftDuration.count() << endl;

    // revert to the original size
    ifftOutput.resize(oldSize);

    // set to stereo
    outputAudio.setNumChannels(2);

    outputAudio.setNumSamplesPerChannel(ifftOutput.size());

    // set the sampleRate to the sample rate we read in the beginning
    outputAudio.setSampleRate(sampleRate);

    // copy the contents of the ifftouput array to the outputAudio object which can than be convert to a wav file
    for (int i = 0; i < ifftOutput.size(); i++)
    {
        for (int channel = 0; channel < outputAudio.getNumChannels(); channel++)
        {
            outputAudio.samples[channel][i] = ifftOutput[i].real();
        }
    }

  outputAudio.save("output.wav", AudioFileFormat::Wave);

    auto stopMain = high_resolution_clock::now();

    auto mainDuration = duration_cast<milliseconds>(stopMain - startMain);

    cout <<"main duration in milliseconds: "<< mainDuration.count() << endl;

    // Open a text file in write mode
    ofstream outfile(arg2 + ".txt");

    // Check if the file stream is open
    if (!outfile.is_open()) {
        cerr << "Failed to open the file." << endl;
        return 1;
    }

    // Iterate over the vector and write each element to the file
    for (int i = 0; i < ifftOutput.size(); i++){
        outfile << ifftOutput[i].real() << endl;
    }

    // Close the file
    outfile.close();

    cout << "Data written to file successfully." << endl;


  Vp vp("VP",argv[1],argv[2]);
  sc_start(10, sc_core::SC_MS);
  //std::printf("Throughput is: %d\n", thr*100);
  return 0;
}