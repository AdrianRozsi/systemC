#include "soft.hpp"
#include <fstream>
#include "AudioFile.h"
#include <complex>
#include <vector>
#include <cstring>
#include <cmath>
#include "typedefs.hpp"
#include <chrono>
#include "utils.hpp"

using namespace std;
using namespace std::chrono;


int thr;
char* g_coeffs;
using Complex = std::complex<double>;
using namespace sc_core;




/*param:
		gauss - array with all values set to zero;
		gauss_freq - array with frequences where we generate our gaussian function
		gauss_amp - amplitude of a gauss function on a given freq
*/

std::vector<double> Soft::generate_gaussian(const std::vector<double>& freq, uint8_t preset)
{
   
  
    	std::vector<double> gauss(freq.size());
	for (int i = 0; i < freq.size(); i++)
	{
		gauss[i] = 0;
	}

	for (int i = 0; i < gauss_freq.size(); i++)
	{
		for (int j = 0; j < freq.size(); j++)
			gauss[j] += gauss_amp[preset][i] * exp(-(pow(freq[j] - gauss_freq[i], 2) / (2 * pow(gauss_freq[i] / 3, 2))));
	}

	return gauss;
}

std::vector<double>Soft::getFreq(int n, int sampleRate)
{
    std::vector<double> freq(n);
     double d   = 1.0 / sampleRate;          // időkülönbség (s)
     double val = 1.0 / (n * d);             // Fs / n  (Hz per bin)
    int N = (n - 1) / 2 + 1;         // utolsó pozitív frekvencia indexe

    for (int i = 0; i <= N; i++)
	{
		freq[i] = i * val;
	}

	int backWards = n / 2;
	for (int i = N + 1; i <= n; i++)
	{
		freq[i] = -backWards * val;
		backWards--;
	}
	return freq;
}



SC_HAS_PROCESS(Soft);
std::vector<num_t> getCoefficients(const char *coefs);

Soft::Soft(sc_core::sc_module_name name, char *coefs, char *samples, uint8_t preset_val)
  : sc_module(name), offset(sc_core::SC_ZERO_TIME), preset(preset_val), samplesPath(samples)
{
 



    SC_THREAD(gen);
    SC_REPORT_INFO("Soft", "Constructed");
}


Soft::~Soft()
{
    

    SC_REPORT_INFO("Soft", "Destroyed.");
}

void Soft::gen() {
    bool loadedOK = audioFile.load(samplesPath);
    if (!loadedOK) {
        std::cout << "Probléma a wav fájl beolvasásával" << std::endl;
        sc_core::sc_stop();
        return;
    }

    oldSize = audioFile.getNumSamplesPerChannel(); // Eredeti méret mentése
    sampleRate = audioFile.getSampleRate();

    // Méret átméretezése a legközelebbi 2 hatványra, 0-ás kitöltéssel
    newSize = resize(audioFile.getNumSamplesPerChannel());
    audioFile.samples[0].resize(newSize);
    for (int i = oldSize; i < newSize; i++)
        audioFile.samples[0][i] = 0;

    // Konvertálás double-ból komplex double-ba
    // A képzetes részeket 0-ra állítjuk
    for (int i = 0; i < audioFile.getNumSamplesPerChannel(); i++)
        audioFileComplex.push_back(audioFile.samples[0][i]);

    const int block_size = SAMPLES_SIZE;
    int numBlocks = (newSize + block_size - 1) / block_size;

    transform = FFT_Soft(audioFileComplex); // FFT végrehajtása

    vector<double> freq(audioFile.getNumSamplesPerChannel());
    vector<complex<double>> filteredTransform(transform);
    vector<double> filteredData(oldSize);

    // A frekvenciák listája
    freq = getFreq(audioFile.getNumSamplesPerChannel(), audioFile.getSampleRate());

    // Gauss függvény generálása
    gauss = generate_gaussian(freq, preset);

    // Gauss szűrés alkalmazása
    for (int i = 0; i < freq.size(); i++)
        filteredTransform[i] += filteredTransform[i] * gauss[i];

    // A jel normalizálása, hogy ne legyen túl nagy erősítés
    filteredTransform = normalizePower(filteredTransform, transform);

    // Kiírás az első 100 FFT és Gauss szűrt eredményről fájlokba
    std::ofstream fftFile("fft_result.txt");
    for (size_t i = 0; i < std::min(100UL, transform.size()); ++i) {
        fftFile << "Index: " << i << " | Real: " << transform[i].real() << " | Imag: " << transform[i].imag() << std::endl;
    }
    fftFile.close();

    std::ofstream gaussFile("gauss_filtered_result.txt");
    for (size_t i = 0; i < std::min(100UL, filteredTransform.size()); ++i) {
        gaussFile << "Index: " << i << " | Real: " << filteredTransform[i].real() << " | Imag: " << filteredTransform[i].imag() << std::endl;
    }
    gaussFile.close();

    // Blokkok feldolgozása és átfedés-adás
    vector<Complex> outputData(newSize, 0.0);  // Az összesített adat
    vector<Complex> overlapBuffer(block_size, 0.0); // Átfedés-tároló puffer

    for (int b = 0; b < numBlocks; ++b) {
        int start_idx = b * block_size;
        size_t end_idx = std::min(static_cast<size_t>(start_idx + block_size), newSize);

        vector<Complex> block(block_size);
        for (int i = 0; i < block_size; ++i) {
            int idx = start_idx + i;
            block[i] = (idx < newSize) ? filteredTransform[idx] : Complex(0.0, 0.0);
        }

        // Blokk küldése BRAM-ba
        for (int i = 0; i < block_size; ++i) {
            write_bram(i * 2, block[i].real(), 1);
            write_bram(i * 2 + 1, block[i].imag(), 1);
        }

        // Várakozás a hardver válaszára
        wait(10.8, sc_core::SC_NS);
        write_hard(ADDR_OFFSET, 0);
        write_hard(ADDR_CMD, 1);

        int attempts = 0;
        while (!read_hard()) {
            wait(100, sc_core::SC_NS);
            if (++attempts > 10000) {
                std::cerr << "[SOFT] Időtúllépés a HARD válaszára várva!" << std::endl;
                sc_stop();
                return;
            }
        }

        // A visszakapott blokkok írása az outputData-ba
        for (int i = 0; i < block_size; ++i) {
            int global_idx = start_idx + i;
            if (global_idx >= newSize) break;

            num_t re = read_bram(i * 2, 2);
            outputData[global_idx] += re; // Átfedés-adás
        }
         outputData.resize(oldSize);
        // Az utolsó blokk átfedésének kezelése
        if (b < numBlocks - 1) {
            // Az overlap-tároló pufferbe beírjuk az aktuális blokk végét
            for (int i = 0; i < overlapBuffer.size(); ++i) {
                overlapBuffer[i] = block[i];
            }
        }
    }

    // Végső adatokat visszaírjuk az audio fájlba
    outputAudio.setNumChannels(2);  // Két csatorna
    outputAudio.setNumSamplesPerChannel(outputData.size());  // Az eredeti méret
    outputAudio.setSampleRate(sampleRate);  // A mintavételi frekvencia

    for (int i = 0; i < oldSize; ++i) {
        for (int channel = 0; channel < outputAudio.getNumChannels(); ++channel) {
            outputAudio.samples[channel][i] = outputData[i].real();  // Csak a valós rész
        }
    }
    

    outputAudio.save("output.wav", AudioFileFormat::Wave);  // Elmentés WAV fájlként
    std::cout << "[SOFT] output.wav fájl generálva." << std::endl;
    sc_core::sc_stop();
}






int Soft::resize(int samples)
{
    int nth_degree[21] = {
        1024, 2048, 4096, 8192, 16384, 32768, 65536,
        131072, 262144, 524288, 1048576, 2097152, 4194304,
        8388608, 16777216, 33554432, 67108864, 134217728,
        268435456, 536870912, 1073741824
    };

    for (int i = 0; i < 21; i++) {
        if (samples < nth_degree[i]){
            return nth_degree[i];
            break;
        }
    }
    //return std::min(samples, 1073741824);
    return samples;
}










std::vector<Complex> Soft::FFT_Soft(const std::vector<std::complex<double>>& in)
{
    std::vector<Complex> x(in);

    //DFT
    unsigned int N = x.size(), k = N, n;
    double thetaT = 3.14159265358979323846264338328L / N;
    Complex phiT = Complex(cos(thetaT), -sin(thetaT)), T;
    while (k > 1)
    {
        n = k;
        k >>= 1;
        phiT = phiT * phiT;
        T = 1.0L;
        for (unsigned int l = 0; l < k; l++)
        {
            for (unsigned int a = l; a < N; a += n)
            {
                unsigned int b = a + k;
                Complex t = x[a] - x[b];
                x[a] += x[b];
                x[b] = t * T;
            }
            T *= phiT;
        }
    }
    // Decimate
    unsigned int m = (unsigned int)log2(N);
    for (unsigned int a = 0; a < N; a++)
    {
        unsigned int b = a;
        // Reverse bits
        b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
        b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
        b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
        b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
        b = ((b >> 16) | (b << 16)) >> (32 - m);
        if (b > a)
		{
			Complex t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
    }

    return x;
}


num_t Soft::read_bram(sc_dt::uint64 addr, unsigned char type)
{
  pl_t pl; //payload
  unsigned char buf[BUFF_SIZE]; //buffer that converts num_t into unsigned_char

uint64_t taddr = 0;
switch(type)
{
  case 0: taddr = VP_ADDR_SAMPLES_BASE + addr * BUFF_SIZE; break;
  case 1: taddr = VP_ADDR_COEFFS_BASE + addr * BUFF_SIZE; break;
  case 2: taddr = VP_ADDR_RESULT_BASE + addr * BUFF_SIZE; break;
  default: SC_REPORT_ERROR("Soft::read_bram", "Invalid BRAM type"); break;
}
  pl.set_address(taddr);
  pl.set_data_length(BUFF_SIZE); //length of buffer
  pl.set_data_ptr(buf);
  pl.set_command( tlm::TLM_READ_COMMAND );
  pl.set_response_status ( tlm::TLM_INCOMPLETE_RESPONSE );
  interconnect_socket->b_transport(pl,offset);
  wait (10.8,sc_core::SC_NS);
  return to_fixed(buf); //unsigned_char to buff;
}




void Soft::write_bram(sc_dt::uint64 addr, num_t val, unsigned char type)
{

  
  pl_t pl;
  unsigned char buf[BUFF_SIZE]; //buffer that converts num_t into unsigned char for ic
  to_uchar(buf,val);
  //sc_dt::uint64 taddr = (addr*BUFF_SIZE) | VP_ADDR_BRAM_BASE;

  uint64_t taddr = 0;
switch(type)
{
  case 0: taddr = VP_ADDR_SAMPLES_BASE + addr * BUFF_SIZE; break;
  case 1: taddr = VP_ADDR_COEFFS_BASE + addr * BUFF_SIZE; break;
  case 2: taddr = VP_ADDR_RESULT_BASE + addr * BUFF_SIZE; break;
  default: SC_REPORT_ERROR("Soft::write_bram", "Invalid BRAM type"); break;
}
 // std::cout << "[SOFT DEBUG] Write COEFFS addr=0x"
   //       << std::hex << taddr << std::dec
     //     << " val = " << val << std::endl;
  pl.set_address(taddr);
  pl.set_data_length(BUFF_SIZE);
  pl.set_data_ptr(buf);
  pl.set_command( tlm::TLM_WRITE_COMMAND );
  pl.set_response_status ( tlm::TLM_INCOMPLETE_RESPONSE );
  interconnect_socket->b_transport(pl,offset);
  wait (10.8,sc_core::SC_NS); 
}




bool Soft::read_hard()
{
    pl_t pl;
    unsigned char buf[BUFF_SIZE];
    pl.set_address(VP_ADDR_HARD_BASE + ADDR_STATUS);
    pl.set_data_length(BUFF_SIZE);
    pl.set_data_ptr(buf);
    pl.set_command(tlm::TLM_READ_COMMAND);
    pl.set_response_status(tlm::TLM_INCOMPLETE_RESPONSE);
    sc_core::sc_time offset = sc_core::SC_ZERO_TIME;
    interconnect_socket->b_transport(pl, offset);
   num_t status_val = to_fixed(buf);
//std::cout << "[SOFT] read_hard status = " << status_val << std::endl;
return (status_val != 0.0);
}





/*
num_t Soft::read_hard(sc_dt::uint64 addr)
{
  pl_t pl;
  unsigned char buf[BUFF_SIZE]; //buffer that converts num_t into unsigned_char
  pl.set_address(VP_ADDR_HARD_BASE + addr);
  pl.set_data_length(BUFF_SIZE);
  pl.set_data_ptr(buf);
  pl.set_command( tlm::TLM_READ_COMMAND );
  pl.set_response_status ( tlm::TLM_INCOMPLETE_RESPONSE );
  sc_core::sc_time offset = sc_core::SC_ZERO_TIME;
  interconnect_socket->b_transport(pl,offset); //trasnfer

  return to_fixed(buf); 

}

*/

void Soft::write_hard(sc_dt::uint64 addr, num_t val)
{
  pl_t pl;
  unsigned char buf[BUFF_SIZE];
  to_uchar (buf, val); //convert val to unsigned char
  pl.set_address(VP_ADDR_HARD_BASE + addr);
  pl.set_data_length(BUFF_SIZE);
  pl.set_data_ptr(buf);
  pl.set_command( tlm::TLM_WRITE_COMMAND );
  pl.set_response_status ( tlm::TLM_INCOMPLETE_RESPONSE );
  interconnect_socket->b_transport(pl,offset); //transfer
}

std::vector<Complex> Soft::normalizePeak(const std::vector<Complex>& fftArr, uint8_t preset)
{
    std::vector<Complex> out(fftArr.size());
    

    int normalizeValue = *max_element(gauss_amp[preset].begin(), gauss_amp[preset].end()) + 1;

    for (int i = 0; i < fftArr.size(); i++){
       out[i].real(fftArr[i].real() / normalizeValue);
    }
    return out;
}

std::vector<Complex> Soft::normalizePower(const std::vector<Complex>& modifiedFftArr, const std::vector<Complex>& fftArr)
{
    double power1 = 0;
    double power2 = 0;
    double multiplier;

    for (int i = 0; i < fftArr.size(); i++)
    { 
           power1 += std::abs(fftArr[i].real());
    }
    for (int i = 0; i < fftArr.size(); i++)
    {
           power2 += std::abs(modifiedFftArr[i].real());
    }
    //std::cout << "[normalizePower] power1 = " << power1 << ", power2 = " << power2 << std::endl;

    if (power2 == 0)
    {
        std::cerr << "[normalizePower] Warning: modified FFT has zero power!\n";
        return modifiedFftArr;
    }

    multiplier = pow(power1 / power2, 2);

    vector<complex<double>> out(fftArr.size());
    for (int i = 0; i < fftArr.size(); i++)
    {
        out[i] = modifiedFftArr[i] * multiplier;
    }
    return out;
}




// Finding nearest number power of 2
/*int Soft::resize(int samples)
{
    int nth_degree[21] = {
        1024, 2048, 4096, 8192, 16384, 32768, 65536,
        131072, 262144, 524288, 1048576, 2097152, 4194304,
        8388608, 16777216, 33554432, 67108864, 134217728,
        268435456, 536870912, 1073741824
    };

    for (int i = 0; i < 21; i++) {
        if (samples < nth_degree[i]){
            return nth_degree[i];
            break;
        }
    }
    return samples;
}*/