#include "soft.hpp"
#include <fstream>
#include "AudioFile.h"  // dodato
#include <complex>      // dodato
#include <vector>       // dodato
#include <cstring>      // Dodato - For strcmp
#include <chrono>  // Hozzáadjuk a chrono könyvtárat
using namespace std::chrono;  // Biztosítjuk, hogy elérhetők legyenek az időzítési funkciók
using Complex = std::complex<double>;


int thr;
char* g_coeffs;       // Coeffs global



const sc_dt::uint64 BASS = 0;
const sc_dt::uint64 MID = 1;
const sc_dt::uint64 TREBLE = 2;



int resize(int samples);

SC_HAS_PROCESS(Soft);
std::vector<num_t> getCoefficients(const char *coefs);

Soft::Soft(sc_core::sc_module_name name, const std::string& wavFile, const std::string& preset)
    : sc_module(name), wavFile(wavFile), preset(preset)
    {
        // Inicializálás
        std::cout << "WAV fájl: " << wavFile << ", preset: " << preset << std::endl;
  
 
    // FFT számítás
    

    // Az eredmények kiírása
   // for (size_t i = 0; i < fft_result.size(); i++) {
     //   std::cout << "Result[" << i << "] = " << fft_result[i].real()
       //           << " + " << fft_result[i].imag() << "i" << std::endl;
   // }

    SC_THREAD(gen); // A feldolgozást elindítjuk egy szálon
    SC_REPORT_INFO("Soft", "Constructed.");
}

int resize(int samples)
  {
	  int nth_degree[21] = {1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824};

	  for (int i = 0; i < 21; i++)
	  {
		  if (samples < nth_degree[i])
		  {
			  return nth_degree[i];
			  break;
		  }
	  }
	  return samples;
  }

Soft::~Soft()
{
  SC_REPORT_INFO("Soft", "Destroyed.");
}

void Soft::gen()
{
    num_t write_val, read_val;
    int t = 0;
  
    int num = (int) std::ceil(count / BLOCK_SIZE);
  
    std::cout << "Number of blocks: " << num << std::endl;
    std::cout << "Coeffs global in gen(): " << g_coeffs << std::endl;
  
    // A koefficienseket a WAV fájlból kinyert értékekből számítjuk
    std::vector<num_t> v_coeffs; // Dinamikusan a WAV-ból jönnek
    std::vector<Complex> complex_coeffs; // A komplex koefficiensek

    // Feltételezzük, hogy az audio fájl már betöltődött, és a minták a complex_coeffs vektorban vannak
    for (size_t i = 0; i < complex_coeffs.size(); ++i) {
        v_coeffs.push_back(static_cast<num_t>(complex_coeffs[i].real())); // Csak a valós részt használjuk, ha nem kell az imaginárius rész
    }

// Inicializáljuk a koefficiens vektort
        

        // Válasszunk presetet
        if (preset == "BASS") {
            v_coeffs = {4, 3, 2, 0, 0, 0, 0, 0, 0, 0}; // BASS koefficiensek
        }
        else if (preset == "MID") {
            v_coeffs = {0, 0, 0, 4, 4, 4, 4, 0, 0, 0}; // MID koefficiensek
        }
        else if (preset == "TREBLE") {
            v_coeffs = {4, 3, 2, 0, 0, 0, 0, 2, 3, 4}; // TREBLE koefficiensek
        }
        else {
            std::cerr << "Érvénytelen preset!" << std::endl;
            return;
        }


    // Koefficiensek írása a BRAM-ba
    for (size_t i = 0; i < v_coeffs.size(); ++i) 
    {
        write_val = v_coeffs[i];
        write_bram(i, write_val, 1);
        std::cout << "Coeff at index " << i << ": " << v_coeffs[i] << std::endl;
    }

    std::cout << "Coeffs written to BRAM." << std::endl;
  
    // Hardware elindítása
    write_hard(ADDR_CMD, 1);  // start hardware

    while (read_hard(ADDR_STATUS) == 0) // wait for hardware to finish
    {
        wait(10.8, sc_core::SC_NS);
    }

    num_t temp_real;
    num_t temp_imag;

    // Eredmények olvasása és kiírása
    for (int i = 0; i < COEFFS_SIZE; i++) {
        // Valós és imaginárius rész olvasása a BRAM-ból
        temp_real = read_bram(i, 2);
        temp_imag = read_bram(i + COEFFS_SIZE, 2);

        std::cout << "Values from HARD[" << i << "] = " << temp_real 
                  << " , " << temp_imag << "i" << std::endl;

        // Eredmények fájlba írása (ha szükséges)
        // output << temp_real << " , " << temp_imag << "i," << std::endl;
    }

    t++;
}


void Soft::processAudio() {
uint8_t preset_val = 0;
if (preset == "bass") {
        preset_val = BASS;
    } else if (preset == "mid") {
        preset_val = MID;
    } else if (preset == "treble") {
        preset_val = TREBLE;
    } else {
        std::cerr << "Hiba: Ismeretlen preset!" << std::endl;
        return;
    }
    
    auto startMain = high_resolution_clock::now();

    AudioFile<double> audioFile;
    vector<complex<double>> audioFileComplex;
    vector<complex<double>> transform;
    size_t oldSize;
    size_t newSize;
    float sampleRate;

    bool loadedOK = audioFile.load(wavFile);
    if (!loadedOK) {
        std::cerr << "Problem with reading the wav file" << std::endl;
        return;
    }

    oldSize = audioFile.getNumSamplesPerChannel();
    sampleRate = audioFile.getSampleRate();
    newSize = resize(audioFile.getNumSamplesPerChannel());
    audioFile.samples[0].resize(newSize, 0);

    for (int i = 0; i < audioFile.getNumSamplesPerChannel(); i++) {
        audioFileComplex.push_back(audioFile.samples[0][i]);
    }

    auto start = high_resolution_clock::now();
    transform = FFT(audioFileComplex);
    auto stop = high_resolution_clock::now();
    std::cout << "FFT duration: " << duration_cast<milliseconds>(stop - start).count() << " ms" << std::endl;

    

    vector<double> freq = getFreq(audioFile.getNumSamplesPerChannel(), sampleRate);
    vector<double> gauss = generate_gaussian(freq, preset_val);
    vector<complex<double>> filteredTransform = transform;

    for (int i = 0; i < freq.size(); i++) {
        filteredTransform[i] += filteredTransform[i] * gauss[i];
    }

    filteredTransform = normalizePower(filteredTransform, transform);

    start = high_resolution_clock::now();
    vector<complex<double>> ifftOutput = IFFT(filteredTransform);
    stop = high_resolution_clock::now();
    std::cout << "IFFT duration: " << duration_cast<milliseconds>(stop - start).count() << " ms" << std::endl;

    if (ifftOutput.empty()) {
        std::cerr << "Error: ifftOutput is empty!" << std::endl;
        return;
    }

    ifftOutput.resize(oldSize);

    AudioFile<double> outputAudio;
    outputAudio.setNumChannels(2);
    outputAudio.setNumSamplesPerChannel(ifftOutput.size());
    outputAudio.setSampleRate(sampleRate);

    for (int i = 0; i < ifftOutput.size(); i++) {
        for (int channel = 0; channel < outputAudio.getNumChannels(); channel++) {
            outputAudio.samples[channel][i] = ifftOutput[i].real();
        }
    }

    if (!outputAudio.save("output.wav", AudioFileFormat::Wave)) {
        std::cerr << "Error saving the WAV file." << std::endl;
        return;
    }

    auto stopMain = high_resolution_clock::now();
    std::cout << "Main duration: " << duration_cast<milliseconds>(stopMain - startMain).count() << " ms" << std::endl;
}






void Soft::FFT_Function()
{
    num_t temp;
    
    // AudioFile objektum létrehozása és a .wav fájl betöltése
    AudioFile<double> audioFile;
    if (!audioFile.load("guitar_mono.wav")) {
        std::cerr << "Hiba a WAV fájl betöltésekor!" << std::endl;
        return;
    }

    // A bal csatornán található minták száma
    size_t numSamples = audioFile.getNumSamplesPerChannel();
    std::cout << "Minták száma: " << numSamples << std::endl;

    // A SAMPLES változó értékét a teljes fájl mintáinak számához állítjuk
    const int SAMPLES = static_cast<int>(numSamples); // Az összes minta feldolgozása

    // A komplex koefficiensek előkészítése a dinamikusan meghatározott minták számával
    std::vector<Complex> complex_coeffs(SAMPLES);

    // Feltételezzük, hogy a WAV fájl sztereó, és a bal csatornát (vagy mono fájl esetén az egyetlen csatornát) szeretnénk használni
    const std::vector<double>& leftChannel = audioFile.samples[0]; // 0 a bal csatorna

    std::cout << "Minták olvasása a WAV fájlból..." << std::endl;

    // A komplex koefficiensek feltöltése a WAV fájlból olvasott valós értékekkel
    for (int i = 0; i < SAMPLES; i++) {
        double temp = leftChannel[i];  // A bal csatornából származó valós érték
        std::cout << "coefsH[" << i << "] = " << temp << std::endl;
        coeffsH[i] = temp;
        complex_coeffs[i] = Complex(temp, 0.0);  // Az imaginárius rész 0
    }

    std::cout << "\nFFT végrehajtása " << SAMPLES << " mintán..." << std::endl;

    // FFT számítása (ez igazodhat a saját FFT függvényedhez)
    std::vector<Complex> fft_result = FFT(complex_coeffs);

    // Az eredmények írása a BRAM-ba
    for (int i = 0; i < SAMPLES; i++) {
        result[i] = static_cast<num_t>(fft_result[i].real());    // Valós rész
        double imag_part = static_cast<num_t>(fft_result[i].imag()); // Imaginárius rész

        std::cout << "Result[" << i << "] = " << result[i]
                  << " , " << imag_part << "i" << std::endl;

        // A valós és imaginárius részek írása a BRAM-ba
        write_bram(i, result[i], 2);
        write_bram(i + SAMPLES, imag_part, 2); // Imaginárius rész külön helyre
    }

    std::cout << "\nFFT számítás kész " << SAMPLES << " mintán! A valós és imaginárius részek írása BRAM-ba.." << std::endl;

    done = true; // A művelet befejezésének jelezése
}


vector<Complex> Soft::FFT(const vector<Complex>& in)
{
  //vector<Complex> x = in; // Copy input data
    vector<Complex> x(in);

	// DFT
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

vector<Complex> Soft::IFFT(const vector<Complex>& in) {
    vector<Complex> x1(in);
    size_t size = in.size();
    
    // Komplex konjugálás
    for (size_t i = 0; i < size; i++) {
        x1[i] = std::conj(x1[i]);
    }

    // FFT meghívása
    x1 = FFT(x1);  // Most már elérhető az FFT függvény

    // Komplex konjugálás ismétlése
    for (size_t i = 0; i < size; i++) {
        x1[i] = std::conj(x1[i]);
    }

    // Normalizálás
    for (size_t i = 0; i < size; i++) {
        x1[i] /= size;
    }

    return x1;
}







num_t Soft::read_bram(sc_dt::uint64 addr, unsigned char type)
{
  pl_t pl; //payload
 unsigned char buf[BUFF_SIZE]; //buffer that converts num_t into unsigned_char

  sc_dt::uint64 taddr= (addr*BUFF_SIZE) | VP_ADDR_BRAM_BASE;

  switch(type)
    {
    case 0:
      taddr |= VP_ADDR_SAMPLES_BASE;
      break;
    case 1:
      taddr |= VP_ADDR_COEFFS_BASE;
      break;
    case 2:
      taddr |= VP_ADDR_RESULT_BASE;
      break;
    default:
      break;
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
  sc_dt::uint64 taddr = (addr*BUFF_SIZE) | VP_ADDR_BRAM_BASE;

  switch(type)
    {
    case 0:
      taddr |= VP_ADDR_SAMPLES_BASE;
      break;
    case 1:
      taddr |= VP_ADDR_COEFFS_BASE;
      break;
    case 2:
      taddr |= VP_ADDR_RESULT_BASE;
      break;
    default:
      break;
    }
  
  pl.set_address(taddr);
  pl.set_data_length(BUFF_SIZE);
  pl.set_data_ptr(buf);
  pl.set_command( tlm::TLM_WRITE_COMMAND );
  pl.set_response_status ( tlm::TLM_INCOMPLETE_RESPONSE );
  interconnect_socket->b_transport(pl,offset);
  wait (10.8,sc_core::SC_NS); 
}

num_t Soft::read_hard(sc_dt::uint64 addr)
{
  pl_t pl;
  unsigned char buf[BUFF_SIZE]; //buffer that converts num_t into unsigned_char
  pl.set_address(addr | VP_ADDR_HARD_BASE);
  pl.set_data_length(BUFF_SIZE);
  pl.set_data_ptr(buf);
  pl.set_command( tlm::TLM_READ_COMMAND );
  pl.set_response_status ( tlm::TLM_INCOMPLETE_RESPONSE );
  sc_core::sc_time offset = sc_core::SC_ZERO_TIME;
  interconnect_socket->b_transport(pl,offset); //trasnfer

  return to_fixed(buf);
}

void Soft::write_hard(sc_dt::uint64 addr, int val)
{
  pl_t pl;
  unsigned char buf[BUFF_SIZE];
  to_uchar (buf, val); //convert val to unsigned char
  pl.set_address(addr | VP_ADDR_HARD_BASE);
  pl.set_data_length(BUFF_SIZE);
  pl.set_data_ptr(buf);
  pl.set_command( tlm::TLM_WRITE_COMMAND );
  pl.set_response_status ( tlm::TLM_INCOMPLETE_RESPONSE );
  interconnect_socket->b_transport(pl,offset); //transfer
}

/* std::vector<num_t> getCoefficients(const char *coefs) {
    std::vector<num_t> coeffs(COEFFS_SIZE, 0); // Initialize with zeros

    if (std::strcmp(coefs, "BASS") == 0) {
        coeffs = {4, 3, 2, 0, 0, 0, 0, 0, 0, 0};
    } else if (std::strcmp(coefs, "TREBLE") == 0) {
        coeffs = {0, 0, 0, 4, 4, 4, 4, 0, 0, 0};
    }

    return coeffs;
}
*/




vector<double> Soft::getFreq(int n, int sampleRate)
{
  vector<double> freq(n);
	double d = 1.0 / sampleRate;
	double val = 1.0 / (n * d);
	int N = (n - 1) / 2 + 1;

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



/*param:
		gauss - array with all values set to zero;
		gauss_freq - array with frequences where we generate our gaussian function
		gauss_amp - amplitude of a gauss function on a given freq
*/


vector<sc_dt::uint64> gauss_freq{31, 62, 125, 250, 500, 1000, 2000, 4000, 8000, 16000};



vector<vector<int>> gauss_amp
{
	{4, 3, 2, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 4, 4, 4, 4, 0, 0, 0},
	{4, 3, 2, 0, 0, 0, 0, 2, 3, 4}
  };

  vector<double> Soft::generate_gaussian(vector<double> freq, uint8_t preset)
{
	vector<double> gauss(freq.size());
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



/*
Peak normalization of a fft Array.

param:
		fftArr - the vector that needs to be normalzied.
		returns vector that is normalized

*/

  vector<Complex> Soft::normalizePeak(vector<Complex> fftArr, uint8_t preset)
{
	vector<Complex> out(fftArr.size());

	int normalizeValue = *max_element(gauss_amp[preset].begin(), gauss_amp[preset].end()) + 1;

	for (int i = 0; i < fftArr.size(); i++)
	{
		out[i].real(fftArr[i].real() / normalizeValue);
	}
	return out;
}

/*
Power normalization:

param:
		modifiedFftArr - the fft vector that we modified by the gauss array
		fftArr - the unmodofied vector
		returns vector that is normalized

*/

vector<Complex> Soft::normalizePower(vector<Complex> modifiedFftArr, vector<Complex> fftArr)
{
	double power1 = 0;
	double power2 = 0;
	double multiplier;

	// calculate the total power of the unmodified signal
	for (int i = 0; i <= fftArr.size(); i++)
	{
		power1 += abs(fftArr[i].real());
	}
	// calculate the total power of the modified signal
	for (int i = 0; i <= fftArr.size(); i++)
	{
		power2 += abs(modifiedFftArr[i].real());
	}

	// the square of the ratios of the 2 powers gives me a multiplier which we will use to normalize our signal
	multiplier = pow(power1 / power2, 2);
	vector<Complex> out(fftArr.size());

	// multiply every element of the array
	for (int i = 0; i < fftArr.size(); i++)
	{
		out[i] = modifiedFftArr[i] * multiplier;
	}

	return out;
}



