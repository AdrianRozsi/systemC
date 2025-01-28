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

int resize(int samples);

SC_HAS_PROCESS(Soft);
std::vector<num_t> getCoefficients(const char *coefs);


Soft::Soft(sc_core::sc_module_name name, char *coefs, char *samples)
  : sc_module(name)
  , offset(sc_core::SC_ZERO_TIME)
{

  coeffsFile.open(coefs);
  samplesFile.open(samples);
  output.open("result.txt");  
  //AudioFile<double> audioFile;  
  samplesFile.open(samples);


  if(!coeffsFile.is_open() || !samplesFile.is_open() )
    SC_REPORT_ERROR("Soft", "Cannot open file.");

  if(samplesFile.is_open())
    {
      while(samplesFile.peek()!=EOF)
	{
	  getline(samplesFile, line);
	  count++;
	}
    }  



  //audioFile.open(samples);
  std::cout<<"arg1: "<<coefs<<std::endl;
  std::cout<<"arg2: "<<samples<<std::endl;

  g_coeffs = coefs;
  std::cout<<"Coeffs global: "<<g_coeffs<<std::endl;


  //Dodato iz main ----
  AudioFile<double> audioFile;              // audiFile (vector)
  vector<complex<double>> audioFileComplex;
  vector<complex<double>> transform;
  size_t oldSize;
  size_t newSize;
  float sampleRate;
  bool loadedOK = audioFile.load(samples); // read file

  if (!loadedOK)
  {
    cout << "problem with reading the wav file";
    cout << "call with -h argument: ./filename -h" << endl;

    
  }
  oldSize = audioFile.getNumSamplesPerChannel(); // save original size
  sampleRate = audioFile.getSampleRate();

  std::cout<<"oldSize: "<<oldSize<<std::endl;
  std::cout<<"sampleRate: "<<sampleRate<<std::endl;

  newSize = resize(audioFile.getNumSamplesPerChannel()); 
  std::cout<<"newSize: "<<newSize<<std::endl;

  audioFile.samples[0].resize(newSize);
  for (int i = oldSize; i < newSize; i++)
    audioFile.samples[0][i] = 0;

  for (int i = 0; i < audioFile.getNumSamplesPerChannel(); i++)
    audioFileComplex.push_back(audioFile.samples[0][i]);
  // Do ovde -----
//////////////////////////////////////////////////////////////////////////////////////////////////////
  output.open("result.txt");
  
  if(!samplesFile.is_open() )
    SC_REPORT_ERROR("Soft", "Cannot open file.");

  if(samplesFile.is_open())
    {
      while(samplesFile.peek()!=EOF)
	    {
	      getline(samplesFile, line);
	      count++;
	    }
      std::cout<<"count: "<<count<<std::endl;
    }
  
  samplesFile.seekg(0, ios::beg);
  thr = 0; 
  SC_THREAD(gen);
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
 
  coeffsFile.close(); 
  samplesFile.close();
 
  SC_REPORT_INFO("Soft", "Destroyed.");
}

void Soft::gen()
{
  num_t write_val, read_val;
  int t = 0;
  
  int num =(int) std::ceil(count / BLOCK_SIZE);
  
  std::cout<<"Number of blocks: "<<num<<std::endl;
  std::cout<<"Coeffs global in gen(): "<<g_coeffs<<std::endl;
  
  // Determine coefficients based on arg1
  std::vector<num_t> v_coeffs(COEFFS_SIZE, 0); // Initialize with zeros
  //std::vector<num_t> v_coeffs(160, 0); // Initialize with zeros
  
  //  v_coeffs = {4, 3, 2, 0, 0, 0, 0, 0, 0, 0};    // Values for BASS
  //  v_coeffs = {0, 0, 0, 4, 4, 4, 4, 0, 0, 0};    // Values for MID
  //  v_coeffs = {4, 3, 2, 0, 0, 0, 0, 2, 3, 4};    // Values for TREBLE
  
  v_coeffs = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};    // HC Values for FFT
  //v_coeffs = {4, 3, 2, 0, 0, 0, 0, 0, 0, 0};    // HC Values for BASS

  /*
  for(int i = 0; i < 160; ++i)
  {
    v_coeffs[i] = i + 2;
  }
   */

  for (size_t i = 0; i < v_coeffs.size(); ++i) 
  {
    write_val = v_coeffs[i];
    write_bram(i, write_val,1);
    std::cout<<"coeff at index : "<< i << ":" <<v_coeffs[i]<<std::endl;
  }
  // ---------------------------------------
    
  std::cout<<"Coeffs written to BRAM. "<<std::endl;
  //wait (21*10.8,sc_core::SC_NS);  
  
  //std::cout<<"Samples written to BRAM " <<std::endl<<std::endl;
   
  write_hard(ADDR_CMD, 1);  //start hardware

  while(read_hard(ADDR_STATUS) == 0) // wait for hardwarw to finish
    {
       wait(10.8, sc_core::SC_NS);
    }

  //std::cout << std::endl << "Hardware finished. " << std::endl<<std::endl;
 
  num_t temp_real;
  num_t temp_imag;

  for (int i = 0; i < COEFFS_SIZE; i++) {
      // Valós rész olvasása BRAM-ból
      temp_real = read_bram(i, 2);
      // Imaginárius rész olvasása BRAM-ból
      temp_imag = read_bram(i + COEFFS_SIZE, 2);

      // Eredmények kiírása konzolra
      std::cout << "Values from HARD[" << i << "] = " << temp_real 
                << " , " << temp_imag << "i" << std::endl;

      // Eredmények kiírása fájlba
      //output << temp_real << " , " << temp_imag << "i," << std::endl;
  }
  t++;
}


void Soft::FFT_Function()
{
 num_t temp;
       const int SAMPLES = 16; // Csak az első 16 minta
    std::vector<Complex> complex_coeffs(SAMPLES);

    // Koeficiensek olvasása BRAM-ból és komplex vektor előkészítése
    for (int i = 0; i < SAMPLES; i++) {
        double temp = read_bram(i, 1); // BRAM olvasás
        std::cout << "coefsH[" << i << "] = " << temp << std::endl;
        coeffsH[i] = temp;
        complex_coeffs[i] = Complex(temp, 0.0); // Valós rész feltöltése, imaginárius 0
    }

    std::cout << "\nPerforming FFT on 16 samples..." << std::endl;

    // FFT kiszámítása csak az első 16 mintán
    std::vector<Complex> fft_result = FFT(complex_coeffs);

    // FFT eredmények BRAM-ba írása (valós és imaginárius részek)
    for (int i = 0; i < SAMPLES; i++) {
        result[i] = static_cast<num_t>(fft_result[i].real());    // Valós rész
        double imag_part = static_cast<num_t>(fft_result[i].imag()); // Imaginárius rész

        // Konzolra kiírás
        std::cout << "Result[" << i << "] = " << result[i]
                  << " , " << imag_part << "i" << std::endl;

        // Valós rész BRAM-ba írása
        write_bram(i, result[i], 2);

        // Imaginárius rész BRAM-ba írása
        write_bram(i + SAMPLES, imag_part, 2); // Imaginárius rész külön helyre
    }

    std::cout << "\nFFT Computation DONE for 16 samples! Writing real and imaginary parts to BRAM.." << std::endl;

    done = true; // Jelzés a művelet befejezéséről
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


const sc_dt::uint64 BASS = 0;
const sc_dt::uint64 MID = 1;
const sc_dt::uint64 TREBLE = 2;

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



