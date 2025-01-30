#include "hard.hpp"
#include <vector>
#include <complex>
#include <cmath>
#include <iostream>
#include "AudioFile.h" 

using Complex = std::complex<double>;
extern int thr;

Hard::Hard(sc_core::sc_module_name name)
  : sc_module(name)
  , start(0)
  , done(0)
  , offset(sc_core::SC_ZERO_TIME)
{
  soft_socket.register_b_transport(this, &Hard::b_transport); 
 
  SC_REPORT_INFO("Hard", "Constructed.");
}

Hard::~Hard()
{
  SC_REPORT_INFO("Hard", "Destroyed.");
}

void Hard::FFT_Function()
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


vector<Complex> Hard::FFT(const vector<Complex>& in)
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

void Hard::b_transport(pl_t &pl, sc_core::sc_time &offset)
{
  tlm::tlm_command cmd = pl.get_command();
  sc_dt::uint64 addr = pl.get_address();
  unsigned int len = pl.get_data_length();
  unsigned char *buf = pl.get_data_ptr();
  pl.set_response_status( tlm::TLM_OK_RESPONSE );

  if(len != BUFF_SIZE)
    {
      pl.set_response_status(tlm::TLM_BURST_ERROR_RESPONSE);
    }

  if(cmd == tlm::TLM_WRITE_COMMAND)
    {
      if(addr == ADDR_CMD)
	{
	  // start = to_fixed(buf);
	  start = true;
	  done = !start;
	  std::cout<<"Call FFT Function()"<<std::endl;
	  FFT_Function();
	}
      else
	{
	  pl.set_response_status(tlm::TLM_ADDRESS_ERROR_RESPONSE);
	}
      
    }
  else if(cmd == tlm::TLM_READ_COMMAND)
    {
      if(addr == ADDR_STATUS)
	{
	  to_uchar(buf, done);
	}
      
      else if(addr >= VP_ADDR_SAMPLES_BASE)
	{
	  to_uchar(buf, result[addr - VP_ADDR_SAMPLES_BASE]);
	  }
      else{
	
	pl.set_response_status(tlm::TLM_ADDRESS_ERROR_RESPONSE);
	
      }

    }
  else{
    
    pl.set_response_status(tlm::TLM_COMMAND_ERROR_RESPONSE);
    
  }
    
   offset += sc_core::sc_time(10, sc_core::SC_NS);
}
    
num_t Hard::read_bram(int addr, unsigned char type)
{
  pl_t pl;
  unsigned char buf[BUFF_SIZE];
  pl.set_address(addr*BUFF_SIZE);
  pl.set_data_length(BUFF_SIZE);
  pl.set_data_ptr(buf);
  pl.set_command( tlm::TLM_READ_COMMAND );
  pl.set_response_status ( tlm::TLM_INCOMPLETE_RESPONSE );
  sc_core::sc_time offset = sc_core::SC_ZERO_TIME;

  switch(type)
    {
    case 0:
      bram_samples_socket->b_transport(pl, offset);
      break;
    case 1:
      bram_coeffs_socket->b_transport(pl, offset);
      break;
    case 2:
      bram_result_socket->b_transport(pl, offset);
      break;
    default:
      break;
    }
  wait (10.8,sc_core::SC_NS);
  return to_fixed(buf);
}

void Hard::write_bram(int addr, num_t val, unsigned char type)
{
  pl_t pl;
  unsigned char buf[BUFF_SIZE];
  to_uchar(buf,val);
  pl.set_address(addr*BUFF_SIZE);
  pl.set_data_length(BUFF_SIZE);
  pl.set_data_ptr(buf);
  pl.set_command( tlm::TLM_WRITE_COMMAND );
  pl.set_response_status ( tlm::TLM_INCOMPLETE_RESPONSE );

  switch(type)
    {
    case 0:
      bram_samples_socket->b_transport(pl, offset);
      break;
    case 1:
      bram_coeffs_socket->b_transport(pl, offset);
      break;
    case 2:
      bram_result_socket->b_transport(pl, offset);
      break;
    default:
      break;
      
    }
  wait (10.8,sc_core::SC_NS);    
  
}