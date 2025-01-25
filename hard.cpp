#include "hard.hpp"
#include <vector>
#include <complex>
#include <cmath>
#include <iostream>


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
    std::vector<std::complex<double>> complex_coeffs(COEFFS_SIZE); // Complex array for FFT

    std::cout << "\nStarting Hard...\n" << std::endl;

    // Read coefficients from BRAM and store in the real part of the complex array
    for (int i = 0; i < COEFFS_SIZE; i++) {
        temp = read_bram(i, 1); // Read from BRAM
        std::cout << "coefsH[" << i << "] = " << temp << std::endl;
        coeffsH[i] = temp;
        complex_coeffs[i] = std::complex<double>(temp, 0.0); // Initialize imaginary part to 0
    }

    std::cout << "\nStarting FFT..." << std::endl;

    // Perform FFT
    std::vector<std::complex<double>> fft_result = FFT(complex_coeffs);

    // Display the FFT results and write them back to BRAM
    //std::cout << "Write FFT results to BRAM..." << std::endl;

    // Imamo problem u vracanju negativnih brojeva, koristimo test vector dok ne resimo problem
    std::vector<num_t> test(COEFFS_SIZE, 0);
    test = {-2, 3, -4, -5, 6, 7, 8, 9, 10, 11, -12, 13, 14, 15, 16, -17};

    for (int i = 0; i < COEFFS_SIZE; i++) {
        result[i] = static_cast<num_t>(fft_result[i].real());
        std::cout << "Result[" << i << "] = " << test[i] << std::endl; //treba result[i]
        write_bram(i,test[i], 2);   //treba result[i]
    }
    std::cout << "\nFFT Computation DONE! \nWriting to BRAM.." << std::endl;
    done = true; // Set the done flag
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