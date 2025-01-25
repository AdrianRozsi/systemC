#include "soft.hpp"
#include <fstream>
#include "AudioFile.h"  // dodato
#include <complex>      // dodato
#include <vector>       // dodato
#include <cstring>      // Dodato - For strcmp
int thr;
char* g_coeffs;       // Coeffs global

int resize(int samples);

SC_HAS_PROCESS(Soft);
std::vector<num_t> getCoefficients(const char *coefs);


Soft::Soft(sc_core::sc_module_name name, char *coefs, char *samples)
  : sc_module(name)
  , offset(sc_core::SC_ZERO_TIME)
{
  //AudioFile<double> audioFile;  
  samplesFile.open(samples);
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
  //coeffsFile.close(); -Nem kell
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
 
  num_t temp;
  //double temp_double = temp; // Implicit conversion
  //double temp_double = static_cast<double>(temp);

  // ---- !!!!!! ----
  // Radi i ovo kako treba samo treba popraviti num_t da uzme negativne vr.
  for(int i = 0; i< COEFFS_SIZE; i++)
    {
        num_t temp = read_bram(i, 2); // Olvasás a BRAM-ból
        std::cout << "Values from HARD[" << i << "] = " << temp << std::endl;
        output << temp << "," << std::endl; // Fájlba írás
    }

    std::cout << "FFT eredmények fájlba mentve: fft_result.txt" << std::endl;
  
  t++;

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
