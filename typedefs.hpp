#ifndef TYPEDEFS_HPP
#define TYPEDEFS_HPP

#include <systemc>
#include <tlm>
#define SC_INCLUDE_FX

typedef float num_t; // lebegőpontos szám használata





typedef tlm::tlm_base_protocol_types::tlm_payload_type pl_t;
typedef tlm::tlm_base_protocol_types::tlm_phase_type ph_t;

#define ADDR_CMD 0x18
#define ADDR_STATUS 0x1c
#define ADDR_RES 0x20
#define ADDR_STATUS_INIT 0x38
#define ADDR_OFFSET 0x30

#define VP_ADDR_BRAM_BASE 0x01000000
#define VP_ADDR_BRAM_L 0x00000000
#define VP_ADDR_BRAM_H 0x0003BFFF

#define VP_ADDR_HARD_BASE 0x02000000
#define VP_ADDR_HARD_L 0x02000000
#define VP_ADDR_HARD_H 0x020000FF

#define VP_ADDR_SAMPLES_BASE 0x00000000
#define VP_ADDR_SAMPLES_L 0x00000000
#define VP_ADDR_SAMPLES_H 0x0000FFFF

#define VP_ADDR_COEFFS_BASE 0x00010000
#define VP_ADDR_COEFFS_L  0x00010000
#define VP_ADDR_COEFFS_H  0x0001FFFF

#define VP_ADDR_RESULT_BASE 0x00020000
#define VP_ADDR_RESULT_L 0x00020000
#define VP_ADDR_RESULT_H 0x0003BFFF

#define SAMPLES_SIZE 1024
#define OUTPUT_SIZE 1024
#define COEFFS_SIZE 1024
#define BLOCK_SIZE 1024
#define BUFFER_LEN (SAMPLES_SIZE -1 + COEFFS_SIZE)
#define MEM_RESERVED 100000
constexpr int BUFF_SIZE = 8;  // double = 8 bájt
#endif // TYPEDEFS_HPP
