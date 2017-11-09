#define DIM 16
#define SBOX_SIZE 4
#define NUM_SBOXES (DIM/SBOX_SIZE)

uint64_t mds[] = {
    0x0000488c,
    0x00002446,
    0x00001223,
    0x0000c11d,
    0x0000c488,
    0x00006244,
    0x00003122,
    0x0000dc11,
    0x00008c48,
    0x00004624,
    0x00002312,
    0x00001dc1,
    0x000088c4,
    0x00004462,
    0x00002231,
    0x000011dc
};