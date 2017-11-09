#define DIM 16
#define SBOX_SIZE 4
#define NUM_SBOXES (DIM/SBOX_SIZE)

uint64_t mds[] = {
    0x00008414,
    0x000042c2,
    0x00002161,
    0x00001c3c,
    0x00004841,
    0x0000242c,
    0x00001216,
    0x0000c1c3,
    0x00002884,
    0x00001442,
    0x0000c221,
    0x0000611c,
    0x00008248,
    0x00004124,
    0x00002c12,
    0x000016c1
};