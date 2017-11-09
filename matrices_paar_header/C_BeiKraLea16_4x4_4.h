#define DIM 16
#define SBOX_SIZE 4
#define NUM_SBOXES (DIM/SBOX_SIZE)

uint64_t mds[] = {
    0x00008b48,
    0x00004924,
    0x00002812,
    0x000014c1,
    0x000088b4,
    0x00004492,
    0x00002281,
    0x0000114c,
    0x0000488b,
    0x00002449,
    0x00001228,
    0x0000c114,
    0x0000b488,
    0x00009244,
    0x00008122,
    0x00004c11
};