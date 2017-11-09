#define DIM 16
#define SBOX_SIZE 4
#define NUM_SBOXES (DIM/SBOX_SIZE)

uint64_t mds[] = {
    0x0000829b,
    0x00004189,
    0x00002c48,
    0x00001624,
    0x000028b9,
    0x00001498,
    0x0000c284,
    0x00006142,
    0x00009b82,
    0x00008941,
    0x0000482c,
    0x00002416,
    0x0000b928,
    0x00009814,
    0x000084c2,
    0x00004261
};