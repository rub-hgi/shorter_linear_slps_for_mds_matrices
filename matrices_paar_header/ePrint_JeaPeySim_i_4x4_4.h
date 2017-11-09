#define DIM 16
#define SBOX_SIZE 4
#define NUM_SBOXES (DIM/SBOX_SIZE)

uint64_t mds[] = {
    0x000048b8,
    0x00002494,
    0x00001282,
    0x0000c141,
    0x0000829b,
    0x00004189,
    0x00002c48,
    0x00001624,
    0x00008f28,
    0x00004b14,
    0x000029c2,
    0x00001861,
    0x00009884,
    0x00008442,
    0x00004221,
    0x0000211c
};