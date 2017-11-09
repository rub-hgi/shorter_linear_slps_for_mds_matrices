#define DIM 16
#define SBOX_SIZE 4
#define NUM_SBOXES (DIM/SBOX_SIZE)

uint64_t mds[] = {
    0x00008419,
    0x000042c8,
    0x00002164,
    0x00001c32,
    0x00004891,
    0x0000248c,
    0x00001246,
    0x0000c123,
    0x00001984,
    0x0000c842,
    0x00006421,
    0x0000321c,
    0x00009148,
    0x00008c24,
    0x00004612,
    0x000023c1
};