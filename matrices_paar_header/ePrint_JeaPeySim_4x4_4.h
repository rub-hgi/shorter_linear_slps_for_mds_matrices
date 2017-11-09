#define DIM 16
#define SBOX_SIZE 4
#define NUM_SBOXES (DIM/SBOX_SIZE)

uint64_t mds[] = {
    0x0000884b,
    0x00004429,
    0x00002218,
    0x000011c4,
    0x000084b8,
    0x00004294,
    0x00002182,
    0x00001c41,
    0x00008b84,
    0x00004942,
    0x00002821,
    0x0000141c,
    0x00004888,
    0x00002444,
    0x00001222,
    0x0000c111
};