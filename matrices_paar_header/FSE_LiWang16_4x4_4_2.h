#define DIM 16
#define SBOX_SIZE 4
#define NUM_SBOXES (DIM/SBOX_SIZE)

uint64_t mds[] = {
    0x00001888,
    0x0000a444,
    0x00008222,
    0x00004111,
    0x00008816,
    0x000044a8,
    0x00002287,
    0x00001142,
    0x00008681,
    0x0000484a,
    0x00002728,
    0x00001214,
    0x00008168,
    0x00004a84,
    0x00002872,
    0x00001421
};
