#define DIM 32
#define SBOX_SIZE 4
#define NUM_SBOXES (DIM/SBOX_SIZE)

uint64_t mds[] = {
    0x24398cde,
    0x12d846a7,
    0xc1a4235f,
    0x6c521deb,
    0x42ced389,
    0x2167ad48,
    0x1c3f5a24,
    0xc6dbe512,
    0x9ed2c834,
    0x87a164d2,
    0x4f5c32a1,
    0x2be6d15c,
    0xd8934e2c,
    0xa48d2716,
    0x524a1fc3,
    0xe125cb6d,
    0x8dec2943,
    0x4a76182d,
    0x25f3c41a,
    0x1ebd62c5,
    0x3c2de498,
    0xd61a7284,
    0xa3c5f142,
    0x5d6ebc21,
    0xe9843dc2,
    0x7842da61,
    0xf421a53c,
    0xb21c5ed6,
    0xc34892ed,
    0x6d24817a,
    0x3a124cf5,
    0xd5c126be
};