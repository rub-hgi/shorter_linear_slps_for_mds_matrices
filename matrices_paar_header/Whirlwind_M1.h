#define DIM 32
#define SBOX_SIZE 4
#define NUM_SBOXES (DIM/SBOX_SIZE)

uint64_t mds[] = {
    0xa72e8cf1,
    0x5f1746bc,
    0xebcf2396,
    0x796b1d83,
    0x7ae2c81f,
    0xf57164cb,
    0xbefc3269,
    0x97b6d138,
    0x2ea7f18c,
    0x175fbc46,
    0xcfeb9623,
    0x6b79831d,
    0xe27a1fc8,
    0x71f5cb64,
    0xfcbe6932,
    0xb69738d1,
    0x8cf1a72e,
    0x46bc5f17,
    0x2396ebcf,
    0x1d83796b,
    0xc81f7ae2,
    0x64cbf571,
    0x3269befc,
    0xd13897b6,
    0xf18c2ea7,
    0xbc46175f,
    0x9623cfeb,
    0x831d6b79,
    0x1fc8e27a,
    0xcb6471f5,
    0x6932fcbe,
    0x38d1b697
};