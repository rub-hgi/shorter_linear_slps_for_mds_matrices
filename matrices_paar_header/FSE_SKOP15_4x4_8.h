#define DIM 32
#define SBOX_SIZE 8
#define NUM_SBOXES (DIM/SBOX_SIZE)

uint64_t mds[] = {
    0x80401090,
    0x4020c080,
    0x20106040,
    0x10c03020,
    0x08040109,
    0x04020c08,
    0x02010604,
    0x010c0302,
    0x40809010,
    0x204080c0,
    0x10204060,
    0xc0102030,
    0x04080901,
    0x0204080c,
    0x01020406,
    0x0c010203,
    0x10908040,
    0xc0804020,
    0x60402010,
    0x302010c0,
    0x01090804,
    0x0c080402,
    0x06040201,
    0x0302010c,
    0x90104080,
    0x80c02040,
    0x40601020,
    0x2030c010,
    0x09010408,
    0x080c0204,
    0x04060102,
    0x02030c01
};