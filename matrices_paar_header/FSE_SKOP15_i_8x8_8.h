#define DIM 64
#define SBOX_SIZE 8
#define NUM_SBOXES (DIM/SBOX_SIZE)

uint64_t mds[] = {
    0x8040c089200ea087,
    0x4020608710075080,
    0x2010308008c02840,
    0x1008184004601420,
    0x8040c2002300a10,
    0x402061001180508,
    0x2010308c30cc104,
    0x1c3c204a206a302,
    0x408089c00e2087a0,
    0x2040876007108050,
    0x10208030c0084028,
    0x810401860042014,
    0x408200c3002100a,
    0x204100618010805,
    0x10208030cc304c1,
    0xc30104c206a202a3,
    0xc0898040a087200e,
    0x6087402050801007,
    0x30802010284008c0,
    0x1840100814200460,
    0xc2008040a100230,
    0x610040205080118,
    0x3080201c104c30c,
    0xc20401c3a302a206,
    0x89c0408087a00e20,
    0x8760204080500710,
    0x803010204028c008,
    0x4018081020146004,
    0x200c0408100a3002,
    0x1006020408051801,
    0x803010204c10cc3,
    0x4c2c30102a306a2,
    0x200ea0878040c089,
    0x1007508040206087,
    0x8c0284020103080,
    0x460142010081840,
    0x2300a1008040c20,
    0x118050804020610,
    0xc30cc10402010308,
    0xa206a30201c3c204,
    0xe2087a0408089c0,
    0x710805020408760,
    0xc008402810208030,
    0x6004201408104018,
    0x3002100a0408200c,
    0x1801080502041006,
    0xcc304c101020803,
    0x6a202a3c30104c2,
    0xa087200ec0898040,
    0x5080100760874020,
    0x284008c030802010,
    0x1420046018401008,
    0xa1002300c200804,
    0x508011806100402,
    0xc104c30c03080201,
    0xa302a206c20401c3,
    0x87a00e2089c04080,
    0x8050071087602040,
    0x4028c00880301020,
    0x2014600440180810,
    0x100a3002200c0408,
    0x805180110060204,
    0x4c10cc308030102,
    0x2a306a204c2c301
};