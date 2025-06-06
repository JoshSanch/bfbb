#include "iMix.h"

#include <dolphin/os.h>

static struct MIXChannel __MIXChannel[64];

static U32 __MIXSoundMode;
static S32 __MIXDvdStreamAttenUser;
static S32 __MIXDvdStreamAttenCurrent;

U16 __MIXVolumeTable[] = {
    0x0,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,
    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,
    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,
    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,
    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,    0x1,
    0x1,    0x1,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,
    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,
    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,    0x2,
    0x2,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,
    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,    0x3,
    0x3,    0x3,    0x4,    0x4,    0x4,    0x4,    0x4,    0x4,    0x4,    0x4,    0x4,    0x4,
    0x4,    0x4,    0x4,    0x4,    0x4,    0x4,    0x4,    0x4,    0x4,    0x5,    0x5,    0x5,
    0x5,    0x5,    0x5,    0x5,    0x5,    0x5,    0x5,    0x5,    0x5,    0x5,    0x5,    0x5,
    0x5,    0x6,    0x6,    0x6,    0x6,    0x6,    0x6,    0x6,    0x6,    0x6,    0x6,    0x6,
    0x6,    0x6,    0x7,    0x7,    0x7,    0x7,    0x7,    0x7,    0x7,    0x7,    0x7,    0x7,
    0x7,    0x7,    0x8,    0x8,    0x8,    0x8,    0x8,    0x8,    0x8,    0x8,    0x8,    0x8,
    0x9,    0x9,    0x9,    0x9,    0x9,    0x9,    0x9,    0x9,    0x9,    0xA,    0xA,    0xA,
    0xA,    0xA,    0xA,    0xA,    0xA,    0xA,    0xB,    0xB,    0xB,    0xB,    0xB,    0xB,
    0xB,    0xC,    0xC,    0xC,    0xC,    0xC,    0xC,    0xC,    0xD,    0xD,    0xD,    0xD,
    0xD,    0xD,    0xD,    0xE,    0xE,    0xE,    0xE,    0xE,    0xE,    0xF,    0xF,    0xF,
    0xF,    0xF,    0x10,   0x10,   0x10,   0x10,   0x10,   0x11,   0x11,   0x11,   0x11,   0x11,
    0x12,   0x12,   0x12,   0x12,   0x12,   0x13,   0x13,   0x13,   0x13,   0x13,   0x14,   0x14,
    0x14,   0x14,   0x15,   0x15,   0x15,   0x15,   0x16,   0x16,   0x16,   0x16,   0x17,   0x17,
    0x17,   0x18,   0x18,   0x18,   0x18,   0x19,   0x19,   0x19,   0x1A,   0x1A,   0x1A,   0x1A,
    0x1B,   0x1B,   0x1B,   0x1C,   0x1C,   0x1C,   0x1D,   0x1D,   0x1D,   0x1E,   0x1E,   0x1E,
    0x1F,   0x1F,   0x20,   0x20,   0x20,   0x21,   0x21,   0x21,   0x22,   0x22,   0x23,   0x23,
    0x23,   0x24,   0x24,   0x25,   0x25,   0x26,   0x26,   0x26,   0x27,   0x27,   0x28,   0x28,
    0x29,   0x29,   0x2A,   0x2A,   0x2B,   0x2B,   0x2C,   0x2C,   0x2D,   0x2D,   0x2E,   0x2E,
    0x2F,   0x2F,   0x30,   0x31,   0x31,   0x32,   0x32,   0x33,   0x33,   0x34,   0x35,   0x35,
    0x36,   0x37,   0x37,   0x38,   0x38,   0x39,   0x3A,   0x3A,   0x3B,   0x3C,   0x3D,   0x3D,
    0x3E,   0x3F,   0x3F,   0x40,   0x41,   0x42,   0x42,   0x43,   0x44,   0x45,   0x46,   0x46,
    0x47,   0x48,   0x49,   0x4A,   0x4B,   0x4B,   0x4C,   0x4D,   0x4E,   0x4F,   0x50,   0x51,
    0x52,   0x53,   0x54,   0x55,   0x56,   0x57,   0x58,   0x59,   0x5A,   0x5B,   0x5C,   0x5D,
    0x5E,   0x5F,   0x60,   0x61,   0x62,   0x64,   0x65,   0x66,   0x67,   0x68,   0x6A,   0x6B,
    0x6C,   0x6D,   0x6F,   0x70,   0x71,   0x72,   0x74,   0x75,   0x76,   0x78,   0x79,   0x7B,
    0x7C,   0x7E,   0x7F,   0x80,   0x82,   0x83,   0x85,   0x87,   0x88,   0x8A,   0x8B,   0x8D,
    0x8F,   0x90,   0x92,   0x94,   0x95,   0x97,   0x99,   0x9B,   0x9C,   0x9E,   0xA0,   0xA2,
    0xA4,   0xA6,   0xA8,   0xAA,   0xAB,   0xAD,   0xAF,   0xB2,   0xB4,   0xB6,   0xB8,   0xBA,
    0xBC,   0xBE,   0xC0,   0xC3,   0xC5,   0xC7,   0xCA,   0xCC,   0xCE,   0xD1,   0xD3,   0xD6,
    0xD8,   0xDB,   0xDD,   0xE0,   0xE2,   0xE5,   0xE7,   0xEA,   0xED,   0xF0,   0xF2,   0xF5,
    0xF8,   0xFB,   0xFE,   0x101,  0x104,  0x107,  0x10A,  0x10D,  0x110,  0x113,  0x116,  0x11A,
    0x11D,  0x120,  0x124,  0x127,  0x12A,  0x12E,  0x131,  0x135,  0x138,  0x13C,  0x140,  0x143,
    0x147,  0x14B,  0x14F,  0x153,  0x157,  0x15B,  0x15F,  0x163,  0x167,  0x16B,  0x16F,  0x173,
    0x178,  0x17C,  0x180,  0x185,  0x189,  0x18E,  0x193,  0x197,  0x19C,  0x1A1,  0x1A6,  0x1AB,
    0x1AF,  0x1B4,  0x1BA,  0x1BF,  0x1C4,  0x1C9,  0x1CE,  0x1D4,  0x1D9,  0x1DF,  0x1E4,  0x1EA,
    0x1EF,  0x1F5,  0x1FB,  0x201,  0x207,  0x20D,  0x213,  0x219,  0x21F,  0x226,  0x22C,  0x232,
    0x239,  0x240,  0x246,  0x24D,  0x254,  0x25B,  0x262,  0x269,  0x270,  0x277,  0x27E,  0x286,
    0x28D,  0x295,  0x29D,  0x2A4,  0x2AC,  0x2B4,  0x2BC,  0x2C4,  0x2CC,  0x2D5,  0x2DD,  0x2E6,
    0x2EE,  0x2F7,  0x300,  0x309,  0x312,  0x31B,  0x324,  0x32D,  0x337,  0x340,  0x34A,  0x354,
    0x35D,  0x367,  0x371,  0x37C,  0x386,  0x390,  0x39B,  0x3A6,  0x3B1,  0x3BB,  0x3C7,  0x3D2,
    0x3DD,  0x3E9,  0x3F4,  0x400,  0x40C,  0x418,  0x424,  0x430,  0x43D,  0x449,  0x456,  0x463,
    0x470,  0x47D,  0x48A,  0x498,  0x4A5,  0x4B3,  0x4C1,  0x4CF,  0x4DD,  0x4EC,  0x4FA,  0x509,
    0x518,  0x527,  0x536,  0x546,  0x555,  0x565,  0x575,  0x586,  0x596,  0x5A6,  0x5B7,  0x5C8,
    0x5D9,  0x5EB,  0x5FC,  0x60E,  0x620,  0x632,  0x644,  0x657,  0x66A,  0x67D,  0x690,  0x6A4,
    0x6B7,  0x6CB,  0x6DF,  0x6F4,  0x708,  0x71D,  0x732,  0x748,  0x75D,  0x773,  0x789,  0x79F,
    0x7B6,  0x7CD,  0x7E4,  0x7FB,  0x813,  0x82B,  0x843,  0x85C,  0x874,  0x88E,  0x8A7,  0x8C1,
    0x8DA,  0x8F5,  0x90F,  0x92A,  0x945,  0x961,  0x97D,  0x999,  0x9B5,  0x9D2,  0x9EF,  0xA0D,
    0xA2A,  0xA48,  0xA67,  0xA86,  0xAA5,  0xAC5,  0xAE5,  0xB05,  0xB25,  0xB47,  0xB68,  0xB8A,
    0xBAC,  0xBCF,  0xBF2,  0xC15,  0xC39,  0xC5D,  0xC82,  0xCA7,  0xCCC,  0xCF2,  0xD19,  0xD3F,
    0xD67,  0xD8E,  0xDB7,  0xDDF,  0xE08,  0xE32,  0xE5C,  0xE87,  0xEB2,  0xEDD,  0xF09,  0xF36,
    0xF63,  0xF91,  0xFBF,  0xFEE,  0x101D, 0x104D, 0x107D, 0x10AE, 0x10DF, 0x1111, 0x1144, 0x1177,
    0x11AB, 0x11DF, 0x1214, 0x124A, 0x1280, 0x12B7, 0x12EE, 0x1326, 0x135F, 0x1399, 0x13D3, 0x140D,
    0x1449, 0x1485, 0x14C2, 0x14FF, 0x153E, 0x157D, 0x15BC, 0x15FD, 0x163E, 0x1680, 0x16C3, 0x1706,
    0x174A, 0x178F, 0x17D5, 0x181C, 0x1863, 0x18AC, 0x18F5, 0x193F, 0x198A, 0x19D5, 0x1A22, 0x1A6F,
    0x1ABE, 0x1B0D, 0x1B5D, 0x1BAE, 0x1C00, 0x1C53, 0x1CA7, 0x1CFC, 0x1D52, 0x1DA9, 0x1E01, 0x1E5A,
    0x1EB4, 0x1F0F, 0x1F6B, 0x1FC8, 0x2026, 0x2086, 0x20E6, 0x2148, 0x21AA, 0x220E, 0x2273, 0x22D9,
    0x2341, 0x23A9, 0x2413, 0x247E, 0x24EA, 0x2557, 0x25C6, 0x2636, 0x26A7, 0x271A, 0x278E, 0x2803,
    0x287A, 0x28F2, 0x296B, 0x29E6, 0x2A62, 0x2AE0, 0x2B5F, 0x2BDF, 0x2C61, 0x2CE5, 0x2D6A, 0x2DF1,
    0x2E79, 0x2F03, 0x2F8E, 0x301B, 0x30AA, 0x313A, 0x31CC, 0x325F, 0x32F5, 0x338C, 0x3425, 0x34BF,
    0x355B, 0x35FA, 0x369A, 0x373C, 0x37DF, 0x3885, 0x392C, 0x39D6, 0x3A81, 0x3B2F, 0x3BDE, 0x3C90,
    0x3D43, 0x3DF9, 0x3EB1, 0x3F6A, 0x4026, 0x40E5, 0x41A5, 0x4268, 0x432C, 0x43F4, 0x44BD, 0x4589,
    0x4657, 0x4727, 0x47FA, 0x48D0, 0x49A8, 0x4A82, 0x4B5F, 0x4C3E, 0x4D20, 0x4E05, 0x4EEC, 0x4FD6,
    0x50C3, 0x51B2, 0x52A4, 0x5399, 0x5491, 0x558C, 0x5689, 0x578A, 0x588D, 0x5994, 0x5A9D, 0x5BAA,
    0x5CBA, 0x5DCD, 0x5EE3, 0x5FFC, 0x6119, 0x6238, 0x635C, 0x6482, 0x65AC, 0x66D9, 0x680A, 0x693F,
    0x6A77, 0x6BB2, 0x6CF2, 0x6E35, 0x6F7B, 0x70C6, 0x7214, 0x7366, 0x74BC, 0x7616, 0x7774, 0x78D6,
    0x7A3D, 0x7BA7, 0x7D16, 0x7E88, 0x7FFF, 0x817B, 0x82FB, 0x847F, 0x8608, 0x8795, 0x8927, 0x8ABE,
    0x8C59, 0x8DF9, 0x8F9E, 0x9148, 0x92F6, 0x94AA, 0x9663, 0x9820, 0x99E3, 0x9BAB, 0x9D79, 0x9F4C,
    0xA124, 0xA302, 0xA4E5, 0xA6CE, 0xA8BC, 0xAAB0, 0xACAA, 0xAEAA, 0xB0B0, 0xB2BC, 0xB4CE, 0xB6E5,
    0xB904, 0xBB28, 0xBD53, 0xBF84, 0xC1BC, 0xC3FA, 0xC63F, 0xC88B, 0xCADD, 0xCD37, 0xCF97, 0xD1FE,
    0xD46D, 0xD6E3, 0xD960, 0xDBE4, 0xDE70, 0xE103, 0xE39E, 0xE641, 0xE8EB, 0xEB9E, 0xEE58, 0xF11B,
    0xF3E6, 0xF6B9, 0xF994, 0xFC78, 0xFF64
};

S32 __MIXPanTable[] = {
    0x0,        0x0,        0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFE, 0xFFFFFFFE, 0xFFFFFFFE,
    0xFFFFFFFD, 0xFFFFFFFD, 0xFFFFFFFC, 0xFFFFFFFC, 0xFFFFFFFC, 0xFFFFFFFB, 0xFFFFFFFB, 0xFFFFFFFB,
    0xFFFFFFFA, 0xFFFFFFFA, 0xFFFFFFF9, 0xFFFFFFF9, 0xFFFFFFF9, 0xFFFFFFF8, 0xFFFFFFF8, 0xFFFFFFF7,
    0xFFFFFFF7, 0xFFFFFFF6, 0xFFFFFFF6, 0xFFFFFFF6, 0xFFFFFFF5, 0xFFFFFFF5, 0xFFFFFFF4, 0xFFFFFFF4,
    0xFFFFFFF3, 0xFFFFFFF3, 0xFFFFFFF2, 0xFFFFFFF2, 0xFFFFFFF2, 0xFFFFFFF1, 0xFFFFFFF1, 0xFFFFFFF0,
    0xFFFFFFF0, 0xFFFFFFEF, 0xFFFFFFEF, 0xFFFFFFEE, 0xFFFFFFEE, 0xFFFFFFED, 0xFFFFFFEC, 0xFFFFFFEC,
    0xFFFFFFEB, 0xFFFFFFEB, 0xFFFFFFEA, 0xFFFFFFEA, 0xFFFFFFE9, 0xFFFFFFE9, 0xFFFFFFE8, 0xFFFFFFE7,
    0xFFFFFFE7, 0xFFFFFFE6, 0xFFFFFFE6, 0xFFFFFFE5, 0xFFFFFFE4, 0xFFFFFFE4, 0xFFFFFFE3, 0xFFFFFFE2,
    0xFFFFFFE2, 0xFFFFFFE1, 0xFFFFFFE0, 0xFFFFFFDF, 0xFFFFFFDF, 0xFFFFFFDE, 0xFFFFFFDD, 0xFFFFFFDC,
    0xFFFFFFDC, 0xFFFFFFDB, 0xFFFFFFDA, 0xFFFFFFD9, 0xFFFFFFD8, 0xFFFFFFD8, 0xFFFFFFD7, 0xFFFFFFD6,
    0xFFFFFFD5, 0xFFFFFFD4, 0xFFFFFFD3, 0xFFFFFFD2, 0xFFFFFFD1, 0xFFFFFFD0, 0xFFFFFFCF, 0xFFFFFFCE,
    0xFFFFFFCD, 0xFFFFFFCC, 0xFFFFFFCA, 0xFFFFFFC9, 0xFFFFFFC8, 0xFFFFFFC7, 0xFFFFFFC5, 0xFFFFFFC4,
    0xFFFFFFC3, 0xFFFFFFC1, 0xFFFFFFC0, 0xFFFFFFBE, 0xFFFFFFBD, 0xFFFFFFBB, 0xFFFFFFB9, 0xFFFFFFB8,
    0xFFFFFFB6, 0xFFFFFFB4, 0xFFFFFFB2, 0xFFFFFFB0, 0xFFFFFFAD, 0xFFFFFFAB, 0xFFFFFFA9, 0xFFFFFFA6,
    0xFFFFFFA3, 0xFFFFFFA0, 0xFFFFFF9D, 0xFFFFFF9A, 0xFFFFFF96, 0xFFFFFF92, 0xFFFFFF8D, 0xFFFFFF88,
    0xFFFFFF82, 0xFFFFFF7B, 0xFFFFFF74, 0xFFFFFF6A, 0xFFFFFF5D, 0xFFFFFF4C, 0xFFFFFF2E, 0xFFFFFC78
};

S16 __MIX_DPL2_front[] = {
    0x0,    0x0,    0x0,    0x0,    0x0,    0x0,    0x0,    0x0,    0x0,    0x0,    0x0,    0x0,
    0x0,    0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFF, 0xFFFE, 0xFFFE, 0xFFFE, 0xFFFE, 0xFFFD, 0xFFFD,
    0xFFFD, 0xFFFC, 0xFFFC, 0xFFFC, 0xFFFB, 0xFFFB, 0xFFFA, 0xFFFA, 0xFFFA, 0xFFF9, 0xFFF9, 0xFFF8,
    0xFFF8, 0xFFF7, 0xFFF7, 0xFFF6, 0xFFF5, 0xFFF5, 0xFFF4, 0xFFF4, 0xFFF3, 0xFFF2, 0xFFF2, 0xFFF1,
    0xFFF0, 0xFFEF, 0xFFEF, 0xFFEE, 0xFFED, 0xFFEC, 0xFFEB, 0xFFEB, 0xFFEA, 0xFFE9, 0xFFE8, 0xFFE7,
    0xFFE6, 0xFFE5, 0xFFE4, 0xFFE3, 0xFFE2, 0xFFE1, 0xFFE0, 0xFFDE, 0xFFDD, 0xFFDC, 0xFFDB, 0xFFDA,
    0xFFD8, 0xFFD7, 0xFFD6, 0xFFD4, 0xFFD3, 0xFFD1, 0xFFD0, 0xFFCE, 0xFFCC, 0xFFCB, 0xFFC9, 0xFFC7,
    0xFFC6, 0xFFC4, 0xFFC2, 0xFFC0, 0xFFBE, 0xFFBC, 0xFFBA, 0xFFB7, 0xFFB5, 0xFFB3, 0xFFB0, 0xFFAE,
    0xFFAB, 0xFFA8, 0xFFA6, 0xFFA3, 0xFFA0, 0xFF9C, 0xFF99, 0xFF96, 0xFF92, 0xFF8E, 0xFF8A, 0xFF86,
    0xFF82, 0xFF7D, 0xFF78, 0xFF73, 0xFF6E, 0xFF68, 0xFF61, 0xFF5A, 0xFF53, 0xFF4B, 0xFF42, 0xFF37,
    0xFF2C, 0xFF1F, 0xFF0F, 0xFEFB, 0xFEE2, 0xFEBF, 0xFE83, 0xFC40
};

S16 __MIX_DPL2_rear[] = {
    0xFFC3, 0xFFC3, 0xFFC4, 0xFFC5, 0xFFC5, 0xFFC6, 0xFFC6, 0xFFC7, 0xFFC8, 0xFFC8, 0xFFC9, 0xFFC9,
    0xFFCA, 0xFFCB, 0xFFCB, 0xFFCC, 0xFFCC, 0xFFCD, 0xFFCE, 0xFFCE, 0xFFCF, 0xFFCF, 0xFFD0, 0xFFD0,
    0xFFD1, 0xFFD1, 0xFFD2, 0xFFD2, 0xFFD3, 0xFFD3, 0xFFD4, 0xFFD4, 0xFFD5, 0xFFD5, 0xFFD6, 0xFFD6,
    0xFFD7, 0xFFD7, 0xFFD8, 0xFFD8, 0xFFD9, 0xFFD9, 0xFFDA, 0xFFDA, 0xFFDA, 0xFFDB, 0xFFDB, 0xFFDC,
    0xFFDC, 0xFFDD, 0xFFDD, 0xFFDD, 0xFFDE, 0xFFDE, 0xFFDF, 0xFFDF, 0xFFE0, 0xFFE0, 0xFFE0, 0xFFE1,
    0xFFE1, 0xFFE1, 0xFFE2, 0xFFE2, 0xFFE3, 0xFFE3, 0xFFE3, 0xFFE4, 0xFFE4, 0xFFE4, 0xFFE5, 0xFFE5,
    0xFFE5, 0xFFE6, 0xFFE6, 0xFFE6, 0xFFE7, 0xFFE7, 0xFFE7, 0xFFE8, 0xFFE8, 0xFFE8, 0xFFE9, 0xFFE9,
    0xFFE9, 0xFFEA, 0xFFEA, 0xFFEA, 0xFFEB, 0xFFEB, 0xFFEB, 0xFFEC, 0xFFEC, 0xFFEC, 0xFFEC, 0xFFED,
    0xFFED, 0xFFED, 0xFFEE, 0xFFEE, 0xFFEE, 0xFFEE, 0xFFEF, 0xFFEF, 0xFFEF, 0xFFEF, 0xFFF0, 0xFFF0,
    0xFFF0, 0xFFF0, 0xFFF1, 0xFFF1, 0xFFF1, 0xFFF1, 0xFFF2, 0xFFF2, 0xFFF2, 0xFFF2, 0xFFF3, 0xFFF3,
    0xFFF3, 0xFFF3, 0xFFF3, 0xFFF4, 0xFFF4, 0xFFF4, 0xFFF4, 0xFFF5
};

static S32 __MIXGetVolume(S32 param_1)
{
    if (param_1 <= -0x388)
    {
        return 0;
    }
    if (0x3c <= param_1)
    {
        return 0xff64;
    }
    return __MIXVolumeTable[param_1 + 0x388];
}

static void __MIXSetPan(struct MIXChannel* param_1)
{
    S32 iVar1 = *(S32*)((S32)param_1 + 0x14);
    S32 iVar3 = *(S32*)((S32)param_1 + 0x18);

    S32 iVar2 = 0x7f - iVar1;
    S32 iVar4 = 0x7f - iVar3;

    if (__MIXSoundMode == 3)
    {
        param_1->data[8] = __MIX_DPL2_front[iVar1];
        param_1->data[9] = __MIX_DPL2_front[iVar2];
        param_1->data[10] = __MIX_DPL2_front[iVar4];
        param_1->data[11] = __MIX_DPL2_front[iVar3];
        param_1->data[12] = __MIX_DPL2_rear[iVar2];
        param_1->data[13] = __MIX_DPL2_rear[iVar1];
    }
    else
    {
        param_1->data[8] = __MIXPanTable[iVar1];
        param_1->data[9] = __MIXPanTable[iVar2];
        param_1->data[10] = __MIXPanTable[iVar4];
        param_1->data[11] = __MIXPanTable[iVar3];
    }
}

static void __MIXResetChannel(struct MIXChannel* param_1)
{
    param_1->data[1] = 0x50000000;
    param_1->data[2] = 0;
    param_1->data[3] = 0xfffffc40;
    param_1->data[4] = 0xfffffc40;
    param_1->data[7] = 0;
    param_1->data[5] = 0x40;
    param_1->data[6] = 0x7f;

    *(U16*)(&param_1->data[23]) = 0;
    *(U16*)(&param_1->data[22]) = 0;
    *(U16*)(&param_1->data[21]) = 0;
    *(U16*)(&param_1->data[20]) = 0;
    *(U16*)(&param_1->data[19]) = 0;
    *(U16*)(&param_1->data[18]) = 0;
    *(U16*)(&param_1->data[17]) = 0;
    *(U16*)(&param_1->data[16]) = 0;
    *(U16*)(&param_1->data[15]) = 0;
    *(U16*)(&param_1->data[14]) = 0;

    __MIXSetPan(param_1);
}

static S32 __MIXClampPan(S32 param)
{
    S32 retval;
    if (param < 0)
    {
        return 0;
    }
    retval = 0x7f;
    if (param <= 0x7f)
    {
        retval = param;
    }
    return retval;
}

void MIXInit()
{
    S32 iVar1 = 0;
    do
    {
        __MIXResetChannel((struct MIXChannel*)&__MIXChannel[iVar1]);
        iVar1++;
    } while (iVar1 < 0x40);
    __MIXDvdStreamAttenCurrent = 0;
    __MIXDvdStreamAttenUser = 0;
    __MIXSoundMode = OSGetSoundMode();
}

void MIXInitChannel(S32 param_1, U32 param_2, S32 param_3, S32 param_4, S32 param_5,
                    S32 param_6, S32 param_7, S32 param_8)
{
    struct MIXChannel* chan = &__MIXChannel[*(S32*)((S32)param_1 + 0x18)];

    S16 sVar1;
    U16 uVar2;
    U16 uVar4;

    chan->data[0] = param_1;
    chan->data[1] = param_2 & 7;
    chan->data[2] = param_3;
    chan->data[3] = param_4;
    chan->data[4] = param_5;
    chan->data[5] = param_6;
    chan->data[6] = param_7;
    chan->data[7] = param_8;

    __MIXSetPan(chan);

    if ((chan->data[1] & 4))
    {
        *(U16*)(&chan->data[0xe]) = 0;
    }
    else
    {
        *(U16*)(&chan->data[0xe]) = __MIXGetVolume(param_3);
    }

    uVar4 = 0;

    switch ((S32)__MIXSoundMode)
    {
    case 0:
        *(U16*)(&chan->data[0xf]) = __MIXGetVolume(chan->data[7] + chan->data[10]);
        *(U16*)(&chan->data[0x10]) = __MIXGetVolume(chan->data[7] + chan->data[10]);
        *(U16*)(&chan->data[0x11]) = __MIXGetVolume(chan->data[7] + chan->data[11]);
        if ((chan->data[1] & 1U) != 0)
        {
            *(U16*)(&chan->data[0x12]) = __MIXGetVolume(chan->data[3] + chan->data[10]);
            *(U16*)(&chan->data[0x13]) = __MIXGetVolume(chan->data[3] + chan->data[10]);
            *(U16*)(&chan->data[0x14]) =
                __MIXGetVolume(chan->data[3] + chan->data[11] - 0x3c);
        }
        else
        {
            *(U16*)(&chan->data[0x12]) =
                __MIXGetVolume(chan->data[7] + chan->data[3] + chan->data[10]);
            *(U16*)(&chan->data[0x13]) =
                __MIXGetVolume(chan->data[7] + chan->data[3] + chan->data[10]);
            *(U16*)(&chan->data[0x14]) =
                __MIXGetVolume(chan->data[7] + chan->data[3] + chan->data[11] - 0x3c);
        }

        if (chan->data[1] & 2)
        {
            *(U16*)(&chan->data[0x15]) = __MIXGetVolume(chan->data[4] + chan->data[10]);
            *(U16*)(&chan->data[0x16]) = __MIXGetVolume(chan->data[4] + chan->data[10]);
            *(U16*)(&chan->data[0x17]) =
                __MIXGetVolume(chan->data[4] + chan->data[11] - 0x3c);
        }
        else
        {
            *(U16*)(&chan->data[0x15]) =
                __MIXGetVolume(chan->data[7] + chan->data[4] + chan->data[10]);
            *(U16*)(&chan->data[0x16]) =
                __MIXGetVolume(chan->data[7] + chan->data[4] + chan->data[10]);
            *(U16*)(&chan->data[0x17]) =
                __MIXGetVolume(chan->data[7] + chan->data[4] + chan->data[11] - 0x3c);
        }

        break;

    case 1:
    case 2:
        *(U16*)(&chan->data[0xf]) =
            __MIXGetVolume(chan->data[7] + chan->data[8] + chan->data[10]);
        *(U16*)(&chan->data[0x10]) =
            __MIXGetVolume(chan->data[7] + chan->data[9] + chan->data[10]);
        *(U16*)(&chan->data[0x11]) = __MIXGetVolume(chan->data[7] + chan->data[11]);

        if ((chan->data[1] & 1) != 0)
        {
            *(U16*)(&chan->data[0x12]) =
                __MIXGetVolume(chan->data[3] + chan->data[8] + chan->data[10]);
            *(U16*)(&chan->data[0x13]) =
                __MIXGetVolume(chan->data[3] + chan->data[9] + chan->data[10]);
            *(U16*)(&chan->data[0x14]) =
                __MIXGetVolume(chan->data[3] + chan->data[11] + -0x3c);
        }
        else
        {
            *(U16*)(&chan->data[0x12]) =
                __MIXGetVolume(chan->data[7] + chan->data[3] + chan->data[8] + chan->data[10]);
            *(U16*)(&chan->data[0x13]) =
                __MIXGetVolume(chan->data[7] + chan->data[3] + chan->data[9] + chan->data[10]);
            *(U16*)(&chan->data[0x14]) =
                __MIXGetVolume(chan->data[7] + chan->data[3] + chan->data[11] + -0x3c);
        }
        if (chan->data[1] & 2U)
        {
            *(U16*)(&chan->data[0x15]) =
                __MIXGetVolume(chan->data[4] + chan->data[8] + chan->data[10]);
            *(U16*)(&chan->data[0x16]) =
                __MIXGetVolume(chan->data[4] + chan->data[9] + chan->data[10]);
            *(U16*)(&chan->data[0x17]) =
                __MIXGetVolume(chan->data[4] + chan->data[11] + -0x3c);
        }
        else
        {
            *(U16*)(&chan->data[0x15]) =
                __MIXGetVolume(chan->data[7] + chan->data[4] + chan->data[8] + chan->data[10]);
            *(U16*)(&chan->data[0x16]) =
                __MIXGetVolume(chan->data[7] + chan->data[4] + chan->data[9] + chan->data[10]);
            *(U16*)(&chan->data[0x17]) =
                __MIXGetVolume(chan->data[7] + chan->data[4] + chan->data[11] + -0x3c);
        }
        break;

    case 3:
        *(U16*)(&chan->data[15]) =
            __MIXGetVolume(chan->data[7] + chan->data[8] + chan->data[10]);
        *(U16*)(&chan->data[16]) =
            __MIXGetVolume(chan->data[7] + chan->data[9] + chan->data[10]);
        *(U16*)(&chan->data[21]) =
            __MIXGetVolume(chan->data[7] + chan->data[12] + chan->data[11]);
        *(U16*)(&chan->data[22]) =
            __MIXGetVolume(chan->data[7] + chan->data[13] + chan->data[11]);

        if ((chan->data[1] & 1) != 0)
        {
            *(U16*)(&chan->data[0x12]) =
                __MIXGetVolume(chan->data[3] + chan->data[8] + chan->data[10]);
            *(U16*)(&chan->data[0x13]) =
                __MIXGetVolume(chan->data[3] + chan->data[9] + chan->data[10]);
            *(U16*)(&chan->data[0x14]) =
                __MIXGetVolume(chan->data[3] + chan->data[12] + chan->data[11]);
            *(U16*)(&chan->data[0x17]) =
                __MIXGetVolume(chan->data[3] + chan->data[13] + chan->data[11]);
        }
        else
        {
            *(U16*)(&chan->data[0x12]) =
                __MIXGetVolume(chan->data[7] + chan->data[3] + chan->data[8] + chan->data[10]);
            *(U16*)(&chan->data[0x13]) =
                __MIXGetVolume(chan->data[7] + chan->data[3] + chan->data[9] + chan->data[10]);
            *(U16*)(&chan->data[0x14]) =
                __MIXGetVolume(chan->data[7] + chan->data[3] + chan->data[12] + chan->data[11]);
            *(U16*)(&chan->data[0x17]) =
                __MIXGetVolume(chan->data[7] + chan->data[3] + chan->data[13] + chan->data[11]);
        }

        uVar4 |= 0x4000;
        break;
    }

    int enabled = OSDisableInterrupts();

    *(U16*)((S32)param_1 + 0x19c) =
        *(U16*)(&chan->data[0xe]);
    *(U16*)((S32)param_1 + 0x19e) = 0;
    if ((*(U16*)((S32)param_1 + 0x14a) = *(U16*)(&chan->data[0xf])))
    {
        uVar4 |= 1;
    }

    *(U16*)((S32)param_1 + 0x14c) = 0;

    if (*(U16*)(param_1 + 0x14e) = *(U16*)(&chan->data[0x10]))
    {
        uVar4 = uVar4 | 2;
    }
    *(U16*)(param_1 + 0x150) = 0;

    if (*(U16*)(param_1 + 0x152) = *(U16*)(&chan->data[0x12]))
    {
        uVar4 = uVar4 | 0x10;
    }
    *(U16*)((S32)param_1 + 0x154) = 0;
    if (*(U16*)((S32)param_1 + 0x156) = *(U16*)(&chan->data[0x13]))
    {
        uVar4 = uVar4 | 0x20;
    }
    *(U16*)(param_1 + 0x158) = 0;
    if (*(U16*)(param_1 + 0x15a) = *(U16*)(&chan->data[0x15]))
    {
        uVar4 = uVar4 | 0x200;
    }
    *(U16*)(param_1 + 0x15c) = 0;
    if (*(U16*)(param_1 + 0x15e) = *(U16*)(&chan->data[0x16]))
    {
        uVar4 = uVar4 | 0x400;
    }
    *(U16*)(param_1 + 0x160) = 0;
    ;
    if (*(U16*)(param_1 + 0x162) = *(U16*)(&chan->data[0x17]))
    {
        uVar4 = uVar4 | 0x1000;
    }
    *(U16*)(param_1 + 0x164) = 0;
    if (*(U16*)((S32)param_1 + 0x166) = *(U16*)(&chan->data[0x11]))
    {
        uVar4 = uVar4 | 4;
    }
    *(U16*)((S32)param_1 + 0x168) = 0;
    if (*(U16*)(param_1 + 0x16a) = *(U16*)(&chan->data[0x14]))
    {
        uVar4 = uVar4 | 0x80;
    }
    *(U16*)((S32)param_1 + 0x16c) = 0;
    *(U16*)((S32)param_1 + 0x144) = uVar4;
    *(U32*)(param_1 + 0x1c) = *(U32*)(param_1 + 0x1c) | 0x212;
    OSRestoreInterrupts(enabled);
}

void MIXReleaseChannel(S32* param_1)
{
    __MIXChannel[*(param_1 + 6)].data[0] = 0;
}

void MIXAdjustInput(S32* param_1, S32 param_2)
{
    S32* handle = &__MIXChannel[*(param_1 + 6)].data[0];
    handle[2] += param_2;
    handle[1] |= 0x10000000;
}

S32 MIXGetInput(S32* param_1)
{
    S32* handle = &__MIXChannel[*(param_1 + 6)].data[0];
    return handle[2];
}

void MIXAdjustPan(S32* param_1, S32 param_2)
{
    struct MIXChannel* chan = &__MIXChannel[*(param_1 + 6)];
    S32* handle = &__MIXChannel[*(param_1 + 6)].data[0];
    handle[5] = __MIXClampPan(handle[5] + param_2);
    __MIXSetPan(chan);
    chan->data[1] |= 0x40000000;
}

S32 MIXGetPan(S32* param_1)
{
    S32* handle = &__MIXChannel[*(param_1 + 6)].data[0];
    return handle[5];
}

void MIXUnMute(S32* param_1)
{
    S32* handle = &__MIXChannel[*(param_1 + 6)].data[0];
    handle[1] &= ~4;
    handle[1] |= 0x10000000;
}

void MIXAdjustFader(S32* param_1, S32 param_2)
{
    S32* handle = &__MIXChannel[*(param_1 + 6)].data[0];
    handle[7] += param_2;
    handle[1] |= 0x40000000;
}

S32 MIXGetFader(S32* param_1)
{
    S32* handle = &__MIXChannel[*(param_1 + 6)].data[0];
    return handle[7];
}