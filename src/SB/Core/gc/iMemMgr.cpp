#include "iMemMgr.h"
#include "iSystem.h"
#include "xMemMgr.h"

#include <types.h>
#include <PowerPC_EABI_Support/MSL_C/MSL_Common/stdlib.h>
#include <dolphin.h>

U32 mem_top_alloc;
U32 mem_base_alloc;
volatile OSHeapHandle the_heap;
OSHeapHandle hs;
OSHeapHandle he;
U32 HeapSize;
extern xMemInfo_tag gMemInfo;
extern unsigned char _stack_addr[];

// Starts going wrong after the if and else statement, everything else before looks fine.
void iMemInit()
{
    OSHeapHandle hi = (OSHeapHandle)OSGetArenaHi();
    he = hi & 0xffffffe0;
    hs = (OSHeapHandle)OSGetArenaLo() + 0x1f & 0xffffffe0;
    the_heap = OSCreateHeap(OSInitAlloc((void*)hs, (void*)(hi & 0xffffffe0), 1), (void*)he);
    OSHeapHandle currHeap = the_heap;
    if (currHeap >= 0)
    {
        OSSetCurrentHeap(currHeap);
    }
    else
    {
        exit(-5);
    }
    gMemInfo.system.addr = 0;
    gMemInfo.system.size = 0x100000;
    gMemInfo.system.flags = 0x20;
    gMemInfo.stack.addr = (U32)&_stack_addr;
    gMemInfo.stack.size = 0xffff8000;
    gMemInfo.stack.flags = gMemInfo.DRAM.flags = 0x820;
    HeapSize = 0x384000;
    gMemInfo.DRAM.addr = (U32)OSAllocFromHeap(__OSCurrHeap, 0x384000);
    gMemInfo.DRAM.size = HeapSize;
    gMemInfo.DRAM.flags = 0x820;
    gMemInfo.SRAM.addr = 0;
    gMemInfo.SRAM.size = 0x200000;
    gMemInfo.SRAM.flags = 0x660;
    mem_top_alloc = gMemInfo.DRAM.addr + HeapSize;
    mem_base_alloc = gMemInfo.DRAM.addr;
}

void iMemExit()
{
    free((void*)gMemInfo.DRAM.addr);
    gMemInfo.DRAM.addr = 0;
}
