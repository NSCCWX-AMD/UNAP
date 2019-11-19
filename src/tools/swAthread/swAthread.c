#include "swAthread.h"
#include <athread.h>

void swacc_init()
{
	athread_init();
}

void swacc_end()
{
	athread_halt();
}

void swacc_init_()
{
	athread_init();
}

void swacc_end_()
{
	athread_halt();
}

