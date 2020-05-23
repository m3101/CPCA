#include <stdlib.h>
#include <stdio.h>
#include "gfx.h"
#include "../src/cpca.h"
#include "../src/linal.h"
#include "../lib/bmpreader.h"

int main()
{
    readbmp();
    gfx_open();
    return 0;
}