// wavefront_reconstruction.cpp : Defines the exported functions for the DLL.
//

#include "pch.h"
#include "framework.h"
#include "wavefront_reconstruction.h"


// This is an example of an exported variable
WAVEFRONTRECONSTRUCTION_API int nwavefrontreconstruction=0;

// This is an example of an exported function.
WAVEFRONTRECONSTRUCTION_API int fnwavefrontreconstruction(void)
{
    return 0;
}

// This is the constructor of a class that has been exported.
Cwavefrontreconstruction::Cwavefrontreconstruction()
{
    return;
}
