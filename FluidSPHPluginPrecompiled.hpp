#pragma once
// We use precompiled headers to support the fastest possible compilation
#define ZeroImportDll
#include <Windows.h>
#include "Zilch.hpp"
#include "ZeroEngine.hpp"

// Declares our Zilch Library (where we get LibraryBuilder from)
// This also handles the plugin initialization
ZilchDeclareStaticLibraryAndPlugin(FluidSPHPluginLibrary, FluidSPHPluginPlugin);

// We also encourage including all files within here, rather than within headers
// If another project needed to include any headers from our project, then they would simply
// include our precompiled header instead (ideally within their own precompiled header)
// This also must means you must order headers in dependency order (who depends on who)
#include "CL/cl.hpp"
#include "util.h"
#include "CLWrapper.h"
#include "FluidSPHPlugin.hpp"
// Auto Includes (used by Visual Studio plugins, do not remove this line)
