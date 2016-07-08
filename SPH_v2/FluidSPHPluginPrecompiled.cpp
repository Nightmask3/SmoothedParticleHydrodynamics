#include "FluidSPHPluginPrecompiled.hpp"

//***************************************************************************
ZilchDefineStaticLibraryAndPlugin(FluidSPHPluginLibrary, FluidSPHPluginPlugin, ZilchDependencyStub(ZeroEngine))
{
  ZilchInitializeType(FluidSPHPlugin);
  ZilchInitializeType(FluidSPHPluginEvent);
  // Auto Initialize (used by Visual Studio plugins, do not remove this line)
}

//***************************************************************************
void FluidSPHPluginPlugin::Initialize(Zilch::BuildEvent* event)
{
  // One time startup logic goes here
  // This runs after our plugin library/reflection is built
  Zilch::Console::WriteLine("FluidSPHPluginPlugin::Initialize");
}

//***************************************************************************
void FluidSPHPluginPlugin::Uninitialize()
{
  // One time shutdown logic goes here
  Zilch::Console::WriteLine("FluidSPHPluginPlugin::Uninitialize");
}