#pragma once

// For more information on binding and using Zilch APIs, visit: http://zilch.digipen.edu/
// For auto binding specifically, visit: http://zilch.digipen.edu/AutomaticBinding.html

// An example component being bound to the engine
class FluidSPHPlugin : public ZeroEngine::ZilchComponent
{
public:
  ZilchDeclareDerivedType(FluidSPHPlugin, ZeroEngine::ZilchComponent);
  
  FluidSPHPlugin();
  ~FluidSPHPlugin();
  
  void Initialize(ZeroEngine::CogInitializer* initializer);
  void OnLogicUpdate(ZeroEngine::UpdateEvent* event);

  /* MEMBER VARIABLES */
  // Number of steps for the fluid solver to iterate
  float NumberofSteps;
  // Size of the smoothing kernel
  float InteractionRadius;
  // The approximate starting density of the fluid volume
  float ReferenceDensity;
  // Physical property based on the material type of the fluid being simulated
  float BulkModulus;
  // Defines resistance of a fluid material to shear forces and changes in velocity
  float Viscosity;
  // Ideal gas constant
  float GasConstant;
  // 2D Vector of external forces applied on fluid
  Zilch::Real2 ExternalForces;
  // Defines how much energy is lost upon contact with the boundary volume
  float Restitution;
  // Half length of the bounding volume along X axis
  float HalfGridBoundsX;
  // Half length of the bounding volume along Y axis
  float HalfGridBoundsY;
  // Name of fluid particle archetype to spawn
  Zero::String FluidParticleArchetypeName;

private:
	void PlaceParticles();
	void NormalizeMass();
	void InitKernel();
	void ComputeInitalDensity();
	void UpdateArray();
	int BoxIndicator(float xVal, float yVal);
	int CircleIndictor(float xVal, float yVal);

	// Array of Handles to the fluid particle gameobjects
	std::vector<Zilch::HandleOf<ZeroEngine::Cog>> FluidParticlesArray;

	std::vector<Particle> HostParticles;
	SystemParameters HostParameters;

	// Fluid Particle properties
	std::vector<cl_float> Position;
	std::vector<cl_float> VelocityFullStep;
	std::vector<cl_float> VelocityHalfStep;
	std::vector<cl_float> Acceleration;
	std::vector<cl_float> Density;

	// Number of particles to initialize  
	int ParticleCount;
	// Mass of each fluid particle - assumed to be 1 in beginning and recalculated
	float ParticleMass = 1.0;
	// Inital density calculation constant kernel term
	float ConstantDensityKernelTerm;
	// Initial density calculation constant sum term
	float ConstantDensitySumTerm;

	// OpenCL wrapper object
	CLWrapper clObject;
	FILE* outputStream;
};

// An example of a custom event that we can send
class FluidSPHPluginEvent : public ZeroEngine::ZilchEvent
{
public:
  ZilchDeclareDerivedType(FluidSPHPluginEvent, ZeroEngine::ZilchEvent);
  
  int mLives;
};
