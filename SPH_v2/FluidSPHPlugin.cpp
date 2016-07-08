#include "FluidSPHPluginPrecompiled.hpp"
#include "SPH_Kernel.cl"
#define PI 3.14159265359f
//***************************************************************************
ZilchDefineType(FluidSPHPlugin, FluidSPHPluginLibrary, builder, type)
{
  // This is required for component binding
  ZilchBindDestructor();
  ZilchBindConstructor();
  ZilchBindMethod(Initialize);
  
  // Note: All event connection methods must be bound
  ZilchBindMethod(OnLogicUpdate);
  ZilchBindFieldProperty(NumberofSteps);
  ZilchBindFieldProperty(InteractionRadius);
  ZilchBindFieldProperty(ReferenceDensity);
  ZilchBindFieldProperty(BulkModulus);
  ZilchBindFieldProperty(Viscosity);
  ZilchBindFieldProperty(ExternalForces);
  ZilchBindFieldProperty(GasConstant);
  ZilchBindFieldProperty(Restitution);
  ZilchBindFieldProperty(HalfGridBoundsX);
  ZilchBindFieldProperty(HalfGridBoundsY);
  ZilchBindFieldProperty(FluidParticleArchetypeName);

  // Using Property at the end is the same as the [Property] attribute
  // You could also use ->AddAttribute after the bind macro
}

//***************************************************************************
FluidSPHPlugin::FluidSPHPlugin()
{
  //Zilch::Console::WriteLine("FluidSPHPlugin::FluidSPHPlugin (Constructor)");
}

//***************************************************************************
FluidSPHPlugin::~FluidSPHPlugin()
{
  //Zilch::Console::WriteLine("FluidSPHPlugin::~FluidSPHPlugin (Destructor)");
}

//***************************************************************************
void FluidSPHPlugin::Initialize(ZeroEngine::CogInitializer* initializer)
{
  //Zilch::Console::WriteLine("FluidSPHPlugin::Initialize");

  HostParameters.Epsilon = 0.0001f;
  HostParameters.DeltaTime = 0.016f;
  HostParameters.Gravity.s[0] = ExternalForces.x;
  HostParameters.Gravity.s[1] = ExternalForces.y;

  HostParameters.WorldSize.s[0] = 2 * HalfGridBoundsX;
  HostParameters.WorldSize.s[1] = 2 * HalfGridBoundsY;

  HostParameters.Gas_Constant = GasConstant;
  HostParameters.KernelLength = InteractionRadius;
  HostParameters.H2 = pow(InteractionRadius, 2);
  HostParameters.Restitution = Restitution;
  HostParameters.Rest_Density = ReferenceDensity;
  HostParameters.Viscosity = Viscosity;
  HostParameters.Surface_Normal = 6.0f;
  HostParameters.Surface_Coefficient = 0.2;

  float H8 = pow(HostParameters.H2, 4);

  // Calculate number of particles from bounding box parameters
  PlaceParticles();

  // Constant density terms calculated with initial default mass (1) first (as they are used in normalizing mass)
  ConstantDensityKernelTerm = 4 * ((ParticleMass / PI) / H8);
  ConstantDensitySumTerm = 4 * ((ParticleMass / PI) / HostParameters.H2);

  NormalizeMass();
  // Assign mass after normalization
  HostParameters.Mass = ParticleMass;

  // load and build kernel program
  clObject.loadProgram(kernel_source);
  clObject.CreateKernel("SPH_kernel");

  void * hostParticlesRawPointer = static_cast<void *>(HostParticles.begin()._Ptr);
  Particle Sample;
  // Request buffer of particles on device memory
  clObject.requestCustomTypeBuffer(hostParticlesRawPointer, ParticleCount, CLBuffer::BufferTypes::READ_WRITE, Sample);

  // Set kernel arguments
  //InitKernel();
  // Write host buffer memory to device memory
  clObject.writeAllBuffers();
 // ZeroConnectThisTo(this->GetSpace(), "LogicUpdate", "OnLogicUpdate"); 
  
  Zilch::Console::Write("Particle Count ->");
  Zilch::Console::WriteLine(ParticleCount);
}

//***************************************************************************
void FluidSPHPlugin::OnLogicUpdate(ZeroEngine::UpdateEvent* event)
{
	for (int i = 0; i < NumberofSteps; ++i)
	{
		//ComputeAcceleration();
		//LeapfrogIntegrator(event->GetDt());
		//clObject.RunKernel(ParticleCount);
	}
	//UpdateArray();
}

void FluidSPHPlugin::ComputeInitalDensity()
{
	// Reset entire density array every time
	for (int i = 0; i < ParticleCount; ++i)
	{
		HostParticles[i].density = 0;
	}
	for (int i = 0; i < ParticleCount; ++i)
	{
		HostParticles[i].density += ConstantDensitySumTerm;
		for (int j = i + 1; j < ParticleCount; ++j)
		{
			// index values for ith element
			float iXPosition = HostParticles[i].position.s[0];
			float iYPosition = HostParticles[i].position.s[1];
			// index values for jth element
			float jXPosition = HostParticles[j].position.s[0];
			float jYPosition = HostParticles[j].position.s[1];

			// Finds change along x and y direction
			float deltaPositionX = iXPosition - jXPosition;
			float deltaPositionY = iYPosition - jYPosition;

			// creates comparision term to be used in kernel
			float r2 = (deltaPositionX * deltaPositionX) + (deltaPositionY * deltaPositionY);
			float z = HostParameters.H2 - r2;

			// If difference between kernel value and gradient values is significant 
			// (i.e. non -zero or a non-negative contribution on the ith particle from jth particle)
			// Use the calculated density to affect both particles
			if (z > 0)
			{
				float rho_ij = ConstantDensityKernelTerm * z * z * z;
				HostParticles[i].density += rho_ij;
				HostParticles[j].density += rho_ij;
			}
		}
	}
}

void FluidSPHPlugin::UpdateArray()
{
	for (int i = 0; i < ParticleCount; ++i)
	{
		Zilch::Real3 newPosition = Zilch::Real3(Position[(4 * i) + 0], Position[(4 * i) + 1], this->GetOwner()->GetTransform()->GetWorldTranslation().z);
		FluidParticlesArray[i].Dereference()->GetTransform()->SetWorldTranslation(newPosition);
	}
}

void FluidSPHPlugin::PlaceParticles()
{
	float interval = InteractionRadius / 1.3f;

	float xMin = this->GetOwner()->GetTransform()->GetWorldTranslation().x - HalfGridBoundsX;
	float xMax = this->GetOwner()->GetTransform()->GetWorldTranslation().x + HalfGridBoundsX;
	float yMin = this->GetOwner()->GetTransform()->GetWorldTranslation().y - HalfGridBoundsY;
	float yMax = this->GetOwner()->GetTransform()->GetWorldTranslation().y + HalfGridBoundsY;

	// Find out how many particles would fit in this region
	for (float x = xMin; x < xMax; x += interval)
	{
		for (float y = yMin; y < yMax; y += interval)
		{
			ParticleCount += BoxIndicator(x, y); // <------ Change this to get different shaped container for fluid volume
		}
	}

	// Initialize host particles
	HostParticles.resize(ParticleCount);

	int counter = 0;
	// Place the particles correspondingly and initialize their velocities
	for (float x = xMin; x < xMax; x += interval)
	{
		for (float y = yMin; y < yMax; y += interval)
		{
			if (BoxIndicator(x, y))       // <------ Change this to get different shaped container for fluid volume
			{
				HostParticles[counter].position.s[0] = x;
				HostParticles[counter].position.s[1] = y;
				++counter;
			}
		}
	}
	// Gets handle to particle archetype
	Zilch::HandleOf<ZeroEngine::Archetype> handle = ZeroEngine::Archetype::Find(this->FluidParticleArchetypeName.c_str());
	// Instantiates the particles and adds them to array
	for (int i = 0; i < ParticleCount; ++i)
	{
		FluidParticlesArray.push_back(this->GetSpace()->CreateAtPosition(handle, Zilch::Real3(HostParticles[i].position.s[0], HostParticles[i].position.s[1], this->GetOwner()->GetTransform()->GetWorldTranslation().z)));
	}
}

void FluidSPHPlugin::NormalizeMass()
{
	ComputeInitalDensity();
	float rho0 = HostParameters.Rest_Density;
	float rho2s = 0.0;
	float rhos = 0.0;
	for (int i = 0; i < ParticleCount; ++i)
	{
		rho2s += pow(HostParticles[i].density, 2);
		rhos += HostParticles[i].density;
	}
	ParticleMass *= ((rho0 * rhos) / rho2s);
}

void FluidSPHPlugin::InitKernel()
{
	// Delta-time
	clObject.SetKernelArgument(0, 0.016f);

	// Epsilon
	clObject.SetKernelArgument(1, 0.0001f);

	// Restitution
	clObject.SetKernelArgument(2, Restitution);

	// Constants
	clObject.SetKernelArgument(3, ConstantDensitySumTerm);
	clObject.SetKernelArgument(4, ConstantDensityKernelTerm);
	clObject.SetKernelArgument(5, HostParameters.H2);
	clObject.SetKernelArgument(6, ReferenceDensity);
	clObject.SetKernelArgument(7, InteractionRadius);
	// Boundary values
	float xMin = this->GetOwner()->GetTransform()->GetWorldTranslation().x - HalfGridBoundsX;
	float xMax = this->GetOwner()->GetTransform()->GetWorldTranslation().x + HalfGridBoundsX;
	float yMin = this->GetOwner()->GetTransform()->GetWorldTranslation().y - HalfGridBoundsY;
	float yMax = this->GetOwner()->GetTransform()->GetWorldTranslation().y + HalfGridBoundsY;
	clObject.SetKernelArgument(11, xMin);
	clObject.SetKernelArgument(12, xMax);
	clObject.SetKernelArgument(13, yMin);
	clObject.SetKernelArgument(14, yMax);

	// Buffer values
	clObject.SetKernelArgument(15, clObject.BufferValue(0));
	clObject.SetKernelArgument(16, clObject.BufferValue(1));
	clObject.SetKernelArgument(17, clObject.BufferValue(2));
	clObject.SetKernelArgument(18, clObject.BufferValue(3));
	clObject.SetKernelArgument(19, clObject.BufferValue(4));

	// Initialize local memory
	clObject.SetKernelArgument(20, cl::Local(ParticleCount));

	clObject.FinishCommandQueue();
}

int FluidSPHPlugin::BoxIndicator(float xVal, float yVal)
{
	return ((std::abs(xVal - this->GetOwner()->GetTransform()->GetWorldTranslation().x) < HalfGridBoundsX) && (std::abs(yVal - this->GetOwner()->GetTransform()->GetWorldTranslation().y) < HalfGridBoundsY));
}

int FluidSPHPlugin::CircleIndictor(float xVal, float yVal)
{
	return 0;
}

//***************************************************************************

//***************************************************************************
ZilchDefineType(FluidSPHPluginEvent, FluidSPHPluginLibrary, builder, type)
{
  // This is required for event binding
  ZilchBindDestructor();
  ZilchBindConstructor();
  ZilchBindFieldProperty(mLives);
}
