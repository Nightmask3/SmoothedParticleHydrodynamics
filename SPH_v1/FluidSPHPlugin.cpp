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
  ZilchBindFieldProperty(GravityMagnitude);
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
  // Close output log stream
 // if(outputStream != NULL)
//	fclose(outputStream);
  //Zilch::Console::WriteLine("FluidSPHPlugin::~FluidSPHPlugin (Destructor)");
}

//***************************************************************************
void FluidSPHPlugin::Initialize(ZeroEngine::CogInitializer* initializer)
{
  //Zilch::Console::WriteLine("FluidSPHPlugin::Initialize");
  
  H2 = (InteractionRadius * InteractionRadius);
  H8 = (H2 * H2) * (H2 * H2);

  // Calculate number of particles from bounding box parameters
  PlaceParticles();
  
  VelocityHalfStep.resize(4 * ParticleCount);
  Acceleration.resize(4 * ParticleCount);
  Density.resize(ParticleCount);

  // Constant density terms calculated with initial default mass (1) first (as they are used in normalizing mass)
  ConstantDensityKernelTerm = 4 * ((ParticleMass / PI) / H8);
  ConstantDensitySumTerm = 4 * ((ParticleMass / PI) / H2);

  NormalizeMass();

  // constant terms calculated once again with new normalized mass
  ConstantDensityKernelTerm = 4 * ((ParticleMass / PI) / H8);
  ConstantDensitySumTerm = 4 * ((ParticleMass / PI) / H2);

  C0 = ((ParticleMass / PI) / (H2 * H2));
  CP = 15 * BulkModulus;
  CV = -40 * Viscosity;

  // load and build kernel program
  clObject.loadProgram(kernel_source);
  clObject.CreateKernel("SPH_kernel");

  // Request buffers
  int sizePos = sizeof(Position);
  clObject.requestFloat4Buffer(&Position[0], ParticleCount, CLBuffer::BufferTypes::READ_WRITE);
  clObject.requestFloat4Buffer(&VelocityHalfStep[0], ParticleCount, CLBuffer::BufferTypes::READ_WRITE);
  clObject.requestFloat4Buffer(&VelocityFullStep[0], ParticleCount, CLBuffer::BufferTypes::READ_WRITE);
  clObject.requestFloat4Buffer(&Acceleration[0], ParticleCount, CLBuffer::BufferTypes::READ_WRITE);
  clObject.requestFloatBuffer(&Density[0], ParticleCount, CLBuffer::BufferTypes::READ_WRITE);

  // Set kernel arguments
  InitKernel();
  // Write host buffer memory to device memory
  clObject.writeAllBuffers();
  ZeroConnectThisTo(this->GetSpace(), "LogicUpdate", "OnLogicUpdate"); 
  
  Zilch::Console::Write("Particle Count ->");
  Zilch::Console::WriteLine(ParticleCount);
  // Redirect stdout to the output log file
  // File pointer to the output log for OpenCL
  //if ((outputStream = freopen("C:\\Users\\sainarayan.n\\Documents\\ZeroProjects\\BackupSPH\\FluidSPH\\Content\\FluidSPHPlugin\\OpenCLOutputLog.txt", "w", stdout)) == NULL)
  //	  exit(-1);
}

//***************************************************************************
void FluidSPHPlugin::OnLogicUpdate(ZeroEngine::UpdateEvent* event)
{
	for (int i = 0; i < NumberofSteps; ++i)
	{
		//ComputeAcceleration();
		//LeapfrogIntegrator(event->GetDt());
		clObject.RunKernel(ParticleCount);
	}
	//ReflectBoundaryConditions();
	UpdateArray();
}

void FluidSPHPlugin::ComputeDensity()
{
	// Reset entire density array every time
	std::fill(Density.begin(), Density.end(), 0.0f);

	for (int i = 0; i < ParticleCount; ++i)
	{
		Density[i] += ConstantDensitySumTerm;
		for (int j = i + 1; j < ParticleCount; ++j)
		{
			// index values for ith element
			float iXPosition = Position[(4 * i) + 0];
			float iYPosition = Position[(4 * i) + 1];
			// index values for jth element
			float jXPosition = Position[(4 * j) + 0];
			float jYPosition = Position[(4 * j) + 1];

			// Finds change along x and y direction
			float deltaPositionX = iXPosition - jXPosition;
			float deltaPositionY = iYPosition - jYPosition;

			// creates comparision term to be used in kernel
			float r2 = (deltaPositionX * deltaPositionX) + (deltaPositionY * deltaPositionY);
			float z = H2 - r2;

			// If difference between kernel value and gradient values is significant 
			// (i.e. non -zero or a non-negative contribution on the ith particle from jth particle)
			// Use the calculated density to affect both particles
			if (z > 0)
			{
				float rho_ij = ConstantDensityKernelTerm * z * z * z;
				Density[i] += rho_ij;
				Density[j] += rho_ij;
			}
		}
	}
}

void FluidSPHPlugin::ComputeAcceleration()
{
	// Compute density
	ComputeDensity();
	// Set accleration parameter due to gravity
	for (int i = 0; i < ParticleCount; ++i)
	{
		Acceleration[(4 * i) + 0] = 0;
		Acceleration[(4 * i) + 1] = -GravityMagnitude;
	}

	// Compute interaction forces
	for (int i = 0; i < ParticleCount; ++i)
	{
		float rhoi = Density[i];
		for (int j = i + 1; j < ParticleCount; ++j)
		{
			// index values for ith element
			float iXPosition = Position[(4 * i) + 0];
			float iYPosition = Position[(4 * i) + 1];

			// index values for jth element
			float jXPosition = Position[(4 * j) + 0];
			float jYPosition = Position[(4 * j) + 1];

			// Finds gradient of position along x and y direction
			float deltaPositionX = iXPosition - jXPosition;
			float deltaPositionY = iYPosition - jYPosition;

			// creates comparision term to be used in kernel
			float r2 = (deltaPositionX * deltaPositionX) + (deltaPositionY * deltaPositionY);
			// if caluculated value is within the squared interaction radius 
			// (i.e. non -zero or a non-negative contribution on the ith particle from jth particle)
			// calculate the effect of these particles on each others accleration
			if (r2 < H2)
			{
				float rhoj = Density[j];
				float q = std::sqrt(r2) / InteractionRadius;
				float u = 1 - q;
				float w0 = C0 * (u / rhoi / rhoj);
				float wP = w0 * CP * (rhoi + rhoj - (2 * ReferenceDensity)) * (u / q);
				float wV = w0 * CV;
				float deltaVelocityX = VelocityFullStep[(4 * i) + 0] - VelocityFullStep[(4 * j) + 0];
				float deltaVelocityY = VelocityFullStep[(4 * i) + 1] - VelocityFullStep[(4 * j) + 1];

				// ith particle
				Acceleration[(4 * i) + 0] += (wP * deltaPositionX) + (wV * deltaVelocityX);
				Acceleration[(4 * i) + 1] += (wP * deltaPositionY) + (wV * deltaVelocityY);
				// jth particle
				Acceleration[(4 * j) + 0] -= (wP * deltaPositionX) + (wV * deltaVelocityX);
				Acceleration[(4 * j) + 1] -= (wP * deltaPositionY) + (wV * deltaVelocityY);
			}
		}
	}
}

void FluidSPHPlugin::DampenReflections(int which, float barrier)
{
	// Ignore degenerate cases
	if (VelocityFullStep[which] == 0.0)
		return;
	// Scale back the distance traveled based on time from collision
	float tbounce = ((Position[which] - barrier) / VelocityFullStep[which]);
	int indexX = 0;
	int indexY = 0;

	// If even, then its an x dimension value
	if (which % 2 == 0)
	{
		indexX = which;
		indexY = indexX + 1;

	}
	// If odd, then its a y dimension value
	else
	{
		indexY = which;
		indexX = indexY - 1;
	}

	// Scale back the distance traveled based on time from collision
	Position[indexX] -= VelocityFullStep[indexX] * (1 - Restitution) * tbounce;
	Position[indexY] -= VelocityFullStep[indexY] * (1 - Restitution) * tbounce;

	// Reflect the position and velocity
	Position[which] = 2 * barrier - Position[which];
	VelocityFullStep[which] = -VelocityFullStep[which];
	VelocityHalfStep[which] = -VelocityHalfStep[which];

	// Damp the velocities
	VelocityFullStep[indexX] *= Restitution; 
	VelocityFullStep[indexY] *= Restitution; 
	
	VelocityHalfStep[indexX] *= Restitution;
	VelocityHalfStep[indexY] *= Restitution;
}

void FluidSPHPlugin::ReflectBoundaryConditions()
{
	// Grid boundary values 
	float xMin = this->GetOwner()->GetTransform()->GetWorldTranslation().x - HalfGridBoundsX;
	float xMax = this->GetOwner()->GetTransform()->GetWorldTranslation().x + HalfGridBoundsX;
	float yMin = this->GetOwner()->GetTransform()->GetWorldTranslation().y - HalfGridBoundsY;
	float yMax = this->GetOwner()->GetTransform()->GetWorldTranslation().y + HalfGridBoundsY;
	for (int i = 0; i < 4 * ParticleCount; i += 4)
	{
		if (Position[i + 1] < yMin)
		{
			DampenReflections(i + 1, yMin);
		}

		if (Position[i + 1] > yMax)
		{
			DampenReflections(i + 1, yMax);
		}

		if (Position[i] < xMin)
		{
			DampenReflections(i, xMin);
		}

		if (Position[i] > xMax)
		{
			DampenReflections(i, xMax);
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

void FluidSPHPlugin::LeapfrogIntegrator(float deltaTime)
{
	for (int i = 0; i < 4 * ParticleCount; ++i)
	{
		VelocityHalfStep[i] = VelocityFullStep[i] + (Acceleration[i] * deltaTime / 2);
	}
	for (int i = 0; i < 4 * ParticleCount; ++i)
	{
		VelocityFullStep[i] += Acceleration[i] * deltaTime;
	}
	for (int i = 0; i < 4 * ParticleCount; ++i)
	{
		Position[i] += VelocityHalfStep[i] * deltaTime;
	}
	ReflectBoundaryConditions();
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

	Position.resize(4 * ParticleCount);
	VelocityFullStep.resize(4 * ParticleCount);
	int counter = 0;
	// Place the particles correspondingly and initialize their velocities
	for (float x = xMin; x < xMax; x += interval)
	{
		for (float y = yMin; y < yMax; y += interval)
		{
			if (BoxIndicator(x, y))       // <------ Change this to get different shaped container for fluid volume
			{
				Position[(4 * counter) + 0] = x;
				Position[(4 * counter) + 1] = y;
				VelocityFullStep[(4 * counter) + 0] = 0.0f;
				VelocityFullStep[(4 * counter) + 1] = 0.0f;
				++counter;
			}
		}
	}
	// Gets handle to particle archetype
	Zilch::HandleOf<ZeroEngine::Archetype> handle = ZeroEngine::Archetype::Find(this->FluidParticleArchetypeName.c_str());
	// Instantiates the particles and adds them to array
	for (int i = 0; i < ParticleCount; ++i)
	{
		FluidParticlesArray.push_back(this->GetSpace()->CreateAtPosition(handle, Zilch::Real3(Position[(4 * i) + 0], Position[(4 * i) + 1], this->GetOwner()->GetTransform()->GetWorldTranslation().z)));
	}
}

void FluidSPHPlugin::NormalizeMass()
{
	ComputeDensity();
	float rho0 = ReferenceDensity;
	float rho2s = 0.0;
	float rhos = 0.0;
	for (int i = 0; i < ParticleCount; ++i)
	{
		rho2s += Density[i] * Density[i];
		rhos += Density[i];
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
	clObject.SetKernelArgument(5, H2);
	clObject.SetKernelArgument(6, ReferenceDensity);
	clObject.SetKernelArgument(7, InteractionRadius);
	clObject.SetKernelArgument(8, C0);
	clObject.SetKernelArgument(9, CP);
	clObject.SetKernelArgument(10, CV);

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
