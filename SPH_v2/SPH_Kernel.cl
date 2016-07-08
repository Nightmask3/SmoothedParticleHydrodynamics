#define STRINGIFY(A) #A
std::string kernel_source = STRINGIFY(

typedef struct Particle
{
	float4 position;
	float4 velocity_full;
	float4 velocity_half;
	float4 acceleration;
	float density;
}Particle;

typedef struct SystemParameters
{
	// System constant parameters
	float DeltaTime;
	float Epsilon;
	float Restitution;
	float2 Gravity;

	// Grid parameters
	float2 WorldSize;
	float CellSize;
	uint2 GridSize;
	uint Total_Cell_Count;

	// Particle parameters
	float Rest_Density;
	float Gas_Constant;
	float Viscosity;
	float Surface_Normal;
	float Surface_Coefficient;
	float Mass;
	float Self_Density;

	// Kernel parameters
	float KernelLength;
	float H2;				// (KernelLength) ^ 2
	float Poly6_Constant;
	float Spiky_Constant;
	float Viscosity_Constant;
	float Poly6_Gradient_Constant;
	float Poly6_Laplacian_Constant;
	// No clue what this one's for
	float Self_Laplacian_Color;

}SystemParameters;

void Density_Calculation(
float ConstantDensitySumTerm,
float ConstantDensityKernelTerm,
float eps,
float H2, 
global float4* position,
global float* density,
local float4* pblock
)
{
    // Id of this work-item in global index space
    int gid = get_global_id(0);
    // Id of this work-item within it's work group
    int tid = get_local_id(0);

    int globalSize = get_global_size(0);
    int localSize = get_local_size(0);
    int numTiles = globalSize/localSize;

    // Zero out the density term of this work-item
    density[gid] = 0.0f;
    density[gid] += ConstantDensitySumTerm;

    float4 thisPosition = position[gid];
    float densityTerm = 0.0f;
	
	// Outer loop iterates over all the work-group blocks
    for(int i = 0; i < numTiles; ++i)
    {
        // Cache the particle position within the work-group
        pblock[tid] = position[(i * localSize) + tid];

        // synchronize to make sure data is available for processing
        barrier(CLK_LOCAL_MEM_FENCE);

        // Inner loop iterates over the work-items in each work-group, after all positions have been cached
        for(int j = 0; j < localSize; ++j)
        {
            float4 otherPosition = pblock[j];

            float4 deltaPosition = thisPosition - otherPosition;

            float r2 = (deltaPosition.x * deltaPosition.x) + (deltaPosition.y * deltaPosition.y) + (deltaPosition.z * deltaPosition.z);
            float z = (H2 - r2) + eps;

            if(z > 0.0f) 
            {
                float rho_ij = ConstantDensityKernelTerm * z * z * z;
                densityTerm += rho_ij;
            }
        }

        // Synchronize so that next tile can be loaded
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    density[gid] += densityTerm;
}

void Acceleration_Calculation(
float eps, 
float ConstantDensitySumTerm, float ConstantDensityKernelTerm,
float H2, float ReferenceDensity, float InteractionRadius,
float C0, float CP, float CV,
global float4* position,
global float4* velocity_full,
global float4* acceleration,
global float* density,
local float4* pblock
)
{
    Density_Calculation(ConstantDensitySumTerm, ConstantDensityKernelTerm, eps, H2, position, density, pblock);

    // Id of this work-item in global index space
    int gid = get_global_id(0);
    // Id of this work-item within it's work group
   // int tid = get_local_id(0);

    int globalSize = get_global_size(0);
    //int localSize = get_local_size(0);
   // int numTiles = globalSize/localSize;

    // Set acceleration parameters
    acceleration[gid].x = 0.0f;
    acceleration[gid].y = -0.08f;

    float4 thisPosition = position[gid];
    float4 thisVelocity = velocity_full[gid];
    float rhoi = density[gid];
    float accelerationTermX = 0.0f;
    float accelerationTermY = 0.0f;

	for (int i = 0; i < globalSize; ++i)
	{
		float4 deltaPosition = thisPosition - position[i];

		float r2 = (deltaPosition.x * deltaPosition.x) + (deltaPosition.y * deltaPosition.y) + (deltaPosition.z * deltaPosition.z);
		float z = (H2 - r2) + eps;
		if (z > 0.0f && r2 != 0.0f) 
		{
			float rhoj = density[i];
			float q = sqrt(r2) / InteractionRadius;
			float u = 1.0f - q;
			float w0 = C0 * (u / rhoi / rhoj);
			float wP = w0 * CP * (rhoi + rhoj - (2 * ReferenceDensity)) * (u / q);
			float wV = w0 * CV;

			float4 deltaVelocity = thisVelocity - velocity_full[i];

			accelerationTermX = (wP * deltaPosition.x) + (wV * deltaVelocity.x);
			accelerationTermY = (wP * deltaPosition.y) + (wV * deltaVelocity.y);

			acceleration[gid].x += accelerationTermX;
			acceleration[gid].y += accelerationTermY;

			acceleration[i].x -= accelerationTermX;
			acceleration[i].y -= accelerationTermY;
		}
	}

  /*  for(int i = 0; i < numTiles; ++i)
    {
        for(int j = 0; j < localSize; ++j)
        {
            float4 otherPosition = pblock[j];
				 
            float4 deltaPosition = thisPosition - otherPosition;

            float r2 = (deltaPosition.x * deltaPosition.x) + (deltaPosition.y * deltaPosition.y) + (deltaPosition.z * deltaPosition.z);

            if(r2 < (H2 + eps) && r2 != 0.0f)
            {
                float rhoj = density[i * localSize + j];
                float q = sqrt(r2) / InteractionRadius;
                float u = 1.0f - q;
                float w0 = C0 * (u / rhoi / rhoj);
                float wP = w0 * CP * (rhoi + rhoj - (2 * ReferenceDensity)) * (u / q);
                float wV = w0 * CV;

                float4 deltaVelocity = thisVelocity - velocity_full[i * localSize + j];

                accelerationTermX = (wP * deltaPosition.x) + (wV * deltaVelocity.x);
                accelerationTermY = (wP * deltaPosition.y) + (wV * deltaVelocity.y);

				acceleration[gid].x += accelerationTermX;
				acceleration[gid].y += accelerationTermY;

				acceleration[i * localSize + j].x -= accelerationTermX;
				acceleration[i * localSize + j].y -= accelerationTermY;
            }
        }
    }*/
	
}

inline void LeapfrogIntegrator(
const float4 dt, 
global float4* position,
global float4* velocity_full,
global float4* velocity_half, 
global float4* acceleration
)
{
    // Id of this work-item in global index space
    int gid = get_global_id(0);
    velocity_half[gid] = velocity_full[gid] + (acceleration[gid] * (dt/2));
    velocity_full[gid] += acceleration[gid] * dt;
    position[gid] += velocity_full[gid] * dt;
}

inline void DampenReflectionsX(
global float4* position,
global float4* velocity_full,
global float4* velocity_half,
float Restitution,
float barrier
)
{
	// Id of this work-item in global index space
	int gid = get_global_id(0);

	// Ignore degenerate cases
	//if (velocity_full[gid].x == 0.0f)
	//	return;

	// find time since crossing barrier
	float tbounce = (position[gid].x - barrier) / velocity_full[gid].x;
	
	// Scale back the distance traveled based on time from collision
	position[gid].x -= velocity_full[gid].x * (1 - Restitution) * tbounce;
	position[gid].y -= velocity_full[gid].y * (1 - Restitution) * tbounce;

	// Reflect the position and velocity according to which axis the penetration occured on
	position[gid].x = 2 * barrier - position[gid].x;
	velocity_full[gid].x = -velocity_full[gid].x;
	velocity_half[gid].x = -velocity_half[gid].x;

	// Damp the velocities
	velocity_full[gid] *= Restitution;
	velocity_half[gid] *= Restitution;
}

inline void DampenReflectionsY(
	global float4* position,
	global float4* velocity_full,
	global float4* velocity_half,
	float Restitution,
	float barrier
)
{
	// Id of this work-item in global index space
	int gid = get_global_id(0);

	// Ignore degenerate cases
	//if (velocity_full[gid].y == 0.0f)
	//	return;

	// find time since crossing barrier
	float tbounce = (position[gid].y - barrier) / velocity_full[gid].y;

	// Scale back the distance traveled based on time from collision
	position[gid].x -= velocity_full[gid].x * (1 - Restitution) * tbounce;
	position[gid].y -= velocity_full[gid].y * (1 - Restitution) * tbounce;

	// Reflect the position and velocity according to which axis the penetration occured on
	position[gid].y = 2 * barrier - position[gid].y;
	velocity_full[gid].y = -velocity_full[gid].y;
	velocity_half[gid].y = -velocity_half[gid].y;
	
	// Damp the velocities
	velocity_full[gid] *= Restitution;
	velocity_half[gid] *= Restitution;
}



void ReflectBoundaryConditions(
	global float4* position,
	global float4* velocity_full,
	global float4* velocity_half,
	float Restitution,
	float xMin,
	float xMax,
	float yMin,
	float yMax
)
{
	// Id of this work-item in global index space
	int gid = get_global_id(0);

	if (position[gid].x < xMin)
		DampenReflectionsX(position, velocity_full, velocity_half, Restitution, xMin);
	else if (position[gid].x > xMax)
		DampenReflectionsX(position, velocity_full, velocity_half, Restitution, xMax);

	if (position[gid].y < yMin)
		DampenReflectionsY(position, velocity_full, velocity_half, Restitution, yMin);
	else if (position[gid].y > yMax)
		DampenReflectionsY(position, velocity_full, velocity_half, Restitution, yMax);

}

void kernel SPH_kernel(
float dt, float eps, float Restitution,
float ConstantDensitySumTerm, float ConstantDensityKernelTerm,
float H2, float ReferenceDensity, float InteractionRadius,
float C0, float CP, float CV,
float xMin, float xMax, float yMin, float yMax,
global float4* position,
global float4* velocity_half,
global float4* velocity_full,
global float4* acceleration,
global float* density,
local float4* pblock
)
{
    const float4 dt4 = (float4)(dt,dt,dt,0.0f);

    Acceleration_Calculation(eps, ConstantDensitySumTerm, ConstantDensityKernelTerm, H2, ReferenceDensity, InteractionRadius, C0, CP, CV,
                             position, velocity_full, acceleration, density, pblock);

   LeapfrogIntegrator(dt4, position, velocity_full, velocity_half, acceleration);

   ReflectBoundaryConditions(position, velocity_full, velocity_half, Restitution, xMin, xMax, yMin, yMax);
}
);