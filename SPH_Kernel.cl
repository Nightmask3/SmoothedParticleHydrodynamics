#define STRINGIFY(A) #A
std::string kernel_source = STRINGIFY(

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
    density[gid] = 0;
    density[gid] += ConstantDensitySumTerm;

    float4 thisPosition = position[gid];
    float densityTerm = 0.0;
    
	/*for (int i = 0; i < globalSize; ++i)
	{
		float4 otherPosition = position[i];

		float4 deltaPosition = thisPosition - otherPosition;

		float r2 = (deltaPosition.x * deltaPosition.x) + (deltaPosition.y * deltaPosition.y) + (deltaPosition.z * deltaPosition.z);
		float z = (H2 - r2) + eps;

		if (z > 0)
		{
			float rho_ij = ConstantDensityKernelTerm * z * z * z;
			density[gid] += rho_ij;
			density[i] += rho_ij;
		}
	}*/
	
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

            if(z > 0) 
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
    int tid = get_local_id(0);

    int globalSize = get_global_size(0);
    int localSize = get_local_size(0);
    int numTiles = globalSize/localSize;

    // Set acceleration parameters
    acceleration[gid].x = 0.0;
    acceleration[gid].y = -0.1;

    float4 thisPosition = position[gid];
    float4 thisVelocity = velocity_full[gid];
    float rhoi = density[gid];
    float accelerationTermX = 0.0;
    float accelerationTermY = 0.0;

	for (int i = 0; i < globalSize; ++i)
	{
		float4 otherPosition = position[i];

		float4 deltaPosition = thisPosition - otherPosition;

		float r2 = (deltaPosition.x * deltaPosition.x) + (deltaPosition.y * deltaPosition.y) + (deltaPosition.z * deltaPosition.z);
		float z = (H2 - r2) + eps;
		if (z > 0 && r2 != 0.0f) 
		{
			float rhoj = density[i];
			float q = sqrt(r2) / InteractionRadius;
			float u = 1 - q;
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
   /* for(int i = 0; i < numTiles; ++i)
    {
        for(int j = 0; j < localSize; ++j)
        {
            float4 otherPosition = pblock[j];
				 
            float4 deltaPosition = thisPosition - otherPosition;

            float r2 = (deltaPosition.x * deltaPosition.x) + (deltaPosition.y * deltaPosition.y) + (deltaPosition.z * deltaPosition.z);

            if(r2 < (H2 + eps))
            {
                float rhoj = density[j];
                float q = sqrt(r2) / InteractionRadius;
                float u = 1 - q;
                float w0 = C0 * (u / rhoi / rhoj);
                float wP = w0 * CP * (rhoi + rhoj - (2 * ReferenceDensity)) * (u / q);
                float wV = w0 * CV;

                float4 deltaVelocity = thisVelocity - velocity_full[j];

                accelerationTermX += (wP * deltaPosition.x) + (wV * deltaVelocity.x);
                accelerationTermY += (wP * deltaPosition.y) + (wV * deltaVelocity.y);
            }
        }
    }*/
	
}

void LeapfrogIntegrator(
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

void ReflectBoundaryConditions(
float xMin,
float xMax,
float yMin,
float yMax
)
{

}

void DampenReflections(
int globalIndex,
float barrier
)
{}

void kernel SPH_kernel(
float dt, float eps,
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

   // LeapfrogIntegrator(dt4, position, velocity_full, velocity_half, acceleration);
}
);