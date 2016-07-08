#pragma once
#include "FluidSPHPluginPrecompiled.hpp"

typedef struct Particle
{
	cl_float4 position;
	cl_float4 velocity_full;
	cl_float4 velocity_half;
	cl_float4 acceleration;
	cl_float density;
public:
	Particle()
	{
		for (int i = 0; i < 4; ++i)
		{
			position.s[i] = 0;
			velocity_full.s[i] = 0;
			velocity_half.s[i] = 0;
			acceleration.s[i] = 0;
		}
		density = 0;
	}
}Particle;


typedef struct SystemParameters
{	
	// System constant parameters
	cl_float DeltaTime;
	cl_float Epsilon;
	cl_float Restitution;
	cl_float2 Gravity;
	
	// Grid parameters
	cl_float2 WorldSize;
	cl_float CellSize;
	cl_uint2 GridSize;
	cl_uint Total_Cell_Count;
	
	// Particle parameters
	cl_float Rest_Density;
	cl_float Gas_Constant;
	cl_float Viscosity;
	/* ? */cl_float Surface_Normal;
	/* ? */cl_float Surface_Coefficient;
	cl_float Mass;
	cl_float Self_Density;

	// Kernel parameters
	cl_float KernelLength;
	cl_float H2;				// (KernelLength) ^ 2
	cl_float Poly6_Constant;
	cl_float Spiky_Constant;
	cl_float Viscosity_Constant;
	cl_float Poly6_Gradient_Constant;
	cl_float Poly6_Laplacian_Constant;
	/* ? */cl_float Self_Laplacian_Color;
}SystemParameters;