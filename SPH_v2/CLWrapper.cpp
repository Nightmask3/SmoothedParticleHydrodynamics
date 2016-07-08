#include "FluidSPHPluginPrecompiled.hpp"

CLWrapper::CLWrapper()
{
	// Initialize OpenCL object and context

	// Get all platforms that are OpenCL compatible
	err_ = cl::Platform::get(&allPlatforms_);
	Zilch::Console::Write("Collect all OpenCL platforms. Status : ");
	Zilch::Console::WriteLine(oclErrorString(err_));

	if (allPlatforms_.size() == 0)
	{
		Zilch::Console::WriteLine("No OpenCL compatible platforms found, check installation!");
		exit(1);
	}

	// Use default platform 
	platform_ = allPlatforms_[0];
	Zilch::Console::Write("Using platform: ");
	Zilch::Console::WriteLine(platform_.getInfo<CL_PLATFORM_NAME>().c_str());

	// Get all devices registered to platform
	err_ = platform_.getDevices(CL_DEVICE_TYPE_ALL, &allDevices_);
	Zilch::Console::Write("Collect all OpenCL devices registered to platform. Status : ");
	Zilch::Console::WriteLine(oclErrorString(err_));

	if (allDevices_.size() == 0) 
	{
		Zilch::Console::WriteLine(" No OpenCL compatible devices found, check installation!");
		exit(1);
	}

	// Use default devices
	device_ = allDevices_[0];

	// Create context
	context_ = cl::Context({ device_ });

	// Create command queue
	queue_ = cl::CommandQueue(context_, device_, 0, &err_);

	Zilch::Console::Write("Creating command queue. -> Status : ");
	Zilch::Console::WriteLine(oclErrorString(err_));

}

CLWrapper::~CLWrapper()
{
}

void CLWrapper::RunKernel(int GlobalWorkSize)
{
	//writeAllBuffers();

	// Changing local work-group size causes issues (Why?)
	err_ = queue_.enqueueNDRangeKernel(kernel_, cl::NullRange, cl::NDRange(GlobalWorkSize), cl::NullRange);
	
	//Zilch::Console::Write("Running kernel. -> Status :");
	//Zilch::Console::WriteLine(oclErrorString(err_));
	err_ = queue_.enqueueReadBuffer(allBuffers_[0].bufferValue, TRUE, 0, allBuffers_[0].bufferSize, allBuffers_[0].dataPointer);
	queue_.finish();
}

void CLWrapper::loadProgram(std::string kernel_source)
{
	int programLength;
	programLength = kernel_source.size();

	// Make program from source
	cl::Program::Sources source(1, std::make_pair(kernel_source.c_str(), programLength));
	program_ = cl::Program(context_, source);

	// Build program
	err_ = program_.build(allDevices_);
	Zilch::Console::Write("Building program... -> Status : ");
	Zilch::Console::WriteLine(oclErrorString(err_));

	if (err_ != CL_SUCCESS)
	{
		Zilch::Console::WriteLine("Building program failed!");
		//exit(1);
	}
}

void CLWrapper::readAllBuffers()
{
	for each (CLBuffer buf in allBuffers_)
	{
		err_ = queue_.enqueueReadBuffer(buf.bufferValue, CL_TRUE, 0, buf.bufferSize, buf.dataPointer, NULL, NULL);
		if (err_ != CL_SUCCESS)
		{
			Zilch::Console::Write("Reading buffers. -> Status : ");
			Zilch::Console::WriteLine(oclErrorString(err_));
		}
	}
}

void CLWrapper::writeAllBuffers()
{
	for each (CLBuffer buf in allBuffers_)
	{
		err_ = queue_.enqueueWriteBuffer(buf.bufferValue, CL_TRUE, 0, buf.bufferSize, buf.dataPointer, NULL, NULL);
		queue_.finish();
		if (err_ != CL_SUCCESS)
		{
			Zilch::Console::Write("Writing buffers. -> Status : ");
			Zilch::Console::WriteLine(oclErrorString(err_));
		}
	}
}

void CLWrapper::requestFloatBuffer(cl_float * ArrayPointer, int BufferSize, CLBuffer::BufferTypes BufferType)
{
	size_t arraySize = BufferSize * sizeof(cl_float);
	CLBuffer newBuffer;
	switch (BufferType)
	{
	case CLBuffer::BufferTypes::READ :
		newBuffer.bufferValue = cl::Buffer(context_, CL_MEM_READ_ONLY, arraySize, NULL, &err_);
		newBuffer.dataPointer = ArrayPointer;
		newBuffer.bufferSize = arraySize;
		break;
	case CLBuffer::BufferTypes::WRITE:
		newBuffer.bufferValue = cl::Buffer(context_, CL_MEM_WRITE_ONLY, arraySize, NULL, &err_);
		newBuffer.dataPointer = ArrayPointer;
		newBuffer.bufferSize = arraySize;
		break;
	case CLBuffer::BufferTypes::READ_WRITE:
		newBuffer.bufferValue = cl::Buffer(context_, CL_MEM_READ_WRITE, arraySize, NULL, &err_);
		newBuffer.dataPointer = ArrayPointer;
		newBuffer.bufferSize = arraySize;
		break;
	}

	allBuffers_.push_back(newBuffer);
}

void CLWrapper::requestFloat4Buffer(cl_float * ArrayPointer, int BufferSize, CLBuffer::BufferTypes BufferType)
{
	size_t arraySize = BufferSize * sizeof(cl_float4);
	CLBuffer newBuffer;
	switch (BufferType)
	{
	case CLBuffer::BufferTypes::READ:
		newBuffer.bufferValue = cl::Buffer(context_, CL_MEM_READ_ONLY, arraySize, NULL, &err_);
		newBuffer.dataPointer = ArrayPointer;
		newBuffer.bufferSize = arraySize;
		break;
	case CLBuffer::BufferTypes::WRITE:
		newBuffer.bufferValue = cl::Buffer(context_, CL_MEM_WRITE_ONLY, arraySize, NULL, &err_);
		newBuffer.dataPointer = ArrayPointer;
		newBuffer.bufferSize = arraySize;
		break;
	case CLBuffer::BufferTypes::READ_WRITE:
		newBuffer.bufferValue = cl::Buffer(context_, CL_MEM_READ_WRITE, arraySize, NULL, &err_);
		newBuffer.dataPointer = ArrayPointer;
		newBuffer.bufferSize = arraySize;
		break;
	}
	if (err_ != CL_SUCCESS)
	{
		Zilch::Console::WriteLine("Error in creating Buffer!");
	}
	allBuffers_.push_back(newBuffer);
}
