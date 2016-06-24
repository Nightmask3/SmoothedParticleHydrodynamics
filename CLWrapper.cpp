#include "FluidSPHPluginPrecompiled.hpp"

CLWrapper::CLWrapper()
{
	// Initialize OpenCL object and context

	// Get all platforms that are OpenCL compatible
	_err = cl::Platform::get(&_allPlatforms);
	Zilch::Console::Write("Collect all OpenCL platforms. Status : ");
	Zilch::Console::WriteLine(oclErrorString(_err));

	if (_allPlatforms.size() == 0)
	{
		Zilch::Console::WriteLine("No OpenCL compatible platforms found, check installation!");
		exit(1);
	}

	// Use default platform 
	_platform = _allPlatforms[0];
	Zilch::Console::Write("Using platform: ");
	Zilch::Console::WriteLine(_platform.getInfo<CL_PLATFORM_NAME>().c_str());

	// Get all devices registered to platform
	_err = _platform.getDevices(CL_DEVICE_TYPE_ALL, &_allDevices);
	Zilch::Console::Write("Collect all OpenCL devices registered to platform. Status : ");
	Zilch::Console::WriteLine(oclErrorString(_err));

	if (_allDevices.size() == 0) 
	{
		Zilch::Console::WriteLine(" No OpenCL compatible devices found, check installation!");
		exit(1);
	}

	// Use default devices
	_device = _allDevices[0];

	// Create context
	_context = cl::Context({ _device });

	// Create command queue
	_queue = cl::CommandQueue(_context, _device, 0, &_err);

	Zilch::Console::Write("Creating command queue. Status : ");
	Zilch::Console::WriteLine(oclErrorString(_err));

}

CLWrapper::~CLWrapper()
{
}

void CLWrapper::loadProgram(std::string kernel_source)
{
	int programLength;
	programLength = kernel_source.size();

	// Make program from source
	cl::Program::Sources source(1, std::make_pair(kernel_source.c_str(), programLength));
	_program = cl::Program(_context, source);

	// Build program
	_err = _program.build(_allDevices);
	Zilch::Console::Write("Building program. Status : ");
	Zilch::Console::WriteLine(oclErrorString(_err));

	if (_err != CL_SUCCESS)
	{
		Zilch::Console::WriteLine("Building program failed!");
		exit(1);
	}
}

void CLWrapper::readAllBuffers()
{
	for each (CLBuffer buf in _allBuffers)
	{
		_err = _queue.enqueueReadBuffer(buf.bufferValue, CL_TRUE, 0, buf.bufferSize, buf.dataPointer, NULL, NULL);
		if (_err != CL_SUCCESS)
		{
			Zilch::Console::WriteLine("Error in reading buffers!");
			exit(1);
		}
	}
}

void CLWrapper::writeAllBuffers()
{
	for each (CLBuffer buf in _allBuffers)
	{
		_err = _queue.enqueueWriteBuffer(buf.bufferValue, CL_TRUE, 0, buf.bufferSize, buf.dataPointer, NULL, NULL);
		if (_err != CL_SUCCESS)
		{
			Zilch::Console::WriteLine("Error in writing buffers!");
			exit(1);
		}
	}
}

void CLWrapper::requestFloatBuffer(float * ArrayPointer, int BufferSize, CLBuffer::BufferTypes BufferType)
{
	float arraySize = BufferSize * sizeof(cl_float);
	CLBuffer newBuffer;
	switch (BufferType)
	{
	case CLBuffer::BufferTypes::READ :
		newBuffer.bufferValue = cl::Buffer(_context, CL_MEM_READ_ONLY, arraySize, NULL, &_err);
		newBuffer.dataPointer = ArrayPointer;
		newBuffer.bufferSize = arraySize;
		break;
	case CLBuffer::BufferTypes::WRITE:
		newBuffer.bufferValue = cl::Buffer(_context, CL_MEM_WRITE_ONLY, arraySize, NULL, &_err);
		newBuffer.dataPointer = ArrayPointer;
		newBuffer.bufferSize = arraySize;
		break;
	case CLBuffer::BufferTypes::READ_WRITE:
		newBuffer.bufferValue = cl::Buffer(_context, CL_MEM_READ_WRITE, arraySize, NULL, &_err);
		newBuffer.dataPointer = ArrayPointer;
		newBuffer.bufferSize = arraySize;
		break;
	}
	_allBuffers.push_back(newBuffer);
}

void CLWrapper::requestFloat4Buffer(float * ArrayPointer, int BufferSize, CLBuffer::BufferTypes BufferType)
{
	float arraySize = BufferSize * sizeof(cl_float4);
	CLBuffer newBuffer;
	switch (BufferType)
	{
	case CLBuffer::BufferTypes::READ:
		newBuffer.bufferValue = cl::Buffer(_context, CL_MEM_READ_ONLY, arraySize, NULL, &_err);
		newBuffer.dataPointer = ArrayPointer;
		newBuffer.bufferSize = arraySize;
		break;
	case CLBuffer::BufferTypes::WRITE:
		newBuffer.bufferValue = cl::Buffer(_context, CL_MEM_WRITE_ONLY, arraySize, NULL, &_err);
		newBuffer.dataPointer = ArrayPointer;
		newBuffer.bufferSize = arraySize;
		break;
	case CLBuffer::BufferTypes::READ_WRITE:
		newBuffer.bufferValue = cl::Buffer(_context, CL_MEM_READ_WRITE, arraySize, NULL, &_err);
		newBuffer.dataPointer = ArrayPointer;
		newBuffer.bufferSize = arraySize;
		break;
	}
	_allBuffers.push_back(newBuffer);
}
