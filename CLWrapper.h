#pragma once
struct CLBuffer
{
	enum BufferTypes
	{
		READ,
		WRITE,
		READ_WRITE
	};
	cl::Buffer bufferValue;
	float bufferSize;
	float * dataPointer;
};

class CLWrapper
{
	CLWrapper();
	~CLWrapper();

	// RunKernel and InitKernel are pure virtual as they are dependent on implementation and must be defined by the client inherited class
	virtual void RunKernel() = 0;
	virtual void InitKernel() = 0;
	// loads, creates and builds kernel program from specified string
	void loadProgram(std::string kernel_source);
	// Bulk Buffer Read/Write methods
	void readAllBuffers();
	void writeAllBuffers();
	// Buffer request methods
	void requestFloatBuffer(float * ArrayPointer, int BufferSize, CLBuffer::BufferTypes BufferType);
	void requestFloat4Buffer(float * ArrayPointer, int BufferSize, CLBuffer::BufferTypes BufferType);
private:
	// Handles to the various OpenCL objects required for execution
	cl::Kernel _kernel;
	cl::Platform _platform;
	cl::Device _device;
	cl::Context _context;
	cl::CommandQueue _queue;
	cl::Program _program;
	std::vector<cl::Platform> _allPlatforms;
	std::vector<cl::Device> _allDevices;

	std::vector<CLBuffer> _allBuffers;

	// debugging variables
	cl_int _err;
};