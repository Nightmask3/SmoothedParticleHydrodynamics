#pragma once
// TODO: Figure out a way to have typed dataPointer for more safety
// Maybe with child classes?
struct CLBuffer
{
	enum BufferTypes
	{
		READ,
		WRITE,
		READ_WRITE
	};
	cl::Buffer bufferValue;
	size_t bufferSize;
	void * dataPointer;
};

class CLWrapper
{
public:
	CLWrapper();
	~CLWrapper();

	// Getters/Setters
	inline cl::Buffer & BufferValue(int index) { return allBuffers_[index].bufferValue; }
	inline void * DataPointer(int index) { return allBuffers_[index].dataPointer; }
	// Utilities
	inline void FinishCommandQueue() { queue_.finish(); }
	inline void CreateKernel(std::string kernelName) { kernel_ = cl::Kernel(program_, kernelName.c_str(), &err_); }

	// Writes buffers, enqueues kernel execution, reads buffers
	void RunKernel(int GlobalWorkSize);
	// Allows for external setting of the arguments for classes without access to the kernel
	template <typename T>
	void SetKernelArgument(int ArgIndex, T Argument);
	// loads, creates and builds kernel program from specified string
	void loadProgram(std::string kernel_source);
	// Bulk Buffer Read/Write methods
	void readAllBuffers();
	void writeAllBuffers();
	// Buffer request methods
	void requestFloatBuffer(cl_float * ArrayPointer, int BufferSize, CLBuffer::BufferTypes BufferType);
	void requestFloat4Buffer(cl_float * ArrayPointer, int BufferSize, CLBuffer::BufferTypes BufferType);
	template<typename T>
	void requestCustomTypeBuffer(void * ArrayPointer, int BufferSize, CLBuffer::BufferTypes BufferType, T SampleObject);
protected:
	// Handles to the various OpenCL objects required for execution
	cl::Kernel kernel_;
	cl::Platform platform_;
	cl::Device device_;
	cl::Context context_;
	cl::CommandQueue queue_;
	cl::Program program_;
	std::vector<cl::Platform> allPlatforms_;
	std::vector<cl::Device> allDevices_;
	std::vector<CLBuffer> allBuffers_;

	// debugging variables
	cl_int err_;
};

template<typename T>
inline void CLWrapper::SetKernelArgument(int ArgIndex, T Argument)
{
	err_ = kernel_.setArg(ArgIndex, Argument);
	Zilch::Console::Write("Setting argument :");
	Zilch::Console::Write(ArgIndex);
	Zilch::Console::Write(" -> Status : ");
	Zilch::Console::WriteLine(oclErrorString(err_));
}

template<typename T>
void CLWrapper::requestCustomTypeBuffer(void * ArrayPointer, int BufferSize, CLBuffer::BufferTypes BufferType, T SampleObject)
{
	size_t arraySize = BufferSize * sizeof(T);
	T * dataPointer = static_cast<T *>(ArrayPointer);
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
}