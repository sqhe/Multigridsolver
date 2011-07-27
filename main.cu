#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cutil.h>

#define BLOCK_NUM 32
#define THREAD_NUM 256

bool InitCUDA(void)
{
	int count = 0;
	int i = 0;

	cudaGetDeviceCount(&count);
	if(count == 0) {
		fprintf(stderr, "There is no device.\n");
		return false;
	}

	for(i = 0; i < count; i++) {
		cudaDeviceProp prop;
		if(cudaGetDeviceProperties(&prop, i) == cudaSuccess) {
			if(prop.major >= 1) {
				break;
			}
		}
	}
	if(i == count) {
		fprintf(stderr, "There is no device supporting CUDA.\n");
		return false;
	}
	cudaSetDevice(i);

	printf("CUDA initialized.\n");
	return true;
}

__global__ static void HelloCUDA(char* result, int num)
{
	extern __shared__ int shared[];
	const int tid = threadIdx.x;
	const int bid = blockIdx.x;
	
	char p_HelloCUDA[] = "Hello CUDA!";
	for(int i = bid*THREAD_NUM + tid; i < num; i+=BLOCK_NUM*THREAD_NUM) {
		result[i] = p_HelloCUDA[i];
	}
}

int main(int argc,char *argv[])
{
	if (!InitCUDA())
	{
		return 0;
	}
	
	char *device_result = 0;
	char host_result[12] = {0};
	
	CUDA_SAFE_CALL(cudaMalloc((void **) &device_result,sizeof(char) *11));
	
	unsigned int timer = 0;
	CUT_SAFE_CALL(cutCreateTimer(&timer));
	CUT_SAFE_CALL(cutStartTimer(timer));
	
	HelloCUDA<<<BLOCK_NUM,THREAD_NUM,0>>>(device_result,11);
	CUT_CHECK_ERROR("Kernel execution failed\n");
	
	CUDA_SAFE_CALL( cudaThreadSynchronize() );
	CUT_SAFE_CALL( cutStopTimer( timer));
	printf("Processing time: %f (ms)\n", cutGetTimerValue( timer));
	CUT_SAFE_CALL( cutDeleteTimer( timer));

	CUDA_SAFE_CALL( cudaMemcpy(&host_result, device_result, sizeof(char) * 11, cudaMemcpyDeviceToHost));
	printf("%s\n", host_result);

	CUDA_SAFE_CALL( cudaFree(device_result));
	CUT_EXIT(argc, argv);	
	
	return 0;
}