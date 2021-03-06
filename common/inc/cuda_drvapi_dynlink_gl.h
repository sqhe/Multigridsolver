/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and 
 * proprietary rights in and to this software and related documentation. 
 * Any use, reproduction, disclosure, or distribution of this software 
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA) 
 * associated with this source code for terms and conditions that govern 
 * your use of this NVIDIA software.
 * 
 */

#ifndef __cuda_drvapi_dynlink_cuda_gl_h__
#define __cuda_drvapi_dynlink_cuda_gl_h__

#ifdef CUDA_INIT_OPENGL

	#ifdef _WIN32
	#  define WINDOWS_LEAN_AND_MEAN
	#  define NOMINMAX
	#  include <windows.h>
	#endif

	// includes, system
	#include <stdlib.h>
	#include <stdio.h>
	#include <string.h>
	#include <math.h>

	// includes, GL
	#include <GL/glew.h>

	#if defined (__APPLE__) || defined(MACOSX)
	#include <GLUT/glut.h>
	#else
	#include <GL/glut.h>
	#endif
		        
		/************************************
		 **
		 **    OpenGL Graphics/Interop
		 **
		 ***********************************/

		// OpenGL/CUDA interop (CUDA 2.0+)
		typedef CUresult CUDAAPI tcuGLCtxCreate( CUcontext *pCtx, unsigned int Flags, CUdevice device );
		typedef CUresult CUDAAPI tcuGraphicsGLRegisterBuffer( CUgraphicsResource *pCudaResource, GLuint buffer, unsigned int Flags );
		typedef CUresult CUDAAPI tcuGraphicsGLRegisterImage( CUgraphicsResource *pCudaResource, GLuint image, GLenum target, unsigned int Flags );

	#ifdef _WIN32
        #include <GL/wglext.h>
		// WIN32
		typedef CUresult CUDAAPI tcuWGLGetDevice( CUdevice *pDevice, HGPUNV hGpu );
	#endif

#endif // CUDA_INIT_OPENGL

#endif // __cuda_drvapi_dynlink_cuda_gl_h__

