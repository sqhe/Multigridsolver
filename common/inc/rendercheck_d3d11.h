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

 
#pragma once

#ifndef _RENDERCHECK_D3D11_H_
#define _RENDERCHECK_D3D11_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <d3d11.h>
#include <d3dx11.h>

class CheckRenderD3D11
{
public:

	CheckRenderD3D11() {}

	static HRESULT ActiveRenderTargetToPPM(ID3D11Device  *pDevice, const char *zFileName);
	static HRESULT ResourceToPPM(ID3D11Device *pDevice, ID3D11Resource *pResource, const char *zFileName);

	static bool PPMvsPPM( const char *src_file, const char *ref_file, const char *exec_path, 
                          const float epsilon, const float threshold = 0.0f );
};

#endif