// Copyright 2013 Intel Corporation
// All Rights Reserved
//
// Permission is granted to use, copy, distribute and prepare derivative works of this
// software for any purpose and without fee, provided, that the above copyright notice
// and this statement appear in all copies.  Intel makes no representations about the
// suitability of this software for any purpose.  THIS SOFTWARE IS PROVIDED "AS IS."
// INTEL SPECIFICALLY DISCLAIMS ALL WARRANTIES, EXPRESS OR IMPLIED, AND ALL LIABILITY,
// INCLUDING CONSEQUENTIAL AND OTHER INDIRECT DAMAGES, FOR THE USE OF THIS SOFTWARE,
// INCLUDING LIABILITY FOR INFRINGEMENT OF ANY PROPRIETARY RIGHTS, AND INCLUDING THE
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  Intel does not
// assume any responsibility for any errors which may appear in this software nor any
// responsibility to update it.
#include "Common.fxh"
#include "CloudsCommon.fxh"
#include "IntelExtensions.hlsl"

#ifndef CLOUD_DENSITY_TEX_DIM
#   define CLOUD_DENSITY_TEX_DIM float2(1024, 1024)
#endif

// Downscale factor for cloud color, transparency and distance buffers
#ifndef BACK_BUFFER_DOWNSCALE_FACTOR
#   define BACK_BUFFER_DOWNSCALE_FACTOR 2
#endif

#ifndef LIGHT_SPACE_PASS
#   define LIGHT_SPACE_PASS 0
#endif

#ifndef VOLUMETRIC_BLENDING
#   define VOLUMETRIC_BLENDING 1
#endif

#if !PS_ORDERING_AVAILABLE
#   undef VOLUMETRIC_BLENDING
#   define VOLUMETRIC_BLENDING 0
#endif

static const float g_fCloudExtinctionCoeff = 100;

// Minimal cloud transparancy not flushed to zero
static const float g_fTransparencyThreshold = 0.01;

// Fraction of the particle cut off distance which serves as
// a transition region from particles to flat clouds
static const float g_fParticleToFlatMorphRatio = 0.2;

static const float g_fTimeScale = 1.f;
static const float2 g_f2CloudDensitySamplingScale = float2(1.f / 200000.f, 1.f / 20000.f);

Texture2DArray<float>  g_tex2DLightSpaceDepthMap_t0 : register( t0 );
Texture2DArray<float>  g_tex2DLiSpCloudTransparency : register( t0 );
Texture2DArray<float2> g_tex2DLiSpCloudMinMaxDepth  : register( t1 );
Texture2D<float>       g_tex2DCloudDensity          : register( t1 );
Texture2D<float3>      g_tex2DWhiteNoise            : register( t3 );
Texture3D<float>       g_tex3DNoise                 : register( t4 );
Texture2D<float>       g_tex2MaxDensityMip          : register( t3 );
StructuredBuffer<uint> g_PackedCellLocations        : register( t0 );
StructuredBuffer<SCloudCellAttribs> g_CloudCells    : register( t2 );
StructuredBuffer<SParticleAttribs>  g_Particles     : register( t3 );
StructuredBuffer<uint>              g_VisibleParticleFlags   : register( t6 );
StructuredBuffer<SCloudParticleLighting> g_ParticlesLighting : register( t7 );
Texture2D<float3>       g_tex2DAmbientSkylight               : register( t7 );
Texture2DArray<float>   g_tex2DLightSpCloudTransparency      : register( t6 );
Texture2DArray<float2>  g_tex2DLightSpCloudMinMaxDepth       : register( t7 );
Texture2D<int>          g_tex2DFirstKnot                     : register( t8 );
StructuredBuffer<SParticleListKnot> g_TiledParticlesList     : register( t9 );
StructuredBuffer<uint>  g_ParticleOrder                      : register( t11 );
Buffer<uint>            g_ValidCellsCounter                  : register( t0 );
StructuredBuffer<uint>  g_ValidCellsUnorderedList            : register( t1 );
StructuredBuffer<uint>  g_ValidParticlesUnorderedList        : register( t1 );
StructuredBuffer<uint>  g_SortedParticlesOrder               : register( t0 );
Texture3D<float2>       g_tex3DParticleDensityLUT    : register( t10 );
Texture3D<float>        g_tex3DSingleScatteringInParticleLUT   : register( t11 );
Texture3D<float>        g_tex3DMultipleScatteringInParticleLUT : register( t12 );

RWStructuredBuffer<SParticleLayer> g_rwbufParticleLayers : register( u3 );
StructuredBuffer<SParticleLayer> g_bufParticleLayers : register( t0 );

StructuredBuffer<float4> g_SunLightAttenuation : register( t4 );

SamplerState samLinearWrap : register( s1 );
SamplerState samPointWrap : register( s2 );

cbuffer cbPostProcessingAttribs : register( b0 )
{
    SGlobalCloudAttribs g_GlobalCloudAttribs;
};

struct SScreenSizeQuadVSOutput
{
    float4 m_f4Pos : SV_Position;
    float2 m_f2PosPS : PosPS; // Position in projection space [-1,1]x[-1,1]
};

// Vertex shader for generating screen-size quad
SScreenSizeQuadVSOutput ScreenSizeQuadVS(in uint VertexId : SV_VertexID)
{
    float4 MinMaxUV = float4(-1, -1, 1, 1);
    
    SScreenSizeQuadVSOutput Verts[4] = 
    {
        {float4(MinMaxUV.xy, 1.0, 1.0), MinMaxUV.xy}, 
        {float4(MinMaxUV.xw, 1.0, 1.0), MinMaxUV.xw},
        {float4(MinMaxUV.zy, 1.0, 1.0), MinMaxUV.zy},
        {float4(MinMaxUV.zw, 1.0, 1.0), MinMaxUV.zw}
    };

    return Verts[VertexId];
}

// This shader computes level 0 of the maximum density mip map
float RenderMaxMipLevel0PS(SScreenSizeQuadVSOutput In) : SV_Target
{
    // Compute maximum density in the 3x3 neighborhood
    // Since dimensions of the source and destination textures are the same,
    // rasterization locations are aligned with source texels
    float2 f2UV = ProjToUV(In.m_f2PosPS.xy);
    float fMaxDensity = 0;
    [unroll]
    for(int i=-1; i <= +1; ++i)
        [unroll]
        for(int j=-1; j <= +1; ++j)
        {
            float fCurrDensity = g_tex2DCloudDensity.SampleLevel(samPointWrap, f2UV, 0, int2(i,j));
            fMaxDensity = max(fMaxDensity, fCurrDensity);
        }
    return fMaxDensity;
}

// This shader renders next coarse level of the maximum density mip map
float RenderCoarseMaxMipLevelPS(SScreenSizeQuadVSOutput In) : SV_Target
{
    // Compute maximum density of 4 finer-level texels
    float2 f2MipSize; 
    float fMipLevels;
    float fSrcMip = g_GlobalCloudAttribs.f4Parameter.x-1;
    g_tex2MaxDensityMip.GetDimensions(fSrcMip, f2MipSize.x, f2MipSize.y, fMipLevels);
    // Note that since dst level is 2 times smaller than the src level, we have to
    // align texture coordinates
    //   _____ _____
    //  |     |     |  1
    //  |_____._____|
    //  |  .  |     |  0
    //  |_____|_____|
    //     0      1
    float2 f2UV = ProjToUV(In.m_f2PosPS.xy) - 0.5f / f2MipSize;
    float fMaxDensity = 0;
    [unroll]
    for(int i=0; i <= +1; ++i)
        [unroll]
        for(int j=0; j <= +1; ++j)
        {
            float fCurrDensity = g_tex2MaxDensityMip.SampleLevel(samPointWrap, f2UV, fSrcMip, int2(i,j));
            fMaxDensity = max(fMaxDensity, fCurrDensity);
        }
    return fMaxDensity;
}

float2 ComputeDensityTexLODsFromUV(in float4 fDeltaUV01)
{
    fDeltaUV01 *= CLOUD_DENSITY_TEX_DIM.xyxy;
    float2 f2UVDeltaLen = float2( length(fDeltaUV01.xy), length(fDeltaUV01.zw) );
    f2UVDeltaLen = max(f2UVDeltaLen,1 );
    return log2( f2UVDeltaLen );
}

float2 ComputeDensityTexLODsFromStep(in float fSamplingStep)
{
    float2 f2dU = fSamplingStep * g_f2CloudDensitySamplingScale * CLOUD_DENSITY_TEX_DIM.xx;
    float2 f2LODs = log2(max(f2dU, 1));
    return f2LODs;
}

float4 GetCloudDensityUV(in float3 CloudPosition, in float fTime)
{
    const float4 f2Offset01 = float4( 0.1*float2(-0.04, +0.01) * fTime, 0.2*float2( 0.01,  0.04) * fTime );
    float4 f2UV01 = CloudPosition.xzxz * g_f2CloudDensitySamplingScale.xxyy + f2Offset01; 
    return f2UV01;
}

float GetCloudDensity(in float4 f4UV01, in float2 f2LODs = float2(0,0))
{
    float fDensity = 
        g_tex2DCloudDensity.SampleLevel(samLinearWrap, f4UV01.xy, f2LODs.x) * 
        g_tex2DCloudDensity.SampleLevel(samLinearWrap, f4UV01.zw, f2LODs.y);

    fDensity = saturate((fDensity-g_GlobalCloudAttribs.fCloudDensityThreshold)/(1-g_GlobalCloudAttribs.fCloudDensityThreshold));

    return fDensity;
}

float GetCloudDensityAutoLOD(in float4 f4UV01)
{
    float fDensity = 
        g_tex2DCloudDensity.Sample(samLinearWrap, f4UV01.xy) * 
        g_tex2DCloudDensity.Sample(samLinearWrap, f4UV01.zw);

    fDensity = saturate((fDensity-g_GlobalCloudAttribs.fCloudDensityThreshold)/(1-g_GlobalCloudAttribs.fCloudDensityThreshold));

    return fDensity;
}

float GetCloudDensity(in float3 CloudPosition, in const float fTime, in float2 f2LODs = float2(0,0))
{
    float4 f4UV01 = GetCloudDensityUV(CloudPosition, fTime);
    return GetCloudDensity(f4UV01, f2LODs);
}

float GetMaxDensity(in float4 f4UV01, in float2 f2LODs = float2(0,0))
{
    float fDensity = 
        g_tex2MaxDensityMip.SampleLevel(samPointWrap, f4UV01.xy, f2LODs.x) * 
        g_tex2MaxDensityMip.SampleLevel(samPointWrap, f4UV01.zw, f2LODs.y);

    fDensity = saturate((fDensity-g_GlobalCloudAttribs.fCloudDensityThreshold)/(1-g_GlobalCloudAttribs.fCloudDensityThreshold));

    return fDensity;
}

float GetMaxDensity(in float3 CloudPosition, in const float fTime, in float2 f2LODs = float2(0,0))
{
    float4 f4UV01 = GetCloudDensityUV(CloudPosition, fTime);
    return GetMaxDensity(f4UV01, f2LODs);
}


// This function performs bilateral upscaling of the cloud color, transparency and distance buffers
void FilterDownscaledCloudBuffers(in float2 f2UV,
                                  in float fDistToCamera,
                                  out float fCloudTransparency,
                                  out float3 f3CloudsColor,
                                  out float fDistToCloud)
{
    fCloudTransparency = 1;
    f3CloudsColor = 0;
    fDistToCloud = +FLT_MAX;

    float2 f2CloudBufferSize = float2(g_GlobalCloudAttribs.fDownscaledBackBufferWidth, g_GlobalCloudAttribs.fDownscaledBackBufferHeight);
    // Get location of the left bottom source texel
    float2 f2SrcIJ = f2UV * f2CloudBufferSize - 0.5;
    float2 f2SrcIJ0 = floor(f2SrcIJ);
    float2 f2UVWeights = f2SrcIJ - f2SrcIJ0;
    // Compute UV coordinates of the gather location
    float2 f2GatherUV = (f2SrcIJ0+1) / f2CloudBufferSize;

    // The values in float4, which Gather() returns are arranged as follows:
    //   _______ _______
    //  |       |       |
    //  |   x   |   y   |
    //  |_______o_______|  o gather location
    //  |       |       |
    //  |   w   |   z   |  
    //  |_______|_______|
    
    // Read data from the source buffers
    float4 f4SrcDistToCloud = g_tex2DScrSpaceCloudMinMaxDist.GatherRed  ( samLinearClamp, f2GatherUV );
    float4 f4SrcCloudTransp = g_tex2DScrSpaceCloudTransparency.GatherRed( samLinearClamp, f2GatherUV );
    float4 f4CloudColorR = g_tex2DScrSpaceCloudColor.GatherRed  ( samLinearClamp, f2GatherUV );
    float4 f4CloudColorG = g_tex2DScrSpaceCloudColor.GatherGreen( samLinearClamp, f2GatherUV );
    float4 f4CloudColorB = g_tex2DScrSpaceCloudColor.GatherBlue ( samLinearClamp, f2GatherUV );
    
    // Compute bilateral weights, start by bilinear:
    float4 f4BilateralWeights = float4(1 - f2UVWeights.x, f2UVWeights.x,   f2UVWeights.x, 1-f2UVWeights.x) * 
                                float4(    f2UVWeights.y, f2UVWeights.y, 1-f2UVWeights.y, 1-f2UVWeights.y);
    
    // Take into account only these source texels, which are closer to the camera than the opaque surface
    //
    //              . .
    //   ---------->.'.      /\
    //   ----->/\   .'.     /  \
    //   -----/--\->.'.    /    \
    //       /    \       /      \
    //
    // To assure smooth filtering at cloud boundaries especially when cloud is in front of opaque surfaces,
    // we also need to account for transparent pixels (identified by dist to clouds == +FLT_MAX)
    //       
    //       Finite +FLT_MAX   
    //        dist 
    //         |     |
    //         |     |
    //      ...V...  |
    //     ''''''''' |
    //               |
    //   ____________|___
    //               |
    //               V 
    //
    f4BilateralWeights.xyzw *= float4( f4SrcDistToCloud.xyzw < fDistToCamera || f4SrcDistToCloud.xyzw ==+FLT_MAX );
    
    float fSumWeight = dot(f4BilateralWeights.xyzw, 1);
    
    if( fSumWeight > 1e-2 )
    {
        f4BilateralWeights /= fSumWeight;
    
        fCloudTransparency = dot(f4SrcCloudTransp, f4BilateralWeights);
    
        f3CloudsColor.r = dot(f4CloudColorR, f4BilateralWeights);
        f3CloudsColor.g = dot(f4CloudColorG, f4BilateralWeights);
        f3CloudsColor.b = dot(f4CloudColorB, f4BilateralWeights);
        
        // Use minimum distance to avoid filtering FLT_MAX
        fDistToCloud = min(f4SrcDistToCloud.x, f4SrcDistToCloud.y);
        fDistToCloud = min(fDistToCloud, f4SrcDistToCloud.z);
        fDistToCloud = min(fDistToCloud, f4SrcDistToCloud.w);
    }
}


// This shader renders flat clouds by sampling cloud density at intersection
// of the view ray with the cloud layer
// If particles are rendred in lower resolution, it also upscales the
// downscaled buffers and combines with the result
void RenderFlatCloudsPS(SScreenSizeQuadVSOutput In,
                          out float fTransparency   : SV_Target0,
                          out float2 f2MinMaxZRange : SV_Target1
                          #if !LIGHT_SPACE_PASS
                            , out float4 f4Color : SV_Target2
                          #endif
                          )
{
    // Load depth from the depth buffer
    float fDepth;
    float3 f3RayStart;
    float2 f2UV = ProjToUV(In.m_f2PosPS.xy);
#if LIGHT_SPACE_PASS
    fDepth = g_tex2DLightSpaceDepthMap_t0.SampleLevel(samLinearClamp, float3(f2UV,g_GlobalCloudAttribs.f4Parameter.x), 0 );
    // For directional light source, we should use position on the near clip plane instead of
    // camera location as a ray start point (use 1.01 to avoid issues when depth == 1)
    float4 f4PosOnNearClipPlaneWS = mul( float4(In.m_f2PosPS.xy,1.01,1), g_CameraAttribs.mViewProjInv );
    f3RayStart = f4PosOnNearClipPlaneWS.xyz/f4PosOnNearClipPlaneWS.w;
#else
    f3RayStart = g_CameraAttribs.f4CameraPos.xyz;
    fDepth = g_tex2DDepthBuffer.Load(int3(In.m_f4Pos.xy,0));
#endif

    // Reconstruct world space position
    float4 f4ReconstructedPosWS = mul( float4(In.m_f2PosPS.xy,fDepth,1), g_CameraAttribs.mViewProjInv );
    float3 f3WorldPos = f4ReconstructedPosWS.xyz / f4ReconstructedPosWS.w;

    // Compute view ray
    float3 f3ViewDir = f3WorldPos - f3RayStart;
    float fDistToCamera = length(f3ViewDir);
    f3ViewDir /= fDistToCamera;
    float fRayLength = fDistToCamera;

    float fTime = g_fTimeScale*g_GlobalCloudAttribs.fTime;

    // Handle the case when the ray does not hit any surface
    // When rendering from light, we do not need to trace the ray
    // further than the far clipping plane because there is nothing 
    // visible there
#if !LIGHT_SPACE_PASS
    if( fDepth < 1e-10 )
        fRayLength = + FLT_MAX;
#endif

    // Compute intersection of the view ray with the Earth and the spherical cloud layer
    float3 f3EarthCentre = float3(0, -g_MediaParams.fEarthRadius, 0);
    float4 f4CloudLayerAndEarthIsecs;
    GetRaySphereIntersection2(f3RayStart.xyz,
                              f3ViewDir,
                              f3EarthCentre,
                              float2(g_MediaParams.fEarthRadius, g_MediaParams.fEarthRadius + g_GlobalCloudAttribs.fCloudAltitude),
                              f4CloudLayerAndEarthIsecs);
    float2 f2EarthIsecs = f4CloudLayerAndEarthIsecs.xy;
    float2 f2CloudLayerIsecs = f4CloudLayerAndEarthIsecs.zw;
    
    bool bIsValid = true;
    // Check if the view ray does not hit the cloud layer
    if( f2CloudLayerIsecs.y < f2CloudLayerIsecs.x )
        bIsValid = false;

    float fFadeOutFactor = 1;
#if LIGHT_SPACE_PASS
    // For light space pass, always use the first intersection with the cloud layer
    float fDistToCloudLayer = f2CloudLayerIsecs.x;
    if(f2EarthIsecs.x < f2EarthIsecs.y)
        fRayLength = min(fRayLength, f2EarthIsecs.x);
#else
    // For the camera space pass, select either first or second intersection
    // If the camera is under the cloud layer, use second intersection
    // If the camera is above the cloud layer, use first intersection
    float fDistToCloudLayer = f2CloudLayerIsecs.x > 0 ? f2CloudLayerIsecs.x : f2CloudLayerIsecs.y;

    // Region [fParticleCutOffDist - fFadeOutDistance, fParticleCutOffDist] servers as transition
    // from particles to flat clouds
    float fFadeOutDistance = g_GlobalCloudAttribs.fParticleCutOffDist * g_fParticleToFlatMorphRatio;
    float fFadeOutStartDistance = g_GlobalCloudAttribs.fParticleCutOffDist - fFadeOutDistance;
    
    if( fDistToCloudLayer < fFadeOutStartDistance )
        bIsValid = false;

    // Compute fade out factor
    fFadeOutFactor = saturate( (fDistToCloudLayer - fFadeOutStartDistance) /  max(fFadeOutDistance,1) );
#endif

    if( fDistToCloudLayer > fRayLength  )
        bIsValid = false;

#if LIGHT_SPACE_PASS || BACK_BUFFER_DOWNSCALE_FACTOR == 1
    if( !bIsValid) discard;
#endif

    float fTotalMass = 0;
    float3 f3CloudLayerIsecPos = 0;
    float fDensity = 0;
    float fCloudPathLen = 0;
    if( bIsValid)
    {
        // Get intersection point
        f3CloudLayerIsecPos = f3RayStart.xyz + f3ViewDir * fDistToCloudLayer;

        // Get cloud density at intersection point
        float4 f4UV01 = GetCloudDensityUV(f3CloudLayerIsecPos, fTime);
        fDensity = GetCloudDensityAutoLOD(f4UV01)* fFadeOutFactor;
        // Fade out clouds when view angle is orthogonal to zenith
        float3 f3ZenithDir = normalize(f3CloudLayerIsecPos - f3EarthCentre);
        float fCosZenithAngle = dot(f3ViewDir, f3ZenithDir);
        fDensity *= abs(fCosZenithAngle)*2;

        fCloudPathLen = g_GlobalCloudAttribs.fCloudThickness;
        fTotalMass = fCloudPathLen * g_GlobalCloudAttribs.fCloudVolumeDensity * fDensity;

        // This helps improve perofrmance by reducing memory bandwidth
        if(fTotalMass < 1e-5)
            bIsValid = false;    
    }

#if LIGHT_SPACE_PASS || BACK_BUFFER_DOWNSCALE_FACTOR == 1
    if(!bIsValid) discard;
#endif

    fTransparency = exp(-g_fCloudExtinctionCoeff*fTotalMass);

#if LIGHT_SPACE_PASS
    // Transform intersection point into light view space
    float4 f4LightSpacePosPS = mul( float4(f3CloudLayerIsecPos,1), g_CameraAttribs.WorldViewProj );
    f2MinMaxZRange.xy = f4LightSpacePosPS.z / f4LightSpacePosPS.w;
#else
    f4Color = 0;
    if( bIsValid )
    {
        float3 f3SunLightExtinction, f3AmbientLight;
        GetSunLightExtinctionAndSkyLight(f3CloudLayerIsecPos, f3SunLightExtinction, f3AmbientLight, g_tex2DOccludedNetDensityToAtmTop, g_tex2DAmbientSkylight);
        f4Color.rgb = fDensity * fCloudPathLen * g_GlobalCloudAttribs.fCloudVolumeDensity * f3SunLightExtinction;
    }
    else
        fDistToCloudLayer = +FLT_MAX;

#   if BACK_BUFFER_DOWNSCALE_FACTOR > 1
        // Upscale buffers
        float fParticleDistToCloud = 0;
        float3 f3ParticleColor = 0;
        float fParticleTransparency = 0;
        FilterDownscaledCloudBuffers( f2UV,
                                      fDistToCamera,
                                      fParticleTransparency,
                                      f3ParticleColor,
                                      fParticleDistToCloud);
        // Combine with flat clouds
        f4Color.rgb = f4Color.rgb * fParticleTransparency + f3ParticleColor;
        fDistToCloudLayer = min(fDistToCloudLayer, fParticleDistToCloud);
        fTransparency *= fParticleTransparency;
        // Save bandwidth by not outputting fully transparent clouds
        if( fTransparency > 1 - 1e-5 )
            discard;
#   endif

    f4Color.a = fTransparency;
    f2MinMaxZRange.xy = fDistToCloudLayer;
#endif
}


technique11 RenderFlatCloudsTech
{
    pass
    {
        SetVertexShader( CompileShader( vs_5_0, ScreenSizeQuadVS() ) );
        SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_5_0, RenderFlatCloudsPS() ) );
    }
}


struct PS_Input
{
    float4 f4Pos : SV_Position;
    nointerpolation uint uiParticleID : PARTICLE_ID;
};

float3 GetParticleScales(in float fSize, in float fNumActiveLayers)
{
    float3 f3Scales = fSize;
    //if( fNumActiveLayers > 1 )
    //    f3Scales.y = max(f3Scales.y, g_GlobalCloudAttribs.fCloudThickness/fNumActiveLayers);
    f3Scales.y = min(f3Scales.y, g_GlobalCloudAttribs.fCloudThickness/2.f);
    return f3Scales;
}

RWStructuredBuffer<SCloudCellAttribs> g_CloudCellsRW        : register( u0 );
AppendStructuredBuffer<uint>          g_ValidCellsAppendBuf : register( u1 );

#ifndef THREAD_GROUP_SIZE
#   define THREAD_GROUP_SIZE 64
#endif

// This shader processes each cell of the cloud grid and outputs all valid cells to the
// append buffer
[numthreads(THREAD_GROUP_SIZE, 1, 1)]
void ProcessCloudGridCS( uint3 Gid  : SV_GroupID, 
                         uint3 GTid : SV_GroupThreadID )
{
    uint uiCellId = Gid.x * THREAD_GROUP_SIZE + GTid.x;
    if( uiCellId >= g_GlobalCloudAttribs.uiNumCells )
        return;

    // Get cell location in the grid
    uint uiPackedLocation = g_PackedCellLocations[uiCellId];
    uint uiCellI, uiCellJ, uiRing, uiLayer;
    UnPackParticleIJRing(uiPackedLocation, uiCellI, uiCellJ, uiRing, uiLayer);

    // Compute cell center world space coordinates
    uint uiRingDimension = g_GlobalCloudAttribs.uiRingDimension;
    const float fRingWorldStep = GetCloudRingWorldStep(uiRing, g_GlobalCloudAttribs);

    //
    // 
    //                                 Camera        
    //                               |<----->|
    //   |   X   |   X   |   X   |   X   |   X   |   X   |   X   |   X   |       CameraI == 4
    //   0  0.5     1.5     2.5     3.5  4  4.5     5.5     6.5     7.5  8       uiRingDimension == 8
    //                                   |
    //                                CameraI
    float fCameraI = floor(g_CameraAttribs.f4CameraPos.x/fRingWorldStep + 0.5);
    float fCameraJ = floor(g_CameraAttribs.f4CameraPos.z/fRingWorldStep + 0.5);

    float3 f3CellCenter;
    f3CellCenter.x = (fCameraI + (float)uiCellI - (uiRingDimension/2) + 0.5) * fRingWorldStep;
    f3CellCenter.z = (fCameraJ + (float)uiCellJ - (uiRingDimension/2) + 0.5) * fRingWorldStep;
    f3CellCenter.y = 0;

    float3 f3EarthCentre = float3(0, -g_MediaParams.fEarthRadius, 0);
    float3 f3DirFromEarthCenter = f3CellCenter - f3EarthCentre;
    float fDistFromCenter = length(f3DirFromEarthCenter);
    f3CellCenter = f3EarthCentre + f3DirFromEarthCenter * ((g_MediaParams.fEarthRadius + g_GlobalCloudAttribs.fCloudAltitude)/fDistFromCenter);

    uint uiNumActiveLayers = GetNumActiveLayers(g_GlobalCloudAttribs.uiMaxLayers, uiRing);

    float fParticleSize = GetParticleSize(fRingWorldStep);
    float3 f3Size = GetParticleScales(fParticleSize, uiNumActiveLayers);
   
    // Construct local frame
    float3 f3Normal = normalize(f3CellCenter.xyz - f3EarthCentre);
    float3 f3Tangent = normalize( cross(f3Normal, float3(0,0,1)) );
    float3 f3Bitangent = normalize( cross(f3Tangent, f3Normal) );

    float fTime = g_fTimeScale*g_GlobalCloudAttribs.fTime;

    // Get cloud density in the cell
    float4 f4UV01 = GetCloudDensityUV(f3CellCenter, fTime);
    float2 f2LODs = ComputeDensityTexLODsFromStep(f3Size.x*2);
    float fMaxDensity = GetMaxDensity( f4UV01, f2LODs );
    
    bool bIsValid = true;
    if( fMaxDensity < 1e-5 )
        bIsValid = false;

    float fDensity = 0;
    float fMorphFadeout = 1;
    if( bIsValid )
    {
        fDensity = saturate(GetCloudDensity( f4UV01, f2LODs ));

        // Compute morph weights for outer and inner boundaries
        {
            float4 f4OuterBoundary = g_CameraAttribs.f4CameraPos.xzxz + float4(-1,-1,+1,+1) * (float)(uiRingDimension/2) * fRingWorldStep;

            //f4OuterBoundary.x                                                  f4OuterBoundary.z
            //      |                                                               |
            //      |       uiRingDimension/2              uiRingDimension/2        |
            //      |<----------------------------->C<----------------------------->|
            //                               |       |
            //   |   X   |   X   |   X   |   X   |   X   |   X   |   X   |   X   |

            float4 f4DistToOuterBnd = float4(1,1,-1,-1)*(f3CellCenter.xzxz  - f4OuterBoundary.xyzw);
            float fMinDist = min(f4DistToOuterBnd.x, f4DistToOuterBnd.y);
            fMinDist = min(fMinDist, f4DistToOuterBnd.z);
            fMinDist = min(fMinDist, f4DistToOuterBnd.w);
            float fOuterMorphRange = g_GlobalCloudAttribs.uiRingExtension * fRingWorldStep;
            float fOuterMorphWeight = saturate( fMinDist / fOuterMorphRange);
            fMorphFadeout *= fOuterMorphWeight;
        }

        if(uiRing > 0)
        {
            float4 f4InnerBoundary = g_CameraAttribs.f4CameraPos.xzxz + float4(-1,-1,+1,+1) * (float)(g_GlobalCloudAttribs.uiInnerRingDim/4 + g_GlobalCloudAttribs.uiRingExtension/2) * fRingWorldStep;

            //               f4InnerBoundary.x                f4InnerBoundary.z
            //                        |                             |
            //                        |                             |
            //                        |<----------->C<------------->|               
            //                               |       |
            //   |   X   |   X   |   X   |   X   |   X   |   X   |   X   |   X   |

            float4 f4DistToInnerBnd = float4(1,1,-1,-1)*(f3CellCenter.xzxz - f4InnerBoundary.xyzw);
            float fMinDist = min(f4DistToInnerBnd.x, f4DistToInnerBnd.y);
            fMinDist = min(fMinDist, f4DistToInnerBnd.z);
            fMinDist = min(fMinDist, f4DistToInnerBnd.w);
            float fInnerMorphRange = g_GlobalCloudAttribs.uiRingExtension/2 * fRingWorldStep;
            float fInnerMorphWeight = 1-saturate( fMinDist / fInnerMorphRange);
            fMorphFadeout *= fInnerMorphWeight;
        }
        
        if( fDensity < 1e-5 )
            bIsValid = false;

        // TODO: perform this check for each particle, not cell:
        float fParticleBoundSphereRadius = length(f3Size);
        if( length(f3CellCenter - g_CameraAttribs.f4CameraPos.xyz) > g_GlobalCloudAttribs.fParticleCutOffDist + fParticleBoundSphereRadius )
            bIsValid = false;
    }
    
    if( bIsValid )
    {
        // If the cell is valid, store the data in the buffer
        g_CloudCellsRW[uiCellId].f3Center = f3CellCenter;
        g_CloudCellsRW[uiCellId].fSize = f3Size.x;
                   
        g_CloudCellsRW[uiCellId].f3Normal.xyz    = f3Normal;
        g_CloudCellsRW[uiCellId].f3Tangent.xyz   = f3Tangent;
        g_CloudCellsRW[uiCellId].f3Bitangent.xyz = f3Bitangent;
        
        g_CloudCellsRW[uiCellId].uiNumActiveLayers = uiNumActiveLayers;

        g_CloudCellsRW[uiCellId].fDensity = fDensity;
        g_CloudCellsRW[uiCellId].fMorphFadeout = fMorphFadeout;
        
        g_CloudCellsRW[uiCellId].uiPackedLocation = uiPackedLocation;

        // Append the cell ID to the list
        g_ValidCellsAppendBuf.Append(uiCellId);
    }
}

RWBuffer<uint> g_DispatchArgsRW : register( u0 );

// This shader compute dispatch arguments for the DispatchIndirect() method
[numthreads(1, 1, 1)]
void ComputeDispatchArgsCS( uint3 Gid  : SV_GroupID, 
                            uint3 GTid : SV_GroupThreadID )
{
    g_DispatchArgsRW[0] = (g_ValidCellsCounter.Load(0) + THREAD_GROUP_SIZE-1) / THREAD_GROUP_SIZE;
}

RWStructuredBuffer<SParticleAttribs> g_CloudParticlesRW        : register( u0 );
AppendStructuredBuffer<uint>         g_ValidParticlesAppendBuf : register( u1 );

// This shader processes all valid cells and for each cell outputs
// appropriate number of particles in this cell
[numthreads(THREAD_GROUP_SIZE, 1, 1)]
void ProcessValidParticlesCS( uint3 Gid  : SV_GroupID, 
                              uint3 GTid : SV_GroupThreadID )
{
    uint uiValidCellNum = Gid.x * THREAD_GROUP_SIZE + GTid.x;
    uint uiNumValidCells = g_ValidCellsCounter.Load(0);
    if( uiValidCellNum >= uiNumValidCells )
        return;

    float fTime = g_fTimeScale*g_GlobalCloudAttribs.fTime;

    uint uiMaxLayers = g_GlobalCloudAttribs.uiMaxLayers;
    // Load valid cell id from the list
    uint uiCellId = g_ValidCellsUnorderedList[uiValidCellNum];
    // Get the cell attributes
    SCloudCellAttribs CellAttrs = g_CloudCells[uiCellId];
    uint uiDstParticleId = uiCellId * uiMaxLayers;
    float3 f3Size = GetParticleScales(CellAttrs.fSize, CellAttrs.uiNumActiveLayers);
    uint uiNumActiveLayers = CellAttrs.uiNumActiveLayers;
    for(uint uiLayer = 0; uiLayer < uiNumActiveLayers; ++uiLayer)
    {
        // Process each layer in the cell
        float3 f3CloudPos = CellAttrs.f3Center;
        //float fLayerAltitude = ((float)uiLayer + 0.5) / (float)uiNumActiveLayers - 0.5;
        //f3CloudPos += CellAttrs.f3Normal.xyz * fLayerAltitude * g_GlobalCloudAttribs.fCloudThickness;
        f3CloudPos += CellAttrs.f3Normal.xyz * (float)uiLayer * f3Size.y * 0.9;
        
        // Apply random displacement to particle pos
        float3 f3Noise = g_tex2DWhiteNoise.SampleLevel(samLinearWrap, (f3CloudPos.yx + f3CloudPos.xz)/1500+0*fTime*5e-4, 0).xyz;
        float3 f3RandomDisplacement = (f3Noise*2-1) * f3Size / 2.f * float3(0.8,0.3,0.8);
        f3CloudPos += f3RandomDisplacement.x * CellAttrs.f3Tangent.xyz + 
                      f3RandomDisplacement.y * CellAttrs.f3Normal.xyz + 
                      f3RandomDisplacement.z * CellAttrs.f3Bitangent.xyz;
        
        float fDensityStep = min(0.05, 1.f/(float)uiNumActiveLayers);
        float fLayerMinDensity = fDensityStep * (float)uiLayer;
        float fLayerMaxDensity = fDensityStep * (float)(uiLayer+1);

        float RelativeDensity = saturate(( CellAttrs.fDensity - fLayerMinDensity ) / ( fLayerMaxDensity - fLayerMinDensity ) );
        float fNoise = dot(f3Noise, 1.f/3.f);
        RelativeDensity *= lerp(1, f3Noise.z, 0.1);
        float fRandomness = lerp(0.5, 0.9, (float)uiLayer / (float)uiNumActiveLayers );
        const float fSizeScale = lerp(1.5, 0.8, (float)uiLayer / (float)uiNumActiveLayers );
        float fParticleSize = CellAttrs.fSize * RelativeDensity * fSizeScale * (1 + (f3Noise.x-0.5) * fRandomness);
        if( fParticleSize > 1e-5 )
        {   
            g_CloudParticlesRW[uiDstParticleId].f3Pos = f3CloudPos;
            g_CloudParticlesRW[uiDstParticleId].fSize = fParticleSize;
            g_CloudParticlesRW[uiDstParticleId].fRndAzimuthBias = f3Noise.y+(f3Noise.x-0.5)*fTime*5e-2;
            g_CloudParticlesRW[uiDstParticleId].fDensity = pow( saturate((RelativeDensity-0.5)*2),5);

            // Append valid particle to the list
            g_ValidParticlesAppendBuf.Append(uiDstParticleId);
        }

        ++uiDstParticleId;
    }
}

RWStructuredBuffer<SCloudParticleLighting> g_ParticlesLightingRW : register(u0);
RWStructuredBuffer<float4>                 g_SunLightAttenuationRW : register(u1);

// This function computes visibility for the particle
bool IsParticleVisibile(in float3 f3Center, in float3 f3Scales, float4 f4ViewFrustumPlanes[6])
{
    float fParticleBoundSphereRadius = length(f3Scales);
    bool bIsVisible = true;
    for(int iPlane = 0; iPlane < 6; ++iPlane)
    {
#if LIGHT_SPACE_PASS
        // Do not clip against far clipping plane for light pass
        if( iPlane == 5 )
            continue;
#endif
        float4 f4CurrPlane = f4ViewFrustumPlanes[iPlane];
#if 1
        // Note that the plane normal is not normalized to 1
        float DMax = dot(f3Center.xyz, f4CurrPlane.xyz) + f4CurrPlane.w + fParticleBoundSphereRadius*length(f4CurrPlane.xyz);
#else
        // This is a bit more accurate but significantly more computationally expensive test
        float DMax = -FLT_MAX;
        for(uint uiCorner=0; uiCorner < 8; ++uiCorner)
        {
            float4 f4CurrCornerWS = ParticleAttrs.f4BoundBoxCornersWS[uiCorner];
            float D = dot( f4CurrCornerWS.xyz, f4CurrPlane.xyz) + f4CurrPlane.w;
            DMax = max(DMax, D);
        }
#endif
        if( DMax < 0 )
        {
            bIsVisible = false;
        }
    }
    return bIsVisible;
}

RWStructuredBuffer<uint> g_ParticlesVisibilityFlagsRW : register(u0);

// This shader computes visibility for each valid particle and packs them in
// as 1-bit flags
[numthreads(THREAD_GROUP_SIZE, 1, 1)]
void ComputeParticlesVisibilityCS(uint3 Gid  : SV_GroupID, 
                                  uint3 GTid : SV_GroupThreadID)
{
    uint uiValidParticleSerialNum = Gid.x * THREAD_GROUP_SIZE + GTid.x;
    uint uiNumValidParticles = g_ValidCellsCounter.Load(0);
    if( uiValidParticleSerialNum >= uiNumValidParticles )
        return;

    uint uiParticleId = g_ValidParticlesUnorderedList[uiValidParticleSerialNum];
    uint uiMaxLayers = g_GlobalCloudAttribs.uiMaxLayers;
    uint uiCellId = uiParticleId/uiMaxLayers;
    SCloudCellAttribs CellAttrs = g_CloudCells[uiCellId];

    float3 f3Size = GetParticleScales(CellAttrs.fSize, CellAttrs.uiNumActiveLayers);
    float3 f3Pos = g_Particles[uiParticleId].f3Pos;


    bool bIsVisible = IsParticleVisibile(f3Pos, f3Size, 
#if LIGHT_SPACE_PASS
                                         g_GlobalCloudAttribs.f4TilingFrustumPlanes
#else
                                         g_CameraAttribs.f4ViewFrustumPlanes
#endif
    );

    // Set appropriate 1-bit flag if the particle is visible
	if( bIsVisible )
		InterlockedOr( g_ParticlesVisibilityFlagsRW[uiParticleId/32], 1 << (uiParticleId&31) );

    ++uiParticleId;
}

void IntersectRayWithParticle(const in SParticleAttribs ParticleAttrs,
                              const in SCloudCellAttribs CellAttrs,
                              const in float3 f3CameraPos, 
                              const in float3 f3ViewRay,
                              out float2 f2RayIsecs,
                              out float3 f3EntryPointUSSpace,
                              out float3 f3ViewRayUSSpace,
                              out float3 f3LightDirUSSpace,
                              out float fDistanceToEntryPoint,
                              out float fDistanceToExitPoint);

float2 ComputeLiSpTransparency(const in SParticleAttribs ParticleAttrs)
{
    float4 f4ParticleCenterPosPS = mul( float4(ParticleAttrs.f3Pos,1), g_GlobalCloudAttribs.mParticleTiling );
    float2 f2ParticleProjSpaceXY = f4ParticleCenterPosPS.xy/f4ParticleCenterPosPS.w;
    if( any( abs(f2ParticleProjSpaceXY) > 1 ) )
        return 1;

    float2 f2ParticleUV = saturate( ProjToUV(f2ParticleProjSpaceXY) );

    float3 f3RayStart = ParticleAttrs.f3Pos.xyz;
    float3 f3RayDir = g_LightAttribs.f4DirOnLight.xyz;

    uint2 uiTileXY = floor(f2ParticleUV * g_GlobalCloudAttribs.uiLiSpFirstListIndTexDim);
    int iCurrKnotInd = g_tex2DFirstKnot.Load( uint3(uiTileXY, 0) );
    float2 f2Transparency = 1;
    
    [loop]
    while( iCurrKnotInd >= 0 )
    {
        SParticleListKnot CurrKnot = g_TiledParticlesList[iCurrKnotInd];
        SParticleAttribs CurrParticleAttrs = g_Particles[CurrKnot.uiParticleID];
        SCloudCellAttribs CurrCellAttrs = g_CloudCells[CurrKnot.uiParticleID / g_GlobalCloudAttribs.uiMaxLayers];

        float2 f2RayIsecs;
        float fDistanceToEntryPoint, fDistanceToExitPoint;
        float3 f3EntryPointUSSpace, f3ViewRayUSSpace, f3LightDirUSSpace;
        IntersectRayWithParticle(CurrParticleAttrs, CurrCellAttrs, f3RayStart,  f3RayDir,
                                 f2RayIsecs, f3EntryPointUSSpace, f3ViewRayUSSpace,
                                 f3LightDirUSSpace,
                                 fDistanceToEntryPoint, fDistanceToExitPoint);
#if PS_ORDERING_AVAILABLE
        if( f2RayIsecs.y > f2RayIsecs.x && fDistanceToExitPoint < 0 )
            break;
#endif
        
        float2 f2CurrTransparency = 1.f;
        if( //CurrKnot.uiParticleID != uiParticleOrder && 
            fDistanceToExitPoint > 0 
            /*f2RayIsecs.y > f2RayIsecs.x && */)
        {
            //float fDensityScale = 1;
            //float fDensity = ParticleAttrs.fDensity * fDensityScale;
            fDistanceToEntryPoint = max(fDistanceToEntryPoint, CurrParticleAttrs.fSize);
            float fCloudMass = max(fDistanceToExitPoint-fDistanceToEntryPoint, 0);//g_GlobalCloudAttribs.fCloudVolumeDensity * fDensity;
            fCloudMass *= CurrCellAttrs.fMorphFadeout * ParticleAttrs.fDensity;
            f2CurrTransparency = exp( -fCloudMass * g_GlobalCloudAttribs.fAttenuationCoeff * float2(0.05,0.025) );
        }

        f2Transparency *= f2CurrTransparency;
        if( all(f2Transparency < g_fTransparencyThreshold) )
            break;

        iCurrKnotInd = CurrKnot.iNextKnotInd;
    }

    //fTransparency = saturate((fTransparency - g_fTransparencyThreshold) / (1-g_fTransparencyThreshold));

    return f2Transparency;
}

// This shader compute lighting for all valid particles
[numthreads(THREAD_GROUP_SIZE, 1, 1)]
void ComputeParticlesLightingCS(uint3 Gid  : SV_GroupID, 
                                uint3 GTid : SV_GroupThreadID)
{
    uint uiParticleSerialNum = Gid.x * THREAD_GROUP_SIZE + GTid.x;
    uint uiNumValidParticles = g_ValidCellsCounter.Load(0);
    if( uiParticleSerialNum >= uiNumValidParticles )
        return;

    uint uiParticleId = g_ValidParticlesUnorderedList[uiParticleSerialNum];

	bool bIsVisible = g_VisibleParticleFlags[uiParticleId/32] & (1 << (uiParticleId&31));
    if( !bIsVisible )
        return;

    uint uiCellId = uiParticleId / g_GlobalCloudAttribs.uiMaxLayers;
    SParticleAttribs ParticleAttrs = g_Particles[uiParticleId];
    
    float3 f3AtmosphereExtinction, f3AmbientSkyLight;
    GetSunLightExtinctionAndSkyLight(ParticleAttrs.f3Pos, f3AtmosphereExtinction, f3AmbientSkyLight, g_tex2DOccludedNetDensityToAtmTop, g_tex2DAmbientSkylight);
    g_ParticlesLightingRW[uiParticleId].f4SunLight = float4(f3AtmosphereExtinction * g_LightAttribs.f4ExtraterrestrialSunColor.rgb, 0);
    g_ParticlesLightingRW[uiParticleId].f4AmbientLight = float4(f3AmbientSkyLight, 1);

    float2 f2CloudLiSpTransparency = ComputeLiSpTransparency(ParticleAttrs);

    SCloudCellAttribs CellAttrs = g_CloudCells[uiCellId];
    uint uiCellX, uiCellY, uiRing, uiLayerUnused;
    UnPackParticleIJRing(CellAttrs.uiPackedLocation, uiCellX, uiCellY, uiRing, uiLayerUnused);
    uint uiLayer = uiParticleId % g_GlobalCloudAttribs.uiMaxLayers;
    uint uiDirectID = 
        (uiRing * g_GlobalCloudAttribs.uiRingDimension * g_GlobalCloudAttribs.uiRingDimension + 
         uiCellX + uiCellY * g_GlobalCloudAttribs.uiRingDimension ) * g_GlobalCloudAttribs.uiMaxLayers + 
         uiLayer;
    g_SunLightAttenuationRW[uiDirectID] = float4(f2CloudLiSpTransparency, 0, ParticleAttrs.fDensity);
}

technique11 ComputeParticlesLightingTech
{
    pass
    {
        SetComputeShader( CompileShader( cs_5_0, ComputeParticlesLightingCS() ) );
    }
}

// This shader compute lighting for all valid particles
[numthreads(THREAD_GROUP_SIZE, 1, 1)]
void SmoothParticlesLightingCS(uint3 Gid  : SV_GroupID, 
                               uint3 GTid : SV_GroupThreadID)
{
    uint uiParticleSerialNum = Gid.x * THREAD_GROUP_SIZE + GTid.x;
    uint uiNumValidParticles = g_ValidCellsCounter.Load(0);
    if( uiParticleSerialNum >= uiNumValidParticles )
        return;

    uint uiParticleId = g_ValidParticlesUnorderedList[uiParticleSerialNum];

	bool bIsVisible = g_VisibleParticleFlags[uiParticleId/32] & (1 << (uiParticleId&31));
    if( !bIsVisible )
        return;

    uint uiCellId = uiParticleId / g_GlobalCloudAttribs.uiMaxLayers;
    SCloudCellAttribs CellAttrs = g_CloudCells[uiCellId];
    uint uiCellX, uiCellY, uiRing, uiLayerUnused;
    UnPackParticleIJRing(CellAttrs.uiPackedLocation, uiCellX, uiCellY, uiRing, uiLayerUnused);
    uint uiLayer = uiParticleId % g_GlobalCloudAttribs.uiMaxLayers;
    uint uiDirectID = 
        (uiRing * g_GlobalCloudAttribs.uiRingDimension * g_GlobalCloudAttribs.uiRingDimension + 
         uiCellX + uiCellY * g_GlobalCloudAttribs.uiRingDimension ) * g_GlobalCloudAttribs.uiMaxLayers + 
         uiLayer;

    float4 f4SmoothedAttenuation = 0;
    float fTotalWegith = 0;
#define CHECK_NEIGHB_PARTICLE(i, j, l)\
    {                                                                                   \
        uint uiSrcParticleID = uiDirectID + ((i) + (j) * g_GlobalCloudAttribs.uiRingDimension) * g_GlobalCloudAttribs.uiMaxLayers + l; \
        float4 f4SrcSunLightAtt = g_SunLightAttenuation[uiSrcParticleID];                   \
        if( f4SrcSunLightAtt.a > 0 )                                                       \
        {                                                                               \
            f4SmoothedAttenuation += float4(f4SrcSunLightAtt.rg*f4SrcSunLightAtt.a, 0, f4SrcSunLightAtt.a);   \
        }                                                                               \
    }

    bool bHasLeftNeighb  = uiCellX >= 1;
    bool bHasRightNeighb = uiCellX < g_GlobalCloudAttribs.uiRingDimension-1;
    bool bHasBottomNeighb= uiCellY >= 1;
    bool bHasTopNeighb   = uiCellY < g_GlobalCloudAttribs.uiRingDimension-1;

    CHECK_NEIGHB_PARTICLE(0,0, 0)
    if(bHasLeftNeighb)  CHECK_NEIGHB_PARTICLE(-1,0,  0)
    if(bHasRightNeighb) CHECK_NEIGHB_PARTICLE(+1,0,  0)
    if(bHasBottomNeighb)CHECK_NEIGHB_PARTICLE(0,-1,  0)
    if(bHasTopNeighb)   CHECK_NEIGHB_PARTICLE(0,+1,  0)
    if( uiLayer < CellAttrs.uiNumActiveLayers-1 )
        CHECK_NEIGHB_PARTICLE(0,0, +1)
    if( uiLayer > 1 )
        CHECK_NEIGHB_PARTICLE(0,0, -1)

    f4SmoothedAttenuation.rgb /= max(f4SmoothedAttenuation.a,1e-5);
    g_ParticlesLightingRW[uiParticleId].f4LightAttenuation = float4(f4SmoothedAttenuation.xy,1,1);
}

struct GS_Input
{
    uint uiParticleID : PARTICLE_ID;
};

GS_Input RenderCloudsVS( uint uiParticleID : PARTICLE_ID )
{
    GS_Input Out = { uiParticleID };
    return Out;
}

#ifndef TILING_MODE
#   define TILING_MODE 0
#endif

// This geometry shader generates either bounding volume (TILING_MODE == 0) or bounding
// sprite (TILING_MODE == 1) for each particle
#if TILING_MODE
[maxvertexcount(4)]
#else
[maxvertexcount(10+4+4)]
#endif
void RenderCloudsGS( point GS_Input In[1], inout TriangleStream<PS_Input> Out )
{	
    uint uiParticleId = In[0].uiParticleID;
    SParticleAttribs ParticleAttrs = g_Particles[uiParticleId];

    // Only visible particles are sent for rendering, so there is no need to
    // test visibility here
    //bool bIsValid = g_VisibleParticleFlags[uiParticleId/32] & (1 << (uiParticleId&31));
    //if( !bIsValid )
    //    return;

    uint uiCellID = uiParticleId / g_GlobalCloudAttribs.uiMaxLayers;
    SCloudCellAttribs CellAttrs = g_CloudCells[uiCellID];

    float fCloudThickness = g_GlobalCloudAttribs.fCloudThickness;
    float3 f3Size = GetParticleScales(ParticleAttrs.fSize, CellAttrs.uiNumActiveLayers);
    
    float3 f3Tangent   = CellAttrs.f3Tangent.xyz;
    float3 f3Normal    = CellAttrs.f3Normal.xyz; 
    float3 f3Bitangent = CellAttrs.f3Bitangent.xyz;
    // Construct particle view-projection matrix
    matrix ParticleObjToWorldSpaceMatr = (matrix)0;
    // Start with rotation:
    ParticleObjToWorldSpaceMatr[0].xyz = f3Tangent.xyz;
    ParticleObjToWorldSpaceMatr[1].xyz = f3Normal.xyz;
    ParticleObjToWorldSpaceMatr[2].xyz = f3Bitangent.xyz;
    // Add translation to particle world position
    ParticleObjToWorldSpaceMatr[3].xyzw = float4(ParticleAttrs.f3Pos,1);

#if TILING_MODE
    float2 f2MinXY = +1;
    float2 f2MaxXY = -1;
    matrix mViewProj = g_GlobalCloudAttribs.mParticleTiling;
#else
    PS_Input Outs[8];
    matrix mViewProj = g_CameraAttribs.WorldViewProj;
#endif
    // Multiply with camera view-proj matrix
    matrix ParticleObjToProjSpaceMatr = mul(ParticleObjToWorldSpaceMatr, mViewProj);

    for(int iCorner = 0; iCorner < 8; ++iCorner)
    {
        float4 f4CurrCornerWorldPos;
        f4CurrCornerWorldPos.xyz = f3Size * float3( (iCorner & 0x01) ? +1 : -1, (iCorner & 0x04) ? +1 : -1, (iCorner & 0x02) ? +1 : -1);
        f4CurrCornerWorldPos.w = 1;

        float4 f4CurrCornerPosPS = mul( f4CurrCornerWorldPos, ParticleObjToProjSpaceMatr );

#if TILING_MODE
        if( f4CurrCornerPosPS.w > 0 )
        {
            float2 f2XY = f4CurrCornerPosPS.xy / f4CurrCornerPosPS.w;
            f2MinXY = min(f2MinXY, f2XY);
            f2MaxXY = max(f2MaxXY, f2XY);
        }
#else
        Outs[iCorner].uiParticleID = uiParticleId;
        Outs[iCorner].f4Pos = f4CurrCornerPosPS;
#endif
    }

#if TILING_MODE
    if( all(f2MinXY < f2MaxXY) )
    {
        // Extend aprite by 0.5 dst pixel size to assure conservative rasterization
        float2 f2TileTexSize = float2(g_GlobalCloudAttribs.fTileTexWidth, g_GlobalCloudAttribs.fTileTexHeight);
        f2MinXY -= 1.0f / f2TileTexSize;
        f2MaxXY += 1.0f / f2TileTexSize;
        float2 f2XY[4] = 
        {
            float2(f2MinXY.x, f2MaxXY.y),
            float2(f2MinXY.x, f2MinXY.y),
            float2(f2MaxXY.x, f2MaxXY.y),
            float2(f2MaxXY.x, f2MinXY.y)
        };
        for(int i = 0; i < 4; ++i)
        {
            PS_Input CurrOut;
            CurrOut.uiParticleID = uiParticleId;
            CurrOut.f4Pos = float4(f2XY[i],1,1);
            Out.Append( CurrOut );
        }
    }
#else
    // Generate bounding box faces
    {
        uint Side[] = {0,4,1,5,3,7,2,6,0,4};
        for(int i = 0; i < 10; ++i)
            Out.Append( Outs[ Side[ i ] ] );
    }
    
    {
        Out.RestartStrip();
        uint uiBottomCap[] = {2,0,3,1};
        for(int i=0; i<4; ++i)
            Out.Append( Outs[ uiBottomCap[ i ] ] );
    }

    {
        Out.RestartStrip();
        uint uiTopCap[] = {4,6,5,7};
        for(int i=0; i<4; ++i)
            Out.Append( Outs[ uiTopCap[ i ] ] );
    }
#endif
}

struct SerializeParticlesGSIn
{
    uint uiGroupNum : GROUP_NUM;
};

SerializeParticlesGSIn SerializeVisibleParticlesVS( uint uiVertId : SV_VertexID )
{
    SerializeParticlesGSIn Out = {uiVertId};
    return Out;
}

struct SerializeParticlesGSOut
{
    uint uiParticleSerailNum : PARTICLE_SERIAL_NUM;
};

#ifndef SERIALIZE_PARTICLE_GS_BATCH_SIZE
#   define SERIALIZE_PARTICLE_GS_BATCH_SIZE 32
#endif

// This shader streams out all the visible particles for rendering with DrawAuto()
// Most important is that particle order is preserved
[maxvertexcount(SERIALIZE_PARTICLE_GS_BATCH_SIZE)]
void SerializeVisibleParticlesGS(point SerializeParticlesGSIn In[1], inout PointStream<SerializeParticlesGSOut> RenderList)
{
    uint uiBatchSize = SERIALIZE_PARTICLE_GS_BATCH_SIZE;
    uint uiGroipNum = In[0].uiGroupNum;
    uint uiParticleSerialNum = uiGroipNum*uiBatchSize;
    for(uint i=0; i < uiBatchSize; ++i, ++uiParticleSerialNum)
    {
        uint uiParticleId = g_SortedParticlesOrder[uiParticleSerialNum];
        bool bIsVisible = g_VisibleParticleFlags[uiParticleId/32] & (1 << (uiParticleId&31));
        if( bIsVisible )
        {
            SerializeParticlesGSOut Particle = {uiParticleId};
            RenderList.Append(Particle);
        }
    }
}

// This helper function computes intersection of the view ray with the particle ellipsoid
void IntersectRayWithParticle(const in SParticleAttribs ParticleAttrs,
                              const in SCloudCellAttribs CellAttrs,
                              const in float3 f3CameraPos, 
                              const in float3 f3ViewRay,
                              out float2 f2RayIsecs,
                              out float3 f3EntryPointUSSpace, // Entry point in Unit Sphere (US) space
                              out float3 f3ViewRayUSSpace,    // View ray direction in Unit Sphere (US) space
                              out float3 f3LightDirUSSpace,   // Light direction in Unit Sphere (US) space
                              out float fDistanceToEntryPoint,
                              out float fDistanceToExitPoint)
{
    // Construct local frame matrix
    float3 f3Normal    = CellAttrs.f3Normal.xyz;
    float3 f3Tangent   = CellAttrs.f3Tangent.xyz;
    float3 f3Bitangent = CellAttrs.f3Bitangent.xyz;
    float3x3 f3x3ObjToWorldSpaceRotation = float3x3(f3Tangent, f3Normal, f3Bitangent);
    // World to obj space is inverse of the obj to world space matrix, which is simply transpose
    // for orthogonal matrix:
    float3x3 f3x3WorldToObjSpaceRotation = transpose(f3x3ObjToWorldSpaceRotation);
    
    // Compute camera location and view direction in particle's object space:
    float3 f3CamPosObjSpace = f3CameraPos - ParticleAttrs.f3Pos;
    f3CamPosObjSpace = mul(f3CamPosObjSpace, f3x3WorldToObjSpaceRotation);
    float3 f3ViewRayObjSpace = mul(f3ViewRay, f3x3WorldToObjSpaceRotation );
    float3 f3LightDirObjSpce = mul(-g_LightAttribs.f4DirOnLight.xyz, f3x3WorldToObjSpaceRotation );

    // Compute scales to transform ellipsoid into the unit sphere:
    float3 f3Scale = 1.f / GetParticleScales(ParticleAttrs.fSize, CellAttrs.uiNumActiveLayers);
    
    float3 f3ScaledCamPosObjSpace;
    f3ScaledCamPosObjSpace  = f3CamPosObjSpace*f3Scale;
    f3ViewRayUSSpace = normalize(f3ViewRayObjSpace*f3Scale);
    f3LightDirUSSpace = normalize(f3LightDirObjSpce*f3Scale);
    // Scale camera pos and view dir in obj space and compute intersection with the unit sphere:
    GetRaySphereIntersection(f3ScaledCamPosObjSpace, f3ViewRayUSSpace, 0, 1.f, f2RayIsecs);

    f3EntryPointUSSpace = f3ScaledCamPosObjSpace + f3ViewRayUSSpace*f2RayIsecs.x;

    fDistanceToEntryPoint = length(f3ViewRayUSSpace/f3Scale) * f2RayIsecs.x;
    fDistanceToExitPoint  = length(f3ViewRayUSSpace/f3Scale) * f2RayIsecs.y;
}

void SwapLayers( inout SParticleLayer Layer0,
                 inout SParticleLayer Layer1 )
{
    SParticleLayer Tmp = Layer0;
    Layer0 = Layer1;
    Layer1 = Tmp;
}

void MergeParticleLayers(in SParticleLayer Layer0,
                         in SParticleLayer Layer1,
                         out SParticleLayer fMergedLayer,
                         out float3 f3ColorToCommit,
                         out float fTranspToCommit,
                         uniform bool BackToFront)
{
    if( Layer0.f2MinMaxDist.x < Layer1.f2MinMaxDist.x )
    {
        SwapLayers( Layer0, Layer1 );
    }

    //
    //       Min1                 Max1
    //        [-------Layer 1------]
    //
    //                  [-------Layer 0------]
    //                 Min0                 Max0
    //
    //     --------------------------------------> Distance
    
    float fMinDist0 = Layer0.f2MinMaxDist.x;
    float fMaxDist0 = Layer0.f2MinMaxDist.y;
    float fMinDist1 = Layer1.f2MinMaxDist.x;
    float fMaxDist1 = Layer1.f2MinMaxDist.y;

    float fIsecLen = min(fMaxDist1, fMaxDist0) - fMinDist0;
    if( fIsecLen <= 0  )
    {
        if( BackToFront )
        {
            fTranspToCommit = exp( -Layer0.fOpticalMass  );
            f3ColorToCommit.rgb = Layer0.f3Color * (1 - fTranspToCommit);
            fMergedLayer = Layer1;
        }
        else
        {
            fTranspToCommit = exp( -Layer1.fOpticalMass );
            f3ColorToCommit.rgb = Layer1.f3Color * (1 - fTranspToCommit);
            fMergedLayer = Layer0;
        }
    }
    else
    {
        float fLayer0Ext = max(fMaxDist0 - fMinDist0, 1e-5);
        float fLayer1Ext = max(fMaxDist1 - fMinDist1, 1e-5);
        float fDensity0 = Layer0.fOpticalMass / fLayer0Ext;
        float fDensity1 = Layer1.fOpticalMass / fLayer1Ext;
        
        float fBackDist = fMaxDist0 - fMaxDist1;
        // fBackDist > 0  (fMaxDist0 > fMaxDist1)
        // ------------------------------------------------------> Distance
        //
        //       Min1                         Max1
        //        [-----------Layer 1----------]
        //
        //                        [-----------Layer 0-----------]
        //                       Min0                          Max0
        //                                      
        //        |               |             |               |
        //             Front             Isec          Back
        //

        // fBackDist < 0 (fMaxDist0 < fMaxDist1)
        // ------------------------------------------------------> Distance
        //
        //       Min1                                                         Max1
        //        [-----------Layer 1------------------------------------------]
        //
        //                        [-----------Layer 0-----------]
        //                       Min0                          Max0
        //                                      
        //        |               |                             |               |
        //             Front                    Isec                   Back
        //

        float fBackDensity = fBackDist > 0 ? fDensity0      : fDensity1;
        float3 f3BackColor = fBackDist > 0 ? Layer0.f3Color : Layer1.f3Color;
        fBackDist = fBackDist > 0 ? fBackDist : -fBackDist;
        float fBackTransparency = exp(-fBackDist*fBackDensity);
        f3BackColor = f3BackColor * (1 - fBackTransparency );
        
        float fIsecTransparency = exp( -(fDensity0 + fDensity1) * fIsecLen );
        float3 f3IsecColor = (Layer0.f3Color * fDensity0 + Layer1.f3Color * fDensity1)/max(fDensity0 + fDensity1, 1e-4) * (1 - fIsecTransparency);
        
        float fFrontDist = fMinDist0 - fMinDist1;
        float fFrontTransparency = exp( -fDensity1 * fFrontDist );
        float3 f3FrontColor = Layer1.f3Color * (1 - fFrontTransparency );
        
        float3 f3Color = 0;
        float fNetTransparency =  1;

        if( BackToFront )
        {
            f3ColorToCommit.rgb = f3BackColor;
            fTranspToCommit = fBackTransparency;
        
            float3 f3Color = f3FrontColor + fFrontTransparency * f3IsecColor;

            float fNetTransparency =  fIsecTransparency * fFrontTransparency;
            fMergedLayer.f3Color = f3Color / max(saturate(1 - fNetTransparency), 1e-10);

            fMergedLayer.f2MinMaxDist.x = fMinDist1;              
            fMergedLayer.f2MinMaxDist.y = min(fMaxDist0, fMaxDist1);
            fMergedLayer.fOpticalMass = fDensity1 * fFrontDist + (fDensity1 + fDensity0) * fIsecLen;
        }
        else
        {
            //if( 1 || ForceMerge )
            //{
                f3ColorToCommit.rgb = 0;
                fTranspToCommit = 1;

                float3 f3Color = f3FrontColor + fFrontTransparency*(f3IsecColor + fIsecTransparency * f3BackColor);

                float fNetTransparency =  fFrontTransparency * fIsecTransparency * fBackTransparency;
                fMergedLayer.f3Color = f3Color / max(saturate(1 - fNetTransparency), 1e-10);

                fMergedLayer.f2MinMaxDist.x = fMinDist1;
                fMergedLayer.f2MinMaxDist.y = max(fMaxDist0, fMaxDist1);
                fMergedLayer.fOpticalMass = Layer0.fOpticalMass + Layer1.fOpticalMass;
            //}
            //else
            //{
            //    f3ColorToCommit.rgb = f3FrontColor;
            //    fTranspToCommit = fFrontTransparency;

            //    float3 f3Color = f3IsecColor + fIsecTransparency * f3BackColor;

            //    float fNetTransparency =  fIsecTransparency * fBackTransparency;
            //    fMergedLayer.f3Color = f3Color / max(saturate(1 - fNetTransparency), 1e-10);

            //    fMergedLayer.f2MinMaxDist.x = fMinDist0;              
            //    fMergedLayer.f2MinMaxDist.y = max(fMaxDist0, fMaxDist1);
            //    fMergedLayer.fOpticalMass = fBackDensity * fBackDist + (fDensity1 + fDensity0) * fIsecLen;
            //}
        }
        // Do not output color
    }
}

// This function computes different attributes of a particle which will be used
// for rendering
void ComputeParticleRenderAttribs(const in SParticleAttribs ParticleAttrs,
                            const in SCloudCellAttribs CellAttrs,
                            in float fTime,
                            in float3 f3CameraPos,
                            in float3 f3ViewRay,
                            in float3 f3EntryPointUSSpace, // Ray entry point in unit sphere (US) space
                            in float3 f3ViewRayUSSpace,    // View direction in unit sphere (US) space
                            in float  fIsecLenUSSpace,     // Length of the intersection of the view ray with the unit sphere
                            in float3 f3LightDirUSSpace,   // Light direction in unit sphere (US) space
                            in float fDistanceToExitPoint,
                            in float fDistanceToEntryPoint,
                            out float fCloudMass,
                            out float fTransparency,
                            uniform in bool bAutoLOD
#if !LIGHT_SPACE_PASS
                            , in SCloudParticleLighting ParticleLighting
                            , out float4 f4Color
#endif
                            )
{
    float3 f3EntryPointWS = f3CameraPos + fDistanceToEntryPoint * f3ViewRay;
    float3 f3ExitPointWS  = f3CameraPos + fDistanceToExitPoint * f3ViewRay;

    float4 f4LUTCoords;
    WorldParamsToOpticalDepthLUTCoords(f3EntryPointUSSpace, f3ViewRayUSSpace, f4LUTCoords);
    // Randomly rotate the sphere
    f4LUTCoords.y += ParticleAttrs.fRndAzimuthBias;
    //float4 f4ddxLUTCoords = ddx(f4LUTCoords);
    //float4 f4ddyLUTCoords = ddy(f4LUTCoords);
    //float4 f4PixelDiff = sqrt(f4ddxLUTCoords*f4ddxLUTCoords + f4ddyLUTCoords*f4ddyLUTCoords) * OPTICAL_DEPTH_LUT_DIM.xyzw;
    //float fMaxDiff = max(f4PixelDiff.x, max(f4PixelDiff.y, max(f4PixelDiff.z,f4PixelDiff.w)));
    float fLOD = 0;//log2(max(fMaxDiff,1));

    float2 f2NormalizedDensityAndDist;
    SAMPLE_4D_LUT(g_tex3DParticleDensityLUT, OPTICAL_DEPTH_LUT_DIM, f4LUTCoords, fLOD, f2NormalizedDensityAndDist);

    float3 f3FirstMatterPointWS = f3CameraPos + (fDistanceToEntryPoint + (fDistanceToExitPoint - fDistanceToEntryPoint) * f2NormalizedDensityAndDist.y) * f3ViewRay;
    float3 f3FirstMatterPointUSSpace = f3EntryPointUSSpace + (fIsecLenUSSpace * f2NormalizedDensityAndDist.y) * f3ViewRayUSSpace;

    float3 f3NoiseSamplingPos = f3FirstMatterPointWS;
    //float3 f3Weights = abs(f3FirstMatterPointUSSpace);
    //f3Weights /= dot(f3Weights, 1);
    float fNoisePeriod = 3412;
    //float fNoise = g_tex2DCloudDensity.Sample(samLinearWrap, (f3NoiseSamplingPos.xz)/fNoisePeriod) * f3Weights.y + 
    //               g_tex2DCloudDensity.Sample(samLinearWrap, (f3NoiseSamplingPos.yz)/fNoisePeriod) * f3Weights.x + 
    //               g_tex2DCloudDensity.Sample(samLinearWrap, (f3NoiseSamplingPos.xy)/fNoisePeriod) * f3Weights.z;
    float fNoise = bAutoLOD ? 
                (g_tex3DNoise.Sample(samLinearWrap, f3NoiseSamplingPos/(fNoisePeriod))*2 + g_tex3DNoise.Sample(samLinearWrap, f3NoiseSamplingPos/(fNoisePeriod/3)))/3 : 
                (g_tex3DNoise.SampleLevel(samLinearWrap, f3NoiseSamplingPos/(fNoisePeriod),0)*2 + g_tex3DNoise.SampleLevel(samLinearWrap, f3NoiseSamplingPos/(fNoisePeriod/3), 0))/3;

    float fNoNoiseY = -0.7;
    float fTransition = min(1-fNoNoiseY, 0.5);
    fNoise = lerp(0.5, fNoise, max(saturate( (f3EntryPointUSSpace.y - fNoNoiseY) / fTransition ), 0.2) );

    //f2NormalizedDensityAndDist.x *= fNoise;

    fCloudMass = f2NormalizedDensityAndDist.x * (fDistanceToExitPoint - fDistanceToEntryPoint);
    float fFadeOutDistance = g_GlobalCloudAttribs.fParticleCutOffDist * g_fParticleToFlatMorphRatio;
    float fFadeOutFactor = saturate( (g_GlobalCloudAttribs.fParticleCutOffDist - fDistanceToEntryPoint) /  max(fFadeOutDistance,1) );
    fCloudMass *= fFadeOutFactor * CellAttrs.fMorphFadeout;
    fCloudMass *= ParticleAttrs.fDensity;// * lerp(fNoise, 1, 0.5);


    //fDistanceToEntryPoint = fDistanceToEntryPoint + (fDistanceToExitPoint - fDistanceToEntryPoint) * f2NormalizedDensityAndDist.y;

    fTransparency = exp( -fCloudMass * g_GlobalCloudAttribs.fAttenuationCoeff );
    
#if !LIGHT_SPACE_PASS
    float fMultipleScatteringDensityScale = fFadeOutFactor * ParticleAttrs.fDensity;
    float fSingleScatteringDensityScale = fMultipleScatteringDensityScale * f2NormalizedDensityAndDist.x;

    float4 f4SingleScatteringLUTCoords   = WorldParamsToParticleScatteringLUT(f3EntryPointUSSpace, f3ViewRayUSSpace, f3LightDirUSSpace, true, fSingleScatteringDensityScale);
    float4 f4MultipleScatteringLUTCoords = WorldParamsToParticleScatteringLUT(f3EntryPointUSSpace, f3ViewRayUSSpace, f3LightDirUSSpace, true, fMultipleScatteringDensityScale);
    float fSingleScattering=0, fMultipleScattering = 0;
    SAMPLE_4D_LUT(g_tex3DSingleScatteringInParticleLUT,   SRF_SCATTERING_IN_PARTICLE_LUT_DIM, f4SingleScatteringLUTCoords,   0, fSingleScattering);
    SAMPLE_4D_LUT(g_tex3DMultipleScatteringInParticleLUT, SRF_SCATTERING_IN_PARTICLE_LUT_DIM, f4MultipleScatteringLUTCoords, 0, fMultipleScattering);
    //float2 f2PrecomputedScattering = g_PrecomputedScatteringInParticle.SampleLevel(samLinearWrap, f3ParticleScatteringLUTCoords, 0);
    
    //float3 _f3EntryPointUSSpace, _f3ViewRayUSSpace, _f3LightDirUSSpace;
    //ParticleScatteringLUTToWorldParams(f3ParticleScatteringLUTCoords, _f3EntryPointUSSpace, _f3ViewRayUSSpace, _f3LightDirUSSpace);
    //float3 _f3ParticleScatteringLUTCoords = WorldParamsToParticleScatteringLUT(_f3EntryPointUSSpace, _f3ViewRayUSSpace, _f3LightDirUSSpace);
    //f4Color.rgb = abs(_f3ParticleScatteringLUTCoords - f3ParticleScatteringLUTCoords)*1e+4;
    float fCosTheta = dot(-f3ViewRayUSSpace, f3LightDirUSSpace);
    float PhaseFunc = HGPhaseFunc(fCosTheta);

    fSingleScattering *= PhaseFunc;
    float fSSSStartRadius = 3.f/5.f;
    float fMultipleScatteringSSSStrength = 0.8f;
    float fFirstMatterPtRadius = length(f3FirstMatterPointUSSpace);
    float fSubSrfScattering = saturate( (1 - fFirstMatterPtRadius)/(1 - fSSSStartRadius) );
    fMultipleScattering *= (1 + fMultipleScatteringSSSStrength * fSubSrfScattering);
    float3 f3Ambient = ParticleLighting.f4AmbientLight.rgb;

    fSingleScattering   *= ParticleLighting.f4LightAttenuation.x ;
    fMultipleScattering *= ParticleLighting.f4LightAttenuation.y * (0.5+1.0*fNoise);

    float fAmbientSSSStrength = (1-fNoise)*0.5;//0.3;
    f3Ambient.rgb *= lerp(1, fSubSrfScattering, fAmbientSSSStrength);
    f4Color.rgb = (f3Ambient.rgb + (fSingleScattering + fMultipleScattering) * ParticleLighting.f4SunLight.rgb ) * PI;
    f4Color.a = fTransparency;
#endif    

}


RWTexture2D<int> g_tex2DFirstKnotRW : register( u0 );
RWStructuredBuffer<SParticleListKnot> g_TiledParticlesListRW : register( u1 );

// This shader performs tiling
// It uses pixel sync extension to add particles to the list in the same order
// they are submitted for rendering
void PerformTilingPS(PS_Input In)
{
    IntelExt_Init();

    uint2 ui2TileXY = In.f4Pos.xy;
    int iNewInd = g_TiledParticlesListRW.IncrementCounter();
    int iOldInd;
    IntelExt_BeginPixelShaderOrdering();
    InterlockedExchange(g_tex2DFirstKnotRW[ui2TileXY], iNewInd, iOldInd);
    g_TiledParticlesListRW[iNewInd].uiParticleID = In.uiParticleID;
    g_TiledParticlesListRW[iNewInd].iNextKnotInd = iOldInd;
}

// This function gets conservative minimal depth (which corresponds to the furthest surface)
// When rendering clouds in lower resolution, it is essentially important to assure that the 
// furthest depth is taken, which is required for proper bilateral upscaling
float GetConservativeScreenDepth(in float2 f2UV)
{
    float fDepth;

#if BACK_BUFFER_DOWNSCALE_FACTOR == 1
    fDepth = g_tex2DDepthBuffer.SampleLevel(samLinearClamp, f2UV, 0 );
#else
    
    float2 f2GatherUV0 = f2UV;
#   if BACK_BUFFER_DOWNSCALE_FACTOR > 2
        float2 f2DepthBufferSize = float2(g_GlobalCloudAttribs.fBackBufferWidth, g_GlobalCloudAttribs.fBackBufferHeight);
        f2GatherUV0 -= (BACK_BUFFER_DOWNSCALE_FACTOR/2.f-1.f)/f2DepthBufferSize.xy;
#   endif

    fDepth = 1;
    [unroll]
    for(int i=0; i < BACK_BUFFER_DOWNSCALE_FACTOR/2; ++i)
        [unroll]
        for(int j=0; j < BACK_BUFFER_DOWNSCALE_FACTOR/2; ++j)
        {
            float4 f4Depths = g_tex2DDepthBuffer.Gather(samLinearClamp, f2GatherUV0, 2*int2(i,j) );
            fDepth = min(fDepth,f4Depths.x);
            fDepth = min(fDepth,f4Depths.y);
            fDepth = min(fDepth,f4Depths.z);
            fDepth = min(fDepth,f4Depths.w);
        }
#endif
    return fDepth;
}

void RenderCloudsPS1(PS_Input In,
                    out float fTransparency : SV_Target0,
                    out float2 f2MinMaxDist : SV_Target1,
                    out float4 f4Color : SV_Target2)

{
    fTransparency = 0.9;
    f2MinMaxDist = 0;
    f4Color = 0.2;
}

// This shader renders particles
void RenderCloudsPS(PS_Input In,
                    out float fTransparency : SV_Target0

#if !VOLUMETRIC_BLENDING || LIGHT_SPACE_PASS
                    , out float fDistToCloud : SV_Target1
#endif

#if !LIGHT_SPACE_PASS
                    , out float4 f4Color : SV_Target2
#endif
                    )
{
#if !LIGHT_SPACE_PASS && PS_ORDERING_AVAILABLE
    IntelExt_Init();
#endif

    SParticleAttribs ParticleAttrs = g_Particles[In.uiParticleID];
    SCloudCellAttribs CellAttribs = g_CloudCells[In.uiParticleID / g_GlobalCloudAttribs.uiMaxLayers];

#if !LIGHT_SPACE_PASS
    SCloudParticleLighting ParticleLighting = g_ParticlesLighting[In.uiParticleID];
#endif
    float fTime = g_fTimeScale*g_GlobalCloudAttribs.fTime;

    float3 f3CameraPos, f3ViewRay;
#if LIGHT_SPACE_PASS
    // For directional light source, we should use position on the near clip plane instead of
    // camera location as a ray start point
    float2 f2PosPS = UVToProj( (In.f4Pos.xy / g_GlobalCloudAttribs.f2LiSpCloudDensityDim.xy) );
    float4 f4PosOnNearClipPlaneWS = mul( float4(f2PosPS.xy,1,1), g_CameraAttribs.mViewProjInv );
    f3CameraPos = f4PosOnNearClipPlaneWS.xyz/f4PosOnNearClipPlaneWS.w;
    
    //f4PosOnNearClipPlaneWS = mul( float4(f2PosPS.xy,1e-4,1), g_CameraAttribs.mViewProjInv );
    //f3CameraPos = f4PosOnNearClipPlaneWS.xyz/f4PosOnNearClipPlaneWS.w;
    float4 f4PosOnFarClipPlaneWS = mul( float4(f2PosPS.xy,0,1), g_CameraAttribs.mViewProjInv );
    f4PosOnFarClipPlaneWS.xyz = f4PosOnFarClipPlaneWS.xyz/f4PosOnFarClipPlaneWS.w;
    f3ViewRay = normalize(f4PosOnFarClipPlaneWS.xyz - f4PosOnNearClipPlaneWS.xyz);
#else
    f3CameraPos = g_CameraAttribs.f4CameraPos.xyz;
    //f3ViewRay = normalize(In.f3ViewRay);
    float2 f2ScreenDim = float2(g_GlobalCloudAttribs.fDownscaledBackBufferWidth, g_GlobalCloudAttribs.fDownscaledBackBufferHeight);
    float2 f2PosPS = UVToProj( In.f4Pos.xy / f2ScreenDim );
    float fDepth = GetConservativeScreenDepth( ProjToUV(f2PosPS.xy) );
    float4 f4ReconstructedPosWS = mul( float4(f2PosPS.xy,fDepth,1.0), g_CameraAttribs.mViewProjInv );
    float3 f3WorldPos = f4ReconstructedPosWS.xyz / f4ReconstructedPosWS.w;

    // Compute view ray
    f3ViewRay = f3WorldPos - f3CameraPos;
    float fRayLength = length(f3ViewRay);
    f3ViewRay /= fRayLength;

#endif

    // Intersect view ray with the particle
    float2 f2RayIsecs;
    float fDistanceToEntryPoint, fDistanceToExitPoint;
    float3 f3EntryPointUSSpace, f3ViewRayUSSpace, f3LightDirUSSpace;
    IntersectRayWithParticle(ParticleAttrs, CellAttribs, f3CameraPos,  f3ViewRay,
                             f2RayIsecs, f3EntryPointUSSpace, f3ViewRayUSSpace,
                             f3LightDirUSSpace,
                             fDistanceToEntryPoint, fDistanceToExitPoint);
   

#if LIGHT_SPACE_PASS
    if( all(f2RayIsecs == NO_INTERSECTIONS) )
        discard;
#else
    if( f2RayIsecs.y < 0 || fRayLength < fDistanceToEntryPoint )
        discard;
    fDistanceToExitPoint = min(fDistanceToExitPoint, fRayLength);
#endif
    float fCloudMass;
    float fIsecLenUSSpace = f2RayIsecs.y - f2RayIsecs.x;
    // Compute particle rendering attributes
    ComputeParticleRenderAttribs(ParticleAttrs, CellAttribs,
                            fTime,
                            f3CameraPos,
                            f3ViewRay,
                            f3EntryPointUSSpace, 
                            f3ViewRayUSSpace,
                            fIsecLenUSSpace,
                            f3LightDirUSSpace,
                            fDistanceToExitPoint,
                            fDistanceToEntryPoint,
                            fCloudMass,
                            fTransparency,
                            true
#if !LIGHT_SPACE_PASS
                            , ParticleLighting
                            , f4Color
#endif
                            );

#if !LIGHT_SPACE_PASS
    #if VOLUMETRIC_BLENDING
        uint2 ui2PixelIJ = In.f4Pos.xy;
        uint uiLayerDataInd = (ui2PixelIJ.x + ui2PixelIJ.y * g_GlobalCloudAttribs.uiDownscaledBackBufferWidth) * NUM_PARTICLE_LAYERS;
        SParticleLayer Layers[NUM_PARTICLE_LAYERS+1];

        Layers[NUM_PARTICLE_LAYERS].f2MinMaxDist = float2(fDistanceToEntryPoint, fDistanceToExitPoint);
        Layers[NUM_PARTICLE_LAYERS].fOpticalMass = fCloudMass * g_GlobalCloudAttribs.fAttenuationCoeff;
        Layers[NUM_PARTICLE_LAYERS].f3Color = f4Color.rgb;

        f4Color = float4(0,0,0,1);
        fTransparency = 1;

        IntelExt_BeginPixelShaderOrdering();

        [unroll]
        for(int iLayer=0; iLayer < NUM_PARTICLE_LAYERS; ++iLayer)
            Layers[iLayer] = g_rwbufParticleLayers[uiLayerDataInd + iLayer];

        // Sort layers
        //for(int i = 0; i < NUM_PARTICLE_LAYERS; ++i)
        //    for(int j = i+1; j < NUM_PARTICLE_LAYERS; ++j)
        //        if( Layers[i].f2MinMaxDist.x < Layers[j].f2MinMaxDist.x )
        //            SwapLayers( Layers[i], Layers[j] );


        // Merge two furthest layers
        SParticleLayer MergedLayer;
        MergeParticleLayers(Layers[0], Layers[1], MergedLayer, f4Color.rgb, fTransparency, true);
        f4Color.a = fTransparency;

        Layers[1] = MergedLayer;

        // Store updated layers
        [unroll]
        for(iLayer = 0; iLayer < NUM_PARTICLE_LAYERS; ++iLayer)
            g_rwbufParticleLayers[uiLayerDataInd + iLayer] = Layers[iLayer+1];
    #else
        f4Color.rgb *= 1-fTransparency;
    #endif
#endif


#if !VOLUMETRIC_BLENDING || LIGHT_SPACE_PASS
    fDistToCloud = fTransparency < 0.9 ? fDistanceToEntryPoint : +FLT_MAX;
#endif
}

void ApplyParticleLayersPS( SScreenSizeQuadVSOutput In,
					        out float fTransparency : SV_Target0,
	                        out float2 f2MinMaxDist : SV_Target1,
					        out float4 f4Color : SV_Target2)
{
    uint2 ui2PixelIJ = In.m_f4Pos.xy;
    uint uiLayerDataInd = (ui2PixelIJ.x + ui2PixelIJ.y * g_GlobalCloudAttribs.uiDownscaledBackBufferWidth) * NUM_PARTICLE_LAYERS;
    SParticleLayer Layers[NUM_PARTICLE_LAYERS];
    [unroll]
    for(int iLayer = 0; iLayer < NUM_PARTICLE_LAYERS; ++iLayer)
        Layers[iLayer] = g_bufParticleLayers[uiLayerDataInd + iLayer];
    
    //for(iLastLayer = NUM_PARTICLE_LAYERS-1; iLastLayer >= 0; --iLastLayer)
    //    if( Layers[iLastLayer].f2MinMaxDist.x < FLT_MAX )
    //        break;

    //if( iLastLayer < 0 )
    //    discard;

    int iLastLayer = NUM_PARTICLE_LAYERS-1;
    SParticleLayer LastLayer = Layers[0];

    f4Color = 0;
    fTransparency = 1;
    [unroll]
    for(iLayer = 1; iLayer <= iLastLayer; ++iLayer)
    {
        float3 f3MergedColor;
        float fMergedTransparency;
        SParticleLayer MergedLayer;
        MergeParticleLayers(LastLayer, Layers[iLayer], MergedLayer, f3MergedColor.rgb, fMergedTransparency, true);

        fTransparency *= fMergedTransparency;
        f4Color.rgb = f4Color.rgb * fMergedTransparency + f3MergedColor;
        LastLayer = MergedLayer;
    }
    float fLastTransparency = exp( -LastLayer.fOpticalMass );
    f4Color.rgb = f4Color.rgb * fLastTransparency + LastLayer.f3Color * (1-fLastTransparency);
    fTransparency *= fLastTransparency;
    f4Color.a = fTransparency;
    f2MinMaxDist = LastLayer.f2MinMaxDist.x;
}

#ifndef TILE_SIZE
#   define TILE_SIZE 8
#endif

float GetLayerIsec(SParticleLayer Layer0, SParticleLayer Layer1)
{
    return min(Layer0.f2MinMaxDist.y, Layer1.f2MinMaxDist.y) - max(Layer0.f2MinMaxDist.x, Layer1.f2MinMaxDist.x);
}

// This shader renders particles using deffered tiled approach
void RenderCloudsTiledPS( SScreenSizeQuadVSOutput In,
					      out float fTransparency : SV_Target0,
	                      out float2 f2MinMaxDist : SV_Target1,
					      out float4 f4Color : SV_Target2)
{
    float3 f3CameraPos = g_CameraAttribs.f4CameraPos.xyz;
    float fTime = g_fTimeScale*g_GlobalCloudAttribs.fTime;

    // Reconstruct world space position
    float fDepth = GetConservativeScreenDepth( ProjToUV(In.m_f2PosPS.xy) );
    float4 f4ReconstructedPosWS = mul( float4(In.m_f2PosPS.xy,fDepth,1), g_CameraAttribs.mViewProjInv );
    float3 f3WorldPos = f4ReconstructedPosWS.xyz / f4ReconstructedPosWS.w;

    // Compute view ray
    float3 f3ViewDir = f3WorldPos - f3CameraPos;
    float fRayLength = length(f3ViewDir);
    f3ViewDir /= fRayLength;
    
    // Get tile location
    uint2 ui2PixelXY = floor(In.m_f4Pos.xy);
    uint2 uiTileXY = ui2PixelXY / TILE_SIZE;
    // Get first tile index
    int iCurrKnotInd = g_tex2DFirstKnot.Load( uint3(uiTileXY, 0) );
    fTransparency = 1;
    f2MinMaxDist = float2(+FLT_MAX, +FLT_MAX);
    f4Color = 0;

    SParticleLayer Layers[NUM_PARTICLE_LAYERS+1];
    for(int i=0; i < NUM_PARTICLE_LAYERS; ++i)
    {
        Layers[i].f3Color = 0;
        Layers[i].f2MinMaxDist = -1;
        Layers[i].fOpticalMass = 0;
    }

    //int iter = 0;
    [loop]
    while( iCurrKnotInd >= 0 )
    {
        //++iter;
        //if(iter > 10)break;

        SParticleListKnot CurrKnot = g_TiledParticlesList[iCurrKnotInd];
        SParticleAttribs ParticleAttrs = g_Particles[CurrKnot.uiParticleID];
        uint uiCellId = CurrKnot.uiParticleID/g_GlobalCloudAttribs.uiMaxLayers;
        SCloudCellAttribs CellAttribs = g_CloudCells[uiCellId];

        float2 f2RayIsecs;
        float fDistanceToEntryPoint, fDistanceToExitPoint;
        float3 f3EntryPointUSSpace, f3ViewRayUSSpace, f3LightDirUSSpace;
        IntersectRayWithParticle(ParticleAttrs, CellAttribs, f3CameraPos,  f3ViewDir,
                                 f2RayIsecs, f3EntryPointUSSpace, f3ViewRayUSSpace,
                                 f3LightDirUSSpace,
                                 fDistanceToEntryPoint, fDistanceToExitPoint);
        
        if( f2RayIsecs.y > 0 && fRayLength > fDistanceToEntryPoint )
        {
            fDistanceToExitPoint = min(fDistanceToExitPoint, fRayLength);

            float fCurrCloudMass;
            float fCurrTransparency = 0.8;
            float4 f4CurrColor = 0.1;
            float fIsecLenUSSpace = f2RayIsecs.y - f2RayIsecs.x;
            ComputeParticleRenderAttribs(ParticleAttrs, CellAttribs,
                                    fTime,
                                    f3CameraPos, f3ViewDir,
                                    f3EntryPointUSSpace, 
                                    f3ViewRayUSSpace,
                                    fIsecLenUSSpace,
                                    f3LightDirUSSpace,
                                    fDistanceToExitPoint,
                                    fDistanceToEntryPoint,
                                    fCurrCloudMass,
                                    fCurrTransparency,
                                    false
        #if !LIGHT_SPACE_PASS
                                    , g_ParticlesLighting[CurrKnot.uiParticleID]
                                    , f4CurrColor
        #endif
                                    );

            float2 f2CurrMinMaxDist = float2(fDistanceToEntryPoint, -fDistanceToExitPoint);
            f2MinMaxDist = min(f2MinMaxDist, f2CurrMinMaxDist);

            Layers[NUM_PARTICLE_LAYERS].f2MinMaxDist = float2(fDistanceToEntryPoint, fDistanceToExitPoint);
            Layers[NUM_PARTICLE_LAYERS].fOpticalMass = fCurrCloudMass * g_GlobalCloudAttribs.fAttenuationCoeff;
            Layers[NUM_PARTICLE_LAYERS].f3Color = f4CurrColor.rgb;
            float fMergedTransparency = 1;
            float4 f4MergedColor = float4(0,0,0,1);
            
            MergeParticleLayers(Layers[1], Layers[0], Layers[0], f4MergedColor.rgb, fMergedTransparency, false);
                        
            [unroll]
            for(int i = 1; i < NUM_PARTICLE_LAYERS; ++i)
                Layers[i] = Layers[i+1];
            
            f4Color += f4MergedColor * fTransparency;
            fTransparency *= fMergedTransparency;

            float fUncommitedLayersOpticalMass = 0;
            for(int i = 0; i < NUM_PARTICLE_LAYERS; ++i)
                fUncommitedLayersOpticalMass += Layers[i].fOpticalMass;
            float fUncommitedLayersTransparency = 1;//exp( -fUncommitedLayersOpticalMass );

            if( fTransparency * fUncommitedLayersTransparency < g_fTransparencyThreshold )
                break;
        }

        iCurrKnotInd = CurrKnot.iNextKnotInd;
    }

    //if( fTransparency > g_fTransparencyThreshold )
    {
        float3 f3MergedColor;
        float fMergedTransparency;
        SParticleLayer LayerToCommit;

        [unroll]
        for(int iLayer = 1; iLayer <= NUM_PARTICLE_LAYERS-1; ++iLayer)
        {
            MergeParticleLayers(Layers[iLayer], Layers[0], Layers[0], f3MergedColor.rgb, fMergedTransparency, false);
            f4Color.rgb += f3MergedColor * fTransparency;
            fTransparency *= fMergedTransparency;
        }
        SParticleLayer LastLayer = Layers[0];

        float fLastTransparency = exp( -LastLayer.fOpticalMass );
        f4Color.rgb += fTransparency * LastLayer.f3Color * (1-fLastTransparency);
        fTransparency *= fLastTransparency;
    }

    if( fTransparency > 1-1e-3 )
        discard;
    
    fTransparency = saturate((fTransparency - g_fTransparencyThreshold) / (1-g_fTransparencyThreshold));

    f4Color.a = fTransparency;
}

technique11 RenderCloudsTech
{
    pass
    {
        SetVertexShader( CompileShader( vs_5_0, RenderCloudsVS() ) );
        SetGeometryShader( CompileShader( gs_5_0, RenderCloudsGS() ) );
        SetPixelShader( CompileShader( ps_5_0, RenderCloudsPS() ) );
    }
}

// This shader combines cloud buffers with the back buffer
float4 CombineWithBackBufferPS(SScreenSizeQuadVSOutput In) : SV_Target
{
    float3 f3BackgroundColor = g_tex2DColorBuffer.Load( int3(In.m_f4Pos.xy,0) ).rgb;
    float3 f3CloudColor = g_tex2DScrSpaceCloudColor.Load( int3(In.m_f4Pos.xy,0) ).rgb;
    float fTransparency = g_tex2DScrSpaceCloudTransparency.Load( int3(In.m_f4Pos.xy,0) );
    return float4(f3BackgroundColor * fTransparency + f3CloudColor, 1);
}

// This shader renders cloud depth into the shadow map
void RenderCloudDepthToShadowMapPS(SScreenSizeQuadVSOutput In,
                                   out float fDepth : SV_Depth,
                                   out float4 f4DummyColor : SV_Target)
{
    // TODO: transparency map / when shadow map dimension ration is <= 1:4, wider filtering is required
    float fLiSpCloudTransparency = g_tex2DLiSpCloudTransparency.Sample(samLinearClamp, float3(ProjToUV(In.m_f2PosPS.xy), g_GlobalCloudAttribs.f4Parameter.x) );
    float fMaxDepth = g_tex2DLiSpCloudMinMaxDepth.Sample(samLinearClamp, float3(ProjToUV(In.m_f2PosPS.xy), g_GlobalCloudAttribs.f4Parameter.x) ).y;
    fDepth = ( fLiSpCloudTransparency < 0.1 ) ? fMaxDepth : 0.f;
    f4DummyColor = 0;
}
