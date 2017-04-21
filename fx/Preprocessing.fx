/////////////////////////////////////////////////////////////////////////////////////////////
// Copyright 2017 Intel Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");// you may not use this file except in compliance with the License.// You may obtain a copy of the License at//// http://www.apache.org/licenses/LICENSE-2.0//// Unless required by applicable law or agreed to in writing, software// distributed under the License is distributed on an "AS IS" BASIS,// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.// See the License for the specific language governing permissions and// limitations under the License.
/////////////////////////////////////////////////////////////////////////////////////////////
#include "Common.fxh"
#include "CloudsCommon.fxh"

#ifndef DENSITY_GENERATION_METHOD
#   define DENSITY_GENERATION_METHOD 0
#endif

Texture3D<float>       g_tex3DNoise                 : register( t0 );
Texture3D<float>       g_tex3DSingleScattering      : register( t0 );
Texture3D<float>       g_tex3DMultipleScattering    : register( t1 );
Texture3D<float>       g_tex3DPrevSctrOrder         : register( t0 );
Texture3D<float>       g_tex3DGatheredScattering    : register( t0 );
Texture3D<float>       g_tex3DSingleScatteringLUT    : register( t0 );
Texture3D<float>       g_tex3DMultipleScatteringLUT  : register( t1 );

cbuffer cbPostProcessingAttribs : register( b0 )
{
    SGlobalCloudAttribs g_GlobalCloudAttribs;
};

SamplerState samLinearWrap : register( s1 );
SamplerState samPointWrap : register( s2 );

#define NUM_INTEGRATION_STEPS 64

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

float GetRandomDensity(in float3 f3Pos, float fStartFreq, int iNumOctaves = 3, float fAmplitudeScale = 0.6)
{
    float fNoise = 0;
    float fAmplitude = 1;
    float fFreq = fStartFreq;
    for(int i=0; i < iNumOctaves; ++i)
    {
        fNoise += (g_tex3DNoise.SampleLevel(samLinearWrap, f3Pos*fFreq, 0) - 0.5) * fAmplitude;
        fFreq *= 1.7;
        fAmplitude *= fAmplitudeScale;
    }
    return fNoise;
}

float GetPyroSphereDensity(float3 f3CurrPos)
{
    float fDistToCenter = length(f3CurrPos);
    float3 f3NormalizedPos = f3CurrPos / fDistToCenter;
    float fNoise = GetRandomDensity(f3NormalizedPos, 0.15, 4, 0.8);
    float fDensity = fDistToCenter + 0.35*fNoise < 0.8 ? 1 : 0;
    return fDensity;
}

float GetMetabolDensity(in float r)
{
    float r2 = r*r;
    float r4 = r2*r2;
    float r6 = r4*r2;
    return saturate(-4.0/9.0 * r6 + 17.0/9.0 * r4 - 22.0/9.0 * r2 + 1);
}

// This shader computes level 0 of the maximum density mip map
float2 PrecomputeOpticalDepthPS(SScreenSizeQuadVSOutput In) : SV_Target
{
    float3 f3NormalizedStartPos, f3RayDir;
    OpticalDepthLUTCoordsToWorldParams( float4(ProjToUV(In.m_f2PosPS), g_GlobalCloudAttribs.f4Parameter.xy), f3NormalizedStartPos, f3RayDir );
    
    // Intersect view ray with the unit sphere:
    float2 f2RayIsecs;
    // f3NormalizedStartPos  is located exactly on the surface; slightly move start pos inside the sphere
    // to avoid precision issues
    GetRaySphereIntersection(f3NormalizedStartPos + f3RayDir*1e-4, f3RayDir, 0, 1.f, f2RayIsecs);
    
    if( f2RayIsecs.x > f2RayIsecs.y )
        return 0;

    float3 f3EndPos = f3NormalizedStartPos + f3RayDir * f2RayIsecs.y;
    float fNumSteps = NUM_INTEGRATION_STEPS;
    float3 f3Step = (f3EndPos - f3NormalizedStartPos) / fNumSteps;
    float fTotalDensity = 0;
    float fDistToFirstMatter = -1;
    for(float fStepNum=0.5; fStepNum < fNumSteps; ++fStepNum)
    {
        float3 f3CurrPos = f3NormalizedStartPos + f3Step * fStepNum;
        
        float fDistToCenter = length(f3CurrPos);
        float fMetabolDensity = GetMetabolDensity(fDistToCenter);
        float fDensity = 1;
#if DENSITY_GENERATION_METHOD == 0
        fDensity = saturate( 1.0*saturate(fMetabolDensity) + 1*pow(fMetabolDensity,0.5)*(GetRandomDensity(f3CurrPos, 0.15, 4, 0.7 )) );
#elif DENSITY_GENERATION_METHOD == 1
        fDensity = 1.0*saturate(fMetabolDensity) + 1.0*pow(fMetabolDensity,0.5)*(GetRandomDensity(f3CurrPos, 0.1,4,0.8)) > 0.1 ? 1 : 0;
        //fDensity = saturate(fMetabolDensity) - 2*pow(fMetabolDensity,0.5)*GetRandomDensity(f3CurrPos, 0.2, 4, 0.7) > 0.05 ? 1 : 0;
#elif DENSITY_GENERATION_METHOD == 2
        fDensity = GetPyroSphereDensity(f3CurrPos);
#endif

        if( fDensity > 0.05 && fDistToFirstMatter < 0 )
            fDistToFirstMatter = fStepNum / fNumSteps;
        
        fTotalDensity += fDensity;
    }

    if( fDistToFirstMatter < 0 )
        fDistToFirstMatter = 1;

    //float fOcclusion = 0;
    //if( !all(f3FirstMatterPoint ==-1) )
    //{
    //    //float3 f3LocalX, f3LocalY, f3LocalZ;
    //    //ConstructLocalFrameXYZ(-normalize(f3FirstMatterPoint), float3(0,0,1), f3LocalX, f3LocalY, f3LocalZ);

    //    float fNumZenithAngles = 8;
    //    float fNumAzimuthAngles = 16;
    //    float fTotalSolidAngle = 0, fOccludedAngle = 0;
    //    for(float fZenith = 0.5; fZenith < fNumZenithAngles; ++fZenith)
    //        for(float fAzimuth = 0.5; fAzimuth < fNumAzimuthAngles; ++fAzimuth)
    //        {
    //            float fZenithAngle = (fZenith/fNumZenithAngles) * PI;
    //            float fAzimuthAngle = (fAzimuth/fNumAzimuthAngles - 0.5) * (2*PI);
    //            float3 f3CurrDir = ZenithAzimuthAngleToDirection(fZenithAngle, fAzimuthAngle);
    //            float3 f3DisplacedPoint = f3FirstMatterPoint + f3CurrDir * 0.05;
    //            float fDensity = GetPyroSphereDensity(f3DisplacedPoint);
    //            float fDiffSolidAngle = sin(fZenithAngle) * (PI / fNumZenithAngles) * (2*PI / fNumAzimuthAngles);
    //            if( fDensity > 0 )
    //                fOccludedAngle+=fDiffSolidAngle;
    //            fTotalSolidAngle += fDiffSolidAngle;
    //        }
    //    fOcclusion = fOccludedAngle / fTotalSolidAngle;
    //}

    return float2(fTotalDensity / fNumSteps, fDistToFirstMatter);
}


float PrecomputeSingleSctrPS(SScreenSizeQuadVSOutput In) : SV_Target
{
    float4 f4LUTCoords = float4(ProjToUV(In.m_f2PosPS), g_GlobalCloudAttribs.f4Parameter.xy);

    float3 f3EntryPointUSSpace, f3ViewRayUSSpace, f3LightDirUSSpace;
    float fDensityScale;
    ParticleScatteringLUTToWorldParams(f4LUTCoords, g_GlobalCloudAttribs.f4Parameter.z, f3EntryPointUSSpace, f3ViewRayUSSpace, f3LightDirUSSpace, false, fDensityScale);

    // Intersect view ray with the unit sphere:
    float2 f2RayIsecs;
    // f3NormalizedStartPos  is located exactly on the surface; slightly move the start pos inside the sphere
    // to avoid precision issues
    float3 f3BiasedEntryPoint = f3EntryPointUSSpace + f3ViewRayUSSpace*1e-4;
    GetRaySphereIntersection(f3BiasedEntryPoint, f3ViewRayUSSpace, 0, 1.f, f2RayIsecs);
    if( f2RayIsecs.y < f2RayIsecs.x )
        return 0;
    float3 f3EndPos = f3BiasedEntryPoint + f3ViewRayUSSpace * f2RayIsecs.y;

    float fNumSteps = NUM_INTEGRATION_STEPS;
    float3 f3Step = (f3EndPos - f3EntryPointUSSpace) / fNumSteps;
    float fStepLen = length(f3Step);
    float fCloudMassToCamera = 0;
    float fParticleRadius = GetParticleSize(GetCloudRingWorldStep(0, g_GlobalCloudAttribs));
    float fInscattering = 0;
    for(float fStepNum=0.5; fStepNum < fNumSteps; ++fStepNum)
    {
        float3 f3CurrPos = f3EntryPointUSSpace + f3Step * fStepNum;
        //float fDistToCenter = length(f3CurrPos);
        //float fMetabolDensity = GetMetabolDensity(fDistToCenter);
        float fDensity = 1;//saturate(fMetabolDensity);
        fDensity *= fDensityScale;

        float fCloudMassToLight = 0;
        GetRaySphereIntersection(f3CurrPos, f3LightDirUSSpace, 0, 1.f, f2RayIsecs);
        if( f2RayIsecs.y > f2RayIsecs.x )
        {
            //float3 f3LightEntryPoint = f3CurrPos + f3LightDirUSSpace * f2RayIsecs.x;
            //const float fNumStepsToLight = 8;
            //float3 f3StepToLight = (f3LightEntryPoint - f3CurrPos) / fNumStepsToLight;
            //for(float fStepNum2=0.5; fStepNum2 < fNumStepsToLight; ++fStepNum2)
            //{
            //    float3 f3CurrPos2 = f3CurrPos + f3StepToLight * fStepNum2;
            //    float fDensity2 = saturate(GetMetabolDensity(f3CurrPos2));
            //    fCloudMassToLight += fDensity2; 
            //}
            //fCloudMassToLight *= length(f3StepToLight) * fParticleRadius;
            //fCloudMassToLight = length(f3LightEntryPoint - f3CurrPos) * fParticleRadius;
            //fCloudMassToLight /= 200;
            fCloudMassToLight = abs(f2RayIsecs.x) * fParticleRadius;
        }

        float fTotalLightAttenuation = exp( -g_GlobalCloudAttribs.fAttenuationCoeff * (fCloudMassToLight + fCloudMassToCamera) );
        
        //float3 f3NormalizedPos = normalize(f3CurrPos);
        //float fNoise = GetRandomDensity(f3NormalizedPos);
        //float fDensity = fDistToCenter + 0.3*fNoise < 0.7 ? 1 : 0;
        fInscattering += fTotalLightAttenuation * fDensity * g_GlobalCloudAttribs.fScatteringCoeff;
        //fInscattering += fLightAttenuation/fNumSteps;//fCloudMassToLight/50000;
        fCloudMassToCamera += fDensity * fStepLen * fParticleRadius;
    }

    return fInscattering * fStepLen * fParticleRadius;
}


float GatherScatteringPS(SScreenSizeQuadVSOutput In) : SV_Target
{
    float4 f4LUTCoords = float4(ProjToUV(In.m_f2PosPS), g_GlobalCloudAttribs.f4Parameter.xy);

    float3 f3PosUSSpace, f3ViewRayUSSpace, f3LightDirUSSpace;
    float fDensityScale;
    ParticleScatteringLUTToWorldParams(f4LUTCoords, g_GlobalCloudAttribs.f4Parameter.z, f3PosUSSpace, f3ViewRayUSSpace, f3LightDirUSSpace, false, fDensityScale);

    float3 f3LocalX, f3LocalY, f3LocalZ;
    ConstructLocalFrameXYZ(-normalize(f3PosUSSpace), f3LightDirUSSpace, f3LocalX, f3LocalY, f3LocalZ);

    float fGatheredScattering = 0;
    float fTotalSolidAngle = 0;
    const float fNumZenithAngles = VOL_SCATTERING_IN_PARTICLE_LUT_DIM.z;
    const float fNumAzimuthAngles = VOL_SCATTERING_IN_PARTICLE_LUT_DIM.y;
    const float fZenithSpan = PI;
    const float fAzimuthSpan = 2*PI;
    for(float ZenithAngleNum = 0.5; ZenithAngleNum < fNumZenithAngles; ++ZenithAngleNum)
        for(float AzimuthAngleNum = 0.5; AzimuthAngleNum < fNumAzimuthAngles; ++AzimuthAngleNum)
        {
            float ZenithAngle = ZenithAngleNum/fNumZenithAngles * fZenithSpan;
            float AzimuthAngle = (AzimuthAngleNum/fNumAzimuthAngles - 0.5) * fAzimuthSpan;
            float3 f3CurrDir = GetDirectionInLocalFrameXYZ(f3LocalX, f3LocalY, f3LocalZ, ZenithAngle, AzimuthAngle);
            float4 f4CurrDirLUTCoords = WorldParamsToParticleScatteringLUT(f3PosUSSpace, f3CurrDir, f3LightDirUSSpace, false, 1);
            float fCurrDirScattering = 0;
            SAMPLE_4D_LUT(g_tex3DPrevSctrOrder, VOL_SCATTERING_IN_PARTICLE_LUT_DIM, f4CurrDirLUTCoords, 0, fCurrDirScattering);
            if( g_GlobalCloudAttribs.f4Parameter.w == 1 )
            {
                fCurrDirScattering *= HGPhaseFunc( dot(-f3CurrDir, f3LightDirUSSpace) );
            }
            fCurrDirScattering *= HGPhaseFunc( dot(f3CurrDir, f3ViewRayUSSpace), 0.7 );

            float fdZenithAngle = fZenithSpan / fNumZenithAngles;
            float fdAzimuthAngle = fAzimuthSpan / fNumAzimuthAngles * sin(ZenithAngle);
            float fDiffSolidAngle = fdZenithAngle * fdAzimuthAngle;
            fTotalSolidAngle += fDiffSolidAngle;
            fGatheredScattering += fCurrDirScattering * fDiffSolidAngle;
        }
    
    // Total solid angle should be 4*PI. Renormalize to fix discretization issues
    fGatheredScattering *= 4*PI / fTotalSolidAngle;

    return fGatheredScattering;
}


float ComputeScatteringOrderPS(SScreenSizeQuadVSOutput In) : SV_Target
{
    float4 f4StartPointLUTCoords = float4(ProjToUV(In.m_f2PosPS), g_GlobalCloudAttribs.f4Parameter.xy);

    float3 f3PosUSSpace, f3ViewRayUSSpace, f3LightDirUSSpace;
    float fDensityScale;
    ParticleScatteringLUTToWorldParams(f4StartPointLUTCoords, g_GlobalCloudAttribs.f4Parameter.z, f3PosUSSpace, f3ViewRayUSSpace, f3LightDirUSSpace, false, fDensityScale);

    // Intersect view ray with the unit sphere:
    float2 f2RayIsecs;
    // f3NormalizedStartPos  is located exactly on the surface; slightly move start pos inside the sphere
    // to avoid precision issues
    float3 f3BiasedPos = f3PosUSSpace + f3ViewRayUSSpace*1e-4;
    GetRaySphereIntersection(f3BiasedPos, f3ViewRayUSSpace, 0, 1.f, f2RayIsecs);
    if( f2RayIsecs.y < f2RayIsecs.x )
        return 0;

    float3 f3EndPos = f3BiasedPos + f3ViewRayUSSpace * f2RayIsecs.y;
    float fNumSteps = max(VOL_SCATTERING_IN_PARTICLE_LUT_DIM.w*2, NUM_INTEGRATION_STEPS)*2;
    float3 f3Step = (f3EndPos - f3PosUSSpace) / fNumSteps;
    float fStepLen = length(f3Step);
    float fCloudMassToCamera = 0;
    float fParticleRadius = GetParticleSize(GetCloudRingWorldStep(0, g_GlobalCloudAttribs));
    float fInscattering = 0;

    float fPrevGatheredSctr = 0;
    SAMPLE_4D_LUT(g_tex3DGatheredScattering, VOL_SCATTERING_IN_PARTICLE_LUT_DIM, f4StartPointLUTCoords, 0, fPrevGatheredSctr);
    float fStartPointDensity = 1;
    fPrevGatheredSctr *= fStartPointDensity; // Light attenuation == 1

    for(float fStepNum=1; fStepNum <= fNumSteps; ++fStepNum)
    {
        float3 f3CurrPos = f3PosUSSpace + f3Step * fStepNum;
        //float fDistToCenter = length(f3CurrPos);
        //float fMetabolDensity = GetMetabolDensity(fDistToCenter);
        float fDensity = 1;//saturate(fMetabolDensity);
        fDensity *= fDensityScale;

        fCloudMassToCamera += fDensity * fStepLen * fParticleRadius;
        float fAttenuationToCamera = exp( -g_GlobalCloudAttribs.fAttenuationCoeff * fCloudMassToCamera );

        float4 f4CurrDirLUTCoords = WorldParamsToParticleScatteringLUT(f3CurrPos, f3ViewRayUSSpace, f3LightDirUSSpace, false, 1);
        float fGatheredScattering = 0;
        SAMPLE_4D_LUT(g_tex3DGatheredScattering, VOL_SCATTERING_IN_PARTICLE_LUT_DIM, f4CurrDirLUTCoords, 0, fGatheredScattering);
        fGatheredScattering *= fAttenuationToCamera * fDensity;

        fInscattering += (fGatheredScattering + fPrevGatheredSctr) /2;
        fPrevGatheredSctr = fGatheredScattering;
    }

    return fInscattering * fStepLen * fParticleRadius * g_GlobalCloudAttribs.fScatteringCoeff;
}


float AccumulateMultipleScattering(SScreenSizeQuadVSOutput In) : SV_Target
{
    float3 f3LUTCoords = float3(ProjToUV(In.m_f2PosPS), g_GlobalCloudAttribs.f4Parameter.x);
    float fMultipleSctr = g_tex3DPrevSctrOrder.SampleLevel(samPointWrap, f3LUTCoords, 0);
    return fMultipleSctr;
}

void RenderScatteringLUTSlicePS(SScreenSizeQuadVSOutput In,
                                out float fSingleSctr : SV_Target0,
                                out float fMultipleSctr : SV_Target1)
{
    // Scattering on the surface of the sphere is stored in the last 4D-slice
    float4 f4LUTCoords = float4(ProjToUV(In.m_f2PosPS), g_GlobalCloudAttribs.f4Parameter.x, 1 - 0.5/VOL_SCATTERING_IN_PARTICLE_LUT_DIM.w);
    // We only need directions into the sphere
    f4LUTCoords.z = (f4LUTCoords.z - 0.5/SRF_SCATTERING_IN_PARTICLE_LUT_DIM.z)/2.f + 0.5/VOL_SCATTERING_IN_PARTICLE_LUT_DIM.z;

    SAMPLE_4D_LUT(g_tex3DSingleScattering,   VOL_SCATTERING_IN_PARTICLE_LUT_DIM, f4LUTCoords, 0, fSingleSctr);
    SAMPLE_4D_LUT(g_tex3DMultipleScattering, VOL_SCATTERING_IN_PARTICLE_LUT_DIM, f4LUTCoords, 0, fMultipleSctr);
}

// This shader computes exiting irradiance (not radiance!) through the specific point on the sphere
float ComputeExitancePS(SScreenSizeQuadVSOutput In) : SV_Target
{
    // Get position of the start point on the sphere
    float fU = ProjToUV(In.m_f2PosPS).x;
    float3 f3PosUSSpace, f3ViewRayUSSpaceUnused, f3LightDirUSSpace;
    float fDensityScale;
    ParticleScatteringLUTToWorldParams(float4(fU,0,0,0), g_GlobalCloudAttribs.f4Parameter.x, f3PosUSSpace, f3ViewRayUSSpaceUnused, f3LightDirUSSpace, true, fDensityScale);

    // Construct local frame for the specified point
    float3 f3LocalX, f3LocalY, f3LocalZ;
    ConstructLocalFrameXYZ(-normalize(f3PosUSSpace), f3LightDirUSSpace, f3LocalX, f3LocalY, f3LocalZ);

    float fGatheredScattering = 0;
    float fTotalSolidAngle = 0;
    const float fNumZenithAngles  = SRF_SCATTERING_IN_PARTICLE_LUT_DIM.z;
    const float fNumAzimuthAngles = SRF_SCATTERING_IN_PARTICLE_LUT_DIM.y;
    const float fZenithSpan = PI/2; // We only need directions on the hemisphere, thus PI/2
    const float fAzimuthSpan = 2 * PI;
    // Go through all directions on the upper hemisphere
    for(float ZenithAngleNum = 0.5; ZenithAngleNum < fNumZenithAngles; ++ZenithAngleNum)
        for(float AzimuthAngleNum = 0.5; AzimuthAngleNum < fNumAzimuthAngles; ++AzimuthAngleNum)
        {
            float ZenithAngle = ZenithAngleNum/fNumZenithAngles * fZenithSpan;
            float AzimuthAngle = (AzimuthAngleNum/fNumAzimuthAngles - 0.5) * fAzimuthSpan;
            
            // Get current direction. We need view direction which is -1 * exiting direction and thus points into the sphere
            float3 f3CurrDir = GetDirectionInLocalFrameXYZ(f3LocalX, f3LocalY, f3LocalZ, ZenithAngle, AzimuthAngle);
            float4 f4CurrDirLUTCoords = WorldParamsToParticleScatteringLUT(f3PosUSSpace, f3CurrDir, f3LightDirUSSpace, true, 1);
            // Load single and multiple scattering for the current direction
            float2 f2CurrDirScattering = 0;
            float fSingleScattering = g_tex3DSingleScatteringLUT.SampleLevel(samLinearClamp, f4CurrDirLUTCoords.xyz, 0);
            // Multiply single scattering with the phase function
            fSingleScattering *= HGPhaseFunc( dot(-f3CurrDir, f3LightDirUSSpace) );
            float fMultScattering = g_tex3DMultipleScatteringLUT.SampleLevel(samLinearClamp, f4CurrDirLUTCoords.xyz, 0);

            float fdZenithAngle  = fZenithSpan / fNumZenithAngles;
            float fdAzimuthAngle = fAzimuthSpan / fNumAzimuthAngles * sin(ZenithAngle);
            float fDiffSolidAngle = fdZenithAngle * fdAzimuthAngle;
            fTotalSolidAngle += fDiffSolidAngle;
            // Since we need irradiance, the single and multiple scattering irradiances must be multiplied by cos(theta)
            fGatheredScattering += (fSingleScattering + fMultScattering) * fDiffSolidAngle * cos(ZenithAngle);
        }
    // Renormalize to upper hemisphere
    fGatheredScattering *= 2*PI / fTotalSolidAngle;

    return fGatheredScattering;
}
