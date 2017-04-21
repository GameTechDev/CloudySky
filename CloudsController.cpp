/////////////////////////////////////////////////////////////////////////////////////////////
// Copyright 2017 Intel Corporation
//
// Licensed under the Apache License, Version 2.0 (the "License");// you may not use this file except in compliance with the License.// You may obtain a copy of the License at//// http://www.apache.org/licenses/LICENSE-2.0//// Unless required by applicable law or agreed to in writing, software// distributed under the License is distributed on an "AS IS" BASIS,// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.// See the License for the specific language governing permissions and// limitations under the License.
/////////////////////////////////////////////////////////////////////////////////////////////
#include "StdAfx.h"
#include "CloudsController.h"
#include "Structures.fxh"
#include "Visibility.h"
#include "ShaderMacroHelper.h"
#include "IGFXExtensionsHelper.h"

// This helper template function traverses the 3D cloud lattice
template<typename CProc>
void TraverseCloudLattice(const UINT iNumRings,
                          const UINT iInnerRingDim,
                          const UINT iRingExtension,
                          const UINT iNumLayers,
                          CProc Proc)
{
    UINT iRingDimension = iInnerRingDim + 2*iRingExtension;
    assert( (iInnerRingDim % 4) == 0 );
    for(int iRing = iNumRings-1; iRing >= 0; --iRing)
    {
        for(UINT iRow = 0; iRow < iRingDimension; ++iRow)
        {
            UINT iFirstQuart = iRingExtension + iInnerRingDim*1/4;
            UINT iThirdQuart = iRingExtension + iInnerRingDim*3/4;
            UINT iStartCol[2] = {0,           iThirdQuart   };
            UINT iEndCol[2]   = {iFirstQuart, iRingDimension};
            if( !(iRing > 0 && iRow >= iFirstQuart && iRow < iThirdQuart) )
                iStartCol[1] = iEndCol[0];

            for(int i=0; i < _countof(iStartCol); ++i)
                for(UINT iCol = iStartCol[i]; iCol < iEndCol[i]; ++iCol)
                    for(UINT iLayer = 0; iLayer < iNumLayers; ++iLayer)
                        Proc(iCol, iRow, iRing, iLayer);
        }
    }
}

CCloudsController::CCloudsController() : 
    m_strEffectPath(L"fx\\Clouds.fx"),
    m_strPreprocessingEffectPath(L"fx\\Preprocessing.fx"),
    m_uiCloudDensityTexWidth(1024), 
    m_uiCloudDensityTexHeight(1024),
    m_bPSOrderingAvailable(false),
    m_f3PrevLightDir(0,0,0)
{

}

CCloudsController::~CCloudsController()
{

}

void RenderQuad(ID3D11DeviceContext *pd3dDeviceCtx, 
                CRenderTechnique &State, 
                int iWidth = 0, int iHeight = 0,
                int iTopLeftX = 0, int iTopLeftY = 0,
                int iNumInstances = 1);

// This method handles resize event
void CCloudsController::OnResize(ID3D11Device *pDevice, 
                                 UINT uiWidth, UINT uiHeight)
{
    m_uiBackBufferWidth  = uiWidth;
    m_uiBackBufferHeight = uiHeight;
    UINT uiDownscaledWidth  = m_uiBackBufferWidth /m_CloudAttribs.uiDownscaleFactor;
    UINT uiDownscaledHeight = m_uiBackBufferHeight/m_CloudAttribs.uiDownscaleFactor;
    m_CloudAttribs.uiBackBufferWidth = m_uiBackBufferWidth;
    m_CloudAttribs.uiBackBufferHeight = m_uiBackBufferHeight;
    m_CloudAttribs.uiDownscaledBackBufferWidth  = uiDownscaledWidth;
    m_CloudAttribs.uiDownscaledBackBufferHeight = uiDownscaledHeight;

    m_CloudAttribs.fBackBufferWidth  = (float)m_uiBackBufferWidth;
    m_CloudAttribs.fBackBufferHeight = (float)m_uiBackBufferHeight;
    m_CloudAttribs.fDownscaledBackBufferWidth  = (float)uiDownscaledWidth;
    m_CloudAttribs.fDownscaledBackBufferHeight = (float)uiDownscaledHeight;

    // Release existing resources
    m_ptex2DScreenCloudColorSRV.Release();
    m_ptex2DScreenCloudColorRTV.Release();
    m_ptex2DScrSpaceCloudTransparencySRV.Release();
    m_ptex2DScrSpaceCloudTransparencyRTV.Release();
    m_ptex2DScrSpaceDistToCloudSRV.Release();
    m_ptex2DScrSpaceDistToCloudRTV.Release();

    m_ptex2DDownscaledScrCloudColorSRV.Release();
    m_ptex2DDownscaledScrCloudColorRTV.Release();
    m_ptex2DDownscaledScrCloudTransparencySRV.Release();
    m_ptex2DDownscaledScrCloudTransparencyRTV.Release();
    m_ptex2DDownscaledScrDistToCloudSRV.Release();
    m_ptex2DDownscaledScrDistToCloudRTV.Release();

    m_pbufParticleLayersSRV.Release();
    m_pbufParticleLayersUAV.Release();
    m_pbufClearParticleLayers.Release();

    // Create screen space cloud color buffer
    D3D11_TEXTURE2D_DESC ScreenCloudColorTexDesc = 
    {
        m_uiBackBufferWidth,                //UINT Width;
        m_uiBackBufferHeight,               //UINT Height;
        1,                                  //UINT MipLevels;
        1,                                  //UINT ArraySize;
        DXGI_FORMAT_R11G11B10_FLOAT,        //DXGI_FORMAT Format;
        {1,0},                              //DXGI_SAMPLE_DESC SampleDesc;
        D3D11_USAGE_DEFAULT,                //D3D11_USAGE Usage;
        D3D11_BIND_RENDER_TARGET | D3D11_BIND_SHADER_RESOURCE,           //UINT BindFlags;
        0,                                  //UINT CPUAccessFlags;
        0,                                  //UINT MiscFlags;
    };

    HRESULT hr;
    CComPtr<ID3D11Texture2D> ptex2DScreenCloudColor;
    V( pDevice->CreateTexture2D(&ScreenCloudColorTexDesc, nullptr, &ptex2DScreenCloudColor) );
    V( pDevice->CreateShaderResourceView(ptex2DScreenCloudColor, nullptr, &m_ptex2DScreenCloudColorSRV) );
    V( pDevice->CreateRenderTargetView(ptex2DScreenCloudColor, nullptr, &m_ptex2DScreenCloudColorRTV) );

    if( m_CloudAttribs.uiDownscaleFactor > 1 )
    {
        // Create downscaled screen space cloud color buffer
        D3D11_TEXTURE2D_DESC DownscaledScreenCloudColorTexDesc = ScreenCloudColorTexDesc;
        DownscaledScreenCloudColorTexDesc.Width  /= m_CloudAttribs.uiDownscaleFactor;
        DownscaledScreenCloudColorTexDesc.Height /= m_CloudAttribs.uiDownscaleFactor;
        CComPtr<ID3D11Texture2D> ptex2DDownscaledScrCloudColor;
        V( pDevice->CreateTexture2D(&DownscaledScreenCloudColorTexDesc, nullptr, &ptex2DDownscaledScrCloudColor) );
        V( pDevice->CreateShaderResourceView(ptex2DDownscaledScrCloudColor, nullptr, &m_ptex2DDownscaledScrCloudColorSRV) );
        V( pDevice->CreateRenderTargetView(ptex2DDownscaledScrCloudColor, nullptr, &m_ptex2DDownscaledScrCloudColorRTV) );
    }

    {
        // Create screen space cloud transparency buffer
        D3D11_TEXTURE2D_DESC ScreenTransparencyTexDesc = ScreenCloudColorTexDesc;
        ScreenTransparencyTexDesc.Format = DXGI_FORMAT_R8_UNORM;
        CComPtr<ID3D11Texture2D> ptex2DScreenTransparency;
        V( pDevice->CreateTexture2D(&ScreenTransparencyTexDesc, nullptr, &ptex2DScreenTransparency) );
        V( pDevice->CreateShaderResourceView(ptex2DScreenTransparency, nullptr, &m_ptex2DScrSpaceCloudTransparencySRV) );
        V( pDevice->CreateRenderTargetView(ptex2DScreenTransparency, nullptr, &m_ptex2DScrSpaceCloudTransparencyRTV) );
        if( m_CloudAttribs.uiDownscaleFactor > 1 )
        {
            // Create downscaled screen space cloud transparency buffer
            ScreenTransparencyTexDesc.Width  /= m_CloudAttribs.uiDownscaleFactor;
            ScreenTransparencyTexDesc.Height /= m_CloudAttribs.uiDownscaleFactor;
            CComPtr<ID3D11Texture2D> ptex2DDownscaledScrTransparency;
            V( pDevice->CreateTexture2D(&ScreenTransparencyTexDesc, nullptr, &ptex2DDownscaledScrTransparency) );
            V( pDevice->CreateShaderResourceView(ptex2DDownscaledScrTransparency, nullptr, &m_ptex2DDownscaledScrCloudTransparencySRV) );
            V( pDevice->CreateRenderTargetView(ptex2DDownscaledScrTransparency, nullptr, &m_ptex2DDownscaledScrCloudTransparencyRTV) );
        }
    }

    {
        // Create screen space distance to cloud buffer
        D3D11_TEXTURE2D_DESC ScreenDistToCloudTexDesc = ScreenCloudColorTexDesc;
        ScreenDistToCloudTexDesc.Format = DXGI_FORMAT_R32_FLOAT; // We need only the closest distance to camera
        CComPtr<ID3D11Texture2D> ptex2DScrSpaceDistToCloud;
        V( pDevice->CreateTexture2D(&ScreenDistToCloudTexDesc, nullptr, &ptex2DScrSpaceDistToCloud) );
        V( pDevice->CreateShaderResourceView(ptex2DScrSpaceDistToCloud, nullptr, &m_ptex2DScrSpaceDistToCloudSRV) );
        V( pDevice->CreateRenderTargetView(ptex2DScrSpaceDistToCloud, nullptr, &m_ptex2DScrSpaceDistToCloudRTV) );
        if( m_CloudAttribs.uiDownscaleFactor > 1 )
        {
            // Create downscaled screen space distance to cloud buffer
            ScreenDistToCloudTexDesc.Width  /= m_CloudAttribs.uiDownscaleFactor;
            ScreenDistToCloudTexDesc.Height /= m_CloudAttribs.uiDownscaleFactor;
            CComPtr<ID3D11Texture2D> ptex2DDownscaledScrDistToCloud;
            V( pDevice->CreateTexture2D(&ScreenDistToCloudTexDesc, nullptr, &ptex2DDownscaledScrDistToCloud) );
            V( pDevice->CreateShaderResourceView(ptex2DDownscaledScrDistToCloud, nullptr, &m_ptex2DDownscaledScrDistToCloudSRV) );
            V( pDevice->CreateRenderTargetView(ptex2DDownscaledScrDistToCloud, nullptr, &m_ptex2DDownscaledScrDistToCloudRTV) );
        }
    }

    if( m_bPSOrderingAvailable )
    {
        int iNumElements = (m_uiBackBufferWidth  / m_CloudAttribs.uiDownscaleFactor) * 
                           (m_uiBackBufferHeight / m_CloudAttribs.uiDownscaleFactor) * 
                           m_CloudAttribs.uiNumParticleLayers;
        D3D11_BUFFER_DESC ParticleLayersBufDesc = 
        {
            iNumElements * sizeof(SParticleLayer), //UINT ByteWidth;
            D3D11_USAGE_DEFAULT,                    //D3D11_USAGE Usage;
            D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS, //UINT BindFlags;
            0,                                      //UINT CPUAccessFlags;
            D3D11_RESOURCE_MISC_BUFFER_STRUCTURED,  //UINT MiscFlags;
            sizeof(SParticleLayer)               //UINT StructureByteStride;
        };
    
        CComPtr<ID3D11Buffer> pbufParticleLayers;
        V(pDevice->CreateBuffer( &ParticleLayersBufDesc, nullptr, &pbufParticleLayers));
        D3D11_SHADER_RESOURCE_VIEW_DESC SRVDesc;
        ZeroMemory(&SRVDesc, sizeof(SRVDesc));
        SRVDesc.Format = DXGI_FORMAT_UNKNOWN;
        SRVDesc.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;
        SRVDesc.Buffer.FirstElement = 0;
        SRVDesc.Buffer.NumElements = iNumElements;
        V(pDevice->CreateShaderResourceView( pbufParticleLayers, &SRVDesc, &m_pbufParticleLayersSRV));
        D3D11_UNORDERED_ACCESS_VIEW_DESC UAVDesc;
        UAVDesc.Format = DXGI_FORMAT_UNKNOWN;
        UAVDesc.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;
        UAVDesc.Buffer.FirstElement = 0;
        UAVDesc.Buffer.NumElements = iNumElements;
        UAVDesc.Buffer.Flags = 0;
        V(pDevice->CreateUnorderedAccessView( pbufParticleLayers, &UAVDesc, &m_pbufParticleLayersUAV));

        std::vector<SParticleLayer> ClearLayers(iNumElements);
        ParticleLayersBufDesc.Usage = D3D11_USAGE_IMMUTABLE;
        ParticleLayersBufDesc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
        D3D11_SUBRESOURCE_DATA InitData = {&ClearLayers[0],0,0};
        V(pDevice->CreateBuffer( &ParticleLayersBufDesc, &InitData, &m_pbufClearParticleLayers));
    }
    // Create buffer storing index of the first particle for each tile in screen space
    m_uiTileTexWidth  = (uiDownscaledWidth  + sm_iTileSize-1)/sm_iTileSize;
    m_uiTileTexHeight = (uiDownscaledHeight + sm_iTileSize-1)/sm_iTileSize;
    D3D11_TEXTURE2D_DESC FirstListIndTexDesc = 
    {
        m_uiTileTexWidth,  //UINT Width;
        m_uiTileTexHeight, //UINT Height;
        1,                                  //UINT MipLevels;
        1,                                  //UINT ArraySize;
        DXGI_FORMAT_R32_SINT,              //DXGI_FORMAT Format;
        {1,0},                              //DXGI_SAMPLE_DESC SampleDesc;
        D3D11_USAGE_DEFAULT,                //D3D11_USAGE Usage;
        D3D11_BIND_UNORDERED_ACCESS | D3D11_BIND_SHADER_RESOURCE, //UINT BindFlags;
        0,                                  //UINT CPUAccessFlags;
        0,                                  //UINT MiscFlags;
    };

    m_ptex2DScrFirstListIndSRV.Release();
    m_ptex2DScrFirstListIndUAV.Release();
    CComPtr<ID3D11Texture2D> ptex2DFirstListInd;
    V( pDevice->CreateTexture2D(&FirstListIndTexDesc, nullptr, &ptex2DFirstListInd) );
    V( pDevice->CreateShaderResourceView(ptex2DFirstListInd, nullptr, &m_ptex2DScrFirstListIndSRV) );
    V( pDevice->CreateUnorderedAccessView(ptex2DFirstListInd, nullptr, &m_ptex2DScrFirstListIndUAV) );

    m_pbufParticleListsBuffSRV.Release();
    m_pbufParticleListsBuffUAV.Release();

    // Create buffer used to store lists of particles after tiling
    int iNumElements = FirstListIndTexDesc.Width * FirstListIndTexDesc.Height * 1000;
    D3D11_BUFFER_DESC ParticleListsBufDesc = 
    {
        iNumElements * sizeof(SParticleListKnot), //UINT ByteWidth;
        D3D11_USAGE_DEFAULT,                    //D3D11_USAGE Usage;
        D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS, //UINT BindFlags;
        0,                                      //UINT CPUAccessFlags;
        D3D11_RESOURCE_MISC_BUFFER_STRUCTURED,  //UINT MiscFlags;
        sizeof(SParticleListKnot)               //UINT StructureByteStride;
    };
    
    CComPtr<ID3D11Buffer> pbufParticleListsBuffDesc;
    V(pDevice->CreateBuffer( &ParticleListsBufDesc, nullptr, &pbufParticleListsBuffDesc));
    D3D11_SHADER_RESOURCE_VIEW_DESC SRVDesc;
    ZeroMemory(&SRVDesc, sizeof(SRVDesc));
    SRVDesc.Format = DXGI_FORMAT_UNKNOWN;
    SRVDesc.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;
    SRVDesc.Buffer.FirstElement = 0;
    SRVDesc.Buffer.NumElements = iNumElements;
    V(pDevice->CreateShaderResourceView( pbufParticleListsBuffDesc, &SRVDesc, &m_pbufParticleListsBuffSRV));
    D3D11_UNORDERED_ACCESS_VIEW_DESC UAVDesc;
    UAVDesc.Format = DXGI_FORMAT_UNKNOWN;
    UAVDesc.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;
    UAVDesc.Buffer.FirstElement = 0;
    UAVDesc.Buffer.NumElements = iNumElements;
    UAVDesc.Buffer.Flags = D3D11_BUFFER_UAV_FLAG_COUNTER;
    V(pDevice->CreateUnorderedAccessView( pbufParticleListsBuffDesc, &UAVDesc, &m_pbufParticleListsBuffUAV));
}

// This method renders maximum density mip map
void CCloudsController::RenderMaxDensityMip( ID3D11Device *pDevice, 
                                             ID3D11DeviceContext *pDeviceContext, 
                                             ID3D11Texture2D *ptex2DMaxDensityMipMap, 
                                             ID3D11Texture2D *ptex2DTmpMaxDensityMipMap, 
                                             const D3D11_TEXTURE2D_DESC &MaxCloudDensityMipDesc )
{
    HRESULT hr;
    CD3DShaderMacroHelper Macros;
    DefineMacros(Macros);
    Macros.Finalize();
     
    // Create techniques
    CRenderTechnique RenderMaxDensityLevel0Tech;
    RenderMaxDensityLevel0Tech.SetDeviceAndContext(pDevice, pDeviceContext);
    RenderMaxDensityLevel0Tech.CreateVGPShadersFromFile(m_strEffectPath, "ScreenSizeQuadVS", nullptr, "RenderMaxMipLevel0PS", Macros);
    RenderMaxDensityLevel0Tech.SetDS( m_pdsDisableDepth );
    RenderMaxDensityLevel0Tech.SetRS( m_prsSolidFillNoCull );
    RenderMaxDensityLevel0Tech.SetBS( m_pbsDefault );

    CRenderTechnique RenderCoarseMaxMipLevelTech;
    RenderCoarseMaxMipLevelTech.SetDeviceAndContext(pDevice, pDeviceContext);
    RenderCoarseMaxMipLevelTech.CreateVGPShadersFromFile(m_strEffectPath, "ScreenSizeQuadVS", nullptr, "RenderCoarseMaxMipLevelPS", Macros);
    RenderCoarseMaxMipLevelTech.SetDS( m_pdsDisableDepth );
    RenderCoarseMaxMipLevelTech.SetRS( m_prsSolidFillNoCull );
    RenderCoarseMaxMipLevelTech.SetBS( m_pbsDefault );

    CComPtr<ID3D11RenderTargetView> pOrigRTV;
    CComPtr<ID3D11DepthStencilView> pOrigDSV;
    pDeviceContext->OMGetRenderTargets(1, &pOrigRTV, &pOrigDSV);

    D3D11_VIEWPORT OrigViewPort;
    UINT iNumOldViewports = 1;
    pDeviceContext->RSGetViewports(&iNumOldViewports, &OrigViewPort);

    UINT uiCurrMipWidth = MaxCloudDensityMipDesc.Width;
    UINT uiCurrMipHeight = MaxCloudDensityMipDesc.Height;
    for(UINT uiMip = 0; uiMip < MaxCloudDensityMipDesc.MipLevels; ++uiMip)
    {
        D3D11_RENDER_TARGET_VIEW_DESC RTVDesc;
        RTVDesc.Format = MaxCloudDensityMipDesc.Format;
        RTVDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE2D;
        RTVDesc.Texture2D.MipSlice = uiMip;
        CComPtr<ID3D11RenderTargetView> ptex2DTmpMaxDensityMipMapRTV;
        V(pDevice->CreateRenderTargetView(ptex2DTmpMaxDensityMipMap, &RTVDesc, &ptex2DTmpMaxDensityMipMapRTV));

        pDeviceContext->OMSetRenderTargets(1, &ptex2DTmpMaxDensityMipMapRTV.p, nullptr);

        ID3D11SamplerState *pSamplers[] = {m_psamPointWrap};
        pDeviceContext->PSSetSamplers(2, _countof(pSamplers), pSamplers);

        m_CloudAttribs.f4Parameter.x = (float)uiMip;
        UpdateConstantBuffer(pDeviceContext, m_pcbGlobalCloudAttribs, &m_CloudAttribs, sizeof(m_CloudAttribs));

        ID3D11Buffer *pCBs[] = {m_pcbGlobalCloudAttribs};
        pDeviceContext->PSSetConstantBuffers(0, _countof(pCBs), pCBs);

        if(uiMip == 0)
        {
            ID3D11ShaderResourceView *pSRVs[] = {m_ptex2DCloudDensitySRV};
            pDeviceContext->PSSetShaderResources(1, _countof(pSRVs), pSRVs);
        }
        else
        {
            ID3D11ShaderResourceView *pSRVs[] = {m_ptex2DMaxDensityMipMapSRV};
            pDeviceContext->PSSetShaderResources(3, _countof(pSRVs), pSRVs);
        }

        RenderQuad(pDeviceContext, uiMip == 0 ? RenderMaxDensityLevel0Tech : RenderCoarseMaxMipLevelTech, uiCurrMipWidth, uiCurrMipHeight);
        
        pDeviceContext->CopySubresourceRegion(ptex2DMaxDensityMipMap, uiMip, 0,0,0, ptex2DTmpMaxDensityMipMap, uiMip, nullptr);

        uiCurrMipWidth /= 2;
        uiCurrMipHeight /= 2;
    }
    assert( uiCurrMipWidth == 0 && uiCurrMipHeight == 0 );

    pDeviceContext->OMSetRenderTargets(1, &pOrigRTV.p, pOrigDSV);
    pDeviceContext->RSSetViewports(iNumOldViewports, &OrigViewPort);
}

// Method creates 2D texture storing first index of the particle list in light space
void CCloudsController::CreateLiSpFirstListIndTex(ID3D11Device *pDevice)
{
    m_ptex2DLiSpFirstListIndUAV.Release();
    m_ptex2DLiSpFirstListIndSRV.Release();

    D3D11_TEXTURE2D_DESC LiSpFirstListIndTexDesc = 
    {
        m_CloudAttribs.uiLiSpFirstListIndTexDim,  //UINT Width;
        m_CloudAttribs.uiLiSpFirstListIndTexDim, //UINT Height;
        1,                                  //UINT MipLevels;
        1,                                  //UINT ArraySize;
        DXGI_FORMAT_R32_SINT,              //DXGI_FORMAT Format;
        {1,0},                              //DXGI_SAMPLE_DESC SampleDesc;
        D3D11_USAGE_DEFAULT,                //D3D11_USAGE Usage;
        D3D11_BIND_UNORDERED_ACCESS | D3D11_BIND_SHADER_RESOURCE, //UINT BindFlags;
        0,                                  //UINT CPUAccessFlags;
        0,                                  //UINT MiscFlags;
    };

    HRESULT hr;
    CComPtr<ID3D11Texture2D> ptex2DFirstListInd;
    V( pDevice->CreateTexture2D(&LiSpFirstListIndTexDesc, nullptr, &ptex2DFirstListInd) );
    V( pDevice->CreateShaderResourceView(ptex2DFirstListInd, nullptr, &m_ptex2DLiSpFirstListIndSRV) );
    V( pDevice->CreateUnorderedAccessView(ptex2DFirstListInd, nullptr, &m_ptex2DLiSpFirstListIndUAV) );
}

// Auxiliary method which creates a buffer and views
HRESULT CCloudsController::CreateBufferAndViews(ID3D11Device *pDevice,
                                                const D3D11_BUFFER_DESC &BuffDesc, 
                                                D3D11_SUBRESOURCE_DATA *pInitData, 
                                                ID3D11Buffer**ppBuffer, 
                                                ID3D11ShaderResourceView **ppSRV, 
                                                ID3D11UnorderedAccessView **ppUAV,
                                                UINT UAVFlags /*= 0*/)
{
    HRESULT hr;

    CComPtr< ID3D11Buffer > ptmpBuffer;
    if( ppBuffer == nullptr )
        ppBuffer = &ptmpBuffer;

    V_RETURN(pDevice->CreateBuffer( &BuffDesc, pInitData, ppBuffer));

    if( ppSRV )
    {
        assert( BuffDesc.MiscFlags & D3D11_RESOURCE_MISC_BUFFER_STRUCTURED );

        D3D11_SHADER_RESOURCE_VIEW_DESC SRVDesc;
        ZeroMemory(&SRVDesc, sizeof(SRVDesc));
        SRVDesc.Format =  DXGI_FORMAT_UNKNOWN;
        SRVDesc.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;
        SRVDesc.Buffer.FirstElement = 0;
        SRVDesc.Buffer.NumElements = BuffDesc.ByteWidth / BuffDesc.StructureByteStride;
        V_RETURN(pDevice->CreateShaderResourceView( *ppBuffer, &SRVDesc, ppSRV));
    }

    if( ppUAV )
    {
        assert( BuffDesc.MiscFlags & D3D11_RESOURCE_MISC_BUFFER_STRUCTURED );

        D3D11_UNORDERED_ACCESS_VIEW_DESC UAVDesc;
        UAVDesc.Format = DXGI_FORMAT_UNKNOWN;
        UAVDesc.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;
        UAVDesc.Buffer.FirstElement = 0;
        UAVDesc.Buffer.NumElements = BuffDesc.ByteWidth / BuffDesc.StructureByteStride;
        UAVDesc.Buffer.Flags = UAVFlags;
        V_RETURN(pDevice->CreateUnorderedAccessView( *ppBuffer, &UAVDesc, ppUAV));
    }

    return S_OK;
}

// Method crates particle buffers
HRESULT CCloudsController::CreateParticleDataBuffer(ID3D11Device *pDevice)
{
    m_CloudAttribs.uiRingDimension = m_CloudAttribs.uiInnerRingDim + m_CloudAttribs.uiRingExtension*2;

    // Populate cell locations array
    m_PackedCellLocations.clear();
    m_PackedCellLocations.reserve(m_CloudAttribs.uiRingDimension * m_CloudAttribs.uiRingDimension * m_CloudAttribs.uiNumRings);
    TraverseCloudLattice(m_CloudAttribs.uiNumRings, m_CloudAttribs.uiInnerRingDim, m_CloudAttribs.uiRingExtension, 1, 
                            [&](UINT i, UINT j, UINT ring, UINT layer)
                            {
                                m_PackedCellLocations.push_back( PackParticleIJRing(i,j,ring, layer) );
                            }
                            );
    m_CloudAttribs.uiNumCells = (UINT)m_PackedCellLocations.size();

    // Populate particle locations array
    m_PackedParticleLocations.clear();
    m_PackedParticleLocations.reserve(m_CloudAttribs.uiNumCells * m_CloudAttribs.uiMaxLayers);
    TraverseCloudLattice(m_CloudAttribs.uiNumRings, m_CloudAttribs.uiInnerRingDim, m_CloudAttribs.uiRingExtension, m_CloudAttribs.uiMaxLayers, 
                            [&](UINT i, UINT j, UINT ring, UINT layer)
                            {
                                m_PackedParticleLocations.push_back( PackParticleIJRing(i,j,ring, layer) );
                            }
                            );
    m_CloudAttribs.uiMaxParticles = (UINT)m_PackedParticleLocations.size();

    // Reserve space for other aux arrays
    m_CamSpaceOrder.reserve( m_PackedParticleLocations.size() );
    m_LightSpaceOrder.reserve( m_PackedParticleLocations.size() );
    m_ParticleDistToCam.resize( m_PackedParticleLocations.size() );
    m_ParticleLightSpaceZ.resize( m_PackedParticleLocations.size() );

    // Populate camera space and light space particle lists
    m_CamSpaceOrder.clear();
    m_LightSpaceOrder.clear();
    for(UINT ind=0; ind < m_PackedParticleLocations.size(); ++ind)
    {
        UINT i,j,ring,layer;
        UnPackParticleIJRing( m_PackedParticleLocations[ind], i,j,ring,layer );
        if( layer < GetNumActiveLayers(m_CloudAttribs.uiMaxLayers, ring) )
        {
            m_CamSpaceOrder.push_back(ind);
            m_LightSpaceOrder.push_back(ind);
        }
    }
    m_RingGridAttribs.resize(m_CloudAttribs.uiNumRings);


    HRESULT hr;

    // Create cloud cell attributes buffer
    {
        D3D11_BUFFER_DESC CloudGridBuffDesc = 
        {
            m_CloudAttribs.uiNumCells * sizeof(SCloudCellAttribs), //UINT ByteWidth;
            D3D11_USAGE_DEFAULT,                    //D3D11_USAGE Usage;
            D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS, //UINT BindFlags;
            0,                                      //UINT CPUAccessFlags;
            D3D11_RESOURCE_MISC_BUFFER_STRUCTURED,  //UINT MiscFlags;
            sizeof(SCloudCellAttribs)  //UINT StructureByteStride;
        };
        m_pbufCloudGridSRV.Release();
        m_pbufCloudGridUAV.Release();
        V( CreateBufferAndViews( pDevice, CloudGridBuffDesc, nullptr, nullptr, &m_pbufCloudGridSRV, &m_pbufCloudGridUAV ) );
    }

    // Create particle attributes buffer
    {
        D3D11_BUFFER_DESC ParticleBuffDesc = 
        {
            m_CloudAttribs.uiMaxParticles * sizeof(SParticleAttribs), //UINT ByteWidth;
            D3D11_USAGE_DEFAULT,                    //D3D11_USAGE Usage;
            D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS, //UINT BindFlags;
            0,                                      //UINT CPUAccessFlags;
            D3D11_RESOURCE_MISC_BUFFER_STRUCTURED,  //UINT MiscFlags;
            sizeof(SParticleAttribs)  //UINT StructureByteStride;
        };
        m_pbufCloudParticlesSRV.Release();
        m_pbufCloudParticlesUAV.Release();
        V( CreateBufferAndViews( pDevice, ParticleBuffDesc, nullptr, nullptr, &m_pbufCloudParticlesSRV, &m_pbufCloudParticlesUAV ) );
    }

    // Create buffer for storing particle visibility flags
    {
	    D3D11_BUFFER_DESC VisibleParticlesFlagsBuffDesc = 
        {
            (m_CloudAttribs.uiMaxParticles + 31)/32 * sizeof(UINT), //UINT ByteWidth;
            D3D11_USAGE_DEFAULT,                    //D3D11_USAGE Usage;
            D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS, //UINT BindFlags;
            0,                                      //UINT CPUAccessFlags;
            D3D11_RESOURCE_MISC_BUFFER_STRUCTURED,  //UINT MiscFlags;
            sizeof(UINT)							//UINT StructureByteStride;
        };
        m_pbufVisibleParticlesFlagsSRV.Release();
        m_pbufVisibleParticlesFlagsUAV.Release();
        V( CreateBufferAndViews( pDevice, VisibleParticlesFlagsBuffDesc, nullptr, nullptr, &m_pbufVisibleParticlesFlagsSRV, &m_pbufVisibleParticlesFlagsUAV) );
    }

    // Create buffer for storing particle lighting info
    {
	    D3D11_BUFFER_DESC LightingBuffDesc = 
        {
            m_CloudAttribs.uiMaxParticles * sizeof(SCloudParticleLighting), //UINT ByteWidth;
            D3D11_USAGE_DEFAULT,                    //D3D11_USAGE Usage;
            D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS, //UINT BindFlags;
            0,                                      //UINT CPUAccessFlags;
            D3D11_RESOURCE_MISC_BUFFER_STRUCTURED,  //UINT MiscFlags;
            sizeof(SCloudParticleLighting)			//UINT StructureByteStride;
        };
        m_pbufParticlesLightingSRV.Release();
        m_pbufParticlesLightingUAV.Release();
        V( CreateBufferAndViews( pDevice, LightingBuffDesc, nullptr, nullptr, &m_pbufParticlesLightingSRV, &m_pbufParticlesLightingUAV) );
    }

    // Create buffer for storing attenuated sun light. Note that this buffer has different addressing scheme to allow direct 
    // access by the particle coordinates
    {
	    D3D11_BUFFER_DESC AttenuatedSunLightBuffDesc = 
        {
            m_CloudAttribs.uiNumRings * m_CloudAttribs.uiRingDimension * m_CloudAttribs.uiRingDimension *  m_CloudAttribs.uiMaxLayers * sizeof(D3DXVECTOR4), //UINT ByteWidth;
            D3D11_USAGE_DEFAULT,                    //D3D11_USAGE Usage;
            D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS, //UINT BindFlags;
            0,                                      //UINT CPUAccessFlags;
            D3D11_RESOURCE_MISC_BUFFER_STRUCTURED,  //UINT MiscFlags;
            sizeof(D3DXVECTOR4)       			    //UINT StructureByteStride;
        };
        m_pbufAttenuatedSunLightSRV.Release();
        m_pbufAttenuatedSunLightUAV.Release();
        V( CreateBufferAndViews( pDevice, AttenuatedSunLightBuffDesc, nullptr, nullptr, &m_pbufAttenuatedSunLightSRV, &m_pbufAttenuatedSunLightUAV) );

        AttenuatedSunLightBuffDesc.Usage = D3D11_USAGE_IMMUTABLE;
        AttenuatedSunLightBuffDesc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
        std::vector<float> MinusOne(AttenuatedSunLightBuffDesc.ByteWidth / sizeof(float), -1);
        D3D11_SUBRESOURCE_DATA InitData = {&MinusOne[0], 0,0};

        m_pbufDefaultAttenuatedSunLight.Release();
        V( CreateBufferAndViews( pDevice, AttenuatedSunLightBuffDesc, &InitData, &m_pbufDefaultAttenuatedSunLight, nullptr, nullptr) );
    }

    // Create buffer for storing cell locations
    {
        D3D11_BUFFER_DESC PackedCellLocationsBuffDesc = 
        {
            m_CloudAttribs.uiNumCells * sizeof(UINT), //UINT ByteWidth;
            D3D11_USAGE_IMMUTABLE,                  //D3D11_USAGE Usage;
            D3D11_BIND_SHADER_RESOURCE,             //UINT BindFlags;
            0,                                      //UINT CPUAccessFlags;
            D3D11_RESOURCE_MISC_BUFFER_STRUCTURED,  //UINT MiscFlags;
            sizeof(UINT)                            //UINT StructureByteStride;
        };
        
        m_pbufPackedCellLocationsSRV.Release();
        D3D11_SUBRESOURCE_DATA InitData = {&m_PackedCellLocations[0], 0, 0};
        CreateBufferAndViews( pDevice, PackedCellLocationsBuffDesc, &InitData, nullptr, &m_pbufPackedCellLocationsSRV, nullptr);
    }

    // Create buffer for storing unordered list of valid cell
    {
	    D3D11_BUFFER_DESC ValidCellsBuffDesc = 
        {
            m_CloudAttribs.uiNumCells * sizeof(UINT),//UINT ByteWidth;
            D3D11_USAGE_DEFAULT,                    //D3D11_USAGE Usage;
            D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS,            //UINT BindFlags;
            0,                                      //UINT CPUAccessFlags;
            D3D11_RESOURCE_MISC_BUFFER_STRUCTURED,  //UINT MiscFlags;
            sizeof(UINT)							//UINT StructureByteStride;
        };
        m_pbufValidCellsUnorderedList.Release();
        m_pbufValidCellsUnorderedListSRV.Release();
        m_pbufValidCellsUnorderedListUAV.Release();
        V( CreateBufferAndViews( pDevice, ValidCellsBuffDesc, nullptr, &m_pbufValidCellsUnorderedList, &m_pbufValidCellsUnorderedListSRV, &m_pbufValidCellsUnorderedListUAV, D3D11_BUFFER_UAV_FLAG_APPEND) );
    }

    // Create buffer for storing unordered list of valid particles
    {
	    D3D11_BUFFER_DESC ValidParticlesBuffDesc = 
        {
            m_CloudAttribs.uiMaxParticles * sizeof(UINT),           //UINT ByteWidth;
            D3D11_USAGE_DEFAULT,                    //D3D11_USAGE Usage;
            D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_UNORDERED_ACCESS,            //UINT BindFlags;
            0,                                      //UINT CPUAccessFlags;
            D3D11_RESOURCE_MISC_BUFFER_STRUCTURED,  //UINT MiscFlags;
            sizeof(UINT)							//UINT StructureByteStride;
        };
        m_pbufValidParticlesUnorderedList.Release();
        m_pbufValidParticlesUnorderedListSRV.Release();
        m_pbufValidParticlesUnorderedListUAV.Release();
        V( CreateBufferAndViews( pDevice, ValidParticlesBuffDesc, nullptr, &m_pbufValidParticlesUnorderedList, &m_pbufValidParticlesUnorderedListSRV, &m_pbufValidParticlesUnorderedListUAV, D3D11_BUFFER_UAV_FLAG_APPEND) );
    }

    // Create buffer for storing sorted index of all particles
    {
	    D3D11_BUFFER_DESC SortedParticlesBuffDesc = 
        {
            m_CloudAttribs.uiMaxParticles * sizeof(UINT),           //UINT ByteWidth;
            D3D11_USAGE_DYNAMIC,                    //D3D11_USAGE Usage;
            D3D11_BIND_SHADER_RESOURCE,             //UINT BindFlags;
            D3D11_CPU_ACCESS_WRITE,                 //UINT CPUAccessFlags;
            D3D11_RESOURCE_MISC_BUFFER_STRUCTURED,  //UINT MiscFlags;
            sizeof(UINT)							//UINT StructureByteStride;
        };
        m_pbufSortedParticlesOrder.Release();
        m_pbufSortedParticlesOrderSRV.Release();
        V( CreateBufferAndViews( pDevice, SortedParticlesBuffDesc, nullptr, &m_pbufSortedParticlesOrder, &m_pbufSortedParticlesOrderSRV) );
    }

    // Create buffer for storing streamed out list of visible particles
    {
	    D3D11_BUFFER_DESC SerializedParticlesBuffDesc = 
        {
            m_CloudAttribs.uiMaxParticles * sizeof(UINT),                           //UINT ByteWidth;
            D3D11_USAGE_DEFAULT,                    //D3D11_USAGE Usage;
            D3D11_BIND_STREAM_OUTPUT|D3D11_BIND_VERTEX_BUFFER|D3D11_BIND_SHADER_RESOURCE,               //UINT BindFlags;
            0,                                      //UINT CPUAccessFlags;
            0,                                      //UINT MiscFlags;
            0                                       //UINT StructureByteStride;
        };

        m_pbufSerializedVisibleParticles.Release();
        V(pDevice->CreateBuffer( &SerializedParticlesBuffDesc, nullptr, &m_pbufSerializedVisibleParticles));

        D3D11_SHADER_RESOURCE_VIEW_DESC SRVDesc;
        ZeroMemory(&SRVDesc, sizeof(SRVDesc));
        SRVDesc.Format =  DXGI_FORMAT_R32_UINT;
        SRVDesc.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;
        SRVDesc.Buffer.ElementOffset = 0;
        SRVDesc.Buffer.ElementWidth = sizeof(UINT);
        m_pbufSerializedVisibleParticlesSRV.Release();
        V(pDevice->CreateShaderResourceView( m_pbufSerializedVisibleParticles, &SRVDesc, &m_pbufSerializedVisibleParticlesSRV));
    }
    return S_OK;
}

HRESULT CCloudsController::PrecomputParticleDensity(ID3D11Device *pDevice, ID3D11DeviceContext *pDeviceContext)
{
    HRESULT hr;
    int iNumStartPosZenithAngles  = m_PrecomputedOpticalDepthTexDim.iNumStartPosZenithAngles;
    int iNumStartPosAzimuthAngles = m_PrecomputedOpticalDepthTexDim.iNumStartPosAzimuthAngles;
    int iNumDirectionZenithAngles = m_PrecomputedOpticalDepthTexDim.iNumDirectionZenithAngles;
    int iNumDirectionAzimuthAngles= m_PrecomputedOpticalDepthTexDim.iNumDirectionAzimuthAngles;

    D3D11_TEXTURE3D_DESC PrecomputedOpticalDepthTexDesc = 
    {
        iNumStartPosZenithAngles,  //UINT Width;
        iNumStartPosAzimuthAngles,  //UINT Height;
        iNumDirectionZenithAngles * iNumDirectionAzimuthAngles,  //UINT Depth;
        5, //UINT MipLevels;
        DXGI_FORMAT_R8G8_UNORM,//DXGI_FORMAT Format;
        D3D11_USAGE_DEFAULT, //D3D11_USAGE Usage;
        D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_RENDER_TARGET,//UINT BindFlags;
        0,//UINT CPUAccessFlags;
        D3D11_RESOURCE_MISC_GENERATE_MIPS //UINT MiscFlags;
    };

    CComPtr<ID3D11Texture3D> ptex3DPrecomputedParticleDensity;
    V_RETURN( pDevice->CreateTexture3D(&PrecomputedOpticalDepthTexDesc, nullptr, &ptex3DPrecomputedParticleDensity));

    m_ptex3DPrecomputedParticleDensitySRV.Release();
    V_RETURN(pDevice->CreateShaderResourceView( ptex3DPrecomputedParticleDensity, nullptr, &m_ptex3DPrecomputedParticleDensitySRV));

    if( !m_ComputeOpticalDepthTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.AddShaderMacro("DENSITY_GENERATION_METHOD", m_CloudAttribs.uiDensityGenerationMethod);
        Macros.Finalize();

        m_ComputeOpticalDepthTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_ComputeOpticalDepthTech.CreateVGPShadersFromFile(m_strPreprocessingEffectPath, "ScreenSizeQuadVS", nullptr, "PrecomputeOpticalDepthPS", Macros);
        m_ComputeOpticalDepthTech.SetDS( m_pdsDisableDepth );
        m_ComputeOpticalDepthTech.SetRS( m_prsSolidFillNoCull );
        m_ComputeOpticalDepthTech.SetBS( m_pbsDefault );
    }

    CComPtr<ID3D11RenderTargetView> pOrigRTV;
    CComPtr<ID3D11DepthStencilView> pOrigDSV;
    pDeviceContext->OMGetRenderTargets(1, &pOrigRTV, &pOrigDSV);
    
    D3D11_VIEWPORT OrigViewPort;
    UINT iNumOldViewports = 1;
    pDeviceContext->RSGetViewports(&iNumOldViewports, &OrigViewPort);

    ID3D11Buffer *pCBs[] = {m_pcbGlobalCloudAttribs/*, RenderAttribs.pcMediaScatteringParams*/};
    pDeviceContext->PSSetConstantBuffers(0, _countof(pCBs), pCBs);

    ID3D11SamplerState *pSamplers[] = {m_psamLinearClamp, m_psamLinearWrap, m_psamPointWrap};
    pDeviceContext->VSSetSamplers(0, _countof(pSamplers), pSamplers);
    pDeviceContext->PSSetSamplers(0, _countof(pSamplers), pSamplers);

    ID3D11ShaderResourceView *pSRVs[] = 
    {
        m_ptex3DNoiseSRV,
    };
    
    for(UINT Slice = 0; Slice < PrecomputedOpticalDepthTexDesc.Depth; ++Slice)
    {
        D3D11_RENDER_TARGET_VIEW_DESC RTVDesc;
        RTVDesc.Format = PrecomputedOpticalDepthTexDesc.Format;
        RTVDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE3D;
        RTVDesc.Texture3D.MipSlice = 0;
        RTVDesc.Texture3D.FirstWSlice = Slice;
        RTVDesc.Texture3D.WSize = 1;

        CComPtr<ID3D11RenderTargetView> pSliceRTV;
        V_RETURN(pDevice->CreateRenderTargetView( ptex3DPrecomputedParticleDensity, &RTVDesc, &pSliceRTV));
        
        UINT uiDirectionZenith = Slice % iNumDirectionZenithAngles;
        UINT uiDirectionAzimuth= Slice / iNumDirectionZenithAngles;
        m_CloudAttribs.f4Parameter.x = ((float)uiDirectionZenith + 0.5f)  / (float)iNumDirectionZenithAngles;
        m_CloudAttribs.f4Parameter.y = ((float)uiDirectionAzimuth + 0.5f) / (float)iNumDirectionAzimuthAngles;
        assert(0 < m_CloudAttribs.f4Parameter.x && m_CloudAttribs.f4Parameter.x < 1);
        assert(0 < m_CloudAttribs.f4Parameter.y && m_CloudAttribs.f4Parameter.y < 1);
        UpdateConstantBuffer(pDeviceContext, m_pcbGlobalCloudAttribs, &m_CloudAttribs, sizeof(m_CloudAttribs));

        pDeviceContext->OMSetRenderTargets(1, &pSliceRTV.p, nullptr);
        
        pDeviceContext->PSSetShaderResources(0, _countof(pSRVs), pSRVs);

        RenderQuad(pDeviceContext, m_ComputeOpticalDepthTech, PrecomputedOpticalDepthTexDesc.Width, PrecomputedOpticalDepthTexDesc.Height);
    }
    pDeviceContext->GenerateMips( m_ptex3DPrecomputedParticleDensitySRV);

    pDeviceContext->OMSetRenderTargets(1, &pOrigRTV.p, pOrigDSV);
    pDeviceContext->RSSetViewports(iNumOldViewports, &OrigViewPort);

    return S_OK;
}

HRESULT CCloudsController::PrecomputeScatteringInParticle(ID3D11Device *pDevice, ID3D11DeviceContext *pDeviceContext)
{
    HRESULT hr;
    
    LPCTSTR SingleSctrTexPath   = L"media\\SingleSctr.dds";
    LPCTSTR MultipleSctrTexPath = L"media\\MultipleSctr.dds";
    HRESULT hr1 = D3DX11CreateShaderResourceViewFromFile(pDevice, SingleSctrTexPath, nullptr, nullptr, &m_ptex3DSingleSctrInParticleLUT_SRV, nullptr);
    HRESULT hr2 = D3DX11CreateShaderResourceViewFromFile(pDevice, MultipleSctrTexPath, nullptr, nullptr, &m_ptex3DMultipleSctrInParticleLUT_SRV, nullptr);
    if( SUCCEEDED(hr1) && SUCCEEDED(hr2) )
        return S_OK;

    D3D11_TEXTURE3D_DESC PrecomputedScatteringTexDesc = 
    {
        m_PrecomputedSctrInParticleLUTDim.iNumStartPosZenithAngles,  //UINT Width;
        m_PrecomputedSctrInParticleLUTDim.iNumViewDirAzimuthAngles,  //UINT Height;
        // We are only interested in rays going into the sphere, which is half of the total number of diections
        m_PrecomputedSctrInParticleLUTDim.iNumViewDirZenithAngles/2 * m_PrecomputedSctrInParticleLUTDim.iNumDensityLevels,  //UINT Depth;
        1, //UINT MipLevels;
        DXGI_FORMAT_R16_FLOAT,//DXGI_FORMAT Format;
        D3D11_USAGE_DEFAULT, //D3D11_USAGE Usage;
        D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_RENDER_TARGET,//UINT BindFlags;
        0,//UINT CPUAccessFlags;
        0 //UINT MiscFlags;
    };

    CComPtr<ID3D11Texture3D> ptex3DSingleScatteringInParticleLUT, ptex3DMultipleScatteringInParticleLUT;
    V_RETURN( pDevice->CreateTexture3D(&PrecomputedScatteringTexDesc, nullptr, &ptex3DSingleScatteringInParticleLUT));
    V_RETURN( pDevice->CreateTexture3D(&PrecomputedScatteringTexDesc, nullptr, &ptex3DMultipleScatteringInParticleLUT));
    m_ptex3DSingleSctrInParticleLUT_SRV.Release();
    m_ptex3DMultipleSctrInParticleLUT_SRV.Release();
    V_RETURN(pDevice->CreateShaderResourceView( ptex3DSingleScatteringInParticleLUT,   nullptr, &m_ptex3DSingleSctrInParticleLUT_SRV));
    V_RETURN(pDevice->CreateShaderResourceView( ptex3DMultipleScatteringInParticleLUT, nullptr, &m_ptex3DMultipleSctrInParticleLUT_SRV));
    
    D3D11_TEXTURE3D_DESC TmpScatteringTexDesc = PrecomputedScatteringTexDesc;
    TmpScatteringTexDesc.Format = DXGI_FORMAT_R32_FLOAT;
    TmpScatteringTexDesc.Depth = m_PrecomputedSctrInParticleLUTDim.iNumViewDirZenithAngles * m_PrecomputedSctrInParticleLUTDim.iNumDistancesFromCenter;

    CComPtr<ID3D11Texture3D> ptex3DSingleSctr, ptex3DGatheredScatteringN, ptex3DSctrOrderN, ptex3DMultipeScattering;
    V_RETURN( pDevice->CreateTexture3D(&TmpScatteringTexDesc, nullptr, &ptex3DSingleSctr));
    V_RETURN( pDevice->CreateTexture3D(&TmpScatteringTexDesc, nullptr, &ptex3DGatheredScatteringN));
    V_RETURN( pDevice->CreateTexture3D(&TmpScatteringTexDesc, nullptr, &ptex3DSctrOrderN));
    V_RETURN( pDevice->CreateTexture3D(&TmpScatteringTexDesc, nullptr, &ptex3DMultipeScattering));

    std::vector< CComPtr<ID3D11RenderTargetView> > ptex3DSingleSctrRTVs(TmpScatteringTexDesc.Depth);
    std::vector< CComPtr<ID3D11RenderTargetView> > ptex3DGatheredScatteringN_RTVs(TmpScatteringTexDesc.Depth);
    std::vector< CComPtr<ID3D11RenderTargetView> > ptex3DSctrOrderN_RTVs(TmpScatteringTexDesc.Depth);
    std::vector< CComPtr<ID3D11RenderTargetView> > ptex3DMultipeScatteringRTVs(TmpScatteringTexDesc.Depth);
    
    CComPtr<ID3D11ShaderResourceView> ptex3DSingleSctrSRV;
    CComPtr<ID3D11ShaderResourceView> ptex3DGatheredScatteringN_SRV;
    CComPtr<ID3D11ShaderResourceView> ptex3DSctrOrderN_SRV;
    CComPtr<ID3D11ShaderResourceView> ptex3DMultipeScatteringSRV;
    V_RETURN(pDevice->CreateShaderResourceView( ptex3DSingleSctr,          nullptr, &ptex3DSingleSctrSRV));
    V_RETURN(pDevice->CreateShaderResourceView( ptex3DGatheredScatteringN, nullptr, &ptex3DGatheredScatteringN_SRV));
    V_RETURN(pDevice->CreateShaderResourceView( ptex3DSctrOrderN,          nullptr, &ptex3DSctrOrderN_SRV));
    V_RETURN(pDevice->CreateShaderResourceView( ptex3DMultipeScattering,   nullptr, &ptex3DMultipeScatteringSRV));

    for(UINT Slice = 0; Slice < TmpScatteringTexDesc.Depth; ++Slice)
    {
        D3D11_RENDER_TARGET_VIEW_DESC RTVDesc;
        RTVDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE3D;
        RTVDesc.Texture3D.MipSlice = 0;
        RTVDesc.Texture3D.FirstWSlice = Slice;
        RTVDesc.Texture3D.WSize = 1;
        RTVDesc.Format = TmpScatteringTexDesc.Format;
        V_RETURN(pDevice->CreateRenderTargetView( ptex3DSingleSctr,          &RTVDesc, &ptex3DSingleSctrRTVs[Slice])          );
        V_RETURN(pDevice->CreateRenderTargetView( ptex3DGatheredScatteringN, &RTVDesc, &ptex3DGatheredScatteringN_RTVs[Slice]));
        V_RETURN(pDevice->CreateRenderTargetView( ptex3DSctrOrderN,          &RTVDesc, &ptex3DSctrOrderN_RTVs[Slice])         );
        V_RETURN(pDevice->CreateRenderTargetView( ptex3DMultipeScattering,   &RTVDesc, &ptex3DMultipeScatteringRTVs[Slice])   );
    }

    
    if( !m_ComputeSingleSctrInParticleTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.Finalize();

        m_ComputeSingleSctrInParticleTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_ComputeSingleSctrInParticleTech.CreateVGPShadersFromFile(m_strPreprocessingEffectPath, "ScreenSizeQuadVS", nullptr, "PrecomputeSingleSctrPS", Macros);
        m_ComputeSingleSctrInParticleTech.SetDS( m_pdsDisableDepth );
        m_ComputeSingleSctrInParticleTech.SetRS( m_prsSolidFillNoCull );
        m_ComputeSingleSctrInParticleTech.SetBS( m_pbsDefault );
    }

    if( !m_RenderScatteringLUTSliceTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.Finalize();

        m_RenderScatteringLUTSliceTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_RenderScatteringLUTSliceTech.CreateVGPShadersFromFile(m_strPreprocessingEffectPath, "ScreenSizeQuadVS", nullptr, "RenderScatteringLUTSlicePS", Macros);
        m_RenderScatteringLUTSliceTech.SetDS( m_pdsDisableDepth );
        m_RenderScatteringLUTSliceTech.SetRS( m_prsSolidFillNoCull );
        m_RenderScatteringLUTSliceTech.SetBS( m_pbsDefault );
    }

    if( !m_GatherPrevSctrOrderTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.Finalize();

        m_GatherPrevSctrOrderTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_GatherPrevSctrOrderTech.CreateVGPShadersFromFile(m_strPreprocessingEffectPath, "ScreenSizeQuadVS", nullptr, "GatherScatteringPS", Macros);
        m_GatherPrevSctrOrderTech.SetDS( m_pdsDisableDepth );
        m_GatherPrevSctrOrderTech.SetRS( m_prsSolidFillNoCull );
        m_GatherPrevSctrOrderTech.SetBS( m_pbsDefault );
    }
     
    if( !m_ComputeScatteringOrderTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.Finalize();

        m_ComputeScatteringOrderTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_ComputeScatteringOrderTech.CreateVGPShadersFromFile(m_strPreprocessingEffectPath, "ScreenSizeQuadVS", nullptr, "ComputeScatteringOrderPS", Macros);
        m_ComputeScatteringOrderTech.SetDS( m_pdsDisableDepth );
        m_ComputeScatteringOrderTech.SetRS( m_prsSolidFillNoCull );
        m_ComputeScatteringOrderTech.SetBS( m_pbsDefault );
    }

    if( !m_AccumulateInscatteringTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.Finalize();

        m_AccumulateInscatteringTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_AccumulateInscatteringTech.CreateVGPShadersFromFile(m_strPreprocessingEffectPath, "ScreenSizeQuadVS", nullptr, "AccumulateMultipleScattering", Macros);
        m_AccumulateInscatteringTech.SetDS( m_pdsDisableDepth );
        m_AccumulateInscatteringTech.SetRS( m_prsSolidFillNoCull );

        D3D11_BLEND_DESC AdditiveBlendStateDesc;
        ZeroMemory(&AdditiveBlendStateDesc, sizeof(AdditiveBlendStateDesc));
        AdditiveBlendStateDesc.IndependentBlendEnable = FALSE;
        for(int i=0; i< _countof(AdditiveBlendStateDesc.RenderTarget); i++)
            AdditiveBlendStateDesc.RenderTarget[i].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;
        AdditiveBlendStateDesc.RenderTarget[0].BlendEnable = TRUE;
        AdditiveBlendStateDesc.RenderTarget[0].BlendOp     = D3D11_BLEND_OP_ADD;
        AdditiveBlendStateDesc.RenderTarget[0].BlendOpAlpha= D3D11_BLEND_OP_ADD;
        AdditiveBlendStateDesc.RenderTarget[0].DestBlend   = D3D11_BLEND_ONE;
        AdditiveBlendStateDesc.RenderTarget[0].DestBlendAlpha= D3D11_BLEND_ONE;
        AdditiveBlendStateDesc.RenderTarget[0].SrcBlend     = D3D11_BLEND_ONE;
        AdditiveBlendStateDesc.RenderTarget[0].SrcBlendAlpha= D3D11_BLEND_ONE;
        CComPtr<ID3D11BlendState> pAdditiveBlendBS;
        V_RETURN( pDevice->CreateBlendState( &AdditiveBlendStateDesc, &pAdditiveBlendBS) );
        m_AccumulateInscatteringTech.SetBS( pAdditiveBlendBS );
    }

    CComPtr<ID3D11RenderTargetView> pOrigRTV;
    CComPtr<ID3D11DepthStencilView> pOrigDSV;
    pDeviceContext->OMGetRenderTargets(1, &pOrigRTV, &pOrigDSV);
    
    D3D11_VIEWPORT OrigViewPort;
    UINT iNumOldViewports = 1;
    pDeviceContext->RSGetViewports(&iNumOldViewports, &OrigViewPort);

    ID3D11Buffer *pCBs[] = {m_pcbGlobalCloudAttribs/*, RenderAttribs.pcMediaScatteringParams*/};
    pDeviceContext->PSSetConstantBuffers(0, _countof(pCBs), pCBs);

    ID3D11SamplerState *pSamplers[] = {m_psamLinearClamp, m_psamLinearWrap, m_psamPointWrap};
    pDeviceContext->VSSetSamplers(0, _countof(pSamplers), pSamplers);
    pDeviceContext->PSSetSamplers(0, _countof(pSamplers), pSamplers);

    for(int iDensityLevel = 0; iDensityLevel < m_PrecomputedSctrInParticleLUTDim.iNumDensityLevels; ++iDensityLevel)
    {
        for(UINT Slice = 0; Slice < TmpScatteringTexDesc.Depth; ++Slice)
        {
            float Zero[4]={0,0,0,0};
            pDeviceContext->ClearRenderTargetView(ptex3DMultipeScatteringRTVs[Slice], Zero);
        }

        // Precompute single scattering
        for(UINT Slice = 0; Slice < TmpScatteringTexDesc.Depth; ++Slice)
        {
            UINT uiViewDirZenith = Slice % m_PrecomputedSctrInParticleLUTDim.iNumViewDirZenithAngles;
            UINT uiDistFromCenter = Slice / m_PrecomputedSctrInParticleLUTDim.iNumViewDirZenithAngles;
            m_CloudAttribs.f4Parameter.x = ((float)uiViewDirZenith + 0.5f)  / (float)m_PrecomputedSctrInParticleLUTDim.iNumViewDirZenithAngles;
            m_CloudAttribs.f4Parameter.y = ((float)uiDistFromCenter + 0.5f) / (float)m_PrecomputedSctrInParticleLUTDim.iNumDistancesFromCenter;
            m_CloudAttribs.f4Parameter.z = ((float)iDensityLevel + 0.5f) / (float)m_PrecomputedSctrInParticleLUTDim.iNumDensityLevels;
            assert(0 < m_CloudAttribs.f4Parameter.x && m_CloudAttribs.f4Parameter.x < 1);
            assert(0 < m_CloudAttribs.f4Parameter.y && m_CloudAttribs.f4Parameter.y < 1);
            assert(0 < m_CloudAttribs.f4Parameter.z && m_CloudAttribs.f4Parameter.z < 1);
            UpdateConstantBuffer(pDeviceContext, m_pcbGlobalCloudAttribs, &m_CloudAttribs, sizeof(m_CloudAttribs));

            ID3D11RenderTargetView *pSliceRTV = ptex3DSingleSctrRTVs[Slice];
            pDeviceContext->OMSetRenderTargets(1, &pSliceRTV, nullptr);
        
            ID3D11ShaderResourceView *pSRVs[] = 
            {
                m_ptex3DNoiseSRV,
            };

            pDeviceContext->PSSetShaderResources(0, _countof(pSRVs), pSRVs);

            RenderQuad(pDeviceContext, m_ComputeSingleSctrInParticleTech, TmpScatteringTexDesc.Width, TmpScatteringTexDesc.Height);
        }
    
        // Number of scattering orders is chosen so as to obtain reasonable exitance through the particle surface
        const int iMaxScatteringOrder = 18;
        for(int iSctrOrder = 1; iSctrOrder < iMaxScatteringOrder; ++iSctrOrder)
        {
            for(int iPass = 0; iPass < 3; ++iPass)
            {
                // Gather scattering of previous order
                for(UINT Slice = 0; Slice < TmpScatteringTexDesc.Depth; ++Slice)
                {
                    if( iPass < 2 )
                    {
                        UINT uiViewDirZenith = Slice % m_PrecomputedSctrInParticleLUTDim.iNumViewDirZenithAngles;
                        UINT uiDistFromCenter = Slice / m_PrecomputedSctrInParticleLUTDim.iNumViewDirZenithAngles;
                        m_CloudAttribs.f4Parameter.x = ((float)uiViewDirZenith + 0.5f)  / (float)m_PrecomputedSctrInParticleLUTDim.iNumViewDirZenithAngles;
                        m_CloudAttribs.f4Parameter.y = ((float)uiDistFromCenter + 0.5f) / (float)m_PrecomputedSctrInParticleLUTDim.iNumDistancesFromCenter;
                        m_CloudAttribs.f4Parameter.z = ((float)iDensityLevel + 0.5f)  / (float)m_PrecomputedSctrInParticleLUTDim.iNumDensityLevels;
                        assert(0 < m_CloudAttribs.f4Parameter.x && m_CloudAttribs.f4Parameter.x < 1);
                        assert(0 < m_CloudAttribs.f4Parameter.y && m_CloudAttribs.f4Parameter.y < 1);
                        assert(0 < m_CloudAttribs.f4Parameter.z && m_CloudAttribs.f4Parameter.z < 1);
                        m_CloudAttribs.f4Parameter.w = (float)iSctrOrder;
                    }
                    else
                    {
                        m_CloudAttribs.f4Parameter.x = ((float)Slice + 0.5f) / (float)TmpScatteringTexDesc.Depth;
                        assert(0 < m_CloudAttribs.f4Parameter.x && m_CloudAttribs.f4Parameter.x < 1);
                    }
                    UpdateConstantBuffer(pDeviceContext, m_pcbGlobalCloudAttribs, &m_CloudAttribs, sizeof(m_CloudAttribs));

                    ID3D11RenderTargetView *pSliceRTV = nullptr;
                    CRenderTechnique *pTechnique = nullptr;
                    ID3D11ShaderResourceView *pSRVs[1] = {nullptr};
                    switch(iPass)
                    {
                        // Gather scattering of previous order
                        case 0: 
                            pSRVs[0] = iSctrOrder > 1 ? ptex3DSctrOrderN_SRV : ptex3DSingleSctrSRV;
                            pSliceRTV = ptex3DGatheredScatteringN_RTVs[Slice];
                            pTechnique = &m_GatherPrevSctrOrderTech;
                        break;

                        // Compute current scattering order
                        case 1: 
                            pSRVs[0] = ptex3DGatheredScatteringN_SRV;
                            pSliceRTV = ptex3DSctrOrderN_RTVs[Slice];
                            pTechnique = &m_ComputeScatteringOrderTech;
                        break;

                        // Accumulate current scattering order
                        case 2: 
                            pSRVs[0] = ptex3DSctrOrderN_SRV;
                            pSliceRTV = ptex3DMultipeScatteringRTVs[Slice];
                            pTechnique = &m_AccumulateInscatteringTech;
                        break;
                    }

                    pDeviceContext->OMSetRenderTargets(1, &pSliceRTV, nullptr);
                    pDeviceContext->PSSetShaderResources(0, _countof(pSRVs), pSRVs);

                    RenderQuad(pDeviceContext, *pTechnique, TmpScatteringTexDesc.Width, TmpScatteringTexDesc.Height);
                }
            }
        }

        // Copy single and multiple scattering to the textures
        UINT uiNumSlicesPerDensityLevel = PrecomputedScatteringTexDesc.Depth / m_PrecomputedSctrInParticleLUTDim.iNumDensityLevels;
        for(UINT Slice = 0; Slice < uiNumSlicesPerDensityLevel; ++Slice)
        {
            D3D11_RENDER_TARGET_VIEW_DESC RTVDesc;
            RTVDesc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE3D;
            RTVDesc.Texture3D.MipSlice = 0;
            RTVDesc.Texture3D.FirstWSlice = Slice + iDensityLevel * uiNumSlicesPerDensityLevel;
            RTVDesc.Texture3D.WSize = 1;
            RTVDesc.Format = PrecomputedScatteringTexDesc.Format;
            CComPtr<ID3D11RenderTargetView> pSingleSctrSliceRTV, pMultSctrSliceRTV;
            V_RETURN(pDevice->CreateRenderTargetView( ptex3DSingleScatteringInParticleLUT, &RTVDesc, &pSingleSctrSliceRTV));
            V_RETURN(pDevice->CreateRenderTargetView( ptex3DMultipleScatteringInParticleLUT, &RTVDesc, &pMultSctrSliceRTV));

            m_CloudAttribs.f4Parameter.x = ((float)Slice + 0.5f)  / (float)uiNumSlicesPerDensityLevel;
            UpdateConstantBuffer(pDeviceContext, m_pcbGlobalCloudAttribs, &m_CloudAttribs, sizeof(m_CloudAttribs));

            ID3D11RenderTargetView *pRTVs[] = {pSingleSctrSliceRTV, pMultSctrSliceRTV};
            pDeviceContext->OMSetRenderTargets(_countof(pRTVs), pRTVs, nullptr);

            ID3D11ShaderResourceView *pSRVs[] = 
            {
                ptex3DSingleSctrSRV,
                ptex3DMultipeScatteringSRV
            };
            pDeviceContext->PSSetShaderResources(0, _countof(pSRVs), pSRVs);

            RenderQuad(pDeviceContext, m_RenderScatteringLUTSliceTech, PrecomputedScatteringTexDesc.Width, PrecomputedScatteringTexDesc.Height);
        }
    }

    D3DX11SaveTextureToFile(pDeviceContext, ptex3DSingleScatteringInParticleLUT, D3DX11_IFF_DDS, SingleSctrTexPath);
    D3DX11SaveTextureToFile(pDeviceContext, ptex3DMultipleScatteringInParticleLUT, D3DX11_IFF_DDS, MultipleSctrTexPath);

    pDeviceContext->OMSetRenderTargets(1, &pOrigRTV.p, pOrigDSV);
    pDeviceContext->RSSetViewports(iNumOldViewports, &OrigViewPort);

    return S_OK;
}

HRESULT CCloudsController::ComputeExitance(ID3D11Device *pDevice, ID3D11DeviceContext *pDeviceContext)
{
    HRESULT hr;

    if( !m_ComputeExitanceTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.Finalize();

        m_ComputeExitanceTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_ComputeExitanceTech.CreateVGPShadersFromFile(m_strPreprocessingEffectPath, "ScreenSizeQuadVS", nullptr, "ComputeExitancePS", Macros);
        m_ComputeExitanceTech.SetDS( m_pdsDisableDepth );
        m_ComputeExitanceTech.SetRS( m_prsSolidFillNoCull );
        m_ComputeExitanceTech.SetBS( m_pbsDefault );
    }

    D3D11_TEXTURE2D_DESC ExitanceTexDesc = 
    {
        m_PrecomputedSctrInParticleLUTDim.iNumStartPosZenithAngles,                //UINT Width;
        1,               //UINT Height;
        1,                                  //UINT MipLevels;
        1,                                  //UINT ArraySize;
        DXGI_FORMAT_R32_FLOAT,              //DXGI_FORMAT Format;
        {1,0},                              //DXGI_SAMPLE_DESC SampleDesc;
        D3D11_USAGE_DEFAULT,                //D3D11_USAGE Usage;
        D3D11_BIND_RENDER_TARGET,           //UINT BindFlags;
        0,                                  //UINT CPUAccessFlags;
        0,                                  //UINT MiscFlags;
    };

    CComPtr<ID3D11Texture2D> ptex2DExitance, ptex2DExitanceStaging;
    CComPtr<ID3D11RenderTargetView> ptex2DExitanceRTV;
    V( pDevice->CreateTexture2D(&ExitanceTexDesc, nullptr, &ptex2DExitance));
    V( pDevice->CreateRenderTargetView(ptex2DExitance, nullptr, &ptex2DExitanceRTV));

    ExitanceTexDesc.BindFlags = 0;
    ExitanceTexDesc.Usage = D3D11_USAGE_STAGING;
    ExitanceTexDesc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
    V( pDevice->CreateTexture2D(&ExitanceTexDesc, nullptr, &ptex2DExitanceStaging));


    CComPtr<ID3D11RenderTargetView> pOrigRTV;
    CComPtr<ID3D11DepthStencilView> pOrigDSV;
    pDeviceContext->OMGetRenderTargets(1, &pOrigRTV, &pOrigDSV);
    
    D3D11_VIEWPORT OrigViewPort;
    UINT iNumOldViewports = 1;
    pDeviceContext->RSGetViewports(&iNumOldViewports, &OrigViewPort);

    ID3D11Buffer *pCBs[] = {m_pcbGlobalCloudAttribs/*, RenderAttribs.pcMediaScatteringParams*/};
    pDeviceContext->PSSetConstantBuffers(0, _countof(pCBs), pCBs);

    ID3D11SamplerState *pSamplers[] = {m_psamLinearClamp, m_psamLinearWrap, m_psamPointWrap};
    pDeviceContext->VSSetSamplers(0, _countof(pSamplers), pSamplers);
    pDeviceContext->PSSetSamplers(0, _countof(pSamplers), pSamplers);

    pDeviceContext->OMSetRenderTargets(1, &ptex2DExitanceRTV.p, nullptr);
    ID3D11ShaderResourceView *pSRVs[] = 
    {
        m_ptex3DSingleSctrInParticleLUT_SRV,
        m_ptex3DMultipleSctrInParticleLUT_SRV
    };
    pDeviceContext->PSSetShaderResources(0, _countof(pSRVs), pSRVs);

    assert(false);
    float fDensityLevel = 0;
    m_CloudAttribs.f4Parameter.z = ((float)fDensityLevel + 0.5f) / (float)m_PrecomputedSctrInParticleLUTDim.iNumDensityLevels;
    UpdateConstantBuffer(pDeviceContext, m_pcbGlobalCloudAttribs, &m_CloudAttribs, sizeof(m_CloudAttribs));

    RenderQuad(pDeviceContext, m_ComputeExitanceTech, ExitanceTexDesc.Width, ExitanceTexDesc.Height);

    pDeviceContext->CopyResource(ptex2DExitanceStaging, ptex2DExitance);
    D3D11_MAPPED_SUBRESOURCE MapData;
    pDeviceContext->Map(ptex2DExitanceStaging, 0, D3D11_MAP_READ, 0, &MapData);
    float *pData = (float*)MapData.pData;
    float fTotalArea = 0;
    float fTotalExitance = 0;
    // Go through all start points on the sphere
    for(UINT iStartPosZenith = 0; iStartPosZenith < ExitanceTexDesc.Width; ++iStartPosZenith)
    {
        float fStartPosZenithAngle = ((float)iStartPosZenith + 0.5f) / (float)ExitanceTexDesc.Width * PI;
        // Exiting irradiance is the same for all points with the fixed zenith angle
        float fdAzimuth = sin(fStartPosZenithAngle) * 2*PI;
        float fdZenith = PI / (float)ExitanceTexDesc.Width;
        // Get area for the current exiting irradiance. Note that since we compute exiting flux, we need
        // area, not solid angle. However, since we assume sphere radius is 1, area is the same as solid angle:
        float fDiffArea = fdAzimuth * fdZenith;
        fTotalArea += fDiffArea;
        // Accumulate exiting flux
        fTotalExitance += pData[iStartPosZenith] * fDiffArea;
    }
    // Renormalize to the area of the unit sphere 4*PI
    fTotalExitance *= 4*PI / fTotalArea;
    pDeviceContext->Unmap(ptex2DExitanceStaging, 0);

    // Total flux of the sun light through the sphere is PI * r^2 * E_Sun
    // Thus for the unit raidus and sun irradiance, the flux will be PI
    // Since almost all sun light is scattered, fTotalExitance should be very close to PI

    pDeviceContext->OMSetRenderTargets(1, &pOrigRTV.p, pOrigDSV);
    pDeviceContext->RSSetViewports(iNumOldViewports, &OrigViewPort);
    
    return S_OK;
}

float CubicInterpolate(float ym1, float y0, float y1, float y2, float x)
{
    float b0 = 0*ym1 + 6*y0 + 0*y1 + 0*y2;
    float b1 =-2*ym1 - 3*y0 + 6*y1 - 1*y2;
    float b2 = 3*ym1 - 6*y0 + 3*y1 + 0*y2;
    float b3 =-1*ym1 + 3*y0 - 3*y1 + 1*y2;
    float x2 = x*x;
    float x3 = x2*x;
    return 1.f/6.f * (b0 + x*b1 + x2*b2 + x3*b3);
}

HRESULT CCloudsController::Create3DNoise(ID3D11Device *pDevice)
{
    HRESULT hr;
    // Create 3D noise
    UINT uiMips = 8;
    UINT uiDim = 1 << (uiMips-1);
    D3D11_TEXTURE3D_DESC NoiseTexDesc = 
    {
        uiDim,  //UINT Width;
        uiDim,  //UINT Height;
        uiDim,  //UINT Depth;
        uiMips, //UINT MipLevels;
        DXGI_FORMAT_R8_UNORM,//DXGI_FORMAT Format;
        D3D11_USAGE_DEFAULT, //D3D11_USAGE Usage;
        D3D11_BIND_SHADER_RESOURCE,//UINT BindFlags;
        0,//UINT CPUAccessFlags;
        0//UINT MiscFlags;
    };
    size_t DataSize = 0;
    for(UINT Mip=0; Mip < uiMips; ++Mip)
        DataSize += (NoiseTexDesc.Width>>Mip) * (NoiseTexDesc.Height>>Mip) * (NoiseTexDesc.Depth>>Mip);
    std::vector<float> NoiseData(DataSize);

    #define NOISE(i,j,k) NoiseData[i + j * NoiseTexDesc.Width + k * (NoiseTexDesc.Width * NoiseTexDesc.Height)]

    // Populate texture with random noise
    UINT InitialStep = 8;
    for(UINT i=0; i < NoiseTexDesc.Width; i+=InitialStep)
        for(UINT j=0; j < NoiseTexDesc.Height; j+=InitialStep)
            for(UINT k=0; k < NoiseTexDesc.Depth; k+=InitialStep)
                NOISE(i,j,k) = (float)rand() / (float)RAND_MAX;

    // Smooth rows
    for(UINT i=0; i < NoiseTexDesc.Width; ++i)
        for(UINT j=0; j < NoiseTexDesc.Height; j+=InitialStep)
            for(UINT k=0; k < NoiseTexDesc.Depth; k+=InitialStep)
            {
                int i0 = (i/InitialStep)*InitialStep;
                int im1 = i0-InitialStep;
                if( im1 < 0 )im1 += NoiseTexDesc.Width;
                int i1 = (i0+InitialStep) % NoiseTexDesc.Width;
                int i2 = (i0+2*InitialStep) % NoiseTexDesc.Width;
                NOISE(i,j,k) = CubicInterpolate( NOISE(im1,j,k), NOISE(i0,j,k), NOISE(i1,j,k), NOISE(i2,j,k), (float)(i-i0) / (float)InitialStep );
            }

    // Smooth columns
    for(UINT i=0; i < NoiseTexDesc.Width; ++i)
        for(UINT j=0; j < NoiseTexDesc.Height; ++j)
            for(UINT k=0; k < NoiseTexDesc.Depth; k+=InitialStep)
            {
                int j0 = (j/InitialStep)*InitialStep;
                int jm1 = j0 - InitialStep;
                if( jm1 < 0 )jm1+=NoiseTexDesc.Height;
                int j1 = (j0+InitialStep) % NoiseTexDesc.Height;
                int j2 = (j0+2*InitialStep) % NoiseTexDesc.Height;
                NOISE(i,j,k) = CubicInterpolate(NOISE(i,jm1,k), NOISE(i,j0,k), NOISE(i,j1,k), NOISE(i,j2,k), (float)(j-j0) / (float)InitialStep);
            }

    // Smooth in depth direction
    for(UINT i=0; i < NoiseTexDesc.Width; ++i)
        for(UINT j=0; j < NoiseTexDesc.Height; ++j)
            for(UINT k=0; k < NoiseTexDesc.Depth; ++k)
            {
                int k0 = (k/InitialStep)*InitialStep;
                int km1 = k0-InitialStep;
                if( km1 < 0 )km1+=NoiseTexDesc.Depth;
                int k1 = (k0+InitialStep) % NoiseTexDesc.Depth;
                int k2 = (k0+2*InitialStep) % NoiseTexDesc.Depth;
                NOISE(i,j,k) = CubicInterpolate(NOISE(i,j,km1), NOISE(i,j,k0), NOISE(i,j,k1), NOISE(i,j,k2), (float)(k-k0) / (float)InitialStep);
            }
    
    // Generate mips
    auto FinerMipIt = NoiseData.begin();
    for(uint Mip = 1; Mip < uiMips; ++Mip)
    {
        UINT uiFinerMipWidth  = NoiseTexDesc.Width  >> (Mip-1);
        UINT uiFinerMipHeight = NoiseTexDesc.Height >> (Mip-1);
        UINT uiFinerMipDepth  = NoiseTexDesc.Depth  >> (Mip-1);

        auto CurrMipIt = FinerMipIt + uiFinerMipWidth * uiFinerMipHeight * uiFinerMipDepth;
        UINT uiMipWidth  = NoiseTexDesc.Width  >> Mip;
        UINT uiMipHeight = NoiseTexDesc.Height >> Mip;
        UINT uiMipDepth  = NoiseTexDesc.Depth  >> Mip;
        for(UINT i=0; i < uiMipWidth; ++i)
            for(UINT j=0; j < uiMipHeight; ++j)
                for(UINT k=0; k < uiMipDepth; ++k)
                {
                    float fVal=0;
                    for(int x=0; x<2;++x)
                        for(int y=0; y<2;++y)
                            for(int z=0; z<2;++z)
                            {
                                fVal += FinerMipIt[(i*2+x) + (j*2 + y) * uiFinerMipWidth + (k*2+z) * (uiFinerMipWidth * uiFinerMipHeight)];
                            }
                    CurrMipIt[i + j * uiMipWidth + k * (uiMipWidth * uiMipHeight)] = fVal / 8.f;
                }
        FinerMipIt = CurrMipIt;
    }
    assert(FinerMipIt+1 == NoiseData.end());

    // Convert to 8-bit
    std::vector<BYTE> NoiseDataR8(NoiseData.size());
    for(auto it=NoiseData.begin(); it != NoiseData.end(); ++it)
        NoiseDataR8[it-NoiseData.begin()] = (BYTE)min(max((int)( *it*255.f), 0),255);

    // Prepare init data
    std::vector<D3D11_SUBRESOURCE_DATA>InitData(uiMips);
    auto CurrMipIt = NoiseDataR8.begin();
    for( UINT Mip = 0; Mip < uiMips; ++Mip )
    {
        UINT uiMipWidth  = NoiseTexDesc.Width  >> Mip;
        UINT uiMipHeight = NoiseTexDesc.Height >> Mip;
        UINT uiMipDepth  = NoiseTexDesc.Depth  >> Mip;
        InitData[Mip].pSysMem = &(*CurrMipIt);
        InitData[Mip].SysMemPitch = uiMipWidth*sizeof(NoiseDataR8[0]);
        InitData[Mip].SysMemSlicePitch = uiMipWidth*uiMipHeight*sizeof(NoiseDataR8[0]);
        CurrMipIt += uiMipWidth * uiMipHeight * uiMipDepth;
    }
    assert(CurrMipIt == NoiseDataR8.end());
    
    // TODO: compress to BC1

    // Create 3D texture
    CComPtr<ID3D11Texture3D> ptex3DNoise;
    V( pDevice->CreateTexture3D(&NoiseTexDesc, &InitData[0], &ptex3DNoise));
    V( pDevice->CreateShaderResourceView(ptex3DNoise, nullptr, &m_ptex3DNoiseSRV));

    return S_OK;
}

HRESULT CCloudsController::OnCreateDevice(ID3D11Device *pDevice, ID3D11DeviceContext *pDeviceContext)
{
    HRESULT hr;

    // Detect and report Intel extensions on this system
    hr = IGFX::Init( pDevice );
	if ( FAILED(hr) )
	{
		//CPUTOSServices::GetOSServices()->OpenMessageBox( _L("Error"), _L("Failed hardware detection initialization: incorrect vendor or device.\n\n") );
	}
    // detect the available extensions
    IGFX::Extensions extensions = IGFX::getAvailableExtensions( pDevice );

    m_bPSOrderingAvailable = extensions.PixelShaderOrdering;  
    
    // Disable the AVSM extension method if the hardware/driver does not support Pixel Shader Ordering feature
    if ( !extensions.PixelShaderOrdering ) 
    {
        CPUTOSServices::GetOSServices()->OpenMessageBox(_L("Pixel Shader Ordering feature not found"), _L("Your hardware or graphics driver does not support the pixel shader ordering feature. Volume-aware blending will be disabled. Please update your driver or run on a system that supports the required feature to see that option."));      
    }

    CreateParticleDataBuffer(pDevice);

    // Create buffer for storing number of valid cells
    {
	    D3D11_BUFFER_DESC ValidCellsCounterBuffDesc = 
        {
            sizeof(UINT)*4,                           //UINT ByteWidth;
            D3D11_USAGE_DEFAULT,                    //D3D11_USAGE Usage;
            D3D11_BIND_SHADER_RESOURCE,             //UINT BindFlags;
            0,                                      //UINT CPUAccessFlags;
            0,                                      //UINT MiscFlags;
            0	            						//UINT StructureByteStride;
        };
        V( CreateBufferAndViews( pDevice, ValidCellsCounterBuffDesc, nullptr, &m_pbufValidCellsCounter) );
        
        D3D11_SHADER_RESOURCE_VIEW_DESC SRVDesc;
        ZeroMemory(&SRVDesc, sizeof(SRVDesc));
        SRVDesc.Format =  DXGI_FORMAT_R32_UINT;
        SRVDesc.ViewDimension = D3D11_SRV_DIMENSION_BUFFER;
        SRVDesc.Buffer.ElementOffset = 0;
        SRVDesc.Buffer.ElementWidth = sizeof(UINT);
        V_RETURN(pDevice->CreateShaderResourceView( m_pbufValidCellsCounter, &SRVDesc, &m_pbufValidCellsCounterSRV));
    }
    
    // Create buffer for storing DispatchIndirect() arguments
    {
        UINT DispatchArgs[] = 
        {
            1, // UINT ThreadGroupCountX
            1, // UINT ThreadGroupCountY
            1, // UINT ThreadGroupCountZ
        };

	    D3D11_BUFFER_DESC DispatchArgsBuffDesc = 
        {
            sizeof(DispatchArgs),                   //UINT ByteWidth;
            D3D11_USAGE_DEFAULT,                    //D3D11_USAGE Usage;
            D3D11_BIND_UNORDERED_ACCESS,            //UINT BindFlags;
            0,                                      //UINT CPUAccessFlags;
            D3D11_RESOURCE_MISC_DRAWINDIRECT_ARGS,  //UINT MiscFlags;
            0                                       //UINT StructureByteStride;
        };

        D3D11_SUBRESOURCE_DATA InitData = {&DispatchArgs, 0, 0};
        V( CreateBufferAndViews( pDevice, DispatchArgsBuffDesc, &InitData, &m_pbufDispatchArgs, nullptr, nullptr) );

        D3D11_UNORDERED_ACCESS_VIEW_DESC UAVDesc;
        UAVDesc.Format = DXGI_FORMAT_R32_UINT;
        UAVDesc.ViewDimension = D3D11_UAV_DIMENSION_BUFFER;
        UAVDesc.Buffer.FirstElement = 0;
        UAVDesc.Buffer.NumElements = _countof(DispatchArgs);
        UAVDesc.Buffer.Flags = 0;
        V_RETURN(pDevice->CreateUnorderedAccessView( m_pbufDispatchArgs, &UAVDesc, &m_pbufDispatchArgsUAV));
    }

    D3D11_BUFFER_DESC GlobalCloudAttribsCBDesc = 
    {
        sizeof(SGlobalCloudAttribs), //UINT ByteWidth;
        D3D11_USAGE_DYNAMIC,         //D3D11_USAGE Usage;
        D3D11_BIND_CONSTANT_BUFFER,  //UINT BindFlags;
        D3D11_CPU_ACCESS_WRITE,      //UINT CPUAccessFlags;
        0,                                      //UINT MiscFlags;
        0                                       //UINT StructureByteStride;
    };
    V(pDevice->CreateBuffer( &GlobalCloudAttribsCBDesc, nullptr, &m_pcbGlobalCloudAttribs));
    
    // Create depth stencil states
    D3D11_DEPTH_STENCIL_DESC EnableDepthTestDSDesc;
    ZeroMemory(&EnableDepthTestDSDesc, sizeof(EnableDepthTestDSDesc));
    EnableDepthTestDSDesc.DepthEnable = TRUE;
    EnableDepthTestDSDesc.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ALL;
    EnableDepthTestDSDesc.DepthFunc = D3D11_COMPARISON_GREATER;
    V( pDevice->CreateDepthStencilState(  &EnableDepthTestDSDesc, &m_pdsEnableDepth) );

    D3D11_DEPTH_STENCIL_DESC DisableDepthTestDSDesc;
    ZeroMemory(&DisableDepthTestDSDesc, sizeof(DisableDepthTestDSDesc));
    DisableDepthTestDSDesc.DepthEnable = FALSE;
    DisableDepthTestDSDesc.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ZERO;
    DisableDepthTestDSDesc.DepthFunc = D3D11_COMPARISON_GREATER;
    V( pDevice->CreateDepthStencilState(  &DisableDepthTestDSDesc, &m_pdsDisableDepth) );
    
    // Create rasterizer states
    D3D11_RASTERIZER_DESC SolidFillCullBackRSDesc;
    ZeroMemory(&SolidFillCullBackRSDesc, sizeof(SolidFillCullBackRSDesc));
    SolidFillCullBackRSDesc.FillMode = D3D11_FILL_SOLID;
    SolidFillCullBackRSDesc.CullMode = D3D11_CULL_FRONT;
    SolidFillCullBackRSDesc.DepthClipEnable = FALSE; // TODO: temp
    V( pDevice->CreateRasterizerState( &SolidFillCullBackRSDesc, &m_prsSolidFillCullFront) );

    D3D11_RASTERIZER_DESC SolidFillNoCullRSDesc;
    ZeroMemory(&SolidFillNoCullRSDesc, sizeof(SolidFillNoCullRSDesc));
    SolidFillNoCullRSDesc.FillMode = D3D11_FILL_SOLID;
    SolidFillNoCullRSDesc.CullMode = D3D11_CULL_NONE;
    SolidFillNoCullRSDesc.DepthClipEnable = TRUE;
    V( pDevice->CreateRasterizerState( &SolidFillNoCullRSDesc, &m_prsSolidFillNoCull) );
   
    // Create default blend state
    D3D11_BLEND_DESC DefaultBlendStateDesc;
    ZeroMemory(&DefaultBlendStateDesc, sizeof(DefaultBlendStateDesc));
    DefaultBlendStateDesc.IndependentBlendEnable = FALSE;
    for(int i=0; i< _countof(DefaultBlendStateDesc.RenderTarget); i++)
        DefaultBlendStateDesc.RenderTarget[i].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;
    V( pDevice->CreateBlendState( &DefaultBlendStateDesc, &m_pbsDefault) );

    // Create blend state for rendering particles
    D3D11_BLEND_DESC AlphaBlendStateDesc;
    ZeroMemory(&AlphaBlendStateDesc, sizeof(AlphaBlendStateDesc));
    AlphaBlendStateDesc.IndependentBlendEnable = TRUE;
    for(int i=0; i< _countof(AlphaBlendStateDesc.RenderTarget); i++)
        AlphaBlendStateDesc.RenderTarget[i].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;
    AlphaBlendStateDesc.RenderTarget[0].BlendEnable = TRUE;
    AlphaBlendStateDesc.RenderTarget[0].BlendOp        = D3D11_BLEND_OP_ADD;
    AlphaBlendStateDesc.RenderTarget[0].SrcBlend       = D3D11_BLEND_ZERO;
    AlphaBlendStateDesc.RenderTarget[0].DestBlend      = D3D11_BLEND_SRC_COLOR;

    AlphaBlendStateDesc.RenderTarget[0].BlendOpAlpha  = D3D11_BLEND_OP_ADD;
    AlphaBlendStateDesc.RenderTarget[0].SrcBlendAlpha = D3D11_BLEND_ZERO;
    AlphaBlendStateDesc.RenderTarget[0].DestBlendAlpha= D3D11_BLEND_SRC_ALPHA;

    AlphaBlendStateDesc.RenderTarget[1].BlendEnable    = TRUE;
    AlphaBlendStateDesc.RenderTarget[1].BlendOp        = D3D11_BLEND_OP_MIN;
    AlphaBlendStateDesc.RenderTarget[1].SrcBlend       = D3D11_BLEND_ONE;
    AlphaBlendStateDesc.RenderTarget[1].DestBlend      = D3D11_BLEND_ONE;
                                     
    AlphaBlendStateDesc.RenderTarget[1].BlendOpAlpha  = D3D11_BLEND_OP_MIN;
    AlphaBlendStateDesc.RenderTarget[1].SrcBlendAlpha = D3D11_BLEND_ONE;
    AlphaBlendStateDesc.RenderTarget[1].DestBlendAlpha= D3D11_BLEND_ONE;

    AlphaBlendStateDesc.RenderTarget[2].BlendEnable = TRUE;
    AlphaBlendStateDesc.RenderTarget[2].BlendOp        = D3D11_BLEND_OP_ADD;
    AlphaBlendStateDesc.RenderTarget[2].SrcBlend       = D3D11_BLEND_ONE;
    AlphaBlendStateDesc.RenderTarget[2].DestBlend      = D3D11_BLEND_SRC_ALPHA;
                                     
    AlphaBlendStateDesc.RenderTarget[2].BlendOpAlpha   = D3D11_BLEND_OP_ADD;
    AlphaBlendStateDesc.RenderTarget[2].SrcBlendAlpha  = D3D11_BLEND_ONE;
    AlphaBlendStateDesc.RenderTarget[2].DestBlendAlpha = D3D11_BLEND_ONE;

    V( pDevice->CreateBlendState( &AlphaBlendStateDesc, &m_pbsRT0MulRT1MinRT2Over) );

    D3D11_SAMPLER_DESC SamLinearWrap = 
    {
        D3D11_FILTER_MIN_MAG_MIP_LINEAR,
        D3D11_TEXTURE_ADDRESS_WRAP,
        D3D11_TEXTURE_ADDRESS_WRAP,
        D3D11_TEXTURE_ADDRESS_WRAP,
        0, //FLOAT MipLODBias;
        0, //UINT MaxAnisotropy;
        D3D11_COMPARISON_NEVER, // D3D11_COMPARISON_FUNC ComparisonFunc;
        {0.f, 0.f, 0.f, 0.f}, //FLOAT BorderColor[ 4 ];
        -FLT_MAX, //FLOAT MinLOD;
        +FLT_MAX //FLOAT MaxLOD;
    };
    V( pDevice->CreateSamplerState( &SamLinearWrap, &m_psamLinearWrap) );

    D3D11_SAMPLER_DESC SamPointWrap = SamLinearWrap;
    SamPointWrap.Filter = D3D11_FILTER_MIN_MAG_MIP_POINT;
    V( pDevice->CreateSamplerState( &SamPointWrap, &m_psamPointWrap) );
    

    SamLinearWrap.AddressU = D3D11_TEXTURE_ADDRESS_CLAMP;
    SamLinearWrap.AddressV = D3D11_TEXTURE_ADDRESS_CLAMP;
    SamLinearWrap.AddressW = D3D11_TEXTURE_ADDRESS_CLAMP;
    V( pDevice->CreateSamplerState( &SamLinearWrap, &m_psamLinearClamp) );

    D3DX11_IMAGE_LOAD_INFO LoadInfo;
    ZeroMemory(&LoadInfo, sizeof(D3DX11_IMAGE_LOAD_INFO));
    LoadInfo.Width          = D3DX11_FROM_FILE;
    LoadInfo.Height         = D3DX11_FROM_FILE;
    LoadInfo.Depth          = D3DX11_FROM_FILE;
    LoadInfo.FirstMipLevel  = D3DX11_FROM_FILE;
    LoadInfo.MipLevels      = D3DX11_DEFAULT;
    LoadInfo.Usage          = D3D11_USAGE_IMMUTABLE;
    LoadInfo.BindFlags      = D3D11_BIND_SHADER_RESOURCE;
    LoadInfo.CpuAccessFlags = 0;
    LoadInfo.MiscFlags      = 0;
    LoadInfo.MipFilter      = D3DX11_FILTER_LINEAR;
    LoadInfo.pSrcInfo       = NULL;
    LoadInfo.Format         = DXGI_FORMAT_BC4_UNORM;
    LoadInfo.Filter         = D3DX11_FILTER_LINEAR;
    
    // Load noise textures. Important to use BC4 compression
    LoadInfo.Format         = DXGI_FORMAT_BC4_UNORM;
    D3DX11CreateShaderResourceViewFromFile(pDevice, L"media\\Noise.png", &LoadInfo, nullptr, &m_ptex2DCloudDensitySRV, nullptr);

    // Noise is not compressed well. Besides, it seems like there are some strange unstable results when using BC1 (?)
    LoadInfo.Format         = DXGI_FORMAT_R8G8B8A8_UNORM;//DXGI_FORMAT_BC1_UNORM;
    D3DX11CreateShaderResourceViewFromFile(pDevice, L"media\\WhiteNoise.png", &LoadInfo, nullptr, &m_ptex2DWhiteNoiseSRV, nullptr);
    
    {
        // Create maximum density mip map
        CComPtr<ID3D11Resource> pCloudDensityRes;
        m_ptex2DCloudDensitySRV->GetResource(&pCloudDensityRes);
        D3D11_TEXTURE2D_DESC CloudDensityTexDesc;
        CComQIPtr<ID3D11Texture2D>(pCloudDensityRes)->GetDesc(&CloudDensityTexDesc);
        m_uiCloudDensityTexWidth = CloudDensityTexDesc.Width;
        m_uiCloudDensityTexHeight = CloudDensityTexDesc.Height;

        D3D11_TEXTURE2D_DESC MaxCloudDensityMipDesc = CloudDensityTexDesc;
        MaxCloudDensityMipDesc.Format = DXGI_FORMAT_R8_UNORM;
        MaxCloudDensityMipDesc.Usage = D3D11_USAGE_DEFAULT;
        MaxCloudDensityMipDesc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
        CComPtr<ID3D11Texture2D> ptex2DMaxDensityMipMap, ptex2DTmpMaxDensityMipMap;
        V(pDevice->CreateTexture2D(&MaxCloudDensityMipDesc, nullptr, &ptex2DMaxDensityMipMap));
        V(pDevice->CreateShaderResourceView(ptex2DMaxDensityMipMap, nullptr, &m_ptex2DMaxDensityMipMapSRV));

        MaxCloudDensityMipDesc.BindFlags = D3D11_BIND_RENDER_TARGET;
        V(pDevice->CreateTexture2D(&MaxCloudDensityMipDesc, nullptr, &ptex2DTmpMaxDensityMipMap));

        RenderMaxDensityMip( pDevice, pDeviceContext, 
                             ptex2DMaxDensityMipMap, ptex2DTmpMaxDensityMipMap, 
                             MaxCloudDensityMipDesc );
    }

    Create3DNoise(pDevice);

    CreateLiSpFirstListIndTex(pDevice);
    
    return S_OK;
}

void CCloudsController::OnDestroyDevice()
{
    m_pcbGlobalCloudAttribs.Release();
    for(int i=0; i < _countof(m_RenderCloudsTech); ++i)
        m_RenderCloudsTech[i].Release();
    for(int i=0; i < _countof(m_RenderFlatCloudsTech); ++i)
        m_RenderFlatCloudsTech[i].Release();
    m_CombineWithBBTech.Release();
    m_RenderCloudDetphToShadowMap.Release();
    m_ProcessCloudGridTech.Release();
    for(int i=0; i < _countof(m_ComputeParticleVisibilityTech); ++i)
        m_ComputeParticleVisibilityTech[i].Release();
    m_ComputeCloudLightingTech.Release();
    m_SmoothParticleLightingTech.Release();
    m_RenderCloudsTiledTech.Release();
    m_PerformTilingTech.Release();
    m_SerializeVisibleParticlesTech.Release();
    m_ProcessParticlesTech.Release();
    m_ComputeDispatchArgsTech.Release();
    m_ComputeOpticalDepthTech.Release();
    m_ApplyParticleLayersTech.Release();
    m_ComputeSingleSctrInParticleTech.Release();
    m_GatherPrevSctrOrderTech.Release();
    m_ComputeScatteringOrderTech.Release();
    m_AccumulateInscatteringTech.Release();
    m_RenderScatteringLUTSliceTech.Release();
    m_ComputeExitanceTech.Release();

    m_pdsEnableDepth.Release();
    m_pdsDisableDepth.Release();
    m_prsSolidFillCullFront.Release();
    m_prsSolidFillNoCull.Release();
    m_pbsDefault.Release();
    m_pbsRT0MulRT1MinRT2Over.Release();
    m_ptex2DCloudDensitySRV.Release();
    m_ptex2DWhiteNoiseSRV.Release();
    m_ptex2DMaxDensityMipMapSRV.Release();
    m_ptex3DNoiseSRV.Release();
    m_psamLinearWrap.Release();
    m_psamPointWrap.Release();
    m_psamLinearClamp.Release();

    m_pbufCloudGridSRV.Release();
    m_pbufCloudGridUAV.Release();
    m_pbufCloudParticlesUAV.Release();
    m_pbufCloudParticlesSRV.Release();
    m_pbufVisibleParticlesFlagsSRV.Release();
    m_pbufVisibleParticlesFlagsUAV.Release();
    m_pbufParticlesLightingSRV.Release();
    m_pbufParticlesLightingUAV.Release();
    m_pbufAttenuatedSunLightSRV.Release();
    m_pbufAttenuatedSunLightUAV.Release();
    m_pbufDefaultAttenuatedSunLight.Release();
    m_pbufValidCellsUnorderedList.Release();
    m_pbufValidCellsUnorderedListUAV.Release();
    m_pbufValidCellsUnorderedListSRV.Release();
    m_pbufValidCellsCounter.Release();
    m_pbufValidCellsCounterSRV.Release();
    m_pbufValidParticlesUnorderedList.Release();
    m_pbufValidParticlesUnorderedListUAV.Release();
    m_pbufValidParticlesUnorderedListSRV.Release();

    m_pbufSortedParticlesOrder.Release();
    m_pbufSortedParticlesOrderSRV.Release();
    
    m_pbufSerializedVisibleParticles.Release();
    m_pbufSerializedVisibleParticlesSRV.Release();
    m_pbufDispatchArgsUAV.Release();
    m_pbufDispatchArgs.Release();
	
    m_pbufPackedCellLocationsSRV.Release();

    m_ptex2DScreenCloudColorSRV.Release();
    m_ptex2DScreenCloudColorRTV.Release();
    m_ptex2DScrSpaceCloudTransparencySRV.Release();
    m_ptex2DScrSpaceCloudTransparencyRTV.Release();
    m_ptex2DScrSpaceDistToCloudSRV.Release();
    m_ptex2DScrSpaceDistToCloudRTV.Release();

    m_ptex2DDownscaledScrCloudColorSRV.Release();
    m_ptex2DDownscaledScrCloudColorRTV.Release();
    m_ptex2DDownscaledScrCloudTransparencySRV.Release();
    m_ptex2DDownscaledScrCloudTransparencyRTV.Release();
    m_ptex2DDownscaledScrDistToCloudSRV.Release();
    m_ptex2DDownscaledScrDistToCloudRTV.Release();

    m_pbufParticleListsBuffUAV.Release();
    m_pbufParticleListsBuffSRV.Release();

    m_ptex2DScrFirstListIndUAV.Release();
    m_ptex2DScrFirstListIndSRV.Release();

    m_ptex2DLiSpFirstListIndUAV.Release();
    m_ptex2DLiSpFirstListIndSRV.Release();

    m_pbufParticleLayersSRV.Release();
    m_pbufParticleLayersUAV.Release();
    m_pbufClearParticleLayers.Release();

    m_ptex3DPrecomputedParticleDensitySRV.Release();
    m_ptex3DSingleSctrInParticleLUT_SRV.Release();
    m_ptex3DMultipleSctrInParticleLUT_SRV.Release();

    
    m_pRenderCloudsInputLayout.Release();
}

// Method computes visibility flags for all valid particles
void CCloudsController::ComputeParticleVisibility(ID3D11Device *pDevice,
                                                  ID3D11DeviceContext *pDeviceContext, 
                                                  ID3D11Buffer *pcbCameraAttribs,
                                                  bool bLightSpacePass)
{
    auto &ComputeParticleVisibilityTech = m_ComputeParticleVisibilityTech[bLightSpacePass ? 1 : 0];
    if( !ComputeParticleVisibilityTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.AddShaderMacro("THREAD_GROUP_SIZE", sm_iCSThreadGroupSize);
        Macros.AddShaderMacro("LIGHT_SPACE_PASS", bLightSpacePass);
        Macros.Finalize();

        ComputeParticleVisibilityTech.SetDeviceAndContext(pDevice, pDeviceContext);
        ComputeParticleVisibilityTech.CreateComputeShaderFromFile(m_strEffectPath, "ComputeParticlesVisibilityCS", Macros);
    }

    ID3D11Buffer *pCBs[] = {m_pcbGlobalCloudAttribs, nullptr, pcbCameraAttribs};
    pDeviceContext->CSSetConstantBuffers(0, _countof(pCBs), pCBs);

    //ID3D11SamplerState *pSamplers[] = {m_psamLinearClamp, m_psamLinearWrap, m_psamPointWrap};
    //pDeviceContext->CSSetSamplers(0, _countof(pSamplers), pSamplers);

    ID3D11ShaderResourceView *pSRVs[] = 
    {
        m_pbufValidCellsCounterSRV,       // Buffer<uint> g_ValidCellsCounter                 : register( t0 );
        m_pbufValidParticlesUnorderedListSRV, // StructuredBuffer<uint> g_ValidParticlesUnorderedList : register( t1 );
        m_pbufCloudGridSRV,               // StructuredBuffer<SCloudCellAttribs> g_CloudCells : register( t2 );
        m_pbufCloudParticlesSRV           // StructuredBuffer<SParticleAttribs> g_Particles : register( t3 );
    };
    pDeviceContext->CSSetShaderResources(0, _countof(pSRVs), pSRVs);

	ID3D11UnorderedAccessView *pUAVs[] = {m_pbufVisibleParticlesFlagsUAV};
	pDeviceContext->CSSetUnorderedAccessViews(0, _countof(pUAVs), pUAVs, NULL);
	UINT Zero[] = {0,0,0,0};
	pDeviceContext->ClearUnorderedAccessViewUint(m_pbufVisibleParticlesFlagsUAV, Zero);

    ComputeParticleVisibilityTech.Apply();
    pDeviceContext->DispatchIndirect(m_pbufDispatchArgs, 0);
    
    memset(pUAVs, 0, sizeof(pUAVs));
    pDeviceContext->CSSetUnorderedAccessViews(0, _countof(pUAVs), pUAVs, NULL);
}

// Method processes all particles and streams out visible ones preserving submit order
void CCloudsController::SerializeVisibileParticles(ID3D11Device *pDevice,
                                                   ID3D11DeviceContext *pDeviceContext,
                                                   const std::vector<UINT> &ParticleOrder)
{
    static const UINT uiBatchSize = 16;
    if( !m_SerializeVisibleParticlesTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.AddShaderMacro("SERIALIZE_PARTICLE_GS_BATCH_SIZE", uiBatchSize);
        Macros.Finalize();

        m_SerializeVisibleParticlesTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_SerializeVisibleParticlesTech.CreateVertexShaderFromFile(m_strEffectPath, "SerializeVisibleParticlesVS", Macros);
        D3D11_SO_DECLARATION_ENTRY SODecl[] = 
        {
            {0, "PARTICLE_SERIAL_NUM", 0, 0, 1, 0}
        };
        UINT Strides[] = {sizeof(UINT)};
        m_SerializeVisibleParticlesTech.CreateGSWithSOFromFile(m_strEffectPath, "SerializeVisibleParticlesGS", Macros, SODecl, _countof(SODecl), Strides, _countof(Strides), D3D11_SO_NO_RASTERIZED_STREAM );
    }

    // Upload sorted particle index
    D3D11_MAPPED_SUBRESOURCE MapData;
    pDeviceContext->Map(m_pbufSortedParticlesOrder, 0, D3D11_MAP_WRITE_DISCARD, 0, &MapData);
    memcpy(MapData.pData, &ParticleOrder[0], sizeof(ParticleOrder[0])*ParticleOrder.size());
    pDeviceContext->Unmap(m_pbufSortedParticlesOrder, 0);

    ID3D11Buffer *pBuffs[] = {m_pbufSerializedVisibleParticles};
    UINT Offsets[_countof(pBuffs)] = { 0 };
    pDeviceContext->SOSetTargets(_countof(pBuffs), pBuffs, Offsets);
    // StructuredBuffer<uint> g_SortedParticlesOrder : register( t0 );
    pDeviceContext->GSSetShaderResources(0, 1, &m_pbufSortedParticlesOrderSRV.p);
    // StructuredBuffer<uint> g_VisibleParticleFlags : register(t6);
    pDeviceContext->GSSetShaderResources(6, 1, &m_pbufVisibleParticlesFlagsSRV.p);
    m_SerializeVisibleParticlesTech.Apply();
    pDeviceContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
    ID3D11Buffer *pDummyVBs[] = { nullptr };
    UINT DummyStrides[_countof(pDummyVBs)] = { 0 };
    UINT DummyOffsets[_countof(pDummyVBs)] = { 0 };
    pDeviceContext->IASetVertexBuffers(0 , 1, pDummyVBs, DummyStrides, DummyOffsets);
    pDeviceContext->IASetInputLayout(nullptr);
    // Round up the number of groups. All particles behind the range are anyway invisible
    pDeviceContext->Draw( ((UINT)ParticleOrder.size() + (uiBatchSize-1))/uiBatchSize, 0);
    memset(pBuffs, 0, sizeof(pBuffs));
    pDeviceContext->SOSetTargets(_countof(pBuffs), pBuffs, Offsets);
}


void CCloudsController::Update( const SGlobalCloudAttribs &NewAttribs,
                                const D3DXVECTOR3 &CameraPos, 
                                const D3DXVECTOR3 &LightDir,
                                ID3D11Device *pDevice, 
                                ID3D11DeviceContext *pDeviceContext, 
                                ID3D11Buffer *pcbCameraAttribs, 
                                ID3D11Buffer *pcbLightAttribs, 
                                ID3D11Buffer *pcMediaScatteringParams)
{
    if(GetAsyncKeyState(VK_F7))
    {
        for(int i=0; i < _countof(m_RenderCloudsTech); ++i)
            m_RenderCloudsTech[i].Release();
        for(int i=0; i < _countof(m_RenderFlatCloudsTech); ++i)
            m_RenderFlatCloudsTech[i].Release();
        m_CombineWithBBTech.Release();
        m_RenderCloudDetphToShadowMap.Release();
        m_ProcessCloudGridTech.Release();
        for(int i=0; i < _countof(m_ComputeParticleVisibilityTech); ++i)
            m_ComputeParticleVisibilityTech[i].Release();
        m_ComputeCloudLightingTech.Release();
        m_SmoothParticleLightingTech.Release();
        m_RenderCloudsTiledTech.Release();
        m_PerformTilingTech.Release();
        m_SerializeVisibleParticlesTech.Release();
        m_ProcessParticlesTech.Release();
        m_ComputeDispatchArgsTech.Release();
        m_ComputeOpticalDepthTech.Release();
        m_ApplyParticleLayersTech.Release();
        m_ComputeSingleSctrInParticleTech.Release();
        m_GatherPrevSctrOrderTech.Release();
        m_ComputeScatteringOrderTech.Release();
        m_AccumulateInscatteringTech.Release();
        m_RenderScatteringLUTSliceTech.Release();
        m_ComputeExitanceTech.Release();
        m_ptex3DPrecomputedParticleDensitySRV.Release();
        m_ptex3DSingleSctrInParticleLUT_SRV.Release();
        m_ptex3DMultipleSctrInParticleLUT_SRV.Release();
    }

    m_CloudAttribs.fCloudDensityThreshold = NewAttribs.fCloudDensityThreshold;
    m_CloudAttribs.fCloudAltitude         = NewAttribs.fCloudAltitude;
    m_CloudAttribs.fCloudThickness        = NewAttribs.fCloudThickness;
    m_CloudAttribs.fParticleCutOffDist    = NewAttribs.fParticleCutOffDist;

    if( m_CloudAttribs.uiNumRings     != NewAttribs.uiNumRings ||
        m_CloudAttribs.uiInnerRingDim != NewAttribs.uiInnerRingDim ||
        m_CloudAttribs.uiMaxLayers    != NewAttribs.uiMaxLayers )
    {
        m_CloudAttribs.uiNumRings     = NewAttribs.uiNumRings;
        m_CloudAttribs.uiInnerRingDim = NewAttribs.uiInnerRingDim;
        m_CloudAttribs.uiMaxLayers    = NewAttribs.uiMaxLayers;
        CreateParticleDataBuffer(pDevice);
        m_f3PrevLightDir = D3DXVECTOR3(0,0,0);
    }

    if( m_CloudAttribs.uiDownscaleFactor != NewAttribs.uiDownscaleFactor )
    {
        m_CloudAttribs.uiDownscaleFactor = NewAttribs.uiDownscaleFactor;
        OnResize(pDevice, m_uiBackBufferWidth, m_uiBackBufferHeight);
        for(int i=0; i < _countof(m_RenderCloudsTech); ++i)
            m_RenderCloudsTech[i].Release();
        for(int i=0; i < _countof(m_RenderFlatCloudsTech); ++i)
            m_RenderFlatCloudsTech[i].Release();
        m_CombineWithBBTech.Release();
    }

    if( m_CloudAttribs.uiNumCascades != NewAttribs.uiNumCascades )
    {
        m_CloudAttribs.uiNumCascades = NewAttribs.uiNumCascades;
        for(int i=0; i < _countof(m_RenderCloudsTech); ++i)
            m_RenderCloudsTech[i].Release();
        for(int i=0; i < _countof(m_RenderFlatCloudsTech); ++i)
            m_RenderFlatCloudsTech[i].Release();
        m_RenderCloudsTiledTech.Release();
    }

    if( m_CloudAttribs.bVolumetricBlending != NewAttribs.bVolumetricBlending )
    {
        m_CloudAttribs.bVolumetricBlending = NewAttribs.bVolumetricBlending;
        m_RenderCloudsTech[0].Release();
    }

    if( m_CloudAttribs.uiDensityGenerationMethod != NewAttribs.uiDensityGenerationMethod )
    {
        m_CloudAttribs.uiDensityGenerationMethod = NewAttribs.uiDensityGenerationMethod;
        m_ptex3DPrecomputedParticleDensitySRV.Release();
        m_ComputeOpticalDepthTech.Release();
    }

    const UINT iNumRings = m_CloudAttribs.uiNumRings;
    const UINT iRingDimension = m_CloudAttribs.uiRingDimension;
    const UINT iMaxLayers = m_CloudAttribs.uiMaxLayers;

    int iHalfRingDim = m_CloudAttribs.uiRingDimension/2;
    const float fLargestRingSize = m_CloudAttribs.fParticleCutOffDist * 2;

    // Compute aux attribs for each ring to simplify coordinate computations
    for(UINT iRing = 0; iRing < m_CloudAttribs.uiNumRings; ++iRing)
    {
        const float fRingWorldStep = fLargestRingSize / (float)((m_CloudAttribs.uiRingDimension) << ((m_CloudAttribs.uiNumRings-1) - iRing));
        float fCameraI = floor(CameraPos.x/fRingWorldStep + 0.5f);
        float fCameraJ = floor(CameraPos.z/fRingWorldStep + 0.5f);
        float fRingStartX = (fCameraI - (float)iHalfRingDim + 0.5f) * fRingWorldStep - CameraPos.x;
        float fRingStartY = (fCameraJ - (float)iHalfRingDim + 0.5f) * fRingWorldStep - CameraPos.z;
        m_RingGridAttribs[iRing] = D3DXVECTOR3(fRingStartX, fRingStartY, fRingWorldStep);
    }
    // Compute distance to camera and light space Z
    auto DistToCamIt = m_ParticleDistToCam.begin();
    auto LiSpZIt = m_ParticleLightSpaceZ.begin();
    for(auto ParticleIt = m_PackedParticleLocations.begin(); ParticleIt != m_PackedParticleLocations.end(); ++ParticleIt, ++DistToCamIt, ++LiSpZIt)
    {
        UINT i,j,ring,layer;
        UnPackParticleIJRing( *ParticleIt, i,j,ring,layer );
        
        UINT uiActiveLayers = GetNumActiveLayers(iMaxLayers, ring);
        if( layer >= uiActiveLayers )
            continue;

        D3DXVECTOR3 &RingGridAttribs = m_RingGridAttribs[ring];
        const float fRingWorldStep = RingGridAttribs.z;
        
        assert( fRingWorldStep == (float)(fLargestRingSize / (float)((m_CloudAttribs.uiRingDimension) << ((m_CloudAttribs.uiNumRings-1) - ring))) );

        D3DXVECTOR3 f3Pos;
        float fCameraI = floor(CameraPos.x/fRingWorldStep + 0.5f);
        float fCameraJ = floor(CameraPos.z/fRingWorldStep + 0.5f);
        assert( RingGridAttribs.x == (float)((fCameraI + (float)((int)0 - iHalfRingDim) + 0.5f) * fRingWorldStep - CameraPos.x) );
        assert( RingGridAttribs.y == (float)((fCameraJ + (float)((int)0 - iHalfRingDim) + 0.5f) * fRingWorldStep - CameraPos.z) );
        //f3Pos.x = (fCameraI + (float)((int)i - iHalfRingDim) + 0.5f) * fRingWorldStep;
        //f3Pos.z = (fCameraJ + (float)((int)j - iHalfRingDim) + 0.5f) * fRingWorldStep;
        f3Pos.x = RingGridAttribs.x + (float)i * fRingWorldStep;
        f3Pos.z = RingGridAttribs.y + (float)j * fRingWorldStep;
        f3Pos.y = m_CloudAttribs.fCloudAltitude + static_cast<float>( (-0.5 + ((float)layer + 0.5f) / uiActiveLayers) * m_CloudAttribs.fCloudThickness) - CameraPos.y;
        float fDistToCam = D3DXVec3Dot(&f3Pos, &f3Pos);
        float fLightSpaceZ = D3DXVec3Dot(&f3Pos, &LightDir);
        *DistToCamIt = fDistToCam;
        *LiSpZIt = fLightSpaceZ;
    }

    {
        // Sort particles
        std::sort( m_CamSpaceOrder.begin(), m_CamSpaceOrder.end(), 
                    [&](UINT uiInd0, UINT uiInd1)
                    {
                        return m_ParticleDistToCam[uiInd0] > m_ParticleDistToCam[uiInd1];
                    }
                  );
    }

    if( m_f3PrevLightDir != LightDir )
    {
        m_f3PrevLightDir = LightDir;
        std::sort( m_LightSpaceOrder.begin(), m_LightSpaceOrder.end(), 
                    [&](UINT uiInd0, UINT uiInd1)
                    {
                        return m_ParticleLightSpaceZ[uiInd0] > m_ParticleLightSpaceZ[uiInd1];
                    }
                  );
    }

    // Process cloud grid
    if( !m_ProcessCloudGridTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.AddShaderMacro("THREAD_GROUP_SIZE", sm_iCSThreadGroupSize);
        Macros.Finalize();

        m_ProcessCloudGridTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_ProcessCloudGridTech.CreateComputeShaderFromFile(m_strEffectPath, "ProcessCloudGridCS", Macros);
    }
    
    UpdateConstantBuffer(pDeviceContext, m_pcbGlobalCloudAttribs, &m_CloudAttribs, sizeof(m_CloudAttribs));

    ID3D11Buffer *pCBs[] = {m_pcbGlobalCloudAttribs, pcMediaScatteringParams, pcbCameraAttribs, pcbLightAttribs};
    pDeviceContext->CSSetConstantBuffers(0, _countof(pCBs), pCBs);

    ID3D11SamplerState *pSamplers[] = {m_psamLinearClamp, m_psamLinearWrap, m_psamPointWrap};
    pDeviceContext->CSSetSamplers(0, _countof(pSamplers), pSamplers);

    ID3D11ShaderResourceView *pSRVs[] = 
    {
        m_pbufPackedCellLocationsSRV, // StructuredBuffer<uint> g_PackedCellLocations : register( t0 );
        m_ptex2DCloudDensitySRV,
        m_ptex3DNoiseSRV,
        m_ptex2DMaxDensityMipMapSRV,
        nullptr,
        nullptr, // Texture2D<float2> g_tex2DOccludedNetDensityToAtmTop : register( t5 );
    };
    pDeviceContext->CSSetShaderResources(0, _countof(pSRVs), pSRVs);
    
    ID3D11UnorderedAccessView *pUAVs[] = {m_pbufCloudGridUAV, m_pbufValidCellsUnorderedListUAV};
    UINT uiZeroCounters[_countof(pUAVs)] =  {0};
    pDeviceContext->CSSetUnorderedAccessViews(0, _countof(pUAVs), pUAVs, uiZeroCounters);

    m_ProcessCloudGridTech.Apply();
    pDeviceContext->Dispatch( (m_CloudAttribs.uiNumCells + (sm_iCSThreadGroupSize-1)) / sm_iCSThreadGroupSize, 1, 1);
    
    memset(pUAVs, 0, sizeof(pUAVs));
    pDeviceContext->CSSetUnorderedAccessViews(0, _countof(pUAVs), pUAVs, nullptr);

    // Compute DispatchIndirect() arguments
    if( !m_ComputeDispatchArgsTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.AddShaderMacro("THREAD_GROUP_SIZE", sm_iCSThreadGroupSize);
        Macros.Finalize();

        m_ComputeDispatchArgsTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_ComputeDispatchArgsTech.CreateComputeShaderFromFile(m_strEffectPath, "ComputeDispatchArgsCS", Macros);
    }

    pDeviceContext->CopyStructureCount(m_pbufValidCellsCounter, 0, m_pbufValidCellsUnorderedListUAV);
    pSRVs[0] = m_pbufValidCellsCounterSRV; // Buffer<uint> g_ValidCellsCounter : register( t0 );
    pDeviceContext->CSSetShaderResources(0, 1, pSRVs);
    pUAVs[0] = m_pbufDispatchArgsUAV;
    pDeviceContext->CSSetUnorderedAccessViews(0, 1, pUAVs, nullptr);
    m_ComputeDispatchArgsTech.Apply();
    pDeviceContext->Dispatch(1,1,1);
    pUAVs[0] = nullptr;
    pDeviceContext->CSSetUnorderedAccessViews(0, 1, pUAVs, nullptr);
    
    // Process all valid cells and generate particles for each cell
    if(!m_ProcessParticlesTech.IsValid())
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.AddShaderMacro("THREAD_GROUP_SIZE", sm_iCSThreadGroupSize);
        Macros.Finalize();

        m_ProcessParticlesTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_ProcessParticlesTech.CreateComputeShaderFromFile(m_strEffectPath, "ProcessValidParticlesCS", Macros);
    }

    pUAVs[0] = m_pbufCloudParticlesUAV;
    pUAVs[1] = m_pbufValidParticlesUnorderedListUAV;
    pDeviceContext->CSSetUnorderedAccessViews(0, 2, pUAVs, uiZeroCounters);
    pSRVs[0] = m_pbufValidCellsCounterSRV;       // Buffer<uint> g_ValidCellsCounter                 : register( t0 );
    pSRVs[1] = m_pbufValidCellsUnorderedListSRV; // StructuredBuffer<uint> g_ValidCellsUnorderedList : register( t1 );
    pSRVs[2] = m_pbufCloudGridSRV;               // StructuredBuffer<SCloudCellAttribs> g_CloudCells : register( t2 );
    pSRVs[3] = m_ptex2DWhiteNoiseSRV;            // Texture2D<float3> g_tex2DWhiteNoise              : register( t3 );
    pDeviceContext->CSSetShaderResources(0, 4, pSRVs);
    m_ProcessParticlesTech.Apply();
    pDeviceContext->DispatchIndirect(m_pbufDispatchArgs, 0);
    memset(pUAVs, 0, sizeof(pUAVs));
    pDeviceContext->CSSetUnorderedAccessViews(0, 2, pUAVs, nullptr);


    pDeviceContext->CopyStructureCount(m_pbufValidCellsCounter, 0, m_pbufValidParticlesUnorderedListUAV);
    pSRVs[0] = m_pbufValidCellsCounterSRV; // Buffer<uint> g_ValidCellsCounter : register( t0 );
    pDeviceContext->CSSetShaderResources(0, 1, pSRVs);
    pUAVs[0] = m_pbufDispatchArgsUAV;
    pDeviceContext->CSSetUnorderedAccessViews(0, 1, pUAVs, nullptr);
    m_ComputeDispatchArgsTech.Apply();
    pDeviceContext->Dispatch(1,1,1);
    pUAVs[0] = nullptr;
    pDeviceContext->CSSetUnorderedAccessViews(0, 1, pUAVs, nullptr);
}

void CCloudsController::PerformTilingInLightSpace(SRenderAttribs &RenderAttribs)
{
    ID3D11DeviceContext *pDeviceContext = RenderAttribs.pDeviceContext;
    ID3D11Device *pDevice = RenderAttribs.pDevice;

    D3DXMATRIX WorldToLightViewSpaceMatr;
    D3DXMatrixTranspose( &WorldToLightViewSpaceMatr, &RenderAttribs.m_pSMAttribs->mWorldToLightViewT );

    D3DXMATRIX mCameraView, mCameraProj, ParticleCutOffProjMatr;
    D3DXMatrixTranspose( &mCameraProj, &RenderAttribs.m_pCameraAttribs->mProjT );
    D3DXMatrixTranspose( &mCameraView, &RenderAttribs.m_pCameraAttribs->mViewT );
    float fMainCamNearPlane = -mCameraProj._43 / mCameraProj._33;
    float fMainCamFarPlane = mCameraProj._33 / (mCameraProj._33-1) * fMainCamNearPlane;
    // Remember that complimentary depth buffering is used
    float fFarPlane  = fMainCamFarPlane;
    float fNearPlane = min(m_CloudAttribs.fParticleCutOffDist, fMainCamNearPlane);
    ParticleCutOffProjMatr = mCameraProj;
    ParticleCutOffProjMatr._33 = fFarPlane / (fFarPlane - fNearPlane);
    ParticleCutOffProjMatr._43 = -fNearPlane * ParticleCutOffProjMatr._33;

    D3DXMATRIX ParticleCutOffFrustumMatrix = mCameraView * ParticleCutOffProjMatr;
    D3DXMATRIX ParticleCutOffFrustumMatrixInv;
    D3DXMatrixInverse(&ParticleCutOffFrustumMatrixInv, nullptr, &ParticleCutOffFrustumMatrix);

    D3DXMATRIX ParticleCutOffProjToLightSpaceMatr = ParticleCutOffFrustumMatrixInv * WorldToLightViewSpaceMatr;

    D3DXVECTOR3 f3MinXYZ(+FLT_MAX,+FLT_MAX,+FLT_MAX), f3MaxXYZ(-FLT_MAX,-FLT_MAX,-FLT_MAX);
    for(int iClipPlaneCorner=0; iClipPlaneCorner < 8; ++iClipPlaneCorner)
    {
        D3DXVECTOR3 f3PlaneCornerProjSpace( (iClipPlaneCorner & 0x01) ? +1.f : - 1.f, 
                                            (iClipPlaneCorner & 0x02) ? +1.f : - 1.f,
                                            // Since we use complimentary depth buffering, 
                                            // far plane has depth 0
                                            (iClipPlaneCorner & 0x04) ? 1.f : 0.f);
        D3DXVECTOR3 f3PlaneCornerLightSpace;
        D3DXVec3TransformCoord(&f3PlaneCornerLightSpace, &f3PlaneCornerProjSpace, &ParticleCutOffProjToLightSpaceMatr);
        D3DXVec3Minimize(&f3MinXYZ, &f3MinXYZ, &f3PlaneCornerLightSpace);
        D3DXVec3Maximize(&f3MaxXYZ, &f3MaxXYZ, &f3PlaneCornerLightSpace);
    }

    D3DXMATRIX ParticleTilingProjMatr;
    D3DXMatrixOrthoOffCenterLH( &ParticleTilingProjMatr, f3MinXYZ.x, f3MaxXYZ.x, f3MinXYZ.y, f3MaxXYZ.y, f3MaxXYZ.z, f3MinXYZ.z);

    D3DXMATRIX ParticleLiSpTilingMatr = WorldToLightViewSpaceMatr * ParticleTilingProjMatr;
    
    ExtractViewFrustumPlanesFromMatrix(ParticleLiSpTilingMatr, (SViewFrustum&)m_CloudAttribs.f4TilingFrustumPlanes);
    UpdateConstantBuffer(pDeviceContext, m_pcbGlobalCloudAttribs, &m_CloudAttribs, sizeof(m_CloudAttribs));

    // Compute visibility to serialize only particles affecting lighting computations
    ComputeParticleVisibility(pDevice, pDeviceContext, nullptr, true);

    SerializeVisibileParticles(pDevice, pDeviceContext, m_LightSpaceOrder);
    
	UINT MinusOne[] = {(UINT)-1, (UINT)-1, (UINT)-1, (UINT)-1};
	pDeviceContext->ClearUnorderedAccessViewUint(m_ptex2DLiSpFirstListIndUAV, MinusOne);

    PerformTiling(RenderAttribs, m_ptex2DLiSpFirstListIndUAV, ParticleLiSpTilingMatr, m_CloudAttribs.uiLiSpFirstListIndTexDim, m_CloudAttribs.uiLiSpFirstListIndTexDim);
}

// Method computes lighting for each valid visible particle
void CCloudsController::ComputeParticleLighting(SRenderAttribs &RenderAttribs)
{
    ID3D11DeviceContext *pDeviceContext = RenderAttribs.pDeviceContext;
    
    if( !m_ComputeCloudLightingTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.AddShaderMacro("THREAD_GROUP_SIZE", sm_iCSThreadGroupSize);
        Macros.Finalize();

        m_ComputeCloudLightingTech.SetDeviceAndContext(RenderAttribs.pDevice, RenderAttribs.pDeviceContext);
        m_ComputeCloudLightingTech.CreateComputeShaderFromFile(m_strEffectPath, "ComputeParticlesLightingCS", Macros);
    }

    if( !m_SmoothParticleLightingTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.AddShaderMacro("THREAD_GROUP_SIZE", sm_iCSThreadGroupSize);
        Macros.Finalize();

        m_SmoothParticleLightingTech.SetDeviceAndContext(RenderAttribs.pDevice, RenderAttribs.pDeviceContext);
        m_SmoothParticleLightingTech.CreateComputeShaderFromFile(m_strEffectPath, "SmoothParticlesLightingCS", Macros);
    }

    ID3D11Buffer *pCBs[] = {m_pcbGlobalCloudAttribs, RenderAttribs.pcMediaScatteringParams, RenderAttribs.pcbCameraAttribs, RenderAttribs.pcbLightAttribs};
    pDeviceContext->CSSetConstantBuffers(0, _countof(pCBs), pCBs);

    ID3D11SamplerState *pSamplers[] = {m_psamLinearClamp, m_psamLinearWrap, m_psamPointWrap};
    pDeviceContext->CSSetSamplers(0, _countof(pSamplers), pSamplers);

    ID3D11ShaderResourceView *pSRVs[] = 
    {
        m_pbufValidCellsCounterSRV,              // Buffer<uint> g_ValidCellsCounter                 : register( t0 );
        m_pbufValidParticlesUnorderedListSRV,    // StructuredBuffer<uint> g_ValidParticlesUnorderedList : register( t1 );
        m_pbufCloudGridSRV,                      // StructuredBuffer<SCloudCellAttribs> g_CloudCells : register( t2 );
        m_pbufCloudParticlesSRV,                 // StructuredBuffer<SParticleAttribs> g_Particles : register( t3 );
        nullptr,
	    RenderAttribs.pPrecomputedNetDensitySRV, // Texture2D<float2> g_tex2DOccludedNetDensityToAtmTop : register( t5 );
        m_pbufVisibleParticlesFlagsSRV,          // StructuredBuffer<uint> g_VisibleParticleFlags : register(t6);
        RenderAttribs.pAmbientSkylightSRV,       // t7
        m_ptex2DLiSpFirstListIndSRV,           // Texture2D<int> g_tex2DFirstKnot : register( t8 );
        m_pbufParticleListsBuffSRV            // StructuredBuffer<SParticleListKnot> g_TiledParticlesList : register( t10 );
	};
    pDeviceContext->CSSetShaderResources(0, _countof(pSRVs), pSRVs);

    CComPtr<ID3D11Resource> pRes;
    m_pbufAttenuatedSunLightSRV->GetResource(&pRes);
    pDeviceContext->CopyResource(pRes, m_pbufDefaultAttenuatedSunLight);

    ID3D11UnorderedAccessView *pUAVs[] = {m_pbufParticlesLightingUAV, m_pbufAttenuatedSunLightUAV};
	pDeviceContext->CSSetUnorderedAccessViews(0, _countof(pUAVs), pUAVs, NULL);

    m_ComputeCloudLightingTech.Apply();
    pDeviceContext->DispatchIndirect(m_pbufDispatchArgs, 0);

    pUAVs[0] = m_pbufParticlesLightingUAV;
    pUAVs[1] = nullptr;
	pDeviceContext->CSSetUnorderedAccessViews(0, _countof(pUAVs), pUAVs, NULL);

    pSRVs[4] = m_pbufAttenuatedSunLightSRV;
    pDeviceContext->CSSetShaderResources(0, _countof(pSRVs), pSRVs);

    m_SmoothParticleLightingTech.Apply();
    pDeviceContext->DispatchIndirect(m_pbufDispatchArgs, 0);

    ID3D11UnorderedAccessView *pDummyUAV[]={NULL, NULL};
    pDeviceContext->CSSetUnorderedAccessViews(0, _countof(pDummyUAV), pDummyUAV, NULL);

}

void CCloudsController :: DefineMacros(class CD3DShaderMacroHelper &Macros)
{
    {
        std::stringstream ss;
        ss<<"float2("<<m_uiCloudDensityTexWidth<<","<<m_uiCloudDensityTexHeight<<")";
        Macros.AddShaderMacro("CLOUD_DENSITY_TEX_DIM", ss.str());
    }
    {
        std::stringstream ss;
        ss<<"float4("<< m_PrecomputedOpticalDepthTexDim.iNumStartPosZenithAngles  <<","
                     << m_PrecomputedOpticalDepthTexDim.iNumStartPosAzimuthAngles <<","
                     << m_PrecomputedOpticalDepthTexDim.iNumDirectionZenithAngles <<","
                     << m_PrecomputedOpticalDepthTexDim.iNumDirectionAzimuthAngles<< ")";
        Macros.AddShaderMacro("OPTICAL_DEPTH_LUT_DIM", ss.str());
    }

    Macros.AddShaderMacro("NUM_PARTICLE_LAYERS", m_CloudAttribs.uiNumParticleLayers);
    Macros.AddShaderMacro("PS_ORDERING_AVAILABLE", m_bPSOrderingAvailable);

    {
        std::stringstream ss;
        ss<<"float4("<< m_PrecomputedSctrInParticleLUTDim.iNumStartPosZenithAngles <<","
                     << m_PrecomputedSctrInParticleLUTDim.iNumViewDirAzimuthAngles <<","
                     << m_PrecomputedSctrInParticleLUTDim.iNumViewDirZenithAngles <<","
                     << m_PrecomputedSctrInParticleLUTDim.iNumDistancesFromCenter << ")";
        Macros.AddShaderMacro("VOL_SCATTERING_IN_PARTICLE_LUT_DIM", ss.str());
    }
    {
        std::stringstream ss;
        ss<<"float4("<< m_PrecomputedSctrInParticleLUTDim.iNumStartPosZenithAngles <<","
                     << m_PrecomputedSctrInParticleLUTDim.iNumViewDirAzimuthAngles <<","
                     << m_PrecomputedSctrInParticleLUTDim.iNumViewDirZenithAngles/2<<","
                     << m_PrecomputedSctrInParticleLUTDim.iNumDensityLevels        <<")";
        Macros.AddShaderMacro("SRF_SCATTERING_IN_PARTICLE_LUT_DIM", ss.str());
    }
    Macros.AddShaderMacro("BACK_BUFFER_DOWNSCALE_FACTOR", m_CloudAttribs.uiDownscaleFactor);
}

// Method performs tiling of visible particles
void CCloudsController::PerformTiling(SRenderAttribs &RenderAttribs,
                                      ID3D11UnorderedAccessView *pFirstListIndUAV,
                                      const D3DXMATRIX &TilingMatrix,
                                      UINT uiUAVWidth, UINT uiUAVHeight)
{
    ID3D11DeviceContext *pDeviceContext = RenderAttribs.pDeviceContext;
    ID3D11Device *pDevice = RenderAttribs.pDevice;

    if( !m_PerformTilingTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.AddShaderMacro("LIGHT_SPACE_PASS", false);
        //Macros.AddShaderMacro("NUM_SHADOW_CASCADES", m_CloudAttribs.uiNumCascades);
        //Macros.AddShaderMacro("BEST_CASCADE_SEARCH", false);
        Macros.AddShaderMacro("TILING_MODE", true);
        Macros.Finalize();

        m_PerformTilingTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_PerformTilingTech.CreateVGPShadersFromFile(m_strEffectPath, "RenderCloudsVS", "RenderCloudsGS", "PerformTilingPS", Macros);
        m_PerformTilingTech.SetDS( m_pdsDisableDepth );
        m_PerformTilingTech.SetRS( m_prsSolidFillNoCull );
        m_PerformTilingTech.SetBS( m_pbsDefault );
    }

    ID3D11RenderTargetView *ppRTVs[3];
    CComPtr<ID3D11RenderTargetView> pOrigRTV[_countof(ppRTVs)];
    CComPtr<ID3D11DepthStencilView> pOrigDSV;
    pDeviceContext->OMGetRenderTargets(_countof(ppRTVs), ppRTVs, &pOrigDSV);
    for(int i=0; i < _countof(ppRTVs); ++i)
        pOrigRTV[i].Attach(ppRTVs[i]);
    
    D3D11_VIEWPORT OrigViewPort;
    UINT iNumOldViewports = 1;
    pDeviceContext->RSGetViewports(&iNumOldViewports, &OrigViewPort);

    D3D11_VIEWPORT NewViewPort = OrigViewPort;
    NewViewPort.Width  = (float)uiUAVWidth;
    NewViewPort.Height = (float)uiUAVHeight;
    pDeviceContext->RSSetViewports(1, &NewViewPort);

    ID3D11UnorderedAccessView *pUAVs[] = {pFirstListIndUAV, m_pbufParticleListsBuffUAV};
    UINT puiInitialCounts[] = {0,0};
    pDeviceContext->OMSetRenderTargetsAndUnorderedAccessViews(0, nullptr, nullptr, 0, _countof(pUAVs), pUAVs, puiInitialCounts);
    UINT MinusOne[] = {(UINT)-1, 0,0,0};
    pDeviceContext->ClearUnorderedAccessViewUint(pFirstListIndUAV, MinusOne);

    ID3D11ShaderResourceView *pSRVs[] = 
    {
        m_pbufCloudGridSRV,                // StructuredBuffer<SCloudCellAttribs> g_CloudCells : register( t2 );
        m_pbufCloudParticlesSRV,           // StructuredBuffer<SParticleAttribs> g_Particles : register( t3 );
        nullptr,
        nullptr,          
        m_pbufVisibleParticlesFlagsSRV
	};

    pDeviceContext->GSSetShaderResources(2, _countof(pSRVs), pSRVs);
    pDeviceContext->PSSetShaderResources(2, _countof(pSRVs), pSRVs);

    m_CloudAttribs.fTileTexWidth     = (float) uiUAVWidth;
    m_CloudAttribs.fTileTexHeight    = (float) uiUAVHeight;
    D3DXMatrixTranspose(&m_CloudAttribs.mParticleTilingT, &TilingMatrix);
    UpdateConstantBuffer(pDeviceContext, m_pcbGlobalCloudAttribs, &m_CloudAttribs, sizeof(m_CloudAttribs));

    ID3D11Buffer *pCBs[] = {m_pcbGlobalCloudAttribs, RenderAttribs.pcMediaScatteringParams, RenderAttribs.pcbCameraAttribs, RenderAttribs.pcbLightAttribs};
    pDeviceContext->VSSetConstantBuffers(0, _countof(pCBs), pCBs);
    pDeviceContext->GSSetConstantBuffers(0, _countof(pCBs), pCBs);
    pDeviceContext->PSSetConstantBuffers(0, _countof(pCBs), pCBs);

    pDeviceContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
    m_PerformTilingTech.Apply();
    UINT Strides[] = {sizeof(UINT)};
    UINT Offsets[] = {0};
    pDeviceContext->IASetVertexBuffers(0, 1, &m_pbufSerializedVisibleParticles.p, Strides, Offsets);
    pDeviceContext->IASetInputLayout(m_pRenderCloudsInputLayout);
    pDeviceContext->DrawAuto();
    ID3D11Buffer *pDummyBuff = {nullptr};
    pDeviceContext->IASetVertexBuffers(0, 1, &pDummyBuff, Strides, Offsets);

    pDeviceContext->OMSetRenderTargets(_countof(ppRTVs), ppRTVs, pOrigDSV);
    pDeviceContext->RSSetViewports(iNumOldViewports, &OrigViewPort);

    UnbindPSResources(pDeviceContext);
    UnbindVSResources(pDeviceContext);
    UnbindGSResources(pDeviceContext);
}

// Renders all visible particles
void CCloudsController::RenderParticles(SRenderAttribs &RenderAttribs)
{
    ID3D11DeviceContext *pDeviceContext = RenderAttribs.pDeviceContext;
    ID3D11Device *pDevice = RenderAttribs.pDevice;

    bool bLightSpacePass = RenderAttribs.bLightSpacePass;

    auto &RenderCloudsTech = m_RenderCloudsTech[bLightSpacePass ? 1 : 0];
    
    if( !RenderCloudsTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.AddShaderMacro("LIGHT_SPACE_PASS", bLightSpacePass);
        Macros.AddShaderMacro("NUM_SHADOW_CASCADES", m_CloudAttribs.uiNumCascades);
        Macros.AddShaderMacro("BEST_CASCADE_SEARCH", false);
        Macros.AddShaderMacro("TILING_MODE", false);
        Macros.AddShaderMacro("VOLUMETRIC_BLENDING", m_bPSOrderingAvailable && m_CloudAttribs.bVolumetricBlending);
        Macros.Finalize();

        RenderCloudsTech.SetDeviceAndContext(pDevice, pDeviceContext);
        RenderCloudsTech.CreateVGPShadersFromFile(m_strEffectPath, "RenderCloudsVS", "RenderCloudsGS", "RenderCloudsPS", Macros);
        RenderCloudsTech.SetDS( m_pdsDisableDepth /*m_pdsEnableDepth*/ );
        RenderCloudsTech.SetRS( m_prsSolidFillCullFront );
        RenderCloudsTech.SetBS( m_pbsRT0MulRT1MinRT2Over );

        if( !m_pRenderCloudsInputLayout )
        {
            // Create vertex input layout
            const D3D11_INPUT_ELEMENT_DESC layout[] =
            {
                { "PARTICLE_ID",  0, DXGI_FORMAT_R32_UINT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0 }
            };

	        auto pVSByteCode = RenderCloudsTech.GetVSByteCode();
            HRESULT hr;
            V( pDevice->CreateInputLayout( layout, ARRAYSIZE( layout ),
                                            pVSByteCode->GetBufferPointer(),
										    pVSByteCode->GetBufferSize(),
                                            &m_pRenderCloudsInputLayout ) );
        }
    }

    if( RenderAttribs.bTiledMode && !m_RenderCloudsTiledTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.AddShaderMacro("LIGHT_SPACE_PASS", false);
        Macros.AddShaderMacro("NUM_SHADOW_CASCADES", m_CloudAttribs.uiNumCascades);
        Macros.AddShaderMacro("BEST_CASCADE_SEARCH", false);
        Macros.AddShaderMacro("TILE_SIZE", sm_iTileSize);
        Macros.Finalize();

        m_RenderCloudsTiledTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_RenderCloudsTiledTech.CreateVGPShadersFromFile(m_strEffectPath, "ScreenSizeQuadVS", nullptr, "RenderCloudsTiledPS", Macros);
        m_RenderCloudsTiledTech.SetDS( m_pdsDisableDepth );
        m_RenderCloudsTiledTech.SetRS( m_prsSolidFillNoCull );
        m_RenderCloudsTiledTech.SetBS( m_pbsRT0MulRT1MinRT2Over );
    }

    SerializeVisibileParticles(pDevice, pDeviceContext, m_CamSpaceOrder);

    ID3D11ShaderResourceView *pSRVs[] = 
    {
        RenderAttribs.pDepthBufferSRV,
        m_ptex2DCloudDensitySRV,
        m_pbufCloudGridSRV,                // StructuredBuffer<SCloudCellAttribs> g_CloudCells : register( t2 );
        m_pbufCloudParticlesSRV,           // StructuredBuffer<SParticleAttribs> g_Particles : register( t3 );
        m_ptex3DNoiseSRV,                  // Texture3D<float> g_tex3DNoise                  : register(t4);
        nullptr,                            // t5
        m_pbufVisibleParticlesFlagsSRV,     // t6
        m_pbufParticlesLightingSRV,         // t7
        m_ptex2DScrFirstListIndSRV,         // t8
        m_pbufParticleListsBuffSRV,         // t9
        m_ptex3DPrecomputedParticleDensitySRV, // t10
        m_ptex3DSingleSctrInParticleLUT_SRV,   // t11
        m_ptex3DMultipleSctrInParticleLUT_SRV  // t12
	};

    if( RenderAttribs.bTiledMode )
    {
        D3DXMATRIX WorldViewProj;
        D3DXMatrixTranspose(&WorldViewProj, &RenderAttribs.m_pCameraAttribs->WorldViewProjT);
        PerformTiling(RenderAttribs, m_ptex2DScrFirstListIndUAV, WorldViewProj, m_uiTileTexWidth, m_uiTileTexHeight);

        pDeviceContext->VSSetShaderResources(0, _countof(pSRVs), pSRVs);
        pDeviceContext->GSSetShaderResources(0, _countof(pSRVs), pSRVs);
        pDeviceContext->PSSetShaderResources(0, _countof(pSRVs), pSRVs);

        RenderQuad(pDeviceContext, m_RenderCloudsTiledTech, m_uiBackBufferWidth/m_CloudAttribs.uiDownscaleFactor, m_uiBackBufferHeight/m_CloudAttribs.uiDownscaleFactor);
    }
    else
    {
        pDeviceContext->VSSetShaderResources(0, _countof(pSRVs), pSRVs);
        pDeviceContext->GSSetShaderResources(0, _countof(pSRVs), pSRVs);
        pDeviceContext->PSSetShaderResources(0, _countof(pSRVs), pSRVs);

        pDeviceContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_POINTLIST);
        RenderCloudsTech.Apply();
        UINT Strides[] = {sizeof(UINT)};
        UINT Offsets[] = {0};
        pDeviceContext->IASetVertexBuffers(0, 1, &m_pbufSerializedVisibleParticles.p, Strides, Offsets);
        pDeviceContext->IASetInputLayout(m_pRenderCloudsInputLayout);
        pDeviceContext->DrawAuto();
    }

    if( !RenderAttribs.bLightSpacePass && m_bPSOrderingAvailable && m_CloudAttribs.bVolumetricBlending )
    {
        if( !m_ApplyParticleLayersTech.IsValid() )
        {
            CD3DShaderMacroHelper Macros;
            DefineMacros(Macros);
            Macros.Finalize();

            m_ApplyParticleLayersTech.SetDeviceAndContext(pDevice, pDeviceContext);
            m_ApplyParticleLayersTech.CreateVGPShadersFromFile(m_strEffectPath, "ScreenSizeQuadVS", nullptr, "ApplyParticleLayersPS", Macros);
            m_ApplyParticleLayersTech.SetDS( m_pdsDisableDepth );
            m_ApplyParticleLayersTech.SetRS( m_prsSolidFillNoCull );
            m_ApplyParticleLayersTech.SetBS( m_pbsRT0MulRT1MinRT2Over );
        }

        // We need to remove UAVs from the pipeline to be able to bind it as shader resource
        ID3D11RenderTargetView *pRTVs[] = {m_ptex2DScrSpaceCloudTransparencyRTV, m_ptex2DScrSpaceDistToCloudRTV, m_ptex2DScreenCloudColorRTV};
        ID3D11RenderTargetView *pDwnsclRTVs[] = {m_ptex2DDownscaledScrCloudTransparencyRTV, m_ptex2DDownscaledScrDistToCloudRTV, m_ptex2DDownscaledScrCloudColorRTV};
        if(m_CloudAttribs.uiDownscaleFactor > 1 )
            pDeviceContext->OMSetRenderTargets(_countof(pDwnsclRTVs), pDwnsclRTVs, nullptr);
        else
            pDeviceContext->OMSetRenderTargets(_countof(pRTVs), pRTVs, nullptr);

        ID3D11ShaderResourceView *pSRVs[] = 
        {
            m_pbufParticleLayersSRV
        };
        pDeviceContext->PSSetShaderResources(0, _countof(pSRVs), pSRVs);

        RenderQuad(pDeviceContext, m_ApplyParticleLayersTech);
    }

    UnbindPSResources(pDeviceContext);
    UnbindVSResources(pDeviceContext);
    UnbindGSResources(pDeviceContext);
}

// Renders flat clouds on a spherical layer
void CCloudsController::RenderFlatClouds(SRenderAttribs &RenderAttribs)
{
    ID3D11DeviceContext *pDeviceContext = RenderAttribs.pDeviceContext;
    ID3D11Device *pDevice = RenderAttribs.pDevice;
    
    bool bLightSpacePass = RenderAttribs.bLightSpacePass;

    auto &RenderFlatCloudsTech = m_RenderFlatCloudsTech[bLightSpacePass ? 1 : 0];
    if(!RenderFlatCloudsTech.IsValid())
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.AddShaderMacro("LIGHT_SPACE_PASS", bLightSpacePass);
        Macros.AddShaderMacro("NUM_SHADOW_CASCADES", m_CloudAttribs.uiNumCascades);
        Macros.AddShaderMacro("BEST_CASCADE_SEARCH", false);
        Macros.Finalize();

        RenderFlatCloudsTech.SetDeviceAndContext(pDevice, pDeviceContext);
        RenderFlatCloudsTech.CreateVGPShadersFromFile(m_strEffectPath, "ScreenSizeQuadVS", nullptr, "RenderFlatCloudsPS", Macros);
        RenderFlatCloudsTech.SetDS( m_pdsDisableDepth );
        RenderFlatCloudsTech.SetRS( m_prsSolidFillNoCull );
        RenderFlatCloudsTech.SetBS( m_pbsDefault );
    }


    ID3D11ShaderResourceView *pSRVs[] = 
    {
        RenderAttribs.pDepthBufferSRV,
        m_ptex2DCloudDensitySRV,
        nullptr,                                 
        m_ptex2DMaxDensityMipMapSRV,             // Texture2D<float> g_tex2MaxDensityMip           : register(t3);
        m_ptex3DNoiseSRV,                        // Texture3D<float> g_tex3DNoise                  : register(t4);
        RenderAttribs.pPrecomputedNetDensitySRV, // Texture2D<float2> g_tex2DOccludedNetDensityToAtmTop : register( t5 );
        nullptr,                                 // t6
        RenderAttribs.pAmbientSkylightSRV        // t7
	};
    pDeviceContext->VSSetShaderResources(0, _countof(pSRVs), pSRVs);
    pDeviceContext->PSSetShaderResources(0, _countof(pSRVs), pSRVs);
        
    if( !bLightSpacePass && m_CloudAttribs.uiDownscaleFactor > 1 )
    {
        ID3D11ShaderResourceView *pSRVs2[] = 
        {
            m_ptex2DDownscaledScrCloudTransparencySRV,
            m_ptex2DDownscaledScrDistToCloudSRV,
            m_ptex2DDownscaledScrCloudColorSRV
        };
        pDeviceContext->PSSetShaderResources(11, _countof(pSRVs2), pSRVs2);
    }

    RenderQuad(pDeviceContext, RenderFlatCloudsTech);

    UnbindPSResources(pDeviceContext);
    UnbindVSResources(pDeviceContext);
    UnbindGSResources(pDeviceContext);
}

// Renders light space density from light
void CCloudsController::RenderLightSpaceDensity(SRenderAttribs &RenderAttribs)
{
    ID3D11DeviceContext *pDeviceContext = RenderAttribs.pDeviceContext;
    ID3D11Device *pDevice = RenderAttribs.pDevice;

    m_CloudAttribs.fTime = RenderAttribs.fCurrTime;
    m_CloudAttribs.f4Parameter.x = (float)RenderAttribs.iCascadeIndex;
    UpdateConstantBuffer(pDeviceContext, m_pcbGlobalCloudAttribs, &m_CloudAttribs, sizeof(m_CloudAttribs));

    ID3D11Buffer *pCBs[] = {m_pcbGlobalCloudAttribs, RenderAttribs.pcMediaScatteringParams, RenderAttribs.pcbCameraAttribs, RenderAttribs.pcbLightAttribs};
    pDeviceContext->VSSetConstantBuffers(0, _countof(pCBs), pCBs);
    pDeviceContext->GSSetConstantBuffers(0, _countof(pCBs), pCBs);
    pDeviceContext->PSSetConstantBuffers(0, _countof(pCBs), pCBs);
    pDeviceContext->CSSetConstantBuffers(0, _countof(pCBs), pCBs);

    ID3D11SamplerState *pSamplers[] = {m_psamLinearClamp, m_psamLinearWrap, m_psamPointWrap};
    pDeviceContext->VSSetSamplers(0, _countof(pSamplers), pSamplers);
    pDeviceContext->GSSetSamplers(0, _countof(pSamplers), pSamplers);
    pDeviceContext->PSSetSamplers(0, _countof(pSamplers), pSamplers);
    pDeviceContext->CSSetSamplers(0, _countof(pSamplers), pSamplers);


    ID3D11RenderTargetView *ppOrigRTVs[2];
    CComPtr<ID3D11DepthStencilView> pOrigDSV;
    pDeviceContext->OMGetRenderTargets(_countof(ppOrigRTVs), ppOrigRTVs, &pOrigDSV);
    CComPtr<ID3D11RenderTargetView> pTransparencyRTV, pMinMaxDepthRTV;
    pTransparencyRTV.Attach(ppOrigRTVs[0]);
    pMinMaxDepthRTV.Attach(ppOrigRTVs[1]);
    
    m_CloudAttribs.f2LiSpCloudDensityDim.x = static_cast<float>(RenderAttribs.uiLiSpCloudDensityDim);
    m_CloudAttribs.f2LiSpCloudDensityDim.y = static_cast<float>(RenderAttribs.uiLiSpCloudDensityDim);

    float fOneMinusEpsilon = 1.f;
    --((int&)fOneMinusEpsilon);
    const float One[4] = {1, 1, 1, fOneMinusEpsilon};
    pDeviceContext->ClearRenderTargetView(pTransparencyRTV, One);

    const float InitialMinMaxDepth[4] = {0, 0, 0, FLT_MIN};
    pDeviceContext->ClearRenderTargetView(pMinMaxDepthRTV, InitialMinMaxDepth);

    RenderAttribs.bLightSpacePass = true;
    RenderFlatClouds(RenderAttribs);
}

// Merges light space distance to cloud with the shadow map
void CCloudsController::MergeLiSpDensityWithShadowMap(SRenderAttribs &RenderAttribs)
{
    ID3D11DeviceContext *pDeviceContext = RenderAttribs.pDeviceContext;
    ID3D11Device *pDevice = RenderAttribs.pDevice;

    if( !m_RenderCloudDetphToShadowMap.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.Finalize();

        m_RenderCloudDetphToShadowMap.SetDeviceAndContext(pDevice, pDeviceContext);
        m_RenderCloudDetphToShadowMap.CreateVGPShadersFromFile(m_strEffectPath, "ScreenSizeQuadVS", nullptr, "RenderCloudDepthToShadowMapPS", Macros);
        m_RenderCloudDetphToShadowMap.SetDS( m_pdsEnableDepth );
        m_RenderCloudDetphToShadowMap.SetRS( m_prsSolidFillNoCull );
        m_RenderCloudDetphToShadowMap.SetBS( m_pbsDefault );
    }

    CComPtr<ID3D11RenderTargetView> pOrigRTV;
    CComPtr<ID3D11DepthStencilView> pOrigDSV;
    pDeviceContext->OMGetRenderTargets(1, &pOrigRTV, &pOrigDSV);
    
    D3D11_VIEWPORT OrigViewPort;
    UINT iNumOldViewports = 1;
    pDeviceContext->RSGetViewports(&iNumOldViewports, &OrigViewPort);

    pDeviceContext->OMSetRenderTargets(0, nullptr, RenderAttribs.pShadowMapDSV);
    ID3D11ShaderResourceView *pSRVs[] = 
    {
        RenderAttribs.pLiSpCloudTransparencySRV,
        RenderAttribs.pLiSpCloudMinMaxDepthSRV
    };
    pDeviceContext->PSSetShaderResources(0, _countof(pSRVs), pSRVs);

    ID3D11SamplerState *pSamplers[] = {m_psamLinearClamp, m_psamLinearWrap, m_psamPointWrap};
    pDeviceContext->VSSetSamplers(0, _countof(pSamplers), pSamplers);
    pDeviceContext->PSSetSamplers(0, _countof(pSamplers), pSamplers);

    m_CloudAttribs.f4Parameter.x = (float)RenderAttribs.iCascadeIndex;
    UpdateConstantBuffer(pDeviceContext, m_pcbGlobalCloudAttribs, &m_CloudAttribs, sizeof(m_CloudAttribs));

    ID3D11Buffer *pCBs[] = {m_pcbGlobalCloudAttribs};
    pDeviceContext->PSSetConstantBuffers(0, _countof(pCBs), pCBs);

    RenderQuad(pDeviceContext, m_RenderCloudDetphToShadowMap);

    pDeviceContext->OMSetRenderTargets(1, &pOrigRTV.p, pOrigDSV);
    pDeviceContext->RSSetViewports(iNumOldViewports, &OrigViewPort);
}

// Renders cloud color, transparency and distance to clouds from camera
void CCloudsController::RenderScreenSpaceDensityAndColor(SRenderAttribs &RenderAttribs)
{
    ID3D11DeviceContext *pDeviceContext = RenderAttribs.pDeviceContext;
    ID3D11Device *pDevice = RenderAttribs.pDevice;

    if( !m_ptex3DPrecomputedParticleDensitySRV )
    {
        HRESULT hr;
        V( PrecomputParticleDensity(pDevice, pDeviceContext) );
    }

    if( !m_ptex3DSingleSctrInParticleLUT_SRV || !m_ptex3DMultipleSctrInParticleLUT_SRV )
    {
        HRESULT hr;
        V( PrecomputeScatteringInParticle(pDevice, pDeviceContext) );

        //V( ComputeExitance(pDevice, pDeviceContext) );
    }

    PerformTilingInLightSpace(RenderAttribs);

    // Compute visibility so that only visible particles are preocessed
    ComputeParticleVisibility(pDevice, pDeviceContext, RenderAttribs.pcbCameraAttribs, RenderAttribs.bLightSpacePass);

    ComputeParticleLighting(RenderAttribs);

    m_CloudAttribs.fTime = RenderAttribs.fCurrTime;
    m_CloudAttribs.f4Parameter.x = (float)RenderAttribs.iCascadeIndex;
    UpdateConstantBuffer(pDeviceContext, m_pcbGlobalCloudAttribs, &m_CloudAttribs, sizeof(m_CloudAttribs));

    ID3D11Buffer *pCBs[] = {m_pcbGlobalCloudAttribs, RenderAttribs.pcMediaScatteringParams, RenderAttribs.pcbCameraAttribs, RenderAttribs.pcbLightAttribs};
    pDeviceContext->VSSetConstantBuffers(0, _countof(pCBs), pCBs);
    pDeviceContext->GSSetConstantBuffers(0, _countof(pCBs), pCBs);
    pDeviceContext->PSSetConstantBuffers(0, _countof(pCBs), pCBs);
    pDeviceContext->CSSetConstantBuffers(0, _countof(pCBs), pCBs);

    ID3D11SamplerState *pSamplers[] = {m_psamLinearClamp, m_psamLinearWrap, m_psamPointWrap};
    pDeviceContext->VSSetSamplers(0, _countof(pSamplers), pSamplers);
    pDeviceContext->GSSetSamplers(0, _countof(pSamplers), pSamplers);
    pDeviceContext->PSSetSamplers(0, _countof(pSamplers), pSamplers);
    pDeviceContext->CSSetSamplers(0, _countof(pSamplers), pSamplers);

    CComPtr<ID3D11RenderTargetView> pOrigRTV;
    CComPtr<ID3D11DepthStencilView> pOrigDSV;
    pDeviceContext->OMGetRenderTargets(1, &pOrigRTV, &pOrigDSV);

    D3D11_VIEWPORT OrigViewPort;
    UINT iNumOldViewports = 1;
    pDeviceContext->RSGetViewports(&iNumOldViewports, &OrigViewPort);

    float Zero[4]={0,0,0,FLT_MIN};
    pDeviceContext->ClearRenderTargetView(m_ptex2DScreenCloudColorRTV, Zero);

    float fOneMinusEpsilon = 1.f;
    --((int&)fOneMinusEpsilon);
    const float One[4] = {1, 1, 1, fOneMinusEpsilon}; // Use 1-Epsilon to block fast clear path
    pDeviceContext->ClearRenderTargetView(m_ptex2DScrSpaceCloudTransparencyRTV, One);

    if( m_bPSOrderingAvailable && m_CloudAttribs.bVolumetricBlending )
    {
        CComPtr<ID3D11Resource> pDstRes;
        m_pbufParticleLayersUAV->GetResource(&pDstRes);
        pDeviceContext->CopyResource(pDstRes, m_pbufClearParticleLayers);
    }
    // With complimentary depth buffer 0 is the far clipping plane
    // TODO: output some distance from shader (or clear with distanc to horizon?). (Do not forget about sample refinement!)
    const float InitialMinMaxZ[4] = {+FLT_MAX, -FLT_MAX, 0, 0};
    pDeviceContext->ClearRenderTargetView(m_ptex2DScrSpaceDistToCloudRTV, InitialMinMaxZ);

    RenderAttribs.bLightSpacePass = false;
    if(m_CloudAttribs.uiDownscaleFactor > 1 )
    {
        D3D11_VIEWPORT NewViewPort = OrigViewPort;
        NewViewPort.Width  = m_CloudAttribs.fDownscaledBackBufferWidth;
        NewViewPort.Height = m_CloudAttribs.fDownscaledBackBufferHeight;
        pDeviceContext->RSSetViewports(1, &NewViewPort);

        pDeviceContext->ClearRenderTargetView(m_ptex2DDownscaledScrCloudColorRTV, Zero);
        pDeviceContext->ClearRenderTargetView(m_ptex2DDownscaledScrCloudTransparencyRTV, One);
        pDeviceContext->ClearRenderTargetView(m_ptex2DDownscaledScrDistToCloudRTV, InitialMinMaxZ);
        ID3D11RenderTargetView *pRTVs[] = {m_ptex2DDownscaledScrCloudTransparencyRTV, m_ptex2DDownscaledScrDistToCloudRTV, m_ptex2DDownscaledScrCloudColorRTV};
        if( m_bPSOrderingAvailable && m_CloudAttribs.bVolumetricBlending )
        {
            ID3D11UnorderedAccessView *pUAVs[] = {m_pbufParticleLayersUAV};
            UINT puiInitialCounts[_countof(pUAVs)] = {0};
            pDeviceContext->OMSetRenderTargetsAndUnorderedAccessViews(_countof(pRTVs), pRTVs, nullptr, 3, _countof(pUAVs), pUAVs, puiInitialCounts);
        }
        else
        {
            pDeviceContext->OMSetRenderTargets(_countof(pRTVs), pRTVs, nullptr);
        }
        RenderParticles(RenderAttribs);
        
        pDeviceContext->RSSetViewports(1, &OrigViewPort);
    }

    ID3D11RenderTargetView *pRTVs[] = {m_ptex2DScrSpaceCloudTransparencyRTV, m_ptex2DScrSpaceDistToCloudRTV, m_ptex2DScreenCloudColorRTV};
    pDeviceContext->OMSetRenderTargets(_countof(pRTVs), pRTVs, nullptr);

    RenderFlatClouds(RenderAttribs);
    if(m_CloudAttribs.uiDownscaleFactor == 1 )
    {
        if( m_bPSOrderingAvailable && m_CloudAttribs.bVolumetricBlending )
        {
            ID3D11RenderTargetView *pRTVs[] = {m_ptex2DScrSpaceCloudTransparencyRTV, nullptr, m_ptex2DScreenCloudColorRTV};
            ID3D11UnorderedAccessView *pUAVs[] = {m_pbufParticleLayersUAV};
            UINT puiInitialCounts[_countof(pUAVs)] = {0};
            pDeviceContext->OMSetRenderTargetsAndUnorderedAccessViews(_countof(pRTVs), pRTVs, nullptr, 3, _countof(pUAVs), pUAVs, puiInitialCounts);
        }
        RenderParticles(RenderAttribs);
    }

    pDeviceContext->OMSetRenderTargets(1, &pOrigRTV.p, pOrigDSV);
    pDeviceContext->RSSetViewports(iNumOldViewports, &OrigViewPort);
}

// Combines cloud color & transparancy with back buffer
void CCloudsController::CombineWithBackBuffer(ID3D11Device *pDevice, 
                                              ID3D11DeviceContext *pDeviceContext, 
                                              ID3D11ShaderResourceView *pDepthBufferSRV,
                                              ID3D11ShaderResourceView *pBackBufferSRV)
{
    if( !m_CombineWithBBTech.IsValid() )
    {
        CD3DShaderMacroHelper Macros;
        DefineMacros(Macros);
        Macros.Finalize();

        m_CombineWithBBTech.SetDeviceAndContext(pDevice, pDeviceContext);
        m_CombineWithBBTech.CreateVGPShadersFromFile(m_strEffectPath, "ScreenSizeQuadVS", nullptr, "CombineWithBackBufferPS", Macros);
        m_CombineWithBBTech.SetDS( m_pdsDisableDepth );
        m_CombineWithBBTech.SetRS( m_prsSolidFillNoCull );
        m_CombineWithBBTech.SetBS( m_pbsDefault );
    }

    ID3D11ShaderResourceView *pSRVs[] = 
    {
        pDepthBufferSRV,
        pBackBufferSRV
    };
    pDeviceContext->PSSetShaderResources(0, _countof(pSRVs), pSRVs);
    
    ID3D11ShaderResourceView *pSRVs2[] = 
    {
        m_ptex2DScrSpaceCloudTransparencySRV,
        m_ptex2DScrSpaceDistToCloudSRV,
        m_ptex2DScreenCloudColorSRV
    };
    pDeviceContext->PSSetShaderResources(11, _countof(pSRVs2), pSRVs2);

    RenderQuad(pDeviceContext, m_CombineWithBBTech);
}
