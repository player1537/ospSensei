// Copyright 2009-2021 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

// plugin_sensei
#include <plugin_sensei/OSPRayStudioVisualization.h>

// stdlib
#include <memory>

// OSPRay Studio
#include "app/ospStudio.h"
#include "app/Plugin.h"

// SENSEI
#include <ConfigurableInTransitDataAdaptor.h>

// VTK
#include <vtkNew.h>

namespace ospray {
namespace plugin_sensei {


struct PluginSENSEI : public Plugin {
  PluginSENSEI()
    : Plugin("SENSEI")
  {}

  void mainMethod(std::shared_ptr<StudioContext> ctx) override;
};

void PluginSENSEI::mainMethod(std::shared_ptr<StudioContext> ctx) {
  using DataAdaptor = sensei::ConfigurableInTransitDataAdaptor;
  using AnalysisAdaptor = plugin_sensei::OSPRayStudioVisualization;

  vtkNew<AnalysisAdaptor> analysisAdaptor;

  if (0 != analysisAdaptor->Initialize()) {
    SENSEI_ERROR("Failed to initialize analysis adaptor")//no semicolon
    return;
  }

  analysisAdaptor->SetRootNode(ctx->frame->child("world").shared_from_this());
  analysisAdaptor->SetScheduler(ctx->scheduler);

  vtkNew<DataAdaptor> dataAdaptor;

  if (0 != dataAdaptor->Initialize("transport.xml")) {
    SENSEI_ERROR("Failed to initialize data adaptor")//no semicolon
    return;
  }

  if (0 != dataAdaptor->OpenStream()) {
    SENSEI_ERROR("Failed to open data adaptor stream")//no semicolon
    return;
  }

  for (;;) {
    if (0 != analysisAdaptor->Execute(dataAdaptor)) {
      SENSEI_ERROR("Failed to execute analysis adaptor")//no semicolon
      return;
    }

    if (0 != dataAdaptor->ReleaseData()) {
      SENSEI_ERROR("Failed to release data adaptor data")//no semicolon
      return;
    }

    if (0 != dataAdaptor->AdvanceStream()) {
      break;
    }
  }

  if (0 != dataAdaptor->CloseStream()) {
    SENSEI_ERROR("Failed to close data adaptor stream")//no semicolon
    return;
  }

  if (0 != dataAdaptor->Finalize()) {
    SENSEI_ERROR("Failed to finalize data adaptor")//no semicolon
    return;
  }

  dataAdaptor.Reset();

  if (0 != analysisAdaptor->Finalize()) {
    SENSEI_ERROR("Failed to finalize analysis adaptor")//no semicolon
    return;
  }

  analysisAdaptor.Reset();
}

extern "C" PLUGIN_INTERFACE Plugin *init_plugin_sensei()
{
  std::cout << "loaded plugin 'SENSEI'!" << std::endl;
  return new PluginSENSEI();
}


} // namespace plugin_sensei
} // namespace ospray
