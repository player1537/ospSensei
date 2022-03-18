// Copyright 2009-2021 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

// plugin_sensei
#include <plugin_sensei/OSPRayStudioVisualization.h>

// stdlib
#include <memory>

// OSPRay Studio
#include "app/Batch.h"
#include "app/ospStudio.h"
#include "app/Plugin.h"

// SENSEI
#include <ConfigurableInTransitDataAdaptor.h>

// VTK
#include <vtkNew.h>
#include <vtkSmartPointer.h>

namespace ospray {
namespace plugin_sensei {


struct PluginSENSEI : public Plugin {
  PluginSENSEI()
    : Plugin("SENSEI")
  {}

  void mainMethod(std::shared_ptr<StudioContext> ctx) override;

private:
  std::string optConfigFilename;

  enum class ShouldContinue {
    YES,
    NO,
  };

  bool first{true};
  void onFirstCall(const std::shared_ptr<StudioContext> &ctx);
  ShouldContinue onSubsequentCalls(const std::shared_ptr<StudioContext> &ctx);
  void onLastCall(const std::shared_ptr<StudioContext> &ctx);

  using DataAdaptor = sensei::ConfigurableInTransitDataAdaptor;
  using AnalysisAdaptor = plugin_sensei::OSPRayStudioVisualization;

  vtkSmartPointer<AnalysisAdaptor> analysisAdaptor;
  vtkSmartPointer<DataAdaptor> dataAdaptor;
};

void PluginSENSEI::mainMethod(std::shared_ptr<StudioContext> ctx) {
  if (first) {
    onFirstCall(ctx);
    first = false;
  }


  ShouldContinue sc = onSubsequentCalls(ctx);
  switch (sc) {
    case ShouldContinue::YES: {
      std::fprintf(stderr, "plugin_sensei: We should continue\n");
      auto batch = std::dynamic_pointer_cast<BatchContext>(ctx);
      batch->shouldContinueRendering = true;
      break;
    }

    case ShouldContinue::NO: {
      std::fprintf(stderr, "plugin_sensei: We should not continue\n");
      onLastCall(ctx);
      break;
    }
  }
}

void PluginSENSEI::onFirstCall(const std::shared_ptr<StudioContext> &ctx) {
  std::fprintf(stderr, "plugin_sensei: First call\n");

  auto &studioCommon = ctx->studioCommon;
  int ac = studioCommon.plugin_argc;
  const char **av = studioCommon.plugin_argv;

  for (int i=0; i<ac; ++i) {
    std::string arg = av[i];
    if (arg == "--plugin:sensei:configFilename") {
      optConfigFilename = av[i + 1];
      ++i;
    }
  }

  int mpiIsInitialized;
  MPI_Initialized(&mpiIsInitialized);
  
  if (!mpiIsInitialized) {
    int argc = 0;
    char *argv[] = { NULL };
    int provided;
    int success = MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);
    if (success != MPI_SUCCESS) {
      SENSEI_ERROR("Error while initializing MPI")//no semicolon
      return;
    }

    if (provided != MPI_THREAD_FUNNELED) {
      SENSEI_ERROR("MPI provided the wrong level of thread support")//no semicolon
      return;
    }
  }

  analysisAdaptor = vtkNew<AnalysisAdaptor>();

  if (0 != analysisAdaptor->Initialize()) {
    SENSEI_ERROR("Failed to initialize analysis adaptor")//no semicolon
    return;
  }

  analysisAdaptor->SetRootNode(ctx->frame->child("world").shared_from_this());
  analysisAdaptor->SetScheduler(ctx->scheduler);

  dataAdaptor = vtkNew<DataAdaptor>();

  if (0 != dataAdaptor->Initialize(optConfigFilename)) {
    SENSEI_ERROR("Failed to initialize data adaptor")//no semicolon
    return;
  }

  if (0 != dataAdaptor->OpenStream()) {
    SENSEI_ERROR("Failed to open data adaptor stream")//no semicolon
    return;
  }
}

PluginSENSEI::ShouldContinue PluginSENSEI::onSubsequentCalls(const std::shared_ptr<StudioContext> &ctx) {
  std::fprintf(stderr, "plugin_sensei: Subsequent call\n");

  if (!analysisAdaptor->Execute(dataAdaptor)) {
    SENSEI_ERROR("Failed to execute analysis adaptor")//no semicolon
    return ShouldContinue::NO;
  }

  if (0 != dataAdaptor->ReleaseData()) {
    SENSEI_ERROR("Failed to release data adaptor data")//no semicolon
    return ShouldContinue::NO;
  }

  if (0 != dataAdaptor->AdvanceStream()) {
    return ShouldContinue::NO;
  }

  return ShouldContinue::YES;
}

void PluginSENSEI::onLastCall(const std::shared_ptr<StudioContext> &ctx) {
  std::fprintf(stderr, "plugin_sensei: Last call\n");

  if (0 != dataAdaptor->CloseStream()) {
    SENSEI_ERROR("Failed to close data adaptor stream")//no semicolon
    return;
  }

  if (0 != dataAdaptor->Finalize()) {
    SENSEI_ERROR("Failed to finalize data adaptor")//no semicolon
    return;
  }

  dataAdaptor = nullptr;

  if (0 != analysisAdaptor->Finalize()) {
    SENSEI_ERROR("Failed to finalize analysis adaptor")//no semicolon
    return;
  }

  analysisAdaptor = nullptr;
}

extern "C" PLUGIN_INTERFACE Plugin *init_plugin_sensei()
{
  std::cout << "loaded plugin 'SENSEI'!" << std::endl;
  return new PluginSENSEI();
}


} // namespace plugin_sensei
} // namespace ospray
