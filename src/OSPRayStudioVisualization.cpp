/**
 *
 */

// self
#include <plugin_sensei/OSPRayStudioVisualization.h>

// stdlib
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// OSPRay Studio
#include "sg/Node.h"

// SENSEI
#include <DataAdaptor.h>
#include <Error.h>

// VTK
#include <vtkDataObject.h>
#include <vtkFloatArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>


namespace ospray {
namespace plugin_sensei {


struct OSPRayStudioVisualization::InternalsType {
  bool Execute(size_t nPoints, float *positions);

  ospray::sg::NodePtr RootNode{};
  ospray::sg::SchedulerPtr Scheduler{};

  ospray::sg::NodePtr Transform{};
  ospray::sg::NodePtr Sphere{};
};


senseiNewMacro(OSPRayStudioVisualization);


int OSPRayStudioVisualization::Initialize() {
  this->Internals = new InternalsType;
  return /*failure=*/0;
}

void OSPRayStudioVisualization::SetRootNode(const ospray::sg::NodePtr &root) {
  this->Internals->RootNode = root;
}

void OSPRayStudioVisualization::SetScheduler(const ospray::sg::SchedulerPtr &scheduler) {
  this->Internals->Scheduler = scheduler;
}

bool OSPRayStudioVisualization::Execute(sensei::DataAdaptor *data) {
  vtkDataObject *mesh;
  if (data->GetMesh("bodies", /*structureOnly=*/false, mesh)) {
    SENSEI_ERROR("Failed to get mesh 'bodies'")//no semicolon
    return /*success=*/false;
  }

  vtkMultiBlockDataSet *multiBlockDataSet;
  multiBlockDataSet = dynamic_cast<vtkMultiBlockDataSet *>(mesh);
  if (!multiBlockDataSet) {
    SENSEI_ERROR("Expected mesh 'bodies' to be a vtkMultiBlockDataSet")//no semicolon
    return /*success=*/false;
  }

  if (multiBlockDataSet->GetNumberOfBlocks() > 1) {
    SENSEI_ERROR("Expected mesh 'bodies' to have only 1 block")//no semicolon
    return /*success=*/false;
  }

  vtkDataObject *dataObject;
  dataObject = multiBlockDataSet->GetBlock(0);

  vtkPolyData *polyData;
  polyData = dynamic_cast<vtkPolyData *>(dataObject);
  if (!polyData) {
    SENSEI_ERROR("Expected mesh 'bodies' first block to be a vtkPolyData")//no semicolon
    return /*success=*/false;
  }

  if (data->AddArray(dataObject, "bodies", vtkDataObject::POINT, "position")) {
    SENSEI_ERROR("Failed to get array 'position' from mesh 'bodies' first block")//no semicolon
    return /*success=*/false;
  }

  vtkPoints *points;
  points = polyData->GetPoints();
  if (!points) {
    SENSEI_ERROR("Expected points from mesh 'bodies' first block to exist")//no semicolon
    return /*success=*/false;
  }

  vtkDataArray *dataArray;
  dataArray = points->GetData();
  if (!dataArray) {
    SENSEI_ERROR("Expected points from mesh 'bodies' first block to have data")//no semicolon
    return /*success=*/false;
  }

  vtkFloatArray *floatArray;
  if (!(floatArray = vtkFloatArray::SafeDownCast(dataArray))) {
    SENSEI_ERROR("Expected array 'position' from mesh 'bodies' first block to be a vtkFloatArray")//no semicolon
    return /*success=*/false;
  }

  size_t nComponents = floatArray->GetNumberOfComponents();
  if (nComponents != 3) {
    SENSEI_ERROR("Expected array 'position' from mesh 'bodies' first block to be a vtkFloatArray with 3 components")//no semicolon
    return /*success=*/false;
  }

  if (!floatArray->HasStandardMemoryLayout()) {
    SENSEI_ERROR("Expected array 'position' from mesh 'bodies' first block to be a vtkFloatArray with 3 components and standard memory layout")//no semicolon
    return /*success=*/false;
  }

  size_t nPoints = floatArray->GetNumberOfTuples();
  float *positions = floatArray->GetPointer(0);

  return /*success=*/this->Internals->Execute(nPoints, positions);
}

bool OSPRayStudioVisualization::InternalsType::Execute(size_t nPoints, float *positions) {
  std::fprintf(stderr, "positions[0] = %+0.2f\n", positions[0]);
  std::fprintf(stderr, "positions[1] = %+0.2f\n", positions[1]);
  std::fprintf(stderr, "positions[2] = %+0.2f\n", positions[2]);

  // OSPData sharedData;
  // sharedData = ospNewSharedData((void *)bodies.pos, OSP_FLOAT, bodies.count);
  // ospCommit(sharedData);

  // OSPData data;
  // data = ospNewData(OSP_VEC3F, bodies.count);
  // ospCopyData(sharedData, data);
  // ospCommit(data);
  // ospRelease(sharedData);

  // OSPGeometry sphere;
  // sphere = ospNewGeometry("sphere");
  // ospSetObject(sphere, "sphere.position", data);
  // ospCommit(sphere);

  return /*success=*/true;
}

int OSPRayStudioVisualization::Finalize() {
  delete this->Internals;

  return 0; // no error
}


} /* namespace plugin_sensei */
} /* namespace ospray */
// vim: ts=2:sts=2
