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
#include <vtkPointData.h>
#include <vtkPolyData.h>


namespace ospray {
namespace plugin_sensei {


struct OSPRayStudioVisualization::InternalsType {
  bool Execute();

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

bool OSPRayStudioVisualization::Execute(sensei::DataAdaptor *data) {
  vtkDataObject *dataObject;
  if (data->GetMesh("bodies", /*structureOnly=*/false, dataObject)) {
    SENSEI_ERROR("Failed to get mesh 'bodies'")//no semicolon
    return /*success=*/false;
  }

  vtkPolyData *polyData;
  polyData = dynamic_cast<vtkPolyData *>(dataObject);
  if (!polyData) {
    SENSEI_ERROR("Expected mesh 'bodies' to be a vtkPolyData")//no semicolon
    return /*success=*/false;
  }

  if (data->AddArray(dataObject, "bodies", vtkDataObject::POINT, "position")) {
    SENSEI_ERROR("Failed to get array 'position' from mesh 'bodies'")//no semicolon
    return /*success=*/false;
  }

  vtkAbstractArray *abstractArray;
  abstractArray = polyData->GetPointData()->GetAbstractArray("position");

  vtkFloatArray *floatArray;
  if (!(floatArray = vtkFloatArray::SafeDownCast(abstractArray))) {
    SENSEI_ERROR("Expected array 'position' from mesh 'bodies' to be a vtkFloatArray")//no semicolon
    return /*success=*/false;
  }

  size_t nComponents = floatArray->GetNumberOfComponents();
  if (nComponents != 3) {
    SENSEI_ERROR("Expected array 'position' from mesh 'bodies' to be a vtkFloatArray with 3 components")//no semicolon
    return /*success=*/false;
  }

  if (!floatArray->HasStandardMemoryLayout()) {
    SENSEI_ERROR("Expected array 'position' from mesh 'bodies' to be a vtkFloatArray with 3 components and standard memory layout")//no semicolon
    return /*success=*/false;
  }

  // size_t nBodies = floatArray->GetNumberOfTuples();

  // Bodies bodies;
  // bodies.pos = (Point *)floatArray->GetPointer(0);

  return /*success=*/Execute();
}

bool OSPRayStudioVisualization::Execute() {
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
