// Copyright 2014 Mark Harris (https://github.com/harrism)
// SPDX-License-Identifier: Apache-2.0
//
// Modifications (2022) by Tanner Hobson (https://github.com/player1537)

// stdlib
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// SENSEI
#include <ConfigurableAnalysis.h>
#include <VTKDataAdaptor.h>

// VTK
#include <vtkDataObject.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>

#define SOFTENING 1e-9f


typedef struct {
  float x, y, z;
} Point;

typedef struct {
  Point *pos;
  Point *vel;
} Bodies;

void randomize(Point *p, size_t n) {
  for (size_t i = 0; i < n; i++) {
    p[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    p[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    p[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
}

void bodyForce(const Bodies &bodies, float dt, size_t n) {
  #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < n; i++) { 
    float Fx = 0.0f;
    float Fy = 0.0f;
    float Fz = 0.0f;

    for (size_t j = 0; j < n; j++) {
      float dx = bodies.pos[j].x - bodies.pos[i].x;
      float dy = bodies.pos[j].y - bodies.pos[i].y;
      float dz = bodies.pos[j].z - bodies.pos[i].z;
      float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * invDist3;
      Fy += dy * invDist3;
      Fz += dz * invDist3;
    }

    bodies.vel[i].x += dt * Fx;
    bodies.vel[i].y += dt * Fy;
    bodies.vel[i].z += dt * Fz;
  }
}


class OSPRayStudioPointsVisualization : public AnalysisAdapter {
public:
  bool Execute(DataAdaptor *data) override;
  int Finalize() override;

private:
};

bool OSPRayStudioPointsVisualization::Execute(DataAdaptor *data) {
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

  return /*success=*/true;
}

int OSPRayStudioPointsVisualization::Finalize() {
  return 0; // no error
}


int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  size_t nBodies = 30000;
  size_t nIters = 10;
  float dt = 0.01f; // time step
  std::string configFilename{"nbody.xml"};
  if (argc > 1 && argv[1][0] != '\0') nBodies = strtoul(argv[1], NULL, 10);
  if (argc > 2 && argv[2][0] != '\0') nIters = strtoul(argv[2], NULL, 10);
  if (argc > 3 && argv[3][0] != '\0') dt = atof(argv[3]);
  if (argc > 4 && argv[4][0] != '\0') configFilename = argv[4];

  Bodies bodies;
  bodies.pos = (Point *)malloc(nBodies * sizeof(*bodies.pos));
  bodies.vel = (Point *)malloc(nBodies * sizeof(*bodies.vel));

  randomize(bodies.pos, nBodies);
  randomize(bodies.vel, nBodies);

  vtkNew<sensei::ConfigurableAnalysis> analysisAdaptor;
  analysisAdaptor->SetCommunicator(MPI_COMM_WORLD);
  analysisAdaptor->Initialize(configFilename);

  for (int iter = 0; iter <= nIters; ++iter) {
    if (iter > 0) {
      bodyForce(bodies, dt, nBodies); // compute interbody forces

      for (int i = 0 ; i < nBodies; i++) { // integrate position
        bodies.pos[i].x += bodies.vel[i].x * dt;
        bodies.pos[i].y += bodies.vel[i].y * dt;
        bodies.pos[i].z += bodies.vel[i].z * dt;
      }
    }

    vtkNew<vtkFloatArray> floatArrayPos;
    floatArrayPos->Initialize();
    floatArrayPos->SetName("position");
    floatArrayPos->SetNumberOfComponents(3);
    floatArrayPos->SetArray((float *)bodies.pos, 3 * nBodies, /*save=*/0);

    vtkNew<vtkFloatArray> floatArrayVel;
    floatArrayVel->Initialize();
    floatArrayVel->SetName("velocity");
    floatArrayVel->SetNumberOfComponents(3);
    floatArrayVel->SetArray((float *)bodies.vel, 3 * nBodies, /*save=*/0);

    vtkNew<vtkPolyData> polyData;
    polyData->Initialize();
    polyData->GetPointData()->AddArray(floatArrayPos);
    polyData->GetPointData()->AddArray(floatArrayVel);

    vtkNew<sensei::VTKDataAdaptor> vtkDataAdaptor;
    vtkDataAdaptor->SetDataTime(iter * dt);
    vtkDataAdaptor->SetDataTimeStep(iter);
    vtkDataAdaptor->SetDataObject("bodies", polyData);

    analysisAdaptor->Execute(vtkDataAdaptor);

    vtkDataAdaptor->ReleaseData();
  }

  analysisAdaptor->Finalize();

  free(bodies.vel);
  free(bodies.pos);

  MPI_Finalize();

  return 0;
}
