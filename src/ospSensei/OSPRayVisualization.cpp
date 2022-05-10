/**
 *
 */

// self
#include <ospSensei/OSPRayVisualization.h>

// stdlib
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// OSPRay
#include <ospray/ospray.h>
#include <ospray/ospray_util.h>

// SENSEI
#include <DataAdaptor.h>
#include <Error.h>

// VTK
#include <vtkDataObject.h>
#include <vtkFloatArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>


namespace ospSensei {

// helper function to write the rendered image as PPM file
static void writePPM(const char *fileName, int size_x, int size_y, const uint32_t *pixel) {
  using namespace std;

  FILE *file = fopen(fileName, "wb");
  if (!file) {
    fprintf(stderr, "fopen('%s', 'wb') failed: %d", fileName, errno);
    return;
  }
  fprintf(file, "P6\n%i %i\n255\n", size_x, size_y);
  unsigned char *out = (unsigned char *)alloca(3 * size_x);
  for (int y = 0; y < size_y; y++) {
    const unsigned char *in =
        (const unsigned char *)&pixel[(size_y - 1 - y) * size_x];
    for (int x = 0; x < size_x; x++) {
      out[3 * x + 0] = in[4 * x + 0];
      out[3 * x + 1] = in[4 * x + 1];
      out[3 * x + 2] = in[4 * x + 2];
    }
    fwrite(out, 3 * size_x, sizeof(char), file);
  }
  fprintf(file, "\n");
  fclose(file);
}


struct OSPRayVisualization::InternalsType {
  bool Execute(MPI_Comm comm, size_t nPoints, float *positions);
  void Finalize();

  int Width{256};
  int Height{256};
  int FrameNumber{0};

  int CommRank{-1};
  int CommSize{-1};

  OSPDevice Device;
  OSPData GeometrySpherePosition;
  OSPGeometry Geometry;
  OSPGeometricModel GeometricModel;
  OSPGroup Group;
  OSPInstance Instance;
  OSPLight Light;
  std::vector<float> WorldRegionValues;
  OSPData WorldRegion;
  OSPWorld World;
  OSPCamera Camera;
  OSPRenderer Renderer;
  OSPFrameBuffer FrameBuffer;
  OSPFuture Future;
  bool HasInitialized{false};
};


senseiNewMacro(OSPRayVisualization);


int OSPRayVisualization::Initialize() {
  this->Internals = new InternalsType;
  return /*failure=*/0;
}

bool OSPRayVisualization::Execute(sensei::DataAdaptor *data) {
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

  if (data->AddArray(mesh, "bodies", vtkDataObject::POINT, "position")) {
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

  return /*success=*/this->Internals->Execute(GetCommunicator(), nPoints, positions);
}

bool OSPRayVisualization::InternalsType::Execute(MPI_Comm comm, size_t nPoints, float *positions) {
  std::fprintf(stderr, "positions[0] = %+0.2f\n", positions[0]);
  std::fprintf(stderr, "positions[1] = %+0.2f\n", positions[1]);
  std::fprintf(stderr, "positions[2] = %+0.2f\n", positions[2]);
  
  if (!HasInitialized) {
    MPI_Comm_rank(comm, &CommRank);
    MPI_Comm_size(comm, &CommSize);

    ospInit(NULL, NULL);

    // ospLoadModule("mpi");

    // Device = ospNewDevice("mpiDistributed");
    // ospDeviceCommit(Device);
    // ospSetCurrentDevice(Device);

    GeometrySpherePosition = nullptr;

    Geometry = ospNewGeometry("sphere");
    ospSetObject(Geometry, "sphere.position", GeometrySpherePosition);
    ospCommit(Geometry);

    GeometricModel = ospNewGeometricModel(NULL);
    ospSetObject(GeometricModel, "geometry", Geometry);
    ospCommit(GeometricModel);

    Group = ospNewGroup();
    ospSetObjectAsData(Group, "geometry", OSP_GEOMETRY, GeometricModel);
    ospCommit(Group);

    Instance = ospNewInstance(NULL);
    ospSetObject(Instance, "group", Group);
    ospCommit(Instance);

    Light = ospNewLight("ambient");
    ospCommit(Light);

    WorldRegionValues.insert(WorldRegionValues.end(), {
      -100.0f, -100.0f, -100.0f,
      +100.0f, +100.0f, +100.0f,
    });

    WorldRegion = ospNewSharedData(WorldRegionValues.data(), OSP_BOX3F, WorldRegionValues.size() / 6);
    ospCommit(WorldRegion);

    World = ospNewWorld();
    ospSetObjectAsData(World, "instance", OSP_INSTANCE, Instance);
    ospSetObjectAsData(World, "light", OSP_LIGHT, Light);
    ospSetObject(World, "region", WorldRegion);
    ospCommit(World);

    Camera = ospNewCamera("perspective");
    ospSetFloat(Camera, "aspect", (float)Width / (float)Height);
    ospSetVec3f(Camera, "position", 100.0f, 0.0f, 0.0f);
    ospSetVec3f(Camera, "direction", -1.0f, 0.0f, 0.0f);
    ospSetVec3f(Camera, "up", 0.0f, 1.0f, 0.0f);
    ospCommit(Camera);

    Renderer = ospNewRenderer("scivis");
    ospSetVec3f(Renderer, "backgroundColor", 0.5f, 0.5f, 0.5f);
    ospCommit(Renderer);

    FrameBuffer = ospNewFrameBuffer(Width, Height, OSP_FB_RGBA8, OSP_FB_COLOR);

    Future = nullptr;

    HasInitialized = true;
  }

  if (GeometrySpherePosition) {
    ospRelease(GeometrySpherePosition);
    GeometrySpherePosition = nullptr;
  }

  GeometrySpherePosition = ospNewSharedData(positions, OSP_VEC3F, nPoints);
  ospCommit(GeometrySpherePosition);
  
  ospSetObject(Geometry, "sphere.position", GeometrySpherePosition);
  ospCommit(Geometry);

  ospSetObject(GeometricModel, "geometry", Geometry);
  ospCommit(GeometricModel);

  ospSetObjectAsData(Group, "geometry", OSP_GEOMETRY, GeometricModel);
  ospCommit(Group);

  ospSetObject(Instance, "group", Group);
  ospCommit(Instance);

  ospSetObjectAsData(World, "instance", OSP_INSTANCE, Instance);
  ospCommit(World);

  ospResetAccumulation(FrameBuffer);
  Future = ospRenderFrame(FrameBuffer, Renderer, Camera, World);
  ospWait(Future, OSP_TASK_FINISHED);

  if (CommRank == 0) {
    std::string filename = std::string("ospSensei.") + std::to_string(FrameNumber) + std::string(".ppm");
    const void *fb = ospMapFrameBuffer(FrameBuffer, OSP_FB_COLOR);
    writePPM(filename.c_str(), Width, Height, static_cast<const uint32_t *>(fb));
    ospUnmapFrameBuffer(fb, FrameBuffer);
  }

  ++FrameNumber;

  return /*success=*/true;
}

void OSPRayVisualization::InternalsType::Finalize() {
  ospRelease(FrameBuffer);
  ospRelease(Renderer);
  ospRelease(Camera);
  ospRelease(World);
  ospRelease(Light);
  ospRelease(Instance);
  ospRelease(Group);
  ospRelease(GeometricModel);
  ospRelease(Geometry);
  if (GeometrySpherePosition) ospRelease(GeometrySpherePosition);
  ospDeviceRelease(Device);

  ospShutdown();
}

int OSPRayVisualization::Finalize() {
  Internals->Finalize();
  delete this->Internals;

  return 0; // no error
}


} /* namespace ospSensei */
// vim: ts=2:sts=2
