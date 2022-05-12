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
  float GeometryRadius{0.1f};

  int FrameNumber{0};
  int CommRank{-1};
  int CommSize{-1};

  OSPDevice Device;
  OSPData GeometrySpherePositionData;
  OSPGeometry Geometry;
  OSPMaterial Material;
  OSPGeometricModel GeometricModel;
  OSPGroup Group;
  OSPInstance Instance;
  OSPLight Light;
  std::vector<float> WorldRegion;
  OSPData WorldRegionData;
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
  MPI_Comm comm;
  comm = GetCommunicator();

  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

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

  if (multiBlockDataSet->GetNumberOfBlocks() != size) {
    SENSEI_ERROR("Expected mesh 'bodies' to have only " << size << " blocks")//no semicolon
    return /*success=*/false;
  }

  vtkDataObject *dataObject;
  dataObject = multiBlockDataSet->GetBlock(rank);
  if (dataObject == nullptr) {
    SENSEI_ERROR("Expected mesh 'bodies' block " << rank << " to be non-null")//no semicolon
    return /*success=*/false;
  }

  vtkPolyData *polyData;
  polyData = dynamic_cast<vtkPolyData *>(dataObject);
  if (!polyData) {
    SENSEI_ERROR("Expected mesh 'bodies' block " << rank << " to be a vtkPolyData")//no semicolon
    return /*success=*/false;
  }

  if (data->AddArray(mesh, "bodies", vtkDataObject::POINT, "position")) {
    SENSEI_ERROR("Failed to get array 'position' from mesh 'bodies' block " << rank)//no semicolon
    return /*success=*/false;
  }

  vtkPoints *points;
  points = polyData->GetPoints();
  if (!points) {
    SENSEI_ERROR("Expected points from mesh 'bodies' block " << rank << " to exist")//no semicolon
    return /*success=*/false;
  }

  vtkDataArray *dataArray;
  dataArray = points->GetData();
  if (!dataArray) {
    SENSEI_ERROR("Expected points from mesh 'bodies' block " << rank << " to have data")//no semicolon
    return /*success=*/false;
  }

  vtkFloatArray *floatArray;
  if (!(floatArray = vtkFloatArray::SafeDownCast(dataArray))) {
    SENSEI_ERROR("Expected array 'position' from mesh 'bodies' block " << rank << " to be a vtkFloatArray")//no semicolon
    return /*success=*/false;
  }

  size_t nComponents = floatArray->GetNumberOfComponents();
  if (nComponents != 3) {
    SENSEI_ERROR("Expected array 'position' from mesh 'bodies' block " << rank << " to be a vtkFloatArray with 3 components")//no semicolon
    return /*success=*/false;
  }

  if (!floatArray->HasStandardMemoryLayout()) {
    SENSEI_ERROR("Expected array 'position' from mesh 'bodies' block " << rank << " to be a vtkFloatArray with 3 components and standard memory layout")//no semicolon
    return /*success=*/false;
  }

  size_t nPoints = floatArray->GetNumberOfTuples();
  float *positions = floatArray->GetPointer(0);

  return /*success=*/this->Internals->Execute(comm, nPoints, positions);
}

bool OSPRayVisualization::InternalsType::Execute(MPI_Comm comm, size_t nPoints, float *positions) {
  float lo[3], hi[3];

  std::fprintf(stderr, "nPoints = %lu\n", nPoints);
  if (nPoints > 0) {
    std::fprintf(stderr, "positions[0] = %+0.2f\n", positions[0]);
    std::fprintf(stderr, "positions[1] = %+0.2f\n", positions[1]);
    std::fprintf(stderr, "positions[2] = %+0.2f\n", positions[2]);

    lo[0] = hi[0] = positions[0*3+0];
    lo[1] = hi[1] = positions[0*3+1];
    lo[2] = hi[2] = positions[0*3+2];
    for (size_t i=1; i<nPoints; ++i) {
      if (positions[i*3+0] < lo[0]) lo[0] = positions[i*3+0];
      if (positions[i*3+1] < lo[1]) lo[1] = positions[i*3+1];
      if (positions[i*3+2] < lo[2]) lo[2] = positions[i*3+2];
      if (positions[i*3+0] > hi[0]) hi[0] = positions[i*3+0];
      if (positions[i*3+1] > hi[1]) hi[1] = positions[i*3+1];
      if (positions[i*3+2] > hi[2]) hi[2] = positions[i*3+2];
    }

    std::fprintf(stderr, "min = %+0.2f, %+0.2f, %+0.2f\n", lo[0], lo[1], lo[2]);
    std::fprintf(stderr, "max = %+0.2f, %+0.2f, %+0.2f\n", hi[0], hi[1], hi[2]);
  }
  
  if (!HasInitialized) {
    MPI_Comm_rank(comm, &CommRank);
    MPI_Comm_size(comm, &CommSize);

    ospLoadModule("mpi");

    Device = ospNewDevice("mpiDistributed");
    ospDeviceCommit(Device);
    ospSetCurrentDevice(Device);

    GeometrySpherePositionData = nullptr;

    Geometry = ospNewGeometry("sphere");
    ospSetObject(Geometry, "sphere.position", GeometrySpherePositionData);
    ospSetFloat(Geometry, "radius", GeometryRadius);
    ospCommit(Geometry);

    Material = ospNewMaterial(nullptr, "obj");
    ospSetVec3f(Material, "kd", (CommRank % 4 == 0 ? 1.0f : 0.0f), (CommRank % 4 == 1 ? 1.0f : 0.0f), (CommRank % 4 == 2 ? 1.0f : 0.0f));
    ospCommit(Material);

    GeometricModel = ospNewGeometricModel(nullptr);
    ospSetObject(GeometricModel, "geometry", Geometry);
    ospSetObject(GeometricModel, "material", Material);
    ospCommit(GeometricModel);

    Group = ospNewGroup();
    ospSetObjectAsData(Group, "geometry", OSP_GEOMETRIC_MODEL, GeometricModel);
    ospCommit(Group);

    Instance = ospNewInstance(nullptr);
    ospSetObject(Instance, "group", Group);
    ospCommit(Instance);

    Light = ospNewLight("ambient");
    ospCommit(Light);

    WorldRegionData = nullptr;

    World = ospNewWorld();
    ospSetObjectAsData(World, "instance", OSP_INSTANCE, Instance);
    ospSetObjectAsData(World, "light", OSP_LIGHT, Light);
    ospSetObject(World, "region", WorldRegionData);
    ospCommit(World);

    Camera = ospNewCamera("perspective");
    ospSetFloat(Camera, "aspect", (float)Width / (float)Height);
    ospSetVec3f(Camera, "position", 5.0f, 0.0f, 0.0f);
    ospSetVec3f(Camera, "direction", -1.0f, 0.0f, 0.0f);
    ospSetVec3f(Camera, "up", 0.0f, 1.0f, 0.0f);
    ospCommit(Camera);

    Renderer = ospNewRenderer("mpiRaycast");
    ospSetVec3f(Renderer, "backgroundColor", 0.5f, 0.5f, 0.5f);
    ospCommit(Renderer);

    FrameBuffer = ospNewFrameBuffer(Width, Height, OSP_FB_SRGBA, OSP_FB_COLOR | OSP_FB_ACCUM);

    Future = nullptr;

    HasInitialized = true;
  }

  if (GeometrySpherePositionData) {
    ospRelease(GeometrySpherePositionData);
    GeometrySpherePositionData = nullptr;
  }

  GeometrySpherePositionData = ospNewSharedData(positions, OSP_VEC3F, nPoints);
  ospCommit(GeometrySpherePositionData);
  
  ospSetObject(Geometry, "sphere.position", GeometrySpherePositionData);
  ospCommit(Geometry);

  ospSetObject(GeometricModel, "geometry", Geometry);
  ospCommit(GeometricModel);

  ospSetObjectAsData(Group, "geometry", OSP_GEOMETRIC_MODEL, GeometricModel);
  ospCommit(Group);

  ospSetObject(Instance, "group", Group);
  ospCommit(Instance);

  WorldRegion.clear();
  WorldRegionData = nullptr;

  if (nPoints > 0) {
    WorldRegion.insert(WorldRegion.end(), {
      lo[0]-GeometryRadius, lo[1]-GeometryRadius, lo[2]-GeometryRadius,
      hi[0]+GeometryRadius, hi[1]+GeometryRadius, hi[2]+GeometryRadius,
    });

    WorldRegionData = ospNewSharedData(WorldRegion.data(), OSP_BOX3F, WorldRegion.size() / 6);
    ospCommit(WorldRegionData);
  }

  ospSetObjectAsData(World, "instance", OSP_INSTANCE, Instance);
  ospSetObject(World, "region", WorldRegionData);
  ospCommit(World);

  ospResetAccumulation(FrameBuffer);
  Future = ospRenderFrame(FrameBuffer, Renderer, Camera, World);
  ospWait(Future, OSP_TASK_FINISHED);
  ospRelease(Future);
  Future = nullptr;

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
  if (GeometrySpherePositionData) ospRelease(GeometrySpherePositionData);
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
