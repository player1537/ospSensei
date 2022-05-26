/**
 *
 */

// stdlib
#include <array>
#include <complex>
#include <cstdint>

// vtk
#include <vtkCellData.h>
#include <vtkDataObject.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkUniformGrid.h>
#include <vtkUniformGridAMR.h>
#include <vtkUnsignedShortArray.h>

// sensei
#include <VTKDataAdaptor.h>

// ospSensei
#include <ospSensei/OSPRayVisualization.h>


struct Mandelbrot {
  using ScalarF = float;
  using ScalarU = uint16_t;
  using BoundsF = std::array<ScalarF, 6>;
  enum Bounds { MinX = 0, MinY, MinZ, MaxX, MaxY, MaxZ };
  using ComplexF = std::complex<ScalarF>;

  enum Debug { OnlyData, OnlyNsteps };

  Mandelbrot(size_t nx_, size_t ny_, size_t nz_, BoundsF bounds_);
  ~Mandelbrot();

  void debug(Debug);
  void step(size_t dt);
  vtkUniformGrid *vtk();

  size_t nx, ny, nz;
  BoundsF bounds;
  ScalarF *data;
  ScalarU *nsteps;
};

Mandelbrot::Mandelbrot(size_t nx_, size_t ny_, size_t nz_, Mandelbrot::BoundsF bounds_)
  : nx(nx_)
  , ny(ny_)
  , nz(nz_)
  , bounds(bounds_)
  , data(new ScalarF[2*nx_*ny_*nz_])
  , nsteps(new ScalarU[nx_*ny_*nz_])
{
  for (size_t zi=0; zi<nz; ++zi) {
    size_t zindex = zi*ny*nx;

    for (size_t yi=0; yi<ny; ++yi) {
      size_t yindex = zindex + yi*nx;

      for (size_t xi=0; xi<nx; ++xi) {
        size_t xindex = yindex + xi;

        data[2*xindex+0] = 0.0f;
        data[2*xindex+1] = 0.0f;
        nsteps[xindex] = 0;
      }
    }
  }
}

Mandelbrot::~Mandelbrot() {
  delete[] nsteps;
  delete[] data;
}

void Mandelbrot::debug(Debug which) {
  for (size_t zi=0; zi<nz; ++zi) {
    size_t zindex = zi*ny*nx;

    std::fprintf(stderr, "[");

    for (size_t yi=0; yi<ny; ++yi) {
      size_t yindex = zindex + yi*nx;

      if (yi == 0) std::fprintf(stderr, " [");
      else std::fprintf(stderr, "  [");

      for (size_t xi=0; xi<nx; ++xi) {
        size_t xindex = yindex + xi;

        if (which == OnlyData) {
          std::fprintf(stderr, " %+0.2f%+0.2fi", data[2*xindex+0], data[2*xindex+1]);
        } else if (which == OnlyNsteps) {
          std::fprintf(stderr, " %03d", nsteps[xindex]);
        }
      }
      
      std::fprintf(stderr, "\n");
    }

    std::fprintf(stderr, "\n");
  }
}

void Mandelbrot::step(size_t dt) {
  for (size_t zi=0; zi<nz; ++zi) {
    ScalarF zratio = (ScalarF)zi / (ScalarF)nz;
    ScalarF z = std::get<MinZ>(bounds) * zratio + std::get<MaxZ>(bounds) * (1.0f - zratio);
    size_t zindex = zi*ny*nx;

    for (size_t yi=0; yi<ny; ++yi) {
      ScalarF yratio = (ScalarF)yi / (ScalarF)ny;
      ScalarF y = std::get<MinY>(bounds) * yratio + std::get<MaxY>(bounds) * (1.0f - yratio);
      size_t yindex = zindex + yi*nx;

      for (size_t xi=0; xi<nx; ++xi) {
        ScalarF xratio = (ScalarF)xi / (ScalarF)nx;
        ScalarF x = std::get<MinX>(bounds) * xratio + std::get<MaxX>(bounds) * (1.0f - xratio);
        size_t xindex = yindex + xi;

        for (size_t ti=0; ti<dt; ++ti) {
          ScalarF xd = data[2*xindex+0];
          ScalarF yd = data[2*xindex+1];

          if (xd*xd + yd*yd >= 2.0) {
            break;
          }

          ComplexF temp = std::pow(ComplexF(xd, yd), z);
          data[2*xindex+0] = temp.real() + x;
          data[2*xindex+1] = temp.imag() + y;
          ++nsteps[xindex];
        }
      }
    }
  }
}

vtkUniformGrid *Mandelbrot::vtk() {
  vtkNew<vtkUnsignedShortArray> unsignedShortArray;
  unsignedShortArray->SetName("nsteps");
  unsignedShortArray->SetArray(nsteps, nx*ny*nz, 1);

  vtkUniformGrid *uniformGrid = vtkUniformGrid::New();
  uniformGrid->SetOrigin(
    std::get<MinX>(bounds),
    std::get<MinY>(bounds),
    std::get<MinZ>(bounds));
  uniformGrid->SetSpacing(
    (std::get<MaxX>(bounds) - std::get<MinX>(bounds)) / (float)nx,
    (std::get<MaxY>(bounds) - std::get<MinY>(bounds)) / (float)ny,
    (std::get<MaxZ>(bounds) - std::get<MinZ>(bounds)) / (float)nz);
  uniformGrid->SetDimensions(nx, ny, nz);
  uniformGrid->GetCellData()->SetScalars(unsignedShortArray);
  uniformGrid->GetCellData()->SetActiveScalars("nsteps");

  return uniformGrid;
}


int main(int argc, char **argv) {
  int provided;
  int success = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  if (success != MPI_SUCCESS) {
    SENSEI_ERROR("Error while initializing MPI")//no semicolon
    return 1;
  }

  if (provided != MPI_THREAD_MULTIPLE) {
    SENSEI_ERROR("MPI provided the wrong level of thread support")//no semicolon
    return 1;
  }

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  size_t opt_nx, opt_ny, opt_nz;
  float opt_minx, opt_miny, opt_minz;
  float opt_maxx, opt_maxy, opt_maxz;
  size_t opt_nsteps;

  opt_nx = 8;
  opt_ny = 8;
  opt_nz = 8;

  opt_minx = -2.0f;
  opt_miny = -2.0f;
  opt_minz = +2.0f;
  opt_maxx = +2.0f;
  opt_maxy = +2.0f;
  opt_maxz = +4.0f;

  opt_nsteps = 8;

#define ARGLOOP \
  if (char *ARGVAL=nullptr) \
    ; \
  else \
    for (int ARGIND=1; ARGIND<argc; ARGVAL=NULL, ++ARGIND) \
      if (0) \
        ;

#define ARG(s) \
      else if (strncmp(argv[ARGIND], s, sizeof(s)) == 0 && ++ARGIND < argc && (ARGVAL = argv[ARGIND], 1))

  ARGLOOP
  ARG("-nx") opt_nx = (size_t)std::stoull(ARGVAL);
  ARG("-ny") opt_ny = (size_t)std::stoull(ARGVAL);
  ARG("-nz") opt_nz = (size_t)std::stoull(ARGVAL);
  ARG("-minx") opt_minx = std::stof(ARGVAL);
  ARG("-miny") opt_miny = std::stof(ARGVAL);
  ARG("-minz") opt_minz = std::stof(ARGVAL);
  ARG("-maxx") opt_maxx = std::stof(ARGVAL);
  ARG("-maxy") opt_maxy = std::stof(ARGVAL);
  ARG("-maxz") opt_maxz = std::stof(ARGVAL);
  ARG("-nsteps") opt_nsteps = (size_t)std::stoull(ARGVAL);

#undef ARG
#undef ARGLOOP

  std::fprintf(stderr, "Configuration:\n");
  std::fprintf(stderr, "  %zu x values in range [%+0.2f, %+0.2f]\n", opt_nx, opt_minx, opt_maxx);
  std::fprintf(stderr, "  %zu y values in range [%+0.2f, %+0.2f]\n", opt_ny, opt_miny, opt_maxy);
  std::fprintf(stderr, "  %zu z values in range [%+0.2f, %+0.2f]\n", opt_nz, opt_minz, opt_maxz);
  std::fprintf(stderr, "  %zu steps\n", opt_nsteps);

  Mandelbrot mandelbrot(opt_nx, opt_ny, opt_nz, {
    opt_minx, opt_miny, opt_minz,
    opt_maxx, opt_maxy, opt_maxz,
  });

  using AnalysisAdaptor = ospSensei::OSPRayVisualization;
  vtkNew<AnalysisAdaptor> analysisAdaptor;
  analysisAdaptor->SetCommunicator(MPI_COMM_WORLD);
  analysisAdaptor->SetMode(AnalysisAdaptor::VOLUME);
  analysisAdaptor->SetMeshName("mandelbrot");
  analysisAdaptor->SetArrayName("nsteps");
  analysisAdaptor->Initialize();

  for (size_t iter=0; iter<1; ++iter) {
    mandelbrot.step(opt_nsteps);

    if (rank == 0) {
      if (opt_nx * opt_ny * opt_nx < 1024UL) {
        mandelbrot.debug(Mandelbrot::OnlyNsteps);
      }
    }

    vtkSmartPointer<vtkUniformGrid> uniformGrid = mandelbrot.vtk();

    vtkNew<vtkUniformGridAMR> uniformGridAMR;
    int blocksPerLevel[] = { size };
    uniformGridAMR->Initialize(1 /* nLevels */, blocksPerLevel);
    uniformGridAMR->SetDataSet(0 /* levels */, rank /* index */, uniformGrid);

    using DataAdaptor = sensei::VTKDataAdaptor;
    vtkNew<DataAdaptor> dataAdaptor;
    dataAdaptor->SetDataTime(iter);
    dataAdaptor->SetDataTimeStep(iter * opt_nsteps);
    dataAdaptor->SetDataObject("mandelbrot", uniformGridAMR);

    analysisAdaptor->Execute(dataAdaptor);

    dataAdaptor->ReleaseData();
  }

  analysisAdaptor->Finalize();
  analysisAdaptor.Reset();

  MPI_Finalize();

  return 0;
}
// vim: ts=2:sts=2
