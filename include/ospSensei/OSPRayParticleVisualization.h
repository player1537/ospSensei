/**
 *
 */

#ifndef __OSPSENSEI_OSPRAY_PARTICLE_VISUALIZATION_H__
#define __OSPSENSEI_OSPRAY_PARTICLE_VISUALIZATION_H__

// SENSEI
#include <AnalysisAdaptor.h>

// VTK
#include <vtkSetGet.h>
#include <vtkObject.h>


//---

namespace ospSensei {

class OSPRayParticleVisualization : public sensei::AnalysisAdaptor {
public:
  senseiBaseTypeMacro(OSPRayParticleVisualization, AnalysisAdaptor);
  static OSPRayParticleVisualization *New();

  virtual void SetWidth(int _arg) {
    this->Width = _arg;
  }

  virtual void SetHeight(int _arg) {
    this->Height = _arg;
  }

  virtual void SetMeshName(const char *_arg) {
    this->MeshName = _arg;
  }

  virtual void SetArrayName(const char *_arg) {
    this->ArrayName = _arg;
  }

  int Initialize();
  bool Execute(sensei::DataAdaptor *) override;
  int Finalize() override;

private:
  int Width;
  int Height;
  std::string MeshName;
  std::string ArrayName;

  struct InternalsType;
  InternalsType *Internals;
};

} /* namespace ospSensei */


//---

#endif /* __OSPSENSEI_OSPRAY_PARTICLE_VISUALIZATION_H__ */
