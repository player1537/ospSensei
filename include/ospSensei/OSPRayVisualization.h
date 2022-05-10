/**
 *
 */

#ifndef __OSPSENSEI_OSPRAYVISUALIZATION_H__
#define __OSPSENSEI_OSPRAYVISUALIZATION_H__


// SENSEI
#include <AnalysisAdaptor.h>


namespace ospSensei {


class OSPRayVisualization : public sensei::AnalysisAdaptor {
public:
  static OSPRayVisualization *New();

  senseiTypeMacro(OSPRayVisualization, AnalysisAdaptor);

  int Initialize();

  bool Execute(sensei::DataAdaptor *) override;
  int Finalize() override;

private:
  bool Execute();

  struct InternalsType;
  InternalsType *Internals;
};


} /* namespace ospSensei */


#endif /* __OSPSENSEI_OSPRAYVISUALIZATION_H__ */
// vim: ts=2:sts=2
