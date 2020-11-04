#ifndef DOUBLE_BUFFER
#define DOUBLE_BUFFER

/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */
#include "grid.hh"
/* -------------------------------------------------------------------------- */

class DoubleBuffer {
public:
  DoubleBuffer(int m, int n);

  Grid & current();
  Grid & old();

  void swap();
private:

  std::unique_ptr<Grid> m_current;
  std::unique_ptr<Grid> m_old;
};

#endif /* DOUBLE_BUFFER */
