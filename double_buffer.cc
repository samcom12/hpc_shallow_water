/* -------------------------------------------------------------------------- */
#include "double_buffer.hh"
#include "grid.hh"
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
DoubleBuffer::DoubleBuffer(int m, int n, MPI_Comm communicator_)
    : m_current(new Grid(m, n, communicator_)), m_old(new Grid(m, n, communicator_)) {}

/* -------------------------------------------------------------------------- */
Grid & DoubleBuffer::current() { return *m_current; }

/* -------------------------------------------------------------------------- */
Grid & DoubleBuffer::old() { return *m_old; }

/* -------------------------------------------------------------------------- */
void DoubleBuffer::swap() {
  m_current.swap(m_old);
}

/* -------------------------------------------------------------------------- */
void DoubleBuffer::resize(int m, int n) {
    m_current->resize(m, n);
    m_old->resize(m, n);
}