/* -------------------------------------------------------------------------- */
#include "grid.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
Grid::Grid(int ny_, int nx_, MPI_Comm communicator_) : ny(ny_),nx(nx_),
    m_storage(ny_ * nx_), communicator(communicator_) { clear(); }

/* -------------------------------------------------------------------------- */
void Grid::read_file(string filename)
{
    int prank, psize;
    MPI_Comm_rank(communicator, &prank);
    MPI_Comm_size(communicator, &psize);
    MPI_File file;
    MPI_Status status;

    int h = ny, start = 0;

    // removing the ghosts from the height
    if (psize > 1) {
        h = (prank == 0 || prank == psize - 1 ? h - 1 : h - 2);
        start = (prank == 0 ? 0 : 1);
    }

    // Gathering the size of every processors, this could be done as in the
    // constructor of the Simulation instead
    std::vector<int> size_per_proc(psize);
    MPI_Allgather(&h, 1, MPI_INT, size_per_proc.data(), 1, MPI_INT, communicator);

    // determining the local offset
    int offset_h = 0;
    for (int i = 0; i < prank; ++i) {
        offset_h += size_per_proc[i];
    }

    MPI_File_open(communicator, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    MPI_File_read_at(file, offset_h * nx * sizeof(double), &m_storage[start * nx], h * nx, MPI_DOUBLE, &status);
    MPI_File_close(&file);
}

/* -------------------------------------------------------------------------- */
void Grid::write_file(string filename)
{
    int prank, psize;
    MPI_Comm_rank(communicator, &prank);
    MPI_Comm_size(communicator, &psize);
    MPI_File file;
    MPI_Status status;

    int h = ny, start = 0;

    // removing the ghosts from the height
    if (psize > 1) {
        h = (prank == 0 || prank == psize - 1 ? h - 1 : h - 2);
        start = (prank == 0 ? 0 : 1);
    }

    // Gathering the size of every processors, this could be done as in the
    // constructor of the Simulation instead
    std::vector<int> size_per_proc(psize);
    MPI_Allgather(&h, 1, MPI_INT, size_per_proc.data(), 1, MPI_INT, communicator);

    // determining the local offset
    int offset_h = 0;
    for (int i = 0; i < prank; ++i) {
        offset_h += size_per_proc[i];
    }

    MPI_File_open(communicator, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
    MPI_File_set_size(file, nx * nx * sizeof(double));   /*in bytes*/
    MPI_File_write_at(file, offset_h * nx * sizeof(double), &m_storage[start * nx], h * nx, MPI_DOUBLE, &status);
    MPI_File_close(&file);
}

/* -------------------------------------------------------------------------- */
void Grid::clear() { std::fill(m_storage.begin(), m_storage.end(), 0.); }

/* -------------------------------------------------------------------------- */
void Grid::resize(int ny_, int nx_) {
    ny = ny_;
    nx = nx_;
    m_storage.resize(ny_ * nx_);
}
/* -------------------------------------------------------------------------- */
int Grid::get_nx() const { return nx; }
int Grid::get_ny() const { return ny; }
