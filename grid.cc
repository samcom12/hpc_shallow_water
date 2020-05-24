/* -------------------------------------------------------------------------- */
#include "grid.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
Grid::Grid(int ny_, int nx_) : ny(ny_),nx(nx_), m_storage(ny_ * nx_) { clear(); }

/* -------------------------------------------------------------------------- */
void Grid::read_file(string filename)
{
    // open file
    ifstream file(filename, ios::in | ios::binary);

    // get length of file:
    file.seekg(0, file.end);
    int length = file.tellg();
    file.seekg(0, file.beg);

    // allocate memory
    auto memblock = new char [length];

    // read data from file
    file.read(memblock, length);
    file.close();

    double* values = (double*)memblock;
    for (auto i = 0; i < ny; i++)
        for (auto j = 0; j < nx; j++)
            m_storage[i * nx + j] = values[i * nx + j];
}

/* -------------------------------------------------------------------------- */
void Grid::write_file(string filename)
{
    // open file
    std::ofstream file(filename, std::ios::out | std::ios::binary);

    // write data to file
    for (auto i = 0; i < ny; i++)
        for (auto j = 0; j < nx; j++)
            file.write((char*)(&m_storage[i * nx + j]), sizeof(double));

    file.close();
}

/* -------------------------------------------------------------------------- */
void Grid::clear() { std::fill(m_storage.begin(), m_storage.end(), 0.); }

/* -------------------------------------------------------------------------- */
int Grid::get_nx() const { return nx; }
int Grid::get_ny() const { return ny; }
