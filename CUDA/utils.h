#ifndef UTILS_H
#define UTILS_H

#include <cstring>
#include <sstream>
#include <fstream>
using namespace std;

void cpy_to(double *target, const double *source, int numElements);
void read_file(string filename, double *array, size_t size);
void write_file(string filename, double *array, size_t size);
void read_data(string filename, double *H, double *HU, double *HV, double *Zdx, double *Zdy, size_t size);
void swap(double *H, double *HU, double *HV, double *Ht, double *HUt, double *HVt);

#endif
