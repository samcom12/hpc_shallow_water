#include "utils.h"

void cpy_to(double *target, const double *source, int numElements){
    for(int i=0; i<numElements; i++){
        target[i] = source[i];
    }
}

void read_file(string filename, double *array, size_t size)
{
    // open file
    ifstream file(filename, ios::in | ios::binary);

    // get length of file:
    file.seekg(0, file.end);
    int length = file.tellg();
    file.seekg(0, file.beg);

    // read data from file
    file.read(reinterpret_cast<char*>(array), length);
    file.close();
}

void write_file(string filename, double *array, size_t size)
{
    // open file
    std::ofstream file(filename, std::ios::out | std::ios::binary);

    // write data to file
    file.write(reinterpret_cast<char*>(&array[0]), size*size*sizeof(double));
    file.close();
}

void read_data(string filename, double *H, double *HU, double *HV, double *Zdx, double *Zdy, size_t size)
{
    read_file(filename + "_h.bin", H, size);
    read_file(filename + "_hu.bin", HU, size);
    read_file(filename + "_hv.bin", HV, size);
    read_file(filename + "_Zdx.bin", Zdx, size);
    read_file(filename + "_Zdy.bin", Zdy, size);
}

void swap(double* &a, double* &b)
{
    double *temp = a;
    a = b;
    b = temp;
}

void swap(double *H, double *HU, double *HV, double *Ht, double *HUt, double *HVt)
{
    swap(H, Ht);
    swap(HU, HUt);
    swap(HV, HVt);
}


