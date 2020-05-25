/* -------------------------------------------------------------------------- */
#include "simulation.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <iostream>
#include <string>
#include <algorithm>
/* -------------------------------------------------------------------------- */
using std::to_string;
using std::max;
using std::cout;
/* -------------------------------------------------------------------------- */
Simulation::Simulation(int nx_, double size_, double tend_, MPI_Comm communicator_)
    : nx(nx_), size(size_), dx(size_ / nx_), Tend(tend_), epsilon(1e-4),
      h(nx_, nx_, communicator_), hu(nx_, nx_, communicator_), hv(nx_, nx_, communicator_),
      zdx(nx_, nx_, communicator_), zdy(nx_, nx_, communicator_),
      communicator(communicator_){

    // retrieving the number of proc and the rank in the proc pool
    MPI_Comm_rank(communicator, &prank);
    MPI_Comm_size(communicator, &psize);

    // computation of the local size of the grid the remainder is spread equally
    // on the first processors
    local_ny = nx / psize + (prank < nx % psize ? 1 : 0);
    local_nx = nx;

    // adding the ghosts lines if needed
    if (psize > 1)
        local_ny += (prank == 0 || prank == psize - 1) ? 1 : 2;

    // computing the offsets of the local grid in the global one
    offset_y = (nx / psize) * prank + (prank < nx % psize ? prank : nx % psize);
    offset_x = 0;

    // resizing the different grids
    h.resize(local_ny, local_nx);
    hu.resize(local_ny, local_nx);
    hv.resize(local_ny, local_nx);
    zdx.resize(local_ny, local_nx);
    zdy.resize(local_ny, local_nx);

    // determining the rank of the neighbors
    north_prank = (prank == 0 ? MPI_PROC_NULL : prank - 1);
    south_prank = (prank == (psize - 1) ? MPI_PROC_NULL : prank + 1);
    std::cout << prank << " " << nx << " " << nx << " "
               << local_ny << " " << local_nx << " " << offset_y << " "
               << offset_x << " " << north_prank << " " << south_prank
               << std::endl;
}

/* -------------------------------------------------------------------------- */
void Simulation::set_initial_conditions(std::string filename) {
    h.old().read_file(filename + "_h.bin");
    hu.old().read_file(filename + "_hu.bin");
    hv.old().read_file(filename + "_hv.bin");

    zdx.read_file(filename + "_Zdx.bin");
    zdy.read_file(filename + "_Zdy.bin");
}

/* -------------------------------------------------------------------------- */
void Simulation::save_results() {
    auto filename = to_string(Tend);
    filename = filename.substr(0, filename.find(".") + 2);
    filename = "./data/Solution_nx" + to_string(nx) + "_" + to_string(size) + "km_T" + filename + "_h.bin";

    Grid & ht = h.old();
    ht.write_file(filename);
}

/* -------------------------------------------------------------------------- */
int Simulation::compute() {
  int nt = 0;

  while (T < (double)Tend){
      compute_step();
      h.swap();
      hu.swap();
      hv.swap();
      ++nt;
  }

  save_results();
  return nt;
}

/* -------------------------------------------------------------------------- */
void Simulation::set_treshold(double epsilon_) { epsilon = epsilon_; }

/* -------------------------------------------------------------------------- */
double Simulation::treshold() const { return epsilon; }

/* -------------------------------------------------------------------------- */
void Simulation::boundary_conditions() {
    int prank, psize;
    MPI_Comm_rank(communicator, &prank);
    MPI_Comm_size(communicator, &psize);


    int height = local_ny;
    int width = local_nx;

    int start = 0;

    // removing the ghosts from the height
    if (psize > 1) {
        height = (prank == 0 || prank == psize - 1 ? height - 1 : height - 2);
        start = (prank == 0 ? 0 : 1);
    }

    Grid & ht = h.old();
    Grid & hvt = hv.old();
    Grid & hut = hu.old();

    for (auto i = start; i < start + height; i++) {
        ht(i, 0) = ht(i, 1);
        ht(i, nx - 1) = ht(i, nx - 2);

        hvt(i, 0) = hvt(i, 1);
        hvt(i, nx - 1) = hvt(i, nx - 2);

        hut(i, 0) = hut(i, 1);
        hut(i, nx - 1) = hut(i, nx - 2);
    }

    if (prank == 0) {
        for (auto i = 0; i < nx; i++) {
            ht(0, i) = ht(1, i);
            hvt(0, i) = hvt(1, i);
            hut(0, i) = hut(1, i);
        }
    }

    if (prank == psize - 1) {
        for (auto i = 0; i < nx; i++) {
            ht(local_ny - 1, i) = ht(local_ny - 2, i);
            hvt(local_ny - 1, i) = hvt(local_ny - 2, i);
            hut(local_ny - 1, i) = hut(local_ny - 2, i);
        }
    }
}

/* -------------------------------------------------------------------------- */
void Simulation::compute_step() {
    double g = 127267.2; // Gravity, 9.82*(3.6)^2*1000 in [km / hr^2]

    int prank, psize;
    MPI_Comm_rank(communicator, &prank);
    MPI_Comm_size(communicator, &psize);

    int height = local_ny;
    int width = local_nx;

    int start = 0;

    // removing the ghosts from the height
    if (psize > 1) {
        height = (prank == 0 || prank == psize - 1 ? height - 1 : height - 2);
        start = (prank == 0 ? 0 : 1);
    }

    // Gathering the size of every processors, this could be done as in the
    // constructor of the Simulation instead
    std::vector<int> size_per_proc(psize);
    MPI_Allgather(&height, 1, MPI_INT, size_per_proc.data(), 1, MPI_INT, communicator);

    // determining the local offset
    int offset_h = 0;
    for (int i = 0; i < prank; ++i) {
        offset_h += size_per_proc[i];
    }

    int total_h = offset_h;
    for (int i = prank; i < psize; ++i) {
        total_h += size_per_proc[i];
    }

    double max_hu, max_hv, mu = 0, dt = 0, C;

    Grid &ho = h.old(), &hvo = hv.old(), &huo = hu.old();
    Grid &ht = h.current(), &hvt = hv.current(), &hut = hu.current();


    // Compute the time-step length
    for (auto i = start; i < start + height; i++) {
        for (auto j = 0; j < nx; j++) {
            max_hu = max(abs(huo(i, j) / ho(i, j) + sqrt(ho(i, j) * g)),
                         abs(huo(i, j) / ho(i, j) - sqrt(ho(i, j) * g)));
            max_hv = max(abs(hvo(i, j) / ho(i, j) + sqrt(ho(i, j) * g)),
                         abs(hvo(i, j) / ho(i, j) - sqrt(ho(i, j) * g)));
            mu = max(mu, sqrt(max_hu * max_hu + max_hv * max_hv));
        }
    }

    std::vector<double> mu_vec(psize);
    MPI_Allgather(&mu, 1, MPI_DOUBLE, mu_vec.data(), 1, MPI_DOUBLE, communicator);
    mu = *std::max_element(mu_vec.begin(), mu_vec.end());
    dt = dx / (sqrt(2) * mu);

    if (T + dt > Tend)
        dt = Tend - T;

    if (prank == 0){
        // Print status
        cout << "Computing T: " << T + dt << ". " << 100 * (T + dt) / Tend << "%\n";
    }

    // Enforce boundary condition
    boundary_conditions();

    // Taking care of communications going up (so receiving from bottom)
    MPI_Sendrecv(&ho(1, 0), nx, MPI_DOUBLE, north_prank, 0,
                 &ho(local_ny - 1, 0), nx, MPI_DOUBLE, south_prank, 0,
                 communicator, MPI_STATUS_IGNORE);

    // Taking care of communications going down (so receiving from top)
    MPI_Sendrecv(&ho(local_ny - 2, 0), nx, MPI_DOUBLE, south_prank, 0,
                 &ho(0, 0), nx, MPI_DOUBLE, north_prank, 0,
                 communicator, MPI_STATUS_IGNORE);

    // Taking care of communications going up (so receiving from bottom)
    MPI_Sendrecv(&huo(1, 0), nx, MPI_DOUBLE, north_prank, 0,
                 &huo(local_ny - 1, 0), nx, MPI_DOUBLE, south_prank, 0,
                 communicator, MPI_STATUS_IGNORE);

    // Taking care of communications going down (so receiving from top)
    MPI_Sendrecv(&huo(local_ny - 2, 0), nx, MPI_DOUBLE, south_prank, 0,
                 &huo(0, 0), nx, MPI_DOUBLE, north_prank, 0,
                 communicator, MPI_STATUS_IGNORE);

    // Taking care of communications going up (so receiving from bottom)
    MPI_Sendrecv(&hvo(1, 0), nx, MPI_DOUBLE, north_prank, 0,
                 &hvo(local_ny - 1, 0), nx, MPI_DOUBLE, south_prank, 0,
                 communicator, MPI_STATUS_IGNORE);

    // Taking care of communications going down (so receiving from top)
    MPI_Sendrecv(&hvo(local_ny - 2, 0), nx, MPI_DOUBLE, south_prank, 0,
                 &hvo(0, 0), nx, MPI_DOUBLE, north_prank, 0,
                 communicator, MPI_STATUS_IGNORE);

    // Compute a time-step
    C = 0.5 * dt / dx;
    for (auto i = 1; i < local_ny - 1; i++)
        for (auto j = 1; j < nx - 1; j++) {
            ht(i, j) = 0.25 * (ho(i, j + 1) + ho(i, j - 1) + ho(i + 1, j) + ho(i - 1, j)) +
                       C * (huo(i - 1, j) - huo(i + 1, j) + hvo(i, j - 1) - hvo(i, j + 1));
        }

    // Taking care of communications going up (so receiving from bottom)
    MPI_Sendrecv(&ht(1, 0), nx, MPI_DOUBLE, north_prank, 0,
                 &ht(local_ny - 1, 0), nx, MPI_DOUBLE, south_prank, 0,
                 communicator, MPI_STATUS_IGNORE);

    // Taking care of communications going down (so receiving from top)
    MPI_Sendrecv(&ht(local_ny - 2, 0), nx, MPI_DOUBLE, south_prank, 0,
                 &ht(0, 0), nx, MPI_DOUBLE, north_prank, 0,
                 communicator, MPI_STATUS_IGNORE);

    for (auto i = 1; i < local_ny - 1; i++)
        for (auto j = 1; j < nx - 1; j++) {
            hut(i, j) = 0.25 * (huo(i, j - 1) + huo(i, j + 1) + huo(i - 1, j) + huo(i + 1, j)) -
                        dt * g * ht(i, j) * zdx(i, j) +
                        C * (huo(i - 1, j) * huo(i - 1, j) / ho(i - 1, j) +
                             0.5 * g * ho(i - 1, j) * ho(i - 1, j) -
                             huo(i + 1, j) * huo(i + 1, j) / ho(i + 1, j) -
                             0.5 * g * ho(i + 1, j) * ho(i + 1, j)) +
                        C * (huo(i, j - 1) * hvo(i, j - 1) / ho(i, j - 1) -
                             huo(i, j + 1) * hvo(i, j + 1) / ho(i, j + 1));
            hvt(i, j) = 0.25 * (hvo(i, j - 1) + hvo(i, j + 1) + hvo(i - 1, j) + hvo(i + 1, j)) -
                        dt * g * ht(i, j) * zdy(i, j) +
                        C * (hvo(i, j - 1) * hvo(i, j - 1) / ho(i, j - 1) +
                             0.5 * g * ho(i, j - 1) * ho(i, j - 1) -
                             hvo(i, j + 1) * hvo(i, j + 1) / ho(i, j + 1) -
                             0.5 * g * ho(i, j + 1) * ho(i, j + 1)) +
                        C * (huo(i - 1, j) * hvo(i - 1, j) / ho(i - 1, j) -
                             huo(i + 1, j) * hvo(i + 1, j) / ho(i + 1, j));
        }

    // Impose tolerances
    for (auto i = start; i < start + height; i++) {
        for (auto j = 1; j < nx - 1; j++) {
            if (ht(i, j) < 0)
                ht(i, j) = epsilon / 10;
            if (ht(i, j) <= epsilon) {
                hut(i, j) = 0;
                hvt(i, j) = 0;
            }
        }
    }

    // Update time T
    T = T + dt;
}