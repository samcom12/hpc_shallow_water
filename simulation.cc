/* -------------------------------------------------------------------------- */
#include "simulation.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <iostream>
#include <string>
/* -------------------------------------------------------------------------- */
using std::to_string, std::max, std::cout;
/* -------------------------------------------------------------------------- */
Simulation::Simulation(int nx_, double size_, double tend_)
    : nx(nx_), size(size_), epsilon(1e-4), dx(size_ / nx_),
      h(nx_, nx_), hu(nx_, nx_), hv(nx_, nx_), zdx(nx_, nx_), zdy(nx_, nx_),
      Tend(tend_) {}

/* -------------------------------------------------------------------------- */
void Simulation::set_initial_conditions() {
    auto filename = to_string(Tend);
    filename = filename.substr(0, filename.find(".") + 2);
    filename = "../data/Data_nx" + to_string(nx) + "_" + to_string(size) + "km_T" + filename;

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
    filename = "../data/Solution_nx" + to_string(nx) + "_" + to_string(size) + "km_T" + filename + "_h.bin";

    h.current().write_file(filename);
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
      save_results();
  }

  return nt;
}

/* -------------------------------------------------------------------------- */
void Simulation::set_treshold(double epsilon_) { epsilon = epsilon_; }

/* -------------------------------------------------------------------------- */
double Simulation::treshold() const { return epsilon; }

/* -------------------------------------------------------------------------- */
void Simulation::boundary_conditions() {
    Grid & ht = h.old();
    Grid & hvt = hv.old();
    Grid & hut = hu.old();

    for (auto i = 0; i < nx; i++) {
        ht(0, i) = ht(1, i);
        ht(nx - 1, i) = ht(nx - 2, i);
        ht(i, 0) = ht(i, 1);
        ht(i, nx - 1) = ht(i, nx - 2);

        hvt(0, i) = hvt(1, i);
        hvt(nx - 1, i) = hvt(nx - 2, i);
        hvt(i, 0) = hvt(i, 1);
        hvt(i, nx - 1) = hvt(i, nx - 2);

        hut(0, i) = hut(1, i);
        hut(nx - 1, i) = hut(nx - 2, i);
        hut(i, 0) = hut(i, 1);
        hut(i, nx - 1) = hut(i, nx - 2);
    }
}

/* -------------------------------------------------------------------------- */
void Simulation::compute_step() {
    double g = 127267.2; // Gravity, 9.82*(3.6)^2*1000 in [km / hr^2]

    double max_hu, max_hv, mu = 0, dt, C;

    Grid &ho = h.old(), &hvo = hv.old(), &huo = hu.old();
    Grid &ht = h.current(), &hvt = hv.current(), &hut = hu.current();

    // Compute the time-step length
    for (auto i = 0; i < nx; i++) {
        for (auto j = 0; j < nx; j++) {
            max_hu = max(abs(huo(i, j) / ho(i, j) + sqrt(ho(i, j) * g)),
                         abs(huo(i, j) / ho(i, j) - sqrt(ho(i, j) * g)));
            max_hv = max(abs(hvo(i, j) / ho(i, j) + sqrt(ho(i, j) * g)),
                         abs(hvo(i, j) / ho(i, j) - sqrt(ho(i, j) * g)));
            mu = max(mu, sqrt(max_hu * max_hu + max_hv * max_hv));
        }
    }

    dt = dx / (sqrt(2) * mu);

    if (T + dt > Tend)
        dt = Tend - T;

    // Print status
    cout << "Computing T: " << T + dt << ". " << 100 * (T + dt) / Tend << "%\n";

    // Enforce boundary condition
    boundary_conditions();

    // Compute a time-step
    C = 0.5 * dt / dx;
    for (auto i = 1; i < nx - 1; i++) {
        for (auto j = 1; j < nx - 1; j++) {
            ht(i, j) = 0.25 * (ho(i, j + 1) + ho(i, j - 1) + ho(i + 1, j) + ho(i - 1, j)) +
                       C * (huo(i, j - 1) - huo(i, j + 1) + hvo(i - 1, j) - huo(i + 1, j));

            hut(i, j) = 0.25 * (huo(i, j - 1) + huo(i, j + 1) + huo(i - 1, j) + huo(i + 1, j)) -
                        dt * g * ht(i, j) * zdx(i, j) +
                        C * (huo(i, j - 1) * huo(i, j - 1) / ho(i, j - 1) +
                             0.5 * g * ho(i, j - 1) * ho(i, j - 1) -
                             huo(i, j + 1) * huo(i, j + 1) / ho(i, j + 1) -
                             0.5 * g * ho(i, j + 1) * ho(i, j + 1)) +
                        C * (huo(i - 1, j) * hvo(i - 1, j) / ho(i - 1, j) -
                             huo(i + 1, j) * hvo(i + 1, j) / ho(i + 1, j));
            hvt(i, j) = 0.25 * (hvo(i, j - 1) + hvo(i, j + 1) + hvo(i - 1, j) + hvo(i + 1, j)) -
                        dt * g * ht(i, j) * zdy(i, j) +
                        C * (hvo(i - 1, j) * hvo(i - 1, j) / ho(i - 1, j) +
                             0.5 * g * ho(i - 1, j) * ho(i - 1, j) -
                             hvo(i + 1, j) * hvo(i + 1, j) / ho(i + 1, j) -
                             0.5 * g * ho(i + 1, j) * ho(i + 1, j)) +
                        C * (huo(i, j - 1) * hvo(i, j - 1) / ho(i, j - 1) -
                             huo(i, j + 1) * hvo(i, j + 1) / ho(i, j + 1));
        }
    }

    // Impose tolerances
    for (auto i = 1; i < nx - 1; i++) {
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