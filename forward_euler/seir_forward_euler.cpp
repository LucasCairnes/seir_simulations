#include <string>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <pybind11/pybind11.h>

namespace py = pybind11; 

double ds_dt(double i, double s, double beta) {
    return -beta * i * s;
}
double de_dt(double i, double s, double beta, double sigma, double e) {
    return beta * i * s - sigma * e;
}
double di_dt(double i, double sigma, double e, double gamma) {
    return sigma * e - gamma * i;
}
double dr_dt(double i, double gamma) {
    return gamma * i;
} 

int seir_forward_euler(double beta, double sigma, double gamma, double s_0, double e_0, double i_0, double r_0, double step, double t_final, const std::string& seir_filename) {
    double t = 0.0, s = s_0, e = e_0, i = i_0, r = r_0;
    
    if (std::abs(s_0 + e_0 + i_0 + r_0 - 1.0) > 0.001) {
        throw std::invalid_argument("SEIR values must add to 1.0");
    }
    if (step > t_final) {
        throw std::invalid_argument("Step must be smaller than t_final.");
    }

    std::vector<double> s_values, e_values, i_values, r_values, t_values;

    s_values.emplace_back(s);  
    e_values.emplace_back(e);
    i_values.emplace_back(i);   
    r_values.emplace_back(r);
    t_values.emplace_back(t);

    while (t < t_final) {
        double ds = ds_dt(i, s, beta);
        double de = de_dt(i, s, beta, sigma, e);
        double di = di_dt(i, sigma, e, gamma);
        double dr = dr_dt(i, gamma); 

        s += ds * step;
        e += de * step;
        i += di * step;
        r += dr * step;
        t += step;

        s_values.emplace_back(s);
        e_values.emplace_back(e);
        i_values.emplace_back(i);   
        r_values.emplace_back(r);
        t_values.emplace_back(t);
    }

    std::ofstream seir_out(seir_filename);
    seir_out << "time,susceptible,exposed,infected,recovered\n";

    for (size_t j = 0; j < t_values.size(); j++) {
        seir_out << t_values[j] << "," << s_values[j] << "," << e_values[j] << "," << i_values[j] << "," << r_values[j] << "\n";
    }
    seir_out.close();

    return 0;
}

PYBIND11_MODULE(seir_forward_euler, m) {
    m.doc() = "A module for solving SEIR ODEs.";
    m.def("forward_euler", &seir_forward_euler, "A function to solve the SEIR ODEs with the forward Euler method.", 
        py::arg("beta"),
        py::arg("sigma"),
        py::arg("gamma"),
        py::arg("s_0"),
        py::arg("e_0"),
        py::arg("i_0"),
        py::arg("r_0"),
        py::arg("step"),
        py::arg("t_final"),
        py::arg("seir_filename") 
    );
}