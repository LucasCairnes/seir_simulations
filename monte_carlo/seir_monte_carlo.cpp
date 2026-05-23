#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <pybind11/pybind11.h>

namespace py = pybind11;

struct Agent {
    int x_pos, y_pos, state;
    Agent(int x, int y, int s)
        : x_pos(x), y_pos(y), state(s) {}
}; 

class System {
    private:
        int length;
        int s, e, i, r;
        double beta, sigma, gamma;

        std::vector<int> lattice; 
        std::vector<Agent> agents;

        std::mt19937 rng {42};
        std::uniform_int_distribution<int> creation_rng;
        std::uniform_int_distribution<int> movement_rng;
        std::uniform_int_distribution<int> direction_rng;
        std::uniform_real_distribution<double> prob_rng;
        
        int get_index(int x_value, int y_value) const;
        int boundary_check(int coord) const;
        
        int random_index();
        int random_movement();
        int random_direction();
        bool random_state(double probability);
        void populate_lattice(int agent_count, double s_0, double e_0, double i_0, double r_0);
        void move_agent(Agent& agent, int x_destination, int y_destination);
        void update_state(Agent& agent); 

    public:
        System(int length, int agent_count, double s_0, double e_0, double i_0, double r_0, double beta, double sigma, double gamma)
            : length(length), creation_rng(0, length - 1), movement_rng(-1, 1), direction_rng(0, 1), prob_rng(0.0, 1.0), beta(beta), sigma(sigma), gamma(gamma) { 
            
            if (std::abs(s_0 + e_0 + i_0 + r_0 - 1.0) > 0.001) {
                throw std::invalid_argument("SEIR values must add to 1.0");
            }
            if (agent_count > length * length) {
                throw std::invalid_argument("Agent count must not exceed the number of tiles.");
            }
            if (length < 0) {
                throw std::invalid_argument("Length must be a positive integer");
            }
            
            lattice.resize(length*length, 0);
            populate_lattice(agent_count, s_0, e_0, i_0, r_0); 
        };
        void run_sim(int MCS, const std::string& seir_filename, const std::string& lattice_filename);
};

int System::random_index() { return creation_rng(rng); }
int System::random_movement() { return movement_rng(rng); }
int System::random_direction() { return direction_rng(rng); }

bool System::random_state(double probability) {
    return (prob_rng(rng) < probability);  
}

int System::get_index(int x_value, int y_value) const {
    return (x_value * length) + y_value;
} 

void System::populate_lattice(int agent_count, double s_0, double e_0, double i_0, double r_0) {
    int x_val, y_val, state;
    s = static_cast<int>(std::floor(agent_count * s_0));
    e = static_cast<int>(std::floor(agent_count * e_0));
    i = static_cast<int>(std::floor(agent_count * i_0));
    r = agent_count - s - e - i; 

    agents.reserve(agent_count); 

    for (int j = 0; j < agent_count; j++) {
        if (j < s) state = 1;
        else if (j < s + e) state = 2;
        else if (j < s + e + i) state = 3;
        else state = 4;
        
        do {
            x_val = random_index();
            y_val = random_index();
        } while (lattice[get_index(x_val, y_val)] != 0); 

        agents.emplace_back(x_val, y_val, state);
        lattice[get_index(x_val, y_val)] = state; 
    }
}

int System::boundary_check(int coord) const {
    if (coord < 0) return length - 1;
    if (coord >= length) return 0;
    return coord; 
}

void System::move_agent(Agent& agent, int x_destination, int y_destination) {
    x_destination = boundary_check(x_destination);
    y_destination = boundary_check(y_destination);

    if (lattice[get_index(x_destination, y_destination)] == 0) {
        lattice[get_index(x_destination, y_destination)] = lattice[get_index(agent.x_pos, agent.y_pos)];
        lattice[get_index(agent.x_pos, agent.y_pos)] = 0; 
    
        agent.x_pos = x_destination;
        agent.y_pos = y_destination;
    }
}

void System::update_state(Agent& agent) {
    int current_state = agent.state;

    if (current_state == 4) return;
    
    bool state_change = false;
    int i_neighbours = 0, agent_x = agent.x_pos, agent_y = agent.y_pos;

    if (current_state == 1) {
        if (lattice[get_index(boundary_check(agent_x + 1), agent_y)] == 3) i_neighbours++;
        if (lattice[get_index(boundary_check(agent_x - 1), agent_y)] == 3) i_neighbours++;
        if (lattice[get_index(boundary_check(agent_x), boundary_check(agent_y + 1))] == 3) i_neighbours++;
        if (lattice[get_index(boundary_check(agent_x), boundary_check(agent_y - 1))] == 3) i_neighbours++; 

        if (i_neighbours == 0) return;

        state_change = random_state(1.0 - std::pow(1.0 - beta, i_neighbours));
    }
    else {
        if (current_state == 2) state_change = random_state(sigma);
        else if (current_state == 3) state_change = random_state(gamma);
    }

    if (state_change) { 
        agent.state = current_state + 1;
        lattice[get_index(agent_x, agent_y)] = current_state + 1;

        if (current_state == 1) {s--; e++;}
        else if (current_state == 2) {e--; i++;}
        else if (current_state == 3) {i--; r++;}
    } 
}

void System::run_sim(int MCS, const std::string& seir_filename, const std::string& lattice_filename) {
    int step = 0;

    if (MCS < 0) {
        throw std::invalid_argument("MCS must be a positive integer");
    }

    std::vector<double> s_values, e_values, i_values, r_values; // Upgraded to double
    std::vector<int> steps;
    std::vector<std::string> lattices;

    s_values.reserve(MCS + 1); 
    e_values.reserve(MCS + 1);
    i_values.reserve(MCS + 1);
    r_values.reserve(MCS + 1);
    steps.reserve(MCS + 1);

    s_values.emplace_back(s);  
    e_values.emplace_back(e);
    i_values.emplace_back(i);   
    r_values.emplace_back(r);
    steps.emplace_back(step);

    std::ostringstream lattice_snapshot;
    lattice_snapshot << step << ",";
    for (int status : lattice) {
        lattice_snapshot << status << " ";
    }
    lattices.emplace_back(lattice_snapshot.str());

    for (int k = 0; k < MCS; k++) {
        std::shuffle(agents.begin(), agents.end(), rng); 
        for (Agent& curr_agent : agents) {
            int movement = random_movement();
            if (movement == 0) continue;

            int direction = random_direction();
            int x_destination, y_destination;

            if (direction == 1) {
                x_destination = curr_agent.x_pos + movement;
                y_destination = curr_agent.y_pos;
            }
            else {
                x_destination = curr_agent.x_pos;
                y_destination = curr_agent.y_pos + movement;
            }

            move_agent(curr_agent, x_destination, y_destination);
        } 

        for (Agent& curr_agent : agents) {
            update_state(curr_agent);
        } 

        step++;

        s_values.emplace_back(s);
        e_values.emplace_back(e);
        i_values.emplace_back(i);   
        r_values.emplace_back(r);
        steps.emplace_back(step);

        std::ostringstream lattice_snapshot_loop;
        lattice_snapshot_loop << step << ",";
        for (int status : lattice) {
            lattice_snapshot_loop << status << " ";
        }
        lattices.emplace_back(lattice_snapshot_loop.str());
    }

    std::ofstream seir_out(seir_filename);
    seir_out << "Monte Carlo step,susceptible,exposed,infected,recovered\n";
    
    for (size_t idx = 0; idx < steps.size(); idx++) {
        seir_out << steps[idx] << "," << s_values[idx] << "," << e_values[idx] << "," << i_values[idx] << "," << r_values[idx] << "\n";
    }
    seir_out.close();

    std::ofstream lattice_out(lattice_filename);
    lattice_out << "Monte Carlo step,lattice\n";
    for (const auto& row : lattices) {
        lattice_out << row << "\n";
    }
    lattice_out.close();
}

PYBIND11_MODULE(seir_monte_carlo, m) {
    m.doc() = "A monte carlo SEIR simulation";

    py::class_<System>(m, "System")
        .def(py::init<int, int, double, double, double, double, double, double, double>(), 
            py::arg("length"),
            py::arg("agent_count"),
            py::arg("s_0"),
            py::arg("e_0"),
            py::arg("i_0"),
            py::arg("r_0"),
            py::arg("beta"),
            py::arg("sigma"),
            py::arg("gamma")  
        )
        .def("run_sim", 
            &System::run_sim,
            py::arg("MCS"),
            py::arg("seir_filename"),
            py::arg("lattice_filename")
        );
}