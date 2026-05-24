# 🦠 C++ SEIR Simulations

![C++](https://img.shields.io/badge/C%2B%2B-00599C?style=flat&logo=c%2B%2B&logoColor=white)
![Python](https://img.shields.io/badge/Python-3776AB?style=flat&logo=python&logoColor=white)
![CMake](https://img.shields.io/badge/CMake-064F8C?style=flat&logo=cmake&logoColor=white)
![Jupyter](https://img.shields.io/badge/Jupyter-F37626?style=flat&logo=jupyter&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-11557c?style=flat&logo=matplotlib&logoColor=white)

A C++ implementation of the SEIR model (Susceptible, Exposed, Infected, Recovered). 

## Project Overview
This project compares two distinct models: a deterministic ODE model solved with the Forward Euler method, and a stochastic, agent-based Monte Carlo simulation on a periodic lattice. The simulation engines for both models are implemented in C++ and, through the use of pybind11, Python is used for analysis and visualisation

**[View My Experimental Findings](https://drive.google.com/file/d/16QFR0962qGyqNGi_RslJRaU5ySu5T3Fz/view?usp=sharing)**

---

## Repository Architecture

**ODE Model** `seir_forward_euler.cpp`
* Solves the SEIR differential equations numerically (Forward Euler)
* Takes initial parameters, sequentially calculates instantaneous rates, and then updates the SEIR value
* Outputs the data to .xlsx files to be analysed by the visualisation script

**Agent Model** `seir_monte_carlo.cpp`
* A discrete model on a periodic lattice which models interactions between individual agents
* Agents are represented by structs and the overall system is encapsulated by the System class
* During each simulation step, agents attempt to random walk
* The disease propagates probabilistically based on neighbouring agents 
* The program outputs .xlsx files for the SEIR values and for snapshots of the lattice

**Visualisation** `visualise_seir.ipynb`
* Imports the pybind11 objects and triggers the simulations with the chosen parameters
* Plots the SEIR trajectories with matplotlib and the lattices with FuncAnimation

---

## Setup

Clone the repository:

* `git clone https://github.com/LucasCairnes/seir_simulations`
* `cd seir_simulations`

Setup a Python environment and install the dependencies:

* `conda create -n "seir_simulation" python=3.13`
* `conda activate seir_simulation`
* `pip install -r requirements.txt`

Compile the C++ code:

* `cmake -S . -B build`
* `cmake --build build`

---

## Execution

* Launch `visualise_seir.ipynb`
* Input your chosen values into the parameters cell
* Execute the cells sequentially
