# QubiSim.jl

[![Build Status](https://github.com/bmorsink/QubiSim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/bmorsink/QubiSim.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://bmorsink.github.io/QubiSim/)
[![Version](https://img.shields.io/github/v/release/bmorsink/QubiSim.jl?label=version)](https://github.com/bmorsink/QubiSim.jl/releases)
[![Downloads](https://img.shields.io/github/downloads/bmorsink/QubiSim.jl/total)](https://github.com/bmorsink/QubiSim.jl/releases)
[![License](https://img.shields.io/github/license/bmorsink/QubiSim.jl)](LICENSE)

A Julia package for simulating quantum algorithms using a Schrödinger-based approach.

QubiSim supports unitary evolution, quantum channels, and measurements on qubits, using both state vectors and density operators. It also includes visualization tools for quantum states and measurement outcomes.

---

## Installation

```julia
using Pkg
Pkg.add("QubiSim")
```

---

## Quick example: Bell state measurement

```julia
using QubiSim

# Create circuit
qc = createQuantumCircuit(2)
hGate!(qc, 1)
cnotGate!(qc, 1, 2)

# Add measurement
measureGate!(qc, [1, 2])

# Compile and run
qp = compileQuantumCircuit(qc)
iqs = createInitialQubitState(vector, [1. 0.; 1. 0.])
qo = runQuantumProgram(qp, iqs, 1000)

# Plot measurement outcomes
probeMeasureOutcome(qo, 4, "Bell state measurement outcomes")
```

This produces the expected correlated outcomes of a Bell state.

---

## What can you do with QubiSim?

* Build quantum circuits with a wide range of gates
* Simulate full quantum state evolution
* Work with both **state vectors** (pure states) and **density operators** (mixed states)
* Model **noise and decoherence** via quantum channels (CPTP maps with Kraus operators)
* Perform **projective (PVM)** and **generalized (POVM)** measurements
* Analyze and visualize results (histograms, Bloch vectors, probabilities, entropy and purity)

---

## Conceptual workflow

QubiSim follows a clear 4-step structure:

1. Define a quantum circuit
2. Compile it into a quantum program
3. Execute the program on an initial state
4. Analyze the results

---

## Design principles

QubiSim balances physical correctness with execution efficiency:

- **Separation of structure and execution**  
  Circuits define *what* is computed; compiled programs define *how* it runs, enabling optimization.

- **Operator-based formalism**  
  State evolution, measurements and quantum channels are modeled via unitary matrices and Kraus operators.

- **Compiled execution model**  
  Circuits are lowered to time-ordered operations, reducing runtime overhead and enabling operation fusion.

- **Lazy, memory-aware construction**  
  Operators are generated on demand, avoiding unnecessary allocations.

- **Type-stable, dispatch-driven design**  
  Type-stable structures and multiple dispatch ensure performance and extensibility.

- **Extensible architecture**  
  New gates, operations and states can be added without modifying the core.

---

## Documentation

QubiSim provides documentation organized into five sections that take you from first use and hands-on learning to conceptual understanding, internals and full reference:

- **Getting Started** — installation and basic usage  
- **Bell State Tutorial** — a step-by-step detailed walkthrough 
- **Guide** — concepts, workflow and core components 
- **Internals** — software architecture and design
- **API Reference** — complete function and type documentation

👉 https://bmorsink.github.io/QubiSim.jl/

---

## Examples included

QubiSim supports simulation of a broad class of quantum algorithms, protocols and experiments, including:

* Bell state experiments
* Quantum teleportation
* Grover’s algorithm
* Quantum phase estimation
* Quantum key distribution
* CHSH game
* Quantum erasure

---

## Status

This package is under active development.

---

## License

MIT License
