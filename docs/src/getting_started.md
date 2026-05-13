# Installation and basic usage  

---

## Installation

Install QubiSim using Julia’s package manager:

```julia
using Pkg
Pkg.add("QubiSim")
```

Load the package:

```julia
using QubiSim
```

---

## Basic workflow

A typical simulation in QubiSim consists of four steps:

1. Define a quantum circuit
2. Compile the circuit into a quantum program
3. Execute the program on an initial state
4. Analyze the results

---

## Minimal example

The following example constructs and executes a simple two-qubit circuit:

```julia
qc = createQuantumCircuit(2)
hGate!(qc, 1)
cnotGate!(qc, 1, 2)
measureGate!(qc, [1, 2])

qp  = compileQuantumCircuit(qc)
iqs  = createInitialQubitState(vector, [1. 0.; 1. 0.])
qo = runQuantumProgram(qp, iqs, 1000)

probeMeasureOutcome(qo, 4, "Bell state measurement outcomes")
```

This produces the expected correlated outcomes of a Bell state.

---

## Next steps

* See the [Bell State Tutorial](tutorials_bell_state.md) for a step-by-step detailed walkthrough
* Explore the [Guide](guide_overview.md) for concepts, workflow and core components
* Explore the [Internals](internals.md) for software architecture and design
* Consult the [API Reference](api.md) for complete function and type documentation
