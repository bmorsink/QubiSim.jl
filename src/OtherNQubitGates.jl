############################################
# OtherNQubitGates

""" 
	bernoulliGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits,
		probability::Float64,
		U::UnitaryOperation; 
		accuracyCheckForUnitarity = 10*eps(1.0)
	) -> Nothing

Adds a **Bernoulli gate** to the `quantumCircuit` on the specified `qubits` with a user-supplied `probability` and  **unitary operation** ``U``.

During quantum program execution, for each shot, the **Bernoulli gate** applies the user-supplied **unitary operation** ``U`` with probability `probability`, while with probability `1-probability` the state is left unchanged.

The function performs an internal unitarity check to ensure that the user-supplied **unitary operation** ``U`` satisfies the condition ``U^†⋅U ≈ I`` within a specified numerical tolerance.

Under the hood, it creates a `BernoulliGate` holding a factory function that can generate the unitary matrix of the **unitary operation**, and adds 
it to the circuit. During compilation, this factory function simply returns the user-suppied **unitary operation** ``U`` as the produced unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices (interpreted using the circuit’s `settings`) on which the gate is applied.
- `probability::Float64` — The probability of applying the unitary operation.
- `U::UnitaryOperation` — A structure encapsulating the unitary matrix ``U`` representation of the **unitary operation**.
- `accuracyCheckForUnitarity::Float64` (optional) — Numerical tolerance for unitarity validation, defaulting to `10 * eps(1.0)`.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- A `String` exception if ``U`` is not unitary within the given tolerance.

# See also
- [`QuantumCircuit`](@ref)
- [`BernoulliGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Supported **unitary operations** include
- Single-qubit gates:
  - [`createSingleQubitOperationU1`](@ref)
  - [`createSingleQubitOperationU2`](@ref)
  - [`createSingleQubitOperationU3`](@ref)
  - [`createSingleQubitOperationX`](@ref)
  - [`createSingleQubitOperationY`](@ref)
  - [`createSingleQubitOperationZ`](@ref)
  - [`createSingleQubitOperationH`](@ref)
  - [`createSingleQubitOperationId`](@ref)
  - [`createSingleQubitOperationRx`](@ref)
  - [`createSingleQubitOperationRy`](@ref)
  - [`createSingleQubitOperationRz`](@ref)
  - [`createSingleQubitOperationRotation`](@ref)
  - [`createSingleQubitOperationT`](@ref)
  - [`createSingleQubitOperationTd`](@ref)
  - [`createSingleQubitOperationS`](@ref)
  - [`createSingleQubitOperationSd`](@ref)
- Double-qubit gates:
  - [`createDoubleQubitOperationCNOT`](@ref)
  - [`createDoubleQubitOperationCNOTReverse`](@ref)
  - [`createDoubleQubitOperationSWAP`](@ref)
  - [`createDoubleQubitOperationPhase`](@ref)
- Multi-qubit gates:
  - [`createNQubitOperationProjection`](@ref)
  - [`createNQubitOperationReflection`](@ref)
  - [`createNQubitOperationExpH`](@ref)
  - [`createNQubitOperationQFT`](@ref)
  - [`createNQubitOperationIQFT`](@ref)
  - [`createNQubitOperationControlledU`](@ref)
  - [`createNQubitOperationQPE`](@ref)
  - [`compileToSingleGate`](@ref)

# Example
Add a **Bernoulli gate** for a Pauli-X operation with probability 0.3 on qubit 1 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> bernoulliGate!(qc, [1], 0.3, createSingleQubitOperationX())
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
BernoulliGate{QubiSim.var"#129#130"{UnitaryOperation}}(QubiSim.var"#129#130"{UnitaryOperation}(UnitaryOperation(ComplexF64[0.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 0.0 + 0.0im])), [1], 0.3)
```
"""
function bernoulliGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, probability::Float64, U::UnitaryOperation; accuracyCheckForUnitarity = 10*eps(1.0))
	if all(abs.(U.U'*U.U-ComplexMatrix(I,size(U.U))).<accuracyCheckForUnitarity)
		addGate!(quantumCircuit, BernoulliGate(
		    () -> U,
    		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
			probability
		))
	else
		throw("U:$U is not unitary")
	end
end

""" 
	stochasticGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits,
		unitaryOperationFactory::Function; 
		accuracyCheckForUnitarity = 10*eps(1.0)
	) -> Nothing

Adds a **Stochastic gate** to `quantumCircuit`, acting on the specified `qubits`.
The gate uses the user-supplied `unitaryOperationFactory` to generate a unitary operation dynamically during program execution.

For each shot, the stochastic gate:

1. Invokes `unitaryOperationFactory` to generate a unitary operation ``U``.
2. Checks whether ``U`` satisfies the unitarity condition ``U^†⋅U ≈ I`` within the specified numerical tolerance.
3. Applies ``U`` to the quantum state.

Under the hood, it creates a `StochasticGate` holding the factory function, and adds it to the circuit.
During compilation, this factory function is passed on to the `StochasticOperation`.
The factory function is invoked during execution, once per shot.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices (interpreted using the circuit’s `settings`) on which the gate is applied.
- `unitaryOperationFactory::Function` — A callable object (typically a closure that captures all the necessary information) that, when invoked, returns a `UnitaryOperation`.
- `accuracyCheckForUnitarity::Float64` (optional) — Numerical tolerance for unitarity validation, defaulting to `10 * eps(1.0)`.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- A `String` exception if ``U`` is not unitary within the given tolerance.

# See also
- [`QuantumCircuit`](@ref)
- [`StochasticGate`](@ref)
- [`StochasticOperation`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Example
Add a **Stochastic gate** for a U₁(λ) unitary operation with random phase rotation angle ``λ`` on qubit 1 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> stochasticGate!(qc, [1], () -> createSingleQubitOperationU1(rand()))
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
StochasticGate{var"#15#16"}(var"#15#16"(), [1], 2.220446049250313e-15)
```
"""
function stochasticGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, unitaryOperationFactory::Function; accuracyCheckForUnitarity = 10*eps(1.0))
	addGate!(quantumCircuit, StochasticGate(
	    unitaryOperationFactory,
   		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
		accuracyCheckForUnitarity
	))
end

#

""" 
	measureGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits; 
		forgetOutcome::Bool = false
	) -> Nothing

	measureGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		sigmas::Sigmas; 
		forgetOutcome::Bool = false
	) -> Nothing

	measureGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		krausOperators::KrausOperators; 
		forgetOutcome::Bool = false
	) -> Nothing

Adds a **measurement gate** to the `quantumCircuit` on the specified `qubits`, with associated measurement direction operators `sigmas` or Kraus operators `krausOperators`. The `forgetOutcome` flag determines whether the measurement outcome is observed or not revealed.

This function performs a quantum measurement based on the provided configuration:
- **Projective measurement (PVM)** is performed when `krausOperators` are not defined.
- **Generalized measurement (POVM)** is performed when `krausOperators` are defined.

Internally, the function constructs a MeasureGate using the given measurement configuration and adds it to the circuit.
During compilation, any **PVM** is automatically converted into a **POVM** by generating the corresponding Kraus operators.

# Arguments
- `quantumCircuit::QuantumCircuit` - The circuit to which the measurement gate is added.
- `qubits::Qubits` - Indices of qubits to be measured, interpreted using the circuit’s `settings`.
- `sigmas::Sigmas` (optional) - Associated measurement direction operators for each qubit, defaulting to ``σ_z``. Used for PVM type measurement.
- `krausOperators::KrausOperators` - A collection of Kraus operators that defines the generalized POVM type measurement.
- `forgetOutcome::Bool` (optional) — If `true`, the measurement result is not revealed (equivalent to taking the partial trace). If `false`, it is observed. Defaulting to `false`.
  - `false` - Simulates measurement with revealing outcome by collapsing the state to one of the possible outcomes.
  - `true` - Simulates measurement without revealing outcome by creating a mixed state over all possible outcomes.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See Also
- [`QuantumCircuit`](@ref)
- [`MeasureGateType`](@ref)
- [`MeasureGate`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`Qubits`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
- [`Settings`](@ref)

# Examples

**Example 1: Basic PVM measurement on qubits 1 and 2 along default direction ``σ_z`` with outcome observed**
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> measureGate!(qc, [1, 2])
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
MeasureGate(QubiSim.PVM, [1, 2], Matrix{ComplexF64}[], nothing, false)
```

**Example 2: PVM measurement along custom directions ``σ_x`` and ``σ_y`` with outcome observed**
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> measureGate!(qc, [1, 2], [sigmaX(), sigmaY()])
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
MeasureGate(QubiSim.PVM, [1, 2], Matrix{ComplexF64}[[0.7071067811865476 + 0.0im 0.7071067811865475 - 0.0im; 0.7071067811865475 + 0.0im -0.7071067811865476 + 0.0im], [0.7071067811865476 + 0.0im 4.329780281177466e-17 - 0.7071067811865475im; 4.329780281177466e-17 + 0.7071067811865475im -0.7071067811865476 + 0.0im]], nothing, false)
```

**Example 3: PVM measurement along direction ``σ_z`` with outcome NOT revealed to the observer**
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> measureGate!(qc, [1, 2], [sigmaZ(), sigmaZ()], forgetOutcome = true)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
MeasureGate(QubiSim.PVM, [1, 2], Matrix{ComplexF64}[[1.0 + 0.0im 0.0 - 0.0im; 0.0 + 0.0im -1.0 + 0.0im], [1.0 + 0.0im 0.0 - 0.0im; 0.0 + 0.0im -1.0 + 0.0im]], nothing, true)
```

**Example 4: POVM measurement using Kraus operators for a PVM measurement along direction ``σ_z`` with outcome observed**
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> measureGate!(qc, [1], generateKrausOperatorsForPVMMeasurement(1, [1], [sigmaZ()]), forgetOutcome = false)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
MeasureGate(QubiSim.POVM, [1], nothing, KrausOperator[KrausOperator(ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], "0"), KrausOperator(ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im], "1")], false)
```
"""
function measureGate!(quantumCircuit::QuantumCircuit, qubits::Qubits; forgetOutcome::Bool = false)
	addGate!(quantumCircuit, MeasureGate(convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits), Sigmas(), forgetOutcome))
end
function measureGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, sigmas::Sigmas; forgetOutcome::Bool = false)
	addGate!(quantumCircuit, MeasureGate(convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits), sigmas, forgetOutcome))
end
function measureGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, krausOperators::KrausOperators; forgetOutcome::Bool = false)
	addGate!(quantumCircuit, MeasureGate(convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits), krausOperators, forgetOutcome))
end

#

""" 
	quantumChannelGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		krausOperators::KrausOperators
	) -> Nothing

Adds a **quantum channel gate** (also known as a *completely positive trace-preserving map*, or CPTP map) to the `quantumCircuit` on the specified `qubits` using the provided `krausOperators`.

This function models a **superoperator** — a general type of quantum transformation that includes noise, decoherence, or collapse due to measurement.

Internally, it creates a `QuantumChannelGate` with the given Kraus operators and adds it to the circuit.

# Arguments
- `quantumCircuit::QuantumCircuit` - The circuit to which the quantum channel gate is added.
- `qubits::Qubits` - Indices of the qubits the gate acts on, interpreted using the circuit’s `settings`.
- `krausOperators::KrausOperators` - A collection of Kraus operators that defines the quantum channel.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See Also
- [`QuantumCircuit`](@ref)
- [`QuantumChannelGate`](@ref)
- [`QuantumChannelOperation`](@ref)
- [`Qubits`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
- [`Settings`](@ref)

# Example
Add a **quantum channel gate** to a single-qubit circuit using Kraus operators for a projective measurement along direction ``σ_z`` with the outcome not revealed:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> quantumChannelGate!(qc, [1], generateKrausOperatorsForPVMMeasurement(1, [1], [sigmaZ()]))
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
QuantumChannelGate([1], Matrix{ComplexF64}[[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im], [0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im]])
```
"""
function quantumChannelGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, krausOperators::KrausOperators)
	barrier!(quantumCircuit)
	addGate!(quantumCircuit, QuantumChannelGate(convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits), krausOperators))
	barrier!(quantumCircuit)
end
