############################################
# KrausOperatorsCreation

function constructKrausOperatorsForPVMMeasurement(numberOfQubits::Int, qubits::Qubits, sigmas::Sigmas)
	list = generateMeasurementIndexMap(numberOfQubits, qubits)
	transform = buildTransformMatrix(numberOfQubits, qubits, sigmas)
	labels = generateMeasurementLabels(numberOfQubits, qubits)

	# Sort list/labels for little-endian order
	ind = sortperm(list[:, 1])
	list = list[ind, :]
	labels = labels[ind]

	return constructKrausOperators(numberOfQubits, list, transform, labels)
end

function generateMeasurementIndexMap(numberOfQubits::Int, qubits::Qubits)
	numMeasured = length(qubits)
	numOutcomes = 2^numMeasured
	list = zeros(Int, numOutcomes, 2^(numberOfQubits - numMeasured))

	for m in 0:numOutcomes-1
		targetBits = reverse(digits(m, base = 2, pad = numMeasured)')
		k = 1
		for n in 0:(2^numberOfQubits - 1)
			fullBits = reverse(digits(n, base = 2, pad = numberOfQubits)')
			if isequal(fullBits[qubits]', targetBits)
				list[m+1, k] = n + 1
				k += 1
			end
		end
	end

	return list
end

function buildTransformMatrix(numberOfQubits::Int, qubits::Qubits, sigmas::Sigmas)
	transform = 1
	for j in 1:numberOfQubits
		ind = indexin(j, qubits)[1]
		V = isnothing(ind) ? Matrix{Float64}(I, 2, 2) : sigmas[ind]
		transform = tensorProduct(transform, V)
	end
	return transform
end

function generateMeasurementLabels(numberOfQubits::Int, qubits::Qubits)
	numMeasured = length(qubits)
	numOutcomes = 2^numMeasured
	labels = String[]

	for r in 1:numOutcomes
		label = []
		for k in 1:numberOfQubits
			push!(label, '*')
		end
		s = reverse(digits(r - 1, base = 2, pad = numMeasured)')
		for (k, qubit) in enumerate(qubits)
			label[qubit] = s[k]
		end
		push!(labels, join(string.(label)))
	end

	return labels
end

function constructKrausOperators(numberOfQubits::Int, list::Matrix{Int}, transform::AbstractMatrix, labels::Vector{String})
	N = 2^numberOfQubits
	krausOperators = KrausOperators()

	for m in axes(list, 1)
		Mm = zeros(N, N)
		for n in axes(list, 2)
			Mm[list[m, n], list[m, n]] = 1
		end
		Mm = transform' * Mm * transform # Mm = T'*|m><m|*T
		Em = Mm # Em = T'*|m><m|*T
		push!(krausOperators, KrausOperator(Em, labels[m]))
	end

	return krausOperators
end

""" 
    generateKrausOperatorsForPVMMeasurement(
        numberOfQubits::Int,
        qubits::Qubits,
        sigmas::Sigmas;
        zeroBasedNumbering::Bool = false,
        indexType::IndexType = circuitIndexBigEndian
    ) -> KrausOperators

Generates Kraus operators for a **projective measurement (PVM)** on a subset of qubits in a quantum circuit.

This function constructs the Kraus operators corresponding to a measurement performed on the specified `qubits` within a system of `numberOfQubits` qubits.

The measurement directions are defined by the list of operators `sigmas`. If `sigmas` is empty, the default ``σ_z`` is used for each qubit.

# Arguments
- `numberOfQubits::Int` - Total number of qubits in the system.
- `qubits::Qubits` - Indices of the qubits to be measured, interpreted using the settings for `zeroBasedNumbering` and `indexType`.
- `sigmas::Sigmas` - Associated measurement direction operators for each qubit. If `sigmas` is empty (`[]`), ``σ_z`` is used for all specified qubits.
- `zeroBasedNumbering::Bool` (optional): 
  - If `true`, qubit indices start from 0.
  - If `false`, qubit indices start from 1. Defaults to `false`.
- `indexType::IndexType` (optional) - Determines qubit indexing convention. Defaults to `circuitIndexBigEndian`.

# Returns
- `KrausOperators` — List of Kraus operators corresponding to the **projective measurement (PVM)**.

# How **PVM** is related to **POVM**
- A **POVM** is a generalized measurement defined by measurement operators ``{Mₖ}`` satisfying ``∑ₖ Mₖ^†·Mₖ = I``. Each operator is given by ``Mₖ = Eₖ^†·Eₖ`` where ``{Eₖ}`` are the Kraus operators.
- A **PVM** is a special case of the **POVM** where each Kraus operator has the form: ``Eₖ = T^†·|k⟩⟨k|·T``. Here, ``|k⟩⟨k|`` denotes the rank-1 projection operator for outcome ``k``, while ``T`` is the tensor product, taken over the measured qubits, of the corresponding **measurement direction operators** ``M(nₖ)`` as defined in `Sigmas`. Then the corresponding measurement operator becomes: ``Mₖ = Eₖ^†·Eₖ = T^†·|k⟩⟨k|·T = Eₖ``.

# See Also
- [`IndexType`](@ref)
- [`Qubits`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`measureGate!`](@ref)
- [`MeasureGate`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`quantumChannelGate!`](@ref)
- [`QuantumChannelGate`](@ref)
- [`QuantumChannelOperation`](@ref)

# Example
We prepare Kraus operators for a **projective measurement (PVM)** on a single qubit, using the measurement operator ``σ_z``, which corresponds to a measurement in the computational basis ``{|0⟩, |1⟩}``.
This results in the following Kraus operators: ``E_0 = |0⟩⟨0| = [1, 0; 0, 0]``, ``E_1 = |1⟩⟨1| = [0, 0; 0, 1]``
```julia-repl
julia> generateKrausOperatorsForPVMMeasurement(1, [1], [sigmaZ()])
2-element Vector{Any}:
 ComplexF64[1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]
 ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im]
```
"""
function generateKrausOperatorsForPVMMeasurement(numberOfQubits::Int, qubits::Qubits, sigmas::Sigmas; zeroBasedNumbering = false, indexType::IndexType = circuitIndexBigEndian)
	# assume bits conforms to Settings
	settings = Settings(zeroBasedNumbering, indexType)
	qubits = convertToByteIndex(settings, numberOfQubits, qubits)
	sigmas = sigmas == [] ? [sigmaZ() for _ in qubits] : sigmas
	krausOperators = constructKrausOperatorsForPVMMeasurement(numberOfQubits, qubits, sigmas)
	return krausOperators
end

""" 
    generateKrausOperatorsForDepolarizingChannel(
		p::Float64
	) -> KrausOperators

Constructs the Kraus operators for a single-qubit **depolarizing quantum channel**.  
This channel applies a random Pauli error ``σ_x``, ``σ_y`` or ``σ_z`` with total probability ``p``, and leaves the state unchanged with probability ``1 - p``.
The depolarizing channel is defined as:

```math
ρ → (1-p)·ρ + (p/3)·(σ_x·ρ·σ_x+σ_y·ρ·σ_y+σ_z·ρ·σ_z)
```

where ``ρ`` is the input density matrix, ``σ_x``, ``σ_y`` and ``σ_z`` are the Pauli matrices and ``p`` is the error probability.
The corresponding Kraus operators are: ``M_I = √(1-p)·I``, ``M_X = √(p/3)·σ_x``, ``M_Y = √(p/3)·σ_y`` and ``M_Z = √(p/3)·σ_z``.

# Effects
- Both populations driven toward ``½``
- Off-diagonals shrink by ``1 - 4p/3``

# Arguments
- `p::Float64` - The probability that a qubit undergoes an error.

# Returns
- `KrausOperators` — List of Kraus operators ``M_I, M_X, M_Y, M_Z`` corresponding to the **depolarizing quantum channel**.

# See Also
- [`KrausOperators`](@ref)
- [`quantumChannelGate!`](@ref)
- [`QuantumChannelGate`](@ref)
- [`QuantumChannelOperation`](@ref)

# Example
Applying a depolarizing channel with error probability `p = 0.75` on a single qubit in superposition:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> quantumChannelGate!(qc, [1], generateKrausOperatorsForDepolarizingChannel(0.75))
julia> qp = compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
julia> iqs = createInitialQubitState(density, [[1., [pi/2 0]]], blochRepresentation=true)
julia> qo = runQuantumProgram(qp, iqs, 1)

julia> # First state: pure superposition (|0⟩ + |1⟩)/√2
julia> round.(real.(qo.output[1].rho), digits=3)
[0.5 0.5; 0.5 0.5]
julia> qo.output[1].entropy
-9.056320547599769e-31
julia> qo.output[1].purityValue
1.0
julia> qo.output[1].purityNorm
1.0

julia> # Second state: completely mixed state (Id/2) due to p = 0.75
julia> round.(real.(qo.output[2].rho), digits=3)
[0.5 0.0; 0.0 0.5]
julia> qo.output[2].entropy
0.6931471805599453
julia> qo.output[2].purityValue
0.5
julia> qo.output[2].purityNorm
0.0
```
"""
function generateKrausOperatorsForDepolarizingChannel(p::Float64)
	krausOperators = KrausOperators()
	push!(krausOperators, KrausOperator(sqrt(1-p)*ComplexMatrix([1. 0; 0. 1.]), "MI"))
	push!(krausOperators, KrausOperator(sqrt(p/3)*pauliX(), "MX"))
	push!(krausOperators, KrausOperator(sqrt(p/3)*pauliY(), "MY"))
	push!(krausOperators, KrausOperator(sqrt(p/3)*pauliZ(), "MZ"))
	return krausOperators
end

""" 
    generateKrausOperatorsForPhaseDampingChannel(
        p::Float64
    ) -> KrausOperators

Constructs the Kraus operators for a single-qubit phase damping quantum channel.
This channel models the loss of quantum coherence (dephasing) without energy dissipation.
With probability ``p``, the off-diagonal elements of the density matrix decay toward zero, while populations remain unchanged.
The phase damping channel acts on a density matrix ``ρ`` as:

```math
ρ = [ρ_{00}, ρ_{01}; ρ_{10}, ρ_{11}] → [ρ_{00}, (1−p)⋅ρ_{01}; (1−p)⋅ρ_{10}, ρ_{11}]
```

where
- the diagonal elements (populations) are unaffected,,
- the off-diagonal elements (coherences) decay linearly with factor ``1-p``.

The corresponding Kraus operators are: ``M_0 = √(1-p)·I`` (induces no decoherence), ``M_1 = √p·|0⟩⟨0|`` (induces decoherence) and ``M_2 = √p·|1⟩⟨1|`` (induces decoherence).

# Effects
- Populations unchanged
- Off-diagonals shrink linearly by ``1-p``

# Arguments
- `p::Float64` - The probability that coherence between ``|0⟩`` and ``|1⟩`` is lost.

# Returns
- `KrausOperators` — List of Kraus operators ``M_0, M_1, M_2`` corresponding to the **phase damping quantum channel**.

# See Also
- [`KrausOperators`](@ref)
- [`quantumChannelGate!`](@ref)
- [`QuantumChannelGate`](@ref)
- [`QuantumChannelOperation`](@ref)

# Example
Applying a phase damping channel with p = 0.5 on a qubit in superposition:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> quantumChannelGate!(qc, [1], generateKrausOperatorsForPhaseDampingChannel(0.5))
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(density, [[1., [pi/2 0]]], blochRepresentation=true)
julia> qo = runQuantumProgram(qp, iqs, 1)

julia> round.(real.(qo.output[2].rho), digits=3)
[0.5 0.25; 0.25 0.5]
```
"""
function generateKrausOperatorsForPhaseDampingChannel(p::Float64)
	krausOperators = KrausOperators()
	push!(krausOperators, KrausOperator(sqrt(1-p)*ComplexMatrix([1. 0; 0. 1.]), "M0"))
	push!(krausOperators, KrausOperator(sqrt(p)*ComplexMatrix([1. 0; 0. 0.]), "M1"))
	push!(krausOperators, KrausOperator(sqrt(p)*ComplexMatrix([0. 0; 0. 1.]), "M2"))
	return krausOperators
end

""" 
    generateKrausOperatorsForAmplitudeDampingChannel(
        p::Float64
    ) -> KrausOperators

Constructs the Kraus operators for a single-qubit amplitude damping quantum channel.
This channel models energy dissipation, where the excited state ``|1⟩`` decays to the ground state ``|0⟩`` with probability ``p``.
The amplitude damping channel acts on a density matrix ``ρ`` as:

```math
ρ = [ρ_{00}, ρ_{01}; ρ_{10}, ρ_{11}] → [ρ_{00}+p⋅ρ_{11}, √(1−p)⋅ρ_{01}; √(1−p)⋅ρ_{10}, (1−p)⋅ρ_{11}]
```

where
- the population of ``|1⟩`` decays into ``|0⟩`` with probability ``p``,
- off-diagonal coherence terms shrink by ``√(1-p)``,
- population of ``|0⟩`` increases accordingly, preserving trace.

The corresponding Kraus operators are: ``M_0 = |0⟩⟨0| + √(1-p)⋅|1⟩⟨1|`` (induces no decay) and ``M_1 = √p⋅|0⟩⟨1|`` (induces decay ``|1⟩ → |0⟩``).

# Effects
- ``ρ_{11} → (1-p)ρ_{11}``, population leaks to ``ρ_{00}``
- Off-diagonals shrink by ``√(1-p)``

# Arguments
- `p::Float64` - The probability of spontaneous decay ``|1⟩ → |0⟩``.

# Returns
- `KrausOperators` — List of Kraus operators ``M_0, M_1`` corresponding to the **amplitude damping quantum channel**.

# See Also
- [`KrausOperators`](@ref)
- [`quantumChannelGate!`](@ref)
- [`QuantumChannelGate`](@ref)
- [`QuantumChannelOperation`](@ref)

# Example
Applying an amplitude damping channel with p = 0.3 on a qubit:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> quantumChannelGate!(qc, [1], generateKrausOperatorsForAmplitudeDampingChannel(0.3))
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(density, [[1., [pi/2 0]]], blochRepresentation=true)
julia> qo = runQuantumProgram(qp, iqs, 1)

julia> round.(real.(qo.output[2].rho), digits=3)
[0.65 0.418; 0.418 0.35]
```
"""
function generateKrausOperatorsForAmplitudeDampingChannel(p::Float64)
	krausOperators = KrausOperators()
	push!(krausOperators, KrausOperator(ComplexMatrix([1. 0; 0. sqrt(1-p)]), "M0"))
	push!(krausOperators, KrausOperator(ComplexMatrix([0. sqrt(p); 0. 0.]), "M1"))
	return krausOperators
end

""" 
    generateKrausOperatorsForBitFlipChannel(
        p::Float64
    ) -> KrausOperators

Constructs the Kraus operators for a single-qubit bit-flip quantum channel.
The bit-flip channel models a probabilistic application of the Pauli-X operation.
With probability `p`, the qubit state is flipped, ``|0⟩ ↔ |1⟩``, while with probability `1-p` the state is left unchanged.

The channel acts on a density matrix ``ρ`` as:

```math
ρ = [ρ_{00}, ρ_{01}; ρ_{10}, ρ_{11}] → (1-p)⋅ρ + p⋅(X⋅ρ⋅X) = [(1-p)⋅ρ_{00}+p⋅ρ_{11}, (1−p)⋅ρ_{01}+p⋅ρ_{10}; (1−p)⋅ρ_{10}+p⋅ρ_{01}, (1−p)⋅ρ_{11}+p⋅ρ_{00}]
```

where ``X`` is the Pauli-X operator.

The corresponding Kraus operators are: ``M_0 = √(1-p)⋅I`` (induces no bit-flip) and ``M_1 = √p⋅X`` (induces bit-flip ``|1⟩ ↔ |0⟩``).

# Effects
- With probability ``1-p``, the qubit is unchanged.
- With probability ``p``, the qubit undergoes a bit-flip ``|0⟩ ↔ |1⟩``.

# Arguments
- `p::Float64` - The probability of spontaneous bit-flip ``|0⟩ ↔ |1⟩``.

# Returns
- `KrausOperators` — List of Kraus operators ``M_0, M_1`` corresponding to the **bit-flip quantum channel**.

# See Also
- [`KrausOperators`](@ref)
- [`quantumChannelGate!`](@ref)
- [`QuantumChannelGate`](@ref)
- [`QuantumChannelOperation`](@ref)

# Example
Applying a bit-flip quantum channel with p = 0.3 on a qubit:
```julia-repl
julia> qc = createQuantumCircuit(1)
julia> quantumChannelGate!(qc, [1], generateKrausOperatorsForBitFlipChannel(0.3))
julia> qp = compileQuantumCircuit(qc)
julia> iqs = createInitialQubitState(density, [[1., [0. 0.]]], blochRepresentation=true)
julia> qo = runQuantumProgram(qp, iqs, 1)

julia> round.(real.(qo.output[2].rho), digits=3)
[0.7 0.0; 0.0 0.3]
```
"""
function generateKrausOperatorsForBitFlipChannel(p::Float64)
	krausOperators = KrausOperators()
	push!(krausOperators, KrausOperator(ComplexMatrix([sqrt(1-p) 0; 0. sqrt(1-p)]), "M0"))
	push!(krausOperators, KrausOperator(ComplexMatrix([0. sqrt(p); sqrt(p) 0.]), "M1"))
	return krausOperators
end

