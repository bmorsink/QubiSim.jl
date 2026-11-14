############################################
# NQubitUnitaryGates

""" 
	projectionGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		v::AbstractVector{<:Number}
	) -> Nothing

Adds a **projection gate** to the `quantumCircuit` on the specified `qubits` that projects onto an arbitrary quantum state vector ``|v⟩``.

The **projection gate** performs a non-unitary operation that projects the system onto the quantum state defined by ``|v⟩``, a complex-valued state vector.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **projection unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createNQubitOperationProjection` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices of length ``M`` (interpreted using the circuit’s `settings`) on which the gate is applied.
- `v::AbstractVector{<:Number}` - A complex-valued vector ``|v⟩`` of length ``2^M`` representing the target quantum state for projection.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- `ErrorException` - If the length of ``|v⟩`` does not match ``2^M``.

# Notes
- The vector ``|v⟩`` should be normalized, i.e., ``‖|v⟩‖ ≈ 1`` to ensure the projector behaves correctly.

# See also
- [`createNQubitOperationProjection`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Example
Add a **projection gate** onto the state vector ``|v⟩ = [1; 0; 0; 0]`` to a 2-qubit quantum circuit on qubits 1 and 2:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> projectionGate!(qc, [1,2], Vector{ComplexF64}([1; 0; 0; 0]))
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#31#32"{Vector{ComplexF64}}}(QubiSim.var"#31#32"{Vector{ComplexF64}}(ComplexF64[1.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im, 0.0 + 0.0im]), [1, 2], "Projection")
```
"""
function projectionGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, v::AbstractVector{<:Number})
	if log2(length(v))==length(qubits)
		addGate!(quantumCircuit, UnitaryGate(
    	    () -> createNQubitOperationProjection(v),
        	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
			"Projection"
    	))
	else
		throw("The length of vector v does not match the specified amount of qubits for the projection gate.")
	end
end

""" 
    createNQubitOperationProjection(
		v::AbstractVector{<:Number}
	) -> UnitaryOperation

Creates a structure encapsulating the matrix representing the **projection operation**.

The **projection operation** projects the system onto the quantum state defined by ``|v⟩``, a complex-valued state vector.

This **projection operation** is not unitary but is returned as a `UnitaryOperation` object for integration into the circuit model.

# Arguments
- `v::AbstractVector{<:Number}` - A complex-valued vector ``|v⟩`` having length equal to a power of two, representing the target quantum state for projection.

# Returns
- `UnitaryOperation` — A structure encapsulating the rank-1 projection matrix ``P = |v⟩⟨v|`` representation of the **projection operation**.

# Notes
- The vector ``|v⟩`` should be normalized, i.e., ``‖|v⟩‖ ≈ 1`` to ensure the projector behaves correctly.
- This operation is not strictly unitary; it is packaged in a UnitaryOperation for compatibility purposes.

# See also
- [`projectionGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the matrix representation of the **projection operation** onto vector ``|v⟩ = [1; 0; 0; 0]``:
```julia-repl
julia> unitaryOperation = createNQubitOperationProjection([1; 0; 0; 0])
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
```
"""
function createNQubitOperationProjection(v::AbstractVector{<:Number})
	UnitaryOperation(v*v')
end

#

""" 
	reflectionGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		v::AbstractVector{<:Number}, 
		theta::Float64
	) -> Nothing

Adds a **reflection gate** to the `quantumCircuit` on the specified `qubit`.

The **reflection gate** performs a unitary operation that reflects the quantum state about the hyperplane orthogonal to an arbitrary quantum state vector ``|v⟩``, scaled by a complex phase factor determined by the angle ``θ``.

This operation implements the generalized reflection operator: ``R(|v⟩, θ) = I - (1 - exp(−i⋅θ))⋅|v⟩⟨v|`` where ``|v⟩`` is a normalized quantum state vector.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **reflection unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createNQubitOperationReflection` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices of length ``M`` (interpreted using the circuit’s `settings`) on which the gate is applied.
- `v::AbstractVector{<:Number}` - A complex-valued vector ``|v⟩`` of length ``2^M`` representing the quantum state vector to reflect about.
- `theta::Float64` - The scaling phase angle ``θ`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- `ErrorException` - If the length of ``|v⟩`` is not equal to ``2^M``.

# Notes
- The vector ``|v⟩`` **must be** normalized, i.e., ``‖|v⟩‖ ≈ 1``.

# See also
- [`createNQubitOperationReflection`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Example
Add a **reflection gate** around vector ``|v⟩ = [1; 0; 0; 0]`` with scaling angle ``θ = π`` to a 2-qubit quantum circuit on qubits 1 and 2:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> reflectionGate!(qc, [1,2], [1;0;0;0], 1.0*pi)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#33#34"{Vector{Int64}, Float64}}(QubiSim.var"#33#34"{Vector{Int64}, Float64}([1, 0, 0, 0], 3.141592653589793), [1, 2], "Reflection")
```
"""
function reflectionGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, v::AbstractVector{<:Number}, theta::Float64)
	if log2(length(v))==length(qubits)
		addGate!(quantumCircuit, UnitaryGate(
    	    () -> createNQubitOperationReflection(v, theta),
        	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
			"Reflection"
    	))
	else
		throw("The length of vector v does not match the specified amount of qubits for the reflection gate.")
	end
end

""" 
    createNQubitOperationReflection(
		v::AbstractVector{<:Number}, 
		theta::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **reflection unitary operation**.

The **reflection unitary operation** reflects the quantum state about the hyperplane orthogonal to an arbitrary quantum state vector ``|v⟩``, scaled by a complex phase factor determined by the angle ``θ``.

This operation implements the generalized reflection operator: ``R(|v⟩, θ) = I - (1 - exp(−i⋅θ))⋅|v⟩⟨v|`` where ``|v⟩`` is a normalized quantum state vector.

# Arguments
- `v::AbstractVector{<:Number}` - A complex-valued vector ``|v⟩`` having length equal to a power of two, representing the quantum state vector to reflect about.
- `theta::Float64` - The scaling phase angle ``θ`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **reflection unitary operation**.

# Notes
- The vector ``|v⟩`` **must be** normalized, i.e., ``‖|v⟩‖ ≈ 1``.

# See also
- [`reflectionGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **reflection unitary operation** about vector ``|v⟩ = [1; 0; 0; 0]`` with scaling angle ``θ = π``:
```julia-repl
julia> unitaryOperation = createNQubitOperationReflection(ComplexVector([1; 0; 0; 0]), 1.0*pi)
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 -1.0-0.0im  0.0-0.0im  0.0-0.0im  0.0-0.0im
  0.0-0.0im  1.0-0.0im  0.0-0.0im  0.0-0.0im
  0.0-0.0im  0.0-0.0im  1.0-0.0im  0.0-0.0im
  0.0-0.0im  0.0-0.0im  0.0-0.0im  1.0-0.0im
```
"""
function createNQubitOperationReflection(v::AbstractVector{<:Number}, theta::Float64)
	UnitaryOperation(ComplexMatrix(I, length(v), length(v))-(1-exp(-im*theta))*v*v')
end

#

""" 
	unitaryUGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		U::UnitaryOperation; 
		accuracyCheckForUnitarity = 10*eps(1.0)
	) -> Nothing

Adds a **unitary gate** to the `quantumCircuit` on the specified `qubits` with a user-supplied **unitary operation** ``U``.

The **unitary gate** performs the user-supplied **unitary operation** ``U``.

The function performs an internal unitarity check to ensure that the user-supplied **unitary operation** ``U`` satisfies the condition ``U^†⋅U ≈ I`` within a specified numerical tolerance.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **unitary operation**, and adds 
it to the circuit. During compilation, this factory function simply returns the user-suppied **unitary operation** `U` as the produced unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices (interpreted using the circuit’s `settings`) on which the gate is applied.
- `U::UnitaryOperation` — A structure encapsulating the unitary matrix ``U`` representation of the **unitary operation**.
- `accuracyCheckForUnitarity::Float64` (optional) — Numerical tolerance for unitarity validation, defaulting to `10 * eps(1.0)`.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- A `String` exception if ``U`` is not unitary within the given tolerance.

# See also
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
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
Add a custom **unitary gate** with a **unitary operation** for a Z-gate on qubit 1 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> unitaryUGate!(qc, [1], createSingleQubitOperationZ())
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#50#51"{UnitaryOperation}}(QubiSim.var"#50#51"{UnitaryOperation}(UnitaryOperation(ComplexF64[1.0 + 0.0im 0.0 - 0.0im; 0.0 + 0.0im -1.0 + 1.2246467991473532e-16im])), [1], "Unitary")
```
"""
function unitaryUGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, U::UnitaryOperation; accuracyCheckForUnitarity = 10*eps(1.0))
	if all(abs.(U.U'*U.U-ComplexMatrix(I,size(U.U))).<accuracyCheckForUnitarity)
		addGate!(quantumCircuit, UnitaryGate(
		    () -> U,
    		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
			"Unitary"
		))
	else
		throw("U:$U is not unitary")
	end
end

#

""" 
	expHGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits, 
		H::ComplexMatrix, 
		t_over_h::Float64; 
		accuracyCheckForHermiticity = 10*eps(1.0))
	) -> Nothing

Adds a **time evolution gate** to the `quantumCircuit` on the specified `qubits` for a Hamiltonian ``H`` and scaled time ``t/ℏ``.

The **time evolution gate** performs a unitary operation that applies the time evolution:
```math
exp(-i·H·t/ℏ)
```

where ``H`` is the system's Hamiltonian, and ``t/ℏ`` is the ratio of evolution time ``t`` over Planck’s constant ``ℏ``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **time evolution unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createNQubitOperationExpH` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices of length ``M`` (interpreted using the circuit’s `settings`) on which the gate is applied.
- `H::ComplexMatrix` — The Hamiltonian matrix ``H`` of size ``2^M × 2^M`` governing the evolution, must be Hermitian.
- `t_over_h::Float64` — Time evolution parameter ``t/ℏ``, expressed in radians.
- `accuracyCheckForHermiticity::Float64` (optional) — Numerical tolerance for Hermiticity validation of ``H``, defaulting to `10 * eps(1.0)`.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- A `String` exception if the dimensions of ``H`` are incompatible with the number of target qubits or fails the Hermiticity check (i.e., ``H^† ≠ H`` within the specified tolerance).

# See also
- [`createNQubitOperationExpH`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Example
Add a **time evolution gate** for the Hamiltonian:

    H = [ 0.460075592255305    -1.110720734539591;
         -1.110720734539591     2.681517061334488 ]

with evolution time ``t/ℏ = 1.0`` on qubit 1 to a 2-qubit circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> expHGate!(qc, [1], ComplexMatrix([0.460075592255305 -1.110720734539591;-1.110720734539591 2.681517061334488]), 1.0)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#53#54"{Matrix{ComplexF64}, Float64}}(QubiSim.var"#53#54"{Matrix{ComplexF64}, Float64}(ComplexF64[0.460075592255305 + 0.0im -1.110720734539591 + 0.0im; -1.110720734539591 + 0.0im 2.681517061334488 + 0.0im], 1.0), [1], "expH")
```
"""
function expHGate!(quantumCircuit::QuantumCircuit, qubits::Qubits, H::ComplexMatrix, t_over_h::Float64; accuracyCheckForHermiticity = 10*eps(1.0))
	if size(H) == (2^length(qubits), 2^length(qubits))
		if all(abs.(H'-H).<accuracyCheckForHermiticity)
			addGate!(quantumCircuit, UnitaryGate(
	    		() -> createNQubitOperationExpH(H, t_over_h),
    			convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
				"expH"
			))
		else
			throw("H is not Hermitian")
		end
	else
		throw("Dimensions of H are incompatible with the number of target qubits")
	end
end

""" 
	createNQubitOperationExpH(
		H::ComplexMatrix, 
		t_over_h::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **time evolution unitary operation**.

The **time evolution unitary operation** applies the time evolution: 
```math
exp(-i·H·t/ℏ)
```

where ``H`` is the system's Hamiltonian, and ``t/ℏ`` is the ratio of evolution time ``t`` over Planck’s constant ``ℏ``.

# Arguments
- `H::ComplexMatrix` — The Hamiltonian matrix ``H`` governing the evolution, must be Hermitian.
- `t_over_h::Float64` — Time evolution parameter ``t/ℏ``, expressed in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **time evolution unitary operation**.

# Notes
- This function supports any ``M``-qubit system as long as the size of ``H`` is ``2^M × 2^M``.
- ``H`` must be Hermitian to ensure the resulting operation is unitary.

# See also
- [`expHGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **time evolution unitary operation** for a single-qubit system with the Hamiltonian:

    H = [ 0.460075592255305    -1.110720734539591;
         -1.110720734539591     2.681517061334488 ]

and evolution time ``t/ℏ = 1.0``:
```julia-repl
julia> unitaryOperation = createNQubitOperationExpH(ComplexMatrix([0.460075592255305  -1.110720734539591;-1.110720734539591   2.681517061334488]), 1.0)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
 0.707107+0.0im   0.707107+0.0im
 0.707107+0.0im  -0.707107+0.0im
```
"""
function createNQubitOperationExpH(H::ComplexMatrix,t_over_h::Float64)
	UnitaryOperation(exp(-im * H * t_over_h))
end

#

""" 
	qftGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits
	) -> Nothing

Adds a **Quantum Fourier Transform (QFT) gate** to the `quantumCircuit` on the specified `qubits`.

The **QFT gate** performs a unitary operation that transforms the basis state ``|x⟩`` into: 
```math
QFT|x⟩ = 1/√N·∑_{y=0}^{N−1} exp(2·π·i·x·y/N)·|y⟩
```

where ``N = 2^M`` with ``M`` equal to the length of `qubits`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **QFT unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createNQubitOperationQFT` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices of length ``M`` (interpreted using the circuit’s `settings`) on which the gate is applied.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createNQubitOperationQFT`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Example
Add a **QFT gate** on qubits 1 and 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> qftGate!(qc, [1,2])
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#55#56"{Vector{Int64}}}(QubiSim.var"#55#56"{Vector{Int64}}([1, 2]), [1, 2], "QFT")
```
"""
function qftGate!(quantumCircuit::QuantumCircuit, qubits::Qubits)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createNQubitOperationQFT(length(qubits)),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
		"QFT"
	))
end

""" 
	createNQubitOperationQFT(
		numberOfQubits::Int
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Quantum Fourier Transform (QFT) unitary operation** on `numberOfQubits` qubits.

The **QFT unitary operation** transforms the basis state ``|x⟩`` into: 
```math
QFT|x⟩ = 1/√N·∑_{y=0}^{N−1} exp(2·π·i·x·y/N)·|y⟩
```

where ``N = 2^M`` with ``M`` equal to `numberOfQubits`.

Under the hood, it builds the full quantum circuit implementing the **QFT** and 
compiles it into a single unitary matrix representing the **QFT unitary operation**.

# Arguments
- `numberOfQubits::Int` — The number of qubits ``M`` over which the **QFT** is applied.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **QFT unitary operation**.

# See also
- [`qftGate!`](@ref)
- [`createNQubitQFTQuantumCircuit`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **QFT unitary operation** acting on 2 qubits:
```julia-repl
julia> unitaryOperation = createNQubitOperationQFT(2)
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 0.5+0.0im   0.5+0.0im   0.5+0.0im   0.5+0.0im
 0.5+0.0im   0.0+0.5im  -0.5+0.0im   0.0-0.5im
 0.5+0.0im  -0.5+0.0im   0.5+0.0im  -0.5+0.0im
 0.5+0.0im   0.0-0.5im  -0.5+0.0im   0.0+0.5im
```
"""
function createNQubitOperationQFT(numberOfQubits::Int)
	quantumCircuit = createNQubitQFTQuantumCircuit(numberOfQubits)
	return compileToSingleGate(quantumCircuit)
end

"""
    createNQubitQFTQuantumCircuit(
		numberOfQubits::Int
	) -> QuantumCircuit{Z, I}

Builds the quantum circuit implementing the **Quantum Fourier Transform (QFT)** on `numberOfQubits` qubits.

The **QFT** transforms the basis state ``|x⟩`` into: 
```math
QFT|x⟩ = 1/√N·∑_{y=0}^{N−1} exp(2·π·i·x·y/N)·|y⟩
```

where ``N = 2^M`` with ``M`` equal to `numberOfQubits`.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `numberOfQubits::Int` — The number of qubits ``M`` over which the **QFT** is applied.

# Returns
- `QuantumCircuit{Z, I}` — A quantum circuit implementing the **QFT**

# See also
- [`createNQubitOperationQFT`](@ref)
- [`QuantumCircuit`](@ref)

# Example
Create the quantum circuit of a 2-qubit **QFT**:
```julia-repl
julia> quantumCircuit = createNQubitQFTQuantumCircuit(2)
julia> quantumCircuit.circuit
4-element Vector{Step}:
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [1], "H")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#47#48"{Float64}}(QubiSim.var"#47#48"{Float64}(1.5707963267948966), [2, 1], "Phase")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [2], "H")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#45#46"}(QubiSim.var"#45#46"(), [1, 2], "SWAP")])
```
"""
function createNQubitQFTQuantumCircuit(numberOfQubits::Int)
	quantumCircuit = createQuantumCircuit(numberOfQubits)
	
	for m in 1:numberOfQubits
        hGate!(quantumCircuit, m)
		for n in 1:(numberOfQubits - m)
            phaseGate!(quantumCircuit, n + m, m, pi / (2 ^ n))
		end
	end
	for m in 1:Int(floor(numberOfQubits / 2))
        swapGate!(quantumCircuit, m, numberOfQubits - m + 1)
	end
	return quantumCircuit
end

#

""" 
	iqftGate!(
		quantumCircuit::QuantumCircuit, 
		qubits::Qubits
	) -> Nothing

Adds an **Inverse Quantum Fourier Transform (IQFT) gate** to the `quantumCircuit` on the specified `qubits`.

The **IQFT gate** performs a unitary operation that transforms the basis state ``|x⟩`` into: 
```math
IQFT|x⟩ = 1/√N·∑_{y=0}^{N−1} exp(-2·π·i·x·y/N)·|y⟩
```

where ``N = 2^M`` with ``M`` equal to the length of `qubits`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **IQFT unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createNQubitOperationIQFT` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubits::Qubits` — A collection of qubit indices of length ``M`` (interpreted using the circuit’s `settings`) on which the gate is applied.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createNQubitOperationIQFT`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Example
Add an **IQFT gate** on qubits 1 and 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc=createQuantumCircuit(2)
julia> iqftGate!(qc, [1,2])
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#57#58"{Vector{Int64}}}(QubiSim.var"#57#58"{Vector{Int64}}([1, 2]), [1, 2], "IQFT")
```
"""
function iqftGate!(quantumCircuit::QuantumCircuit, qubits::Qubits)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createNQubitOperationIQFT(length(qubits)),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, qubits),
		"IQFT"
	))
end

""" 
    createNQubitOperationIQFT(
		numberOfQubits::Int
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Inverse Quantum Fourier Transform (IQFT) unitary operation** on `numberOfQubits` qubits.

The **IQFT unitary operation** transforms the basis state ``|x⟩`` into:
```math
IQFT|x⟩ = 1/√N·∑_{y=0}^{N−1} exp(-2·π·i·x·y/N)·|y⟩
```

where ``N = 2^M`` with ``M`` equal to `numberOfQubits`.

Under the hood, it builds the full quantum circuit implementing the **IQFT** and 
compiles it into a single unitary matrix representing the **IQFT unitary operation**.

# Arguments
- `numberOfQubits::Int` — The number of qubits ``M`` over which the **IQFT** is applied.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **IQFT unitary operation**.

# See also
- [`iqftGate!`](@ref)
- [`createNQubitIQFTQuantumCircuit`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **IQFT unitary operation** acting on 2 qubits:
```julia-repl
julia> unitaryOperation=createNQubitOperationIQFT(2)
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 0.5+0.0im   0.5+0.0im   0.5+0.0im   0.5+0.0im
 0.5+0.0im   0.0-0.5im  -0.5+0.0im   0.0+0.5im
 0.5+0.0im  -0.5+0.0im   0.5+0.0im  -0.5+0.0im
 0.5+0.0im   0.0+0.5im  -0.5+0.0im   0.0-0.5im
```
"""
function createNQubitOperationIQFT(numberOfQubits::Int)
	quantumCircuit = createNQubitIQFTQuantumCircuit(numberOfQubits)
	return compileToSingleGate(quantumCircuit)
end

""" 
    createNQubitIQFTQuantumCircuit(
		numberOfQubits::Int
	) -> QuantumCircuit{Z, I}

Builds the quantum circuit implementing the **Inverse Quantum Fourier Transform (IQFT)** on `numberOfQubits` qubits.

The **IQFT** transforms the basis state ``|x⟩`` into: 
```math
IQFT|x⟩ = 1/√N·∑_{y=0}^{N−1} exp(-2·π·i·x·y/N)·|y⟩
```

where ``N = 2^M`` with ``M`` equal to `numberOfQubits`.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `numberOfQubits::Int` — The number of qubits ``M`` over which the **IQFT** is applied.

# Returns
- `QuantumCircuit{Z, I}` — A quantum circuit implementing the **IQFT**

# See also
- [`createNQubitOperationIQFT`](@ref)
- [`QuantumCircuit`](@ref)

# Example
Create the quantum circuit of a 2-qubit **IQFT**:
```julia-repl
julia> quantumCircuit = createNQubitIQFTQuantumCircuit(2)
julia> quantumCircuit.circuit
4-element Vector{Step}:
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#45#46"}(QubiSim.var"#45#46"(), [1, 2], "SWAP")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [2], "H")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#47#48"{Float64}}(QubiSim.var"#47#48"{Float64}(-1.5707963267948966), [2, 1], "Phase")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [1], "H")])
```
"""
function createNQubitIQFTQuantumCircuit(numberOfQubits::Int)
	quantumCircuit = createQuantumCircuit(numberOfQubits)

	for m in 1:Int(floor(numberOfQubits / 2))
       swapGate!(quantumCircuit, m, numberOfQubits - m + 1)
	end
	for m in numberOfQubits:-1:1
		for n in 1:(numberOfQubits - m)
            phaseGate!(quantumCircuit, n + m, m, -pi / (2 ^ n))
		end
        hGate!(quantumCircuit, m)
	end
	return quantumCircuit
end

#

""" 
	controlledUGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubits::Qubits, 
		targetQubits::Qubits, 
		UTarget::UnitaryOperation
	) -> Nothing

Adds a **controlled unitary gate** to the `quantumCircuit` on the specified qubits `controlQubits` and `targetQubits` and unitary operation `UTarget`.

The **controlled unitary gate** performs a unitary operation that applies the unitary operation `UTarget` on the `targetQubits` if the `controlQubits` are **all** in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubits`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createNQubitOperationControlledU` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubits::Qubits` — Control qubit indices (interpreted via the circuit's `settings`).
- `targetQubits::Qubits` — Target qubit indices (interpreted via the circuit's `settings`).
- `UTarget::UnitaryOperation` — The unitary operation to be applied.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- A `String` exception if the dimensions of `UTarget` are incompatible with the number of `targetQubits`.

# See also
- [`createNQubitOperationControlledU`](@ref)
- [`UnitaryOperation`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`Qubits`](@ref)

# Example
Add a **controlled Pauli-Z gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledUGate!(qc, [1], [2], createSingleQubitOperationZ())
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}}(QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}([1], [2], UnitaryOperation(ComplexF64[1.0 + 0.0im 0.0 - 0.0im; 0.0 + 0.0im -1.0 + 1.2246467991473532e-16im])), [1, 2], "Controlled-Unitary")
```
"""
function controlledUGate!(quantumCircuit::QuantumCircuit, controlQubits::Qubits, targetQubits::Qubits, UTarget::UnitaryOperation)
	if size(UTarget.U) == (2^length(targetQubits), 2^length(targetQubits))
		addGate!(quantumCircuit, UnitaryGate(
			() -> createNQubitOperationControlledU(length(controlQubits), length(targetQubits), UTarget),
			convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, union(controlQubits, targetQubits)),
			"Controlled-Unitary"
		))
	else
		throw("Dimensions of UTarget are incompatible with the number of target qubits")
	end
end

""" 
    createNQubitOperationControlledU(
		numberOfControlQubits::Int, 
		numberOfTargetQubits::Int, 
		UTarget::UnitaryOperation
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled unitary operation**.

The **controlled unitary operation** applies the unitary operation `UTarget` on the target qubits if the control qubits are **all** in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubits.

# Arguments
- `numberOfControlQubits::Int` — Number of control qubits.
- `numberOfTargetQubits::Int` — Number of target qubits.
- `UTarget::UnitaryOperation` — The unitary operation to be applied.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled unitary operation**.

# See also
- [`controlledUGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **controlled Pauli-Z unitary operation** with one control and one target qubit:
```julia-repl
julia> unitaryOperation = createNQubitOperationControlledU(1, 1, createSingleQubitOperationZ())
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im   0.0-0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  -1.0+0.0im
```
"""
function createNQubitOperationControlledU(numberOfControlQubits::Int, numberOfTargetQubits::Int, UTarget::UnitaryOperation)
	totalNumberOfQubits = numberOfControlQubits+numberOfTargetQubits
	controlledU = createNQubitOperationId(totalNumberOfQubits)
    controlledU.U[2^totalNumberOfQubits-size(UTarget.U,1)+1:end, 
                  2^totalNumberOfQubits-size(UTarget.U,2)+1:end] = UTarget.U

	return controlledU
end

#

""" 
	qpeGate!(
		quantumCircuit::QuantumCircuit, 
		estimationQubits::Qubits,
		targetQubits::Qubits,
		U::UnitaryOperation;
		accuracyCheckForUnitarity = 10 * eps(1.0)
	) -> Nothing

Adds a **Quantum Phase Estimation (QPE) gate** to `quantumCircuit`, acting on the specified `estimationQubits` and `targetQubits`, for the given unitary operation ``U``.

The **QPE gate** estimates the eigenphase(s) of the unitary operation ``U``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **QPE unitary operation**, and adds
it to the circuit. During compilation, this factory function calls `createNQubitOperationQPE` to produce the unitary matrix.

## Conceptual Overview
- The **QPE algorithm** estimates the eigenphase ``φ`` of a unitary operation ``U`` such that ``U|ψ⟩ = exp(2·π·i·φ)·|ψ⟩``. It operates on two quantum registers:
  - **Estimation register (`estimationQubits`)** — A sequence of qubits that stores the binary representation of the estimated eigenphase ``φ``.
  - **Target register (`targetQubits`)** — The qubits on which the unitary operation ``U`` acts. This register must be initialized in an eigenstate ``|ψ⟩`` of ``U`` (or a superposition of eigenstates) to yield meaningful phase information.
- The **QPE algorithm**:
  - Applies Hadamard operations on the estimation register to bring them in a superposition, 
  - Executes a specific sequence of **controlled applications of ``U^{2^k}``** from the estimation onto the target register, 
  - Applies an **inverse Quantum Fourier Transform (IQFT)** on the estimation register to recover the phase estimate.

## Details
- The unitary operation ``U`` must act on a number of qubits equal to the length of `targetQubits`.
- A numerical check ensures that ``U`` satisfies unitarity: ``U^†·U ≈ I`` within the specified numerical tolerance.
- If the matrix dimension of ``U`` does not match the number of `targetQubits`, or ``U`` fails the unitarity check, an exception is thrown.
- The estimation and target qubit indexing is interpreted using the circuit’s `settings`.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `estimationQubits::Qubits` — Indices of qubits forming the **estimation register**, used to encode the measured phase in binary form.
- `targetQubits::Qubits` — Indices of qubits forming the **target register**, on which the unitary operation ``U`` acts.
- `U::UnitaryOperation` — The unitary operation ``U`` whose eigenphase(s) will be estimated.
- `accuracyCheckForUnitarity::Float64` (optional) — Numerical tolerance for verifying the unitarity of ``U``, defaulting to `10 * eps(1.0)`.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# Throws
- A `String` exception if the dimensions of ``U`` are incompatible with the number of `targetQubits` or fails the unitarity check.

# See also
- [`createNQubitOperationQPE`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`UnitaryOperation`](@ref)
- [`Qubits`](@ref)
- [`Settings`](@ref)

# Example
Use a **QPE gate** to estimate the eigenphase of a **phase gate** with 4-qubit accuracy:
```julia-repl
julia> accuracy = 4
julia> phasePart = 11/15
julia> theta = 2*pi*phasePart
julia> qc = createQuantumCircuit(accuracy + 1)
julia> xGate!(qc, accuracy + 1) # to create eigenstate |1> with eigenphase theta from initial state |0>
julia> qpeGate!(qc, Vector(1:accuracy), [accuracy + 1], createSingleQubitOperationU1(theta))
julia> measureGate!(qc, Vector(1:accuracy))
julia> qp = compileQuantumCircuit(qc; optimizeNumberOfSteps=true)
julia> iqs = createInitialQubitState(vector, repeat([1. 0.], outer=[accuracy + 1, 1]))
julia> qo = runQuantumProgram(qp, iqs, 1000)
julia> phaseEstimate = 2 * pi * findmax(countmap([shot.measured.outcome for shot in qo.output[3, :]]))[2] / (2^accuracy)
julia> probeMeasureOutcome(qo, 3, "Estimated phase of a "*string(round(theta, digits=3))*" [rad] phase gate: "*string(round(phaseEstimate, digits=3))*" [rad]")
```
"""
function qpeGate!(quantumCircuit::QuantumCircuit, estimationQubits::Qubits, targetQubits::Qubits, U::UnitaryOperation; accuracyCheckForUnitarity = 10*eps(1.0))
	if size(U.U) == (2^length(targetQubits), 2^length(targetQubits))
		if all(abs.(U.U'*U.U-ComplexMatrix(I,size(U.U))).<accuracyCheckForUnitarity)
			addGate!(quantumCircuit, UnitaryGate(
				() -> createNQubitOperationQPE(length(estimationQubits), length(targetQubits), U, quantumCircuit.settings),
				convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, union(estimationQubits, targetQubits)),
				"QPE"
			))
		else
			throw("U:$U is not unitary")
		end
	else
		throw("Dimensions of U are incompatible with the number of target qubits")
	end
end

""" 
	createNQubitOperationQPE(
		numberOfEstimationQubits::Int, 
		numberOfTargetQubits::Int, 
		U::UnitaryOperation, 
		settings::Settings
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Quantum Phase Estimation (QPE) unitary operation** acting on the specified `numberOfEstimationQubits` and `numberOfTargetQubits`, for the given unitary operation ``U`` using qubit indexing `settings`.

The **QPE unitary operation** estimates the eigenphase(s) of the unitary operation ``U``.

Under the hood, it builds the full quantum circuit implementing the **QPE algorithm** and compiles it into a single unitary matrix representing the **QPE unitary operation**.

## Conceptual Overview
- The **QPE algorithm** estimates the eigenphase ``φ`` of a unitary operation ``U`` such that ``U|ψ⟩ = exp(2·π·i·φ)·|ψ⟩``. It operates on two quantum registers:
  - **Estimation register (size `numberOfEstimationQubits`)** — A sequence of qubits that stores the binary representation of the estimated eigenphase ``φ``.
  - **Target register (size `numberOfTargetQubits`)** — The qubits on which the unitary operation ``U`` acts. This register must be initialized in an eigenstate ``|ψ⟩`` of ``U`` (or a superposition of eigenstates) to yield meaningful phase information.
- The **QPE algorithm**:
  - Applies Hadamard operations on the estimation register to bring them in a superposition, 
  - Executes a specific sequence of **controlled applications of ``U^{2^k}``** from the estimation onto the target register, 
  - Applies an **inverse Quantum Fourier Transform (IQFT)** on the estimation register to recover the phase estimate.

## Details
- The unitary operation ``U`` must act on exactly `numberOfTargetQubits` qubits.
- The resulting unitary operation acts on `numberOfEstimationQubits + numberOfTargetQubits` qubits:
  - The **first `numberOfEstimationQubits`** represent the phase estimation register.
  - The **remaining `numberOfTargetQubits`** represent the target register.
- The qubit indexing is determined by the provided `settings`.

# Arguments
- `numberOfEstimationQubits::Int` — The number of qubits forming the **estimation register**, used to encode the measured phase in binary form.
- `numberOfTargetQubits::Int` — The number of qubits forming the **target register**, on which the unitary operation ``U`` acts.
- `U::UnitaryOperation` — The unitary operation ``U`` whose eigenphase(s) will be estimated.
- `settings::Settings` — Qubit indexing settings.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **QPE unitary operation**.

# See also
- [`qpeGate!`](@ref)
- [`createNQubitQPEQuantumCircuit`](@ref)
- [`UnitaryOperation`](@ref)
- [`Settings`](@ref)

# Example
Extract the unitary matrix representation of the **QPE unitary operation** that estimates the eigenphase of a **T-dagger gate** with single-qubit accuracy:
```julia-repl
julia> unitaryOperation = createNQubitOperationQPE(1, 1, createSingleQubitOperationTd(), createSettings(indexType=circuitIndexBigEndian))
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
         1.0-6.12323e-17im       0.0+0.0im       1.01465e-17+6.12323e-17im       0.0+0.0im
         0.0+0.0im          0.853553-0.353553im          0.0+0.0im          0.146447+0.353553im
 1.01465e-17+6.12323e-17im       0.0+0.0im               1.0-1.83697e-16im       0.0+0.0im
         0.0+0.0im          0.146447+0.353553im          0.0+0.0im          0.853553-0.353553im
```
"""
function createNQubitOperationQPE(numberOfEstimationQubits::Int, numberOfTargetQubits::Int, U::UnitaryOperation, settings::Settings)
	quantumCircuit = createNQubitQPEQuantumCircuit(numberOfEstimationQubits, numberOfTargetQubits, U, settings)
	return compileToSingleGate(quantumCircuit)
end

""" 
    createNQubitQPEQuantumCircuit(
		numberOfEstimationQubits::Int, 
		numberOfTargetQubits::Int, 
		U::UnitaryOperation, 
		settings::Settings
	) -> QuantumCircuit{Z, I}

Constructs a **Quantum Phase Estimation (QPE) circuit** that estimates the eigenphase(s) of the unitary operation ``U``.  

## Conceptual Overview
- The **QPE algorithm** estimates the eigenphase ``φ`` of a unitary operation ``U`` such that ``U|ψ⟩ = exp(2·π·i·φ)·|ψ⟩``. It operates on two quantum registers:
  - **Estimation register (size `numberOfEstimationQubits`)** — A sequence of qubits that stores the binary representation of the estimated eigenphase ``φ``.
  - **Target register (size `numberOfTargetQubits`)** — The qubits on which the unitary operation ``U`` acts. This register must be initialized in an eigenstate ``|ψ⟩`` of ``U`` (or a superposition of eigenstates) to yield meaningful phase information.
- The **QPE algorithm**:
  - Applies Hadamard operations on the estimation register to bring them in a superposition, 
  - Executes a specific sequence of **controlled applications of ``U^{2^k}``** from the estimation onto the target register, 
  - Applies an **inverse Quantum Fourier Transform (IQFT)** on the estimation register to recover the phase estimate.

## Details
- The unitary operation ``U`` must act on exactly `numberOfTargetQubits` qubits.
- The total number of qubits in the circuit is `numberOfEstimationQubits + numberOfTargetQubits`:
  - The **first `numberOfEstimationQubits`** represent the phase estimation register.
  - The **remaining `numberOfTargetQubits`** represent the target register.
- The qubit indexing is determined by the provided `settings`.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `numberOfEstimationQubits::Int` — The number of qubits forming the **estimation register**, used to encode the measured phase in binary form.
- `numberOfTargetQubits::Int` — The number of qubits forming the **target register**, on which the unitary operation ``U`` acts.
- `U::UnitaryOperation` — The unitary operation ``U`` whose eigenphase(s) will be estimated.
- `settings::Settings` — Qubit indexing settings.

# Returns
- `QuantumCircuit{Z, I}` — A quantum circuit implementing the **QPE**

# See also
- [`createNQubitOperationQPE`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryOperation`](@ref)
- [`Settings`](@ref)

# Example
Create the quantum circuit of a **QPE** that estimates the eigenphase of a **T-dagger gate** with 2-qubit accuracy:
```julia-repl
julia> quantumCircuit = createNQubitQPEQuantumCircuit(2, 1, createSingleQubitOperationTd(), createSettings(indexType=circuitIndexBigEndian))
julia> quantumCircuit.circuit
5-element Vector{Step}:
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [1], "H"), UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [2], "H")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}}(QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}([2], [3], UnitaryOperation(ComplexF64[1.0 + 0.0im -0.0 + 0.0im; 0.0 + 0.0im 0.7071067811865476 - 0.7071067811865475im])), [2, 3], "Controlled-Unitary")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}}(QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}([1], [3], UnitaryOperation(ComplexF64[1.0 + 0.0im -0.0 + 0.0im; 0.0 + 0.0im 0.7071067811865476 - 0.7071067811865475im])), [1, 3], "Controlled-Unitary")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}}(QubiSim.var"#59#60"{Vector{Int64}, Vector{Int64}, UnitaryOperation}([1], [3], UnitaryOperation(ComplexF64[1.0 + 0.0im -0.0 + 0.0im; 0.0 + 0.0im 0.7071067811865476 - 0.7071067811865475im])), [1, 3], "Controlled-Unitary")])
 UnitaryStep(UnitaryGate[UnitaryGate{QubiSim.var"#57#58"{Vector{Int64}}}(QubiSim.var"#57#58"{Vector{Int64}}([1, 2]), [1, 2], "IQFT")])
```
"""
function createNQubitQPEQuantumCircuit(numberOfEstimationQubits::Int, numberOfTargetQubits::Int, U::UnitaryOperation, settings::Settings)
	numberOfQubits = numberOfEstimationQubits + numberOfTargetQubits
	estimationQubits = Vector(1:numberOfEstimationQubits)
	targetQubits = Vector((numberOfEstimationQubits + 1):numberOfQubits)

	quantumCircuit = createQuantumCircuit(numberOfQubits)
	for m in estimationQubits
        hGate!(quantumCircuit, m)
	end
	estimationQubitsIterator = getEstimationQubitOrder(settings, estimationQubits)

    for m in estimationQubitsIterator
		powerOfU = getPowerOfU(settings, m, numberOfEstimationQubits)
        for _ in 1:powerOfU
            controlledUGate!(quantumCircuit, [m], targetQubits, U)
        end
    end
    iqftGate!(quantumCircuit, getIQFTOrder(settings, estimationQubits))
    return quantumCircuit
end

getEstimationQubitOrder(::Settings{Z, I}, qubits) where {Z, I} = qubits
getEstimationQubitOrder(::Settings{Z, circuitIndexBigEndian}, qubits) where {Z} = reverse(qubits)

getPowerOfU(::Settings{Z, I}, m, numberOfEstimationQubits) where {Z, I} = 2^(m - 1)
getPowerOfU(::Settings{Z, circuitIndexBigEndian}, m, numberOfEstimationQubits) where {Z} = 2^(numberOfEstimationQubits - m)

getIQFTOrder(::Settings{Z, I}, qubits) where {Z, I} = reverse(qubits)
getIQFTOrder(::Settings{Z, circuitIndexBigEndian}, qubits) where {Z} = qubits
