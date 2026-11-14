############################################
# DoubleQubitUnitaryGates

""" 
	cnotGate!(
		quantumCircuit::QuantumCircuit, 
		cqubit::Int, 
		tqubit::Int
	) -> Nothing

Adds a **CNOT gate** to the `quantumCircuit` on the specified control qubit `cqubit` and target qubit `tqubit`.

The **CNOT gate** (also known as the Controlled-NOT gate) performs a unitary operation that flips the target qubit (``|0⟩ ↔ |1⟩``) *if and only if* the control qubit is in the state ``|1⟩``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **CNOT unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationCNOT` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `cqubit::Int` — Index of the control qubit (interpreted using the circuit’s `settings`).
- `tqubit::Int` — Index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationCNOT`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **CNOT gate** to a 3-qubit quantum circuit on control qubit 1 and target qubit 2:
# Example
Add a CNOT gate  to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> cnotGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#41#42"}(QubiSim.var"#41#42"(), [1, 2], "CNOT")
```
"""
function cnotGate!(quantumCircuit::QuantumCircuit, cqubit::Int, tqubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
	    () -> createDoubleQubitOperationCNOT(),
    	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [cqubit,tqubit]),
		"CNOT"
	))
end

""" 
    createDoubleQubitOperationCNOT() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **CNOT unitary operation**.

The **CNOT unitary operation** flips the target qubit (``|0⟩ ↔ |1⟩``) *if and only if* the control qubit is in the state ``|1⟩``.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **CNOT unitary operation**.

# See also
- [`cnotGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **CNOT unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationCNOT()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
```
"""
function createDoubleQubitOperationCNOT()
	UnitaryOperation([1. 0. 0. 0.;0. 1. 0. 0.;0. 0. 0. 1.;0. 0. 1. 0.])
end

#

""" 
	cnotReverseGate!(
		quantumCircuit::QuantumCircuit, 
		cqubit::Int, 
		tqubit::Int
	) -> Nothing

Adds a **reversed CNOT gate** to the `quantumCircuit` on the specified control qubit `cqubit` and target qubit `tqubit`.

The **reversed CNOT gate** (also known as a target-controlled NOT gate) performs a unitary operation that flips the control qubit (``|0⟩ ↔ |1⟩``) *if and only if* the target qubit is in the state ``|1⟩``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **reversed CNOT unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationCNOTReverse` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `cqubit::Int` — Index of the control qubit (interpreted using the circuit’s `settings`).
- `tqubit::Int` — Index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationCNOTReverse`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **reversed CNOT gate** to a 3-qubit quantum circuit on control qubit 1 and target qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> cnotReverseGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#43#44"}(QubiSim.var"#43#44"(), [1, 2], "CNOT-Reverse")
```
"""
function cnotReverseGate!(quantumCircuit::QuantumCircuit, cqubit::Int, tqubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
	    () -> createDoubleQubitOperationCNOTReverse(),
    	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [cqubit,tqubit]),
		"CNOT-Reverse"
	))
end

""" 
    createDoubleQubitOperationCNOTReverse() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **reversed CNOT unitary operation**.

The **reversed CNOT unitary operation** flips the control qubit (``|0⟩ ↔ |1⟩``) *if and only if* the target qubit is in the state ``|1⟩``.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **reversed CNOT unitary operation**.

# See also
- [`cnotReverseGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **reversed CNOT unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationCNOTReverse()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
```
"""
function createDoubleQubitOperationCNOTReverse()
    UnitaryOperation([1. 0. 0. 0.;0. 0. 0. 1.;0. 0. 1. 0.;0. 1. 0. 0.])
end

#

""" 
	swapGate!(
		quantumCircuit::QuantumCircuit, 
		qubit1::Int, 
		qubit2::Int
	) -> Nothing

Adds a **SWAP gate** to the `quantumCircuit` on the specified qubits `qubit1` and `qubit2`.

The **SWAP gate** performs a unitary operation that fully exchanges the basis states between two qubits: ``|a⟩⊗|b⟩ → |b⟩⊗|a⟩``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **SWAP unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationSWAP` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit1::Int` — Index of the first qubit (interpreted using the circuit’s `settings`).
- `qubit2::Int` — Index of the second qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationSWAP`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **SWAP gate** to a 3-qubit quantum circuit on qubits 1 and 2:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> swapGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#45#46"}(QubiSim.var"#45#46"(), [1, 2], "SWAP")
```
"""
function swapGate!(quantumCircuit::QuantumCircuit, qubit1::Int, qubit2::Int)
	addGate!(quantumCircuit, UnitaryGate(
	    () -> createDoubleQubitOperationSWAP(),
    	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit1,qubit2]),
		"SWAP"
	))
end

""" 
    createDoubleQubitOperationSWAP() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **SWAP unitary operation**.

The **SWAP unitary operation** fully exchanges the basis states between two qubits: ``|a⟩⊗|b⟩ → |b⟩⊗|a⟩``.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **SWAP unitary operation**.

# See also
- [`swapGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **SWAP unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationSWAP()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
```
"""
function createDoubleQubitOperationSWAP()
    UnitaryOperation([1. 0. 0. 0.;0. 0. 1. 0.;0. 1. 0. 0.;0. 0. 0. 1.])
end

#

""" 
	phaseGate!(
		quantumCircuit::QuantumCircuit, 
		qubit1::Int, 
		qubit2::Int, 
		phi::Float64
	) -> Nothing

Adds a **controlled phase gate** to the `quantumCircuit` on the specified qubits `qubit1` and `qubit2`, with a phase angle ``ϕ``.

The **controlled phase gate** performs a unitary operation that applies a conditional phase rotation ``|11⟩ ↦ exp(i⋅ϕ)⋅|11⟩`` while leaving the other computational basis states unchanged.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled phase unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationPhase` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit1::Int` — Index of the first qubit (interpreted using the circuit’s `settings`).
- `qubit2::Int` — Index of the second qubit (interpreted using the circuit’s `settings`).
- `phi::Float64` — Phase angle ``ϕ`` in radians (i.e., the amount of phase applied to the ``|11⟩`` state).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationPhase`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **controlled phase gate** with phase angle ``ϕ = 0.5π`` on qubits 1 and 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> phaseGate!(qc, 1, 2, 0.5*pi)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#47#48"{Float64}}(QubiSim.var"#47#48"{Float64}(1.5707963267948966), [1, 2], "Phase")
```
"""
function phaseGate!(quantumCircuit::QuantumCircuit, qubit1::Int, qubit2::Int, phi::Float64)
	addGate!(quantumCircuit, UnitaryGate(
	    () -> createDoubleQubitOperationPhase(phi),
    	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit1,qubit2]),
		"Phase"
	))
end

""" 
    createDoubleQubitOperationPhase(
		phi::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled phase unitary operation**.

The **controlled phase unitary operation** applies a conditional phase rotation ``|11⟩ ↦ exp(i⋅ϕ)⋅|11⟩`` while leaving the other computational basis states unchanged.

# Arguments
- `phi::Float64` — Phase angle ``ϕ`` in radians (i.e., the amount of phase applied to the ``|11⟩`` state).

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled phase unitary operation**.

# See also
- [`phaseGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **controlled phase unitary operation** with phase angle ``ϕ = 0.5π``:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationPhase(0.5*pi)
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+1.0im
```
"""
function createDoubleQubitOperationPhase(phi::Float64)
	UnitaryOperation([1. 0. 0. 0.;0. 1. 0. 0.;0. 0. 1. 0.;0. 0. 0. exp(im*phi)])
end

#

""" 
	controlledHGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled Hadamard gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled Hadamard gate** performs a unitary operation that applies the Hadamard transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled Hadamard operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledH` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledH`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`hGate!`](@ref)

# Example
Add a **controlled Hadamard gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledHGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#61#62"}(QubiSim.var"#61#62"(), [1, 2], "Controlled-H")
```
"""
function controlledHGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledH(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-H"
	))
end

""" 
    createDoubleQubitOperationControlledH() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled Hadamard operation**.

The **controlled Hadamard operation** applies the Hadamard transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled Hadamard operation**.

# See also
- [`controlledHGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationH`](@ref)

# Example
Extract the unitary matrix representation of the **controlled Hadamard unitary operation**:
julia> unitaryOperation = createDoubleQubitOperationControlledH()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im       0.0+0.0im        0.0+0.0im
 0.0+0.0im  1.0+0.0im       0.0+0.0im        0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.707107+0.0im   0.707107+0.0im
 0.0+0.0im  0.0+0.0im  0.707107+0.0im  -0.707107+0.0im
```
"""
function createDoubleQubitOperationControlledH()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationH()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledXGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled Pauli-X gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled Pauli-X gate** performs a unitary operation that applies the Pauli-X transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled Pauli-X operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledX` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledX`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`xGate!`](@ref)

# Example
Add a **controlled Pauli-X gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledXGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#63#64"}(QubiSim.var"#63#64"(), [1, 2], "Controlled-X")
```
"""
function controlledXGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledX(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-X"
	))
end

""" 
    createDoubleQubitOperationControlledX() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled Pauli-X operation**.

The **controlled Pauli-X operation** applies the Pauli-X transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled Pauli-X operation**.

# See also
- [`controlledXGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationX`](@ref)

# Example
Extract the unitary matrix representation of the **controlled Pauli-X unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledX()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  1.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
```
"""
function createDoubleQubitOperationControlledX()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationX()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledYGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled Pauli-Y gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled Pauli-Y gate** performs a unitary operation that applies the Pauli-Y transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled Pauli-Y operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledY` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledY`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`hGate!`](@ref)

# Example
Add a **controlled Pauli-Y gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledYGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#65#66"}(QubiSim.var"#65#66"(), [1, 2], "Controlled-Y")
```
"""
function controlledYGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledY(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-Y"
	))
end

""" 
    createDoubleQubitOperationControlledY() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled Pauli-Y operation**.

The **controlled Pauli-Y operation** applies the Pauli-Y transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled Pauli-Y operation**.

# See also
- [`controlledYGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationY`](@ref)

# Example
Extract the unitary matrix representation of the **controlled Pauli-Y unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledY()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0-1.0im
 0.0+0.0im  0.0+0.0im  0.0+1.0im  0.0+0.0im
```
"""
function createDoubleQubitOperationControlledY()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationY()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledZGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled Pauli-Z gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled Pauli-Z gate** performs a unitary operation that applies the Pauli-Z transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled Pauli-Z operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledZ` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledZ`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`zGate!`](@ref)

# Example
Add a **controlled Pauli-Z gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledZGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#67#68"}(QubiSim.var"#67#68"(), [1, 2], "Controlled-Z")
```
"""
function controlledZGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledZ(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-Z"
	))
end

""" 
    createDoubleQubitOperationControlledZ() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled Pauli-Z operation**.

The **controlled Pauli-Z operation** applies the Pauli-Z transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled Pauli-Z operation**.

# See also
- [`controlledZGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationZ`](@ref)

# Example
Extract the unitary matrix representation of the **controlled Pauli-Z unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledZ()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im   0.0-0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  -1.0+0.0im
```
"""
function createDoubleQubitOperationControlledZ()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationZ()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledTGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled T gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled T gate** performs a unitary operation that applies the T transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled T operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledT` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledT`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`tGate!`](@ref)

# Example
Add a **controlled T gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledTGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#69#70"}(QubiSim.var"#69#70"(), [1, 2], "Controlled-T")
```
"""
function controlledTGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledT(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-T"
	))
end

""" 
    createDoubleQubitOperationControlledT() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled T operation**.

The **controlled T operation** applies the T transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled T operation**.

# See also
- [`controlledTGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationT`](@ref)

# Example
Extract the unitary matrix representation of the **controlled T unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledT()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im       0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im       0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im      -0.0-0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.707107+0.707107im
```
"""
function createDoubleQubitOperationControlledT()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationT()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledTdGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled T-dagger gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled T-dagger gate** performs a unitary operation that applies the T-dagger transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled T-dagger operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledTd` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledTd`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`tdGate!`](@ref)

# Example
Add a **controlled T-dagger gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledTdGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#71#72"}(QubiSim.var"#71#72"(), [1, 2], "Controlled-Td")
```
"""
function controlledTdGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledTd(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-Td"
	))
end

""" 
    createDoubleQubitOperationControlledTd() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled T-dagger operation**.

The **controlled T-dagger operation** applies the T-dagger transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled T-dagger operation**.

# See also
- [`controlledTdGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationTd`](@ref)

# Example
Extract the unitary matrix representation of the **controlled T-dagger unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledTd()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im       0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im       0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im      -0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.707107-0.707107im
```
"""
function createDoubleQubitOperationControlledTd()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationTd()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledSGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled S gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled S gate** performs a unitary operation that applies the S transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled S operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledS` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledS`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`sGate!`](@ref)

# Example
Add a **controlled S gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledSGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#73#74"}(QubiSim.var"#73#74"(), [1, 2], "Controlled-S")
```
"""
function controlledSGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledS(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-S"
	))
end

""" 
    createDoubleQubitOperationControlledS() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled S operation**.

The **controlled S operation** applies the S transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled S operation**.

# See also
- [`controlledSGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationS`](@ref)

# Example
Extract the unitary matrix representation of the **controlled S unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledS()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im   0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  -0.0-0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im   0.0+1.0im
```
"""
function createDoubleQubitOperationControlledS()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationS()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end

#

""" 
	controlledSdGate!(
		quantumCircuit::QuantumCircuit, 
		controlQubit::Int, 
		targetQubit::Int
	) -> Nothing

Adds a **controlled S-dagger gate** to the `quantumCircuit` on the specified qubits `controlQubit` and `targetQubit`.

The **controlled S-dagger gate** performs a unitary operation that applies the S-dagger transform on the `targetQubit` if the `controlQubit` is in the ``|1⟩`` state.
Otherwise it acts as the identity unitary operation on the `targetQubit`.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **controlled S-dagger operation**, and adds 
it to the circuit. During compilation, this factory function calls `createDoubleQubitOperationControlledSd` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `controlQubit::Int` — The index of the control qubit (interpreted using the circuit’s `settings`).
- `targetQubit::Int` — The index of the target qubit (interpreted using the circuit’s `settings`).

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createDoubleQubitOperationControlledSd`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)
- [`sdGate!`](@ref)

# Example
Add a **controlled S-dagger gate** on control qubit 1 and target qubit 2 to a 2-qubit quantum circuit:
```julia-repl
julia> qc = createQuantumCircuit(2)
julia> controlledSdGate!(qc, 1, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#75#76"}(QubiSim.var"#75#76"(), [1, 2], "Controlled-Sd")
```
"""
function controlledSdGate!(quantumCircuit::QuantumCircuit, controlQubit::Int, targetQubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
		() -> createDoubleQubitOperationControlledSd(),
		convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [controlQubit, targetQubit]),
		"Controlled-Sd"
	))
end

""" 
    createDoubleQubitOperationControlledSd() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **controlled S-dagger operation**.

The **controlled S-dagger operation** applies the S-dagger transform on the target qubit if the control qubit is in the ``|1⟩`` state. 
Otherwise it acts as the identity unitary operation on the target qubit.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **controlled S-dagger operation**.

# See also
- [`controlledSdGate!`](@ref)
- [`UnitaryOperation`](@ref)
- [`createSingleQubitOperationSd`](@ref)

# Example
Extract the unitary matrix representation of the **controlled S-dagger unitary operation**:
```julia-repl
julia> unitaryOperation = createDoubleQubitOperationControlledSd()
julia> unitaryOperation.U
4×4 Matrix{ComplexF64}:
 1.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  1.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0-1.0im
```
"""
function createDoubleQubitOperationControlledSd()
	numberOfControlQubits = 1
	numberOfTargetQubits = 1
	UTarget = createSingleQubitOperationSd()
	createNQubitOperationControlledU(numberOfControlQubits, numberOfTargetQubits, UTarget)
end
