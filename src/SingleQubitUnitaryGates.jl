############################################
# SingleQubitUnitaryGates

""" 
	xGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **Pauli-X gate** to the `quantumCircuit` on the specified `qubit`.

The **Pauli-X gate** (also known as a quantum NOT gate) performs a unitary operation that flips the state of the qubit from ``|0⟩`` to ``|1⟩`` and vice versa.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **Pauli-X unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationX` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationX`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **Pauli-X gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc=createQuantumCircuit(3)
julia> xGate!(qc,2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#5#6"}(QubiSim.var"#5#6"(), [2], "X")
```
"""
function xGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationX(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"X"
    ))
end

""" 
    createSingleQubitOperationX() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Pauli-X unitary operation**.

The **Pauli-X unitary operation** flips the state of the qubit from ``|0⟩`` to ``|1⟩`` and vice versa.

This unitary operation is equivalent to applying the ``U_3(π, 0, π)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **Pauli-X unitary operation**.

# See also
- [`xGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **Pauli-X unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationX()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
 0.0+0.0im  1.0+0.0im
 1.0+0.0im  0.0+0.0im
```
"""
function createSingleQubitOperationX()
    createSingleQubitOperationU3(pi/1., 0., pi/1.)
end

#

""" 
	yGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **Pauli-Y gate** to the `quantumCircuit` on the specified `qubit`.

The **Pauli-Y gate** performs a unitary operation that introduces both a bit and phase flip on a single qubit.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **Pauli-Y unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationY` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationY`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **Pauli-Y gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> yGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#7#8"}(QubiSim.var"#7#8"(), [2], "Y")
```
"""
function yGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationY(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"Y"
    ))
end

""" 
    createSingleQubitOperationY() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Pauli-Y unitary operation**.

The **Pauli-Y unitary operation** flips both the bit and phase of a single qubit.

This unitary operation is equivalent to applying the ``U_3(π, π/2, π/2)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **Pauli-Y unitary operation**.

# See also
- [`yGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **Pauli-Y unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationY()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0-1.0im
 0.0+1.0im  0.0+0.0im
```
"""
function createSingleQubitOperationY()
    createSingleQubitOperationU3(pi/1., pi/2., pi/2.)
end

#

""" 
	zGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **Pauli-Z gate** to the `quantumCircuit` on the specified `qubit`.

The **Pauli-Z gate** performs a unitary operation that flips the phase of a single qubit.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **Pauli-Z unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationZ` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationZ`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **Pauli-Z gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> zGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#9#10"}(QubiSim.var"#9#10"(), [2], "Z")
```
"""
function zGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationZ(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"Z"
    ))
end

""" 
    createSingleQubitOperationZ() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Pauli-Z unitary operation**.

The **Pauli-Z unitary operation** flips the phase of a single qubit.

This unitary operation is equivalent to applying the ``U_1(π)`` or ``U_3(0, 0, π)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **Pauli-Z unitary operation**.

# See also
- [`zGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **Pauli-Z unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationZ()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
 1.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0+0.0im
```
"""
function createSingleQubitOperationZ()
    createSingleQubitOperationU1(pi/1.)
end

#

""" 
	hGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **Hadamard gate** to the `quantumCircuit` on the specified `qubit`.

The **Hadamard gate** performs a unitary operation that creates a superposition of the ``|0⟩`` and ``|1⟩`` states. It performs both a **bit and phase rotation**, mapping computational basis states as follows:

- ``|0⟩ → (|0⟩ + |1⟩) / √2``
- ``|1⟩ → (|0⟩ - |1⟩) / √2``

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **Hadamard unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationH` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationH`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **Hadamard gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> hGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#11#12"}(QubiSim.var"#11#12"(), [2], "H")
```
"""
function hGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationH(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"H"
    ))
end

""" 
    createSingleQubitOperationH() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **Hadamard unitary operation**.

The **Hadamard unitary operation** creates a superposition of the ``|0⟩`` and ``|1⟩`` states. It performs both a **bit and phase rotation**, mapping computational basis states as follows:

- ``|0⟩ → (|0⟩ + |1⟩) / √2``
- ``|1⟩ → (|0⟩ - |1⟩) / √2``

This unitary operation is equivalent to applying the ``U_2(0, π)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **Hadamard unitary operation**.

# See also
- [`hGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **Hadamard unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationH()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
0.707107+0.0im   0.707107+0.0im
0.707107+0.0im  -0.707107+0.0im
```
"""
function createSingleQubitOperationH()
    createSingleQubitOperationU2(0., pi/1.)
end

#

""" 
	idGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **identity gate** to the `quantumCircuit` on the specified `qubit`.

The **identity gate** performs a unitary operation that leaves the qubit state unchanged.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **identity unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationId` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationId`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **identity gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> idGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#13#14"}(QubiSim.var"#13#14"(), [2], "Id")
```
"""
function idGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationId(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"Id"
    ))
end

""" 
    createSingleQubitOperationId() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **identity unitary operation**.

The **identity unitary operation** leaves the qubit state unchanged.

This unitary operation is equivalent to applying the ``U_1(0)`` or ``U_3(0, 0, 0)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **identity unitary operation**.

# See also
- [`idGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **identity unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationId()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
1.0+0.0im  0.0+0.0im
0.0+0.0im  1.0+0.0im
```
"""
function createSingleQubitOperationId()
    createSingleQubitOperationU1(0.)
end

#

""" 
	rxGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		theta::Float64
	)	-> Nothing

Adds an **RX rotation gate** to the `quantumCircuit` on the specified `qubit`.

The **RX rotation gate** performs a unitary operation that rotates the qubit state around the **X-axis** of the Bloch sphere by an angle ``θ``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **RX rotation unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationRx` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `theta::Float64` — Rotation angle ``θ`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationRx`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **RX rotation gate** by an angle ``θ = 0.5``[rad] to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> rxGate!(qc, 2, 0.5)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#15#16"{Float64}}(QubiSim.var"#15#16"{Float64}(0.5), [2], "RX")
```
"""
function rxGate!(quantumCircuit::QuantumCircuit, qubit::Int, theta::Float64)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationRx(theta),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"RX"
    ))
end

""" 
    createSingleQubitOperationRx(
		theta::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **RX rotation unitary operation**.

The **RX rotation unitary operation** rotates the qubit state around the **X-axis** of the Bloch sphere by an angle ``θ``.

# Arguments
- `theta::Float64` — Rotation angle ``θ`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **RX rotation unitary operation**.

# See also
- [`rxGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **RX rotation unitary operation** by an angle ``θ = 0.5``[rad]:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationRx(0.5)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
0.968912-0.0im  0.0-0.247404im
0.0-0.247404im  0.968912+0.0im
```
"""
function createSingleQubitOperationRx(theta::Float64)
    # createSingleQubitOperationU3(theta, -pi/2., pi/2.)
	UnitaryOperation(exp(-im * (theta/2.)*createSingleQubitMeasureOperator(pi/2, 0.)))
end

#

""" 
	ryGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		theta::Float64
	)	-> Nothing

Adds an **RY rotation gate** to the `quantumCircuit` on the specified `qubit`.

The **RY rotation gate** performs a unitary operation that rotates the qubit state around the **Y-axis** of the Bloch sphere by an angle ``θ``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **RY rotation unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationRy` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `theta::Float64` — Rotation angle ``θ`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationRy`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **RY rotation gate** by an angle ``θ = 0.5``[rad] to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> ryGate!(qc, 2, 0.5)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#17#18"{Float64}}(QubiSim.var"#17#18"{Float64}(0.5), [2], "RY")
```
"""
function ryGate!(quantumCircuit::QuantumCircuit, qubit::Int, theta::Float64)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationRy(theta),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"RY"
    ))
end

""" 
    createSingleQubitOperationRy(
		theta::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **RY rotation unitary operation**.

The **RY rotation unitary operation** rotates the qubit state around the **Y-axis** of the Bloch sphere by an angle ``θ``.

# Arguments
- `theta::Float64` — Rotation angle ``θ`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **RY rotation unitary operation**.

# See also
- [`ryGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **RY rotation unitary operation** by an angle ``θ = 0.5``[rad]:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationRy(0.5)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
0.968912+0.0im  -0.247404+0.0im
0.247404+0.0im   0.968912+0.0im
```
"""
function createSingleQubitOperationRy(theta::Float64)
    # createSingleQubitOperationU3(theta, 0., 0.)
	UnitaryOperation(exp(-im * (theta/2.)*createSingleQubitMeasureOperator(pi/2, pi/2)))
end

#

""" 
	rzGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		theta::Float64
	)	-> Nothing

Adds an **RZ rotation gate** to the `quantumCircuit` on the specified `qubit`.

The **RZ rotation gate** performs a unitary operation that rotates the qubit state around the **Z-axis** of the Bloch sphere by an angle ``θ``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **RZ rotation unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationRz` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `theta::Float64` — Rotation angle ``θ`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationRz`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **RZ rotation gate** by an angle ``θ = 0.5``[rad] to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> rzGate!(qc, 2, 0.5)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#19#20"{Float64}}(QubiSim.var"#19#20"{Float64}(0.5), [2], "RZ")
```
"""
function rzGate!(quantumCircuit::QuantumCircuit, qubit::Int, theta::Float64)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationRz(theta),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"RZ"
    ))
end

""" 
    createSingleQubitOperationRz(
		theta::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **RZ rotation unitary operation**.

The **RZ rotation unitary operation** rotates the qubit state around the **Z-axis** of the Bloch sphere by an angle ``θ``.

# Arguments
- `theta::Float64` — Rotation angle ``θ`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **RZ rotation unitary operation**.

# See also
- [`rzGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **RZ rotation unitary operation** by an angle ``θ = 0.5``[rad]:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationRz(0.5)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
0.968912-0.247404im       0.0+0.0im
0.0+0.0im       0.968912+0.247404im
```
"""
function createSingleQubitOperationRz(theta::Float64)
    # createSingleQubitOperationU3(0., 0., theta) # seems to be wrong
	UnitaryOperation(exp(-im * (theta/2.)*createSingleQubitMeasureOperator(0., 0.)))
end

#

""" 
	rotationGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		theta::Float64, 
		phi::Float64, 
		alpha::Float64
	)	-> Nothing

Adds a **rotation gate** around an arbitrary axis to the `quantumCircuit` on the specified `qubit`.

The **rotation gate** around an arbitrary axis performs a unitary operation that rotates the qubit state around the arbitrary axis by an angle ``α``.

The	arbitrary axis lies in the Bloch sphere and is specified by the polar angles ``θ`` and ``φ``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **rotation unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationRotation` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `theta::Float64` — Rotation angle ``θ`` in radians.
- `phi::Float64` — Rotation angle ``φ`` in radians.
- `alpha::Float64` — Rotation angle ``α`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationRotation`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **rotation gate** around the arbitrary axis with polar angles ``θ=0.2``[rad] and ``φ=0.3``[rad] about angle ``α=0.5``[rad] to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> rotationGate!(qc, 2, 0.2, 0.3, 0.5)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#21#22"{Float64, Float64, Float64}}(QubiSim.var"#21#22"{Float64, Float64, Float64}(0.2, 0.3, 0.5), [2], "Rotation")
```
"""
function rotationGate!(quantumCircuit::QuantumCircuit, qubit::Int, theta::Float64, phi::Float64, alpha::Float64)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationRotation(theta, phi, alpha),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"Rotation"
    ))
end

""" 
    createSingleQubitOperationRotation(
		theta::Float64,
		phi::Float64,
		alpha::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **rotation unitary operation** around an arbitrary axis.

The **rotation unitary operation** rotates the qubit state around the arbitrary axis by an angle ``α``.

The	arbitrary axis lies in the Bloch sphere and is specified by the polar angles ``θ`` and ``φ``.

# Arguments
- `theta::Float64` — Rotation angle ``θ`` in radians.
- `phi::Float64` — Rotation angle ``φ`` in radians.
- `alpha::Float64` — Rotation angle ``α`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **rotation unitary operation** around an arbitrary axis.

# See also
- [`rotationGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **rotation unitary operation** around an arbitrary axis with polar angles ``θ=0.2``[rad] and ``φ=0.3``[rad] about angle ``α=0.5``[rad]:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationRotation(0.2, 0.3, 0.5)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
  0.968912-0.242472im   -0.0145253-0.0469563im
 0.0145253-0.0469563im    0.968912+0.242472im
```
"""
function createSingleQubitOperationRotation(theta::Float64, phi::Float64, alpha::Float64)
	UnitaryOperation(exp(-im * (alpha/2.)*createSingleQubitMeasureOperator(theta, phi)))
end

#

""" 
	tGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **T gate** to the `quantumCircuit` on the specified `qubit`.

The **T gate** (also known as the ``π/8`` gate) performs a unitary operation that rotates the phase of the qubit state by ``π/4`` radians about the Z-axis of the Bloch sphere.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **T unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationT` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationT`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **T gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> tGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#23#24"}(QubiSim.var"#23#24"(), [2], "T")
```
"""
function tGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationT(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"T"
    ))
end

""" 
    createSingleQubitOperationT() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **T unitary operation**.

The **T unitary operation** rotates the phase of the qubit state by ``π/4`` radians about the Z-axis of the Bloch sphere.

This unitary operation is equivalent to applying the ``U_1(π/4)`` or ``U_3(0, 0, π/4)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **T unitary operation**.

# See also
- [`tGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **T unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationT()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
1.0+0.0im  0.0+0.0im
0.0+0.0im  0.707107+0.707107im
```
"""
function createSingleQubitOperationT()
    createSingleQubitOperationU1(pi/4.)
end

#

""" 
	tdGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds a **T-dagger gate** to the `quantumCircuit` on the specified `qubit`.

The **T-dagger gate** (also known as the inverse T gate or the ``-π/8`` gate) performs a unitary operation that rotates the phase of the qubit state by ``-π/4`` radians about the Z-axis of the Bloch sphere.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **T-dagger unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationTd` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationTd`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **T-dagger gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> tdGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#25#26"}(QubiSim.var"#25#26"(), [2], "Td")
```
"""
function tdGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationTd(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"Td"
    ))
end

""" 
    createSingleQubitOperationTd() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **T-dagger unitary operation**.

The **T-dagger unitary operation** rotates the phase of the qubit state by ``-π/4`` radians about the Z-axis of the Bloch sphere.

This unitary operation is equivalent to applying the ``U_1(-π/4)`` or ``U_3(0, 0, -π/4)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **T-dagger unitary operation**.

# See also
- [`tdGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **T-dagger unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationTd()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
1.0+0.0im  0.0+0.0im
0.0+0.0im  0.707107-0.707107im
```
"""
function createSingleQubitOperationTd()
    createSingleQubitOperationU1(-pi/4.)
end

#

""" 
	sGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds an **S gate** to the `quantumCircuit` on the specified `qubit`.

The **S gate** (also known as the Phase gate) performs a unitary operation that rotates the phase of the qubit state by ``π/2`` radians about the Z-axis of the Bloch sphere.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **S unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationS` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationS`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add an **S gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> sGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#27#28"}(QubiSim.var"#27#28"(), [2], "S")
```
"""
function sGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationS(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"S"
    ))
end

""" 
    createSingleQubitOperationS() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **S unitary operation**.

The **S unitary operation** rotates the phase of the qubit state by ``π/2`` radians about the Z-axis of the Bloch sphere.

It is the square root of the **Z unitary operation** (``S² = Z``).

This unitary operation is equivalent to applying the ``U_1(π/2)`` or ``U_3(0, 0, π/2)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **S unitary operation**.

# See also
- [`sGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **S unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationS()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
1.0+0.0im  0.0+0.0im
0.0+0.0im  0.0+1.0im
```
"""
function createSingleQubitOperationS()
    createSingleQubitOperationU1(pi/2.)
end

#

""" 
	sdGate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int
	) -> Nothing

Adds an **S-dagger gate** to the `quantumCircuit` on the specified `qubit`.

The **S-dagger gate** (also known as the inverse S or inverse Phase gate) performs a unitary operation that rotates the phase of the qubit state by ``-π/2`` radians about the Z-axis of the Bloch sphere.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **S-dagger unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationSd` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationSd`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add an **S-dagger gate** to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> sdGate!(qc, 2)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#29#30"}(QubiSim.var"#29#30"(), [2], "Sd")
```
"""
function sdGate!(quantumCircuit::QuantumCircuit, qubit::Int)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationSd(),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"Sd"
    ))
end

""" 
    createSingleQubitOperationSd() -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **S-dagger unitary operation**.

The **S-dagger unitary operation** rotates the phase of the qubit state by ``-π/2`` radians about the Z-axis of the Bloch sphere.

It is the square root of the **Z unitary operation** (``S^†² = Z``).

This unitary operation is equivalent to applying the ``U_1(-π/2)`` or ``U_3(0, 0, -π/2)`` operation from the universal unitary operation set.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **S-dagger unitary operation**.

# See also
- [`sdGate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **S-dagger unitary operation**:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationSd()
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
1.0+0.0im  0.0+0.0im
0.0+0.0im  0.0-1.0im
```
"""
function createSingleQubitOperationSd()
	createSingleQubitOperationU1(-pi/2.)
end

#

""" 
	u1Gate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		lambda::Float64
	) -> Nothing

Adds a **U₁(λ) gate** to the `quantumCircuit` on the specified `qubit` with phase rotation angle ``λ``.

The **U₁(λ) gate** performs a unitary operation that applies a phase rotation by angle ``λ`` around the Z-axis of the Bloch sphere: ``U₁(λ) = U₃(0, 0, λ) = [1, 0; 0, exp(i⋅λ)]``.

This unitary operation is equivalent to `RZ(λ)` up to a global phase.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **U₁(λ) unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationU1` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `lambda::Float64` - The phase rotation angle ``λ`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationU1`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **``U₁(λ)`` gate** with phase rotation ``λ = π`` to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> u1Gate!(qc, 2, 1.0*pi)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#35#36"{Float64}}(QubiSim.var"#35#36"{Float64}(3.141592653589793), [2], "U1")
```
"""
function u1Gate!(quantumCircuit::QuantumCircuit, qubit::Int, lambda::Float64)
	addGate!(quantumCircuit, UnitaryGate(
	    () -> createSingleQubitOperationU1(lambda),
    	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"U1"
	))
end

""" 
    createSingleQubitOperationU1(
		lambda::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **U₁(λ) unitary operation**.

The **U₁(λ) unitary operation** applies a phase rotation by angle ``λ`` around the Z-axis of the Bloch sphere: ``U₁(λ) = U₃(0, 0, λ) = [1, 0; 0, exp(i⋅λ)]``.

This unitary operation is equivalent to `RZ(λ)` up to a global phase.

# Arguments
- `lambda::Float64` - The phase rotation angle ``λ`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **U₁(λ) unitary operation**.

# See also
- [`u1Gate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **U₁(λ) unitary operation** with phase rotation ``λ = π``:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationU1(1.0*pi)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
 1.0+0.0im   0.0-0.0im
 0.0+0.0im  -1.0+0.0im
```
"""
function createSingleQubitOperationU1(lambda::Float64)
    createSingleQubitOperationU(0., 0., lambda)
end

#

""" 
	u2Gate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		phi::Float64, 
		lambda::Float64
	) -> Nothing

Adds a **U₂(φ, λ) gate** to the `quantumCircuit` on the specified `qubit` with phase rotation angles ``φ`` and ``λ``.

The **U₂(φ, λ) gate** performs a unitary operation that applies a double phase rotation defined by: ``U₂(φ, λ) = U₃(π/2, φ, λ) = 1/√2⋅[1, -exp(i⋅λ); exp(i⋅φ), exp(i⋅(φ+λ))]``.

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **U₂(φ, λ) unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationU2` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `phi::Float64` -  The first phase rotation angle ``φ`` in radians.
- `lambda::Float64` - The second phase rotation angle ``λ`` in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationU2`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **U₂(φ, λ) gate** with phase rotations ``φ = 0.5π`` and ``λ = π`` to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> u2Gate!(qc, 2, 0.5*pi, 1.0*pi)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#37#38"{Float64, Float64}}(QubiSim.var"#37#38"{Float64, Float64}(1.5707963267948966, 3.141592653589793), [2], "U2")
```
"""
function u2Gate!(quantumCircuit::QuantumCircuit, qubit::Int, phi::Float64, lambda::Float64)
	addGate!(quantumCircuit, UnitaryGate(
	    () -> createSingleQubitOperationU2(phi, lambda),
    	convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"U2"
	))
end

""" 
    createSingleQubitOperationU2(
		phi::Float64, 
		lambda::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **U₂(φ, λ) unitary operation**.

The **U₂(φ, λ) unitary operation** applies a double phase rotation defined by: ``U₂(φ, λ) = U₃(π/2, φ, λ) = 1/√2⋅[1, -exp(i⋅λ); exp(i⋅φ), exp(i⋅(φ+λ))]``.

# Arguments
- `phi::Float64` - The first phase rotation angle ``φ`` in radians.
- `lambda::Float64` - The second phase rotation angle ``λ`` in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **U₂(φ, λ) unitary operation**.

# See also
- [`u2Gate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **U₂(φ, λ) unitary operation** with phase rotations ``φ = 0.5π`` and ``λ = π``:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationU2(0.5*pi, 1.0*pi)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
    0.707107+0.0im  0.707107+0.0im
    0.0+0.707107im  0.0-0.707107im
```
"""
function createSingleQubitOperationU2(phi::Float64, lambda::Float64)
    createSingleQubitOperationU(pi/2,phi,lambda)
end

#

""" 
	u3Gate!(
		quantumCircuit::QuantumCircuit, 
		qubit::Int, 
		theta::Float64, 
		phi::Float64, 
		lambda::Float64
	) -> Nothing

Adds a **U₃(θ, φ, λ) gate** to the `quantumCircuit` on the specified `qubit` with phase rotation angles ``θ``, ``φ`` and ``λ``.

The **U₃(θ, φ, λ) gate** performs the most general unitary operation that applies a triple phase rotation defined by: ``U₃(θ, φ, λ) = [cos(θ/2), -exp(i⋅λ)⋅sin(θ/2); exp(i⋅φ)⋅sin(θ/2), exp(i⋅(φ+λ))⋅cos(θ/2)]``

Special cases:
- ``U₃(0, 0, 0)`` = Identity
- ``U₃(π/2, φ, λ) = U₂(φ, λ)``
- ``U₃(0, 0, λ) = U₁(λ)``

Under the hood, it creates a `UnitaryGate` holding a factory function that can generate the unitary matrix of the **U₃(θ, φ, λ) unitary operation**, and adds 
it to the circuit. During compilation, this factory function calls `createSingleQubitOperationU3` to produce the unitary matrix.

# Arguments
- `quantumCircuit::QuantumCircuit` — The quantum circuit to which the gate is added.
- `qubit::Int` — The index of the qubit (interpreted using the circuit’s `settings`) where the gate is added.
- `theta::Float64` - Rotation angle ``θ`` around the Y-axis in radians.
- `phi::Float64` - Phase angle ``φ`` for pre-Z rotation in radians.
- `lambda::Float64` - Phase angle ``λ`` for post-Z rotation in radians.

# Returns
- `Nothing` — This function mutates the input `quantumCircuit` by adding the gate in-place.

# See also
- [`createSingleQubitOperationU3`](@ref)
- [`QuantumCircuit`](@ref)
- [`UnitaryGate`](@ref)
- [`Settings`](@ref)

# Example
Add a **U₃(θ, φ, λ) gate** with phase rotations ``θ = 0.25π``, ``φ = 0.5π`` and ``λ = π`` to a 3-qubit quantum circuit on qubit 2:
```julia-repl
julia> qc = createQuantumCircuit(3)
julia> u3Gate!(qc, 2, 0.25*pi, 0.5*pi, 1.0*pi)
julia> getGate(qc, getStep(qc, 1), 1) # first step, first gate
UnitaryGate{QubiSim.var"#39#40"{Float64, Float64, Float64}}(QubiSim.var"#39#40"{Float64, Float64, Float64}(0.7853981633974483, 1.5707963267948966, 3.141592653589793), [2], "U3")
```
"""
function u3Gate!(quantumCircuit::QuantumCircuit, qubit::Int, theta::Float64, phi::Float64, lambda::Float64)
	addGate!(quantumCircuit, UnitaryGate(
        () -> createSingleQubitOperationU3(theta, phi, lambda),
        convertToByteIndex(quantumCircuit.settings, quantumCircuit.numberOfQubits, [qubit]),
		"U3"
    ))
end

""" 
    createSingleQubitOperationU3(
		theta::Float64, 
		phi::Float64, 
		lambda::Float64
	) -> UnitaryOperation

Creates a structure encapsulating the unitary matrix representing the **U₃(θ, φ, λ) unitary operation**.

The **U₃(θ, φ, λ) unitary operation** is the most general unitary operation that applies a triple phase rotation defined by: ``U₃(θ, φ, λ) = [cos(θ/2), -exp(i⋅λ)⋅sin(θ/2); exp(i⋅φ)⋅sin(θ/2), exp(i⋅(φ+λ))⋅cos(θ/2)]``

Special cases:
- ``U₃(0, 0, 0)`` = Identity
- ``U₃(π/2, φ, λ) = U₂(φ, λ)``
- ``U₃(0, 0, λ) = U₁(λ)``

# Arguments
- `theta::Float64` - Rotation angle ``θ`` around the Y-axis in radians.
- `phi::Float64` - Phase angle ``φ`` for pre-Z rotation in radians.
- `lambda::Float64` - Phase angle ``λ`` for post-Z rotation in radians.

# Returns
- `UnitaryOperation` — A structure encapsulating the unitary matrix representation of the **U₃(θ, φ, λ) unitary operation**.

# See also
- [`u3Gate!`](@ref)
- [`UnitaryOperation`](@ref)

# Example
Extract the unitary matrix representation of the **U₃(θ, φ, λ) unitary operation** with phase rotations ``θ = 0.25π``, ``φ = 0.5π`` and ``λ = π``:
```julia-repl
julia> unitaryOperation = createSingleQubitOperationU3(0.25*pi, 0.5*pi, 1.0*pi)
julia> unitaryOperation.U
2×2 Matrix{ComplexF64}:
    0.92388+0.0im   0.382683+0.0im
    0.0+0.382683im  0.0-0.92388im
```
"""
function createSingleQubitOperationU3(theta::Float64, phi::Float64, lambda::Float64)
    createSingleQubitOperationU(theta, phi, lambda)
end

#
