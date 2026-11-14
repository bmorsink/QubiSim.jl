############################################
# MeasurementDirectionOperators

""" 
    sigmaX() -> ComplexMatrix

Creates the **measurement direction operator** corresponding to the X-axis on the Bloch sphere.

This operator defines a projective measurement along the X direction and is a special case of the general Bloch-sphere measurement operator expressed in terms of spherical polar angles ``θ`` and ``φ``:
```math
    M(θ, φ) = cos(θ/2) ⋅ σ_z + sin(θ/2) ⋅ (cos(φ) ⋅ σ_x + sin(φ) ⋅ σ_y)
```

where ``σ_x``, ``σ_y``, and ``σ_z`` are the Pauli matrices associated with the X, Y, and Z axes, respectively.

Specifically for `sigmaX()`, the measurement direction corresponds to spherical polar angles ``θ = π/2`` and ``φ = 0``, aligning with the positive X-axis.

# Returns
- `ComplexMatrix` - A 2×2 Hermitian matrix representing the **measurement direction operator** along the X-axis.

# Calling structure
```julia-repl
sigmaX()
    └── sigmaN(theta::Float64, phi::Float64)
        └── sigmaN(n::Vector{Float64})
            └── createSingleQubitMeasureOperator(n::Vector{Float64})
```

# See Also
- [`pauliX`](@ref)
- [`pauliY`](@ref)
- [`pauliZ`](@ref)
- [`sigmaY`](@ref)
- [`sigmaZ`](@ref)
- [`sigmaN`](@ref)
- [`measureGate!`](@ref)
- [`MeasureGateType`](@ref)
- [`MeasureGate`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
- [`Measured`](@ref)

# Example
```julia-repl
julia> sigmaX()
2×2 Matrix{ComplexF64}:
 0.707107+0.0im   0.707107-0.0im
 0.707107+0.0im  -0.707107+0.0im
```
"""
function sigmaX()
	return sigmaN(pi/2, 0.)
	# sx=pauliX() # [0 1;1 0]
	# V=eigvecs(sx)
	# D=eigvals(sx)
	# ind=sortperm(D,rev=true)
	# D=D[ind]
	# V=V[:,ind]
	# V[:,2]=-V[:,2]
	# return V
end

""" 
    sigmaY() -> ComplexMatrix

Creates the **measurement direction operator** corresponding to the Y-axis on the Bloch sphere.

This operator defines a projective measurement along the Y direction and is a special case of the general Bloch-sphere measurement operator expressed in terms of spherical polar angles ``θ`` and ``φ``:
```math
    M(θ, φ) = cos(θ/2) ⋅ σ_z + sin(θ/2) ⋅ (cos(φ) ⋅ σ_x + sin(φ) ⋅ σ_y)
```

where ``σ_x``, ``σ_y``, and ``σ_z`` are the Pauli matrices associated with the X, Y, and Z axes, respectively.

Specifically for `sigmaY()`, the measurement direction corresponds to spherical polar angles ``θ = π/2`` and ``φ = π/2``, aligning with the positive Z-axis.

# Returns
- `ComplexMatrix` - A 2×2 Hermitian matrix representing the **measurement direction operator** along the Y-axis.

# Calling structure
```julia-repl
sigmaY()
    └── sigmaN(theta::Float64, phi::Float64)
        └── sigmaN(n::Vector{Float64})
            └── createSingleQubitMeasureOperator(n::Vector{Float64})
```

# See Also
- [`pauliX`](@ref)
- [`pauliY`](@ref)
- [`pauliZ`](@ref)
- [`sigmaX`](@ref)
- [`sigmaZ`](@ref)
- [`sigmaN`](@ref)
- [`measureGate!`](@ref)
- [`MeasureGateType`](@ref)
- [`MeasureGate`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
- [`Measured`](@ref)

# Example
```julia-repl
julia> sigmaY()
2×2 Matrix{ComplexF64}:
    0.707107+0.0im   0.0-0.707107im
    0.0+0.707107im  -0.707107+0.0im
```
"""
function sigmaY()
	return sigmaN(pi/2, pi/2)
	# sy=pauliY() # [0 -1im;1im 0]
	# V=eigvecs(sy)
	# D=eigvals(sy)
	# ind=sortperm(D,rev=true)
	# D=D[ind]
	# V=V[:,ind]
	# V=1im*V
	# return V
end

""" 
    sigmaZ() -> ComplexMatrix

Creates the **measurement direction operator** corresponding to the Z-axis on the Bloch sphere.

This operator defines a projective measurement along the Z direction and is a special case of the general Bloch-sphere measurement operator expressed in terms of spherical polar angles ``θ`` and ``φ``:
```math
    M(θ, φ) = cos(θ/2) ⋅ σ_z + sin(θ/2) ⋅ (cos(φ) ⋅ σ_x + sin(φ) ⋅ σ_y)
```

where ``σ_x``, ``σ_y``, and ``σ_z`` are the Pauli matrices associated with the X, Y, and Z axes, respectively.

Specifically for `sigmaZ()`, the measurement direction corresponds to spherical polar angles ``θ = 0`` and ``φ = 0``, aligning with the positive Z-axis.

# Returns
- `ComplexMatrix` - A 2×2 Hermitian matrix representing the **measurement direction operator** along the Z-axis.

# Calling structure
```julia-repl
sigmaZ()
    └── sigmaN(theta::Float64, phi::Float64)
        └── sigmaN(n::Vector{Float64})
            └── createSingleQubitMeasureOperator(n::Vector{Float64})
```

# See Also
- [`pauliX`](@ref)
- [`pauliY`](@ref)
- [`pauliZ`](@ref)
- [`sigmaX`](@ref)
- [`sigmaY`](@ref)
- [`sigmaN`](@ref)
- [`measureGate!`](@ref)
- [`MeasureGateType`](@ref)
- [`MeasureGate`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
- [`Measured`](@ref)

# Example
```julia-repl
julia> sigmaZ()
2×2 Matrix{ComplexF64}:
 1.0+0.0im   0.0-0.0im
 0.0+0.0im  -1.0+0.0im
```
"""
function sigmaZ()
	return sigmaN(0., 0.)
	# sz=pauliZ() # [1 0;0 -1]
	# V=eigvecs(sz)
	# D=eigvals(sz)
	# ind=sortperm(D,rev=true)
	# D=D[ind]
	# V=V[:,ind]
	# return V
	# return sz
end

""" 
    sigmaN(
		n::Vector{Float64}
	) -> ComplexMatrix
    
	sigmaN(
		theta::Float64, 
		phi::Float64
	) -> ComplexMatrix

Creates the **measurement direction operator** corresponding to an arbitrary direction on the Bloch sphere, specified either by a 3D unit vector ``n`` or by spherical polar angles ``θ`` and ``φ``.

Each **measurement direction operator** can be parameterized by a 3D unit vector ``n = [n_x, n_y, n_z]`` on the Bloch sphere:
```math
    M(n) = n ⋅ σ = n_x ⋅ σ_x + n_y ⋅ σ_y + n_z ⋅ σ_z
```

with ``σ = [σ_x, σ_y, σ_z]`` the vector of Pauli matrices associated with the X, Y, and Z axes, respectively.

Alternatively, if the unit vector ``n`` is expressed in terms of spherical polar angles ``θ`` and ``φ``, the operator takes the form:
```math
    M(θ, φ) = cos(θ/2) ⋅ σ_z + sin(θ/2) ⋅ (cos(φ) ⋅ σ_x + sin(φ) ⋅ σ_y)
```

# Arguments
- `n::Vector{Float64}` - A 3-dimensional unit vector ``n = [n_x, n_y, n_z]`` specifying a direction on the Bloch sphere.
- `theta::Float64` - Polar angle ``θ`` in radians (angle from the Z-axis).
- `phi::Float64` - Azimuthal angle ``φ`` in radians (angle from the X-axis in the XY-plane).

# Returns
- `ComplexMatrix` - A 2×2 Hermitian matrix representing the **measurement direction operator** for the specified direction.

# Calling structure
```julia-repl
sigmaN(theta::Float64, phi::Float64)
    └── sigmaN(n::Vector{Float64})
        └── createSingleQubitMeasureOperator(n::Vector{Float64})
```

# See Also
- [`pauliX`](@ref)
- [`pauliY`](@ref)
- [`pauliZ`](@ref)
- [`sigmaX`](@ref)
- [`sigmaY`](@ref)
- [`sigmaZ`](@ref)
- [`measureGate!`](@ref)
- [`MeasureGateType`](@ref)
- [`MeasureGate`](@ref)
- [`MeasureOperation`](@ref)
- [`MeasureAndForgetOperation`](@ref)
- [`Sigmas`](@ref)
- [`KrausOperators`](@ref)
- [`generateKrausOperatorsForPVMMeasurement`](@ref)
- [`Measured`](@ref)

# Examples

**Example 1: Using a unit vector direction**

```julia-repl
julia> sigmaN(vec([1 1 1])/sqrt(3))
2×2 Matrix{ComplexF64}:
 0.57735+0.0im       0.57735-0.57735im
 0.57735+0.57735im  -0.57735+0.0im
```

**Example 2: Using spherical coordinates ``θ = π/3``, ``φ = π/4``**

```julia-repl
julia> sigmaN(pi/3., pi/4.)
2×2 Matrix{ComplexF64}:
 0.866025+0.0im        0.353553-0.353553im
 0.353553+0.353553im  -0.866025+0.0im
```
"""
function sigmaN(theta::Float64, phi::Float64)
	return sigmaN(vec([sin(theta/2)*cos(phi) sin(theta/2)*sin(phi) cos(theta/2)]))
end

function sigmaN(n::AbstractVector{<:Number})
	if length(n) == 3
		return createSingleQubitMeasureOperator(Float64.(n))
	else
		throw("length of vector n is not equal to three")
	end
	# V=eigvecs(sn)
	# D=eigvals(sn)
	# ind=sortperm(D,rev=true)
	# D=D[ind]
	# V=V[:,ind]
	# return V
end

# # Calling structure
# ```julia-repl
# createSingleQubitMeasureOperator(theta::Float64, phi::Float64)
#     └── createSingleQubitMeasureOperator(n::Vector{Float64})
#         └── [matrix construction]
# ```

function createSingleQubitMeasureOperator(theta::Float64, phi::Float64)
	return createSingleQubitMeasureOperator(vec([sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]))
end

function createSingleQubitMeasureOperator(n::Vector{Float64})
	# nsigma = a*sigma with sigma the 3 pauli matrices
	# rho = n(1)*[[0,1];[1,0]]+n(2)*[[0,-i];[i,0]]+n(3)*[[1,0];[0,-1]]
	return [n[3] n[1]-im*n[2]; n[1]+im*n[2] -n[3]]
end

""" 
    pauliX() -> ComplexMatrix

# Returns
- `ComplexMatrix` — 2×2 unitary matrix representation of the **Pauli-X unitary operation**.

The **Pauli-X unitary operation** flips the bit of a single qubit from ``|0⟩`` to ``|1⟩`` and vice versa.

# Example
```julia-repl
julia> pauliX()
2×2 Matrix{ComplexF64}:
 0.0+0.0im  1.0+0.0im
 1.0+0.0im  0.0+0.0im
```
"""
function pauliX()
	return ComplexMatrix([0 1;1 0])
end

""" 
    pauliY() -> ComplexMatrix

# Returns
- `ComplexMatrix` — 2×2 unitary matrix representation of the **Pauli-Y unitary operation**.

The **Pauli-Y unitary operation** flips both the bit and the phase of a single qubit.

# Example
```julia-repl
julia> pauliY()
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0-1.0im
 0.0+1.0im  0.0+0.0im
```
"""
function pauliY()
	return ComplexMatrix([0 -1im;1im 0])
end

""" 
    pauliZ() -> ComplexMatrix

# Returns
- `ComplexMatrix` — 2×2 unitary matrix representation of the **Pauli-Z unitary operation**.

The **Pauli-Z unitary operation** flips the phase of a single qubit.

# Example
```julia-repl
julia> pauliZ()
2×2 Matrix{ComplexF64}:
 1.0+0.0im   0.0+0.0im
 0.0+0.0im  -1.0+0.0im
```
"""
function pauliZ()
	return ComplexMatrix([1 0;0 -1])
end
