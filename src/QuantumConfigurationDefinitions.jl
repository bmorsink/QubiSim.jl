############################################
# QuantumConfigurationDefinitions

""" 
	@enum IndexType circuitIndexBigEndian circuitIndexLittleEndian byteIndex

Enumeration defining the qubit indexing strategy used in a quantum circuit.

# Variants

- `circuitIndexBigEndian`: 
  - Qubits are indexed top-to-bottom in the circuit view (1 to n).
  - This corresponds to left-to-right (1 to n) in the byte representation.

- `circuitIndexLittleEndian`: 
  - Qubits are indexed top-to-bottom (1 to n) in the circuit view.
  - This corresponds to right-to-left (n to 1) in the byte representation.

- `byteIndex`: 
  - Qubits are indexed left-to-right (1 to n) directly in byte form.

# Indexing Illustration

## Byte Index

`byteIndex` — left(1)-to-right(n) in byte-form
```julia-repl
nonzero based index, zero based index
|1>|2>...|n>         |0>|1>...|n-1> 
```

## Big Endian (Circuit View)

`circuitIndexBigEndian` — top(1)-to-bottom(n) in circuit-form corresponds to left(1)-to-right(n) in byte-form\n
```julia-repl
nonzero based index, zero based index
|1>|2>...|n> = |1>   |0>|1>...|n-1> = |0>
               |2>                    |1>
               .                      .
               |n>                    |n-1>
```

## Little Endian (Circuit View)

`circuitIndexLittleEndian` — top(1)-to-bottom(n) in circuit-form corresponds to left(n)-to-right(1) in byte-form
```julia-repl
nonzero based index, zero based index
|n>...|2>|1> = |1>   |n-1>...|1>|0> = |0>
               |2>                    |1>
               .                      .
               |n>                    |n-1>
```

# See also
- [`Settings`](@ref)
"""
@enum IndexType circuitIndexBigEndian circuitIndexLittleEndian byteIndex

""" 
	struct Settings{Z, I}

Compile-time structure defining **global indexing settings** for qubits in a quantum circuit.  
This type encodes its parameters as type-level values for efficient specialization and dispatch.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Constructors
- `Settings(zeroBasedNumbering::Bool, indexType::IndexType)` - Creates a specialized Settings{Z, I} instance, where both parameters `Z` and `I` are encoded at compile time.

# Accessors
- `zeroBasedNumbering(settings::Settings) -> Bool` - Returns whether qubit indices start from 0 or 1:
  - If `true`, qubit indices start from 0.
  - If `false`, qubit indices start from 1.

- `indexType(settings::Settings) -> IndexType` - Returns the index convention used:
  - Determines how qubits are mapped between circuit view and byte representation. 
  - Can be `circuitIndexBigEndian`, `circuitIndexLittleEndian`, or `byteIndex`.

# Indexing Illustration

## Byte Index

`byteIndex` — left(1)-to-right(n) in byte-form
```julia-repl
nonzero based index, zero based index
|1>|2>...|n>         |0>|1>...|n-1> 
```

## Big Endian (Circuit View)

`circuitIndexBigEndian` — top(1)-to-bottom(n) in circuit-form corresponds to left(1)-to-right(n) in byte-form\n
```julia-repl
nonzero based index, zero based index
|1>|2>...|n> = |1>   |0>|1>...|n-1> = |0>
               |2>                    |1>
               .                      .
               |n>                    |n-1>
```

## Little Endian (Circuit View)

`circuitIndexLittleEndian` — top(1)-to-bottom(n) in circuit-form corresponds to left(n)-to-right(1) in byte-form
```julia-repl
nonzero based index, zero based index
|n>...|2>|1> = |1>   |n-1>...|1>|0> = |0>
               |2>                    |1>
               .                      .
               |n>                    |n-1>
```

# See also
- [`IndexType`](@ref)

# Example
```julia-repl
julia> s = Settings(true, circuitIndexBigEndian)
Settings{true, circuitIndexBigEndian}()
julia> zeroBasedNumbering(s)
true
julia> indexType(s)
circuitIndexBigEndian::IndexType = 0
```
"""
struct Settings{Z, I} end
Settings(zeroBasedNumbering::Bool, indexType::IndexType) = Settings{zeroBasedNumbering, indexType}()
zeroBasedNumbering(::Settings{Z, I}) where {Z, I} = Z
indexType(::Settings{Z, I}) where {Z, I} = I

# """
#     convertToOneBasedNumbering(
# 		settings::Settings, 
# 		index
# 	) -> Int

# Converts an index to **1-based numbering** based on the provided settings.

# This function is **compile-time specialized** on the `zeroBasedNumbering` parameter of `settings`, ensuring optimized and branch-free code generation.

# # Behavior
# - If `settings` was constructed with `zeroBasedNumbering = true`, the function **adds 1** to the given `index`.
# - If `zeroBasedNumbering = false`, the `index` is **returned unchanged**.

# # See also
# - [`Settings`](@ref)
# """
convertToOneBasedNumbering(::Settings{true, I}, index) where {I} = index.+1
convertToOneBasedNumbering(::Settings{false, I}, index) where {I} = index

# """
#     convertFromOneToSettingsBasedNumbering(
# 		settings::Settings, 
# 		index
# 	) -> Int

# Converts a **1-based index** to the numbering convention specified in the given `settings`.

# This function is **compile-time specialized** on the `zeroBasedNumbering` parameter of `settings`, ensuring optimized and branch-free code generation.

# # Behavior
# - If `settings` was constructed with `zeroBasedNumbering = true`, the function **subtracts 1** from the given `index`.
# - If `zeroBasedNumbering = false`, the `index` is **returned unchanged**.

# # See also
# - [`Settings`](@ref)
# """
convertFromOneToSettingsBasedNumbering(::Settings{true, I}, index) where {I} = index.-1
convertFromOneToSettingsBasedNumbering(::Settings{false, I}, index) where {I} = index

# """
#     convertToByteIndex(
# 		settings::Settings, 
# 		numberOfQubits::Int, 
# 		index
# 	) -> Int

# Converts a **qubit index** to its corresponding **byte index**, taking into account both the **endianness** and the **indexing base** defined in the given `settings`.

# This function is **compile-time specialized** on the `indexType` and `zeroBasedNumbering` parameters of `settings`, ensuring optimized and branch-free code generation.

# # Behavior
# - If `settings` was constructed with `indexType = circuitIndexBigEndian` or `indexType = byteIndex`, the function returns the same (but one-based converted) `index`.
# - If `indexType = circuitIndexLittleEndian`, the **flipped** one-based converted `index` relative to the `numberOfQubits` is returned.

# # See also
# - [`Settings`](@ref)
# ```
# """
convertToByteIndex(settings::Settings{Z, I}, numberOfQubits::Int, index) where {Z, I} = convertToOneBasedNumbering(settings, index)
convertToByteIndex(settings::Settings{Z, circuitIndexLittleEndian}, numberOfQubits::Int, index) where {Z} = numberOfQubits.-convertToOneBasedNumbering(settings, index).+1

#

""" 
    createSettings(;
		zeroBasedNumbering = false, 
		indexType = circuitIndexBigEndian
	) -> Settings{Z, I}

Creates a compile-time structure defining **global indexing settings** for qubits in a quantum circuit.  
This type encodes its parameters as type-level values for efficient specialization and dispatch.

# Type Parameters
- `Z`: Boolean type parameter indicating whether **zero-based numbering** is used.
- `I`: `IndexType` parameter specifying the **index mapping convention**.

# Arguments
- `zeroBasedNumbering::Bool`: 
  - If `true`, qubit indices start from 0.
  - If `false`, qubit indices start from 1.

- `indexType::IndexType`: 
  - Determines how qubits are mapped between circuit view and byte representation. 
  - Can be `circuitIndexBigEndian`, `circuitIndexLittleEndian`, or `byteIndex`.

# Returns
- `Settings{Z, I}` — Creates a specialized instance, where both parameters `Z` and `I` are encoded at compile time.

# See Also
- [`IndexType`](@ref)
- [`Settings`](@ref)

# Indexing Illustration

## Byte Index

`byteIndex` — left(1)-to-right(n) in byte-form
```julia-repl
nonzero based index, zero based index
|1>|2>...|n>         |0>|1>...|n-1> 
```

## Big Endian (Circuit View)

`circuitIndexBigEndian` — top(1)-to-bottom(n) in circuit-form corresponds to left(1)-to-right(n) in byte-form\n
```julia-repl
nonzero based index, zero based index
|1>|2>...|n> = |1>   |0>|1>...|n-1> = |0>
               |2>                    |1>
               .                      .
               |n>                    |n-1>
```

## Little Endian (Circuit View)

`circuitIndexLittleEndian` — top(1)-to-bottom(n) in circuit-form corresponds to left(n)-to-right(1) in byte-form
```julia-repl
nonzero based index, zero based index
|n>...|2>|1> = |1>   |n-1>...|1>|0> = |0>
               |2>                    |1>
               .                      .
               |n>                    |n-1>
```

# Example
Create the default indexing settings (zeroBasedNumbering = false, indexType = circuitIndexBigEndian)
```julia-repl
julia> settings = createSettings()
Settings{false, circuitIndexBigEndian}()
```
"""
function createSettings(;zeroBasedNumbering = false, indexType = circuitIndexBigEndian)
	Settings(zeroBasedNumbering, indexType)
end
