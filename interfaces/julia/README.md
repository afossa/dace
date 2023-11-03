# DACE Julia Interface

This interface is generated using the Julia package
[CxxWrap.jl](https://github.com/JuliaInterop/CxxWrap.jl).

For example, paste the following into *script.jl*

```julia
module DACE
    using CxxWrap

    @wrapmodule(() -> "interfaces/cxx/libdace.so", :define_julia_module)
end

# initialise DACE for 20th-order computations in 1 variable
DACE.init(10, 1)

# initialise x as DA
x = DACE.DA(1)

# compute y = sin(x)
y = sin(x)

# print x and y to screen
println("x")
print(x)
println("y = sin(x)")
print(y)
```

and run it with something like:

```
julia script.jl
```

Note that *CxxWrap* will need to be installed first.
