#include <iostream>
#include <jlcxx/jlcxx.hpp>
#include <dace/dace.h>


JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    using namespace DACE;

    // add DA static methods separately
    mod.method("init", [](const unsigned int ord, const unsigned int nvar) {
        DA::init(ord, nvar);
    });
    mod.method("getMaxOrder", []() { return DA::getMaxOrder(); });
    mod.method("getMaxVariables", []() { return DA::getMaxVariables(); });
    mod.method("getMaxMonomials", []() { return DA::getMaxMonomials(); });

    // add the DA object
    mod.add_type<DA>("DA")
        .constructor<const double>()
        .constructor<const unsigned int, const double>()
        .constructor<const int, const double>()
        .method("toString", &DA::toString);

    // TODO: add finaliser(s)??? is it necessary?

    // adding methods to Base
    mod.set_override_module(jl_base_module);
    // maths functions
    mod.method("sin", [](DA& da) { return da.sin(); });
    mod.method("cos", [](DA& da) { return da.cos(); });
    mod.method("tan", [](DA& da) { return da.tan(); });
    // displaying
//    mod.method("show", )
    mod.method("print", [](DA& da) { std::cout << da.toString(); });
    // end adding methods to base
    mod.unset_override_module();
}
