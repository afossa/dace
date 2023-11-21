#include "dace/AlgebraicVector.h"
#include "dace/DA.h"
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

    // add the Monomial object
    mod.add_type<Monomial>("Monomial")
        .method("order", &Monomial::order)
        .method("toString", &Monomial::toString);

    // add the DA object
    mod.add_type<DA>("DA")
        .constructor<const double>()
        .constructor<const unsigned int, const double>()
        .constructor<const int, const double>()
        .method("getCoefficient", &DA::getCoefficient)
        .method("toString", &DA::toString);

    // TODO: add finaliser(s)??? is it necessary?

    // adding DA methods to Base
    mod.set_override_module(jl_base_module);
    // operators
    mod.method("+", [](const DA& da1, const DA& da2) { return da1 + da2; });
    mod.method("+", [](const DA& da, const double c) { return da + c; });
    mod.method("+", [](const double c, const DA& da) { return c + da; });
    mod.method("-", [](const DA& da1, const DA& da2) { return da1 - da2; });
    mod.method("-", [](const DA& da, const double c) { return da - c; });
    mod.method("-", [](const double c, const DA& da) { return c - da; });
    // maths functions
    mod.method("sin", [](const DA& da) { return da.sin(); });
    mod.method("cos", [](const DA& da) { return da.cos(); });
    mod.method("tan", [](const DA& da) { return da.tan(); });
    // displaying
//    mod.method("show", )  // TODO: how to do show (then don't need print/println below)
    mod.method("print", [](DA& da) { std::cout << da.toString(); });
    mod.method("println", [](DA& da) { std::cout << da.toString(); });
    // end adding methods to base
    mod.unset_override_module();


    // adding AlgebraicVector
    mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("AlgebraicVector")
            .apply<AlgebraicVector<DA>, AlgebraicVector<double>>([](auto wrapped) {
        typedef typename decltype(wrapped)::type WrappedT;
        typedef typename WrappedT::value_type ScalarT;  // AlgebraicVector inherits from std::vector which sets value_type

        wrapped.template constructor<const size_t>();

        wrapped.method("sin", [](const WrappedT& avec) { return sin(avec); });
        wrapped.method("cos", [](const WrappedT& avec) { return cos(avec); });
        wrapped.method("tan", [](const WrappedT& avec) { return tan(avec); });

        // add methods to Base
        wrapped.module().set_override_module(jl_base_module);
        // TODO: add show instead of print(ln)
        wrapped.module().method("print", [](const WrappedT& avec) { std::cout << avec; });
        wrapped.module().method("println", [](const WrappedT& avec) { std::cout << avec; });
        // implementing the indexing interface
        wrapped.module().method("getindex", [](const WrappedT& avec, const int i) { return avec[i-1]; });
        wrapped.module().method("setindex!", [](WrappedT& avec, const ScalarT& v, const int i) { avec[i-1] = v; });
        wrapped.module().method("firstindex", [](const WrappedT& avec) { return 1; });  // TODO: do we want to index from 1?
        wrapped.module().method("lastindex", [](const WrappedT& avec) { return avec.size(); });  // TODO: do we want to index from 1?
        // stop adding methods to base
        wrapped.module().unset_override_module();
    });
}
