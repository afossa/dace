#include "dace/AlgebraicVector.h"
#include "dace/DA.h"
#include <iostream>
#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>
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

    mod.set_override_module(jl_base_module);
    mod.method("print", [](const Monomial& m) { std::cout << m.toString(); });
    mod.method("println", [](const Monomial& m) { std::cout << m.toString(); });
    mod.unset_override_module();

    jlcxx::stl::apply_stl<Monomial>(mod);


    // add the DA object
    mod.add_type<DA>("DA")
        .constructor<>()
        .constructor<const double>()
        .constructor<const unsigned int, const double>()
        .constructor<const int, const double>()
        .method("getCoefficient", &DA::getCoefficient)
        .method("multiplyMonomials", &DA::multiplyMonomials)
        .method("sqr", &DA::sqr)
        .method("sqrt", &DA::sqrt)
        .method("toString", &DA::toString);

    // TODO: add finaliser(s)??? is it necessary?

    mod.method("getMonomials", [](const DA& da)->std::vector<Monomial> { return da.getMonomials(); });
    mod.method("deriv", [](const DA& da, const unsigned int i) { return da.deriv(i); });
    mod.method("deriv", [](const DA& da, const std::vector<unsigned int> ind) { return da.deriv(ind); });
    mod.method("integ", [](const DA& da, const unsigned int i) { return da.integ(i); });
    mod.method("integ", [](const DA& da, const std::vector<unsigned int> ind) { return da.integ(ind); });

    // adding DA methods to Base
    mod.set_override_module(jl_base_module);
    // operators
    mod.method("+", [](const DA& da1, const DA& da2) { return da1 + da2; });
    mod.method("+", [](const DA& da, const double c) { return da + c; });
    mod.method("+", [](const double c, const DA& da) { return c + da; });
    mod.method("-", [](const DA& da1, const DA& da2) { return da1 - da2; });
    mod.method("-", [](const DA& da, const double c) { return da - c; });
    mod.method("-", [](const double c, const DA& da) { return c - da; });
    mod.method("*", [](const DA& da1, const DA& da2) { return da1 * da2; });
    mod.method("*", [](const DA& da, const double c) { return da * c; });
    mod.method("*", [](const double c, const DA& da) { return c * da; });
    mod.method("/", [](const DA& da1, const DA& da2) { return da1 / da2; });
    mod.method("/", [](const DA& da, const double c) { return da / c; });
    mod.method("/", [](const double c, const DA& da) { return c / da; });
    mod.method("^", [](const DA& da, const jlcxx::StrictlyTypedNumber<int> p) { return da.pow(p.value); });
    mod.method("^", [](const DA& da, const jlcxx::StrictlyTypedNumber<double> p) { return da.pow(p.value); });
    // maths functions
    mod.method("sin", [](const DA& da) { return da.sin(); });
    mod.method("cos", [](const DA& da) { return da.cos(); });
    mod.method("tan", [](const DA& da) { return da.tan(); });
    mod.method("asin", [](const DA& da) { return da.asin(); });
    mod.method("acos", [](const DA& da) { return da.acos(); });
    mod.method("atan", [](const DA& da) { return da.atan(); });
    mod.method("sinh", [](const DA& da) { return da.sinh(); });
    mod.method("cosh", [](const DA& da) { return da.cosh(); });
    mod.method("tanh", [](const DA& da) { return da.tanh(); });
    mod.method("asinh", [](const DA& da) { return da.asinh(); });
    mod.method("acosh", [](const DA& da) { return da.acosh(); });
    mod.method("atanh", [](const DA& da) { return da.atanh(); });
    mod.method("exp", [](const DA& da) { return exp(da); });
    mod.method("log", [](const DA& da) { return da.log(); });
    mod.method("log10", [](const DA& da) { return da.log10(); });
    mod.method("log2", [](const DA& da) { return da.log2(); });
    mod.method("sqrt", [](const DA& da) { return sqrt(da); });
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

        //wrapped.method("sin", [](const WrappedT& avec) { return sin(avec); });
        //wrapped.method("cos", [](const WrappedT& avec) { return cos(avec); });
        //wrapped.method("tan", [](const WrappedT& avec) { return tan(avec); });

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
        // maths functions
        wrapped.module().method("sin", [](const WrappedT& avec) { return sin(avec); });
        wrapped.module().method("cos", [](const WrappedT& avec) { return cos(avec); });
        wrapped.module().method("tan", [](const WrappedT& avec) { return tan(avec); });
        // operators
        // TODO: ...
        // stop adding methods to base
        wrapped.module().unset_override_module();
    });

    mod.method("trim", [](const AlgebraicVector<DA>& da, const unsigned int min, const unsigned int max) { return da.trim(min, max); });

    mod.method("gradient", [](const DA& da)->AlgebraicVector<DA> { return da.gradient(); });


    // adding compiledDA
    mod.add_type<compiledDA>("compiledDA")
        .constructor<const DA&>()
        .method("getDim", &compiledDA::getDim)
        .method("getOrd", &compiledDA::getOrd)
        .method("getVars", &compiledDA::getVars)
        .method("getTerms", &compiledDA::getTerms);


    // adding AlgebraicMatrix
    mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("AlgebraicMatrix")
            .apply<AlgebraicMatrix<DA>, AlgebraicMatrix<double>>([](auto wrapped) {
        typedef typename decltype(wrapped)::type WrappedT;
        //typedef typename WrappedT::value_type ScalarT;

        wrapped.template constructor<const int>();
        wrapped.template constructor<const int, const int>();
        //wrapped.template constructor<const int, const int, const ScalarT&>();

        // add methods to base
        wrapped.module().set_override_module(jl_base_module);
        // TODO: add show instead of print/println
        wrapped.module().method("print", [](const WrappedT& amat) { std::cout << amat; });
        wrapped.module().method("println", [](const WrappedT& amat) { std::cout << amat; });
        // implementing the indexing interface
        wrapped.module().method("getindex", [](const WrappedT& amat, const int irow, const int icol)->const auto& { return amat.at(irow, icol); });
        // stop adding methods to base
        wrapped.module().unset_override_module();
    });
}
