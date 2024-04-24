#include "dace/AlgebraicVector.h"
#include "dace/DA.h"
#include <iostream>
#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>
#include <jlcxx/tuple.hpp>
#include <dace/dace.h>

// map trivial layouts directly, see "Breaking changes in v0.9" in CxxWrap README
// template<> struct jlcxx::IsMirroredType<DACE::Interval> : std::false_type { };

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    using namespace DACE;

    // add DA static methods separately
    mod.method("init", [](const unsigned int ord, const unsigned int nvar) {
            DA::init(ord, nvar);
        });
    mod.method("getMaxOrder", []()->int64_t { return DA::getMaxOrder(); });
    mod.method("getMaxVariables", []()->int64_t { return DA::getMaxVariables(); });
    mod.method("getMaxMonomials", []()->int64_t { return DA::getMaxMonomials(); });
    mod.method("setEps", [](const double eps) { return DA::setEps(eps); });
    mod.method("getEps", []() { return DA::getEps(); });
    mod.method("getEpsMac", []() { return DA::getEpsMac(); });
    mod.method("setTO", [](const unsigned int ot)->int64_t { return DA::setTO(ot); });
    mod.method("getTO", []()->int64_t { return DA::getTO(); });
    mod.method("pushTO", [](const unsigned int ot) { DA::pushTO(ot); });
    mod.method("popTO", []() { DA::popTO(); });

    // add the Monomial object
    mod.add_type<Monomial>("Monomial")
        .method("order", &Monomial::order)
        .method("toString", &Monomial::toString);

    mod.set_override_module(jl_base_module);
    mod.method("print", [](const Monomial& m) { std::cout << m.toString(); });
    mod.method("println", [](const Monomial& m) { std::cout << m.toString(); });
    mod.unset_override_module();

    jlcxx::stl::apply_stl<Monomial>(mod);

    // map the Interval object to a predefined Julia structure
    mod.map_type<Interval>("Interval");
    jlcxx::stl::apply_stl<Interval>(mod);

    // add the DA object
    mod.add_type<DA>("DA", jlcxx::julia_type("Real", "Base"))
        .constructor<>()
        .constructor<const double>()
        .constructor<const int, const double>()
        .method("getCoefficient", &DA::getCoefficient)
        .method("getCoefficient", [](const DA& da, jlcxx::ArrayRef<unsigned int> jj) {
                std::vector<unsigned int> jjvec(jj.begin(), jj.end());
                return da.getCoefficient(jjvec);
            })
        .method("multiplyMonomials", &DA::multiplyMonomials)
        .method("sqr", &DA::sqr)
        .method("cons", &DA::cons)
        .method("toString", &DA::toString);

    // TODO: add finaliser(s)??? is it necessary?

    jlcxx::stl::apply_stl<DA>(mod);

    mod.method("getMonomials", [](const DA& da)->std::vector<Monomial> { return da.getMonomials(); });
    mod.method("deriv", [](const DA& da, const unsigned int i) { return da.deriv(i); });
    mod.method("deriv", [](const DA& da, const std::vector<unsigned int> ind) { return da.deriv(ind); });
    mod.method("integ", [](const DA& da, const unsigned int i) { return da.integ(i); });
    mod.method("integ", [](const DA& da, const std::vector<unsigned int> ind) { return da.integ(ind); });
    mod.method("PsiFunction", [](const DA& da, const unsigned int n) { return da.PsiFunction(n); });
    mod.method("isrt", [](const DA& da) { return da.isrt(); });
    mod.method("icrt", [](const DA& da) { return da.icrt(); });
    mod.method("root", [](const DA& da, const int p) { return da.root(p); });
    mod.method("root", [](const DA& da) { return da.root(2); });
    mod.method("mod", [](const DA& da, const double p) { return da.mod(p); });
    mod.method("trim", [](const DA& da, const unsigned int min, const unsigned int max) { return da.trim(min, max); });
    mod.method("divide", [](const DA& da, const unsigned int var, const unsigned int p) { return da.divide(var, p); });
    mod.method("erf", [](const DA& da) { return da.erf(); });
    mod.method("erfc", [](const DA& da) { return da.erfc(); });
    mod.method("besselj", [](const int n, const DA& da) { return da.BesselJFunction(n); });
    mod.method("bessely", [](const int n, const DA& da) { return da.BesselYFunction(n); });
    mod.method("besseli", [](const int n, const DA& da) { return da.BesselIFunction(n); });
    mod.method("besselk", [](const int n, const DA& da) { return da.BesselKFunction(n); });
    mod.method("gamma", [](const DA& da) { return da.GammaFunction(); });
    mod.method("loggamma", [](const DA& da) { return da.LogGammaFunction(); });
    mod.method("powi", [](const DA& da, const int p) { return da.pow(p); });
    mod.method("powd", [](const DA& da, const double p) { return da.pow(p); });

    // norm and estimation routines
    mod.method("abs", [](const DA& da) { return da.abs(); });
    mod.method("norm", [](const DA& da, const unsigned int p) { return da.norm(p); });
    mod.method("orderNorm", [](const DA& da, const unsigned int v, const unsigned int p) { return da.orderNorm(v, p); });
    mod.method("estimNorm", [](const DA& da, const unsigned int v, const unsigned int p, const unsigned int o) { return da.estimNorm(v, p, o); });
    mod.method("bound", [](const DA& da) { return da.bound(); });

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
    mod.method("-", [](const DA& da) { return -da; });
    // maths functions
    mod.method("sin", [](const DA& da) { return da.sin(); });
    mod.method("cos", [](const DA& da) { return da.cos(); });
    mod.method("tan", [](const DA& da) { return da.tan(); });
    mod.method("asin", [](const DA& da) { return da.asin(); });
    mod.method("acos", [](const DA& da) { return da.acos(); });
    mod.method("atan", [](const DA& da) { return da.atan(); });
    mod.method("atan", [](const DA& da1, const DA& da2) { return da1.atan2(da2); });
    mod.method("sinh", [](const DA& da) { return da.sinh(); });
    mod.method("cosh", [](const DA& da) { return da.cosh(); });
    mod.method("tanh", [](const DA& da) { return da.tanh(); });
    mod.method("asinh", [](const DA& da) { return da.asinh(); });
    mod.method("acosh", [](const DA& da) { return da.acosh(); });
    mod.method("atanh", [](const DA& da) { return da.atanh(); });
    mod.method("exp", [](const DA& da) { return da.exp(); });
    mod.method("log", [](const DA& da) { return da.log(); });
    mod.method("log", [](const double b, const DA& da) { return da.logb(b); });
    mod.method("log10", [](const DA& da) { return da.log10(); });
    mod.method("log2", [](const DA& da) { return da.log2(); });
    mod.method("sqrt", [](const DA& da) { return da.sqrt(); });
    mod.method("cbrt", [](const DA& da) { return da.cbrt(); });
    mod.method("hypot", [](const DA& da1, const DA& da2) { return da1.hypot(da2); });
    mod.method("inv", [](const DA& da) { return da.minv(); });
    mod.method("round", [](const DA& da) { return da.round(); });
    mod.method("trunc", [](const DA& da) { return da.trunc(); });
    // end adding methods to base
    mod.unset_override_module();


    // adding AlgebraicVector
    mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("AlgebraicVector", jlcxx::julia_type("AbstractVector", "Base"))
            .apply<AlgebraicVector<DA>, AlgebraicVector<double>>([](auto wrapped) {
        typedef typename decltype(wrapped)::type WrappedT;
        typedef typename WrappedT::value_type ScalarT;  // AlgebraicVector inherits from std::vector which sets value_type

        wrapped.template constructor<const size_t>();                // constructor with size
        wrapped.template constructor<const std::vector<ScalarT>&>(); // copy constructor

        wrapped.method("toString", [](const WrappedT& avec) { return avec.toString(); });
        wrapped.method("sqr", [](const WrappedT& avec) { return avec.sqr(); });

        // add methods to Base
        wrapped.module().set_override_module(jl_base_module);
        // implementing the AbstractArray interface
        wrapped.module().method("size", [](const WrappedT& avec) { return std::make_tuple(avec.size()); });
        wrapped.module().method("getindex", [](const WrappedT& avec, const int i) { return avec.at(i-1); });
        wrapped.module().method("setindex!", [](WrappedT& avec, const ScalarT& v, const int i) { avec.at(i-1) = v; });
        wrapped.module().method("firstindex", [](const WrappedT& avec) { return 1; });
        wrapped.module().method("lastindex", [](const WrappedT& avec) { return avec.size(); });
        // maths functions
        wrapped.module().method("sin", [](const WrappedT& avec) { return sin(avec); });
        wrapped.module().method("cos", [](const WrappedT& avec) { return cos(avec); });
        wrapped.module().method("tan", [](const WrappedT& avec) { return tan(avec); });
        // operators
        // TODO: ...
        // stop adding methods to base
        wrapped.module().unset_override_module();

        // TODO how can we use templates to wrap single functions?
        wrapped.method("eval", [](const DA& da, const WrappedT& args) { return da.eval(args); });
        wrapped.method("eval", [](const AlgebraicVector<DA>& vec, const WrappedT& args) { return vec.eval(args); });
    });

    mod.method("trim", [](const AlgebraicVector<DA>& da, const unsigned int min, const unsigned int max) { return da.trim(min, max); });

    mod.method("gradient", [](const DA& da)->AlgebraicVector<DA> { return da.gradient(); });
    mod.method("deriv", [](const AlgebraicVector<DA>& vec, const unsigned int p)->AlgebraicVector<DA> { return vec.deriv(p); });
    mod.method("integ", [](const AlgebraicVector<DA>& obj, const unsigned int p) { return obj.integ(p); });
    mod.method("linear", [](const DA& da)->AlgebraicVector<double> { return da.linear(); });
    mod.method("invert", [](const AlgebraicVector<DA>& vec) { return vec.invert(); });
    mod.method("cons", [](const AlgebraicVector<DA>& vec)->AlgebraicVector<double> { return vec.cons(); });
    // mod.method("eval", [](const AlgebraicVector<DA>& obj, AlgebraicVector<DA>& args) { return obj.eval(args); });

    // add AlgebraicVector methods to Base
    mod.set_override_module(jl_base_module);
    mod.method("+", [](const AlgebraicVector<DA>& vec1, const AlgebraicVector<DA>& vec2) { return vec1 + vec2; });
    mod.method("+", [](const AlgebraicVector<DA>& vec, const double scalar) { return vec + scalar; });
    mod.method("+", [](const double scalar, const AlgebraicVector<DA>& vec) { return vec + scalar; });
    mod.method("+", [](const AlgebraicVector<DA>& vec, const DA& scalar) { return vec + scalar; });
    mod.method("+", [](const DA& scalar, const AlgebraicVector<DA>& vec) { return vec + scalar; });
    mod.method("-", [](const AlgebraicVector<DA>& vec1, const AlgebraicVector<DA>& vec2) { return vec1 - vec2; });
    mod.method("-", [](const AlgebraicVector<DA>& vec, const double scalar) { return vec - scalar; });
    mod.method("-", [](const double scalar, const AlgebraicVector<DA>& vec) { return scalar - vec; });
    mod.method("-", [](const AlgebraicVector<DA>& vec, const DA& scalar) { return vec - scalar; });
    mod.method("-", [](const DA& scalar, const AlgebraicVector<DA>& vec) { return scalar - vec; });
    mod.method("*", [](const AlgebraicVector<DA>& vec1, const AlgebraicVector<DA>& vec2) { return vec1 * vec2; });
    mod.method("*", [](const AlgebraicVector<DA>& vec, const double scalar) { return vec * scalar; });
    mod.method("*", [](const double scalar, const AlgebraicVector<DA>& vec) { return scalar * vec; });
    mod.method("*", [](const AlgebraicVector<DA>& vec, const DA& scalar) { return vec * scalar; });
    mod.method("*", [](const DA& scalar, const AlgebraicVector<DA>& vec) { return scalar * vec; });
    mod.method("/", [](const AlgebraicVector<DA>& vec1, const AlgebraicVector<DA>& vec2) { return vec1 / vec2; });
    mod.method("/", [](const AlgebraicVector<DA>& vec, const double scalar) { return vec / scalar; });
    mod.method("/", [](const double scalar, const AlgebraicVector<DA>& vec) { return scalar / vec; });
    mod.method("/", [](const AlgebraicVector<DA>& vec, const DA& scalar) { return vec / scalar; });
    mod.method("/", [](const DA& scalar, const AlgebraicVector<DA>& vec) { return scalar / vec; });
    // end adding methods to base
    mod.unset_override_module();

    // adding compiledDA
    mod.add_type<compiledDA>("compiledDA")
        .constructor<const DA&>()
        .constructor<std::vector<DA>&>()
        .constructor<AlgebraicVector<DA>&>() // TODO how to leverage inheritance here?
        .method("getDim", &compiledDA::getDim)
        .method("getOrd", &compiledDA::getOrd)
        .method("getVars", &compiledDA::getVars)
        .method("getTerms", &compiledDA::getTerms);

    // DA polynomial evaluation routines
    mod.method("compile", [](const DA& da) { return da.compile(); });

    // TODO how to avoid duplicate code?
    mod.method("eval", [](const compiledDA& cda, const AlgebraicVector<DA>& args) { return cda.eval(args); });
    mod.method("eval", [](const compiledDA& cda, const AlgebraicVector<double>& args) { return cda.eval(args); });
    mod.method("eval", [](const compiledDA& cda, const std::vector<DA>& args) { return cda.eval(args); });
    mod.method("eval", [](const compiledDA& cda, const std::vector<double>& args) { return cda.eval(args); });
    mod.method("eval", [](const compiledDA& cda, std::vector<double>& args, std::vector<double>& res) { cda.eval(args, res); });

    mod.method("evalScalar", [](const DA& da, const double arg) { return da.evalScalar(arg); });

    // adding AlgebraicMatrix
    mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("AlgebraicMatrix")
            .apply<AlgebraicMatrix<DA>, AlgebraicMatrix<double>>([](auto wrapped) {
        typedef typename decltype(wrapped)::type WrappedT;
        // TODO: how to get the scalar type (i.e. DA or double)

        wrapped.template constructor<const int>();
        wrapped.template constructor<const int, const int>();
        //wrapped.template constructor<const int, const int, const ScalarT&>();

        // add methods to base
        wrapped.module().set_override_module(jl_base_module);
        // implementing the abstract array interface
        wrapped.module().method("getindex", [](const WrappedT& amat, const int irow, const int icol)->const auto& { return amat.at(irow, icol); });
        wrapped.module().method("size", [](const WrappedT& amat) { return std::make_tuple(amat.nrows(), amat.ncols()); });
        wrapped.module().method("length", [](const WrappedT& amat) { return amat.size(); });
        // stop adding methods to base
        wrapped.module().unset_override_module();
    });
}
