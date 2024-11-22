#include <iostream>
#include <stdexcept>
#include <cstdlib>
#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>
#include <jlcxx/tuple.hpp>
#include <dace/dace.h>

// map trivial layouts directly, see https://github.com/JuliaInterop/CxxWrap.jl?tab=readme-ov-file#breaking-changes-in-v09
// template<> struct jlcxx::IsMirroredType<DACE::Interval> : std::false_type { };

// macros for evaluation routines
#define EVAL(T, U) \
    mod.method("eval", [](const T& obj, const U& args) { return obj.eval(args); }, \
            "Evaluation of `arg1` with a vector of arguments, `arg2`.");
#define EVAL_COMPILED(T) \
    mod.method("eval", [](const compiledDA& cda, const T& args, T& res) { cda.eval(args, res); }, \
            "Evaluate the compiled polynomial, `arg1`, with a vector of any arithmetic type, `arg2`, and return vector of results, `arg3`.");
#define EVAL_SCALAR(T, U) \
    mod.method("evalScalar", [](const T& obj, const U& arg) { return obj.evalScalar(arg); }, \
            "Evaluation of `arg1` with a single arithmetic type argument, `arg2`.");


// override the exit command (mainly so that the code doesn't exit when init isn't called first)
// not sure how portable or safe this is so has to be enabled in cmake with -DCUSTOM_EXIT=ON
#ifdef CUSTOM_EXIT
void exit(int code) {
    // TODO: should treat code==0 differently?
    throw std::runtime_error("Caught exit call (code=" + std::to_string(code) + ")");
}
#endif


JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    using namespace DACE;

    /***********************************************************************************
    *   Static methods
    ************************************************************************************/
    mod.method("init", [](const unsigned int ord, const unsigned int nvar) {
            DA::init(ord, nvar);
        }, "Initialize the DACE control arrays and set the maximum order, `arg1`, and the maximum number of variables, `arg2`.\n\n"
           "Note: must be called before any other DA routine can be used!");
    mod.method("getMaxOrder", []()->int64_t { return DA::getMaxOrder(); },
            "Return the maximum order currently set for the computation, or zero if undefined.");
    mod.method("getMaxVariables", []()->int64_t { return DA::getMaxVariables(); },
            "Return the maximum number of variables set for the computations, or zero if undefined.");
    mod.method("getMaxMonomials", []()->int64_t { return DA::getMaxMonomials(); },
            "Return the maximum number of monomials available with the order and number of variables specified, or zero if undefined.");
    mod.method("setEps", [](const double eps) { return DA::setEps(eps); },
            "Set the cutoff value eps to `arg1` and return the previous value, or zero if undefined.");
    mod.method("getEps", []() { return DA::getEps(); },
            "Return the cutoff value eps currently set for computations, or zero if undefined.");
    mod.method("getEpsMac", []() { return DA::getEpsMac(); },
            "Return the machine epsilon (pessimistic estimate), or zero if undefined.");
    mod.method("setTO", [](const unsigned int ot)->int64_t { return DA::setTO(ot); },
            "Set the truncation order to `arg1` and return the previous value, or zero if undefined.\n\n"
            "All terms larger than the truncation order are discarded in subsequent operations.");
    mod.method("getTO", []()->int64_t { return DA::getTO(); },
            "Return the truncation order currently set for computations, or zero if undefined.\n\n"
            "All terms larger than the truncation order are discarded in subsequent operations.");
    mod.method("pushTO", [](const unsigned int ot) { DA::pushTO(ot); },
            "Set a new truncation order (`arg1`), saving the previous one on the truncation order stack.\n\n"
            "All terms larger than the truncation order are discarded in subsequent operations.");
    mod.method("popTO", []() { DA::popTO(); },
            "Restore the previous truncation order from the truncation order stack.\n\n"
            "All terms larger than the truncation order are discarded in subsequent operations.");


    /***********************************************************************************
    *   Monomial object
    ************************************************************************************/
    mod.add_type<Monomial>("Monomial")
        .method("order", &Monomial::order)
        .method("toString", &Monomial::toString);

    // override methods in Base
    mod.set_override_module(jl_base_module);
    mod.method("print", [](const Monomial& m) { std::cout << m.toString(); });
    mod.method("println", [](const Monomial& m) { std::cout << m.toString(); });
    mod.unset_override_module();

    // treat Monomial as an element of STL containers
    // jlcxx::stl::apply_stl<Monomial>(mod);


    /***********************************************************************************
    *   Interval object
    ************************************************************************************/
    mod.map_type<Interval>("Interval");

    // treat Interval as an element of STL containers
    // jlcxx::stl::apply_stl<Interval>(mod);


    /***********************************************************************************
    *   DA object
    ************************************************************************************/
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

    // treat DA as an element of STL containers
    // jlcxx::stl::apply_stl<DA>(mod);

    // DA specific methods
    mod.method("getMonomials", [](const DA& da)->std::vector<Monomial> { return da.getMonomials(); },
        "Get vector of all non-zero Monomials for DA `arg1`");
    mod.method("deriv", [](const DA& da, const unsigned int i) { return da.deriv(i); },
        "Derivative of DA `arg1` with respect to given variable, `arg2`");
    mod.method("deriv", [](const DA& da, const std::vector<unsigned int> ind) { return da.deriv(ind); },
        "Derivative of DA `arg1` with respect to given variables, `arg2`");
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
    // mod.method("abs", [](const DA& da) { return da.abs(); });
    mod.method("norm", [](const DA& da, const unsigned int p) { return da.norm(p); });
    mod.method("orderNorm", [](const DA& da, const unsigned int v, const unsigned int p) { return da.orderNorm(v, p); });
    mod.method("estimNorm", [](const DA& da, const unsigned int v, const unsigned int p, const unsigned int o) { return da.estimNorm(v, p, o); });
    mod.method("bound", [](const DA& da) { return da.bound(); });

    // polynomial evaluation routines
    EVAL(DA, std::vector<DA>);
    EVAL(DA, std::vector<double>);

    EVAL_SCALAR(DA, double);
    EVAL_SCALAR(DA, DA);

    // override methods in Base
    mod.set_override_module(jl_base_module);
    // arithmetic operators
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
    // comparison operators
    mod.method("==", [](const DA &da1, const DA &da2) { return da1 == da2; });
    mod.method("==", [](const DA &da, const double c) { return da == c; });
    mod.method("==", [](const double c, const DA &da) { return c == da; });
    mod.method("!=", [](const DA &da1, const DA &da2) { return da1 != da2; });
    mod.method("!=", [](const DA &da, const double c) { return da != c; });
    mod.method("!=", [](const double c, const DA &da) { return c != da; });
    mod.method("<", [](const DA &da1, const DA &da2) { return da1 < da2; });
    mod.method("<", [](const DA &da, const double c) { return da < c; });
    mod.method("<", [](const double c, const DA &da) { return c < da; });
    mod.method(">", [](const DA &da1, const DA &da2) { return da1 > da2; });
    mod.method(">", [](const DA &da, const double c) { return da > c; });
    mod.method(">", [](const double c, const DA &da) { return c > da; });
    mod.method("<=", [](const DA &da1, const DA &da2) { return da1 <= da2; });
    mod.method("<=", [](const DA &da, const double c) { return da <= c; });
    mod.method("<=", [](const double c, const DA &da) { return c <= da; });
    mod.method(">=", [](const DA &da1, const DA &da2) { return da1 >= da2; });
    mod.method(">=", [](const DA &da, const double c) { return da >= c; });
    mod.method(">=", [](const double c, const DA &da) { return c >= da; });
    // math functions
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
    mod.unset_override_module();

    // static factory routines
    mod.method("random", [](const double cm) { return DA::random(cm); });
    mod.method("identity", [](const unsigned int var) { return DA::identity(var); });


    /***********************************************************************************
    *   AlgebraicVector object
    ************************************************************************************/
    mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("AlgebraicVector", jlcxx::julia_type("AbstractVector", "Base"))
            .apply<AlgebraicVector<DA>, AlgebraicVector<double>>([](auto wrapped) {
        typedef typename decltype(wrapped)::type WrappedT;
        typedef typename WrappedT::value_type ScalarT;  // AlgebraicVector inherits from std::vector which sets value_type

        wrapped.template constructor<const size_t>();                // constructor with size
        wrapped.template constructor<const std::vector<ScalarT>&>(); // copy constructor

        wrapped.method("toString", [](const WrappedT& avec) { return avec.toString(); });
        wrapped.method("sqr", [](const WrappedT& avec) { return avec.sqr(); });

        // override methods in Base
        wrapped.module().set_override_module(jl_base_module);
        // implementing the AbstractArray interface
        wrapped.module().method("size", [](const WrappedT& avec)->std::tuple<int64_t> { return std::make_tuple(avec.size()); });
        wrapped.module().method("length", [](const WrappedT& avec)->int64_t { return avec.size(); });
        wrapped.module().method("getindex", [](const WrappedT& avec, const int i) { return avec.at(i-1); });
        wrapped.module().method("setindex!", [](WrappedT& avec, const ScalarT& v, const int i) { avec.at(i-1) = v; });
        wrapped.module().method("firstindex", [](const WrappedT& avec)->int64_t { return 1; });
        wrapped.module().method("lastindex", [](const WrappedT& avec)->int64_t { return avec.size(); });
        // maths functions
        wrapped.module().method("sin", [](const WrappedT& avec) { return sin(avec); });
        wrapped.module().method("cos", [](const WrappedT& avec) { return cos(avec); });
        wrapped.module().method("tan", [](const WrappedT& avec) { return tan(avec); });
        // operators
        // TODO: ...
        wrapped.module().unset_override_module();
    });

    // AlgebraicVector specific methods
    mod.method("trim", [](const AlgebraicVector<DA>& da, const unsigned int min, const unsigned int max) { return da.trim(min, max); });
    mod.method("gradient", [](const DA& da)->AlgebraicVector<DA> { return da.gradient(); });
    mod.method("deriv", [](const AlgebraicVector<DA>& vec, const unsigned int p)->AlgebraicVector<DA> { return vec.deriv(p); });
    mod.method("integ", [](const AlgebraicVector<DA>& obj, const unsigned int p) { return obj.integ(p); });
    mod.method("linear", [](const DA& da)->AlgebraicVector<double> { return da.linear(); });
    mod.method("invert", [](const AlgebraicVector<DA>& vec) { return vec.invert(); });
    mod.method("cons", [](const AlgebraicVector<DA>& vec)->AlgebraicVector<double> { return vec.cons(); });

    // polynomial evaluation routines
    EVAL(DA, AlgebraicVector<DA>);
    EVAL(DA, AlgebraicVector<double>);
    EVAL(AlgebraicVector<DA>, AlgebraicVector<DA>);
    EVAL(AlgebraicVector<DA>, AlgebraicVector<double>);
    EVAL(AlgebraicVector<DA>, std::vector<DA>);
    EVAL(AlgebraicVector<DA>, std::vector<double>);

    EVAL_SCALAR(AlgebraicVector<DA>, double);
    EVAL_SCALAR(AlgebraicVector<DA>, DA);

    // override methods in Base
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
    mod.unset_override_module();


    /***********************************************************************************
    *   compiledDA object
    ************************************************************************************/
    mod.add_type<compiledDA>("compiledDA")
        .constructor<const DA&>()
        .constructor<std::vector<DA>&>()
        .constructor<AlgebraicVector<DA>&>() // TODO how to leverage inheritance here?
        .method("getDim", &compiledDA::getDim)
        .method("getOrd", &compiledDA::getOrd)
        .method("getVars", &compiledDA::getVars)
        .method("getTerms", &compiledDA::getTerms);

    mod.method("compile", [](const DA& da) { return da.compile(); });
    mod.method("compile", [](const std::vector<DA>& vec) { return compiledDA(vec); });
    mod.method("compile", [](const AlgebraicVector<DA>& vec) { return vec.compile(); });

    // polynomial evaluation routines
    EVAL(compiledDA, AlgebraicVector<DA>);
    EVAL(compiledDA, AlgebraicVector<double>);
    EVAL(compiledDA, std::vector<DA>);
    EVAL(compiledDA, std::vector<double>);

    EVAL_COMPILED(AlgebraicVector<DA>);
    EVAL_COMPILED(AlgebraicVector<double>);
    EVAL_COMPILED(std::vector<DA>);
    EVAL_COMPILED(std::vector<double>);

    EVAL_SCALAR(compiledDA, double);
    EVAL_SCALAR(compiledDA, DA);


    /***********************************************************************************
    *   AlgebraicMatrix object
    ************************************************************************************/
    mod.add_type<jlcxx::Parametric<jlcxx::TypeVar<1>>>("AlgebraicMatrix", jlcxx::julia_type("AbstractMatrix", "Base"))
            .apply<AlgebraicMatrix<DA>, AlgebraicMatrix<double>>([](auto wrapped) {
        typedef typename decltype(wrapped)::type WrappedT;
        typedef typename WrappedT::value_type ScalarT; // AlgebraicMatrix encapsulates a std::vector which sets value_type

        wrapped.template constructor<const int>();
        wrapped.template constructor<const int, const int>();
        wrapped.template constructor<const int, const int, const ScalarT&>();

        // pretty print
        wrapped.method("toString", [](const WrappedT& amat) { return amat.toString(); });

        // override methods in Base
        wrapped.module().set_override_module(jl_base_module);
        // implementing the abstract array interface
        // TODO check why the matrix bounds check is not working
        wrapped.module().method("size", [](const WrappedT& amat)->std::tuple<int64_t, int64_t> { return std::make_tuple(amat.nrows(), amat.ncols()); });
        wrapped.module().method("length", [](const WrappedT& amat)->int64_t { return amat.size(); });
        wrapped.module().method("getindex", [](const WrappedT& amat, const int irow, const int icol) { return amat.at(irow-1, icol-1); });
        wrapped.module().method("setindex!", [](WrappedT& amat, const ScalarT& val, const int irow, const int icol) { amat.at(irow-1, icol-1) = val; });
        wrapped.module().method("firstindex", [](const WrappedT& amat)->int64_t { return 1; });
        wrapped.module().method("lastindex", [](const WrappedT& amat)->int64_t { return amat.size(); });
        wrapped.module().unset_override_module();

        // matrix specific methods
        wrapped.method("transpose", [](const WrappedT& amat) { return amat.transpose(); });
        wrapped.method("det", [](const WrappedT& amat) { return amat.det(); });
        wrapped.method("inv", [](const WrappedT& amat) { return amat.inv(); });
        wrapped.method("frobenius", [](const WrappedT& amat) { return amat.frobenius(); });
#ifdef WITH_EIGEN
        wrapped.module().method("eigh", [](const WrappedT& amat) {
            auto [val, vec] = amat.eigh();
            return std::make_tuple(val, vec);
        });
#endif
    });

    mod.method("cons", [](const AlgebraicMatrix<DA>& mat)->AlgebraicMatrix<double> { return mat.cons(); });

    // Jacobian of a DA object (mimics the behavior of AlgebraicVecotr::jacobian() when called on a vector of size 1)
    mod.method("jacobian", [](const DA& da) { return da.jacobian(); });

    // Jacobian and linear part of an AlgebraicVector (requires definition of AlgebraicMatrix)
    mod.method("jacobian", [](const AlgebraicVector<DA>& vec) { return vec.jacobian(); });
    mod.method("linear", [](const AlgebraicVector<DA>& vec) { return vec.linear(); });

    // Hessian matrix of a DA object and of the elements of an AlgebraicVector
    mod.method("hessian", [](const DA& da) { return da.hessian(); });
    mod.method("hessian", [](const AlgebraicVector<DA>& vec) { return vec.hessian(); },
        "Returns a vector of Hessians where the `i`th entry is the Hessian of the `i`th element of `arg1`.");

    // temporary workaround for https://github.com/JuliaInterop/libcxxwrap-julia/issues/173
    /*
    mod.method("hess_vec", [](const AlgebraicVector<DA>& vec) {
        std::vector<AlgebraicMatrix<DA>> hess = vec.hessian();
        jlcxx::Array<AlgebraicMatrix<DA>> out{};
        for (size_t i = 0; i < vec.size(); ++i) out.push_back(hess[i]);
        return out;
    }, "Returns an array of Hessians where the `i`th entry is the Hessian of the `i`th element of `arg1`.");
    */

    /***********************************************************************************
    *   Free functions
    ************************************************************************************/
    mod.method("getMultiIndices", [](const unsigned int no, const unsigned int nv) { return getMultiIndices(no, nv); },
        "Get all multi-indices of order `arg1` in `arg2` variables.");
    mod.method("getRawMoments", [](const DA& mgf, const unsigned int no) {
                auto [mi, rm] = getRawMoments(mgf, no);
                return std::make_tuple(mi, rm);
        }, "Compute the raw moments up to order `arg2` given the moment generating function `arg1`.");
    mod.method("getCentralMoments", [](const DA& mgf, const unsigned int no) {
                auto [mi, cm] = getCentralMoments(mgf, no);
                return std::make_tuple(mi, cm);
        }, "Compute the central moments up to order `arg2` given the moment generating function `arg1`.");
}
