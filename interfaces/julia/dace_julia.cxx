#include <string>
#include <jlcxx/jlcxx.hpp>
#include <dace/dace.h>

using namespace DACE;


std::string greet() {
    return "Hello, World!";
}


JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {
    mod.method("greet", &greet);
}
