set(SRCJULIAIF dace_julia.cxx)

add_library(dacejulia OBJECT ${SRCJULIAIF})
set_target_properties(dacejulia PROPERTIES POSITION_INDEPENDENT_CODE True)
target_link_libraries(dacejulia PUBLIC JlCxx::cxxwrap_julia)
if(WIN32)
    target_compile_definitions(dacejulia PRIVATE "DACE_API=__declspec(dllexport)")
endif(WIN32)
