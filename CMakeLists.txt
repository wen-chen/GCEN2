cmake_minimum_required(VERSION 3.18)
project(gcen2_core CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Python + pybind11 ──────────────────────────────────────────────────────
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# ── Compiler flags ─────────────────────────────────────────────────────────
if(MSVC)
    add_compile_options(/W3 /O2)
else()
    # -ffp-contract=off everywhere: no a*b+c fused into an FMA, so results do
    # not depend on the target (AArch64 GCC and -march=native contract by
    # default) and the C++ paths stay bit-identical to their numpy versions
    add_compile_options(-Wall -Wextra -O3 -ffp-contract=off)
    # -march=native gives ~10-15% speedup but produces non-portable binaries;
    # disable for wheel builds, enable for local dev via:
    #   cmake -DGCEN2_NATIVE=ON ..
    option(GCEN2_NATIVE "Enable -march=native" OFF)
    if(GCEN2_NATIVE)
        add_compile_options(-march=native)
    endif()
endif()

# ── Include paths ──────────────────────────────────────────────────────────
include_directories(
    ${CMAKE_SOURCE_DIR}/src          # correlation.hpp, enrichment.hpp, ...
    ${CMAKE_SOURCE_DIR}/src/util     # parallel.hpp
)

# ── Single shared extension: _gcen_core ───────────────────────────────────
pybind11_add_module(_gcen_core
    # pybind11 entry point
    src/bindings/module.cpp

    # Per-submodule bindings
    src/bindings/correlation.cpp
    src/bindings/enrichment.cpp
    src/bindings/module_identify.cpp
    src/bindings/network.cpp
    src/bindings/rwr.cpp
    src/bindings/tom_modules.cpp

    # Implementation files
    src/build_file.cpp
    src/correlation.cpp
    src/edge_shuffle.cpp
    src/enrichment.cpp
    src/module_identify.cpp
    src/network_stats.cpp
    src/edge_files.cpp
    src/edge_table.cpp
    src/rwr.cpp
    src/table_writer.cpp
    src/layer_file.cpp
    src/layer_fuse.cpp
    src/network_stream.cpp
    src/tom.cpp
    src/tom_gemm.cpp
    src/hclust.cpp
)

# TOM products: an AVX build of the kernel, picked at run time when the CPU
# has AVX (same results as the baseline kernel: see src/tom_gemm.hpp). Only
# for an x86 target: with CMAKE_OSX_ARCHITECTURES it must be x86_64 alone (a
# universal2 build compiles every file for arm64 too, and uses the baseline
# kernel: SSE2 / NEON per slice).
set(GCEN2_AVX_KERNEL OFF)
if(NOT MSVC)
    if(CMAKE_OSX_ARCHITECTURES)
        if(CMAKE_OSX_ARCHITECTURES STREQUAL "x86_64")
            set(GCEN2_AVX_KERNEL ON)
        endif()
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64|i[3-6]86)$")
        set(GCEN2_AVX_KERNEL ON)
    endif()
endif()
if(GCEN2_AVX_KERNEL)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-mavx GCEN2_HAVE_MAVX)
    if(NOT GCEN2_HAVE_MAVX)
        set(GCEN2_AVX_KERNEL OFF)
    endif()
endif()
message(STATUS "GCEN2: AVX kernel for the TOM products: ${GCEN2_AVX_KERNEL}")
if(GCEN2_AVX_KERNEL)
    target_sources(_gcen_core PRIVATE src/tom_gemm_avx.cpp)
    set_source_files_properties(src/tom_gemm_avx.cpp PROPERTIES COMPILE_OPTIONS
                                "-mavx")
    set_source_files_properties(src/tom_gemm.cpp PROPERTIES COMPILE_DEFINITIONS GCEN2_AVX_KERNEL)
endif()

target_include_directories(_gcen_core PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/util
)

# ── Install ────────────────────────────────────────────────────────────────
install(TARGETS _gcen_core DESTINATION gcen2/core COMPONENT python)
