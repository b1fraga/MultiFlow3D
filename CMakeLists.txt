cmake_minimum_required(VERSION 3.19)
project(M3D_v2 LANGUAGES Fortran)
include(GNUInstallDirs)


# -----------------------
# Options
# -----------------------
option(M3D_v2_USE_JSON "Enable JSON support" OFF)
option(M3D_v2_USE_HDF5 "Enable HDF5 support" OFF)
option(M3D_v2_BUILD_TESTS "Build tests" OFF)

# -----------------------
# Compiler flags
# -----------------------
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

add_compile_options(-O2 -g -cpp)
if(M3D_v2_USE_JSON)
    add_compile_definitions(USE_JSON=1)
else()
    add_compile_definitions(USE_JSON=0)
endif()

if(M3D_v2_USE_HDF5)
    add_compile_definitions(USE_HDF5=1)
else()
    add_compile_definitions(USE_HDF5=0)
endif()

# -----------------------
# Find dependencies
# -----------------------
find_package(OpenMP REQUIRED Fortran)

if(M3D_v2_USE_JSON)
    find_package(jsonfortran-gnu REQUIRED CONFIG)
endif()

if(M3D_v2_USE_HDF5)
    find_package(HDF5 REQUIRED COMPONENTS Fortran)
endif()

if(M3D_v2_BUILD_TESTS)
    find_package(test-drive REQUIRED)
endif()

# -----------------------
# Sources
# -----------------------
set(SOURCES
    src/alloc_dom.f90
    src/alloc_pt.f90
    src/averaging.f90
    src/bounds.f90
    src/bounds_keps.f90
    src/bounds_lsm.f90
    src/checkdt.f90
    src/collision.f90
    src/convection.f90
    src/delta_func.f90
    src/diffusion.f90
    src/eddyvis_1eqn.f90
    src/eddyvis_keps.f90
    src/eddyvis_smag.f90
    src/eddyvis_wale.f90
    src/energy.f90
    src/exchange.f90
    src/exchange_bc.f90
    src/exchange_bcphi.f90
    src/exchange_phi.f90
    src/exchangep.f90
    src/exchangepp.f90
    src/exchangesca.f90
    src/exchangeu.f90
    src/exchangev.f90
    src/exchangew.f90
    src/flosol.f90
    src/imb.f90
    src/initial.f90
    src/init_particle.f90
    src/io.f90
    src/localparameters.f90
    src/log_law.f90
    src/LPT.f90
    src/lsm.f90
    src/mgsolver.f90
    src/module_LSM.f90
    src/module_mpi.f90
    src/module_multidata.f90
    src/module_vars.f90
    src/module_vars_pt.f90
    src/MPI_pt.f90
    src/multiflow3d_sem.f90
    src/newsolv_mg.f90
    src/post.f90
    src/press.f90
    src/roughness_function.f90
    src/rungek.f90
    src/sediment.f90
    src/shapes.f90
    src/sipsol.f90
    src/timesig.f90
    src/wall_function.f90
    src/weno.f90
    app/fdstag.f90
)

set_source_files_properties(
    src/multiflow3d_sem.f90
    src/io.f90
    src/MPI_pt.f90
    app/fdstag.f90
    PROPERTIES COMPILE_OPTIONS "-Wall;-Wextra"
)

if(M3D_v2_USE_JSON)
    list(APPEND SOURCES
        src/json_io.f90
    )
endif()

if(M3D_v2_USE_HDF5)
    list(APPEND SOURCES src/hdf5_io.f90)
endif()

# -----------------------
# Main executable
# -----------------------
add_executable(M3D_v2.exe ${SOURCES})

target_include_directories(M3D_v2.exe PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/app
)

target_link_libraries(M3D_v2.exe PRIVATE
    OpenMP::OpenMP_Fortran
)

if(M3D_v2_USE_JSON)
    target_link_libraries(M3D_v2.exe PRIVATE jsonfortran-gnu::jsonfortran-static)
endif()

if(M3D_v2_USE_HDF5)
    target_link_libraries(M3D_v2.exe PRIVATE hdf5::hdf5_fortran)
endif()

install(TARGETS M3D_v2.exe
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)


if(M3D_v2_BUILD_TESTS)
    # -----------------------
    # Tests
    # -----------------------
    set(TEST_SOURCES tests/main.f90)

    if(M3D_v2_USE_JSON)    
    	list(APPEND TEST_SOURCES
            src/module_multidata.f90
            src/module_vars.f90
            src/module_mpi.f90
            src/io.f90
            src/json_io.f90
	    tests/test_json_io.f90
	    tests/test_io.f90
        )
    endif()

    if(M3D_v2_USE_HDF5)
        list(APPEND TEST_SOURCES src/hdf5_io.f90 tests/test_hdf5_io.f90)
    endif()

    add_executable(M3D_v2-tests.exe ${TEST_SOURCES})

    add_dependencies(M3D_v2-tests.exe M3D_v2.exe)

    target_link_libraries(M3D_v2-tests.exe PRIVATE
        OpenMP::OpenMP_Fortran
        test-drive::test-drive
    )

    if(M3D_v2_USE_JSON)
        target_link_libraries(M3D_v2-tests.exe PRIVATE
            jsonfortran-gnu::jsonfortran-static
        )
    endif()

    if(M3D_v2_USE_HDF5)
        target_link_libraries(M3D_v2-tests.exe PRIVATE
            hdf5::hdf5_fortran
        )
    endif()

    # -----------------------
    # Test data copy
    # -----------------------
    if(M3D_v2_USE_JSON)
        add_custom_command(TARGET M3D_v2-tests.exe POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/tests/test_file.json
                $<TARGET_FILE_DIR:M3D_v2-tests.exe>/
            COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/tests/test_io.json
                $<TARGET_FILE_DIR:M3D_v2-tests.exe>/
        )
    endif()

    # -----------------------
    # Run tests
    # -----------------------
    add_custom_target(check-M3D_v2
       COMMAND $<TARGET_FILE:M3D_v2-tests.exe>
       DEPENDS M3D_v2-tests.exe
       WORKING_DIRECTORY $<TARGET_FILE_DIR:M3D_v2-tests.exe>
    )
endif()
