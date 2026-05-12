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
    alloc_dom.f90
    alloc_pt.f90
    averaging.f90
    bounds.f90
    bounds_keps.f90
    bounds_lsm.f90
    checkdt.f90
    collision.f90
    convection.f90
    delta_func.f90
    diffusion.f90
    eddyvis_1eqn.f90
    eddyvis_keps.f90
    eddyvis_smag.f90
    eddyvis_wale.f90
    energy.f90
    exchange.f90
    exchange_bc.f90
    exchange_bcphi.f90
    exchange_phi.f90
    exchangep.f90
    exchangepp.f90
    exchangesca.f90
    exchangeu.f90
    exchangev.f90
    exchangew.f90
    fdstag.f90
    flosol.f90
    imb.f90
    initial.f90
    init_particle.f90
    localparameters.f90
    log_law.f90
    LPT.f90
    lsm.f90
    mgsolver.f90
    module_LSM.f90
    module_mpi.f90
    module_multidata.f90
    module_vars.f90
    module_vars_pt.f90
    MPI_pt.f90
    multiflow3d_sem.f90
    newsolv_mg.f90
    post.f90
    press.f90
    roughness_function.f90
    rungek.f90
    sediment.f90
    shapes.f90
    sipsol.f90
    timesig.f90
    wall_function.f90
    weno.f90
)

set_source_files_properties(multiflow3d_sem.f90
    PROPERTIES COMPILE_OPTIONS "-Wall;-Wextra"
)

if(M3D_v2_USE_JSON)
    list(APPEND SOURCES
        io.f90
        json_io.f90
    )
endif()

if(M3D_v2_USE_HDF5)
    list(APPEND SOURCES hdf5_io.f90)
endif()

# -----------------------
# Main executable
# -----------------------
add_executable(M3D_v2.exe ${SOURCES})

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
            io.f90
            json_io.f90
	    tests/test_json_io.f90
	    tests/test_io.f90
        )
    endif()

    if(M3D_v2_USE_HDF5)
        list(APPEND TEST_SOURCES hdf5_io.f90 tests/test_hdf5_io.f90)
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
