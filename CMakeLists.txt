# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required(VERSION 3.5)

# Set the project name
project (orbit-propagator)

# Create a SOURCES variable with a link to all source files to compile
set(SOURCES
    src/Mainvinti.c
    src/Vinti6.c
)

# Add an executable with the above sources
add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME}
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(orbit-propagator PUBLIC m)

include(CTest)

add_executable(test_vinti6
    src/Vinti6.c
    Unity/src/unity.c
    test_vinti6.c
    src/anp.c
    src/bpn2xy.c
    src/c2ixys.c
    src/cal2jd.c
    src/cp.c
    src/cr.c
    src/d2dtf.c
    src/d2tf.c
    src/dat.c
    src/dtf2d.c
    src/era00.c
    src/fad03.c
    src/fae03.c
    src/faf03.c
    src/fal03.c
    src/falp03.c
    src/faom03.c
    src/fapa03.c
    src/fave03.c
    src/fw2m.c
    src/ir.c
    src/jd2cal.c
    src/nut06a.c
    src/pfw06.c
    src/pnm06a.c
    src/pom00.c
    src/rx.c
    src/rxp.c
    src/rxr.c
    src/ry.c
    src/rz.c
    src/s06.c
    src/sp00.c
    src/taitt.c
    src/taiut1.c
    src/tr.c
    src/utctai.c
    src/utcut1.c
    src/xys06a.c
    src/obl06.c
    src/nut00a.c
    src/fama03.c
    src/faju03.c
    src/fasa03.c
    src/faur03.c
    src/fame03.c
    src/coord_conversion.c
)

target_include_directories(test_vinti6
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(test_vinti6 PUBLIC m)

add_test(test_vinti6
    test_vinti6)
