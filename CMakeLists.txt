# Set the minimum version of CMake that can be used
# To find the cmake version run
# $ cmake --version
cmake_minimum_required(VERSION 3.5)

# Set the project name
project (VintiCode)

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

target_link_libraries(VintiCode PUBLIC m)

include(CTest)

add_executable(test_vinti6
    src/Vinti6.c
    Unity/src/unity.c
    test_vinti6.c
)

target_include_directories(test_vinti6
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(test_vinti6 PUBLIC m)

add_test(test_vinti6
    test_vinti6)
