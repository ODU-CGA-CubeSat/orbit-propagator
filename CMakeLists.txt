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
