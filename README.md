# VintiCode

## Cloning the repo

Clone the https://github.com/odu-cga-cubesat/VintiCode.git[VintiCode] repo. Don't forget to use `--recurse-submodules` flag, or else you won't pull down some of the code needed to run unit tests.

```bash
git clone --recurse-submodules https://github.com/odu-cga-cubesat/VintiCode.git
cd VintiCode
```

Note: If you accidentally cloned without using `--recurse-submodules`, you can run `git submodule update --init --recursive` to pull down submodules needed to run unit tests.

## Build with CMake

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

### Running main executable

Create an `inputStateVect.txt` in the `build/` directory, example below (ECI)

```
4063.75
0
5134.54
0
7.826
0
```

In the `build/` directory, run main executable

```bash
./VintiCode
```

This outputs a new state vector to `outputStateVect.txt` (ECI cartesian and mean elements)

#### State Vector Example

The rows in `inputStateVect.txt` and `outputStateVect.txt` correspond to [ECI](https://en.wikipedia.org/wiki/Earth-centered_inertial) state vector which gives position and velocity in cartesian coordinates.
```
x  =      2328.9659400000   km
y  =     -5995.2160000000   km
z  =      1719.9789400000   km
xd =         2.9111011300   km/s
yd =         -.9816405300   km/s
zd =        -7.0904992200   km/s
```

### Running unit tests

In the `build/` directory, run

```bash
ctest -V
```
