# Orbit Propagator

C library for performing orbit propagation, based on Vinti6 algorithm.

## Cloning the repo

Clone the [Orbit Propagator](https://github.com/odu-cga-cubesat/orbit-propagator.git) repo. Don't forget to use `--recurse-submodules` flag, or else you won't pull down some of the code needed to run unit tests.

```bash
git clone --recurse-submodules git@github.com:ODU-CGA-CubeSat/orbit-propogator.git
cd orbit-propagator
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

## Research Papers for Sealion Mission Architecture

### Recommended doctools

It is recommended you install the following doctools for generating PDF or HTML documents.

* [JabRef v5.6](https://github.com/JabRef/jabref/releases/tag/v5.6).
* [asciidoctor-bibtex](https://github.com/asciidoctor/asciidoctor-bibtex#install).

Alternatively, you can install & run the [SeaLion Mission Workspace (Kasm Image)](https://github.com/ODU-CGA-CubeSat/kasm-sealion-workspace)

### Generating PDF or HTML documents

For PDF documents:

```
asciidoctor research/OrbitPropagation.adoc -o dist/OrbitPropagation.pdf -r asciidoctor-pdf -r asciidoctor-diagram -r asciidoctor-bibtex -b pdf
````

For HTML documents:


```
asciidoctor research/OrbitPropagation.adoc -o dist/OrbitPropagation.html -r asciidoctor-diagram -r asciidoctor-bibtex
```

Once you run this step, you can locally view the documentation by opening `dist/abstract.html` in a web browser or by opening `dist/abstract.pdf` in a pdf viewer.
