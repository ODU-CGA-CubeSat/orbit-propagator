# VintiCode

## Build with CMake

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

### Input File

Create `inputStateVect.txt` in `build/` directory, example below (cartesian)

```
4063.75
0
5134.54
0
7.826
0
```

### Running

```bash
./VintiCode
```

Outputs new state vector to `outputStateVect.txt` (ECI cartesian and mean elements)
