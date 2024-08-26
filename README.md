# Realtime Compilation for Continuous Angle QEC Architectures

Compiler for the paper "Realtime Compilation for Continuous Angle QEC Architectures".

## Requirements
Boost is required for the compiler. You can install it with the following command for Debian-based systems:
```bash
sudo apt-get install libboost-all-dev
```
For MacOS, you can install it with the following command (assuming you have homebrew installed):
```bash
brew install boost
```

## Compilation
To compile the compiler, run the following command:
```bash
make all
```
The executable will be compiled as `sim` in the base directory.

## Usage
To run the compiler, use the following command:
```bash
./sim <config_file>
```

## Configuration File
The configuration file is used to pass the parameters to the compiler and the simulator. The configuration file has the following format:
```
<parameter_name> = <parameter_value>
```

### Supported Parameters and Descriptions
- `input_dir`: The input directory containing the circuits to be executed. (Default: `input_dir`)
- `input_file`: Regular expression for the circuits to be executed in the input directory. (Default: `gate_.*`)
- `output_dir`: The output file where the simulation results will be stored. (Default: `output_dir`)
- `code_distance`: The code distance of the surface code.
- `physical_qubit_error_rate`: The error rate of the physical qubits.
- `debug`: If set to `1`, the compiler will print debug information. (Default: `0`)
- `number_of_runs`: The number of times to run each circuit (on different seeds). (Default: `100`)
- `compiler`: The compiler to use. Supported values are `dynamic`, `static` and `autobraid`. (Default: `dynamic`)
- `compression_factor`: The compression factor for the surface code grid. Supported range is between `0-1`. (Default: `0`)
- `mst_computation_frequency`: The frequency (in cycles) of the MST computation. (Default: `100`)
