# svxyz

`svxyz` is a suite of command-line tools designed to handle and process XYZ data for scientific and computational tasks. This package provides three powerful scripts: `txyz`, `dxyz`, and `pxyz`, each tailored for specific use cases like filtering, visualization, and data extraction.

---

## Features

- **`txyz`**:  
  Processes `vasprun.xml/OUTCAR` files to filter atomic configurations based on energy, force, stress, and frame ranges. Supports flexible JSON configuration for custom workflows.

- **`xyz2pos`**:
  Extract specific configuration POSCAR* from `.xyz` file.
  
- **`dxyz`**:  
  Extracts energy, forces, and stresses from XYZ files and writes them to separate `.dat` files for further analysis.

- **`pxyz`**:  
  Visualizes data from `.dat` or `infile.dat` formats with dynamic plotting capabilities, including KDE distributions and projections.

---

## Installation

### From Source

Clone the repository and install the package:

```bash
git clone https://github.com/yourusername/svxyz.git
cd svxyz
pip install .

