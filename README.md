# svxyz

`svxyz` is a versatile suite of command-line tools designed to handle and process atomic configuration data for scientific and computational tasks. It provides an intuitive interface for tasks such as filtering, extracting, visualizing, analyzing symmetry, and converting file formats.

---

## Features

### **`txyz`**
Processes `vasprun.xml/OUTCAR` files to filter atomic configurations based on:
- Energy
- Force
- Stress
- Frame ranges
- Volume
- Pressure
- Virial tensors

Supports flexible JSON configurations (`txyz.json`) for custom workflows. Automatically generates configurations if not present.

---

### **`dxyz`**
Extracts data from `.xyz` files, including:
- Energy
- Atomic forces
- Stress tensors
- Volume, pressure, and temperature information

Outputs separate `.dat` files (`E.dat`, `F.dat`, etc.) for detailed analysis.

---

### **`pxyz`**
Visualizes data from `.dat` or `infile.dat` formats, featuring:
- KDE distributions
- Projections
- Customizable plotting layouts

Supports JSON-based configurations for plot styling and layout (`pxyz.json`).

---

### **`asefmt`**
Converts file formats for atomic structures using `ASE`, `Pymatgen`, and `Ovito`. Features include:
- Supports a wide range of formats like `POSCAR`, `CIF`, `XYZ`, `LAMMPS-data`, and more.
- Configurable JSON file (`ccfmt.json`) for default settings.
- Customizable replication and charge settings for specific formats.

---

### **`analpos`**
Analyzes atomic configurations for:
- Pairwise atomic distances
- Classification by element pairs
- Symmetry detection using tolerance
- Crystal system, space group, and point group information

Outputs results to distance and symmetry files (`distance_POSCAR.dat`, `sym_POSCAR.dat`).

---

### **`xyz2pos`**
Extracts specific configurations from `.xyz` files and writes them as `POSCAR` files for further processing.

---

## Installation

### From Source
Clone the repository and install the package:
```bash
git clone https://github.com/Superionichuan/svxyz.git
cd svxyz
pip install .

