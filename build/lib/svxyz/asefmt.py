import argparse
import json
import os
from ase.io import read as ase_read, write as ase_write

# 默认配置
DEFAULT_CONFIG = {
    "input_file": "POSCAR",
    "output_file": "output.cif",
    "input_format": None,  # ASE 自动识别
    "output_format": None,  # ASE 自动识别
    "_comment": {
        "input_format": "支持的输入文件格式包括但不限于：POSCAR, CIF, XYZ, PDB, OUTCAR, JSON, VASPRUN.xml, EXTXYZ, LAMMPS-data, LAMMPS-dump, ARC, GEN, XSF, CUBE, TRAJ, GPW, FDF, GRO, PRMTOP, GIN, MOP 等。",
        "output_format": "支持的输出文件格式包括但不限于：POSCAR, CIF, XYZ, PDB, OUTCAR, JSON, VASPRUN.xml, EXTXYZ, LAMMPS-data, LAMMPS-dump, ARC, GEN, XSF, CUBE, TRAJ, POV, GPW, FDF, GRO, PRMTOP, GIN, MOP 等。",
        "xyz": "XYZ 格式，常见的轻量级分子结构格式",
        "vasp": "VASP 的 POSCAR/CONTCAR 格式",
        "cif": "Crystallographic Information File，晶体信息文件",
        "pdb": "蛋白质结构文件，广泛应用于生物化学",
        "xsf": "XCrySDen 结构文件，用于晶体和分子建模",
        "cube": "Gaussian 立方网格格式，用于量子化学分析",
        "traj": "ASE 的轨迹格式，用于保存模拟轨迹",
        "lammps-data": "LAMMPS 的 Data 格式，用于分子动力学模拟",
        "lammps-dump": "LAMMPS 的 Dump 格式，用于轨迹文件保存",
        "outcar": "VASP 的 OUTCAR 文件，包含模拟结果详细信息",
        "json": "通用 JSON 格式，用于结构和属性描述",
        "gpw": "GPAW 的二进制文件格式，用于电子结构计算",
        "fdf": "FHI-aims 或 Siesta 的输入文件格式",
        "gen": "GULP 或 DFTB+ 的 GEN 文件，用于定义系统结构",
        "arc": "Materials Studio 的 Arc 格式，用于分子模拟",
        "extxyz": "扩展的 XYZ 格式，支持原子属性",
        "prmtop": "AMBER 的拓扑文件格式，用于分子动力学模拟",
        "gro": "GROMACS 的结构文件格式",
        "gin": "GULP 的输入文件",
        "mop": "MOPAC 的输入文件格式",
        "pov": "POV-Ray 格式，用于图形渲染",
        "vasprun.xml": "VASP 的 Vasprun 文件，用于存储模拟详细结果",
        "xdatcar": "VASP 的 XDATCAR 文件，用于轨迹保存",
        "doscar": "VASP 的 DOSCAR 文件，用于密度状态分析",
        "eigenval": "VASP 的 EIGENVAL 文件，用于电子结构分析"
    }
}

def create_default_config(file_name="asefmt.json"):
    """生成默认 JSON 配置文件"""
    with open(file_name, "w", encoding="utf-8") as f:
        json.dump(DEFAULT_CONFIG, f, indent=4, ensure_ascii=False)
    print(f"Default configuration created: {file_name}")

def load_or_create_config(file_name="asefmt.json"):
    """加载配置文件，如果不存在则生成并退出"""
    if not os.path.exists(file_name):
        create_default_config(file_name)
        print(f"JSON configuration file not found. Generated a default configuration. Exiting now.")
        exit(0)
    with open(file_name, "r", encoding="utf-8") as f:
        return json.load(f)

def convert_with_ase(input_file, output_file, input_format=None, output_format=None):
    """使用 ASE 转换文件格式"""
    try:
        atoms = ase_read(input_file, format=input_format)
        ase_write(output_file, atoms, format=output_format)
        print(f"Converted {input_file} -> {output_file} using ASE successfully!")
    except Exception as e:
        print(f"Error with ASE: {e}")

def main():
    # 命令行参数解析
    parser = argparse.ArgumentParser(description="Convert material structure files using ASE.")
    parser.add_argument("--i", required=True, help="Input file path")
    parser.add_argument("--o", required=True, help="Output file path")
    parser.add_argument("--ifmt", default=None, help="Input file format (optional)")
    parser.add_argument("--ofmt", default=None, help="Output file format (optional)")
    args = parser.parse_args()

    # 加载配置
    config_file = "asefmt.json"
    config = load_or_create_config(config_file)

    # 使用命令行参数覆盖配置
    input_file = args.i or config["input_file"]
    output_file = args.o or config["output_file"]
    input_format = args.ifmt or config["input_format"]
    output_format = args.ofmt or config["output_format"]

    # 转换文件格式
    convert_with_ase(input_file, output_file, input_format, output_format)

if __name__ == "__main__":
    main()

