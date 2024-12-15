import argparse
import os
import json
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# 默认配置
DEFAULT_CONFIG = {
    "tol": 0.01,  # 对称性分析容差
    "input_file": "POSCAR"
}

def load_or_create_config(config_file="analpos.json"):
    """加载配置文件，如果不存在则生成默认配置"""
    if not os.path.exists(config_file):
        with open(config_file, "w") as f:
            json.dump(DEFAULT_CONFIG, f, indent=4)
        print(f"Default configuration created: {config_file}")
        return DEFAULT_CONFIG
    else:
        with open(config_file, "r") as f:
            return json.load(f)

def save_config(config, config_file="analpos.json"):
    """保存更新后的配置到 JSON 文件"""
    with open(config_file, "w") as f:
        json.dump(config, f, indent=4)
    print(f"Configuration updated: {config_file}")

def read_poscar(file_path):
    """读取 POSCAR 文件，返回 pymatgen 的 Structure 对象"""
    try:
        structure = Structure.from_file(file_path)
    except Exception as e:
        raise ValueError(f"Failed to read POSCAR file: {e}")
    return structure

def calculate_distances(structure):
    """计算所有原子对的距离，并按元素对分类"""
    positions = structure.cart_coords
    species = structure.species
    n_atoms = len(positions)
    distances = {}

    # 遍历所有原子对
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):  # 只考虑 i < j
            distance = np.linalg.norm(positions[i] - positions[j])
            pair = f"{species[i]}-{species[j]}"
            atom_pair = f"{species[i]}{i+1}-{species[j]}{j+1}"
            if pair not in distances:
                distances[pair] = []
            distances[pair].append((distance, atom_pair))  # 添加距离和原子对信息

    # 对每一类距离排序
    for pair in distances:
        distances[pair].sort(key=lambda x: x[0])

    return distances
def save_distances(distances, output_file):
    """将分类后的最小距离和原子对保存到文件"""
    with open(output_file, 'w') as f:
        # 写入标题行
        f.write(f"{'Rank':<8}{'Pair':<12}{'Min Distance (Å)':<20}{'Atom Pair':<20}\n")
        f.write("=" * 60 + "\n")
        summary = summarize_distances(distances)
        for rank, (pair, stats) in enumerate(summary, start=1):
            f.write(f"{rank:<8}{pair:<12}{stats['min']:<20.4f}{stats['atom_pair']:<20}\n")



def summarize_distances(distances):
    """生成每类距离的最小值及对应的原子对"""
    summary = []
    for key, value in distances.items():
        min_distance = value[0][0]
        min_atom_pair = value[0][1]
        summary.append((key, {"min": min_distance, "atom_pair": min_atom_pair}))
    # 按最小距离排序
    return sorted(summary, key=lambda x: x[1]["min"])

def analyze_symmetry(structure, tol):
    """使用 pymatgen 分析晶体对称性"""
    try:
        analyzer = SpacegroupAnalyzer(structure, symprec=tol)
        space_group_symbol = analyzer.get_space_group_symbol()
        space_group_number = analyzer.get_space_group_number()
        crystal_system = analyzer.get_crystal_system()
        point_group = analyzer.get_point_group_symbol()

        return {
            "space_group_symbol": space_group_symbol,
            "space_group_number": space_group_number,
            "crystal_system": crystal_system,
            "point_group": point_group
        }
    except Exception as e:
        raise ValueError(f"Failed to analyze symmetry: {e}")

def save_symmetry(symmetry_info, distances, output_file):
    """保存对称性和最小距离统计到文件"""
    with open(output_file, 'w') as f:
        # 写入对称性信息
        f.write("Symmetry Information:\n")
        f.write(f"space_group_symbol: {symmetry_info['space_group_symbol']}\n")
        f.write(f"space_group_number: {symmetry_info['space_group_number']}\n")
        f.write(f"crystal_system: {symmetry_info['crystal_system']}\n")
        f.write(f"point_group: {symmetry_info['point_group']}\n\n")

        # 写入最小距离统计
        f.write(f"{'Rank':<8}{'Pair':<12}{'Min Distance (Å)':<20}{'Atom Pair':<20}\n")
        f.write("=" * 60 + "\n")
        summary = summarize_distances(distances)
        for rank, (pair, stats) in enumerate(summary, start=1):
            f.write(f"{rank:<8}{pair:<12}{stats['min']:<20.4f}{stats['atom_pair']:<20}\n")

def display_summary(symmetry_info, distances):
    """在屏幕上显示对称性和最小距离统计"""
    print("\nSymmetry Information:")
    print(f"  space_group_symbol: {symmetry_info['space_group_symbol']}")
    print(f"  space_group_number: {symmetry_info['space_group_number']}")
    print(f"  crystal_system: {symmetry_info['crystal_system']}")
    print(f"  point_group: {symmetry_info['point_group']}")

    print("\nTotal Distances:")
    summary = summarize_distances(distances)
    for rank, (pair, stats) in enumerate(summary, start=1):
        print(f"  {rank:<3} {pair:<12} Min: {stats['min']:<10.4f} Pair: {stats['atom_pair']}")

def main():
    parser = argparse.ArgumentParser(description="Analyze POSCAR for distances and symmetry.")
    parser.add_argument("--f", help="Input POSCAR file")
    parser.add_argument("--tol", type=float, help="Tolerance for symmetry analysis")
    args = parser.parse_args()

    # 加载或生成配置文件
    config = load_or_create_config()

    # 覆盖配置文件参数
    if args.f:
        config["input_file"] = args.f
    if args.tol is not None:
        config["tol"] = args.tol

    # 保存更新后的配置
    save_config(config)

    input_file = config["input_file"]
    tol = config["tol"]

    # 文件名定义
    base_name = os.path.basename(input_file)
    distance_file = f"distance_{base_name}.dat"
    symmetry_file = f"sym_{base_name}.dat"

    print(f"Reading POSCAR file: {input_file}...")
    structure = read_poscar(input_file)

    print("Calculating distances...")
    distances = calculate_distances(structure)

    print(f"Saving distances to {distance_file}...")
    save_distances(distances, distance_file)

    print("Analyzing symmetry...")
    symmetry_info = analyze_symmetry(structure, tol)

    print(f"Saving symmetry information to {symmetry_file}...")
    save_symmetry(symmetry_info, distances, symmetry_file)

    print("Displaying summary...")
    display_summary(symmetry_info, distances)

    print("Done!")

if __name__ == "__main__":
    main()

