import os
import json
import numpy as np
from ase.io import read, write

# 常量: 从 eV/Å³ 转换为 GPa
EV_A3_TO_GPA = 160.21766208

# 默认配置文件内容
DEFAULT_CONFIG = {
    "input_file": "vasprun.xml",      # 默认输入文件名
    "input_format": "vasp-xml",      # 默认输入文件格式
    "output_file": "filtered_output.xyz",  # 默认输出文件名
    "skip": {
        "on": 1,
        "count": 500
    },
    "frame_range": [
        None,
        None
    ],
    "energy_range": [
        None,
        None
    ],
    "max_atomic_force_range": [
        None,
        None
    ],
    "stress_filters": {
        "S_xx": [
            None,
            None
        ],
        "S_yy": [
            None,
            None
        ],
        "S_zz": [
            None,
            None
        ],
        "S_yz": [
            None,
            None
        ],
        "S_xz": [
            None,
            None
        ],
        "S_xy": [
            None,
            None
        ]
    },
    "atomic_filters": ["O", "Fe"],  # 仅筛选含氧和铁的帧
    "stress_unit": "GPa",
    "show_summary": True
}

def create_default_config(file_name="txyz.json"):
    """生成默认的 JSON 配置文件"""
    with open(file_name, "w") as f:
        json.dump(DEFAULT_CONFIG, f, indent=4)
    print(f"Default configuration created: {file_name}")

def load_or_create_config(file_name="txyz.json"):
    """加载配置文件，如果不存在则生成并退出"""
    if not os.path.exists(file_name):
        create_default_config(file_name)
        print("JSON configuration file not found. Generated a default configuration. Exiting now.")
        exit(0)  # 生成 JSON 文件后退出程序
    with open(file_name, "r") as f:
        return json.load(f)

def filter_atoms(atoms_list, config):
    """按配置筛选帧"""
    skip_count = config["skip"]["count"] if config["skip"]["on"] else 0
    frame_range = config["frame_range"]
    energy_range = config["energy_range"]
    max_atomic_force_range = config["max_atomic_force_range"]
    stress_filters = config["stress_filters"]
    atomic_filters = config.get("atomic_filters", None)  # 获取原子种类筛选

    filtered_atoms = []
    for i, atoms in enumerate(atoms_list):
        # 跳过前 skip_count 帧
        if i < skip_count:
            continue

        # 按帧范围过滤
        if frame_range[0] is not None and i < frame_range[0]:
            continue
        if frame_range[1] is not None and i > frame_range[1]:
            continue

        # 获取帧的能量
        energy = atoms.get_potential_energy()

        # 判断能量下限
        if energy_range[0] is not None and energy < energy_range[0]:
            continue

        # 判断能量上限
        if energy_range[1] is not None and energy > energy_range[1]:
            continue

        # 获取帧的力
        max_force = atoms.get_forces().max()
        if max_atomic_force_range[0] is not None and max_force < max_atomic_force_range[0]:
            continue
        if max_atomic_force_range[1] is not None and max_force > max_atomic_force_range[1]:
            continue

        # 获取帧的应力
        stress = atoms.get_stress(voigt=True) * EV_A3_TO_GPA  # 转换为 GPa
        stress_filters_passed = True
        for key, (lower, upper) in stress_filters.items():
            idx = ["S_xx", "S_yy", "S_zz", "S_yz", "S_xz", "S_xy"].index(key)
            value = -stress[idx]  # 取负值
            if lower is not None and value < lower:
                stress_filters_passed = False
                break
            if upper is not None and value > upper:
                stress_filters_passed = False
                break

        if not stress_filters_passed:
            continue

        # 按原子种类筛选
        if atomic_filters:
            symbols = atoms.get_chemical_symbols()
            if not any(el in atomic_filters for el in symbols):
                continue

        # 如果所有条件通过，添加到结果中
        filtered_atoms.append(atoms)

    return filtered_atoms

def main():
    # 加载或生成配置文件
    config = load_or_create_config()

    # 提取输入和输出文件信息
    input_file = config.get("input_file", "vasprun.xml")
    input_format = config.get("input_format", "vasp-xml")
    output_file = config.get("output_file", "filtered_output.xyz")

    # 读取输入数据
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found.")
        return

    print(f"Reading frames from {input_file} with format {input_format}...")
    atoms_list = read(input_file, index=':', format=input_format)

    # 筛选帧
    print("Filtering frames based on configuration...")
    filtered_atoms = filter_atoms(atoms_list, config)

    # 保存筛选结果
    write(output_file, filtered_atoms)
    print(f"Filtered {len(filtered_atoms)} frames saved to {output_file}")

    # 打印摘要
    if config["show_summary"]:
        print(f"\nSummary:")
        print(f"  Total frames: {len(atoms_list)}")
        print(f"  Skipped frames: {config['skip']['count']}")
        print(f"  Filtered frames: {len(filtered_atoms)}")

if __name__ == "__main__":
    main()

