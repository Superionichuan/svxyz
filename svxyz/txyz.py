import os
import re
import json
import numpy as np
from glob import glob
from ase.io import read, write

# 常量: 从 eV/Å³ 转换为 GPa
EV_A3_TO_GPA = 160.21766208

# 默认配置文件内容
DEFAULT_CONFIG = {
    "input_files": ["vasprun.xml"],
    "input_format": "vasp-xml",
    "output_file": "filtered_output.xyz",
    "skip": {"on": 1, "count": 500},
    "frame_range": [None, None],
    "energy_range": [None, None],
    "max_atomic_force_range": [None, None],
    "pressure_range": [None, None],
    "volume_range": [None, None],
    "temperature_range": [None, None],
    "virial_filters": {
        "V_xx": [None, None],
        "V_yy": [None, None],
        "V_zz": [None, None],
        "V_yz": [None, None],
        "V_xz": [None, None],
        "V_xy": [None, None]
    },
    "stress_filters": {
        "S_xx": [None, None],
        "S_yy": [None, None],
        "S_zz": [None, None],
        "S_yz": [None, None],
        "S_xz": [None, None],
        "S_xy": [None, None]
    },
    "atomic_filters": None,
    "min_distance": None,  # 最小距离 (Å)
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
        exit(0)
    with open(file_name, "r") as f:
        return json.load(f)

# 从文件中提取温度
def extract_temperatures_from_file(file_path, pattern):
    temperatures = []
    if os.path.exists(file_path):
        with open(file_path, "r") as file:
            for line in file:
                match = re.search(pattern, line)
                if match:
                    temperatures.append(float(match.group(1)))
    return temperatures

# 获取温度信息
def get_temperatures(tb_path, outcar_path, n_frames):
    tb_pattern = r"([\d\.]+)"  # TB.dat 格式：每行一个温度值
    outcar_pattern = r"kin\. lattice\s+EKIN_LAT=.*\(temperature\s+([\d\.]+)\s+K\)"  # OUTCAR 格式

    temperatures = extract_temperatures_from_file(tb_path, tb_pattern)
    if not temperatures:
        temperatures = extract_temperatures_from_file(outcar_path, outcar_pattern)
    if not temperatures:
        temperatures = [0.0] * n_frames
    return temperatures

# 从 ST.dat 提取应力
def extract_stress_from_file(file_path):
    stresses = []
    if os.path.exists(file_path):
        with open(file_path, "r") as file:
            for line in file:
                tokens = line.split()
                if len(tokens) >= 7:
                    try:
                        stress_values = list(map(float, tokens[1:7]))
                        reordered_stress = [
                            stress_values[0],  # XX
                            stress_values[1],  # YY
                            stress_values[2],  # ZZ
                            stress_values[4],  # YZ
                            stress_values[5],  # XZ
                            stress_values[3],  # XY
                        ]
                        stresses.append(np.array(reordered_stress))  # 单位为 GPa
                    except ValueError:
                        continue
    return stresses

# 从 OUTCAR 提取应力
def extract_fstress_from_outcar(outcar_path):
    stresses = []
    if os.path.exists(outcar_path):
        with open(outcar_path, "r") as f:
            for line in f:
                if "Total+kin." in line:
                    stress_values = list(map(float, line.split()[1:7]))
                    reordered_stress = [
                        stress_values[0],  # XX
                        stress_values[1],  # YY
                        stress_values[2],  # ZZ
                        stress_values[4],  # YZ
                        stress_values[5],  # XZ
                        stress_values[3],  # XY
                    ]
                    stresses.append(np.array(reordered_stress) / 10.0)  # 转换为 GPa
    return stresses

# 获取应力信息
def get_fstress(st_path, outcar_path, atoms_list):
    if os.path.exists(st_path):
        stresses = extract_stress_from_file(st_path)
    elif os.path.exists(outcar_path):
        stresses = extract_fstress_from_outcar(outcar_path)
    else:
        stresses = [
            atoms.get_stress(voigt=True) * -EV_A3_TO_GPA for atoms in atoms_list
        ]
    return stresses

# 计算压力 (GPa)
def calculate_pressure(stress):
    return np.mean(stress[:3])  # 压力为主应力的平均值


# 计算最小原子间距及其原子对
def calculate_min_distance(atoms):
    positions = atoms.get_positions()
    n_atoms = len(positions)
    min_distance = float("inf")
    min_pair = None
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            distance = np.linalg.norm(positions[i] - positions[j])
            if distance < min_distance:
                min_distance = distance
                min_pair = (i+1, j+1)
    return min_distance, min_pair

# 筛选并添加温度、压力、体积、应力和 virial 信息
def filter_atoms(atoms_list, stresses, temperatures, config):
    skip_count = config["skip"]["count"] if config["skip"].get("on", 0) else 0
    frame_range = config["frame_range"]
    energy_range = config["energy_range"]
    max_atomic_force_range = config["max_atomic_force_range"]
    pressure_range = config["pressure_range"]
    volume_range = config["volume_range"]
    temperature_range = config["temperature_range"]
    min_distance_threshold = config.get("min_distance")

    filtered_atoms = []
    min_distance_violations = []
    for i, atoms in enumerate(atoms_list):
        if i < skip_count:
            continue

        if frame_range[0] is not None and i < frame_range[0]:
            continue
        if frame_range[1] is not None and i > frame_range[1]:
            continue

        energy = atoms.get_potential_energy()
        if energy_range[0] is not None and energy < energy_range[0]:
            continue
        if energy_range[1] is not None and energy > energy_range[1]:
            continue

        max_force = atoms.get_forces().max()
        if max_atomic_force_range[0] is not None and max_force < max_atomic_force_range[0]:
            continue
        if max_atomic_force_range[1] is not None and max_force > max_atomic_force_range[1]:
            continue

        volume = atoms.get_volume()
        if volume_range[0] is not None and volume < volume_range[0]:
            continue
        if volume_range[1] is not None and volume > volume_range[1]:
            continue

        stress = stresses[i] if stresses and i < len(stresses) else None
        pressure = calculate_pressure(stress) if stress is not None else None
        if pressure_range[0] is not None and pressure < pressure_range[0]:
            continue
        if pressure_range[1] is not None and pressure > pressure_range[1]:
            continue

        temperature = temperatures[i]
        if temperature_range[0] is not None and temperature < temperature_range[0]:
            continue
        if temperature_range[1] is not None and temperature > temperature_range[1]:
            continue

        min_distance, min_pair = calculate_min_distance(atoms)
        if min_distance_threshold is not None and min_distance < min_distance_threshold:
            min_distance_violations.append(f"Frame {i}, Atoms {min_pair}, Distance {min_distance:.4f}")
            continue

        atoms.info["temperature"] = f"{temperature:.2f}"
        atoms.info["volume"] = f"{volume:.4f}"
        atoms.info["pressure"] = f"{pressure:.4f}" if pressure is not None else "N/A"
        atoms.info["mindistance"] = f"{min_distance:.4f} {min_pair[0]}-{min_pair[1]}"
        if stress is not None:
            atoms.info["fstress"] = ", ".join(f"{s:.4f}" for s in stress)

        filtered_atoms.append(atoms)

    return filtered_atoms, skip_count, min_distance_violations

# 主函数
def main():
    config = load_or_create_config()

    input_files = config.get("input_files", ["vasprun.xml"])
    input_format = config.get("input_format", "vasp-xml")
    output_file = config.get("output_file", "filtered_output.xyz")

    all_files = []
    for pattern in input_files:
        all_files.extend(glob(pattern))

    if not all_files:
        print(f"Error: No input files found matching patterns: {input_files}")
        return

    filtered_atoms_total = []
    total_frames = 0
    min_distance_violations = []
    for input_file in all_files:
        print(f"Reading frames from {input_file} with format {input_format}...")
        try:
            atoms_list = read(input_file, index=":", format=input_format)
        except Exception as e:
            print(f"Error reading file {input_file} with format {input_format}: {e}")
            continue

        total_frames += len(atoms_list)

        tb_path = os.path.join(os.path.dirname(input_file), "TB.dat")
        outcar_path = os.path.join(os.path.dirname(input_file), "OUTCAR")
        st_path = os.path.join(os.path.dirname(input_file), "ST.dat")

        temperatures = get_temperatures(tb_path, outcar_path, len(atoms_list))
        stresses = get_fstress(st_path, outcar_path, atoms_list)

        print(f"Filtering frames from {input_file}...")
        filtered_atoms, skip_count, violations = filter_atoms(atoms_list, stresses, temperatures, config)
        filtered_atoms_total.extend(filtered_atoms)
        min_distance_violations.extend(violations)

    write(output_file, filtered_atoms_total)
    print(f"Filtered {len(filtered_atoms_total)} frames saved to {output_file}")

    if config["show_summary"]:
        print(f"\nSummary:")
        print(f"  Total frames: {total_frames}")
        print(f"  Skipped frames: {skip_count}")
        print(f"  Filtered frames: {len(filtered_atoms_total)}")
        if min_distance_violations:
            print("\nMinimum distance violations:")
            for violation in min_distance_violations:
                print(f"  {violation}")


if __name__ == "__main__":
    main()

