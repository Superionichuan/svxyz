import sys
import os
import json
import argparse
from ase.io import read
import numpy as np

def extract_properties_with_ids(xyz_file):
    """
    从包含信息的 XYZ 文件中提取能量、力、virial、应力、温度、压力和体积数据，
    分别写入 E.dat, F.dat, virial.dat, stress.dat 和 extra_info.dat。

    Args:
        xyz_file (str): 输入的 XYZ 文件路径。
    """
    try:
        # 加载所有帧
        atoms_list = read(xyz_file, index=":")

        # 检查是否存在 `mindistance` 信息
        has_mindistance = any("mindistance" in atoms.info for atoms in atoms_list)

        # 打开输出文件
        with open("E.dat", "w") as e_file, \
             open("F.dat", "w") as f_file, \
             open("virial.dat", "w") as vir_file, \
             open("stress.dat", "w") as stress_file, \
             open("extra_info.dat", "w") as extra_file:

            # 写入标题行
            e_file.write(f"{'Energy(eV)':<20}{'System_ID':<10}\n")
            f_file.write(f"{'Max_Force':<20}{'Mean_Force':<20}{'System_ID':<10}\n")
            vir_file.write(f"{'V_xx':<18}{'V_yy':<18}{'V_zz':<18}{'V_yz':<18}{'V_xz':<18}{'V_xy':<18}{'System_ID':<10}\n")
            stress_file.write(f"{'S_xx':<18}{'S_yy':<18}{'S_zz':<18}{'S_yz':<18}{'S_xz':<18}{'S_xy':<18}{'System_ID':<10}\n")

            # 根据是否有 `mindistance` 信息写入标题行
            if has_mindistance:
                extra_file.write(f"{'Temperature(K)':<16}{'Pressure(GPa)':<15}{'Volume(Å³)':<15}{'Min_Distance(Å)':<18}{'Atom_Pair':<15}{'System_ID':<10}\n")
            else:
                extra_file.write(f"{'Temperature(K)':<16}{'Pressure(GPa)':<15}{'Volume(Å³)':<15}{'System_ID':<10}\n")

            # 遍历每一帧
            for i, atoms in enumerate(atoms_list):
                # 提取能量
                energy = atoms.get_potential_energy()
                if energy is None:
                    raise ValueError(f"Frame {i} is missing 'energy' in atoms.info.")
                e_file.write(f"{energy:<20.6f}{i:<10}\n")

                # 提取力范数
                forces = atoms.get_forces()
                mean_force_norm = np.linalg.norm(forces, axis=1).mean()  # 平均力范数
                max_force_norm = np.linalg.norm(forces, axis=1).max()  # 最大力范数
                f_file.write(f"{max_force_norm:<20.6f}{mean_force_norm:<20.6f}{i:<10}\n")

                # 提取应力
                stress = atoms.info.get("fstress", None)
                if stress is None:
                    raise ValueError(f"Frame {i} is missing 'fstress' in atoms.info.")

                # 处理应力（fstress）
                if isinstance(stress, str):
                    stress_values = list(map(float, stress.split(", ")))
                elif isinstance(stress, (list, np.ndarray)):
                    stress_values = list(map(float, stress))
                else:
                    raise TypeError(f"Unsupported type for 'fstress': {type(stress)} in frame {i}")

                # 写入应力（fstress）
                stress_file.write(f"{stress_values[0]:<18.6f}{stress_values[1]:<18.6f}{stress_values[2]:<18.6f}"
                                  f"{stress_values[3]:<18.6f}{stress_values[4]:<18.6f}{stress_values[5]:<18.6f}{i:<10}\n")

                # 计算 virial
                stress_voigt = atoms.get_stress(voigt=True)  # ASE 应力（不含电子动能项）
                volume = atoms.get_volume()
                virial = [-s * volume for s in stress_voigt]  # Virial = -stress * volume

                # 写入 virial
                vir_file.write(f"{virial[0]:<18.6f}{virial[1]:<18.6f}{virial[2]:<18.6f}"
                               f"{virial[3]:<18.6f}{virial[4]:<18.6f}{virial[5]:<18.6f} {i:<10}\n")

                # 提取额外信息：温度、压力和体积
                temperature = atoms.info.get("temperature", "N/A")
                pressure = atoms.info.get("pressure", "N/A")

                # 写入 `extra_info.dat`，根据是否有 `mindistance` 动态调整列
                if has_mindistance:
                    mindistance = atoms.info.get("mindistance", "N/A")
                    min_pair = atoms.info.get("min_pair", "N/A")
                    extra_file.write(f"{temperature:<16}{pressure:<15}{volume:<15.6f}{mindistance:<18}{min_pair:<15}{i:<10}\n")
                else:
                    extra_file.write(f"{temperature:<16}{pressure:<15}{volume:<15.6f}{i:<10}\n")

        print("数据已成功写入 E.dat, F.dat, virial.dat, stress.dat 和 extra_info.dat。")

    except Exception as e:
        print(f"发生错误: {e}")

def main():
    parser = argparse.ArgumentParser(
        description=(
            "dxyz: A tool to extract energy, forces, virial, stress, and additional properties from XYZ files.\n\n"
            "Usage:\n"
            "  dxyz [INPUT_FILE] \n  Note: If the input_file is not specified in the command line, the script will attempt to load it from txyz.json. If txyz.json is not found, a default txyz.json will be generated. Command line arguments will override settings in the JSON file.\n\n"
            "Options:\n"
            "  -h, --help        Show this help message and exit.\n\n"
            "Description:\n"
            "  This tool processes an XYZ file with embedded property information (e.g., energy, forces, virial, stress,\n"
            "  temperature, pressure, volume). It writes extracted data into separate `.dat` files:\n"
            "    - E.dat: Energy and system ID.\n"
            "    - F.dat: Maximum and mean atomic force norms with system ID.\n"
            "    - virial.dat: Virial stress tensor components with system ID.\n"
            "    - stress.dat: Stress tensor components with system ID.\n"
            "    - extra_info.dat: Additional properties such as temperature, pressure, volume, and optionally, the\n"
            "                      minimum interatomic distance and associated atomic pair.\n\n"
            "Example:\n"
            "  dxyz sample.xyz\n"
            "  This will process 'sample.xyz' and generate the output `.dat` files in the current directory.\n\n"
            "Notes:\n"
            "  If the input file contains 'mindistance' information in the `info` attribute of each frame,\n"
            "  this tool will include additional columns for minimum interatomic distances and atomic pairs\n"
            "  in the `extra_info.dat` file.\n\n"
            "JSON Configuration:\n"
            "  The tool supports optional configuration via a `dxyz.json` file. If the file is absent,\n"
            "  it will be created with default values. The input XYZ file can also be specified in the JSON file.\n\n"
            "Configuration Keys:\n"
            "  - xyz_file: Path to the input XYZ file."
        ),
        formatter_class=argparse.RawTextHelpFormatter  # 保证换行符正常显示
    )
    
    parser.add_argument(
        "input_file",
        help="Input XYZ file path."
    )   
        
    args = parser.parse_args()
        
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1) 

    # 默认 JSON 文件名
    config_file = "dxyz.json"

    # 如果有命令行参数，则保存到 JSON 文件
    if len(sys.argv) > 1:
        xyz_file = sys.argv[1]

        # 保存参数到 JSON 文件
        with open(config_file, "w") as f:
            json.dump({"xyz_file": xyz_file}, f, indent=4)
        print(f"参数已保存到 {config_file}。")

        # 执行提取
        extract_properties_with_ids(xyz_file)

    # 如果没有命令行参数，则读取 JSON 文件
    else:
        if os.path.exists(config_file):
            with open(config_file, "r") as f:
                config = json.load(f)
                xyz_file = config.get("xyz_file")

                if xyz_file:
                    # 执行提取
                    extract_properties_with_ids(xyz_file)
                else:
                    print("配置文件中没有找到 'xyz_file' 参数。")
        else:
            print(f"请提供输入文件名，或确保 {config_file} 存在并包含正确参数。")

if __name__ == "__main__":
    main()

