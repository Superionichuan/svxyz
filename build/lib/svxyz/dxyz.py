import sys
import os
import json
from ase.io import read
import numpy as np

def extract_properties_with_ids(xyz_file):
    """
    从包含信息的 XYZ 文件中提取能量、力、virial、应力、温度、压力、体积和最小距离数据，
    分别写入 E.dat, F.dat, virial.dat, stress.dat 和 extra_info.dat。

    Args:
        xyz_file (str): 输入的 XYZ 文件路径。
    """
    try:
        # 加载所有帧
        atoms_list = read(xyz_file, index=":")

        # 打开输出文件
        with open("E.dat", "w") as e_file, \
             open("F.dat", "w") as f_file, \
             open("virial.dat", "w") as vir_file, \
             open("stress.dat", "w") as stress_file, \
             open("extra_info.dat", "w") as extra_file:

            # 写入标题行
            e_file.write("Energy(eV) System_ID\n")
            f_file.write("Max_atomic_force_norm Mean_atomic_force_norm System_ID\n")
            vir_file.write("V_xx V_yy V_zz V_yz V_xz V_xy System_ID\n")
            stress_file.write("S_xx S_yy S_zz S_yz S_xz S_xy System_ID\n")
            extra_file.write("Temperature(K) Pressure(GPa) Volume(Å³) Min_Distance(Å) Atom_Pair System_ID\n")

            # 遍历每一帧
            for i, atoms in enumerate(atoms_list):
                # 提取能量
                energy = atoms.info.get("energy", None)
                if energy is None:
                    raise ValueError(f"Frame {i} is missing 'energy' in atoms.info.")
                e_file.write(f"{energy:.6f} {i}\n")

                # 提取力范数
                forces = atoms.get_forces()
                mean_force_norm = np.linalg.norm(forces, axis=1).mean()  # 平均力范数
                max_force_norm = np.linalg.norm(forces, axis=1).max()  # 最大力范数
                f_file.write(f"{max_force_norm:.6f} {mean_force_norm:.6f} {i}\n")

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
                stress_file.write(f"{stress_values[0]:.6f} {stress_values[1]:.6f} {stress_values[2]:.6f} "
                                  f"{stress_values[3]:.6f} {stress_values[4]:.6f} {stress_values[5]:.6f} {i}\n")

                # 计算 virial
                stress_voigt = atoms.get_stress(voigt=True)  # ASE 应力（不含电子动能项）
                volume = atoms.get_volume()
                virial = [-s * volume for s in stress_voigt]  # Virial = -stress * volume

                # 写入 virial
                vir_file.write(f"{virial[0]:.6f} {virial[1]:.6f} {virial[2]:.6f} "
                               f"{virial[3]:.6f} {virial[4]:.6f} {virial[5]:.6f} {i}\n")

                # 提取额外信息：温度、压力、体积和最小距离
                temperature = atoms.info.get("temperature", "N/A")
                pressure = atoms.info.get("pressure", "N/A")
                mindistance = atoms.info.get("mindistance", "N/A")  # 新增

                # 默认值
                min_distance, atom_pair = "N/A", "N/A"
                if isinstance(mindistance, str):
                    parts = mindistance.split()
                    if len(parts) == 2:
                        min_distance = parts[0]
                        atom_pair = parts[1]

                # 写入额外信息
                extra_file.write(f"{temperature} {pressure} {volume:.6f} {min_distance} {atom_pair} {i}\n")

        print("数据已成功写入 E.dat, F.dat, virial.dat, stress.dat 和 extra_info.dat。")

    except Exception as e:
        print(f"发生错误: {e}")

def main():
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

