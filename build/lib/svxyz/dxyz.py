import sys
import os
import json
from ase.io import read
import numpy as np

# eV/Å³ 转 GPa 的转换因子
EV_A3_TO_GPA = 160.21766208

def extract_properties_with_ids(xyz_file):
    """
    从 XYZ 文件中提取能量、力、virial 和应力数据，分别写入 E.dat, F.dat, virial.dat, 和 stress.dat。
    应力值转换为 GPa 并取相反值。

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
             open("stress.dat", "w") as stress_file:

            # 写入标题行
            e_file.write("Energy(eV) System_ID\n")
            f_file.write("Max_atomic_force_norm Mean_atomic_force_norm System_ID\n")
            vir_file.write("$\\tau_{xx}$ $\\tau_{yy}$ $\\tau_{zz}$ $\\tau_{yz}$ $\\tau_{xz}$ $\\tau_{xy} System_ID\n")
            stress_file.write("$\\sigma_{xx}$ $\\sigma_{yy}$ $\\sigma_{zz}$ $\\sigma_{yz} $\\sigma_{xz} $\\sigma_{xy} System_ID\n")

            # 遍历每一帧
            for i, atoms in enumerate(atoms_list):
                # 提取能量
                energy = atoms.get_potential_energy()
                e_file.write(f"{energy:.6f} {i}\n")

                # 提取力范数
                forces = atoms.get_forces()
                mean_force_norm = np.linalg.norm(forces, axis=1).mean()  # 每个原子的力范数的平均值
                max_force_norm = np.linalg.norm(forces, axis=1).max()  # 每个原子的力范数的最大值
                f_file.write(f"{max_force_norm:.6f} {mean_force_norm:.6f} {i}\n")

                # 提取 virial 和应力六分量
                stress = atoms.get_stress(voigt=True)  # Voigt 6 分量表示
                volume = atoms.get_volume()  # ASE 中的体积
                virial = -stress * volume  # Virial 张量的六分量形式

                # 转换 stress 为 GPa 并取负值
                stress_gpa = -stress * EV_A3_TO_GPA

                # 写入 virial 六分量
                vir_file.write(f"{virial[0]:.6f} {virial[1]:.6f} {virial[2]:.6f} {virial[3]:.6f} {virial[4]:.6f} {virial[5]:.6f} {i}\n")

                # 写入应力六分量（单位 GPa，取负值）
                stress_file.write(f"{stress_gpa[0]:.6f} {stress_gpa[1]:.6f} {stress_gpa[2]:.6f} {stress_gpa[3]:.6f} {stress_gpa[4]:.6f} {stress_gpa[5]:.6f} {i}\n")

        print("数据已成功写入 E.dat, F.dat, virial.dat 和 stress.dat。")

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

