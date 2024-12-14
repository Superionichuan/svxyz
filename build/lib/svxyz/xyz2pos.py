import sys
import os
import json
from ase.io import read, write

# 配置文件路径
CONFIG_FILE = "./xyz2pos.json"

def get_last_used_file():
    """
    从配置文件中获取上次使用的文件名。
    """
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, "r") as f:
            try:
                config = json.load(f)
                return config.get("last_used_file", "")
            except json.JSONDecodeError:
                return ""
    return ""

def save_last_used_file(filename):
    """
    将当前使用的文件名保存到配置文件。
    """
    config = {"last_used_file": filename}
    with open(CONFIG_FILE, "w") as f:
        json.dump(config, f)

def main():
    """
    Convert a specific frame in an XYZ file to POSCAR format.

    Usage:
        python xyz_to_poscar.py [frame_index]

    Args:
        frame_index (int): Index of the frame to extract (0-based).
    """
    # 获取上一次使用的文件名
    last_used_file = get_last_used_file()

    # 如果没有指定输入文件且没有记录文件，提示用户输入文件名
    if len(sys.argv) == 1 and not last_used_file:
        print("Error: No input file specified and no previous file found.")
        print("Usage: python xyz_to_poscar.py [input_xyz_file] <frame_index>")
        sys.exit(1)

    # 如果只提供了帧号，则使用记录的文件名
    if len(sys.argv) == 2:
        xyz_file = last_used_file
        frame_index = int(sys.argv[1])
    # 如果提供了文件名和帧号
    elif len(sys.argv) > 2:
        xyz_file = sys.argv[1]
        frame_index = int(sys.argv[2])
        save_last_used_file(xyz_file)  # 保存当前使用的文件名
    else:
        print("Usage: python xyz_to_poscar.py [frame_index]")
        sys.exit(1)

    # 自动生成输出文件名
    poscar_file = f"POSCAR{frame_index}"

    try:
        # 读取指定帧
        atoms = read(xyz_file, index=frame_index)

        # 写入到 POSCAR 格式
        write(poscar_file, atoms, format='vasp')
        print(f"Frame {frame_index} from {xyz_file} has been written to {poscar_file}.")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

