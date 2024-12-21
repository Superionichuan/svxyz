import os
import sys
import json
from ase.io import read, write

# 配置文件路径
CONFIG_FILE = "./xyz2pos.json"

# 默认支持的格式
SUPPORTED_FORMATS = {
    "input_formats": [
        "xyz", "vasp-xml", "vasp-out", "lammps-data", "pdb", "cube",
        "cif", "xsf", "json", "traj", "dcd",
        "proteindatabank", "ase", "extxyz", "mdl", "povray",
        "espresso-in", "espresso-out", "qe-in", "qe-out", "gaussian"
    ],
    "output_formats": [
        "vasp", "xyz", "lammps-data", "cif", "pdb",
        "xsf", "traj", "json", "extxyz", "proteindatabank",
        "povray", "espresso-in", "espresso-out", "qe-in", "qe-out",
        "cube", "dcd", "mdl", "gaussian", "aims"
    ]
}

DEFAULT_CONFIG = {
    "last_used_file": None,
    "input_format": None,  # 默认不指定格式
    "output_format": "vasp",  # 默认输出格式为 VASP 格式
    "frame_index": 0  # 默认帧索引
}

def create_default_config():
    """生成默认的 JSON 配置文件"""
    with open(CONFIG_FILE, "w") as f:
        json.dump(DEFAULT_CONFIG, f, indent=4)
    print(f"Default configuration created: {CONFIG_FILE}")

def load_or_create_config():
    """加载配置文件，如果不存在则生成并退出"""
    if not os.path.exists(CONFIG_FILE):
        create_default_config()
        print("JSON configuration file not found. Generated a default configuration. Exiting now.")
        sys.exit(0)
    with open(CONFIG_FILE, "r") as f:
        return json.load(f)

def update_config(config):
    """更新配置文件"""
    with open(CONFIG_FILE, "w") as f:
        json.dump(config, f, indent=4)

def show_help():
    """显示帮助信息"""
    print("""
Usage: xyz2pos [OPTIONS] [input_file] <frame_index>

Convert a specific frame from an input file to a specified format.

Options:
  --ifmt FORMAT        Input file format. Supported formats: {}
  --ofmt FORMAT        Output file format. Supported formats: {}
  -h, --help           Show this help message and exit.

Arguments:
  [input_file]         Path to the input file. If not specified, the last used file from the JSON configuration will be used.
  <frame_index>        Index of the frame to extract (0-based). This is required.
    """.format(
        ", ".join(SUPPORTED_FORMATS["input_formats"]), ", ".join(SUPPORTED_FORMATS["output_formats"])
    ))

def parse_arguments(args, config):
    """解析命令行参数并更新配置"""
    if "-h" in args or "--help" in args:
        show_help()
        sys.exit(0)

    if "--ifmt" in args:
        input_format = args[args.index("--ifmt") + 1]
        if input_format not in SUPPORTED_FORMATS["input_formats"]:
            print(f"Error: Unsupported input format '{input_format}'.")
            sys.exit(1)
        config["input_format"] = input_format

    if "--ofmt" in args:
        output_format = args[args.index("--ofmt") + 1]
        if output_format not in SUPPORTED_FORMATS["output_formats"]:
            print(f"Error: Unsupported output format '{output_format}'.")
            sys.exit(1)
        config["output_format"] = output_format

    if len(args) >= 2:
        config["last_used_file"] = args[-2]
        try:
            config["frame_index"] = int(args[-1])
        except ValueError:
            print(f"Error: Invalid frame index '{args[-1]}'. It must be an integer.")
            sys.exit(1)

    update_config(config)

def main():
    """主函数"""
    config = load_or_create_config()
    args = sys.argv[1:]

    if args:
        parse_arguments(args, config)

    input_file = config.get("last_used_file")
    input_format = config.get("input_format")
    output_format = config.get("output_format")
    frame_index = config.get("frame_index", 0)

    if not input_file:
        print("Error: No input file specified.")
        sys.exit(1)

    # 自动生成输出文件名
    output_file = f"POSCAR{frame_index}" if output_format == "vasp" else f"output_frame{frame_index}.{output_format}"

    print(f"Reading file: {input_file}")
    try:
        if input_format:
            atoms = read(input_file, index=frame_index, format=input_format)
        else:
            atoms = read(input_file, index=frame_index)

        write(output_file, atoms, format=output_format)
        print(f"Frame {frame_index} from {input_file} has been written to {output_file}.")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

