import matplotlib.pyplot as plt
import numpy as np
import json
import os
import sys
from scipy.stats import gaussian_kde
import mplcursors  # 用于动态提示

# 单位映射和列标识
UNITS = {
    "E": "Energy (eV)",
    "F": "Force Norm (eV/Å)",
    "virial": "Virial (eV)",
    "stress": "Stress Components (GPa)"
}

LABELS = {
    "E": ["Energy"],
    "F": ["||F|| max","||F|| mean"],
    "virial": ["$\\tau_{xx}$", "$\\tau_{yy}$", "$\\tau_{zz}$", "$\\tau_{yz}$", "$\\tau_{xz}$", "$\\tau_{xy}$"],
    "stress": ["$\\sigma_{xx}$", "$\\sigma_{yy}$", "$\\sigma_{zz}$", "$\\sigma_{yz}$", "$\\sigma_{xz}$", "$\\sigma_{xy}$"]
}

# 默认配置
DEFAULT_CONFIG = {
    "ylabel": "$data$",
    "show": True,
    "output_file": "output.pdf",
    "x1min": None,
    "x1max": None,
    "y1min": None,
    "y1max": None,
    "x2min": None,
    "x2max": None,
    "y2min": None,
    "y2max": None,
    "cmap": "viridis"
}

# 应力转换因子
STRESS_CONVERSION_FACTOR = 160.21766  # eV/Å³ to GPa

def load_data(file, indices, data_type):
    """
    加载数据，根据文件类型进行处理。

    Args:
        file (str): 数据文件路径。
        indices (list): 选择的曲线索引。
        data_type (str): 数据类型。

    Returns:
        tuple: 数据矩阵 (list of curves), System_ID 列表, 和曲线标题（可选）。
    """
    if data_type == "infile":
        # 读取 infile.dat 文件
        with open(file, "r") as f:
            data = [line.split() for line in f.readlines()]
        
        # 提取标题和数值部分
        #titles = [row[0] for row in data]  # 第一列为标题
        titles = [row[0].replace("\\\\", "\\") for row in data] 
        values = [list(map(float, row[1:])) for row in data]  # 其余列为数值数据
        
        # 检查索引范围
        if max(indices) >= len(values):
            raise IndexError(f"Index {max(indices)} is out of bounds for data with {len(values)} curves.")

        # 根据 indices 选择曲线数据
        curves = [values[i] for i in indices]
        selected_titles = [titles[i] for i in indices]  # 同时选择对应的标题
        
        # 对齐曲线长度（填充 NaN 以对齐）
        max_length = max(len(curve) for curve in curves)
        aligned_curves = [curve + [np.nan] * (max_length - len(curve)) for curve in curves]
        
        # 定义 System_ID 列表
        system_ids = list(range(max_length))
        
        # 返回对齐的曲线、System_ID 和选择的标题
        return aligned_curves, system_ids, selected_titles

    else:
        # 处理其他数据类型（如 E.dat 等）
        data = np.loadtxt(file, skiprows=1)
        curves = [data[:, i] for i in indices]
        system_ids = data[:, -1].astype(int)
        titles = None
        return curves, system_ids, titles




def plot_distribution_and_projection(curves, system_ids, ylabel, labels, titles=None, config=None,  indices=None):
    """
    绘制分布图和投影图。

    Args:
        curves (list): 数据曲线列表，每个子列表表示一条曲线。
        system_ids (list): System_ID 列表。
        ylabel (str): Y轴标签。
        labels (list): 每列的标识信息。
        titles (list, optional): 曲线标题，仅在 infile 模式下使用。
        config (dict): 配置字典。
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    # 左图：分布曲线
    ax1 = axes[0]
    for i, curve in enumerate(curves):
        curve_array = np.array(curve)
        curve_array = curve_array[~np.isnan(curve_array)]

        kde = gaussian_kde(curve_array)
        x = np.linspace(curve_array.min(), curve_array.max(), 1000)

        mean = np.mean(curve_array)
        peak_x = x[np.argmax(kde(x))]
        std_dev = np.std(curve_array)

        color = f"C{i}" if len(curves) > 1 else "blue"
        #legend_label = titles[i] if titles else labels[i]
        #legend_label = titles[i] if titles else labels[indices[i]]
        # 修正 legend_label 的生成逻辑
        if titles:  # 如果是 infile 模式
            legend_label = titles[i]
        elif indices:  # 如果是其他模式且 indices 存在
            legend_label = labels[indices[i]]
        else:  # 默认情况
            legend_label = labels[i]



        line, = ax1.plot(x, kde(x), label=legend_label, color=color, linewidth=1.5)

        cursor = mplcursors.cursor(line, hover=True)
        cursor.connect(
            "add",
            lambda sel, lbl=legend_label, mean=mean, peak_x=peak_x, std_dev=std_dev: (
                sel.annotation.set(
                    text=f"{lbl}\n$\\mu={mean:.2f}$\n$x_m={peak_x:.2f}$\n$\\sigma={std_dev:.2f}$",
                    color="black"
                ),
                sel.annotation.get_bbox_patch().set_edgecolor("none"),
                sel.annotation.get_bbox_patch().set_facecolor("none"),
                sel.annotation.get_bbox_patch().set_alpha(0)
            )
        )

        ax1.scatter(
            curve_array, kde(curve_array),
            c=curve_array, cmap=config.get("cmap", "viridis"),
            edgecolor="black" if len(curves) == 1 else color,
            alpha=1.0,
            label=f"{legend_label} Points"
        )

    if config["x1min"] or config["x1max"]:
        ax1.set_xlim(config["x1min"], config["x1max"])
    if config["y1min"] or config["y1max"]:
        ax1.set_ylim(config["y1min"], config["y1max"])

    ax1.set_xlabel(ylabel)
    ax1.set_ylabel("Density")
    ax1.set_title(f"{ylabel} Distribution")
    ax1.legend()

    # 右图：投影图
    ax2 = axes[1]
    for i, curve in enumerate(curves):
        curve_array = np.array(curve)
        curve_array = curve_array[~np.isnan(curve_array)]

        scatter = ax2.scatter(
            system_ids[:len(curve_array)], curve_array,
            c=curve_array, cmap=config.get("cmap", "viridis"),
            edgecolor="black" if len(curves) == 1 else f"C{i}",
            alpha=1.0
        )

        # 为右图的点添加动态标签
        cursor = mplcursors.cursor(scatter, hover=True)
        cursor.connect(
            "add",
            lambda sel, lbl=(titles[i] if titles else labels[i]): (
                sel.annotation.set(
                    text=f"{lbl}\nSystem ID: {system_ids[int(sel.index)]}\nValue: {sel.target[1]:.2f}",
                    color="black"
                ),
                sel.annotation.get_bbox_patch().set_edgecolor("none"),  # 边框无色
                sel.annotation.get_bbox_patch().set_facecolor("none"),  # 背景无色
                sel.annotation.get_bbox_patch().set_alpha(0)           # 背景透明
            )
        )




    if config["x2min"] or config["x2max"]:
        ax2.set_xlim(config["x2min"], config["x2max"])
    if config["y2min"] or config["y2max"]:
        ax2.set_ylim(config["y2min"], config["y2max"])

    ax2.set_xlabel("System ID")
    ax2.set_ylabel(ylabel)
    ax2.set_title(f"{ylabel} Projection")
    if len(curves) == 1:
        plt.colorbar(scatter, ax=ax2, label=ylabel)

    plt.tight_layout()
    if config["show"]:
        plt.show()
    else:
        plt.savefig(config["output_file"])
        print(f"Saved plot to {config['output_file']}")


def main():
    json_file = "pxyz.json"
    if not os.path.exists(json_file):
        with open(json_file, "w") as f:
            json.dump(DEFAULT_CONFIG, f, indent=4)
        print(f"No configuration file found. A default configuration has been created as {json_file}.")
        sys.exit(0)

    with open(json_file, "r") as f:
        config = json.load(f)

    if len(sys.argv) < 3:
        print("Usage: python pxyz.py <data_type> <indices...>")
        sys.exit(1)

    data_type = sys.argv[1]
    indices = list(map(int, sys.argv[2:]))

    data_files = {
        "E": "E.dat",
        "F": "F.dat",
        "virial": "virial.dat",
        "stress": "stress.dat",
        "infile": "infile.dat"
    }

    if data_type not in data_files:
        print(f"Unsupported data type: {data_type}")
        sys.exit(1)

    file = data_files[data_type]
    ylabel = UNITS.get(data_type, config["ylabel"])
    labels = LABELS.get(data_type, ["Column"])

    curves, system_ids, titles = load_data(file, indices, data_type)
    plot_distribution_and_projection(curves, system_ids, ylabel, labels, titles, config,  indices)


if __name__ == "__main__":
    main()

