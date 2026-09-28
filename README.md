# zeolite-structure-tools

把原来的单文件分子筛结构处理脚本拆成可维护的 Python 项目。原脚本共 395 行，包含 ARC 输出、超胞、距离矩阵、分子旋转、NH3/NH4/CH3OH 构造、AlPO/吸附结构生成等逻辑。fileciteturn5file0L22-L37

## 主要改造
- 几何算法 → `geometry.py`
- 分子模板与构造 → `molecules.py`
- 结构替换/拼接 → `structure.py`
- ARC 输入输出 → `io.py` + `arc_writer.py`
- 具体科研流程 → `operations/`
- 命令行入口 → `cli/`
- 单元测试 → `tests/`
- 不再使用 `DataFrame.append()`；统一 `pd.concat()`
- `xyz2rsf()` 改为 `atan2`，避免原代码大量手工象限判断。原函数位于旧脚本 103--124 行。fileciteturn5file0L103-L124
- 不再把 `/home10/msc/bin/...` 写死在业务代码里；路径通过参数/配置传入。旧脚本中确实存在固定路径，例如周期表和分子模板路径。fileciteturn5file0L20-L20
- 随机数改成 `numpy.random.Generator`，支持 seed，便于科研结果复现。

## 安装
```bash
pip install -e .
```

## 运行
```bash
python -m zeolite_structure_tools.cli.main --input-dir acid_sites --output-dir output --molecule-name h2o
```
