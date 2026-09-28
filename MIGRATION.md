# 迁移建议

先把旧脚本复制为 `legacy.py` 保留原始版本，然后逐个函数迁移。推荐顺序：

1. `WriteArc` → `arc_writer.write_arc`
2. `distanceMat/findneighbor` → `geometry`
3. `rotation* / xyz2rsf / supercell` → `geometry`
4. `addH/addNH3/addNH4/addCH3OH/getnormalizedmolecule` → `molecules`
5. `replace` → `structure.replace_atoms`
6. `addAdsorbedNH3` → `operations.add_adsorbate`
7. 最后再迁移 `getBAS` 等专项流程。

这样可以保证每次只改一块，并能用测试验证。
