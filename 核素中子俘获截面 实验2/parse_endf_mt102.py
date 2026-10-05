# -*- coding: utf-8 -*-
"""
parse_endf_mt102.py — 从 ENDF-6 文件提取 MF=3 / MT=102 (n,γ) 平滑截面。
用法: python3 parse_endf_mt102.py <endf文件> <输出txt>
ENDF-6 记录: 每行 80 字符 = 66 位数据 + MAT(4) MF(2) MT(3) NS(4)。
MT=102 数据以 LIST 格式 (E,σ) 成对出现。
"""
import sys, re
import numpy as np

def parse_mt102(fname):
    with open(fname, 'r', errors='replace') as f:
        lines = f.read().split('\n')
    rec = re.compile(r'^(\d{4})\s*(\d)\s*(\d{3})\s*(\d+)\s*$')
    endf_num = re.compile(r'^([-+]?\d*\.?\d+)([+-]\d+)$')
    def endf_float(tok):
        """ENDF Fortran 科学计数(如 8.000000+3) 转 python float。"""
        tok = tok.strip()
        m = endf_num.match(tok)
        return float(m.group(1) + 'e' + m.group(2)) if m else float(tok)
    E, S = [], []
    for ln in lines:
        if len(ln) < 70:
            continue
        body = ln[:66]
        tail = ln[66:].strip()
        m = rec.match(tail)
        if not m:
            continue
        mat, mf, mt, ns = m.groups()
        if int(mf) != 3 or int(mt) != 102:
            continue
        if int(ns) < 5:                # 跳过 TAB1 头部(ZA/INT/断点行)
            continue
        # 数据行：66 位内为 6 个 11 位浮点
        try:
            vals = [endf_float(body[11*i:11*(i+1)]) for i in range(6)]
        except ValueError:
            continue
        for j in range(0, 6, 2):
            e, s = vals[j], vals[j+1]
            if e > 0 and s > 0:          # 过滤 TAB1 断点行(σ=0)
                E.append(e); S.append(s)
    if not E:
        raise ValueError(f"{fname}: MF=3/MT=102 数据解析失败")
    return np.array(E), np.array(S)

if __name__ == "__main__":
    E, S = parse_mt102(sys.argv[1])
    np.savetxt(sys.argv[2], np.column_stack([E, S]), header="E[eV]  sigma[barn]", fmt="%.8g")
    print(f"{sys.argv[1]}: {len(E)} 个点, E∈[{E.min():.2e},{E.max():.2e}] eV, σ∈[{S.min():.3g},{S.max():.3g}] b")
