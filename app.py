import streamlit as st
import matplotlib.pyplot as plt
import itertools
import random

# Hàm chuyển số sang mã Gray
def gray_code(n):
    return n ^ (n >> 1)

# Thuật toán Quine–McCluskey đơn giản
def qm_minimize(minterms, dont_cares, num_vars):
    terms = sorted(set(minterms) | set(dont_cares))
    groups = {}
    for t in terms:
        ones = bin(t).count('1')
        groups.setdefault(ones, set()).add((t,))

    prime_implicants = set()

    while groups:
        new_groups = {}
        marked = set()
        all_group_keys = sorted(groups.keys())
        for i in range(len(all_group_keys) - 1):
            g1, g2 = groups[all_group_keys[i]], groups[all_group_keys[i+1]]
            for term1 in g1:
                for term2 in g2:
                    diff = term1[0] ^ term2[0]
                    if bin(diff).count('1') == 1:
                        merged = tuple(sorted(set(term1) | set(term2)))
                        new_groups.setdefault(all_group_keys[i], set()).add(merged)
                        marked.add(term1)
                        marked.add(term2)
        for group in groups.values():
            for term in group:
                if term not in marked:
                    prime_implicants.add(term)
        groups = new_groups

    essential = set()
    for m in minterms:
        covers = [pi for pi in prime_implicants if m in pi]
        if len(covers) == 1:
            essential.add(covers[0])
    return list(essential)

# Vẽ K-map với nhóm tô màu nền
def draw_kmap(num_vars, minterms, groups):
    if num_vars <= 3:
        rows = 1 << (num_vars // 2)
        cols = 1 << (num_vars - num_vars // 2)
    else:
        rows = 1 << (num_vars // 2)
        cols = 1 << (num_vars - num_vars // 2)

    fig, ax = plt.subplots()
    ax.set_xticks(range(cols+1))
    ax.set_yticks(range(rows+1))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.grid(True)

    # Map giá trị vào ô K-map
    cell_map = {}
    for r in range(rows):
        for c in range(cols):
            row_gray = gray_code(r)
            col_gray = gray_code(c)
            val = (row_gray << (num_vars - num_vars // 2)) | col_gray
            cell_map[(r, c)] = val
            if val in minterms:
                ax.text(c+0.5, r+0.5, '1', va='center', ha='center', fontsize=14)

    # Tô màu nhóm
    colors = ['#ffcccc','#ccffcc','#ccccff','#ffffcc','#ccffff','#ffccff']
    for g_idx, group in enumerate(groups):
        color = colors[g_idx % len(colors)]
        for cell, val in cell_map.items():
            if val in group:
                ax.add_patch(plt.Rectangle((cell[1], cell[0]), 1, 1, color=color, alpha=0.3))

    ax.invert_yaxis()
    return fig

# ============================
st.title("K-map Minimizer 2–6 biến")

num_vars = st.slider("Số biến", 2, 6, 3)
minterms_input = st.text_input("Minterms (cách nhau bằng dấu ,)", "1,3,7")
dont_cares_input = st.text_input("Don't care terms", "")

try:
    minterms = [int(x) for x in minterms_input.split(",") if x.strip() != ""]
    dont_cares = [int(x) for x in dont_cares_input.split(",") if x.strip() != ""]
except:
    st.error("Nhập sai định dạng!")
    st.stop()

if st.button("Tối giản và Vẽ K-map"):
    groups = qm_minimize(minterms, dont_cares, num_vars)
    st.write("Nhóm tối giản:", groups)
    fig = draw_kmap(num_vars, minterms, groups)
    st.pyplot(fig)
