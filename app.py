import streamlit as st
import matplotlib.pyplot as plt
import re

# ====== Hàm hỗ trợ ======
def gray_code(n):
    return n ^ (n >> 1)

def num_to_binary_list(num, bits):
    """Trả về list bit (int) từ MSB -> LSB"""
    return [(num >> i) & 1 for i in reversed(range(bits))]

def term_to_mask(group, num_vars):
    """
    group: tuple hoặc list chứa minterm numbers (có thể nested)
    Trả về list các ký tự '0','1','-'
    """
    # Flatten group về list số nguyên
    flat_nums = []
    for item in group:
        if isinstance(item, (list, tuple, set)):
            flat_nums.extend(list(item))
        else:
            flat_nums.append(item)
    nums = sorted(set(int(x) for x in flat_nums))

    # bắt đầu từ số đầu tiên, chuyển thành chuỗi '0'/'1'
    bits = [str(b) for b in num_to_binary_list(nums[0], num_vars)]
    for other in nums[1:]:
        bits2 = [str(b) for b in num_to_binary_list(other, num_vars)]
        bits = ['-' if b1 != b2 else b1 for b1, b2 in zip(bits, bits2)]
    return bits


def term_to_expression(bits):
    """
    bits: list của '0','1','-'
    Trả về chuỗi biểu thức dạng AB'C...
    """
    variables = [chr(ord('A') + i) for i in range(len(bits))]
    expr_parts = []
    for bit, var in zip(bits, variables):
        if bit == '1':
            expr_parts.append(var)
        elif bit == '0':
            expr_parts.append(var + "'")
        # '-' bỏ qua
    return ''.join(expr_parts) if expr_parts else '1'

def expr_list_to_latex(expr_list):
    """
    Chuyển list các biểu thức (ví dụ ["A B'","C'D"]) sang LaTeX:
    - Thay X' -> X^{\\prime}
    - Nối bằng dấu +
    """
    latex_terms = [re.sub(r"([A-Z])'", r"\1^{\\prime}", term) for term in expr_list]
    return " + ".join(latex_terms)

# ====== Quine–McCluskey tối giản (đơn giản, trả về essential implicants) ======
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
            g1 = groups[all_group_keys[i]]
            g2 = groups[all_group_keys[i + 1]]
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

    # tìm essential implicants (rất cơ bản)
    essential = []
    for m in minterms:
        covers = [pi for pi in prime_implicants if m in pi]
        if len(covers) == 1:
            essential.append(covers[0])

    # Nếu không có essential nào (hoặc thiếu che phủ), trả về prime_implicants để người dùng thấy
    if not essential and prime_implicants:
        return list(sorted(prime_implicants))
    return list(sorted(essential))

# ====== Vẽ K-map ======
def draw_kmap(num_vars, minterms, groups, colors):
    """
    groups: list of iterables (groups), enumerate sẽ bắt đầu từ 0 -> color index tương ứng
    """
    rows = 1 << (num_vars // 2)
    cols = 1 << (num_vars - num_vars // 2)

    fig, ax = plt.subplots(figsize=(cols * 0.6, rows * 0.6))
    ax.set_xticks(range(cols + 1))
    ax.set_yticks(range(rows + 1))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.grid(True, linewidth=0.6)

    cell_map = {}
    for r in range(rows):
        for c in range(cols):
            row_gray = gray_code(r)
            col_gray = gray_code(c)
            val = (row_gray << (num_vars - num_vars // 2)) | col_gray
            cell_map[(r, c)] = val
            if val in minterms:
                ax.text(c + 0.5, r + 0.5, '1', va='center', ha='center', fontsize=12)

    # Tô màu nhóm (mỗi ô nếu thuộc nhóm sẽ được tô)
    for g_idx, group in enumerate(groups):
        color = colors[g_idx % len(colors)]
        for cell, val in cell_map.items():
            if val in group:
                ax.add_patch(plt.Rectangle((cell[1], cell[0]), 1, 1, color=color, alpha=0.35, linewidth=0))

    ax.invert_yaxis()
    ax.set_aspect('equal')
    plt.tight_layout()
    return fig

# ============================ Streamlit UI ============================
st.set_page_config(layout="wide")
st.title("K-map Minimizer 2–6 biến (sửa lỗi & hiển thị chi tiết nhóm)")

# Input
num

