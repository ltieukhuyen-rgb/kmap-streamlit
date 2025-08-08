# app.py
import streamlit as st
import matplotlib.pyplot as plt
import re

# ====== Hàm hỗ trợ ======
def gray_code(n):
    return n ^ (n >> 1)

def num_to_binary_list(num, bits):
    """Trả về list bit (int) từ MSB -> LSB"""
    return [(num >> i) & 1 for i in reversed(range(bits))]

def flatten_group(group):
    """Flatten một group có thể chứa int hoặc tuple/list/set lồng nhau -> trả về list số nguyên"""
    nums = []
    for item in group:
        if isinstance(item, (list, tuple, set)):
            for sub in item:
                if isinstance(sub, (list, tuple, set)):
                    nums.extend(flatten_group(sub))
                else:
                    nums.append(int(sub))
        else:
            nums.append(int(item))
    return nums

def term_to_mask(group, num_vars):
    """
    group: tuple hoặc iterable các minterm số (ints) hoặc nested.
    Trả về list các ký tự '0','1','-' (chuỗi) biểu diễn mask.
    """
    nums = sorted(set(flatten_group(group)))
    if not nums:
        return ['-'] * num_vars
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
    return ''.join(expr_parts) if expr_parts else '1'

def expr_list_to_latex(expr_list):
    """
    Chuyển danh sách biểu thức (A B'...) sang LaTeX-friendly: A^{\prime} thay cho A'
    """
    latex_terms = [re.sub(r"([A-Z])'", r"\1^{\\prime}", term) for term in expr_list]
    return " + ".join(latex_terms)

# ====== Quine–McCluskey tối giản (triển khai đơn giản) ======
def qm_minimize(minterms, dont_cares, num_vars):
    """
    Trả về danh sách prime/essential implicants ở dạng tuple các minterm.
    (Lưu ý: đây là phiên bản giản lược; không có Petrick method đầy đủ.)
    """
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

    # Tìm essential implicants (cách cơ bản: implicant nào là cover duy nhất cho một minterm)
    essential = []
    for m in minterms:
        covers = [pi for pi in prime_implicants if m in pi]
        if len(covers) == 1:
            if covers[0] not in essential:
                essential.append(covers[0])

    # Nếu không có essential, trả về prime implicants để người dùng thấy
    if not essential and prime_implicants:
        return sorted(prime_implicants)
    return sorted(essential)

# ====== Vẽ K-map ======
def draw_kmap(num_vars, minterms, groups, colors):
    rows = 1 << (num_vars // 2)
    cols = 1 << (num_vars - num_vars // 2)

    fig, ax = plt.subplots(figsize=(cols * 0.9, rows * 0.9))
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
        covered = set(flatten_group(group))
        for cell, val in cell_map.items():
            if val in covered:
                ax.add_patch(plt.Rectangle((cell[1], cell[0]), 1, 1, color=color, alpha=0.35, linewidth=0))

    ax.invert_yaxis()
    ax.set_aspect('equal')
    plt.tight_layout()
    return fig

# ============================ Streamlit UI ============================
st.set_page_config(layout="wide")
st.title("K-map Minimizer 2–6 biến — Hiển thị chi tiết nhóm")

# Sidebar input
num_vars = st.sidebar.slider("Số biến", 2, 6, 4)
minterms_input = st.sidebar.text_input("Minterms (cách nhau bằng dấu ,)", "1,3,7,11,15")
dont_cares_input = st.sidebar.text_input("Don't care terms", "2,5,6,9,10")

# Parse
try:
    minterms = [int(x) for x in minterms_input.split(",") if x.strip() != ""]
    dont_cares = [int(x) for x in dont_cares_input.split(",") if x.strip() != ""]
except ValueError:
    st.error("Nhập sai định dạng! Hãy chỉ nhập số nguyên, cách nhau bởi dấu phẩy.")
    st.stop()

if st.button("Tối giản và Vẽ K-map"):
    # Hiển thị dữ liệu đầu vào
    st.subheader("Dữ liệu đầu vào")
    c1, c2, c3 = st.columns([1, 2, 2])
    c1.write("**Số biến**")
    c1.write(num_vars)
    c2.write("**Minterms**")
    c2.write(minterms if minterms else "—")
    c3.write("**Don't cares**")
    c3.write(dont_cares if dont_cares else "—")

    # Tối giản
    groups = qm_minimize(minterms, dont_cares, num_vars)

    # Màu cho nhóm
    colors = ['#ffcccc', '#ccffcc', '#ccccff', '#ffffcc', '#ccffff', '#ffccff',
              '#ffd9b3', '#e6ccff', '#d9ffcc', '#ffb3b3']

    # Hiển thị nhóm chi tiết
    st.subheader("Các nhóm tối giản (chi tiết)")
    if not groups:
        st.info("Không tìm thấy nhóm tối giản (có thể input rỗng hoặc thuật toán không tìm được implicant).")
    else:
        expr_list = []
        for idx, g in enumerate(groups, 1):
            mask = term_to_mask(g, num_vars)
            mask_str = ''.join(mask)
            expr = term_to_expression(mask)
            expr_list.append(expr)
            covered = sorted(set(flatten_group(g)))
            color = colors[(idx - 1) % len(colors)]

            # Hiển thị 1 khối màu cho mỗi nhóm, bao gồm ô bao phủ, mask, biểu thức
            st.markdown(
                f"<div style='background-color:{color};padding:8px;border-radius:6px;margin-bottom:6px;'>"
                f"<b>Nhóm {idx}:</b> <code>{covered}</code> &nbsp;|&nbsp; "
                f"<b>Mask:</b> <code>{mask_str}</code> &nbsp;|&nbsp; "
                f"<b>Biểu thức:</b> <code>{expr}</code>"
                f"</div>",
                unsafe_allow_html=True
            )

        # Biểu thức SOP (LaTeX)
        st.subheader("Biểu thức tối giản (SOP)")
        latex_sop = expr_list_to_latex(expr_list) if expr_list else "0"
        st.latex(latex_sop)

    # Vẽ K-map bên cạnh chi tiết (dùng 2 cột)
    st.subheader("K-map (ô có '1' = minterms; màu = nhóm)")
    fig = draw_kmap(num_vars, minterms, groups, colors)
    st.pyplot(fig)
