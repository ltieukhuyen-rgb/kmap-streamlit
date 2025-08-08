# app.py
import streamlit as st
import matplotlib.pyplot as plt
import re

# ---------------- Helpers ----------------
def int_to_bin_str(n, bits):
    return format(n, '0{}b'.format(bits))

def gray_code(n):
    return n ^ (n >> 1)

def mask_to_expr(mask):
    """mask: string like '1-0' -> expr 'A C' with primes"""
    vars_ = [chr(ord('A') + i) for i in range(len(mask))]
    parts = []
    for ch, v in zip(mask, vars_):
        if ch == '1':
            parts.append(v)
        elif ch == '0':
            parts.append(v + "'")
    return ''.join(parts) if parts else '1'

def mask_literal_count(mask):
    return sum(1 for ch in mask if ch != '-')

# ---------------- Quine-McCluskey + Petrick ----------------
def generate_prime_implicants(minterms, dont_cares, num_vars):
    """
    Return list of dicts: {'mask': mask_str, 'covers': set_of_ints}
    """
    # initial implicants from minterms + don't cares
    all_terms = sorted(set(minterms) | set(dont_cares))
    implicants = []
    for t in all_terms:
        mask = int_to_bin_str(t, num_vars)
        implicants.append({'mask': mask, 'covers': {t}})

    prime_implicants = []
    # iterative combining
    while True:
        # group by count of ones (ignoring '-')
        groups = {}
        for imp in implicants:
            ones = imp['mask'].count('1')
            groups.setdefault(ones, []).append(imp)

        new_map = {}  # mask_str -> covers set
        combined_flags = set()  # masks that got combined this round

        group_keys = sorted(groups.keys())
        for i in range(len(group_keys) - 1):
            g1 = groups[group_keys[i]]
            g2 = groups[group_keys[i + 1]]
            for a in g1:
                for b in g2:
                    m1 = a['mask']; m2 = b['mask']
                    # can combine if differ in exactly one position and positions with '-' must match
                    diffs = 0
                    ok = True
                    for x, y in zip(m1, m2):
                        if x != y:
                            if x == '-' or y == '-':
                                ok = False
                                break
                            diffs += 1
                            if diffs > 1:
                                ok = False
                                break
                    if not ok or diffs != 1:
                        continue
                    # combine
                    new_mask = ''.join(x if x == y else '-' for x, y in zip(m1, m2))
                    new_covers = a['covers'] | b['covers']
                    if new_mask in new_map:
                        new_map[new_mask] |= new_covers
                    else:
                        new_map[new_mask] = set(new_covers)
                    combined_flags.add(m1)
                    combined_flags.add(m2)

        # add uncombined implicants to prime_implicants
        for imp in implicants:
            if imp['mask'] not in combined_flags:
                # avoid duplicates by mask
                exists = next((p for p in prime_implicants if p['mask'] == imp['mask']), None)
                if exists is None:
                    prime_implicants.append({'mask': imp['mask'], 'covers': set(imp['covers'])})
                else:
                    exists['covers'] |= imp['covers']

        # prepare next round
        if not new_map:
            break
        implicants = [{'mask': m, 'covers': c} for m, c in new_map.items()]

    # deduplicate prime_implicants masks and unify covers
    uniq = {}
    for p in prime_implicants:
        if p['mask'] in uniq:
            uniq[p['mask']] |= p['covers']
        else:
            uniq[p['mask']] = set(p['covers'])
    result = [{'mask': m, 'covers': uniq[m], 'literals': mask_literal_count(m)} for m in sorted(uniq.keys())]
    return result

def petrick_method(chart, prime_info):
    """
    chart: dict minterm -> list of prime indices that cover it
    prime_info: list of prime implicant dicts (for weights)
    Returns: set of chosen prime indices (optimal by cardinality then literal count)
    """
    # Start product as list containing empty set
    product = [frozenset()]
    # for each minterm, multiply by sum(implicant indices covering that minterm)
    for m, covers in chart.items():
        if not covers:
            # impossible cover
            continue
        new_product = set()
        for term in product:
            for idx in covers:
                new_term = term | {idx}
                new_product.add(frozenset(new_term))
        # reduce by removing supersets
        # convert to list and sort by size then remove supersets
        terms = sorted(new_product, key=lambda s: (len(s), sum(prime_info[i]['literals'] for i in s)))
        reduced = []
        for t in terms:
            if not any(t > r for r in reduced):  # if t not a strict superset of any already kept
                reduced.append(t)
        product = reduced

    # product is list of frozensets (possible covers). Choose minimal by size then by literals weight.
    if not product:
        return set()
    min_size = min(len(s) for s in product)
    candidates = [s for s in product if len(s) == min_size]
    # tie-breaker: minimal literal count
    best = min(candidates, key=lambda s: sum(prime_info[i]['literals'] for i in s))
    return set(best)

def solve_minimization(minterms, dont_cares, num_vars):
    # 1) generate prime implicants
    prime_info = generate_prime_implicants(minterms, dont_cares, num_vars)
    if not prime_info:
        return {'prime_info': [], 'selected_idx': set(), 'final_exprs': [], 'essential_idx': set(), 'chart': {}}

    # 2) build prime implicant chart only for real minterms (not don't cares)
    chart = {}
    for m in minterms:
        covers = []
        for idx, p in enumerate(prime_info):
            if m in p['covers']:
                covers.append(idx)
        chart[m] = covers

    # 3) find essential implicants
    essential_idx = set()
    covered_minterms = set()
    for m, covers in chart.items():
        if len(covers) == 1:
            essential_idx.add(covers[0])
    # mark covered minterms by essentials
    for idx in list(essential_idx):
        covered_minterms |= set(prime_info[idx]['covers']) & set(minterms)

    # 4) remaining minterms
    remaining = [m for m in minterms if m not in covered_minterms]

    selected_idx = set(essential_idx)
    if remaining:
        # build reduced chart for remaining minterms
        reduced_chart = {m: chart[m] for m in remaining}
        # apply petrick
        petrick_choice = petrick_method(reduced_chart, prime_info)
        selected_idx |= set(petrick_choice)
    # final expressions
    final_exprs = []
    for idx in sorted(selected_idx):
        final_exprs.append({'mask': prime_info[idx]['mask'], 'expr': mask_to_expr(prime_info[idx]['mask']), 'covers': sorted(prime_info[idx]['covers']), 'idx': idx})
    return {
        'prime_info': prime_info,
        'selected_idx': selected_idx,
        'final_exprs': final_exprs,
        'essential_idx': essential_idx,
        'chart': chart
    }

# ---------------- K-map drawing ----------------
def draw_kmap(num_vars, minterms, dont_cares, selected_groups, colors):
    # xác định số hàng, cột
    rows = 1 << (num_vars // 2)
    cols = 1 << (num_vars - num_vars // 2)

    # danh sách Gray code cho hàng và cột
    row_gray_list = [gray_code(i) for i in range(rows)]
    col_gray_list = [gray_code(i) for i in range(cols)]

    fig, ax = plt.subplots(figsize=(cols * 0.9, rows * 0.9))
    ax.set_xticks(range(cols + 1))
    ax.set_yticks(range(rows + 1))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.grid(True, linewidth=0.6)

    # map từ tọa độ -> minterm đúng thứ tự
    cell_map = {}
    for r in range(rows):
        for c in range(cols):
            row_bits = row_gray_list[r]
            col_bits = col_gray_list[c]
            # hàng = MSB, cột = LSB
            val = (row_bits << (num_vars - num_vars // 2)) | col_bits
            cell_map[(r, c)] = val
            if val in minterms:
                ax.text(c + 0.5, r + 0.5, '1', va='center', ha='center', fontsize=12)
            elif val in dont_cares:
                ax.text(c + 0.5, r + 0.5, 'X', va='center', ha='center', fontsize=12, color='blue')

    # Tô màu nhóm đã chọn
    for g_idx, grp in enumerate(selected_groups):
        color = colors[g_idx % len(colors)]
        covered = set(grp['covers'])
        for cell, val in cell_map.items():
            if val in covered:
                ax.add_patch(plt.Rectangle((cell[1], cell[0]), 1, 1, color=color, alpha=0.35, linewidth=0))

    ax.invert_yaxis()
    ax.set_aspect('equal')
    plt.tight_layout()
    return fig


# ---------------- Streamlit UI ----------------
st.set_page_config(layout="wide")
st.title("K-map Minimizer — QM + Petrick (hiển thị chi tiết nhóm)")

# Sidebar inputs
num_vars = st.sidebar.slider("Số biến", 2, 6, 4)
minterms_input = st.sidebar.text_input("Minterms (cách nhau bằng dấu ,)", "1,3,7,11,15")
dont_cares_input = st.sidebar.text_input("Don't care terms", "2,5,6,9,10")

# parse inputs
def parse_list(s):
    items = []
    for v in s.split(","):
        v = v.strip()
        if v == "":
            continue
        try:
            items.append(int(v))
        except:
            st.error(f"Không thể parse số: {v}")
            st.stop()
    return items

minterms = parse_list(minterms_input)
dont_cares = parse_list(dont_cares_input)

# validate range and overlap
max_val = (1 << num_vars) - 1
for v in minterms + dont_cares:
    if v < 0 or v > max_val:
        st.error(f"Các giá trị phải nằm trong [0, {max_val}] cho {num_vars} biến. Giá trị {v} không hợp lệ.")
        st.stop()

# remove any dont_cares that are actually minterms
dont_cares = [d for d in dont_cares if d not in minterms]

if st.button("Tối giản và Vẽ K-map"):
    st.subheader("Dữ liệu đầu vào")
    c1, c2, c3 = st.columns([1, 2, 2])
    c1.write("**Số biến**"); c1.write(num_vars)
    c2.write("**Minterms**"); c2.write(minterms if minterms else "—")
    c3.write("**Don't cares**"); c3.write(dont_cares if dont_cares else "—")

    # run solver
    res = solve_minimization(minterms, dont_cares, num_vars)
    prime_info = res['prime_info']
    selected_idx = res['selected_idx']
    final_exprs = res['final_exprs']
    essential_idx = res['essential_idx']
    chart = res['chart']

    colors = ['#ffcccc', '#ccffcc', '#ccccff', '#ffffcc', '#ccffff', '#ffccff',
              '#ffd9b3', '#e6ccff', '#d9ffcc', '#ffb3b3']

    st.subheader("Prime implicants (tất cả)")
    if not prime_info:
        st.write("Không có prime implicant (có thể input rỗng).")
    else:
        # table-like display
        for idx, p in enumerate(prime_info):
            marker = ""
            if idx in essential_idx:
                marker = " (essential)"
            st.markdown(
                f"<div style='padding:6px;border-radius:6px;background:#f7f7f7;margin-bottom:4px;'>"
                f"<b>Index {idx}</b>{marker} &nbsp;|&nbsp; Mask: <code>{p['mask']}</code> &nbsp;|&nbsp; Covers: <code>{sorted(p['covers'])}</code> &nbsp;|&nbsp; Literals: {p['literals']}"
                f"</div>",
                unsafe_allow_html=True
            )

    # show prime implicant chart
    st.subheader("Prime implicant chart (minterm -> prime indices)")
    if chart:
        for m, covers in sorted(chart.items()):
            st.write(f"{m} -> {covers}")
    else:
        st.write("—")

    # show selected implicants
    st.subheader("Implicants được chọn (cuối cùng)")
    if not final_exprs:
        st.info("Không có implicant được chọn (hàm có thể là 0).")
    else:
        # display in color blocks in order of selected_idx (consistent ordering)
        sel_list = sorted(list(selected_idx))
        selected_groups = []
        for order, idx in enumerate(sel_list):
            p = prime_info[idx]
            block = {'mask': p['mask'], 'covers': sorted(p['covers']), 'expr': mask_to_expr(p['mask']), 'idx': idx}
            selected_groups.append(block)
            color = colors[order % len(colors)]
            essential_mark = " (essential)" if idx in essential_idx else ""
            st.markdown(
                f"<div style='background-color:{color};padding:8px;border-radius:6px;margin-bottom:6px;'>"
                f"<b>Selected {order+1} (prime idx {idx}){essential_mark}:</b> &nbsp; Mask: <code>{p['mask']}</code> &nbsp;|&nbsp; Covers: <code>{sorted(p['covers'])}</code> &nbsp;|&nbsp; Expr: <b>{mask_to_expr(p['mask'])}</b>"
                f"</div>",
                unsafe_allow_html=True
            )

        # final SOP (LaTeX)
        st.subheader("Biểu thức tối giản (SOP) - kết quả cuối")
        exprs = [block['expr'] for block in selected_groups]
        # convert A' to A^{\prime} for LaTeX
        latex_exprs = [re.sub(r"([A-Z])'", r"\1^{\\prime}", e) for e in exprs]
        if latex_exprs:
            st.latex(" + ".join(latex_exprs))
        else:
            st.write("0 (hàm luôn bằng 0)")

        # draw K-map using selected groups
        st.subheader("K-map (ô có '1' là minterms; màu = implicant đã chọn)")
        fig = draw_kmap(num_vars, minterms, selected_groups, colors)
        st.pyplot(fig)

