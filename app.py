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
    all_terms = sorted(set(minterms) | set(dont_cares))
    implicants = [{'mask': int_to_bin_str(t, num_vars), 'covers': {t}} for t in all_terms]

    prime_implicants = []
    while True:
        groups = {}
        for imp in implicants:
            ones = imp['mask'].count('1')
            groups.setdefault(ones, []).append(imp)

        new_map = {}
        combined_flags = set()

        group_keys = sorted(groups.keys())
        for i in range(len(group_keys) - 1):
            for a in groups[group_keys[i]]:
                for b in groups[group_keys[i + 1]]:
                    m1, m2 = a['mask'], b['mask']
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
                    new_mask = ''.join(x if x == y else '-' for x, y in zip(m1, m2))
                    new_covers = a['covers'] | b['covers']
                    new_map.setdefault(new_mask, set()).update(new_covers)
                    combined_flags.add(m1)
                    combined_flags.add(m2)

        for imp in implicants:
            if imp['mask'] not in combined_flags:
                exists = next((p for p in prime_implicants if p['mask'] == imp['mask']), None)
                if exists is None:
                    prime_implicants.append({'mask': imp['mask'], 'covers': set(imp['covers'])})
                else:
                    exists['covers'].update(imp['covers'])

        if not new_map:
            break
        implicants = [{'mask': m, 'covers': c} for m, c in new_map.items()]

    uniq = {}
    for p in prime_implicants:
        uniq.setdefault(p['mask'], set()).update(p['covers'])
    return [{'mask': m, 'covers': uniq[m], 'literals': mask_literal_count(m)} for m in sorted(uniq.keys())]

def petrick_method(chart, prime_info):
    product = [frozenset()]
    for covers in chart.values():
        if not covers:
            continue
        new_product = set()
        for term in product:
            for idx in covers:
                new_term = term | {idx}
                new_product.add(frozenset(new_term))
        terms = sorted(new_product, key=lambda s: (len(s), sum(prime_info[i]['literals'] for i in s)))
        reduced = []
        for t in terms:
            if not any(t > r for r in reduced):
                reduced.append(t)
        product = reduced

    if not product:
        return set()
    min_size = min(len(s) for s in product)
    candidates = [s for s in product if len(s) == min_size]
    return min(candidates, key=lambda s: sum(prime_info[i]['literals'] for i in s))

def solve_minimization(minterms, dont_cares, num_vars):
    prime_info = generate_prime_implicants(minterms, dont_cares, num_vars)
    chart = {m: [idx for idx, p in enumerate(prime_info) if m in p['covers']] for m in minterms}

    essential_idx = {covers[0] for m, covers in chart.items() if len(covers) == 1}
    covered_minterms = set()
    for idx in essential_idx:
        covered_minterms |= (prime_info[idx]['covers'] & set(minterms))

    remaining = [m for m in minterms if m not in covered_minterms]
    selected_idx = set(essential_idx)
    if remaining:
        reduced_chart = {m: chart[m] for m in remaining}
        selected_idx |= petrick_method(reduced_chart, prime_info)

    final_exprs = [{'mask': prime_info[idx]['mask'],
                    'expr': mask_to_expr(prime_info[idx]['mask']),
                    'covers': sorted(prime_info[idx]['covers']),
                    'idx': idx}
                   for idx in sorted(selected_idx)]
    return prime_info, selected_idx, final_exprs, essential_idx

# ---------------- K-map drawing ----------------
def draw_kmap(num_vars, minterms, dont_cares, selected_groups, colors):
    rows = 1 << (num_vars // 2)
    cols = 1 << (num_vars - num_vars // 2)
    row_gray_list = [gray_code(i) for i in range(rows)]
    col_gray_list = [gray_code(i) for i in range(cols)]

    fig, ax = plt.subplots(figsize=(cols * 0.9, rows * 0.9))
    ax.set_xticks(range(cols + 1))
    ax.set_yticks(range(rows + 1))
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.grid(True, linewidth=0.6)

    cell_map = {}
    for r in range(rows):
        for c in range(cols):
            row_bits = row_gray_list[r]
            col_bits = col_gray_list[c]
            val = (row_bits << (num_vars - num_vars // 2)) | col_bits
            cell_map[(r, c)] = val
            if val in minterms:
                ax.text(c + 0.5, r + 0.5, '1', va='center', ha='center', fontsize=12)
            elif val in dont_cares:
                ax.text(c + 0.5, r + 0.5, 'X', va='center', ha='center', fontsize=12, color='blue')

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
st.title("K-map Minimizer — QM + Petrick (đúng mapping & hiển thị X)")

num_vars = st.sidebar.slider("Số biến", 2, 6, 4)
minterms_input = st.sidebar.text_input("Minterms (cách nhau bằng dấu ,)", "1,3,7,11,15")
dont_cares_input = st.sidebar.text_input("Don't care terms", "2,5,6,9,10")

def parse_list(s):
    return [int(x.strip()) for x in s.split(",") if x.strip()]

minterms = parse_list(minterms_input)
dont_cares = parse_list(dont_cares_input)

max_val = (1 << num_vars) - 1
for v in minterms + dont_cares:
    if v < 0 or v > max_val:
        st.error(f"Giá trị {v} không hợp lệ cho {num_vars} biến.")
        st.stop()

dont_cares = [d for d in dont_cares if d not in minterms]

if st.button("Tối giản và Vẽ K-map"):
    st.subheader("Dữ liệu đầu vào")
    c1, c2, c3 = st.columns([1, 2, 2])
    c1.write("**Số biến**"); c1.write(num_vars)
    c2.write("**Minterms**"); c2.write(minterms if minterms else "—")
    c3.write("**Don't cares**"); c3.write(dont_cares if dont_cares else "—")

    prime_info, selected_idx, final_exprs, essential_idx = solve_minimization(minterms, dont_cares, num_vars)

    colors = ['#ffcccc', '#ccffcc', '#ccccff', '#ffffcc', '#ccffff', '#ffccff',
              '#ffd9b3', '#e6ccff', '#d9ffcc', '#ffb3b3']

    st.subheader("Prime implicants (tất cả)")
    for idx, p in enumerate(prime_info):
        marker = " (essential)" if idx in essential_idx else ""
        st.markdown(
            f"<div style='padding:6px;border-radius:6px;background:#f7f7f7;margin-bottom:4px;'>"
            f"<b>Index {idx}</b>{marker} &nbsp;|&nbsp; Mask: <code>{p['mask']}</code> "
            f"&nbsp;|&nbsp; Covers: <code>{sorted(p['covers'])}</code> "
            f"&nbsp;|&nbsp; Literals: {p['literals']}</div>",
            unsafe_allow_html=True
        )

    st.subheader("Implicants được chọn")
    if not final_exprs:
        st.info("Không có implicant được chọn.")
    else:
        selected_groups = []
        for order, block in enumerate(final_exprs):
            selected_groups.append(block)
            color = colors[order % len(colors)]
            essential_mark = " (essential)" if block['idx'] in essential_idx else ""
            st.markdown(
                f"<div style='background-color:{color};padding:8px;border-radius:6px;margin-bottom:6px;'>"
                f"<b>Selected {order+1}{essential_mark}:</b> Mask: <code>{block['mask']}</code> "
                f"| Covers: <code>{block['covers']}</code> "
                f"| Expr: <b>{block['expr']}</b></div>",
                unsafe_allow_html=True
            )

        st.subheader("Biểu thức tối giản (SOP)")
        exprs = [block['expr'] for block in selected_groups]
        latex_exprs = [re.sub(r"([A-Z])'", r"\1^{\\prime}", e) for e in exprs]
        st.latex(" + ".join(latex_exprs))

        st.subheader("K-map (1 = minterm, X = don't care)")
        fig = draw_kmap(num_vars, minterms, dont_cares, selected_groups, colors)
        st.pyplot(fig)
