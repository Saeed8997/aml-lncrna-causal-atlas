"""
dashboard/app.py - AML Causal lncRNA-mRNA Regulatory Atlas Web App (< 80 lines)
Nature-Grade Bioinformatics Web Platform with D3.js Causal Simulator & Multi-Drug Gauges.
"""
import os
import sys
import streamlit as st
import streamlit.components.v1 as components

root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
dash_dir = os.path.dirname(__file__)
for p in (root_dir, dash_dir):
    if p not in sys.path:
        sys.path.insert(0, p)

try:
    from dashboard.data_bridge import build_network_payload, load_table
    from dashboard.network_canvas import render_d3_causal_network
    from dashboard.analytics_tabs import render_multi_drug_tab, render_novelty_and_rewiring_tab, render_crispri_gallery_tab
except ImportError:
    from data_bridge import build_network_payload, load_table
    from network_canvas import render_d3_causal_network
    from analytics_tabs import render_multi_drug_tab, render_novelty_and_rewiring_tab, render_crispri_gallery_tab

st.set_page_config(page_title="AML Causal lncRNA Atlas", page_icon="🧬", layout="wide", initial_sidebar_state="collapsed")

st.markdown("""
<style>
  .block-container { padding-top: 1rem; padding-bottom: 2rem; max-width: 98%; }
  .hero-hdr { display: flex; justify-content: space-between; align-items: flex-end; margin-bottom: 12px; }
  .hero-title { font-size: 1.65rem; font-weight: 800; background: linear-gradient(90deg, #38bdf8, #818cf8, #c084fc); -webkit-background-clip: text; -webkit-text-fill-color: transparent; margin: 0; }
  .hero-sub { color: #94a3b8; font-size: 0.85rem; margin-top: 2px; }
  .badge-tag { background: rgba(99,102,241,0.15); border: 1px solid #6366f1; color: #a5b4fc; font-size: 0.72rem; padding: 2px 8px; border-radius: 9999px; font-weight: 600; }
  .kpi-row { display: grid; grid-template-columns: repeat(5, 1fr); gap: 10px; margin-bottom: 14px; }
  .kpi-card { background: #0c1222; border: 1px solid #1e293b; border-radius: 8px; padding: 9px 13px; box-shadow: 0 4px 15px rgba(0,0,0,0.3); transition: transform 0.2s; }
  .kpi-card:hover { transform: translateY(-2px); border-color: #38bdf8; }
  .kpi-val { font-size: 1.3rem; font-weight: 700; color: #38bdf8; margin: 0; }
  .kpi-lbl { font-size: 0.72rem; color: #94a3b8; text-transform: uppercase; letter-spacing: 0.4px; }
</style>
<div class="hero-hdr">
  <div>
    <h1 class="hero-title">🧬 Causal lncRNA–mRNA Regulatory Atlas of Acute Myeloid Leukemia</h1>
    <div class="hero-sub">Multi-Resolution Constraint-Based Discovery on Single-Cell Metacell Atlases with In-Silico do-Calculus & CRISPRi Validation</div>
  </div>
  <div style="text-align:right;">
    <span class="badge-tag">CELLxGENE 748,679 Cells</span>
    <span class="badge-tag">2,476 Metacells</span>
    <span class="badge-tag">Version 8 Production</span>
  </div>
</div>
<div class="kpi-row">
  <div class="kpi-card"><div class="kpi-val">748,679</div><div class="kpi-lbl">Cells Profiled (CELLxGENE Atlas)</div></div>
  <div class="kpi-card"><div class="kpi-val">16,782</div><div class="kpi-lbl">Causal Edges (9,020 Directed)</div></div>
  <div class="kpi-card"><div class="kpi-val">516</div><div class="kpi-lbl">Tier-1 Novel lncRNAs (0 AML Papers)</div></div>
  <div class="kpi-card"><div class="kpi-val">10,625</div><div class="kpi-lbl">LSC vs. Monocytic Rewiring Edges</div></div>
  <div class="kpi-card"><div class="kpi-val">5 Drug Axes</div><div class="kpi-lbl">Multi-Drug Sensitivity Profiles</div></div>
</div>
""", unsafe_allow_html=True)

@st.cache_data
def get_cached_payload():
    return build_network_payload()

payload = get_cached_payload()
components.html(render_d3_causal_network(payload), height=705, scrolling=False)

tab_drug, tab_rewire, tab_crispr = st.tabs([
    "🌊 3-Hop Causal Cascades & Multi-Drug Matrix",
    "💎 Tier-1 Novel Targets & State Rewiring",
    "🔬 CRISPRi Functional Validation & Publication Gallery"
])

with tab_drug:
    render_multi_drug_tab(load_table("ko_effects.csv"), load_table("novel_lncrna_candidates.csv"))

with tab_rewire:
    render_novelty_and_rewiring_tab(load_table("novel_lncrna_candidates.csv"), load_table("state_rewiring.csv"))

with tab_crispr:
    render_crispri_gallery_tab(load_table("crispri_validation.csv"))
