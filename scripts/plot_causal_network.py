"""
scripts/plot_causal_network.py - Publication-Grade Causal Network Connection Diagram
Renders the discovered causal connections between Master Regulator lncRNAs, TFs, and Effectors.
"""
import os
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np

def generate_causal_network_figure():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    edges_path = os.path.join(base_dir, "results", "tables", "stable_edges.csv")
    mrs_path = os.path.join(base_dir, "results", "tables", "mrs_rankings.csv")
    out_path = os.path.join(base_dir, "results", "figures", "causal_network_connections.png")

    df_edges = pd.read_csv(edges_path)
    df_mrs = pd.read_csv(mrs_path).set_index("symbol")

    # Define key nodes
    lncrnas = [
        "LINC00662", "LINC00539", "LINC01003", "STX18-AS1",
        "LINC00667", "CARD8-AS1", "MCM3AP-AS1", "DLGAP1-AS1"
    ]
    tfs = ["HOXA9", "MEIS1", "MYC", "RUNX1", "SPI1"]
    effectors = ["BCL2", "MCL1", "BAX", "PMAIP1", "CD34", "FLT3"]
    all_nodes = lncrnas + tfs + effectors

    # Filter edges between selected nodes
    sub_edges = df_edges[
        df_edges["source_symbol"].isin(all_nodes) & 
        df_edges["target_symbol"].isin(all_nodes)
    ].copy()

    # Create directed graph
    G = nx.DiGraph()
    for n in all_nodes:
        biotype = "lncRNA" if n in lncrnas else ("TF" if n in tfs else "Effector")
        mrs = df_mrs.loc[n, "mrs_score"] if n in df_mrs.index else 0.5
        G.add_node(n, biotype=biotype, mrs=mrs)

    for _, r in sub_edges.iterrows():
        G.add_edge(r["source_symbol"], r["target_symbol"], weight=float(r["bootstrap_freq"]))

    # Setup multipartite / layered hierarchical layout
    pos = {}
    # Left column: lncRNAs
    y_lnc = np.linspace(0.90, 0.10, len(lncrnas))
    for i, n in enumerate(lncrnas):
        pos[n] = np.array([0.12, y_lnc[i]])

    # Middle column: Transcription Factors
    y_tf = np.linspace(0.80, 0.20, len(tfs))
    for i, n in enumerate(tfs):
        pos[n] = np.array([0.50, y_tf[i]])

    # Right column: Effectors
    y_eff = np.linspace(0.85, 0.15, len(effectors))
    for i, n in enumerate(effectors):
        pos[n] = np.array([0.88, y_eff[i]])

    # Setup Plot
    fig, ax = plt.subplots(figsize=(15, 10), facecolor="#0B0F19")
    ax.set_facecolor("#0B0F19")

    # Draw edges with varying alpha & curve
    key_chains = [
        ("LINC00667", "MYC"), ("MYC", "BCL2"), ("MYC", "BAX"),
        ("LINC00539", "MEIS1"), ("MEIS1", "BCL2"), ("MEIS1", "HOXA9"),
        ("CARD8-AS1", "HOXA9"), ("HOXA9", "CD34"), ("LINC01003", "RUNX1")
    ]

    for u, v, d in G.edges(data=True):
        is_key = (u, v) in key_chains
        w = d.get("weight", 0.8)
        color = "#00F0FF" if is_key else ("#4B5563" if w < 0.85 else "#6366F1")
        alpha = 0.85 if is_key else 0.35
        lw = 2.4 if is_key else 1.0
        rad = 0.08 if is_key else 0.04
        
        arrow = nx.draw_networkx_edges(
            G, pos, edgelist=[(u, v)], ax=ax,
            edge_color=color, alpha=alpha, width=lw,
            arrows=True, arrowstyle="-|>", arrowsize=14,
            connectionstyle=f"arc3,rad={rad}",
            min_source_margin=18, min_target_margin=18
        )

    # Node styling
    node_colors = {
        "lncRNA": "#7928CA",   # Neon violet
        "TF": "#F59E0B",       # Amber gold
        "Effector": "#EF4444"  # Rose red
    }
    node_borders = {
        "lncRNA": "#C084FC",
        "TF": "#FDE68A",
        "Effector": "#FCA5A5"
    }

    for biotype, color in node_colors.items():
        nlist = [n for n, d in G.nodes(data=True) if d["biotype"] == biotype]
        sizes = [1800 if biotype == "TF" else (1600 if biotype == "lncRNA" else 1500) for _ in nlist]
        nx.draw_networkx_nodes(
            G, pos, nodelist=nlist, ax=ax,
            node_color=color, node_size=sizes,
            edgecolors=node_borders[biotype], linewidths=2.0, alpha=0.95
        )

    # Node labels
    labels = {n: n for n in G.nodes()}
    nx.draw_networkx_labels(
        G, pos, labels=labels, ax=ax,
        font_size=9.5, font_weight="bold", font_color="#FFFFFF", font_family="sans-serif"
    )

    # Subtitles for the three functional layers
    ax.text(0.12, 0.97, "Master lncRNAs\n(Pearlian Drivers)", color="#C084FC",
            fontsize=12, fontweight="bold", ha="center", va="center")
    ax.text(0.50, 0.97, "Transcription Factors\n(Lineage Regulators)", color="#FDE68A",
            fontsize=12, fontweight="bold", ha="center", va="center")
    ax.text(0.88, 0.97, "Downstream Effectors\n(Apoptosis & Stemness)", color="#FCA5A5",
            fontsize=12, fontweight="bold", ha="center", va="center")

    # Add custom legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='lncRNA Master Driver (Tier 1/2)',
               markerfacecolor='#7928CA', markeredgecolor='#C084FC', markersize=11),
        Line2D([0], [0], marker='o', color='w', label='Transcription Factor Hub',
               markerfacecolor='#F59E0B', markeredgecolor='#FDE68A', markersize=11),
        Line2D([0], [0], marker='o', color='w', label='Effector / Resistance Node',
               markerfacecolor='#EF4444', markeredgecolor='#FCA5A5', markersize=11),
        Line2D([0], [0], color='#00F0FF', lw=2.5, label='High-Impact Causal Path (do-knockout validated)'),
        Line2D([0], [0], color='#6366F1', lw=1.2, alpha=0.6, label='Discovered Causal Regulation (DAG Edge)')
    ]
    leg = ax.legend(handles=legend_elements, loc="lower center", bbox_to_anchor=(0.5, -0.05),
                    ncol=3, frameon=True, facecolor="#111827", edgecolor="#374151",
                    fontsize=9.5, labelcolor="#E5E7EB")

    ax.set_title("Causal lncRNA–mRNA Regulatory Atlas of Acute Myeloid Leukemia (AML)\nDiscovered Hierarchical Directed Connections (PC-Algorithm & RobustFisherZ)",
                 color="#F9FAFB", fontsize=15, fontweight="bold", pad=25)

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.08, 1.04)
    ax.axis("off")

    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    plt.savefig(out_path, dpi=300, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close()
    print(f"Network figure saved successfully to: {out_path}")

if __name__ == "__main__":
    generate_causal_network_figure()
