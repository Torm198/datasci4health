import pandas as pd
import numpy as np
import networkx as nx
import community as community_louvain
import re
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# =========================
# CONFIG
# =========================
LOG2FC_THRESHOLD = 1
TOP_GENES = 50

BETA = 6
WEIGHT_THRESHOLD = 0.1
MAX_EDGES = 5000

#CONTROLE
RUN_EXPRESSION = True
RUN_NETWORK = False
RUN_ADVANCED = False
EXPORT_NETWORK = False


# =========================
# 0. LEITURA ROBUSTA
# =========================
def safe_read_csv(path, sep="\t", compression="infer"):
    encodings = ["utf-8", "latin-1", "ISO-8859-1"]

    for enc in encodings:
        try:
            df = pd.read_csv(
                path,
                sep=sep,
                compression=compression,
                encoding=enc,
                on_bad_lines="skip"
            )
            print(f"Arquivo lido com encoding: {enc}")
            return df
        except:
            continue

    raise ValueError(f"Erro ao ler arquivo: {path}")


# =========================
# 1. EXPRESSÃO
# =========================
def compute_log2fc(df, control_cols, treated_cols):
    df["mean_control"] = df[control_cols].mean(axis=1)
    df["mean_treated"] = df[treated_cols].mean(axis=1)

    df["log2FC"] = np.log2(
        (df["mean_treated"] + 1) /
        (df["mean_control"] + 1)
    )
    return df


def get_top_genes(results):
    results = results[results["padj"] < 0.05]
    results = results[abs(results["log2FoldChange"]) > LOG2FC_THRESHOLD]
    results = results.sort_values(by="log2FoldChange", ascending=False)
    return results.head(TOP_GENES)


def export_expression(df, prefix):
    script_dir = os.path.dirname(__file__)
    path_to_save = os.path.join(script_dir, f"{prefix}_top_genes.csv")
    df[["log2FoldChange", "padj"]].to_csv(path_to_save, index=True)


# =========================
# 2. REDE (WGCNA)
# =========================
def compute_correlation(df, expr_cols):
    expr = df[expr_cols]
    expr.index = df.index  # assuming df.index is GeneName
    return expr.T.corr(method="spearman")


def build_edges(corr):
    edges = []

    for i in range(len(corr)):
        for j in range(i + 1, len(corr)):
            weight = abs(corr.iloc[i, j]) ** BETA

            if weight > WEIGHT_THRESHOLD:
                edges.append((corr.index[i], corr.index[j], weight))

    return sorted(edges, key=lambda x: x[2], reverse=True)[:MAX_EDGES]


def build_graph(genes, edges):
    G = nx.Graph()
    G.add_nodes_from(genes)
    G.add_weighted_edges_from(edges)
    return G


# =========================
# 3. ANÁLISE AVANÇADA
# =========================
def compute_advanced(G):
    centrality = nx.degree_centrality(G)
    partition = community_louvain.best_partition(G)
    return centrality, partition


# =========================
# 4. EXPORT REDE
# =========================
def export_network(top_genes, edges, centrality, partition, prefix):
    nodes = top_genes[["log2FoldChange"]].copy()
    nodes["id"] = nodes.index

    if centrality:
        nodes["degree"] = nodes["id"].map(centrality)

    if partition:
        nodes["community"] = nodes["id"].map(partition)

    nodes.to_csv(f"{prefix}_nodes.csv", index=False)

    edges_df = pd.DataFrame(edges, columns=["source", "target", "weight"])
    edges_df.to_csv(f"{prefix}_edges.csv", index=False)


# =========================
# PIPELINE
# =========================
def run_analysis(df, gene_col, control_cols, treated_cols, prefix):

    df = df.dropna(subset=[gene_col])
    df = df.drop_duplicates(subset=[gene_col])
    df = df.rename(columns={gene_col: "GeneName"}).set_index("GeneName")

    # Run DESeq2
    counts = df[control_cols + treated_cols].T
    metadata = pd.DataFrame({
        "condition": ["control"] * len(control_cols) + ["treated"] * len(treated_cols)
    }, index=control_cols + treated_cols)

    dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors="condition")
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=["condition", "treated", "control"], n_cpus=1)
    stat_res.summary()
    results = stat_res.results_df

    # Get top genes
    top_genes = get_top_genes(results)

    print(f"\n[{prefix}] Top genes:")
    print(top_genes[["log2FoldChange", "padj"]].head(10))

    if RUN_EXPRESSION:
        export_expression(top_genes, prefix)

    if not RUN_NETWORK:
        return

    # REDE
    expr_cols = control_cols + treated_cols
    df_top = df.loc[top_genes.index]
    corr = compute_correlation(df_top, expr_cols)
    edges = build_edges(corr)
    G = build_graph(top_genes.index, edges)

    print(f"[{prefix}] Nós:", G.number_of_nodes())
    print(f"[{prefix}] Arestas:", G.number_of_edges())
    print(f"[{prefix}] Densidade:", round(nx.density(G), 4))

    centrality = None
    partition = None

    if RUN_ADVANCED:
        centrality, partition = compute_advanced(G)

        top_hubs = sorted(
            centrality.items(),
            key=lambda x: x[1],
            reverse=True
        )[:10]

        print(f"\n[{prefix}] Top hubs:")
        for g, v in top_hubs:
            print(g, v)

        print(f"[{prefix}] Comunidades:", len(set(partition.values())))

    if EXPORT_NETWORK:
        export_network(top_genes, edges, centrality, partition, prefix)

if __name__ == "__main__":
    script_dir = os.path.dirname(__file__)
    
    # =========================
    # 1. KLEBSIELLA
    # =========================
    df_kleb = safe_read_csv(os.path.join(script_dir, "GSE307523_gene_fpkm.txt.gz"))

    control_kleb = ["NK01067-1_Count","NK01067-2_Count","NK01067-3_Count"]
    treated_kleb = ["NK01067-MEM-1_Count","NK01067-MEM-2_Count","NK01067-MEM-3_Count"]

    run_analysis(df_kleb, "GeneName", control_kleb, treated_kleb, "kleb")


    # =========================
    # 2. ACINETOBACTER
    # =========================
    df_aci = safe_read_csv(os.path.join(script_dir, "GSE190441_Complete_Raw_gene_counts_matrix.txt.gz"))

    cols_3h = [c for c in df_aci.columns if c.startswith("3h_")]
    df_aci = df_aci[["Name"] + cols_3h]

    samples = {}
    for col in cols_3h:
        base = re.sub(r"_R[12]", "", col)
        samples.setdefault(base, []).append(col)

    df_sum = pd.DataFrame()
    df_sum["Name"] = df_aci["Name"]

    for sample, cols in samples.items():
        df_sum[sample] = df_aci[cols].sum(axis=1)

    control_aci = [c for c in df_sum.columns if "no_mero" in c]
    treated_aci = [c for c in df_sum.columns if "mero" in c and "no" not in c]

    run_analysis(df_sum, "Name", control_aci, treated_aci, "aci")


    # =========================
    # 3. PSEUDOMONAS
    # =========================
    df_pa = pd.read_excel(os.path.join(script_dir, "GSE167137_P_aeruginosa_count_data.xlsx"))

    control_pa = [c for c in df_pa.columns if "PA-0M-" in c]
    treated_pa = [c for c in df_pa.columns if "PA-5M-" in c]

    run_analysis(df_pa, "Unnamed: 0", control_pa, treated_pa, "pa")