import re
import time
import json
import urllib.request
import urllib.error
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict

DEBUG = True

def log(msg):
    if DEBUG:
        print(f"[DEBUG] {time.strftime('%H:%M:%S')} - {msg}")

# =========================================================
# CONFIG
# =========================================================
LFC_THRESHOLD   = 1.0
KEGG_CACHE_FILE = "kegg_ko_pathway_cache.json"
ORA_PVAL_CUTOFF = 0.05

SKIP_PATHWAYS = {
    'Metabolic pathways',
    'Biosynthesis of secondary metabolites',
    'Carbon metabolism',
    'Biosynthesis of amino acids',
}

ORGANISMS = {
    "pseu": {"gff": "pseu_gem/genomic.gff", "ko": "pseu_annotation.txt", "geo": "pseudo.xlsx"},
    "aci":  {"gff": "aci_gem/genomic.gff",  "ko": "aci_annotation.txt",  "geo": "acineto.txt.gz"},
    "kleb": {"gff": "kleb_gem/genomic.gff", "ko": "kleb_annotation.txt", "geo": "kleb.txt.gz"},
}

# =========================================================
# STEP 1: Cálculo de LFC (log2 simples)
# =========================================================

def _lfc_from_counts(df_counts, ctrl_cols, treat_cols, gene_col):
    """LFC = mean(log2(tratamento+1)) - mean(log2(controle+1))"""
    ctrl  = df_counts[ctrl_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    treat = df_counts[treat_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    lfc   = np.log2(treat + 1).mean(axis=1) - np.log2(ctrl + 1).mean(axis=1)
    result = pd.DataFrame({
        "gene_id_raw": df_counts[gene_col].astype(str).str.strip(),
        "LFC": lfc.values,
    })
    return result[result["LFC"] != 0].reset_index(drop=True)


def calc_lfc_pseudo(geo_path):
    log("[pseu] Lendo pseudo.xlsx e calculando LFC")
    df = pd.read_excel(geo_path)
    df.columns = [c.strip() for c in df.columns]
    if df.columns[0].startswith("Unnamed"):
        df = df.rename(columns={df.columns[0]: "Gene ID"})
    ctrl  = [c for c in df.columns if "0M" in c]
    treat = [c for c in df.columns if "5M" in c]
    log(f"[pseu] Controle ({len(ctrl)}): {ctrl}")
    log(f"[pseu] Tratamento ({len(treat)}): {treat}")
    return _lfc_from_counts(df, ctrl, treat, "Gene ID")


def calc_lfc_acineto(geo_path):
    log("[aci] Lendo acineto.txt.gz e calculando LFC")
    df = pd.read_csv(geo_path, sep="\t", compression="gzip")
    df.columns = [c.strip() for c in df.columns]
    gene_col = "Name"
    cols = [c for c in df.columns if c.startswith("3h_") or c == gene_col]
    df = df[cols].copy()
    pat = re.compile(r"^(3h_.+?)_(\d+)_(R1|R2)$")
    totals = {}
    for col in df.columns:
        if col == gene_col:
            continue
        m = pat.match(col)
        if m:
            key = f"{m.group(1)}_{m.group(2)}_Total"
            vals = df[col].apply(pd.to_numeric, errors="coerce").fillna(0).values
            totals[key] = totals.get(key, np.zeros(len(df))) + vals
    total_df = pd.DataFrame(totals)
    total_df.insert(0, gene_col, df[gene_col].values)
    log(f"[aci] Colunas após soma R1+R2: {list(total_df.columns)}")
    ctrl  = [c for c in total_df.columns if "no_mero" in c]
    treat = [c for c in total_df.columns if "mero" in c and "no_mero" not in c]
    log(f"[aci] Controle ({len(ctrl)}): {ctrl}")
    log(f"[aci] Tratamento ({len(treat)}): {treat}")
    return _lfc_from_counts(total_df, ctrl, treat, gene_col)


def calc_lfc_kleb(geo_path):
    log("[kleb] Lendo kleb.txt.gz e calculando LFC")
    df = pd.read_csv(geo_path, sep="\t", compression="gzip",
                     encoding="utf-8", encoding_errors="replace")
    df.columns = [c.strip() for c in df.columns]
    counts = [c for c in df.columns if c.endswith("_Count")]
    ctrl   = [c for c in counts if "NK01067" in c and "MEM" not in c]
    treat  = [c for c in counts if "NK01067" in c and "MEM" in c]
    log(f"[kleb] Controle ({len(ctrl)}): {ctrl}")
    log(f"[kleb] Tratamento ({len(treat)}): {treat}")
    lfc_df = _lfc_from_counts(df, ctrl, treat, "GeneId")
    if "GeneName" in df.columns:
        gene_names = df["GeneName"].astype(str).str.strip()
        lfc_df["gene_name_alt"] = gene_names.iloc[lfc_df.index].values
    return lfc_df


LFC_CALCULATORS = {
    "pseu": calc_lfc_pseudo,
    "aci":  calc_lfc_acineto,
    "kleb": calc_lfc_kleb,
}

def parse_gff(gff_path, species):
    log(f"[{species}] Lendo GFF: {gff_path}")

    def get_attr(attrs, key):
        m = re.search(rf"(?:^|;){key}=([^;]+)", attrs)
        return m.group(1) if m else ""

    gene_info = {}
    records   = []

    with open(gff_path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            if parts[2] == "gene":
                gid = get_attr(attrs, "ID")
                if gid:
                    gene_info[gid] = {
                        "locus_tag": get_attr(attrs, "locus_tag"),
                        "gene_name": get_attr(attrs, "gene"),
                    }
            elif parts[2] == "CDS":
                pid = get_attr(attrs, "protein_id")
                if not pid:
                    continue
                lt = get_attr(attrs, "locus_tag")
                gn = get_attr(attrs, "gene")
                gi = get_attr(attrs, "ID")
                if not lt or not gn:
                    parent = get_attr(attrs, "Parent")
                    if parent in gene_info:
                        lt = lt or gene_info[parent]["locus_tag"]
                        gn = gn or gene_info[parent]["gene_name"]
                records.append({"protein_id": pid, "locus_tag": lt,
                                 "gene_name": gn, "gene_id": gi})

    df = pd.DataFrame(records)
    log(f"[{species}] {len(df)} CDS extraídas do GFF")
    return df

def parse_ko(ko_path, species):
    log(f"[{species}] Lendo anotações KO: {ko_path}")
    df = pd.read_csv(ko_path, sep=r"\s+", header=None,
                     names=["protein_id", "KO"], comment="#", engine="python"
                     ).dropna(subset=["protein_id", "KO"])
    log(f"[{species}] {len(df)} anotações KO carregadas")
    return df

def integrate_organism(org_id, paths):
    lfc_df = LFC_CALCULATORS[org_id](paths["geo"])
    log(f"[{org_id}] {len(lfc_df)} genes com LFC calculado")

    gff    = parse_gff(paths["gff"], org_id)
    ko     = parse_ko(paths["ko"], org_id)
    merged = gff.merge(ko, on="protein_id", how="left")

    lfc_by_id  = dict(zip(lfc_df["gene_id_raw"], lfc_df["LFC"]))
    lfc_by_alt = {}
    if "gene_name_alt" in lfc_df.columns:
        lfc_by_alt = dict(zip(lfc_df["gene_name_alt"], lfc_df["LFC"]))

    def lookup_lfc(row):
        for key in [row["gene_name"], row["locus_tag"], row["protein_id"]]:
            if key and key in lfc_by_id:
                return lfc_by_id[key]
        for key in [row["gene_name"], row["locus_tag"]]:
            if key and key in lfc_by_alt:
                return lfc_by_alt[key]
        return None

    merged["LFC"] = merged.apply(lookup_lfc, axis=1)
    merged.insert(0, "organism", org_id)

    n_ko  = merged["KO"].notna().sum()
    n_lfc = merged["LFC"].notna().sum()
    log(f"[{org_id}] Integração: {n_ko} com KO | {n_lfc} com LFC")

    if n_lfc == 0:
        log(f"[{org_id}] AVISO: 0 matches de LFC!")
        log(f"  GEO (5 ex): {lfc_df['gene_id_raw'].head(5).tolist()}")
        log(f"  GFF locus_tag (5 ex): {merged['locus_tag'].dropna().head(5).tolist()}")
        log(f"  GFF gene_name (5 ex): {merged['gene_name'].dropna().head(5).tolist()}")

    ko_to_genes = defaultdict(list)
    for _, row in merged.dropna(subset=["KO"]).iterrows():
        ko_to_genes[row["KO"]].append(row["protein_id"])

    return merged, ko_to_genes


def build_orthology_matrix(all_species_kos):
    log("Construindo matriz de ortologia")
    all_kos = sorted(set(ko for d in all_species_kos.values() for ko in d))
    data = []
    for ko in all_kos:
        row = {"KO_ID": ko}
        for sp, sp_dict in all_species_kos.items():
            row[sp] = ";".join(sp_dict.get(ko, [])) or "-"
        data.append(row)
    return pd.DataFrame(data)

def build_core_de(all_integrated, species_list, threshold=LFC_THRESHOLD):
    log(f"Calculando core ortólogo DE (threshold |LFC| >= {threshold})")
    ko_de_per_sp = {}
    for org_id, df in all_integrated.items():
        ko_map = defaultdict(list)
        for _, row in df.iterrows():
            if pd.isna(row["LFC"]) or pd.isna(row["KO"]):
                continue
            if abs(row["LFC"]) >= threshold:
                name = row["gene_name"] if row.get("gene_name") else row.get("locus_tag", row["protein_id"])
                ko_map[row["KO"]].append((name, row["LFC"], row["protein_id"]))
        ko_de_per_sp[org_id] = ko_map
        log(f"  [{org_id}] {len(ko_map)} KOs com |LFC| >= {threshold}")

    sets = [set(d.keys()) for d in ko_de_per_sp.values()]
    core_kos = set.intersection(*sets) if sets else set()
    log(f"  Core ortólogo DE: {len(core_kos)} KOs compartilhados entre todas as espécies")

    rows = []
    for ko in sorted(core_kos):
        row = {"KO": ko}
        for sp in species_list:
            if sp not in ko_de_per_sp:
                continue
            hits = ko_de_per_sp[sp].get(ko, [])
            row[f"{sp}_gene"]       = ";".join(g for g, _, _p in hits)
            row[f"{sp}_LFC"]        = ";".join(f"{l:.4f}" for _, l, _p in hits)
            row[f"{sp}_protein_id"] = ";".join(p for _, _, p in hits)
        rows.append(row)

    return pd.DataFrame(rows)

def build_cytoscape_network(ortho_df, all_integrated, species_list):
    log("Gerando rede para Cytoscape")

    lfc_lu  = {}
    name_lu = {}
    for org_id, df in all_integrated.items():
        for _, row in df.iterrows():
            pid = row["protein_id"]
            lfc_lu[(org_id, pid)]  = row["LFC"]
            name_lu[(org_id, pid)] = row["gene_name"] if row.get("gene_name") else row.get("locus_tag", "")

    sps   = [s for s in species_list if s in ortho_df.columns]
    edges = []
    nodes = []

    for _, row in ortho_df.iterrows():
        ko = row["KO_ID"]
        for sp in sps:
            if row[sp] == "-":
                continue
            for pid in row[sp].split(";"):
                lfc = lfc_lu.get((sp, pid))
                nodes.append({
                    "ID":        pid,
                    "gene_name": name_lu.get((sp, pid), pid),
                    "Especie":   sp,
                    "KO":        ko,
                    "LFC":       lfc if lfc is not None else "",
                    "DE":        "yes" if lfc is not None else "no",
                    "node_type": "gene",
                })
        sps_ok = [sp for sp in sps if row[sp] != "-"]
        for i in range(len(sps_ok)):
            for j in range(i + 1, len(sps_ok)):
                si, sj = sps_ok[i], sps_ok[j]
                for gi in row[si].split(";"):
                    for gj in row[sj].split(";"):
                        edges.append({
                            "source": gi, "target": gj, "KO": ko,
                            "type": f"{si}_{sj}",
                            "LFC_source": lfc_lu.get((si, gi)),
                            "LFC_target": lfc_lu.get((sj, gj)),
                            "top_pathway": "",
                        })

    edges_df = pd.DataFrame(edges)
    nodes_df = pd.DataFrame(nodes).drop_duplicates(subset=["ID"])
    log(f"Rede: {len(edges_df)} arestas | {len(nodes_df)} nós")
    return edges_df, nodes_df

def download_kegg_ko_pathway_map(cache_file=KEGG_CACHE_FILE):
    cache_path = Path(cache_file)
    if cache_path.exists():
        log(f"KEGG cache encontrado: {cache_file} — carregando localmente")
        with open(cache_path) as f:
            data = json.load(f)
        return data["ko_to_pathways"], data["pathway_to_info"]

    log("KEGG cache nao encontrado — iniciando download completo (unica vez)...")

    def kegg_get(url, retries=3, delay=1.0):
        for attempt in range(retries):
            try:
                with urllib.request.urlopen(url, timeout=30) as r:
                    return r.read().decode("utf-8")
            except Exception as e:
                if attempt < retries - 1:
                    time.sleep(delay)
                else:
                    raise e

    log("  Baixando lista de pathways KEGG...")
    raw = kegg_get("https://rest.kegg.jp/list/pathway")
    pathway_ids = []
    pathway_to_info = {}
    for line in raw.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        pid  = parts[0].replace("path:", "").strip()
        name = parts[1].strip() if len(parts) > 1 else pid
        pathway_ids.append(pid)
        pathway_to_info[pid] = {"name": name, "kos": []}
    log(f"  {len(pathway_ids)} pathways encontrados")

    ko_to_pathways = defaultdict(list)
    total = len(pathway_ids)
    for i in range(0, total, 10):
        batch = pathway_ids[i:i+10]
        url = "https://rest.kegg.jp/link/ko/" + "+".join(batch)
        try:
            raw = kegg_get(url)
            for line in raw.strip().split("\n"):
                if not line or "\t" not in line:
                    continue
                pid_raw, ko_raw = line.split("\t")
                pid = pid_raw.replace("path:", "").strip()
                ko  = ko_raw.replace("ko:", "").strip()
                if pid in pathway_to_info:
                    pathway_to_info[pid]["kos"].append(ko)
                ko_to_pathways[ko].append(pid)
        except Exception as e:
            log(f"  AVISO: erro no lote {i}: {e}")
        if (i // 10) % 20 == 0:
            log(f"  Progresso: {min(i+10, total)}/{total} pathways")

    pathway_to_info = {p: v for p, v in pathway_to_info.items() if v["kos"]}
    ko_to_pathways  = dict(ko_to_pathways)
    with open(cache_path, "w") as f:
        json.dump({"ko_to_pathways": ko_to_pathways, "pathway_to_info": pathway_to_info}, f)
    log(f"  Cache salvo: {len(ko_to_pathways)} KOs | {len(pathway_to_info)} pathways")
    return ko_to_pathways, pathway_to_info


def run_ora(query_kos, background_kos, ko_to_pathways, pathway_to_info,
            pval_cutoff=ORA_PVAL_CUTOFF):
    q = set(query_kos)
    b = set(background_kos)
    results = []
    for pid, info in pathway_to_info.items():
        pkos = set(info["kos"])
        a = len(q & pkos)
        if a == 0:
            continue
        b2 = len(b & pkos) - a
        c  = len(q) - a
        d  = len(b) - a - b2 - c
        _, pval = fisher_exact([[a, b2], [c, d]], alternative="greater")
        results.append({
            "pathway_id":   pid,
            "pathway_name": info["name"],
            "n_query":      a,
            "n_background": a + b2,
            "n_total_ko":   len(pkos),
            "pval":         pval,
            "query_kos":    ";".join(sorted(q & pkos)),
        })
    if not results:
        log("ORA: nenhum pathway com overlap")
        return pd.DataFrame()
    df = pd.DataFrame(results)
    _, padj, _, _ = multipletests(df["pval"], method="fdr_bh")
    df["pval_adj"] = padj
    df = df[df["pval_adj"] <= pval_cutoff].sort_values("pval_adj")
    log(f"ORA: {len(df)} pathways significativos (padj <= {pval_cutoff})")
    return df.reset_index(drop=True)


def annotate_with_pathway(df, ora_df, ko_col="KO", fill="Other", top_n=5):
    top = ora_df[~ora_df["pathway_name"].isin(SKIP_PATHWAYS)].head(top_n)
    mapping = {}
    for _, row in top.iterrows():
        for ko in row["query_kos"].split(";"):
            ko = ko.strip()
            if ko and ko not in mapping:
                mapping[ko] = row["pathway_name"]
    df = df.copy()
    df["top_pathway"] = df[ko_col].map(mapping).fillna(fill)
    return df


def build_pathway_star_edges(nodes_df, pathway_to_info, ora_df, top_n=5):
    top = ora_df[~ora_df["pathway_name"].isin(SKIP_PATHWAYS)].head(top_n)
    edges = []
    hub_nodes = []
    for _, pw_row in top.iterrows():
        pid          = pw_row["pathway_id"]
        pathway_name = pw_row["pathway_name"]
        hub_id       = f"HUB_{pid}"
        if pid not in pathway_to_info:
            continue
        pathway_kos = set(pathway_to_info[pid]["kos"])
        hub_nodes.append({
            "ID": hub_id, "gene_name": pathway_name,
            "Especie": "pathway_hub", "KO": pid,
            "LFC": "", "DE": "no",
            "top_pathway": pathway_name, "node_type": "pathway_hub",
        })
        for _, gr in nodes_df[nodes_df["KO"].isin(pathway_kos)].iterrows():
            edges.append({
                "source": gr["ID"], "target": hub_id,
                "KO": gr["KO"], "type": "gene_to_pathway",
                "top_pathway": pathway_name,
                "LFC_source": gr.get("LFC", ""), "LFC_target": "",
            })
    log(f"Estrela de pathway: {len(hub_nodes)} hubs | {len(edges)} arestas gene->pathway")
    return pd.DataFrame(edges), pd.DataFrame(hub_nodes)


def plot_enrichment(ora_df, output_file="pathway_enrichment_plot.png", top_n=15):
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
    except ImportError:
        log("matplotlib nao instalado (pip install matplotlib)"); return
    df = ora_df[~ora_df["pathway_name"].isin(SKIP_PATHWAYS)].head(top_n).copy()
    df = df.sort_values("pval_adj", ascending=False)
    df["x"]     = -np.log10(df["pval_adj"].clip(lower=1e-300))
    df["label"] = df["pathway_name"].str[:50]
    norm   = plt.Normalize(df["n_query"].min(), df["n_query"].max())
    colors = cm.YlOrRd(norm(df["n_query"].values))
    fig, ax = plt.subplots(figsize=(10, max(5, len(df) * 0.45)))
    ax.barh(df["label"], df["x"], color=colors, edgecolor="grey", linewidth=0.5)
    ax.axvline(-np.log10(0.05), color="red", linestyle="--", linewidth=1, label="padj=0.05")
    ax.set_xlabel("-log10(p-value ajustado)", fontsize=11)
    ax.set_title("Enriquecimento de Pathways KEGG\nCore ortologs DE - todas as especies", fontsize=12)
    ax.legend(fontsize=9)
    sm = cm.ScalarMappable(cmap="YlOrRd", norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax).set_label("N genes DE no pathway", fontsize=9)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"Gráfico salvo: {output_file}")

def main():
    log("===== INICIANDO PIPELINE DE ORTOLOGIA =====")

    all_integrated  = {}
    all_species_kos = {}

    for org_id, paths in ORGANISMS.items():
        missing = [k for k, v in paths.items() if not Path(v).exists()]
        if missing:
            log(f"[{org_id}] AVISO: arquivos ausentes {missing} — pulando.")
            continue
        df_int, ko_genes = integrate_organism(org_id, paths)
        all_integrated[org_id]  = df_int
        all_species_kos[org_id] = ko_genes

    if not all_integrated:
        log("Nenhum organismo processado."); return

    # Saída 1: Tabela integrada
    integrated_all = pd.concat(all_integrated.values(), ignore_index=True)
    col_order = ["organism", "gene_id", "gene_name", "locus_tag", "protein_id", "KO", "LFC"]
    integrated_all[[c for c in col_order if c in integrated_all.columns]]\
        .to_csv("orthologs_integrated.tsv", sep="\t", index=False)
    log(f"Tabela integrada: orthologs_integrated.tsv ({len(integrated_all)} linhas)")

    # Saída 2: Matriz de ortologia
    ortho_df = build_orthology_matrix(all_species_kos)
    ortho_df.to_csv("dicionario_ortologia_kegg.tsv", sep="\t", index=False)
    log(f"Dicionário: dicionario_ortologia_kegg.tsv ({len(ortho_df)} grupos KO)")

    # Saída 3: Rede Cytoscape
    edges_df, nodes_df = build_cytoscape_network(ortho_df, all_integrated, list(ORGANISMS.keys()))
    edges_df.to_csv("rede_ortologia_edges.tsv", sep="\t", index=False)
    nodes_df.to_csv("nodes_metadata.tsv",       sep="\t", index=False)

    # Saída 4: Core ortólogo DE
    core_df = build_core_de(all_integrated, list(ORGANISMS.keys()))
    core_df.to_csv("core_ortologs_DE.tsv", sep="\t", index=False)
    log(f"Core DE: core_ortologs_DE.tsv ({len(core_df)} KOs)")

    core_ids = set()
    for sp in ORGANISMS:
        col = f"{sp}_protein_id"
        if col in core_df.columns:
            for pids in core_df[col].dropna():
                core_ids.update(p.strip() for p in pids.split(";") if p.strip())
    with open("core_node_ids.txt", "w") as f:
        f.write("\n".join(sorted(core_ids)))
    log(f"IDs para Cytoscape: core_node_ids.txt ({len(core_ids)} proteínas)")

    if not edges_df.empty and len(core_ids) > 0:
        core_edges = edges_df[edges_df["source"].isin(core_ids) & edges_df["target"].isin(core_ids)]
        core_nodes = nodes_df[nodes_df["ID"].isin(core_ids)]
        core_edges.to_csv("core_edges.tsv", sep="\t", index=False)
        core_nodes.to_csv("core_nodes.tsv", sep="\t", index=False)
        log(f"Subgrafo core: {len(core_edges)} arestas | {len(core_nodes)} nós")
    else:
        log("AVISO: core vazio — core_edges/nodes nao criados")
        core_nodes = pd.DataFrame()
        core_edges = pd.DataFrame()

    # Saída 5: KEGG + ORA
    log("Iniciando análise de enriquecimento KEGG...")
    ko_to_pathways, pathway_to_info = download_kegg_ko_pathway_map()

    background_kos = set(integrated_all["KO"].dropna().unique())

    # Query: uniao de todos os KOs DE de qualquer especie
    query_kos = set()
    for org_id, df in all_integrated.items():
        query_kos.update(df[df["LFC"].abs() >= LFC_THRESHOLD]["KO"].dropna().unique())
    log(f"ORA query (uniao): {len(query_kos)} KOs DE")

    if not query_kos:
        log("Nenhum KO DE encontrado — ORA ignorada"); return

    ora_df = run_ora(query_kos, background_kos, ko_to_pathways, pathway_to_info)

    if ora_df.empty:
        log("Nenhum pathway significativo — tente reduzir LFC_THRESHOLD"); return

    ora_df.to_csv("pathway_enrichment.tsv", sep="\t", index=False)
    log(f"Enriquecimento: pathway_enrichment.tsv ({len(ora_df)} pathways)")

    # Anota rede completa
    edges_df = annotate_with_pathway(edges_df, ora_df)
    nodes_df = annotate_with_pathway(nodes_df, ora_df, fill="")
    edges_df.to_csv("rede_ortologia_edges.tsv", sep="\t", index=False)
    nodes_df.to_csv("nodes_metadata.tsv",       sep="\t", index=False)

    # Subgrafo pathway + LFC
    de_ids = set(nodes_df[
        (nodes_df["LFC"].apply(lambda x: pd.to_numeric(x, errors="coerce")).abs() >= LFC_THRESHOLD) &
        (nodes_df["top_pathway"] != "")
    ]["ID"])

    pw_edges = edges_df[
        (edges_df["top_pathway"] != "Other") &
        (edges_df["source"].isin(de_ids) | edges_df["target"].isin(de_ids))
    ].copy()
    pw_node_ids = set(pw_edges["source"]) | set(pw_edges["target"])
    pw_nodes    = nodes_df[nodes_df["ID"].isin(pw_node_ids)].copy()
    pw_edges.to_csv("pathway_edges.tsv", sep="\t", index=False)
    pw_nodes.to_csv("pathway_nodes.tsv", sep="\t", index=False)
    log(f"Subgrafo pathway+LFC: {len(pw_edges)} arestas | {len(pw_nodes)} nós")

    # Estrela de pathway (hubs)
    star_edges, hub_nodes = build_pathway_star_edges(pw_nodes, pathway_to_info, ora_df)

    if not star_edges.empty:
        combined_nodes = pd.concat([pw_nodes, hub_nodes], ignore_index=True)
        combined_edges = pd.concat([pw_edges, star_edges],  ignore_index=True)
        combined_edges.to_csv("combined_edges.tsv", sep="\t", index=False)
        combined_nodes.to_csv("combined_nodes.tsv", sep="\t", index=False)
        hub_nodes.to_csv("pathway_hub_nodes.tsv",   sep="\t", index=False)
        log(f"Rede combinada: {len(combined_edges)} arestas | {len(combined_nodes)} nós")

    plot_enrichment(ora_df)

    log("===== PIPELINE CONCLUÍDO =====")
    print("\nArquivos gerados:")
    print("  orthologs_integrated.tsv       <- tabela completa por gene")
    print("  dicionario_ortologia_kegg.tsv  <- matriz KO x espécie")
    print("  rede_ortologia_edges.tsv       <- rede completa")
    print("  nodes_metadata.tsv             <- nós completos")
    print("  core_ortologs_DE.tsv           <- KOs compartilhados nas 3 espécies")
    print("  core_edges/nodes.tsv           <- subgrafo do core")
    print("  pathway_enrichment.tsv         <- pathways enriquecidos (ORA)")
    print("  pathway_enrichment_plot.png    <- gráfico de barras")
    print("  pathway_edges/nodes.tsv        <- subgrafo por pathway")
    print("  combined_edges/nodes.tsv       <- rede final para Cytoscape")
    print("  pathway_hub_nodes.tsv          <- nós hub de pathway")


if __name__ == "__main__":
    main()