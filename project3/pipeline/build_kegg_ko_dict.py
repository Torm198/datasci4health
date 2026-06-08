import re
import time
import json
import urllib.request
import urllib.error
import urllib.parse
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

SPECIES_TAXID = {
    "pseu":  208964,  # Pseudomonas aeruginosa PAO1
    "aci":   400667,  # Acinetobacter baumannii ATCC 17978
    "kleb":  272620,  # Klebsiella pneumoniae ATCC 13883
}

STRING_MIN_SCORE = 400  # medium confidence

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

# =========================================================
# STEP 2: Parse GFF
# =========================================================

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

# =========================================================
# STEP 3: Parse KO
# =========================================================

def parse_ko(ko_path, species):
    log(f"[{species}] Lendo anotações KO: {ko_path}")
    df = pd.read_csv(ko_path, sep=r"\s+", header=None,
                     names=["protein_id", "KO"], comment="#", engine="python"
                     ).dropna(subset=["protein_id", "KO"])
    log(f"[{species}] {len(df)} anotações KO carregadas")
    return df

# =========================================================
# STEP 4: Integração por organismo
# =========================================================

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

# =========================================================
# STEP 5: Matriz de ortologia
# =========================================================

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

# =========================================================
# STEP 5b: Core ortólogo DE
# =========================================================

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

# =========================================================
# STEP 6: Rede Cytoscape
# =========================================================

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

# =========================================================
# STEP 7: KEGG download + ORA
# =========================================================

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

# =========================================================
# STEP 8: Anotação + Hubs + Plot
# =========================================================

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

# =========================================================
# REDE STRING: PPI por espécie com ortho_score e string_score
# =========================================================

def build_string_network(all_integrated, edges_df, combined_nodes, min_score=STRING_MIN_SCORE):
    """
    Gera 3 pares de arquivos separados por espécie:
      string_pseu_edges.tsv / string_pseu_nodes.tsv
      string_aci_edges.tsv  / string_aci_nodes.tsv
      string_kleb_edges.tsv / string_kleb_nodes.tsv

    Usa exatamente os mesmos genes do combined_nodes (sem hubs de pathway).
    Colunas dos nodes: ID, gene_name, Especie, KO, LFC, ortho_score, string_score
    Colunas das edges: source, target, Especie, string_interaction_score
    """
    log("Iniciando rede STRING...")

    # Filtra pelos IDs do combined_nodes excluindo hubs
    gene_ids = set(
        combined_nodes[combined_nodes["node_type"] != "pathway_hub"]["ID"].tolist()
    )
    log(f"  Genes do combined_nodes (sem hubs): {len(gene_ids)}")

    # Monta tabela base a partir do all_integrated filtrado pelos IDs do combined
    base_rows = []
    for org_id, df in all_integrated.items():
        for _, row in df.iterrows():
            if row["protein_id"] not in gene_ids:
                continue
            gname = row.get("gene_name", "") or row.get("locus_tag", row["protein_id"])
            if not gname:
                gname = row["protein_id"]
            base_rows.append({
                "ID":        row["protein_id"],
                "gene_name": gname,
                "Especie":   org_id,
                "KO":        row.get("KO", ""),
                "LFC":       row["LFC"],
            })
    base_df = pd.DataFrame(base_rows).drop_duplicates(subset=["ID"])
    log(f"  Base STRING após filtro: {len(base_df)} genes")

    # ortho_score: grau de conexões cross-espécie na rede de ortologia
    ortho_edges = edges_df[edges_df["type"] != "gene_to_pathway"] if not edges_df.empty else pd.DataFrame()
    degree = {}
    for _, row in ortho_edges.iterrows():
        degree[row["source"]] = degree.get(row["source"], 0) + 1
        degree[row["target"]] = degree.get(row["target"], 0) + 1
    base_df["ortho_score"] = base_df["ID"].map(degree).fillna(0).astype(int)
    log(f"  ortho_score: max={base_df['ortho_score'].max()} | com score>0: {(base_df['ortho_score']>0).sum()}")

    # Busca STRING e gera arquivos separados por espécie
    for species, taxid in SPECIES_TAXID.items():
        sp_df = base_df[base_df["Especie"] == species].copy()
        if sp_df.empty:
            log(f"  [{species}] nenhum gene — pulando")
            continue

        gene_names = [g for g in sp_df["gene_name"].dropna().unique() if g]
        if not gene_names:
            continue

        log(f"  [{species}] buscando {len(gene_names)} genes no STRING (taxid={taxid})...")

        str_edges = []
        string_degree = {}

        for i in range(0, len(gene_names), 100):
            batch = gene_names[i:i+100]
            identifiers = "%0d".join(batch)
            url = "https://string-db.org/api/json/network"
            params = urllib.parse.urlencode({
                "identifiers":     identifiers,
                "species":         taxid,
                "required_score":  min_score,
                "caller_identity": "ortologia_pipeline",
            }).encode("utf-8")
            try:
                req = urllib.request.Request(url, data=params)
                with urllib.request.urlopen(req, timeout=30) as r:
                    data = json.loads(r.read().decode("utf-8"))
                for inter in data:
                    gA = inter.get("preferredName_A", "")
                    gB = inter.get("preferredName_B", "")
                    sc = inter.get("score", 0)
                    if gA and gB:
                        str_edges.append({
                            "source":                   gA,
                            "target":                   gB,
                            "Especie":                  species,
                            "string_interaction_score": round(sc, 3),
                        })
                        string_degree[gA] = string_degree.get(gA, 0) + 1
                        string_degree[gB] = string_degree.get(gB, 0) + 1
            except Exception as e:
                log(f"  AVISO STRING [{species}] lote {i}: {e}")
            time.sleep(0.5)

        sp_df["string_score"] = sp_df["gene_name"].map(string_degree).fillna(0).astype(int)
        sp_df["ortho_score"]  = sp_df["ID"].map(degree).fillna(0).astype(int)

        # Normaliza scores dentro da espécie (0→1)
        def norm_col(series):
            mn, mx = series.min(), series.max()
            return (series - mn) / (mx - mn) if mx != mn else series * 0.0

        ortho_n  = norm_col(sp_df["ortho_score"])
        string_n = norm_col(sp_df["string_score"])
        lfc_abs  = sp_df["LFC"].abs()

        # impact_score = (ortho_norm + string_norm) × |LFC|
        sp_df["impact_score"] = ((ortho_n + string_n) * lfc_abs).round(4)

        str_edges_df = pd.DataFrame(str_edges)

        sp_df.to_csv(f"string_{species}_nodes.tsv", sep="\t", index=False)
        if not str_edges_df.empty:
            str_edges_df.to_csv(f"string_{species}_edges.tsv", sep="\t", index=False)
        log(f"  [{species}] {len(str_edges_df)} arestas | {len(sp_df)} nós | string_score max={sp_df['string_score'].max()}")


# =========================================================
# PLOT: Top 5 genes por espécie pelo impact_score
# =========================================================

def plot_impact_scores():
    """
    Bubble chart: X=LFC, Y=string_score, tamanho=impact_score
    Todos circulos, labels ajustados para evitar sobreposicao.
    """
    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patheffects as pe
    except ImportError:
        log("matplotlib nao instalado"); return

    try:
        from adjustText import adjust_text
        HAS_ADJUSTTEXT = True
    except ImportError:
        HAS_ADJUSTTEXT = False

    SPECIES_CONFIG = {
        "pseu":  {"label": "Pseudomonas aeruginosa",  "color": "#2E86C1"},
        "aci":   {"label": "Acinetobacter baumannii", "color": "#E67E22"},
        "kleb":  {"label": "Klebsiella pneumoniae",   "color": "#27AE60"},
    }

    fig, ax = plt.subplots(figsize=(13, 9))
    all_top = []
    texts   = []

    for species, cfg in SPECIES_CONFIG.items():
        fpath = f"string_{species}_nodes.tsv"
        if not Path(fpath).exists():
            log(f"  AVISO: {fpath} nao encontrado — pulando")
            continue

        df = pd.read_csv(fpath, sep="\t")
        needed = ["LFC", "string_score", "impact_score"]
        if not all(c in df.columns for c in needed):
            log(f"  AVISO: colunas ausentes em {fpath}")
            continue

        for col in needed:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        df = df.dropna(subset=needed)

        # Normaliza ortho + string para tamanho do ponto (20 a 400)
        def norm01(s):
            mn, mx = s.min(), s.max()
            return (s - mn) / (mx - mn) if mx != mn else s * 0.0

        # Normaliza X e Y dentro da espécie (0→1)
        ortho_n  = norm01(df["ortho_score"])
        string_n = norm01(df["string_score"])

        lfc_abs   = df["LFC"].abs()
        lfc_min   = lfc_abs.min()
        lfc_max   = lfc_abs.max()
        lfc_range = lfc_max - lfc_min if lfc_max != lfc_min else 1
        sizes     = 20 + ((lfc_abs - lfc_min) / lfc_range) * 400

        # Todos os genes — pequenos e transparentes
        ax.scatter(
            ortho_n, string_n,
            s=sizes, color=cfg["color"], alpha=0.18,
            marker="o", zorder=2,
        )

        # Top 5 por impact_score
        top5      = df.nlargest(5, "impact_score")
        ortho_top = norm01(top5["ortho_score"])
        string_top = norm01(top5["string_score"])
        lfc_top   = top5["LFC"].abs()
        sizes_top = 20 + ((lfc_top - lfc_min) / lfc_range) * 400

        ax.scatter(
            ortho_top, string_top,
            s=sizes_top, color=cfg["color"], alpha=0.92,
            marker="o", zorder=4,
            label=cfg["label"],
            edgecolors="white", linewidths=1.2,
        )

        for i, (_, row) in enumerate(top5.iterrows()):
            gname = str(row.get("gene_name", row.get("ID", "")))
            ox = ortho_top.iloc[i]
            sy = string_top.iloc[i]
            t = ax.text(
                ox, sy,
                gname,
                fontsize=9, fontweight="bold", color=cfg["color"],
                zorder=6,
                path_effects=[pe.withStroke(linewidth=3, foreground="white")],
            )
            texts.append(t)
        all_top.append(top5)

    # Ajusta labels automaticamente se adjustText disponivel
    if HAS_ADJUSTTEXT and texts:
        adjust_text(
            texts, ax=ax,
            arrowprops=dict(arrowstyle="-", color="grey", lw=0.6, alpha=0.5),
            expand_points=(1.8, 1.8),
            expand_text=(1.4, 1.4),
        )
    else:
        offsets = [(0.15, 1.5),(-0.15, 1.5),(0.15,-1.5),(-0.15,-1.5),
                   (0.25, 0),(-0.25, 0),(0.1, 2.5),(-0.1,-2.5),(0.2,1.8),(-0.2,1.8)]
        for i, t in enumerate(texts):
            x, y = t.get_position()
            off  = offsets[i % len(offsets)]
            t.set_position((x + off[0], y + off[1]))

    ax.set_xlabel("Ortho Score normalizado (conexoes cross-especie)", fontsize=13)
    ax.set_ylabel("String Score normalizado (interacoes PPI)", fontsize=13)
    ax.set_title(
        "Genes mais relevantes por especie\n"
        "tamanho = |LFC|  |  cor = especie  |  top 5 destacados por impact score",
        fontsize=13
    )
    ax.legend(fontsize=10, framealpha=0.9, markerscale=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig("impact_score_plot.png", dpi=180, bbox_inches="tight")
    plt.close()
    log("Plot salvo: impact_score_plot.png")

    if all_top:
        top_df = pd.concat(all_top, ignore_index=True)
        cols = ["Especie", "gene_name", "KO", "LFC", "ortho_score", "string_score", "impact_score"]
        top_df[[c for c in cols if c in top_df.columns]]\
            .sort_values(["Especie", "impact_score"], ascending=[True, False])\
            .to_csv("top5_impact_genes.tsv", sep="\t", index=False)
        log("Tabela salva: top5_impact_genes.tsv")

# MAIN
# =========================================================

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

    edges_df = annotate_with_pathway(edges_df, ora_df)
    nodes_df = annotate_with_pathway(nodes_df, ora_df, fill="")
    edges_df.to_csv("rede_ortologia_edges.tsv", sep="\t", index=False)
    nodes_df.to_csv("nodes_metadata.tsv",       sep="\t", index=False)

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

    star_edges, hub_nodes = build_pathway_star_edges(pw_nodes, pathway_to_info, ora_df)

    if not star_edges.empty:
        combined_nodes = pd.concat([pw_nodes, hub_nodes], ignore_index=True)
        combined_edges = pd.concat([pw_edges, star_edges], ignore_index=True)
        combined_edges.to_csv("combined_edges.tsv", sep="\t", index=False)
        combined_nodes.to_csv("combined_nodes.tsv", sep="\t", index=False)
        hub_nodes.to_csv("pathway_hub_nodes.tsv",   sep="\t", index=False)
        log(f"Rede combinada: {len(combined_edges)} arestas | {len(combined_nodes)} nós")

    # Saída 6: Rede STRING (PPI por espécie com ortho_score e string_score)
    build_string_network(all_integrated, edges_df, combined_nodes)
    plot_impact_scores()

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
    print("  string_pseu_edges/nodes.tsv    <- rede PPI Pseudomonas")
    print("  string_aci_edges/nodes.tsv     <- rede PPI Acinetobacter")
    print("  string_kleb_edges/nodes.tsv    <- rede PPI Klebsiella")
    print("  impact_score_plot.png          <- top 5 genes por espécie")
    print("  top5_impact_genes.tsv          <- tabela top 5 genes")


if __name__ == "__main__":
    main()