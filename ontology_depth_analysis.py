# go_analysis.py

import obonet
import networkx as nx
import pandas as pd
import mygene
from itertools import combinations

def load_go_ontology(obo_file="go-basic.obo"):
    print("Loading GO ontology...")
    graph = obonet.read_obo(obo_file)
    print(f"Loaded {len(graph)} terms.")
    return graph

def compute_namespace_depths(graph):
    """Compute GO term depth from each root (BP, CC, MF) by reverse traversal."""
    roots = {
        "BP": "GO:0008150",
        "CC": "GO:0005575",
        "MF": "GO:0003674",
    }

    namespace_depths = {}
    for ns, root in roots.items():
        visited = {}
        queue = [(root, 0)]
        while queue:
            term, depth = queue.pop(0)
            if term in visited:
                continue
            visited[term] = depth
            for child in graph.predecessors(term):  # reverse traversal
                queue.append((child, depth + 1))
        namespace_depths[ns] = visited
    return namespace_depths

def get_namespace_map(graph):
    return {
        term: data.get("namespace")
        for term, data in graph.nodes(data=True)
    }

def get_go_term_info(term_id, graph, namespace_depths):
    node = graph.nodes.get(term_id, {})
    annotated_ns = node.get("namespace")
    depths = {
        ns: namespace_depths[ns].get(term_id)
        for ns in namespace_depths
        if term_id in namespace_depths[ns]
    }
    return {
        "term_id": term_id,
        "name": node.get("name"),
        "annotated_namespace": annotated_ns,
        "depths": depths
    }

def shared_go_term_with_min_depth(gene_a, gene_b, go_map, namespace_depths, namespace="CC", min_depth=3):
    terms_a = go_map.get(gene_a, {}).get(namespace, set())
    terms_b = go_map.get(gene_b, {}).get(namespace, set())
    shared_terms = terms_a & terms_b
    for term in shared_terms:
        depth = namespace_depths.get(namespace, {}).get(term, 0)
        if depth >= min_depth:
            return True
    return False

def find_deepest_shared_go_match(gene_a, gene_b, go_map, namespace_depths, namespace="BP", max_depth=15):
    last_matching_depth = None
    for depth in range(max_depth + 1):
        if shared_go_term_with_min_depth(gene_a, gene_b, go_map, namespace_depths, namespace, depth):
            last_matching_depth = depth
        else:
            break
    return last_matching_depth


mg = mygene.MyGeneInfo()

def fetch_go_terms_for_gene(gene_symbol):
    result = mg.query(gene_symbol, scopes="symbol", species="human", fields="go")
    bp_terms, cc_terms = set(), set()
    if result and "hits" in result:
        for hit in result["hits"]:
            go_data = hit.get("go", {})
            for entry in go_data.get("BP", []):
                bp_terms.add(entry["id"])
            for entry in go_data.get("CC", []):
                cc_terms.add(entry["id"])
    return {"BP": bp_terms, "CC": cc_terms}

def build_go_map(gene_list):
    return {gene: fetch_go_terms_for_gene(gene) for gene in gene_list}


def analyze_gene_set_go_depths(
    gene_list,
    namespace_depths,
    namespace="BP",
    max_depth=15,
    go_map_builder=build_go_map,
    go_match_function=find_deepest_shared_go_match
):
    gene_list = sorted(list(gene_list))
    go_map = go_map_builder(gene_list)
    results = []
    for gene_a, gene_b in combinations(gene_list, 2):
        depth = go_match_function(gene_a, gene_b, go_map, namespace_depths, namespace, max_depth)
        results.append((gene_a, gene_b, depth))
    return pd.DataFrame(results, columns=["Gene_a", "Gene_b", f"deepest_{namespace.lower()}_depth"])

if __name__ == "__main__":
    
    go_graph = load_go_ontology()
    namespace_depths = compute_namespace_depths(go_graph)

    gene_sets = {
        "proteasome_genes": {
        "PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7",
        "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7"},
        "EIF3_complex": {"EIF3D", "EIF3F", "EIF3H", "EIF3CL"},
        "COPI_complex": {"COPE", "COPB1", "COPA"},
        "Mediator_complex": {"MED19", "MED7", "MED14", "MED24"},
        "IFT_complex": {"IFT122", "WDR19", "IFT43", "IFT140"},
        "CCT_complex": {"TCP1", "CCT3", "CCT5", "CCT8"},
    }

    all_results = []
    for label, genes in gene_sets.items():
        # Compute both BP and CC depths
        df_bp = analyze_gene_set_go_depths(
            gene_list=genes,
            namespace_depths=namespace_depths,
            namespace="BP",
            max_depth=15
        )

        df_cc = analyze_gene_set_go_depths(
            gene_list=genes,
            namespace_depths=namespace_depths,
            namespace="CC",
            max_depth=15
        )

        # Merge on Gene_a and Gene_b
        df_merged = pd.merge(df_bp, df_cc, on=["Gene_a", "Gene_b"])
        df_merged["complex"] = label

        all_results.append(df_merged)

    final_df = pd.concat(all_results, ignore_index=True)

    with pd.option_context('display.max_rows', None):
        print(final_df)

    print("\nSummary by complex (BP):")
    print(final_df.groupby("complex")["deepest_bp_depth"].value_counts())
    print(final_df.groupby("complex")["deepest_bp_depth"].describe())

    print("\nSummary by complex (CC):")
    print(final_df.groupby("complex")["deepest_cc_depth"].value_counts())
    print(final_df.groupby("complex")["deepest_cc_depth"].describe())


    """
    go_graph = load_go_ontology()
    namespace_depths = compute_namespace_depths(go_graph)

    # Proteasome genes
    proteasome_genes = {
        "PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7",
        "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7"
    }

    df_proteasome = analyze_gene_set_go_depths(
        gene_list=proteasome_genes,
        namespace_depths=namespace_depths,
        namespace="BP"
    )
    print(df_proteasome)
    print(df_proteasome['deepest_bp_depth'].value_counts())

    # Example for EMC genes
    EMC_genes = {"EMC2", "EMC4", "EMC8", "EMC10"}
    df_emc = analyze_gene_set_go_depths(EMC_genes, namespace_depths, namespace="BP")
    print(df_emc)

    # Example GO term depth lookup
    go_terms = ['GO:0010595', 'GO:0001525', 'GO:0045050']
    for term in go_terms:
        print(f"{term} → BP depth:", namespace_depths["BP"].get(term, "Not found"))
    """
