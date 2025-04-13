import gzip

def parse_gaf(file_path, namespace_filter="biological_process"):
    gene2go = {}
    go2genes = {}

    aspect_map = {
        "P": "biological_process",
        "F": "molecular_function",
        "C": "cellular_component"
    }

    open_func = gzip.open if file_path.endswith(".gz") else open

    with open_func(file_path, "rt") as f:
        for line in f:
            if line.startswith("!"):
                continue  # Skip comments

            fields = line.strip().split("\t")
            if len(fields) < 15:
                continue

            #print(fields)
            gene = fields[2]
            #print(gene)
            go_term = fields[4]
            aspect_code = fields[8]
            namespace = aspect_map.get(aspect_code, None)

            if namespace_filter and namespace != namespace_filter:
                continue

            # Update gene2go
            gene2go.setdefault(gene, set()).add(go_term)

            # Update go2genes
            go2genes.setdefault(go_term, set()).add(gene)

    return gene2go, go2genes

gaf_path = "goa_human.gaf.gz"  # or .gz
gene2go, go2genes = parse_gaf(gaf_path, namespace_filter="biological_process")

print('gene 2 go')
print(gene2go)

print(' ')
print('go 2 gene')
# print(go2genes)
