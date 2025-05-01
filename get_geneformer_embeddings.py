import pandas as pd
import mygene
import pickle
import sys

genelist = sys.argv[1]
emb_csv = sys.argv[2]

def get_gformer_dict(genelist, emb_csv):
    mg = mygene.MyGeneInfo()
    entrez_ids = pd.read_csv(genelist, header=None)[0].to_list()
    gf=pd.read_csv(emb_csv, header=None)
    
    gnames = []
    
    result = mg.querymany(entrez_ids, scopes="entrezgene", fields="symbol", species="human")
    
    entrez_to_gene = {}
    for r in result:
        if 'symbol' in r:
            entrez_to_gene[r['query']] = r['symbol']
    
    for entrez_id, gene_name in entrez_to_gene.items():
        print(f"Entrez ID: {entrez_id}, Gene Symbol: {gene_name}")
        gnames.append(gene_name)
        
    gformer=dict()
    for _, row in gf.iterrows():
        gene = gnames[_]
        r = row.to_list()
        gformer[gene] = r

    return gformer

def save_dict_to_pkl(dictionary, output_file):
    """
    Saves a dictionary to a pickle (.pkl) file.

    Parameters:
    - dictionary: dict
        The dictionary you want to save.
    - output_file: str
        The path to the output .pkl file.
    """
    with open(output_file, 'wb') as f:
        pickle.dump(dictionary, f, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == "__main__":
    print('Getting embedding dictionary')
    embedding_dict = get_gformer_dict(genelist, emb_csv)
    print('Writing dictionary to pickle')
    save_dict_to_pkl(embedding_dict, f'{emb_csv[:-4]}.pkl')
    