import pandas as pd
import requests
import argparse
import time
import os

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--table", required=True)
    parser.add_argument("-o", "--output", default="final_names_report.txt")
    return parser.parse_args()

pdb_mapping_cache = {}
uniprot_name_cache = {}

def get_all_mappings_for_pdb(pdb_id):
    pdb_id = pdb_id.lower()
    if pdb_id in pdb_mapping_cache:
        return pdb_mapping_cache[pdb_id]
    
    # Используем PDBe API (EBI) — он лучше переваривает рибосомы
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}"
    mapping = {}
    try:
        r = requests.get(url, timeout=15)
        if r.status_code == 200:
            data = r.json().get(pdb_id, {}).get('UniProt', {})
            for acc, details in data.items():
                for m in details.get('mappings', []):
                    entity_id = str(m.get('entity_id'))
                    mapping[entity_id] = acc
            pdb_mapping_cache[pdb_id] = mapping
            return mapping
    except:
        pass
    
    pdb_mapping_cache[pdb_id] = {}
    return {}

def get_uniprot_name(acc):
    if not acc or acc == "No Mapping": return None
    if acc in uniprot_name_cache: return uniprot_name_cache[acc]
    
    url = f"https://rest.uniprot.org/uniprotkb/{acc}?fields=protein_name"
    try:
        r = requests.get(url, timeout=10)
        if r.status_code == 200:
            data = r.json()
            name = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value')
            if not name:
                subs = data.get('proteinDescription', {}).get('submissionNames', [])
                if subs: name = subs[0].get('fullName', {}).get('value')
            
            uniprot_name_cache[acc] = name if name else "Unknown Protein"
            return uniprot_name_cache[acc]
    except:
        pass
    return "Name not found"

def extract_pdb_info(text):
    """
    Разбирает строку типа 'Toxoplasma_gondii_mito_9I05_17|Chain'
    Берет строго последние два элемента перед пайпом.
    """
    # 1. Отсекаем всё, что после '|'
    clean_text = str(text).split('|')[0]
    # 2. Дробим по подчеркиваниям
    parts = clean_text.split('_')
    
    if len(parts) >= 2:
        pdb_id = parts[-2]      # Предположительно 9I05
        entity_id = parts[-1]   # Предположительно 17
        
        # Проверяем, что PDB похож на правду (4 символа), а Entity — это число
        if len(pdb_id) == 4 and entity_id.isdigit():
            return pdb_id.upper(), entity_id
            
    return None, None

def main():
    args = parse_args()
    if not os.path.exists(args.table):
        print("Файл не найден!"); return

    df = pd.read_csv(args.table, sep='\t')
    if df.shape[1] < 2:
        df = pd.read_csv(args.table, sep=None, engine='python')
    
    df.columns = [c.strip() for c in df.columns]
    
    print(f"[*] Читаю таблицу. Колонки: {list(df.columns)}")

    with open(args.output, "w", encoding="utf-8") as out:
        for idx, row in df.iterrows():
            node = row.get('Species Tree Node', 'N/A')
            out.write(f"\nNODE: {node}\n")
            print(f"[*] Строка {idx+1}: Node {node}")
            
            for col in ['Genes 1', 'Genes 2']:
                out.write(f"{col}:\n")
                genes_val = str(row.get(col, ''))
                if genes_val == 'nan' or not genes_val.strip(): continue
                
                genes_list = genes_val.split(',')
                for g in genes_list:
                    g = g.strip()
                    pdb_id, ent_id = extract_pdb_info(g)
                    
                    if pdb_id:
                        full_map = get_all_mappings_for_pdb(pdb_id)
                        acc = full_map.get(ent_id)
                        
                        if acc:
                            name = get_uniprot_name(acc)
                            out.write(f"  - {acc}: {name} (от {g})\n")
                            print(f"    [+] Нашел: {pdb_id}_{ent_id} -> {acc}")
                        else:
                            out.write(f"  - NO UNIPROT: {g}\n")
                    else:
                        out.write(f"  - PARSE ERROR: {g}\n")
            
            out.write("-" * 60 + "\n")
            time.sleep(0.05)

    print(f"\n[+] Всё готово. Результаты в {args.output}")

if __name__ == "__main__":
    main()
