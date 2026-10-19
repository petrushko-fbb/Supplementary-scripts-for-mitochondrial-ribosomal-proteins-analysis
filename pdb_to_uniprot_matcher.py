import argparse
import time
import pandas as pd
import requests


def parse_args():
    parser = argparse.ArgumentParser(description="Map PDB entity IDs to UniProt accessions via RCSB PDB API.")
    parser.add_argument("-t", "--table", required=True, help="Path to the input table file")
    parser.add_argument("-o", "--output", required=True, help="Path to the output TSV file")
    return parser.parse_args()


def get_uniprot_from_pdb(pdb_id, entity_id):
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id.upper()}/{entity_id}"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            uniprot_ids = data.get("rcsb_polymer_entity_container_identifiers", {}).get("uniprot_ids", [])
            if uniprot_ids:
                return uniprot_ids[0]
    except requests.RequestException:
        pass
    return None


def main():
    args = parse_args()

    try:
        df = pd.read_csv(args.table, sep=None, engine="python")
    except Exception as e:
        print(f"Error reading input table: {e}")
        return

    cache = {}
    print("[*] Starting mapping via RCSB PDB API...")

    for idx, row in df.iterrows():
        for col in df.columns[1:]:
            val = str(row[col])
            if val == "nan" or not val.strip():
                continue

            parts = [p.strip() for p in val.split(",")]
            new_parts = []

            for p in parts:
                raw_id = p.split("|")[0].strip()

                if "_" in raw_id:
                    if raw_id not in cache:
                        sub_parts = raw_id.split("_")
                        pdb_id = sub_parts[0]
                        entity_id = sub_parts[1]

                        if entity_id.isdigit():
                            print(f"  [?] Querying PDB: {pdb_id}_{entity_id}...", end=" ", flush=True)
                            u_id = get_uniprot_from_pdb(pdb_id, entity_id)
                            cache[raw_id] = u_id
                            print(f"Found: {u_id}" if u_id else "Not found")
                        else:
                            cache[raw_id] = None

                    found_acc = cache.get(raw_id)
                    new_parts.append(found_acc if found_acc else raw_id)
                else:
                    new_parts.append(raw_id)

            df.at[idx, col] = ", ".join(new_parts)

        time.sleep(0.1)

    df.to_csv(args.output, sep="\t", index=False)
    print(f"[+] Result saved to: {args.output}")


if __name__ == "__main__":
    main()
