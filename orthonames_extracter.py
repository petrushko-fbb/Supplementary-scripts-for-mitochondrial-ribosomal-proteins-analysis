import argparse
import os
import pandas as pd
from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(description="Extract full FASTA sequence descriptions by Orthogroup.")
    parser.add_argument("-t", "--table", required=True, help="Path to the input table file (TSV/CSV)")
    parser.add_argument("-d", "--dir", required=True, help="Path to the directory containing FASTA files")
    parser.add_argument("-o", "--output", required=True, help="Path to the output text file")
    return parser.parse_args()


def load_full_names(folder_path):
    name_map = {}
    print(f"[*] Scanning FASTA files in {folder_path}...")

    for file in os.listdir(folder_path):
        if file.endswith((".fasta", ".fa", ".faa")):
            path = os.path.join(folder_path, file)
            for record in SeqIO.parse(path, "fasta"):
                name_map[record.id] = record.description

    print(f"[+] Loaded {len(name_map)} unique sequence headers.")
    return name_map


def main():
    args = parse_args()

    try:
        df = pd.read_csv(args.table, sep=None, engine="python")
    except Exception as e:
        print(f"Error reading input table: {e}")
        return

    full_names_dict = load_full_names(args.dir)

    with open(args.output, "w", encoding="utf-8") as f:
        for _, row in df.iterrows():
            og_id = str(row.iloc[0])

            ids_in_row = []
            for col in df.columns[1:]:
                cell = str(row[col])
                if cell != "nan" and cell.strip():
                    sub_ids = [i.strip() for i in cell.split(",")]
                    ids_in_row.extend(sub_ids)

            if ids_in_row:
                f.write(f"{og_id}\n")
                for q_id in ids_in_row:
                    full_name = full_names_dict.get(q_id, f"{q_id} (NOT FOUND IN FILES)")
                    f.write(f"{full_name}\n")
                f.write("\n")

    print(f"[+] Full names saved to: {args.output}")


if __name__ == "__main__":
    main()
