import sys
import csv
# Input/output files
input_csv = "/work/rr151/SamAIM2neuro_data/subtype_lib.csv"
sam_file = "/work/rr151/SamAIM2neuro_data/SJR_subtype_lib_final.sam"
gtf_file = "/work/rr151/SamAIM2neuro_data/gencode.v47.annotation.gtf"
output_file = "/work/rr151/SamAIM2neuro_data/SJR_subtype_lib_final_metadata.tsv"
output_fasta = '/work/rr151/SamAIM2neuro_data/SJR_subtype_lib_final.fasta'
# Step 1: Load guide sequences
guides = {}
with open(input_csv, 'r', newline='') as csvfile, open(output_fasta, 'w') as fastafile:
    reader = csv.DictReader(csvfile)
    reader.fieldnames = [field.strip().replace('\ufeff', '') for field in reader.fieldnames]
    for row in reader:
        gene = row["gene"].strip()
        spacer = row["protospacer"].strip()
        guides[gene] = spacer
        fastafile.write(f'>{gene}\n{spacer}\n')
# Step 2: Parse GTF file to build TSS lookup
tss_lookup = {}  # key = gene_name.upper(), value = (chr, tss_start, tss_end, gene_id)
with open(gtf_file, "r") as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        chrom, _, feature_type, start, end, _, strand, _, attributes = parts
        if feature_type != "gene":
            continue
        # Extract gene_id and gene_name from attributes field
        attr_dict = {}
        for attr in attributes.strip().split(";"):
            if attr.strip():
                key_value = attr.strip().split(" ", 1)
                if len(key_value) == 2:
                    key, value = key_value
                    attr_dict[key] = value.replace('"', '').strip()
        gene_id = attr_dict.get("gene_id", "")
        gene_name = attr_dict.get("gene_name", "").upper()
        if not gene_id or not gene_name:
            continue
        # Determine TSS
        start, end = int(start), int(end)
        if strand == "+":
            tss_start, tss_end = start, start + 1
        else:
            tss_start, tss_end = end - 1, end
        tss_lookup[gene_name] = (chrom, tss_start, tss_end, gene_id)
# Step 3: Parse SAM file and build metadata
metadata = {}
with open(sam_file, "r") as f:
    for line in f:
        if line.startswith("@"):
            continue
        parts = line.strip().split("\t")
        guide_id = parts[0]
        guide_chr = parts[2]
        guide_start = int(parts[3])
        flag = int(parts[1])
        strand = "-" if flag & 16 else "+"
        seq = parts[9]
        pam = seq[-3:] if len(seq) >= 3 else ""
        guide_end = guide_start + len(seq) - 1
        gene_symbol = guide_id.rsplit('.', 1)[0].upper()
        tss = tss_lookup.get(gene_symbol, ("", "", "", ""))
        metadata[guide_id] = {
            "guide_id": guide_id,
            "spacer": guides.get(guide_id, ""),
            "targeting": "TRUE",
            "type": "targeting",
            "guide_chr": guide_chr,
            "guide_start": guide_start,
            "guide_end": guide_end,
            "strand": strand,
            "pam": pam,
            "genomic_element": "tss",
            "intended_target_name": tss[3],  # Ensembl gene ID
            "intended_target_chr": tss[0],
            "intended_target_start": tss[1],
            "intended_target_end": tss[2],
            "putative_target_genes": "",
            "reporter": "",
            "imperfect": ""
        }
# Step 4: Write output
headers = [
    "guide_id", "spacer", "targeting", "type", "guide_chr",
    "guide_start", "guide_end", "strand", "pam", "genomic_element",
    "intended_target_name", "intended_target_chr", "intended_target_start",
    "intended_target_end", "putative_target_genes", "reporter", "imperfect"
]
with open(output_file, "w", newline="") as out_f:
    writer = csv.DictWriter(out_f, fieldnames=headers, delimiter="\t")
    writer.writeheader()
    for guide_id in guides:
        gene_symbol = guide_id.rsplit('.', 1)[0].upper()
        tss = tss_lookup.get(gene_symbol, ("", "", "", ""))
        row = metadata.get(guide_id, {
            "guide_id": guide_id,
            "spacer": guides[guide_id],
            "targeting": "TRUE",
            "type": "targeting",
            "guide_chr": "",
            "guide_start": "",
            "guide_end": "",
            "strand": "",
            "pam": "",
            "genomic_element": "tss",
            "intended_target_name": tss[3],
            "intended_target_chr": tss[0],
            "intended_target_start": tss[1],
            "intended_target_end": tss[2],
            "putative_target_genes": "",
            "reporter": "",
            "imperfect": ""
        })
        writer.writerow(row)
