def parse_fasta(filename):
    seqs = {}
    name = None
    chunks = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name:
            seqs[name] = "".join(chunks)
    return seqs

fasta_file = "C:\\bio_project\\test.fa.txt"

fasta = parse_fasta(fasta_file)

print("✅ Sequences in file:", list(fasta.keys()))
print("✅ Length of chr22:", len(fasta['chr22']))


import pandas as pd
# Step 1: Create a tiny VCF file
vcf_content = """##fileformat=VCFv4.2
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO
22  5   .   A   G   60  PASS    DP=20;AF=0.5
22  10  .   C   T   50  PASS    DP=30;AF=0.2
22  15  .   G   A   15  q10     DP=5;AF=0.1
"""

with open("C:/bio_project/test.vcf", "w") as f:
    f.write(vcf_content)

print("Test VCF file created at C:/bio_project/test.vcf")

def parse_vcf(filename):
    headers = []
    rows = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("##"):   # skip metadata
                continue
            if line.startswith("#CHROM"):
                headers = line.strip().split()
                continue
            parts = line.strip().split()
            rows.append(parts)
    return pd.DataFrame(rows, columns=headers)
# Step 3: Parse the tiny VCF
vcf_file = "C:/bio_project/test.vcf"
vcf_df = parse_vcf(vcf_file)

# Step 4: Show the data
print("\n VCF contents (as DataFrame):")
print(vcf_df)

#to parse DP and AF into seperate columns
def parse_info(info_str):
    fields = {}
    for item in info_str.split(";"):
        if "=" in item:
            k,v = item.split("=")
            fields[k] = v
    return fields

vcf_df["DP"] = vcf_df["INFO"].apply(lambda x: int(parse_info(x).get("DP", 0)))
vcf_df["AF"] = vcf_df["INFO"].apply(lambda x: float(parse_info(x).get("AF", 0)))
print(vcf_df[["POS","REF","ALT","DP","AF"]])


#to filter variants
filtered_vcf = vcf_df[(vcf_df["QUAL"].astype(float) >= 30) & 
                      (vcf_df["DP"] >= 10) & 
                      (vcf_df["FILTER"] == "PASS")]
print(filtered_vcf)


#variant classification
def classify_variant(ref, alt):
    if len(ref)==1 and len(alt)==1:
        if {ref,alt} in [{"A","G"},{"C","T"}]:
            return "Transition"
        else:
            return "Transversion"
    return "INDEL"

vcf_df["Type"] = vcf_df.apply(lambda row: classify_variant(row["REF"], row["ALT"]), axis=1)
print(vcf_df[["POS","REF","ALT","Type"]])

#Plot visualization by using matplotlib
import matplotlib.pyplot as plt

vcf_df["QUAL"].astype(float).hist(bins=10)
plt.title("QUAL Distribution")
plt.show()

#Validation by linking variants to FASTA
chrom_seq = fasta["chr22"]   # from FASTA
def check_ref(row):
    pos = int(row["POS"]) - 1   # FASTA is 0-based
    return chrom_seq[pos] == row["REF"]

vcf_df["Ref_match"] = vcf_df.apply(check_ref, axis=1)
print(vcf_df[["POS","REF","Ref_match"]])

