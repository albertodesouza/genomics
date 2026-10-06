import random
chrom_lengths = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895
}
random.seed(42)
with open("random_regions.bed", "w") as f:
    for chrom, length in chrom_lengths.items():
        starts = sorted(random.sample(range(1, length - 1000), 100))
        for s in starts:
            f.write(f"{chrom}\t{s}\t{s+1000}\n")
