"""Per-haplotype FASTA loading and indel-drift-corrected local coordinates."""

import gzip

_HAPLOTYPE_DRIFT_CACHE = {}


def load_haplotype_fasta(dataset_dir, sample_id, gene, haplotype):
    path = dataset_dir / "individuals" / sample_id / "windows" / gene / f"{sample_id}.{haplotype}.window.fixed.fa"
    lines = path.read_text().splitlines()
    return "".join(line for line in lines if not line.startswith(">"))


def _haplotype_indel_events(dataset_dir, sample_id, gene, haplotype):
    '''(pos_1based, length_delta) for every length-changing variant this haplotype carries in the
    gene's window VCF, sorted by position -- used to replay bcftools consensus's position-shifting
    effect. Handles multi-allelic sites (picks the ALT the phased GT actually points to) and
    symbolic structural-variant alleles: the dataset's consensus-ready VCF is filtered with
    `bcftools view -e \'ALT~"<" && ALT!="<DEL>" && ALT!="<NON_REF>"\'`, so <DEL> entries (using
    their SVLEN/END, not a literal string length) DO reach bcftools consensus and must be counted
    -- e.g. HG02461/TYRP1 carries a real 2,797 bp <DEL> on H1 upstream of the TSS, which a naive
    len(alt)-len(ref) on the literal "<DEL>" text would undercount by ~2.8 kb. Every other symbolic
    type (<DUP>, <INS>, <INS:ME:ALU>, ...) is dropped by that same filter and never reaches the
    consensus sequence, so it's skipped here too.'''
    key = (sample_id, gene, haplotype)
    if key in _HAPLOTYPE_DRIFT_CACHE:
        return _HAPLOTYPE_DRIFT_CACHE[key]

    vcf_path = dataset_dir / "individuals" / sample_id / "windows" / gene / f"{sample_id}.window.vcf.gz"
    hap_index = 0 if haplotype == "H1" else 1
    events = []
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            pos = int(cols[1])
            ref = cols[3]
            alt_alleles = cols[4].split(",")
            info = cols[7]
            gt = cols[9].split(":")[0]
            sep = "|" if "|" in gt else "/"
            alleles = gt.split(sep)
            if len(alleles) != 2:
                continue
            allele_str = alleles[hap_index]
            if allele_str in ("0", "."):
                continue
            allele_num = int(allele_str)
            if not (1 <= allele_num <= len(alt_alleles)):
                continue
            alt = alt_alleles[allele_num - 1]

            if alt.startswith("<"):
                if alt != "<DEL>":
                    continue  # filtered out before bcftools consensus -- never shifts the sequence
                delta = None
                for field in info.split(";"):
                    if field.startswith("SVLEN="):
                        delta = int(field.split("=", 1)[1].split(",")[0])
                        break
                if delta is None:
                    end = None
                    for field in info.split(";"):
                        if field.startswith("END="):
                            end = int(field.split("=", 1)[1])
                            break
                    if end is None:
                        continue
                    delta = -(end - pos)
            else:
                delta = len(alt) - len(ref)

            if delta != 0:
                events.append((pos, delta))
    events.sort()
    _HAPLOTYPE_DRIFT_CACHE[key] = events
    return events


def haplotype_local_idx(dataset_dir, sample_id, gene, haplotype, start_1based, target_pos_1based):
    '''Reference position target_pos_1based -> 0-based local index into this haplotype's own
    fixed FASTA, correcting for this haplotype's own indels between the window start and the
    target.'''
    events = _haplotype_indel_events(dataset_dir, sample_id, gene, haplotype)
    drift = sum(delta for pos, delta in events if pos < target_pos_1based)
    return (target_pos_1based - start_1based) + drift
