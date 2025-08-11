import sys

def parse_fasta(fasta_path):

    """
	Parse a FASTA file into a dict: {chr: seq}
	"""

    with open(fasta_path) as f:
        contents = [line.rstrip('\n') for line in f]

    chrs = contents[::2]
    seqs = contents[1::2]

    fasta = {}

    for i, (chr_, seq) in enumerate(zip(chrs, seqs)):
        if not chr_.startswith(">") or seq.startswith(">"):
            line = i * 2
            raise ValueError(
				f"Wrong FASTA contents in lines {line}-{line + 1}:\n"
				f"            {line}: Awaited chr name ('>chr...'), got \t{chr_[:10]}...\n"
				f"            {line+1}: Awaited chr sequence ('string'), got \t{seq[:10]}...")
            fasta[chr_.lstrip(">")] = seq.upper()
    return fasta

def parse_spike_table(spike_table_path, parced_fasta):

	"""
	Parse the spike table into a dict of mutation dicts 
	{type: 
		{chr: chr_name:chr, 
		 pos : position:int, 
		 abundance : abundance:float, 
		 extra : extra:chr}}.
	"""

	spikes = []
	with open(spike_table_path) as f:
		for i, (spike) in enumerate(f):
			line = spike.strip().split("\t")
			if not line[0] in parced_fasta.keys():
				raise ValueError(f"In line {i} of spikes.tsv: {line} is not a chromosome name present in input FASTA")
			chr_n, pos, abundance, typ, *extra = line[:5]
			spikes.append({
				"chr": chr_n,
				"pos": int(pos),
				"type": typ,
				"abundance": float(abundance),
				"extra": extra[0] if extra else ""})
	return spikes

def build_spike_dict(parsed_fasta, parced_spikes):

    """
	Build a dict: {chr: {pos: [spike_dict, ...]}} for fast lookup.
	"""

    spike_dict = {}

    for chr_ in parced_fasta.keys:
        spike_dict[chr_]
    return spike_dict

def apply_spikes(seq, spikes_at_pos):
    """Apply spikes at each position to the sequence."""
    # This is the core algorithm: for now, just a placeholder
    # Should handle SNP, INS, DEL, MNV, skip redundant
    # For Cython, use for-loops and indexing, avoid dynamic Python features
    # seq is a list of chars!
    new_seq = seq[:]  # Copy
    # ... implement mutation logic here ...
    return new_seq

def main(fasta_path, spike_table_path, out_prefix):
    seqs = parse_fasta(fasta_path)
    spikes = parse_spike_table(spike_table_path)
    spike_dict = build_spike_dict(spikes)

    # For each chromosome
    for chr_n, seq in seqs.items():
        mutated_seqs = []
        # Determine how many times to spike each (abundance)
        # For each abundance, make a new variant
        if chr_n in spike_dict:
            for spike_list in spike_dict[chr_n].values():
                # For each spike, apply mutation as needed
                # For now just pass
                pass

        # Output mutated sequences to a new FASTA (for now, just write original)
        with open(f"{out_prefix}_{chr_n}.fa", "w") as f:
            f.write(f">{chr_n}\n")
            f.write("".join(seq) + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} <reference.fasta> <spike_table.tsv> <output_prefix>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
