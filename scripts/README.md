# NCBI sequence fetcher

This repository now includes a small utility for downloading sequences from NCBI using common or scientific names and gene symbols.

## Requirements

- Python 3.8+ (standard library only; no external dependencies).
- An internet connection.
- An optional NCBI API key if you need higher rate limits.

## Usage

The script lives at `scripts/ncbi_sequence_fetcher.py`. Run it from the repository root:

```bash
python scripts/ncbi_sequence_fetcher.py \
  --species "Homo sapiens" "house mouse" \
  --genes BRCA1 \
  --retmax 3 \
  --outdir sequence_outputs
```

Key options:

- `--species`: One or more common or scientific names (quote multi-word names).
- `--genes`: One or more gene names or symbols.
- `--db`: NCBI database to query (`nuccore` or `protein`). Defaults to `nuccore`.
- `--retmax`: Maximum IDs to fetch per species/gene combination. Defaults to 5.
- `--rettype`/`--retmode`: Passed to `efetch` (defaults `fasta`/`text`).
- `--email`: Contact email sent to NCBI (recommended).
- `--api-key`: NCBI API key (optional but useful for higher throughput).
- `--delay`: Delay between requests to respect NCBI limits (default `0.34` seconds).
- `--flank-bp`: Fetch genomic context upstream/downstream of the gene in base pairs (default `100000` bp).
- `--neighbor-genes`: Include N upstream/downstream genes to define the context window (default `0`).

Outputs are written to the `sequence_outputs` directory by default, with one FASTA file per species/gene combination.
