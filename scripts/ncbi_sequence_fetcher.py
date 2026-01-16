#!/usr/bin/env python3
"""
Fetch nucleotide or protein sequences from NCBI using common or scientific names and gene symbols.

Examples
--------
Fetch three FASTA sequences for BRCA1 from human and house mouse nucleotide databases:

    python scripts/ncbi_sequence_fetcher.py \\
        --species "Homo sapiens" "house mouse" \\
        --genes BRCA1 \\
        --retmax 3

Fetch protein sequences with an API key and custom output directory:

    python scripts/ncbi_sequence_fetcher.py \\
        --species "Homo sapiens" \\
        --genes TP53 \\
        --db protein \\
        --api-key YOUR_NCBI_KEY \\
        --outdir outputs
"""

from __future__ import annotations

import argparse
import re
import sys
import time
import xml.etree.ElementTree as ET
from typing import Iterable, List, Optional, Tuple
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL_NAME = "ncbi-sequence-fetcher"
DEFAULT_FLANK_BP = 100_000


def slugify(value: str) -> str:
    """Convert a string into a safe filename fragment."""
    normalized = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    return normalized or "item"


def perform_get(endpoint: str, params: dict, delay: float, retries: int) -> str:
    """Perform a GET request against an NCBI endpoint with retries."""
    query = urlencode(params)
    url = f"{NCBI_BASE}{endpoint}?{query}"
    last_error: Optional[BaseException] = None

    for attempt in range(1, retries + 1):
        if delay > 0:
            time.sleep(delay)
        try:
            with urlopen(url) as response:
                return response.read().decode("utf-8")
        except (HTTPError, URLError) as exc:  # noqa: PERF203
            last_error = exc
            if attempt < retries:
                continue
            break

    raise RuntimeError(f"Failed to fetch from NCBI endpoint {endpoint}: {last_error}")


def parse_id_list(xml_text: str) -> List[str]:
    """Extract a list of IDs from an Entrez E-Utilities XML response."""
    root = ET.fromstring(xml_text)
    return [elem.text for elem in root.findall(".//IdList/Id") if elem.text]


def perform_json_get(endpoint: str, params: dict, delay: float, retries: int) -> dict:
    """Perform a GET request and parse JSON."""
    params = {**params, "retmode": "json"}
    raw = perform_get(endpoint, params, delay=delay, retries=retries)
    try:
        import json
    except ImportError as exc:  # pragma: no cover - stdlib guard
        raise RuntimeError("Failed to import json module") from exc
    return json.loads(raw)


def fetch_tax_id(
    organism_name: str, *, email: Optional[str], api_key: Optional[str], delay: float, retries: int
) -> Optional[str]:
    """Look up a taxonomy ID from a common or scientific name."""
    params = {
        "db": "taxonomy",
        "term": f"{organism_name}[All Names]",
        "retmode": "xml",
        "tool": TOOL_NAME,
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    response = perform_get("esearch.fcgi", params, delay=delay, retries=retries)
    ids = parse_id_list(response)
    return ids[0] if ids else None


def fetch_sequence_ids(
    tax_id: str,
    gene: str,
    *,
    db: str,
    retmax: int,
    email: Optional[str],
    api_key: Optional[str],
    delay: float,
    retries: int,
) -> List[str]:
    """Retrieve Entrez IDs for a given gene and organism."""
    query = f"{gene}[Gene] AND txid{tax_id}[Organism:exp]"
    params = {
        "db": db,
        "term": query,
        "retmode": "xml",
        "retmax": str(retmax),
        "tool": TOOL_NAME,
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    response = perform_get("esearch.fcgi", params, delay=delay, retries=retries)
    return parse_id_list(response)


def fetch_gene_info(
    tax_id: str,
    gene: str,
    *,
    email: Optional[str],
    api_key: Optional[str],
    delay: float,
    retries: int,
) -> Optional[Tuple[str, str, int, int]]:
    """
    Fetch gene ID and genomic coordinates (accession, start, stop).

    Returns a tuple of (gene_id, chr_accver, start, stop) or None if not found.
    """
    search_params = {
        "db": "gene",
        "term": f"{gene}[Gene] AND txid{tax_id}[Organism:exp]",
        "tool": TOOL_NAME,
    }
    if email:
        search_params["email"] = email
    if api_key:
        search_params["api_key"] = api_key

    search_json = perform_json_get("esearch.fcgi", search_params, delay=delay, retries=retries)
    ids = search_json.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None

    gene_id = ids[0]
    summary_params = {"db": "gene", "id": gene_id, "tool": TOOL_NAME}
    if email:
        summary_params["email"] = email
    if api_key:
        summary_params["api_key"] = api_key

    summary_json = perform_json_get("esummary.fcgi", summary_params, delay=delay, retries=retries)
    result = summary_json.get("result", {})
    doc = result.get(gene_id, {})
    genomic_info = doc.get("genomicinfo", [])
    if not genomic_info:
        return None

    locus = genomic_info[0]
    return (
        gene_id,
        locus.get("chraccver"),
        int(locus.get("chrstart")),
        int(locus.get("chrstop")),
    )


def find_neighbor_gene_ids(
    tax_id: str,
    chromosome: str,
    center_start: int,
    center_stop: int,
    neighbors: int,
    *,
    email: Optional[str],
    api_key: Optional[str],
    delay: float,
    retries: int,
) -> List[Tuple[str, int, int]]:
    """Find neighbor gene IDs and coordinates around a locus."""
    if neighbors <= 0:
        return []

    span = max(center_stop - center_start, 1)
    window = span * 10 + DEFAULT_FLANK_BP
    region_start = max(0, center_start - window)
    region_stop = center_stop + window

    term = (
        f"{chromosome}[Chromosome] AND {region_start}:{region_stop}[chrpos] "
        f"AND txid{tax_id}[Organism:exp]"
    )
    search_params = {"db": "gene", "term": term, "retmax": 500, "tool": TOOL_NAME}
    if email:
        search_params["email"] = email
    if api_key:
        search_params["api_key"] = api_key

    search_json = perform_json_get("esearch.fcgi", search_params, delay=delay, retries=retries)
    idlist = search_json.get("esearchresult", {}).get("idlist", [])
    if not idlist:
        return []

    summary_params = {"db": "gene", "id": ",".join(idlist), "tool": TOOL_NAME}
    if email:
        summary_params["email"] = email
    if api_key:
        summary_params["api_key"] = api_key
    summary_json = perform_json_get("esummary.fcgi", summary_params, delay=delay, retries=retries)
    result = summary_json.get("result", {})

    entries: List[Tuple[str, int, int]] = []
    for gid in idlist:
        info = result.get(gid, {})
        genomic_info = info.get("genomicinfo", [])
        if not genomic_info:
            continue
        locus = genomic_info[0]
        entries.append((gid, int(locus.get("chrstart")), int(locus.get("chrstop"))))

    entries.sort(key=lambda item: item[1])
    return entries


def fetch_sequences(
    ids: Iterable[str],
    *,
    db: str,
    rettype: str,
    retmode: str,
    email: Optional[str],
    api_key: Optional[str],
    delay: float,
    retries: int,
) -> str:
    """Fetch sequences for a list of IDs from Entrez."""
    id_str = ",".join(ids)
    params = {
        "db": db,
        "id": id_str,
        "rettype": rettype,
        "retmode": retmode,
        "tool": TOOL_NAME,
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    return perform_get("efetch.fcgi", params, delay=delay, retries=retries)


def fetch_genomic_slice(
    accession: str,
    start: int,
    stop: int,
    *,
    email: Optional[str],
    api_key: Optional[str],
    delay: float,
    retries: int,
) -> str:
    """Fetch a genomic slice from nuccore."""
    params = {
        "db": "nuccore",
        "id": accession,
        "seq_start": start + 1,  # NCBI is 1-based inclusive
        "seq_stop": stop + 1,
        "rettype": "fasta",
        "retmode": "text",
        "tool": TOOL_NAME,
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    return perform_get("efetch.fcgi", params, delay=delay, retries=retries)


def fasta_length(fasta_text: str) -> int:
    """Compute sequence length from FASTA text."""
    return sum(len(line.strip()) for line in fasta_text.splitlines() if line and not line.startswith(">"))


def save_sequences(content: str, outdir: str, species: str, gene: str, db: str, suffix: str = "") -> Path:
    """Write retrieved sequences to disk and return the path."""
    Path(outdir).mkdir(parents=True, exist_ok=True)
    suffix_part = f"_{suffix}" if suffix else ""
    filename = f"{slugify(species)}_{slugify(gene)}_{db}{suffix_part}.fasta"
    path = Path(outdir) / filename
    path.write_text(content, encoding="utf-8")
    return path


def process_request(
    species: str,
    gene: str,
    args: argparse.Namespace,
) -> Tuple[Optional[str], Optional[str]]:
    """Resolve taxonomy, fetch IDs, retrieve sequences, and save them."""
    tax_id = fetch_tax_id(
        species,
        email=args.email,
        api_key=args.api_key,
        delay=args.delay,
        retries=args.retries,
    )
    if not tax_id:
        print(f"[WARN] No taxonomy ID found for '{species}'.", file=sys.stderr)
        return None, None

    sequence_ids = fetch_sequence_ids(
        tax_id,
        gene,
        db=args.db,
        retmax=args.retmax,
        email=args.email,
        api_key=args.api_key,
        delay=args.delay,
        retries=args.retries,
    )
    output_path = None
    if sequence_ids:
        sequence_data = fetch_sequences(
            sequence_ids,
            db=args.db,
            rettype=args.rettype,
            retmode=args.retmode,
            email=args.email,
            api_key=args.api_key,
            delay=args.delay,
            retries=args.retries,
        )
        output_path = save_sequences(sequence_data, args.outdir, species, gene, args.db)
        seq_len = fasta_length(sequence_data)
        print(f"[INFO] Primary sequences length: {seq_len} bp")
    else:
        print(
            f"[WARN] No {args.db} IDs found for gene '{gene}' in '{species}' (taxid {tax_id}).",
            file=sys.stderr,
        )

    if args.flank_bp or args.neighbor_genes:
        gene_info = fetch_gene_info(
            tax_id,
            gene,
            email=args.email,
            api_key=args.api_key,
            delay=args.delay,
            retries=args.retries,
        )
        if not gene_info:
            print(
                f"[WARN] Could not resolve genomic coordinates for '{gene}' in '{species}'.",
                file=sys.stderr,
            )
            return tax_id, output_path

        gene_id, chr_acc, start, stop = gene_info
        region_start = max(0, start - args.flank_bp)
        region_stop = stop + args.flank_bp

        if args.neighbor_genes:
            neighbors = find_neighbor_gene_ids(
                tax_id,
                chr_acc,
                start,
                stop,
                args.neighbor_genes,
                email=args.email,
                api_key=args.api_key,
                delay=args.delay,
                retries=args.retries,
            )
            indices = [i for i, item in enumerate(neighbors) if item[0] == gene_id]
            if indices:
                idx = indices[0]
                lower = max(0, idx - args.neighbor_genes)
                upper = min(len(neighbors) - 1, idx + args.neighbor_genes)
                region_start = min(region_start, neighbors[lower][1])
                region_stop = max(region_stop, neighbors[upper][2])
            else:
                print(
                    f"[WARN] Could not find neighboring genes for '{gene}' in '{species}'.",
                    file=sys.stderr,
                )

        try:
            genomic_fasta = fetch_genomic_slice(
                chr_acc,
                region_start,
                region_stop,
                email=args.email,
                api_key=args.api_key,
                delay=args.delay,
                retries=args.retries,
            )
            context_path = save_sequences(
                genomic_fasta,
                args.outdir,
                species,
                gene,
                "genomic_context",
                suffix=f"{region_start}_{region_stop}",
            )
            seq_len = fasta_length(genomic_fasta)
            print(
                f"[OK] Saved genomic context for {species} / {gene} ({seq_len} bp) to {context_path}"
            )
        except RuntimeError as exc:
            print(f"[WARN] Failed to fetch genomic context: {exc}", file=sys.stderr)

    return tax_id, output_path


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="Fetch sequences from NCBI by organism (common or scientific name) and gene name."
    )
    parser.add_argument(
        "-s",
        "--species",
        nargs="+",
        required=True,
        help="One or more common or scientific organism names (quote multi-word names).",
    )
    parser.add_argument(
        "-g",
        "--genes",
        nargs="+",
        required=True,
        help="One or more gene names or symbols to search for.",
    )
    parser.add_argument(
        "--db",
        default="nuccore",
        choices=["nuccore", "protein"],
        help="Entrez database to query.",
    )
    parser.add_argument(
        "--rettype",
        default="fasta",
        help="Entrez rettype passed to efetch (default: fasta).",
    )
    parser.add_argument(
        "--retmode",
        default="text",
        help="Entrez retmode passed to efetch (default: text).",
    )
    parser.add_argument(
        "--retmax",
        type=int,
        default=5,
        help="Maximum number of IDs to fetch per organism/gene combination.",
    )
    parser.add_argument(
        "--flank-bp",
        type=int,
        default=DEFAULT_FLANK_BP,
        help="Genomic flanking size (in bp) upstream and downstream around the gene (default: 100000).",
    )
    parser.add_argument(
        "--neighbor-genes",
        type=int,
        default=0,
        help="Number of upstream/downstream genes to include on either side (default: 0).",
    )
    parser.add_argument(
        "--outdir",
        default="sequence_outputs",
        help="Directory where FASTA files will be written.",
    )
    parser.add_argument(
        "--email",
        help="Contact email to include in NCBI requests (recommended by NCBI).",
    )
    parser.add_argument(
        "--api-key",
        dest="api_key",
        help="Optional NCBI API key for higher rate limits.",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.34,
        help="Delay (in seconds) between requests to respect NCBI rate limits.",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=3,
        help="Number of retries for failed network requests.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    """Entrypoint for the CLI."""
    args = parse_args(argv)
    results = []

    for species in args.species:
        for gene in args.genes:
            try:
                tax_id, path = process_request(species, gene, args)
            except RuntimeError as exc:
                print(f"[ERROR] Failed to process {species} / {gene}: {exc}", file=sys.stderr)
                continue

            results.append((species, gene, tax_id, path))
            if path:
                print(f"[OK] Saved {args.db} sequences for {species} / {gene} to {path}")

    if not results:
        return 1

    unresolved = [item for item in results if not item[3]]
    if unresolved:
        print(
            f"[INFO] Completed with warnings: {len(unresolved)} combinations had no sequences.",
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
