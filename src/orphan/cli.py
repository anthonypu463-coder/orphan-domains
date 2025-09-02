import logging, shutil, subprocess
from pathlib import Path
from typing import Optional

import typer, yaml

from .logging_utils import configure_logging
from .config import ProjectPaths
from .ids import IdAllocator

app = typer.Typer(add_completion=False)

@app.command()
def env_check(
    q: bool = typer.Option(False, "-q", "--quiet", help="Reduce log volume (WARNING+)."),
    log_format: str = typer.Option("json", "--log-format", help="json or text"),
    log_level: str = typer.Option("INFO", "--log-level"),
):
    """Report availability and versions of required binaries."""
    configure_logging(level=log_level, json_format=(log_format == "json"), quiet=q)
    log = logging.getLogger("orphan.env")
    paths = ProjectPaths.from_env()
    log.info("paths", extra={
        "data_root": str(paths.data_root),
        "afdb": str(paths.afdb),
        "pfam": str(paths.pfam),
        "outputs": str(paths.outputs),
        "docs": str(paths.docs),
    })
    for tool, version_flag in [("mmseqs", "--version"), ("hmmscan", "-h"), ("foldseek", "--version")]:
        path = shutil.which(tool)
        if not path:
            log.warning("tool_not_found", extra={"tool": tool})
            continue
        try:
            proc = subprocess.run([tool, version_flag], capture_output=True, text=True, check=False)
            first_line = (proc.stdout or proc.stderr).splitlines()[0].strip()
        except Exception as e:
            first_line = f"unavailable ({e.__class__.__name__})"
        log.info("tool", extra={"tool": tool, "path": path, "version": first_line})

@app.command()
def init_versions(
    pfam_version: str = typer.Option("37.4", help="Pfam-A version string."),
    pfam_families: str = typer.Option("~24736", help="Approximate number of families."),
    afdb_release: str = typer.Option("AFDB Swiss-Prot snapshot (2022)", help="AlphaFold DB release description."),
    afdb_count: str = typer.Option("~542000", help="Approximate number of AFDB Swiss-Prot models."),
    uniprot_ref: str = typer.Option("Swiss-Prot (2022)", help="Reference for Swiss-Prot snapshot."),
    force: bool = typer.Option(False, "--force", help="Overwrite if file exists."),
):
    """
    Create docs/dataset_versions.yaml capturing the external dataset baselines.
    Values are editable; this command is idempotent unless --force is provided.
    """
    configure_logging()
    paths = ProjectPaths.from_env()
    paths.docs.mkdir(parents=True, exist_ok=True)
    yml = paths.docs / "dataset_versions.yaml"
    if yml.exists() and not force:
        typer.echo(f"{yml} already exists. Use --force to overwrite.")
        raise typer.Exit(code=0)
    data = {
        "afdb": {"release": afdb_release, "entries": afdb_count},
        "pfam": {"version": pfam_version, "families": pfam_families},
        "uniprot": {"reference": uniprot_ref},
        "notes": "Counts are approximate; pin exact dates/DOIs when available.",
    }
    with yml.open("w") as fh:
        yaml.safe_dump(data, fh, sort_keys=False)
    typer.echo(f"Wrote {yml}")

@app.command()
def show_versions():
    """Print dataset_versions.yaml."""
    paths = ProjectPaths.from_env()
    yml = paths.docs / "dataset_versions.yaml"
    if not yml.exists():
        raise typer.Exit(code=1)
    typer.echo(yml.read_text())

def _infer_dataset_tag(paths: ProjectPaths) -> str:
    yml = paths.docs / "dataset_versions.yaml"
    if not yml.exists():
        return "unversioned"
    data = yaml.safe_load(yml.read_text()) or {}
    pfam_v = (data.get("pfam") or {}).get("version", "NA")
    afdb_r = (data.get("afdb") or {}).get("release", "NA")
    return f"afdb:{afdb_r} | pfam:{pfam_v}"

@app.command()
def mint_domain_id(
    uniprot: str = typer.Option(..., "--uniprot", help="UniProt accession for source protein."),
    start: int = typer.Option(..., "--start", help="Inclusive residue start (1-based)."),
    end: int = typer.Option(..., "--end", help="Inclusive residue end (1-based)."),
    sequence: Optional[str] = typer.Option(None, "--sequence", help="Optional AA sequence for MD5 fingerprint."),
    dataset_tag: Optional[str] = typer.Option(None, "--dataset-tag", help="Override dataset tag; default inferred from dataset_versions.yaml."),
):
    """
    Mint a new Orphan Domain ID (OD#######) and append metadata to registry TSV.
    Registry path: data/outputs/registry/domain_registry.tsv
    """
    configure_logging()
    paths = ProjectPaths.from_env()
    allocator = IdAllocator(paths)
    tag = dataset_tag or _infer_dataset_tag(paths)
    rec = allocator.mint(uniprot_id=uniprot, start=start, end=end, sequence=sequence, dataset_tag=tag)
    typer.echo(f"{rec.domain_id}\t{rec.uniprot_id}\t{rec.start}-{rec.end}\tlen={rec.length}\tmd5={rec.seq_md5}\t{rec.dataset_tag}")
    
if __name__ == "__main__":
    app()
