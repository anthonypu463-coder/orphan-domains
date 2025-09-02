from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
import csv, fcntl, hashlib, json, logging, time

from .config import ProjectPaths

ID_PREFIX = "OD"
ID_WIDTH = 7

def format_domain_id(n: int) -> str:
    return f"{ID_PREFIX}{n:0{ID_WIDTH}d}"

def parse_domain_id(s: str) -> int:
    if not (s.startswith(ID_PREFIX) and s[len(ID_PREFIX):].isdigit()):
        raise ValueError(f"Invalid domain id: {s}")
    return int(s[len(ID_PREFIX):])

@dataclass(frozen=True)
class DomainRecord:
    domain_id: str
    uniprot_id: str
    start: int
    end: int
    seq_md5: str
    length: int
    dataset_tag: str
    created_at: str  # ISO 8601

class IdAllocator:
    """
    File-backed monotonic counter with advisory locking.
    Ensures stable, unique IDs across sessions.
    """
    def __init__(self, paths: ProjectPaths):
        self.paths = paths
        self.reg_dir = paths.outputs / "registry"
        self.reg_dir.mkdir(parents=True, exist_ok=True)
        self.counter_file = self.reg_dir / "domain_id_counter.json"
        self.registry_tsv = self.reg_dir / "domain_registry.tsv"

    def _reserve_next_int(self) -> int:
        self.counter_file.touch(exist_ok=True)
        with self.counter_file.open("r+") as fh:
            fcntl.flock(fh, fcntl.LOCK_EX)
            try:
                try:
                    data = json.load(fh)
                except json.JSONDecodeError:
                    data = {"next": 1}
                n = int(data.get("next", 1))
                data["next"] = n + 1
                fh.seek(0); fh.truncate(); json.dump(data, fh)
            finally:
                fcntl.flock(fh, fcntl.LOCK_UN)
        return n

    def mint(self, *, uniprot_id: str, start: int, end: int,
             sequence: Optional[str], dataset_tag: str) -> DomainRecord:
        if start < 1 or end < start:
            raise ValueError("Invalid coordinates: start must be ≥1 and ≤ end.")
        seq_md5 = hashlib.md5((sequence or "").encode("utf-8")).hexdigest()
        domain_int = self._reserve_next_int()
        domain_id = format_domain_id(domain_int)
        rec = DomainRecord(
            domain_id=domain_id,
            uniprot_id=uniprot_id,
            start=start,
            end=end,
            seq_md5=seq_md5,
            length=(end - start + 1),
            dataset_tag=dataset_tag,
            created_at=time.strftime("%Y-%m-%dT%H:%M:%S"),
        )
        header = list(DomainRecord.__annotations__.keys())
        new_file = not self.registry_tsv.exists()
        with self.registry_tsv.open("a", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
            if new_file:
                w.writeheader()
            w.writerow(rec.__dict__)
        logging.getLogger("orphan.ids").info("minted_domain_id", extra={"domain_id": rec.domain_id})
        return rec
