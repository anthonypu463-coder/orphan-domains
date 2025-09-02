from dataclasses import dataclass
from pathlib import Path
import os

@dataclass(frozen=True)
class ProjectPaths:
    root: Path
    data_root: Path
    afdb: Path
    pfam: Path
    outputs: Path
    docs: Path

    @classmethod
    def from_env(cls, root: str | None = None) -> "ProjectPaths":
        root_path = Path(root or os.getenv("ORPHAN_ROOT", Path.cwd()))
        data_root = Path(os.getenv("DATA_ROOT", root_path / "data"))
        return cls(
            root=root_path,
            data_root=data_root,
            afdb=data_root / "afdb",
            pfam=data_root / "pfam",
            outputs=data_root / "outputs",
            docs=root_path / "docs",
        )
