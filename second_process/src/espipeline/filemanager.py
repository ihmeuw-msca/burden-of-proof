"""
File manager: organize file paths
"""
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

import pandas as pd


@dataclass
class FileManager:
    """
    Manager files for all risk outcome pair in the folder
    """

    i_path: Path
    o_path: Path
    pairs: List[str] = field(init=False)
    pair_paths: Dict[str, Path] = field(init=False)
    pair_types: List[str] = field(init=False)
    pairs_by_type: Dict[str, List[str]] = field(init=False)

    def __post_init__(self):
        self.i_path = Path(self.i_path)
        self.o_path = Path(self.o_path)
        if not self.i_path.exists():
            raise FileNotFoundError(str(self.i_path))
        if not self.o_path.exists():
            self.o_path.mkdir()

        self.pairs = [path.name
                      for path in self.i_path.iterdir() if path.is_dir()]
        self.pair_paths = {pair: self.i_path / pair for pair in self.pairs}
        self.pairs_by_type = defaultdict(list)
        self.sort_pairs()
        self.pair_types = list(self.pairs_by_type.keys())

    def sort_pairs(self):
        for pair, path in self.pair_paths.items():
            meta = path / "risk_cause_metadata.csv"
            if not meta.exists():
                raise FileNotFoundError(str(meta))
            df_meta = pd.read_csv(meta)
            self.pairs_by_type[df_meta.risk_type[0]].append(pair)

    def __getitem__(self, pair: str) -> Path:
        return self.pair_paths[pair]

    def __repr__(self) -> str:
        return f"{type(self).__name__}(num_pairs={len(self.pairs)})"
