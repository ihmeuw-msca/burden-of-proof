"""
Process Class
"""
from pathlib import Path
import shutil

import pickle5 as pkl
import pandas as pd

from espipeline.utils import list_all_files


class Process:
    """
    Process function
    """

    def __init__(self, i_path: Path, o_path: Path):
        self.i_path = Path(i_path)
        self.o_path = Path(o_path)
        self.name = i_path.name

        # load all files
        for f in list_all_files(self.i_path):
            stem, suffix = f.stem, f.suffix
            if suffix == ".pkl":
                setattr(self, stem, pkl.load(open(f, "rb")))
            elif suffix == ".csv":
                setattr(self, stem, pd.read_csv(f))

    def run(self):
        if self.o_path.exists():
            shutil.rmtree(self.o_path)
        shutil.copytree(self.i_path, self.o_path)
