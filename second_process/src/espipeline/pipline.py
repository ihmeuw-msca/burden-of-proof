"""
Pipeline Class
"""
from typing import Callable, List

import numpy as np

from espipeline.filemanager import FileManager
from espipeline.utils import list_all_files


class Pipeline:
    """
    Main Pipeline class
    """

    def __init__(self,
                 name: str,
                 fm: FileManager,
                 process_constructor: Callable,
                 pairs: List[str] = None):
        self.name = name
        self.pairs = fm.pairs_by_type[self.name] if pairs is None else pairs
        self.i_pair_paths = {
            pair: fm.pair_paths[pair]
            for pair in self.pairs
        }
        self.o_pair_paths = {
            pair: fm.o_path / pair
            for pair in self.pairs
        }
        self.process_constructor = process_constructor

    @property
    def num_pairs(self) -> int:
        return len(self.pairs)

    def get_process(self, pair: str) -> "Process":
        return self.process_constructor(self.i_pair_paths[pair],
                                        self.o_pair_paths[pair])

    def run(self):
        for pair in self.pairs:
            i_path = self.i_pair_paths[pair]
            o_path = self.o_pair_paths[pair]
            i_time = max(f.stat().st_mtime for f in list_all_files(i_path))
            if (not o_path.exists() or o_path.stat().st_size == 0):
                o_time = -np.inf
            else:
                o_time = max(f.stat().st_mtime for f in list_all_files(o_path))
            if i_time > o_time:
                print(pair)
                process = self.get_process(pair)
                process.run()

    def __repr__(self) -> str:
        return (f"{type(self).__name__}(name={self.name}, "
                f"num_pairs={self.num_pairs})")
