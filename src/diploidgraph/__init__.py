"""Diploid sequence graphs from two-terminal DAG sequence graphs."""

from .decompose import DecompNode, Decomposition, decompose
from .diploid import DiploidGraph, NonSeriesParallelError, build_diploid, naive_pair_count
from .seqgraph import NotTwoTerminalDAG, SequenceGraph

__all__ = [
    "SequenceGraph",
    "NotTwoTerminalDAG",
    "DecompNode",
    "Decomposition",
    "decompose",
    "DiploidGraph",
    "NonSeriesParallelError",
    "build_diploid",
    "naive_pair_count",
]
__version__ = "0.1.0"
