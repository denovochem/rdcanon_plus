from importlib.metadata import PackageNotFoundError, version

from rdkit import RDLogger

from rdcanon.main import (
    Graph,
    canon_reaction_smarts,
    canon_smarts,
    debug,
    gen_canon_repl_dict,
    random_smarts,
)
from rdcanon.pubchem_prims import prims
from rdcanon.token_parser import order_token_canon

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

__all__ = [
    "canon_smarts",
    "canon_reaction_smarts",
    "random_smarts",
    "debug",
    "gen_canon_repl_dict",
    "order_token_canon",
    "Graph",
    "prims",
]

try:
    __version__ = version("rdcanon_plus")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
