"""Global configuration singleton for GCEN2 runtime settings."""

from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class GCENConfig:
    """Runtime configuration for GCEN2.

    This is a simple dataclass; instantiate it and pass to functions, or
    modify the module-level ``default_config`` to change global defaults.
    """

    n_threads: int = field(default_factory=lambda: int(os.getenv("GCEN2_THREADS", "2")))
    random_seed: int = 42
    verbosity: int = 1  # 0=silent, 1=progress bars, 2=debug
    temp_dir: Optional[str] = None
    figure_dpi: int = 300
    figure_format: str = "pdf"  # pdf | png | svg

    def get_temp_dir(self) -> str:
        if self.temp_dir is None:
            return tempfile.gettempdir()
        return self.temp_dir


#: Module-level default configuration used when no config is explicitly passed.
default_config = GCENConfig()
