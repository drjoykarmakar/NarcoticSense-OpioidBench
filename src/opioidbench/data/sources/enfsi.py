"""ENFSI-DWG EI-MS source adapter metadata."""

from opioidbench.data.sources.base import DataSource


class ENFSIDWGMassSpecSource(DataSource):
    """Adapter metadata for locally supplied ENFSI-DWG EI-MS library files."""

    def supported_suffixes(self) -> tuple[str, ...]:
        return (".jdx", ".dx", ".zip")
