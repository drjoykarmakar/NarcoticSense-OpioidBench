"""SWGDRUG EI-MS source adapter."""

from opioidbench.data.sources.base import DataSource


class SWGDRUGMassSpecSource(DataSource):
    """Adapter metadata for SWGDRUG EI-MS JCAMP/HPJ inputs."""

    def supported_suffixes(self) -> tuple[str, ...]:
        return (".jdx", ".dx", ".hpj", ".txt")
