"""One variant in one cell, as read from the input."""

from dataclasses import dataclass, field


@dataclass
class Variant:
    cell: str
    chrom: str
    pos: str
    ref: str
    alt: str
    qual: str = ""
    info: str = ""
    format: str = ""
    sample: str = ""

    # which database columns the input actually had, not what they contained
    databases: set = field(default_factory=set)
    dbsnp: str = ""
    gnomad_af: list = field(default_factory=list)
    rediportal: str = ""
    cosmic: str = ""

    features: dict = field(default_factory=dict)

    @property
    def key(self):
        return self.chrom, self.pos, self.ref, self.alt
