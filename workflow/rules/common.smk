import os


def pfx(path, suffix=None):
    """Return the fileset *prefix* a tool wants from a Snakemake-tracked path.

    plink/plink2/GenotypeHarmonizer/etc. take a basename that implies sibling
    files, while Snakemake tracks individual files. With `suffix` (the dotted
    extension, e.g. ".pgen" or ".dose.vcf.gz") it strips exactly that suffix --
    matching the shell `${v%.suffix}` idiom and handling multi-part suffixes.
    Without `suffix` it drops only the final extension.
    """
    s = str(path)
    if suffix is None:
        return os.path.splitext(s)[0]
    assert s.endswith(suffix), f"{repr(s)} does not end with {repr(suffix)}"
    return s[: -len(suffix)]
