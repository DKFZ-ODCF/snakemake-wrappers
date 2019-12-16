__author__ = "Johannes Köster"
__copyright__ = "Copyright 2019, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"
__contributor__= "Johannes Köster, Fritjof Lammers"

from urllib import request
from subprocess import run
from os.path import basename

species = snakemake.params.species.lower()
release = snakemake.params.release
fmt = snakemake.params.fmt
build = snakemake.params.build


def checksum():
    """ Retrieve CHECKSUMS file from ENSEMBL FTP directory.
    and compares to locally computed checksum of downloaded data.
    Returns true on match, exits (code 1) on failure."""

    cksum_url = "{baseurl}/CHECKSUMS".format(baseurl=url.rsplit("/", 1)[0])

    try:
        response = request.urlopen(cksum_url)
    except:
        print("Error: Could not retrieve CHECKSUMS %s" % cksum_url)

    lines = response.read().decode("UTF-8").strip().split("\n")
    for line in lines:
        fields = line.strip().split()
        cksum = int(fields[0])
        filename = fields[2]
        if filename == basename(snakemake.output[0]):
            cksum_local = int(run(["sum", snakemake.output[0]], capture_output=True).stdout.strip().split()[0])
            if cksum_local == cksum:
                print("CHECKSUM OK: %s" % snakemake.output[0])
                return True
            else:
                print("CHECKSUM FAILED: %s" % snakemake.output[0])
                exit(1)
        else:
            continue

suffix = ""
if fmt == "gtf":
    suffix = "gtf.gz"
elif fmt == "gff3":
    suffix = "gff3.gz"


with open(snakemake.output[0], "wb") as out:
    url = "ftp://ftp.ensembl.org/pub/release-{release}/{fmt}/{species}/{species_cap}.{build}.{release}.{suffix}".format(
            release=release,
            build=build,
            species=species,
            fmt=fmt,
            species_cap=species.capitalize(),
            suffix=suffix)
    try:
        r = request.urlopen(url)
    except:
        print("Error: could not retrieve %s" % url)

    out.write(r.read())
    print(url)


checksum()



