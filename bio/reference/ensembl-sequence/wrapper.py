__author__ = "Johannes Köster"
__copyright__ = "Copyright 2019, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from urllib import request
from subprocess import run
from os.path import basename

species = snakemake.params.species.lower()
release = snakemake.params.release
build = snakemake.params.build

suffixes = ""
datatype = snakemake.params.get("datatype", "")


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


if datatype == "dna":
    suffixes = ["dna.primary_assembly.fa.gz", "dna.toplevel.fa.gz"]
elif datatype == "cdna":
    suffixes = ["cdna.all.fa.gz"]
elif datatype == "cds":
    suffixes = ["cds.all.fa.gz"]
elif datatype == "ncrna":
    suffixes = ["ncrna.fa.gz"]
elif datatype == "pep":
    suffixes = ["pep.all.fa.gz"]
else:
    raise ValueError("invalid datatype, must be one of dna, cdna, cds, ncrna, pep")

success = False

with open(snakemake.output[0], "wb") as out:
    for suffix in suffixes:
        url = "ftp://ftp.ensembl.org/pub/release-{release}/fasta/{species}/{datatype}/{species_cap}.{build}.{suffix}".format(
            release=release,
            species=species,
            datatype=datatype,
            build=build,
            suffix=suffix,
            species_cap=species.capitalize())
        try:
            response = request.urlopen(url)
            print(f"Trying to download: {url}")
        except:
            continue

    out.write(response.read())
    success = True
    print(url)


if not success:
    raise ValueError(
        "Requested sequence does not seem to exist on ensembl FTP servers (url {})".format(
            url
        )
    )

checksum()

