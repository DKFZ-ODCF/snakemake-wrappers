__author__ = "Johannes Köster"
__copyright__ = "Copyright 2019, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

from urllib import request
from io import StringIO
from subprocess import run
from os.path import basename


species = snakemake.params.species.lower()
release = int(snakemake.params.release)
fmt = snakemake.params.fmt
build = snakemake.params.build
flavor = snakemake.params.get("flavor", "")

branch = ""
if release >= 81 and build == "GRCh37":
    # use the special grch37 branch for new releases
    branch = "grch37/"

if flavor:
    flavor += "."

log = snakemake.log_fmt_shell(stdout=False, stderr=True)


def checksum():
    lines = r.read().decode("UTF-8").strip().split("\n")
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

r = StringIO()

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
    success = True
    print(url)


cksum_url = "{baseurl}/CHECKSUMS".format(baseurl=url.rsplit("/", 1)[0])

try:
    r = request.urlopen(cksum_url)
except:
    print("Error: Could not retrieve CHECKSUMS %s" % cksum_url)

checksum()