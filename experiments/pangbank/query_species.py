"""Look up PanGBank pangenomes (all collections and releases) for GTDB species and their genera.
usage: query_species.py species.txt species_hits.tsv genus_hits.tsv"""
import json, sys, urllib.parse, urllib.request
API = "https://pangbank-api.genoscope.cns.fr"
def get(path, **q):
    with urllib.request.urlopen(f"{API}{path}?{urllib.parse.urlencode(q)}", timeout=120) as r:
        return json.load(r)
def row(h):
    cr = h["collection_release"]
    taxa = {t["rank"]: t["name"] for t in cr.get("taxonomy", {}).get("taxa", [])} if False else {}
    return [cr.get("collection_name"), cr["version"], cr["taxonomy_source"]["version"], h["id"], h["name"],
            h["genome_count"], h["family_count"], h["persistent_family_count"], h["shell_family_count"],
            h["cloud_family_count"], h["module_count"], round(h["mean_completeness"], 1), h["file_md5sum"]]
hdr = "collection\trelease\tgtdb\tpangenome_id\tname\tgenomes\tfamilies\tpersistent\tshell\tcloud\tmodules\tmean_compl\tmd5"
so, go = open(sys.argv[2], "w"), open(sys.argv[3], "w")
so.write("species\t" + hdr + "\n"); go.write("query_genus\t" + hdr + "\n")
seen_genus = set()
for sp in filter(None, map(str.strip, open(sys.argv[1]))):
    genus = sp.split(" ")[0]
    hits = get("/pangenomes/", taxon_name=sp, substring_taxon_match="true", only_latest_release="false", limit=100)
    hits = [h for h in hits if h["name"] == "s__" + sp.replace(" ", "_")]
    for h in hits: so.write("\t".join(map(str, [sp] + row(h))) + "\n")
    if not hits: so.write(f"{sp}\tNA\n")
    if genus in seen_genus: continue
    seen_genus.add(genus)
    gh = get("/pangenomes/", taxon_name="g__" + genus, substring_taxon_match="false", only_latest_release="false", limit=100)
    for h in gh: go.write("\t".join(map(str, [genus] + row(h))) + "\n")
    if not gh: go.write(f"{genus}\tNA\n")
